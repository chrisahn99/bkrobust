"""Driver for the Axis Robustness survival-curve sweep.

Generates the flip-arm and tiered-arm instance populations, samples both
corruption processes, and writes the four output files under
``results/axis_robustness/`` described in the task brief and
``results/axis_robustness/PREREGISTRATION.md``. See
:mod:`bkrobust.robustness.survival` for the corruption/scoring logic; this
module is orchestration and I/O only.

Usage::

    PYTHONPATH=src python3 -m bkrobust.robustness.run_survival --mode smoke
    PYTHONPATH=src python3 -m bkrobust.robustness.run_survival --mode full
    PYTHONPATH=src python3 -m bkrobust.robustness.run_survival --mode determinism --out-dir /tmp/x

No global RNG is touched anywhere in this module either -- every draw is
inside :mod:`bkrobust.robustness.survival`, which threads explicit
``np.random.Generator`` objects throughout.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

from bkrobust.core.resultsio import ResultWriter, write_manifest
from bkrobust.robustness import survival as sv

SAMPLES_FIELDS = [
    "instance_id",
    "arm",
    "gate",
    "coverage",
    "base_wrongness",
    "n_tiers",
    "d",
    "rep",
    "flip_rate",
    "corruption_rate",
    "seed",
    "status",
    "survived",
    "symdiff_proxy_not_distance",
]

CURVES_FIELDS = [
    "instance_id",
    "arm",
    "gate",
    "component_size",
    "separation",
    "coverage",
    "base_wrongness",
    "n_tiers",
    "d",
    "n_eval",
    "n_contradictory",
    "S",
    "S_contra_as_fail",
]

INSTANCES_FIELDS = [
    "instance_id",
    "base_instance_id",
    "arm",
    "gate",
    "component_size",
    "separation",
    "coverage",
    "base_wrongness",
    "n_tiers",
    "seed",
    "x",
    "y",
    "n_k",
    "shd_truth",
    "shd_cpdag",
    "undirected_fraction",
    "z_size",
    "z_star",
    "s0_ok",
    "r_val",
    "r_status",
    "r_method",
    "r_oracle",
    "r_exact",
    "r_assumes",
    "r_witness",
    "r_search_seconds",
    "r_ladder_seconds",
    "r_total_seconds",
    "wall_until_timeout_s",
    "r_stat_elements_visited",
    "r_stat_closures",
    "r_stat_validity_checks",
    "r_stat_extensions_enumerated",
    "AUC_frac",
    "AUC_abs",
    "n_depths_evaluated",
    "n_depths_all_contradictory",
]


def _parse_float_list(s: str) -> list[float]:
    return [float(x) for x in s.split(",") if x.strip() != ""]


def _parse_int_list(s: str) -> list[int]:
    return [int(x) for x in s.split(",") if x.strip() != ""]


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode", choices=["full", "smoke", "determinism"], default="full")
    p.add_argument("--out-dir", default="results/axis_robustness")
    p.add_argument("--reps", type=int, default=200)
    p.add_argument("--instances-per-cell", type=int, default=10)
    p.add_argument("--seed-budget", type=int, default=150)
    p.add_argument("--time-limit-s", type=float, default=60.0)
    p.add_argument("--component-sizes", default="6,8,10,12")
    p.add_argument("--coverages", default="0.5,1.0")
    p.add_argument("--base-wrongness", default="0.0,0.10,0.25")
    p.add_argument("--n-tiers", default="2,3,4")
    p.add_argument(
        "--corruption-rates",
        default="0.0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50",
    )
    p.add_argument("--skip-flip", action="store_true")
    p.add_argument("--skip-tiered", action="store_true")
    p.add_argument("--resume", action="store_true")
    return p


def _fill(row: dict[str, Any], fields: list[str]) -> dict[str, Any]:
    """Project ``row`` onto exactly ``fields``, filling any gap with ``""``."""
    return {k: row.get(k, "") for k in fields}


def run_flip_cell(
    component_size: int,
    separation: int,
    coverage: float,
    base_wrongness: float,
    *,
    instances_per_cell: int,
    seed_budget: int,
    reps: int,
    time_limit_s: float,
    samples_w: ResultWriter,
    curves_w: ResultWriter,
    instances_w: ResultWriter,
    exclusions: dict[str, int],
    s0_failures: list[str],
    curves_for_monotonicity: dict[Any, dict[int, dict[str, Any]]],
) -> int:
    """Fill one flip-arm cell up to quota; returns count accepted."""
    n_accepted = 0
    for seed in range(seed_budget):
        if n_accepted >= instances_per_cell:
            break
        inst, reason = sv.build_flip_instance(component_size, separation, coverage, base_wrongness, seed)
        if inst is None:
            exclusions[reason] = exclusions.get(reason, 0) + 1
            continue
        n_accepted += 1

        s0_ok = sv.is_gac_valid_mpdag(inst["g0"], inst["x"], inst["y"], inst["z_star"])
        if not s0_ok:
            s0_failures.append(inst["instance_id"])

        n_k = inst["n_k"]
        samples: list[dict[str, Any]] = []
        for d in range(1, n_k + 1):
            for rep in range(reps):
                row = sv.sample_flip_state(inst, d, rep)
                samples.append(row)
                samples_w.write(_fill(row, SAMPLES_FIELDS))

        curve = sv.build_curve(samples)
        curves_for_monotonicity[inst["instance_id"]] = curve
        n_depths_all_contra = 0
        for d, stats in curve.items():
            crow = {
                "instance_id": inst["instance_id"],
                "arm": "flip",
                "gate": sv.GATE_NAME,
                "component_size": component_size,
                "separation": separation,
                "coverage": coverage,
                "base_wrongness": base_wrongness,
                "n_tiers": "",
                "d": d,
                **stats,
            }
            curves_w.write(_fill(crow, CURVES_FIELDS))
            if stats["S"] == "":
                n_depths_all_contra += 1

        rv = sv.compute_r_val(inst["cpdag"], inst["g0"], inst["x"], inst["y"], inst["z_star"], time_limit_s=time_limit_s)
        irow = {
            "instance_id": inst["instance_id"],
            "base_instance_id": inst["base_instance_id"],
            "arm": "flip",
            "gate": sv.GATE_NAME,
            "component_size": component_size,
            "separation": separation,
            "coverage": coverage,
            "base_wrongness": base_wrongness,
            "n_tiers": "",
            "seed": seed,
            "x": inst["x"],
            "y": inst["y"],
            "n_k": n_k,
            "shd_truth": inst["shd_truth"],
            "shd_cpdag": inst["shd_cpdag"],
            "undirected_fraction": inst["undirected_fraction"],
            "z_size": inst["z_size"],
            "z_star": ",".join(sorted(inst["z_star"])),
            "s0_ok": s0_ok,
            **rv,
            "AUC_frac": sv.auc_frac(curve, n_k),
            "AUC_abs": sv.auc_abs(curve),
            "n_depths_evaluated": len(curve),
            "n_depths_all_contradictory": n_depths_all_contra,
        }
        instances_w.write(_fill(irow, INSTANCES_FIELDS))
    return n_accepted


def run_tiered_cell(
    component_size: int,
    separation: int,
    n_tiers: int,
    *,
    instances_per_cell: int,
    seed_budget: int,
    reps: int,
    corruption_rates: list[float],
    time_limit_s: float,
    samples_w: ResultWriter,
    curves_w: ResultWriter,
    instances_w: ResultWriter,
    exclusions: dict[str, int],
    s0_failures: list[str],
    curves_for_monotonicity: dict[Any, dict[int, dict[str, Any]]],
) -> int:
    """Fill one tiered-arm cell up to quota; returns count accepted."""
    n_accepted = 0
    for seed in range(seed_budget):
        if n_accepted >= instances_per_cell:
            break
        inst, reason = sv.build_tiered_instance(component_size, separation, seed, n_tiers)
        if inst is None:
            exclusions[reason] = exclusions.get(reason, 0) + 1
            continue
        n_accepted += 1

        s0_ok = sv.is_gac_valid_mpdag(inst["g0"], inst["x"], inst["y"], inst["z_star"])
        if not s0_ok:
            s0_failures.append(inst["instance_id"])

        n_k = inst["n_k"]
        samples: list[dict[str, Any]] = []
        for rate in corruption_rates:
            for rep in range(reps):
                row = sv.sample_tiered_state(inst, rate, rep)
                samples.append(row)
                samples_w.write(_fill(row, SAMPLES_FIELDS))

        curve = sv.build_curve(samples)
        curves_for_monotonicity[inst["instance_id"]] = curve
        n_depths_all_contra = 0
        for d, stats in curve.items():
            crow = {
                "instance_id": inst["instance_id"],
                "arm": "tiered",
                "gate": sv.GATE_NAME,
                "component_size": component_size,
                "separation": separation,
                "coverage": "",
                "base_wrongness": "",
                "n_tiers": n_tiers,
                "d": d,
                **stats,
            }
            curves_w.write(_fill(crow, CURVES_FIELDS))
            if stats["S"] == "":
                n_depths_all_contra += 1

        rv = sv.compute_r_val(inst["cpdag"], inst["g0"], inst["x"], inst["y"], inst["z_star"], time_limit_s=time_limit_s)
        irow = {
            "instance_id": inst["instance_id"],
            "base_instance_id": inst["base_instance_id"],
            "arm": "tiered",
            "gate": sv.GATE_NAME,
            "component_size": component_size,
            "separation": separation,
            "coverage": "",
            "base_wrongness": "",
            "n_tiers": n_tiers,
            "seed": seed,
            "x": inst["x"],
            "y": inst["y"],
            "n_k": n_k,
            "shd_truth": inst["shd_truth"],
            "shd_cpdag": inst["shd_cpdag"],
            "undirected_fraction": inst["undirected_fraction"],
            "z_size": inst["z_size"],
            "z_star": ",".join(sorted(inst["z_star"])),
            "s0_ok": s0_ok,
            **rv,
            "AUC_frac": sv.auc_frac(curve, n_k),
            "AUC_abs": sv.auc_abs(curve),
            "n_depths_evaluated": len(curve),
            "n_depths_all_contradictory": n_depths_all_contra,
        }
        instances_w.write(_fill(irow, INSTANCES_FIELDS))
    return n_accepted


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    component_sizes = _parse_int_list(args.component_sizes)
    coverages = _parse_float_list(args.coverages)
    base_wrongness = _parse_float_list(args.base_wrongness)
    n_tiers_list = _parse_int_list(args.n_tiers)
    corruption_rates = _parse_float_list(args.corruption_rates)

    if args.mode == "determinism":
        component_sizes = [6, 8]
        coverages = [1.0]
        base_wrongness = [0.10]
        n_tiers_list = [2]
        corruption_rates = [0.0, 0.25]
        args.instances_per_cell = 1
        args.reps = 5
    elif args.mode == "smoke":
        component_sizes = component_sizes[:1]
        coverages = coverages[:1]
        base_wrongness = base_wrongness[:1]
        n_tiers_list = n_tiers_list[:1]
        args.instances_per_cell = 1

    samples_path = out_dir / "survival_samples.csv"
    curves_path = out_dir / "survival_curves.csv"
    instances_path = out_dir / "survival_instances.csv"

    exclusions: dict[str, int] = {}
    s0_failures: list[str] = []
    curves_for_monotonicity: dict[Any, dict[int, dict[str, Any]]] = {}
    cell_counts: dict[str, int] = {}

    t_start = time.perf_counter()

    with ResultWriter(samples_path, resume=args.resume) as samples_w, ResultWriter(
        curves_path, resume=args.resume
    ) as curves_w, ResultWriter(instances_path, resume=args.resume) as instances_w:
        if not args.skip_flip:
            for c in sorted(component_sizes):
                seps = sv.separations_for(c)
                for s in sorted(set(seps)):
                    for cov in sorted(coverages):
                        for b in sorted(base_wrongness):
                            n = run_flip_cell(
                                c,
                                s,
                                cov,
                                b,
                                instances_per_cell=args.instances_per_cell,
                                seed_budget=args.seed_budget,
                                reps=args.reps,
                                time_limit_s=args.time_limit_s,
                                samples_w=samples_w,
                                curves_w=curves_w,
                                instances_w=instances_w,
                                exclusions=exclusions,
                                s0_failures=s0_failures,
                                curves_for_monotonicity=curves_for_monotonicity,
                            )
                            cell_counts[f"flip:c{c}_s{s}_cov{cov}_b{b}"] = n

        if not args.skip_tiered:
            for c in sorted(component_sizes):
                seps = sv.separations_for(c)
                for s in sorted(set(seps)):
                    for nt in sorted(n_tiers_list):
                        n = run_tiered_cell(
                            c,
                            s,
                            nt,
                            instances_per_cell=args.instances_per_cell,
                            seed_budget=args.seed_budget,
                            reps=args.reps,
                            corruption_rates=corruption_rates,
                            time_limit_s=args.time_limit_s,
                            samples_w=samples_w,
                            curves_w=curves_w,
                            instances_w=instances_w,
                            exclusions=exclusions,
                            s0_failures=s0_failures,
                            curves_for_monotonicity=curves_for_monotonicity,
                        )
                        cell_counts[f"tiered:c{c}_s{s}_nt{nt}"] = n

    elapsed = time.perf_counter() - t_start
    frac_non_mono = sv.fraction_non_monotonic(curves_for_monotonicity)

    if args.mode != "determinism":
        grid = {
            "component_sizes": component_sizes,
            "coverages": coverages,
            "base_wrongness": base_wrongness,
            "n_tiers": n_tiers_list,
            "corruption_rates": corruption_rates,
            "reps": args.reps,
            "instances_per_cell": args.instances_per_cell,
            "seed_budget": args.seed_budget,
            "time_limit_s": args.time_limit_s,
            "mode": args.mode,
        }
        write_manifest(out_dir, seed=0, grid=grid, extra={"gate": sv.GATE_NAME})

        summary = {
            "elapsed_seconds": elapsed,
            "cell_counts": cell_counts,
            "exclusions": exclusions,
            "n_s0_failures": len(s0_failures),
            "s0_failures": s0_failures,
            "fraction_non_monotonic": frac_non_mono,
            "n_instances_checked_for_monotonicity": sum(
                1
                for c in curves_for_monotonicity.values()
                if len([1 for row in c.values() if row["S"] != ""]) >= 2
            ),
        }
        (out_dir / f"_run_summary_{args.mode}.json").write_text(json.dumps(summary, indent=2, default=str))
        print(json.dumps(summary, indent=2, default=str))

    print(f"elapsed_seconds={elapsed:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
