"""Driver for the v2 "unconfounded" tiered-arm sweep (separation not targeted).

Companion to ``results/axis_robustness_p6/run_survival_p6.py``, restricted to
the tiered arm and swapping instance/query generation for
:mod:`bkrobust.robustness.unconfounded_generator` -- see that module's
docstring for exactly what "not targeted" means here. Everything downstream
of instance construction is reused, unmodified, from
:mod:`bkrobust.robustness.survival` and
:mod:`bkrobust.robustness.survival_p6`: the corruption sampler
(``sample_tiered_state_p6``), the curve builder (``build_curve``), the AUC
endpoint (``auc_frac``/``auc_abs``), ``compute_r_val``, and the new baseline
predictors (``instance_predictors_p6``, though ``separation`` /
``separation_status`` are overridden with this run's own measured values
rather than the p6 module's spec-derived ones -- see ``run_cell`` below).

Grid: ``component_size in {6, 8, 10, 12}`` x ``n_tiers in {2, 3, 4}`` -- 12
cells, no ``separation`` dimension (nothing to target). Same ``reps=1000``,
same ``corruption_rate`` grid (11 points, 0.0..0.5 step 0.05), same
``fast_gate`` admission, same predictor columns as ``run_survival_p6.py``'s
tiered arm.

Sharded by ``component_size`` (4 shards), each an independent OS process
writing its own CSVs/manifest/log under ``--out-dir``, resumable via
``--resume`` (``ResultWriter`` is append-only; killing a shard loses at most
the row in flight, and ``--resume`` skips seeds already written for cells
that reached quota -- see ``run_cell``). ``--merge`` concatenates shards
deterministically, exactly as ``run_survival_p6.py --merge`` does.

Usage::

    PYTHONPATH=src .venv/bin/python -m bkrobust.robustness.run_unconfounded_v2 \\
        --out-dir results/axis_robustness_unconfounded_v2/shards/c06 \\
        --shard-name c06 --component-sizes 6 --reps 1000 --instances-per-cell 40

    PYTHONPATH=src .venv/bin/python -m bkrobust.robustness.run_unconfounded_v2 \\
        --merge --out-dir results/axis_robustness_unconfounded_v2 \\
        --shards results/axis_robustness_unconfounded_v2/shards/c06,...
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path
from typing import Any

from bkrobust.core.resultsio import ResultWriter, write_manifest
from bkrobust.robustness import survival as sv
from bkrobust.robustness import survival_p6 as sv6
from bkrobust.robustness.unconfounded_generator import build_unconfounded_tiered_instance

SAMPLES_FIELDS = sv6.SAMPLES_FIELDS_P6
CURVES_FIELDS = sv6.CURVES_FIELDS_P6
INSTANCES_FIELDS = [
    *sv6.INSTANCES_FIELDS_P6,
    "anchor",
    "n_admissible_queries",
]


def _parse_int_list(s: str) -> list[int]:
    return [int(x) for x in s.split(",") if x.strip() != ""]


def _parse_float_list(s: str) -> list[float]:
    return [float(x) for x in s.split(",") if x.strip() != ""]


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode", choices=["full", "smoke", "pilot"], default="full")
    p.add_argument("--out-dir", default="results/axis_robustness_unconfounded_v2")
    p.add_argument("--shard-name", default=None)
    p.add_argument("--log-dir", default=None, help="Defaults to <out-dir's parent's parent>/_logs.")
    p.add_argument("--reps", type=int, default=1000)
    p.add_argument("--instances-per-cell", type=int, default=40)
    p.add_argument("--seed-budget", type=int, default=200)
    p.add_argument("--time-limit-s", type=float, default=60.0)
    p.add_argument("--component-sizes", default="6,8,10,12")
    p.add_argument("--n-tiers", default="2,3,4")
    p.add_argument(
        "--corruption-rates",
        default="0.0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50",
    )
    p.add_argument("--triangle-prob", type=float, default=0.3)
    p.add_argument("--resume", action="store_true")
    p.add_argument("--merge", action="store_true")
    p.add_argument("--shards", default="")
    return p


def _fill(row: dict[str, Any], fields: list[str]) -> dict[str, Any]:
    return {k: row.get(k, "") for k in fields}


class _Logger:
    def __init__(self, path: Path, shard_name: str) -> None:
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.shard_name = shard_name
        self._fh = self.path.open("a")

    def log(self, msg: str) -> None:
        line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] [{self.shard_name}] {msg}"
        print(line, flush=True)
        self._fh.write(line + "\n")
        self._fh.flush()

    def close(self) -> None:
        self._fh.close()


def run_cell(
    component_size: int,
    n_tiers: int,
    *,
    instances_per_cell: int,
    seed_budget: int,
    reps: int,
    corruption_rates: list[float],
    time_limit_s: float,
    triangle_prob: float,
    samples_w: ResultWriter,
    curves_w: ResultWriter,
    instances_w: ResultWriter,
    exclusions: dict[str, int],
    s0_failures: list[str],
    curves_for_monotonicity: dict[Any, dict[int, dict[str, Any]]],
    denom_totals: dict[str, int],
    separation_status_totals: dict[str, int],
    already_done_instance_ids: set[str],
    logger: _Logger,
) -> int:
    """Fill one (component_size, n_tiers) cell up to quota; returns count accepted.

    ``already_done_instance_ids`` supports ``--resume``: instances whose id
    (deterministic function of ``component_size``, ``seed``, the drawn ``x``,
    and ``n_tiers`` -- see ``unconfounded_generator.build_unconfounded_tiered_instance``)
    already appears in a previously-written ``survival_instances.csv`` are
    skipped without recomputation, but still count toward quota.
    """
    n_accepted = 0
    n_samples_cell = 0
    for seed in range(seed_budget):
        if n_accepted >= instances_per_cell:
            break
        inst, reason = build_unconfounded_tiered_instance(
            component_size, seed, n_tiers, triangle_prob=triangle_prob
        )
        if inst is None:
            exclusions[reason] = exclusions.get(reason, 0) + 1
            continue
        if inst["instance_id"] in already_done_instance_ids:
            n_accepted += 1
            continue
        n_accepted += 1

        separation_status_totals[inst["separation_status"]] = (
            separation_status_totals.get(inst["separation_status"], 0) + 1
        )

        s0_ok = sv.is_gac_valid_mpdag(inst["g0"], inst["x"], inst["y"], inst["z_star"])
        if not s0_ok:
            s0_failures.append(inst["instance_id"])

        n_k = inst["n_k"]
        samples: list[dict[str, Any]] = []
        for rate in corruption_rates:
            for rep in range(reps):
                row = sv6.sample_tiered_state_p6(inst, rate, rep)
                samples.append(row)
                samples_w.write(_fill(row, SAMPLES_FIELDS))
                n_samples_cell += 1
                denom_totals["n_claims_attempted"] += row["n_claims_attempted"]
                denom_totals["n_claims_reversed"] += row["n_claims_reversed"]
                denom_totals["corruption_rejected"] += row["corruption_rejected"]
                denom_totals["corruption_accepted"] += row["corruption_accepted"]
                denom_totals["closure_inert"] += row["closure_inert"]
                denom_totals["n_rows"] += 1

        curve = sv.build_curve(samples)
        curves_for_monotonicity[inst["instance_id"]] = curve
        n_depths_all_contra = 0
        for d, stats in curve.items():
            crow = {
                "instance_id": inst["instance_id"],
                "arm": "tiered_unconfounded",
                "gate": sv.GATE_NAME,
                "component_size": component_size,
                "separation": inst["separation"] if inst["separation"] is not None else "",
                "coverage": "",
                "base_wrongness": "",
                "n_tiers": n_tiers,
                "d": d,
                **stats,
            }
            curves_w.write(_fill(crow, CURVES_FIELDS))
            if stats["S"] == "":
                n_depths_all_contra += 1

        rv = sv.compute_r_val(
            inst["cpdag"], inst["g0"], inst["x"], inst["y"], inst["z_star"], time_limit_s=time_limit_s
        )
        # instance_predictors_p6 assumes instance["separation"] is the
        # spec-targeted value (always "measured" there, by construction of
        # component_generator's admission gate -- see survival_p6.py). This
        # run has no spec-targeted separation, so separation /
        # separation_status are taken from this module's own measured value
        # instead of instance_predictors_p6's copy-through.
        preds = sv6.instance_predictors_p6(inst)
        preds["separation"] = inst["separation"] if inst["separation"] is not None else ""
        preds["separation_status"] = inst["separation_status"]
        irow = {
            "instance_id": inst["instance_id"],
            "base_instance_id": inst["base_instance_id"],
            "arm": "tiered_unconfounded",
            "gate": sv.GATE_NAME,
            "component_size": component_size,
            "coverage": "",
            "base_wrongness": "",
            "n_tiers": n_tiers,
            "seed": seed,
            "x": inst["x"],
            "y": inst["y"],
            "anchor": inst["anchor"],
            "n_admissible_queries": inst["n_admissible_queries"],
            "n_k": n_k,
            "shd_truth": inst["shd_truth"],
            "shd_cpdag": inst["shd_cpdag"],
            "undirected_fraction": inst["undirected_fraction"],
            "z_size": inst["z_size"],
            "z_star": ",".join(sorted(inst["z_star"])),
            "s0_ok": s0_ok,
            **preds,
            **rv,
            "AUC_frac": sv.auc_frac(curve, n_k),
            "AUC_abs": sv.auc_abs(curve),
            "n_depths_evaluated": len(curve),
            "n_depths_all_contradictory": n_depths_all_contra,
        }
        instances_w.write(_fill(irow, INSTANCES_FIELDS))
    logger.log(
        f"tiered_unconfounded cell done: c={component_size} nt={n_tiers} "
        f"accepted={n_accepted}/{instances_per_cell} samples={n_samples_cell}"
    )
    return n_accepted


def run(args: argparse.Namespace) -> int:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    shard_name = args.shard_name or out_dir.name
    log_dir = Path(args.log_dir) if args.log_dir else out_dir.parent.parent / "_logs"
    logger = _Logger(log_dir / f"{shard_name}.log", shard_name)

    component_sizes = _parse_int_list(args.component_sizes)
    n_tiers_list = _parse_int_list(args.n_tiers)
    corruption_rates = _parse_float_list(args.corruption_rates)

    if args.mode == "smoke":
        component_sizes = component_sizes[:1]
        n_tiers_list = n_tiers_list[:1]
        args.instances_per_cell = 1
    elif args.mode == "pilot":
        args.instances_per_cell = min(args.instances_per_cell, 3)

    logger.log(
        f"starting: mode={args.mode} reps={args.reps} instances_per_cell={args.instances_per_cell} "
        f"seed_budget={args.seed_budget} component_sizes={component_sizes} n_tiers={n_tiers_list} "
        f"triangle_prob={args.triangle_prob} resume={args.resume}"
    )

    samples_path = out_dir / "survival_samples.csv"
    curves_path = out_dir / "survival_curves.csv"
    instances_path = out_dir / "survival_instances.csv"

    already_done_instance_ids: set[str] = set()
    if args.resume and instances_path.exists():
        with instances_path.open(newline="") as fh:
            for row in csv.DictReader(fh):
                iid = row.get("instance_id")
                if iid:
                    already_done_instance_ids.add(iid)
        logger.log(f"resume: found {len(already_done_instance_ids)} already-written instances")

    exclusions: dict[str, int] = {}
    s0_failures: list[str] = []
    curves_for_monotonicity: dict[Any, dict[int, dict[str, Any]]] = {}
    cell_counts: dict[str, int] = {}
    separation_status_totals: dict[str, int] = {}
    denom_totals: dict[str, int] = {
        "n_claims_attempted": 0,
        "n_claims_reversed": 0,
        "corruption_rejected": 0,
        "corruption_accepted": 0,
        "closure_inert": 0,
        "n_rows": 0,
    }

    t_start = time.perf_counter()
    n_cells_total = len(set(component_sizes)) * len(set(n_tiers_list))
    n_cells_done = 0

    with ResultWriter(samples_path, resume=args.resume) as samples_w, ResultWriter(
        curves_path, resume=args.resume
    ) as curves_w, ResultWriter(instances_path, resume=args.resume) as instances_w:
        for c in sorted(set(component_sizes)):
            for nt in sorted(set(n_tiers_list)):
                t_cell = time.perf_counter()
                n = run_cell(
                    c,
                    nt,
                    instances_per_cell=args.instances_per_cell,
                    seed_budget=args.seed_budget,
                    reps=args.reps,
                    corruption_rates=corruption_rates,
                    time_limit_s=args.time_limit_s,
                    triangle_prob=args.triangle_prob,
                    samples_w=samples_w,
                    curves_w=curves_w,
                    instances_w=instances_w,
                    exclusions=exclusions,
                    s0_failures=s0_failures,
                    curves_for_monotonicity=curves_for_monotonicity,
                    denom_totals=denom_totals,
                    separation_status_totals=separation_status_totals,
                    already_done_instance_ids=already_done_instance_ids,
                    logger=logger,
                )
                cell_counts[f"tiered_unconfounded:c{c}_nt{nt}"] = n
                n_cells_done += 1
                elapsed = time.perf_counter() - t_start
                logger.log(
                    f"progress: {n_cells_done}/{n_cells_total} cells, "
                    f"cell_wall={time.perf_counter() - t_cell:.1f}s, "
                    f"total_elapsed={elapsed:.1f}s, total_rows={denom_totals['n_rows']}"
                )

    elapsed = time.perf_counter() - t_start
    frac_non_mono = sv.fraction_non_monotonic(curves_for_monotonicity)

    grid = {
        "component_sizes": component_sizes,
        "n_tiers": n_tiers_list,
        "corruption_rates": corruption_rates,
        "reps": args.reps,
        "instances_per_cell": args.instances_per_cell,
        "seed_budget": args.seed_budget,
        "time_limit_s": args.time_limit_s,
        "triangle_prob": args.triangle_prob,
        "mode": args.mode,
        "shard_name": shard_name,
    }
    write_manifest(
        out_dir,
        seed=0,
        grid=grid,
        extra={
            "gate": sv.GATE_NAME,
            "claim_radius_max_depth": sv6.CLAIM_RADIUS_MAX_DEPTH,
            "design": (
                "separation is NOT a targeted spec parameter here; the query "
                "(treatment) is drawn uniformly among admissible candidates "
                "after the graph is built. See unconfounded_generator.py "
                "module docstring for the exact rule."
            ),
        },
    )

    summary = {
        "shard_name": shard_name,
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
        "corruption_denominator_totals": denom_totals,
        "separation_status_totals": separation_status_totals,
    }
    (out_dir / f"_run_summary_{shard_name}.json").write_text(json.dumps(summary, indent=2, default=str))
    logger.log(f"finished: elapsed_seconds={elapsed:.2f} cells={n_cells_done}/{n_cells_total}")
    logger.log(json.dumps(summary, default=str))
    logger.close()
    print(f"elapsed_seconds={elapsed:.2f}")
    return 0


# --- merge mode (identical scheme to run_survival_p6.py --merge) -------------


def _read_csv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as fh:
        r = csv.DictReader(fh)
        fields = r.fieldnames or []
        rows = list(r)
    return fields, rows


def _sort_key_samples(row: dict[str, str]):
    d = row.get("d", "")
    corr = row.get("corruption_rate", "")
    rep = row.get("rep", "")
    return (
        row.get("instance_id", ""),
        float(d) if d != "" else -1.0,
        float(corr) if corr != "" else -1.0,
        int(rep) if rep != "" else -1,
    )


def _sort_key_curves(row: dict[str, str]):
    d = row.get("d", "")
    return (row.get("instance_id", ""), float(d) if d != "" else -1.0)


def _sort_key_instances(row: dict[str, str]):
    return (row.get("instance_id", ""),)


def merge(args: argparse.Namespace) -> int:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    shard_dirs = [Path(s) for s in args.shards.split(",") if s.strip()]
    if not shard_dirs:
        print("no --shards given", file=sys.stderr)
        return 1

    merged: dict[str, list[dict[str, str]]] = {"samples": [], "curves": [], "instances": []}
    fields: dict[str, list[str]] = {}
    shard_manifests = []
    shard_summaries = []
    for sd in shard_dirs:
        for key, fname in (
            ("samples", "survival_samples.csv"),
            ("curves", "survival_curves.csv"),
            ("instances", "survival_instances.csv"),
        ):
            p = sd / fname
            if not p.exists():
                continue
            f, rows = _read_csv_rows(p)
            if key in fields and fields[key] != f:
                raise ValueError(f"schema mismatch in {p}: {f} != {fields[key]}")
            fields[key] = f
            merged[key].extend(rows)
        manifest_p = sd / "manifest.json"
        if manifest_p.exists():
            shard_manifests.append(json.loads(manifest_p.read_text()))
        for summary_p in sorted(sd.glob("_run_summary_*.json")):
            shard_summaries.append(json.loads(summary_p.read_text()))

    merged["samples"].sort(key=_sort_key_samples)
    merged["curves"].sort(key=_sort_key_curves)
    merged["instances"].sort(key=_sort_key_instances)

    for key, fname, gz in (
        ("samples", "survival_samples.csv.gz", True),
        ("curves", "survival_curves.csv", False),
        ("instances", "survival_instances.csv", False),
    ):
        if key not in fields:
            continue
        out_path = out_dir / fname
        if gz:
            import gzip

            with gzip.open(out_path, "wt", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=fields[key])
                w.writeheader()
                for row in merged[key]:
                    w.writerow(row)
        else:
            with out_path.open("w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=fields[key])
                w.writeheader()
                for row in merged[key]:
                    w.writerow(row)
        print(f"wrote {out_path}: {len(merged[key])} rows")

    total_denom = {
        "n_claims_attempted": 0,
        "n_claims_reversed": 0,
        "corruption_rejected": 0,
        "corruption_accepted": 0,
        "closure_inert": 0,
        "n_rows": 0,
    }
    total_elapsed = 0.0
    all_s0_failures: list[str] = []
    separation_status_totals: dict[str, int] = {}
    for s in shard_summaries:
        total_elapsed += s.get("elapsed_seconds", 0.0)
        all_s0_failures.extend(s.get("s0_failures", []))
        dt = s.get("corruption_denominator_totals", {})
        for k in total_denom:
            total_denom[k] += dt.get(k, 0)
        for k, v in s.get("separation_status_totals", {}).items():
            separation_status_totals[k] = separation_status_totals.get(k, 0) + v

    manifest_payload = {
        "merged_from_shards": [str(s) for s in shard_dirs],
        "shard_manifests": shard_manifests,
        "shard_summaries_elapsed_seconds": {
            s.get("shard_name", "?"): s.get("elapsed_seconds") for s in shard_summaries
        },
        "sum_of_shard_elapsed_seconds_cpu": total_elapsed,
        "total_s0_failures": len(all_s0_failures),
        "s0_failures": all_s0_failures,
        "corruption_denominator_totals": total_denom,
        "separation_status_totals": separation_status_totals,
        "n_samples_rows": len(merged["samples"]),
        "n_curves_rows": len(merged["curves"]),
        "n_instances_rows": len(merged["instances"]),
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest_payload, indent=2, default=str))
    print(json.dumps(manifest_payload, indent=2, default=str)[:4000])
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if args.merge:
        return merge(args)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
