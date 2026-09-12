"""Driver for the null-stratum isolation run (flip, coverage=0.5,
base_wrongness=0.25). See :mod:`bkrobust.robustness.nullcell` for the design
rationale and every stochastic/statistical primitive; this module is
orchestration and I/O only, mirroring the shape of ``run_survival.py`` (never
imported or invoked -- this module writes only ``null_*``-prefixed files
under ``results/axis_robustness/``, never touching any session-7 file).

Usage::

    PYTHONPATH=src python3 -m bkrobust.robustness.run_nullcell --mode full
    PYTHONPATH=src python3 -m bkrobust.robustness.run_nullcell --mode determinism

``--mode determinism`` builds a small 3-instance subset and prints a sha256
of its deterministic columns (wall-clock columns excluded) to stdout, for
comparison across ``PYTHONHASHSEED`` values (K4) -- it writes no files.

No global RNG is touched anywhere in this module or in ``nullcell.py``, and
neither module imports or calls any other sweep driver in this package, the
exponential-enumeration adjustment-set search, or its gate wrapper.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import tempfile
import time
from pathlib import Path
from typing import Any

from bkrobust.core.resultsio import ResultWriter, write_manifest
from bkrobust.robustness import nullcell as nc
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
    "component_size",
    "separation",
    "coverage",
    "base_wrongness",
    "n_draws",
    "d",
    "n_eval",
    "n_contradictory",
    "contradiction_rate",
    "S",
    "S_contra_as_fail",
]

INSTANCES_FIELDS = [
    "instance_id", "base_instance_id", "arm", "gate", "component_size", "separation",
    "coverage", "base_wrongness", "seed", "x", "y", "n_k", "shd_truth", "shd_cpdag",
    "undirected_fraction", "z_size", "z_star", "s0_ok",
    "r_val", "r_status", "r_method", "r_oracle", "r_exact", "r_assumes", "r_witness",
    "r_search_seconds", "r_ladder_seconds", "r_total_seconds", "wall_until_timeout_s",
    "r_stat_elements_visited", "r_stat_closures", "r_stat_validity_checks",
    "r_stat_extensions_enumerated",
    "AUC_frac_n200", "AUC_frac_n1000", "AUC_frac_usable_n1000",
    "n_depths_evaluated_n200", "n_depths_evaluated_n1000",
    "n_usable_depths_n200", "n_usable_depths_n1000",
    "se_auc_frac_n200", "se_auc_frac_n1000",
]

TAU_FIELDS = [
    "check", "endpoint", "n_draws", "n_total", "n", "n_excluded_unreached",
    "n_excluded_nan", "tau_b", "p_value", "ci_lo_2p5", "ci_hi_97p5", "n_boot_nan",
]


def _fill(row: dict[str, Any], fields: list[str]) -> dict[str, Any]:
    return {k: row.get(k, "") for k in fields}


# --- determinism mode (K4) ------------------------------------------------------

DETERMINISM_DEPTH_REPS = 50


def _canonical_determinism_payload() -> str:
    """3 instances from the first cell, sampled at a small N, rendered as a
    canonical string over deterministic columns only (wall-clock columns --
    r_search_seconds, r_ladder_seconds, r_total_seconds, wall_until_timeout_s,
    r_stat_* which can differ by machine load even at fixed inputs on some
    platforms -- are excluded on purpose, noted at the call site).
    """
    c, s = nc.cells()[0]
    exclusions: dict[str, int] = {}
    instances = nc.build_cell_instances(c, s, exclusions=exclusions)[:3]
    lines = []
    for inst in instances:
        rows = nc.sample_instance(inst, n_draws=DETERMINISM_DEPTH_REPS)
        rv = sv.compute_r_val(inst["cpdag"], inst["g0"], inst["x"], inst["y"], inst["z_star"])
        lines.append(
            f"instance_id={inst['instance_id']} n_k={inst['n_k']} "
            f"z_star={','.join(sorted(inst['z_star']))} r_val={rv['r_val']} "
            f"r_status={rv['r_status']} shd_truth={inst['shd_truth']}"
        )
        for row in sorted(rows, key=lambda r: (r["d"], r["rep"])):
            lines.append(
                f"  d={row['d']} rep={row['rep']} seed={row['seed']} status={row['status']} "
                f"survived={row['survived']} symdiff={row['symdiff_proxy_not_distance']}"
            )
    return "\n".join(lines)


def run_determinism_mode() -> int:
    payload = _canonical_determinism_payload()
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    print(f"n_lines_hashed={len(payload.splitlines())}")
    print(f"sha256={digest}")
    return 0


# --- full mode -------------------------------------------------------------------


def run_full_mode(out_dir: Path) -> int:
    samples_path = out_dir / "null_samples.csv"
    curves_path = out_dir / "null_curves.csv"
    instances_path = out_dir / "null_instances.csv"
    tau_path = out_dir / "null_tau_comparison.csv"

    t_start = time.perf_counter()

    print("[1/5] building instance population (fixed grid: component_sizes="
          f"{nc.COMPONENT_SIZES}, coverage={nc.COVERAGE}, base_wrongness={nc.BASE_WRONGNESS}, "
          f"instances_per_cell={nc.INSTANCES_PER_CELL}, seed_budget={nc.SEED_BUDGET}) ...")
    instances, build_exclusions = nc.build_all_instances()
    print(f"    accepted {len(instances)} instances across {len(nc.cells())} cells")
    print(f"    exclusions: {build_exclusions}")

    print("[2/5] K2 check: S(0) = 1.0 for every instance ...")
    s0_failures = nc.check_s0(instances)
    if s0_failures:
        print(f"    STOP: {len(s0_failures)} instance(s) fail S(0)=1.0: {s0_failures}")
        return 1
    print(f"    OK: all {len(instances)} instances pass (Z* valid at their own G0)")

    print("[3/5] K1 timing: one instance at N=1000 ...")
    t0 = time.perf_counter()
    probe_rows = nc.sample_instance(instances[0], n_draws=nc.N_DRAWS_FULL)
    t1 = time.perf_counter()
    probe_seconds = t1 - t0
    total_nk = sum(inst["n_k"] for inst in instances)
    projected_seconds = probe_seconds / instances[0]["n_k"] * total_nk
    print(
        f"    instance {instances[0]['instance_id']} (n_k={instances[0]['n_k']}): "
        f"{probe_seconds:.3f}s for {len(probe_rows)} draws "
        f"({probe_seconds / len(probe_rows) * 1000:.4f} ms/draw)"
    )
    print(
        f"    projection: sum(n_k)={total_nk} across {len(instances)} instances -> "
        f"~{projected_seconds:.1f}s for all sampling (excludes instance-build and r_val time, "
        "both separately measured as negligible below)"
    )

    print("[4/5] sampling all instances at N=1000, writing incrementally ...")
    tau_rows_by_check: dict[str, dict[str, list[float]]] = {
        "check1_n200_reproduction": {"r_val": [], "endpoint": []},
        "check2_n1000_full": {"r_val": [], "endpoint": []},
        "check3_n1000_usable": {"r_val": [], "endpoint": []},
    }
    contra_by_depth_n200: dict[int, list[int]] = {}
    contra_by_depth_n1000: dict[int, list[int]] = {}

    with ResultWriter(samples_path) as samples_w, ResultWriter(curves_path) as curves_w, ResultWriter(
        instances_path
    ) as instances_w:
        for i, inst in enumerate(instances):
            rows = probe_rows if i == 0 else nc.sample_instance(inst, n_draws=nc.N_DRAWS_FULL)
            for row in rows:
                samples_w.write(_fill(row, SAMPLES_FIELDS))

            curve_1000 = nc.curve_from_rows(rows)
            curve_200 = nc.curve_from_rows(rows, rep_limit=nc.N_DRAWS_CHECK)

            for n_draws, curve, agg in (
                (200, curve_200, contra_by_depth_n200),
                (1000, curve_1000, contra_by_depth_n1000),
            ):
                rates = nc.contradiction_rate_by_depth(curve)
                for d, stats in curve.items():
                    crow = {
                        "instance_id": inst["instance_id"],
                        "component_size": inst["component_size"],
                        "separation": inst["separation"],
                        "coverage": inst["coverage"],
                        "base_wrongness": inst["base_wrongness"],
                        "n_draws": n_draws,
                        "d": d,
                        "contradiction_rate": rates[d],
                        **stats,
                    }
                    curves_w.write(_fill(crow, CURVES_FIELDS))
                    bucket = agg.setdefault(d, [0, 0])  # [n_eval, n_contradictory]
                    bucket[0] += stats["n_eval"]
                    bucket[1] += stats["n_contradictory"]

            n_k = inst["n_k"]
            auc_200 = sv.auc_frac(curve_200, n_k)
            auc_1000 = sv.auc_frac(curve_1000, n_k)
            auc_usable_1000 = nc.auc_frac_usable(curve_1000, n_k)
            se_200 = nc.auc_frac_se(curve_200, n_k)
            se_1000 = nc.auc_frac_se(curve_1000, n_k)
            n_usable_200 = nc.n_usable_depths(curve_200)
            n_usable_1000 = nc.n_usable_depths(curve_1000)

            rv = sv.compute_r_val(inst["cpdag"], inst["g0"], inst["x"], inst["y"], inst["z_star"])

            irow = {
                "instance_id": inst["instance_id"],
                "base_instance_id": inst["base_instance_id"],
                "arm": "flip",
                "gate": sv.GATE_NAME,
                "component_size": inst["component_size"],
                "separation": inst["separation"],
                "coverage": inst["coverage"],
                "base_wrongness": inst["base_wrongness"],
                "seed": inst["seed"],
                "x": inst["x"],
                "y": inst["y"],
                "n_k": n_k,
                "shd_truth": inst["shd_truth"],
                "shd_cpdag": inst["shd_cpdag"],
                "undirected_fraction": inst["undirected_fraction"],
                "z_size": inst["z_size"],
                "z_star": ",".join(sorted(inst["z_star"])),
                "s0_ok": True,
                **rv,
                "AUC_frac_n200": auc_200,
                "AUC_frac_n1000": auc_1000,
                "AUC_frac_usable_n1000": auc_usable_1000,
                "n_depths_evaluated_n200": len(curve_200),
                "n_depths_evaluated_n1000": len(curve_1000),
                "n_usable_depths_n200": n_usable_200,
                "n_usable_depths_n1000": n_usable_1000,
                "se_auc_frac_n200": se_200,
                "se_auc_frac_n1000": se_1000,
            }
            instances_w.write(_fill(irow, INSTANCES_FIELDS))

            r_val = rv["r_val"] if rv["r_status"] == "ok" else float("nan")
            tau_rows_by_check["check1_n200_reproduction"]["r_val"].append(r_val)
            tau_rows_by_check["check1_n200_reproduction"]["endpoint"].append(
                float(auc_200) if auc_200 != "" else float("nan")
            )
            tau_rows_by_check["check2_n1000_full"]["r_val"].append(r_val)
            tau_rows_by_check["check2_n1000_full"]["endpoint"].append(
                float(auc_1000) if auc_1000 != "" else float("nan")
            )
            tau_rows_by_check["check3_n1000_usable"]["r_val"].append(r_val)
            tau_rows_by_check["check3_n1000_usable"]["endpoint"].append(
                float(auc_usable_1000) if auc_usable_1000 != "" else float("nan")
            )

            if (i + 1) % 25 == 0 or (i + 1) == len(instances):
                print(f"    {i + 1}/{len(instances)} instances done "
                      f"({time.perf_counter() - t_start:.1f}s elapsed)")

    sampling_elapsed = time.perf_counter() - t_start
    print(f"    sampling phase complete in {sampling_elapsed:.1f}s")

    print("[5/5] computing the three tau_b comparisons and writing outputs ...")
    r_assumes = "Conjecture 2 (hence Anti-Exchange Case B, verified not proved)"
    tau_out_rows = []
    endpoint_names = {
        "check1_n200_reproduction": "AUC_frac (first 200 of 1000 draws)",
        "check2_n1000_full": "AUC_frac (full 1000 draws)",
        "check3_n1000_usable": "AUC_frac_usable, n_eval>=30 (1000 draws)",
    }
    n_draws_names = {"check1_n200_reproduction": 200, "check2_n1000_full": 1000, "check3_n1000_usable": 1000}
    for check, data in tau_rows_by_check.items():
        result = nc.tau_b_stratum(data["r_val"], data["endpoint"])
        tau_out_rows.append(
            _fill(
                {
                    "check": check,
                    "endpoint": endpoint_names[check],
                    "n_draws": n_draws_names[check],
                    **result,
                },
                TAU_FIELDS,
            )
        )
        print(f"    {check}: tau_b={result['tau_b']:.4f} "
              f"CI=[{result['ci_lo_2p5']:.4f}, {result['ci_hi_97p5']:.4f}] n={result['n']}")

    with ResultWriter(tau_path) as tau_w:
        for row in tau_out_rows:
            tau_w.write(row)

    # gzip the samples file if it exceeds 20MB, per the output spec.
    samples_size = samples_path.stat().st_size
    final_samples_path = samples_path
    if samples_size > 20 * 1024 * 1024:
        gz_path = samples_path.with_suffix(".csv.gz")
        with samples_path.open("rb") as f_in, gzip.open(gz_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        samples_path.unlink()
        final_samples_path = gz_path
        print(f"    compressed {samples_size / 1e6:.1f}MB -> {gz_path.stat().st_size / 1e6:.1f}MB ({gz_path.name})")

    contradiction_summary = {
        "n200": {d: (b[1] / (b[0] + b[1]) if (b[0] + b[1]) else None) for d, b in sorted(contra_by_depth_n200.items())},
        "n1000": {d: (b[1] / (b[0] + b[1]) if (b[0] + b[1]) else None) for d, b in sorted(contra_by_depth_n1000.items())},
    }

    grid = {
        "component_sizes": list(nc.COMPONENT_SIZES),
        "coverage": nc.COVERAGE,
        "base_wrongness": nc.BASE_WRONGNESS,
        "instances_per_cell": nc.INSTANCES_PER_CELL,
        "seed_budget": nc.SEED_BUDGET,
        "n_draws_check": nc.N_DRAWS_CHECK,
        "n_draws_full": nc.N_DRAWS_FULL,
        "usable_min_n_eval": nc.USABLE_MIN_N_EVAL,
        "n_boot": nc.N_BOOT,
        "boot_seed": nc.BOOT_SEED,
        "mode": "full",
    }
    # write_manifest() unconditionally targets "<dir>/manifest.json" -- and
    # results/axis_robustness/manifest.json is session 7's committed,
    # read-only manifest. Writing there (even transiently, even followed by
    # a rename) would clobber it in the working tree. So write_manifest is
    # pointed at a throwaway temp directory instead, and only *that* file's
    # bytes are copied into out_dir/null_manifest.json -- a brand-new path,
    # never an existing one.
    with tempfile.TemporaryDirectory() as tmp:
        manifest_path = write_manifest(
            tmp,
            seed=0,
            grid=grid,
            extra={
                "gate": sv.GATE_NAME,
                "output_prefix": "null_",
                "n_instances": len(instances),
                "build_exclusions": build_exclusions,
                "s0_failures": s0_failures,
                "contradiction_rate_by_depth": contradiction_summary,
                "elapsed_seconds": time.perf_counter() - t_start,
                "samples_file": final_samples_path.name,
            },
        )
        null_manifest_path = out_dir / "null_manifest.json"
        if null_manifest_path.exists():
            raise FileExistsError(f"{null_manifest_path} exists; refusing to overwrite")
        shutil.copyfile(manifest_path, null_manifest_path)

    print(f"elapsed_seconds={time.perf_counter() - t_start:.2f}")
    return 0


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode", choices=["full", "determinism"], default="full")
    p.add_argument("--out-dir", default="results/axis_robustness")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if args.mode == "determinism":
        return run_determinism_mode()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    return run_full_mode(out_dir)


if __name__ == "__main__":
    raise SystemExit(main())
