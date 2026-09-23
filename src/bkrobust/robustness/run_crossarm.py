"""Driver for the matched cross-arm corruption sweep (Axis Robustness).

Builds matched instances (shared ``cpdag``, ``K_ref``, ``G0``, ``Z*``), then
corrupts each one both ways -- ``tiered`` and ``flip`` -- scored on the
shared, unified intensity axis. See :mod:`bkrobust.robustness.crossarm` for
the corruption/scoring logic; this module is orchestration and I/O only, and
never touches a global RNG.

Usage::

    PYTHONPATH=src python3 -m bkrobust.robustness.run_crossarm sweep \\
        --out-dir results/axis_robustness --component-sizes 6 \\
        --separations 1,3,5 --n-tiers 2,3,4 --instances-per-cell 5 \\
        --reps 200 --resume

    PYTHONPATH=src python3 -m bkrobust.robustness.run_crossarm finalize \\
        --out-dir results/axis_robustness

    PYTHONPATH=src python3 -m bkrobust.robustness.run_crossarm time-one \\
        --component-size 6 --n-tiers 2

Run in chunks (e.g. one Bash call per component size) with ``sweep
--resume`` across calls to stay under any single-command time budget --
every writer is :class:`~bkrobust.core.resultsio.ResultWriter`, which
flushes per row, so a chunk that stops partway loses at most the row it was
computing. ``finalize`` is idempotent-safe to call once, after every sweep
chunk is done: it reads back ``xarm_bins.csv`` (a file this run created, not
a pre-existing one) to compute the paired cross-arm test and writes the
remaining new (``xarm_``-prefixed) output files.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from bkrobust.core.resultsio import ResultWriter, write_manifest
from bkrobust.gac import is_gac_valid_mpdag
from bkrobust.robustness import crossarm as cx

SAMPLES_FIELDS = [
    "instance_id",
    "arm",
    "gate",
    "component_size",
    "separation",
    "n_tiers",
    "grid_point",
    "flip_rate",
    "corruption_rate",
    "d",
    "rep",
    "seed",
    "status",
    "survived",
    "symdiff",
    "n_dir_g0",
    "intensity",
    "intensity_bin",
]

BINS_FIELDS = [
    "instance_id",
    "arm",
    "gate",
    "component_size",
    "separation",
    "n_tiers",
    "bin",
    "n_samples",
    "n_contradictory",
    "n_noncontra",
    "contradiction_rate",
    "S",
    "S_contra_as_fail",
]

INSTANCES_FIELDS = [
    "instance_id",
    "base_instance_id",
    "gate",
    "component_size",
    "separation",
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
    "n_dir_g0",
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
    "n_grid_points_tiered",
    "n_grid_points_unresolved_tiered",
    "n_grid_points_flip",
    "n_grid_points_unresolved_flip",
    "contradiction_rate_overall_tiered",
    "contradiction_rate_overall_flip",
    "AUC_intensity_usable_tiered",
    "AUC_intensity_usable_flip",
]

PAIRED_TEST_FIELDS = [
    "bin",
    "n_instances_matched",
    "mean_diff_S_tiered_minus_flip",
    "median_diff",
    "n_pos",
    "n_neg",
    "n_zero",
    "sign_test_p",
    "bootstrap_mean",
    "bootstrap_ci_lo",
    "bootstrap_ci_hi",
    "n_bootstrap",
]


def _fill(row: dict[str, Any], fields: list[str]) -> dict[str, Any]:
    return {k: row.get(k, "") for k in fields}


def _parse_int_list(s: str) -> list[int]:
    return [int(x) for x in s.split(",") if x.strip() != ""]


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    sweep = sub.add_parser("sweep", help="Build matched instances and sample both arms.")
    sweep.add_argument("--out-dir", default="results/axis_robustness")
    sweep.add_argument("--reps", type=int, default=cx.N_REPS_DEFAULT)
    sweep.add_argument("--instances-per-cell", type=int, default=20)
    sweep.add_argument("--seed-budget", type=int, default=200)
    sweep.add_argument("--time-limit-s", type=float, default=60.0)
    sweep.add_argument("--component-sizes", default="6,8,10")
    sweep.add_argument("--separations", default="", help="Comma list; empty = all achievable per size.")
    sweep.add_argument("--n-tiers", default="2,3,4")
    sweep.add_argument("--resume", action="store_true")
    sweep.add_argument("--skip-r-val", action="store_true", help="Skip breakdown_radius (time-saving).")

    fin = sub.add_parser("finalize", help="Compute the paired test and write manifest/notes.")
    fin.add_argument("--out-dir", default="results/axis_robustness")

    t1 = sub.add_parser("time-one", help="Time exactly one matched instance end-to-end.")
    t1.add_argument("--component-size", type=int, default=6)
    t1.add_argument("--separation", type=int, default=None)
    t1.add_argument("--n-tiers", type=int, default=2)
    t1.add_argument("--seed-budget", type=int, default=10)
    t1.add_argument("--reps", type=int, default=cx.N_REPS_DEFAULT)
    t1.add_argument("--time-limit-s", type=float, default=60.0)

    return p


# --- sweep ---------------------------------------------------------------------


def run_cell(
    component_size: int,
    separation: int,
    n_tiers: int,
    *,
    instances_per_cell: int,
    seed_budget: int,
    reps: int,
    time_limit_s: float,
    skip_r_val: bool,
    samples_w: ResultWriter,
    bins_w: ResultWriter,
    instances_w: ResultWriter,
    exclusions: dict[str, int],
    s0_failures: list[str],
) -> int:
    """Fill one (component_size, separation, n_tiers) cell up to quota."""
    n_accepted = 0
    for seed in range(seed_budget):
        if n_accepted >= instances_per_cell:
            break
        inst, reason = cx.build_matched_instance(component_size, separation, seed, n_tiers)
        if inst is None:
            exclusions[reason] = exclusions.get(reason, 0) + 1
            continue
        n_accepted += 1

        s0_ok = is_gac_valid_mpdag(inst["g0"], inst["x"], inst["y"], inst["z_star"])
        if not s0_ok:
            s0_failures.append(inst["instance_id"])

        samples_tiered: list[dict[str, Any]] = []
        for rate in cx.TIERED_RATE_GRID:
            for rep in range(reps):
                row = cx.sample_tiered_state(inst, rate, rep)
                samples_tiered.append(row)
                samples_w.write(_fill(row, SAMPLES_FIELDS))

        n_k = inst["n_k"]
        samples_flip: list[dict[str, Any]] = []
        for d in range(1, n_k + 1):
            for rep in range(reps):
                row = cx.sample_flip_state(inst, d, rep)
                samples_flip.append(row)
                samples_w.write(_fill(row, SAMPLES_FIELDS))

        bins_tiered = cx.build_intensity_bins(samples_tiered)
        bins_flip = cx.build_intensity_bins(samples_flip)

        for b, stats in bins_tiered.items():
            brow = {
                "instance_id": inst["instance_id"],
                "arm": "tiered",
                "gate": cx.GATE_NAME,
                "component_size": component_size,
                "separation": separation,
                "n_tiers": n_tiers,
                "bin": b,
                **stats,
            }
            bins_w.write(_fill(brow, BINS_FIELDS))
        for b, stats in bins_flip.items():
            brow = {
                "instance_id": inst["instance_id"],
                "arm": "flip",
                "gate": cx.GATE_NAME,
                "component_size": component_size,
                "separation": separation,
                "n_tiers": n_tiers,
                "bin": b,
                **stats,
            }
            bins_w.write(_fill(brow, BINS_FIELDS))

        n_unresolved_tiered = 1 if "unresolved_all_contradictory" in bins_tiered else 0
        n_unresolved_flip = 1 if "unresolved_all_contradictory" in bins_flip else 0
        # Count of *grid points* folded into "unresolved", not just whether
        # the sentinel bin exists -- reconstructed from raw samples so the
        # exclusion count is exact, not just a yes/no flag.
        unresolved_gp_tiered = sum(
            1
            for rate in cx.TIERED_RATE_GRID
            if not any(r["status"] == "ok" and r["grid_point"] == rate for r in samples_tiered)
        )
        unresolved_gp_flip = sum(
            1
            for d in range(1, n_k + 1)
            if not any(r["status"] == "ok" and r["d"] == d for r in samples_flip)
        )

        rv: dict[str, Any] = {}
        if not skip_r_val:
            rv = cx.compute_r_val(inst["cpdag"], inst["g0"], inst["x"], inst["y"], inst["z_star"], time_limit_s=time_limit_s)

        irow = {
            "instance_id": inst["instance_id"],
            "base_instance_id": inst["base_instance_id"],
            "gate": cx.GATE_NAME,
            "component_size": component_size,
            "separation": separation,
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
            "n_dir_g0": inst["n_dir_g0"],
            "s0_ok": s0_ok,
            **rv,
            "n_grid_points_tiered": len(cx.TIERED_RATE_GRID),
            "n_grid_points_unresolved_tiered": unresolved_gp_tiered,
            "n_grid_points_flip": n_k,
            "n_grid_points_unresolved_flip": unresolved_gp_flip,
            "contradiction_rate_overall_tiered": cx.contradiction_rate_overall(samples_tiered),
            "contradiction_rate_overall_flip": cx.contradiction_rate_overall(samples_flip),
            "AUC_intensity_usable_tiered": cx.auc_intensity_usable(bins_tiered),
            "AUC_intensity_usable_flip": cx.auc_intensity_usable(bins_flip),
        }
        instances_w.write(_fill(irow, INSTANCES_FIELDS))
    return n_accepted


def cmd_sweep(args: argparse.Namespace) -> int:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    component_sizes = _parse_int_list(args.component_sizes)
    n_tiers_list = _parse_int_list(args.n_tiers)
    sep_override = _parse_int_list(args.separations) if args.separations else None

    samples_path = out_dir / "xarm_samples.csv"
    bins_path = out_dir / "xarm_bins.csv"
    instances_path = out_dir / "xarm_instances.csv"

    exclusions: dict[str, int] = {}
    s0_failures: list[str] = []
    cell_counts: dict[str, int] = {}

    t_start = time.perf_counter()
    with ResultWriter(samples_path, resume=args.resume) as samples_w, ResultWriter(
        bins_path, resume=args.resume
    ) as bins_w, ResultWriter(instances_path, resume=args.resume) as instances_w:
        for c in sorted(component_sizes):
            seps = sorted(sep_override) if sep_override else sorted(cx.achievable_separations(c))
            for s in seps:
                for nt in sorted(n_tiers_list):
                    n = run_cell(
                        c,
                        s,
                        nt,
                        instances_per_cell=args.instances_per_cell,
                        seed_budget=args.seed_budget,
                        reps=args.reps,
                        time_limit_s=args.time_limit_s,
                        skip_r_val=args.skip_r_val,
                        samples_w=samples_w,
                        bins_w=bins_w,
                        instances_w=instances_w,
                        exclusions=exclusions,
                        s0_failures=s0_failures,
                    )
                    cell_counts[f"c{c}_s{s}_nt{nt}"] = n
                    print(f"cell c{c}_s{s}_nt{nt}: accepted={n}", flush=True)

    elapsed = time.perf_counter() - t_start
    summary = {
        "elapsed_seconds": elapsed,
        "cell_counts": cell_counts,
        "exclusions": exclusions,
        "n_s0_failures": len(s0_failures),
        "s0_failures": s0_failures,
    }
    # NOTE: time.perf_counter() is monotonic but its epoch is arbitrary and
    # can coincide across separate process runs -- int(t_start) is NOT a
    # safe uniqueness key across chunked invocations. Use wall-clock time
    # (fine-grained) plus the PID so summaries from different chunks never
    # collide and overwrite each other.
    import os as _os

    tag = f"{int(time.time() * 1000)}_{_os.getpid()}"
    (out_dir / f"_xarm_sweep_summary_{tag}.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))
    return 0


# --- finalize --------------------------------------------------------------------


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def _num_or_blank(s: str) -> Any:
    if s == "":
        return ""
    try:
        if "." in s or "e" in s.lower():
            return float(s)
        return int(s)
    except ValueError:
        return s


def paired_test_from_bins(bin_rows: list[dict[str, Any]], *, min_n: int = cx.MIN_USABLE_N, n_boot: int = 2000) -> dict[Any, dict[str, Any]]:
    """Cross-arm paired test per shared intensity bin, restricted to matched cells.

    A cell (``instance_id``, ``bin``) counts only when *both* arms have
    ``n_noncontra >= min_n`` there. Uses a two-sided sign test (binomial,
    ``scipy.stats.binomtest``) and a paired bootstrap over instances, both
    seeded deterministically via :func:`bkrobust.robustness.crossarm.derived_seed`.
    """
    from scipy.stats import binomtest

    by_key: dict[tuple[str, Any], dict[str, tuple[float, int]]] = defaultdict(dict)
    for row in bin_rows:
        b = row["bin"]
        if not isinstance(b, (int, float)):
            continue
        if row["S"] == "" or row["n_noncontra"] == "":
            continue
        by_key[(row["instance_id"], b)][row["arm"]] = (float(row["S"]), int(row["n_noncontra"]))

    by_bin: dict[Any, list[tuple[str, float, float]]] = defaultdict(list)
    for (inst_id, b), arms in by_key.items():
        if "tiered" not in arms or "flip" not in arms:
            continue
        s_t, n_t = arms["tiered"]
        s_f, n_f = arms["flip"]
        if n_t < min_n or n_f < min_n:
            continue
        by_bin[b].append((inst_id, s_t, s_f))

    results: dict[Any, dict[str, Any]] = {}
    for b, rows in sorted(by_bin.items()):
        diffs = [s_t - s_f for (_, s_t, s_f) in rows]
        n = len(diffs)
        n_pos = sum(1 for d in diffs if d > 0)
        n_neg = sum(1 for d in diffs if d < 0)
        n_zero = n - n_pos - n_neg
        n_nonzero = n_pos + n_neg
        if n_nonzero > 0:
            sign_p = binomtest(min(n_pos, n_neg), n_nonzero, 0.5, alternative="two-sided").pvalue
        else:
            sign_p = ""
        if n > 0:
            rng = np.random.default_rng(cx.derived_seed("paired_bootstrap", b))
            boots = []
            for _ in range(n_boot):
                idx = rng.integers(0, n, size=n)
                boots.append(sum(diffs[i] for i in idx) / n)
            boots.sort()
            lo = boots[max(0, int(0.025 * n_boot) - 1)]
            hi = boots[min(n_boot - 1, int(0.975 * n_boot))]
            boot_mean = sum(boots) / len(boots)
            sorted_diffs = sorted(diffs)
            median_diff = sorted_diffs[n // 2] if n % 2 == 1 else (sorted_diffs[n // 2 - 1] + sorted_diffs[n // 2]) / 2
        else:
            lo = hi = boot_mean = median_diff = ""
        results[b] = {
            "n_instances_matched": n,
            "mean_diff_S_tiered_minus_flip": (sum(diffs) / n) if n else "",
            "median_diff": median_diff,
            "n_pos": n_pos,
            "n_neg": n_neg,
            "n_zero": n_zero,
            "sign_test_p": sign_p,
            "bootstrap_mean": boot_mean,
            "bootstrap_ci_lo": lo,
            "bootstrap_ci_hi": hi,
            "n_bootstrap": n_boot if n > 0 else 0,
        }
    return results


def cmd_finalize(args: argparse.Namespace) -> int:
    out_dir = Path(args.out_dir)
    bins_path = out_dir / "xarm_bins.csv"
    instances_path = out_dir / "xarm_instances.csv"
    samples_path = out_dir / "xarm_samples.csv"
    paired_path = out_dir / "xarm_paired_test.csv"

    bin_rows_raw = _read_csv_rows(bins_path)
    bin_rows: list[dict[str, Any]] = []
    for r in bin_rows_raw:
        row = dict(r)
        row["bin"] = _num_or_blank(r["bin"]) if r["bin"] != "unresolved_all_contradictory" else r["bin"]
        row["n_noncontra"] = _num_or_blank(r["n_noncontra"])
        row["S"] = _num_or_blank(r["S"])
        bin_rows.append(row)

    paired = paired_test_from_bins(bin_rows)
    if paired:
        with ResultWriter(paired_path, resume=False) as w:
            for b, stats in sorted(paired.items()):
                w.write(_fill({"bin": b, **stats}, PAIRED_TEST_FIELDS))
    else:
        # No bin anywhere had >= min_n non-contradictory samples in BOTH
        # arms for the same instance -- a real, reportable outcome (small
        # smoke runs, or a design where the arms genuinely never overlap on
        # the intensity axis), not an error. Still write a header-only file
        # rather than silently omitting the output.
        with paired_path.open("w", newline="") as fh:
            csv.DictWriter(fh, fieldnames=PAIRED_TEST_FIELDS).writeheader()

    # Achieved-intensity histogram per arm, pooled across all instances,
    # from the bins table (non-contradictory counts only; sentinel
    # "unresolved_all_contradictory" cells reported separately).
    hist: dict[str, dict[Any, int]] = {"tiered": defaultdict(int), "flip": defaultdict(int)}
    unresolved_cells: dict[str, int] = {"tiered": 0, "flip": 0}
    for row in bin_rows:
        arm = row["arm"]
        if row["bin"] == "unresolved_all_contradictory":
            unresolved_cells[arm] += 1
            continue
        hist[arm][row["bin"]] += int(row["n_noncontra"]) if row["n_noncontra"] != "" else 0

    instances_rows = _read_csv_rows(instances_path)
    n_instances = len(instances_rows)
    n_unreached = sum(1 for r in instances_rows if r.get("r_status") == "unreached")
    n_timeout = sum(1 for r in instances_rows if r.get("r_status") == "timeout")
    n_s0_bad = sum(1 for r in instances_rows if r.get("s0_ok") not in ("True", "1", True))

    samples_size = samples_path.stat().st_size if samples_path.exists() else 0
    gz_path = None
    if samples_size > 20 * 1024 * 1024:
        gz_path = out_dir / "xarm_samples.csv.gz"
        with samples_path.open("rb") as f_in, gzip.open(gz_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        samples_path.unlink()

    finalize_summary = {
        "n_instances": n_instances,
        "n_s0_failures": n_s0_bad,
        "n_r_val_unreached": n_unreached,
        "n_r_val_timeout": n_timeout,
        "achieved_intensity_histogram": {
            arm: {str(k): v for k, v in sorted(h.items(), key=lambda kv: kv[0])} for arm, h in hist.items()
        },
        "n_bin_cells_unresolved_all_contradictory": unresolved_cells,
        "paired_test_bins_with_data": sorted(str(b) for b in paired if paired[b]["n_instances_matched"] > 0),
        "samples_compressed": gz_path.name if gz_path else None,
    }
    (out_dir / "_xarm_finalize_summary.json").write_text(json.dumps(finalize_summary, indent=2, default=str))
    print(json.dumps(finalize_summary, indent=2, default=str))

    grid = {
        "component_sizes_seen": sorted({int(r["component_size"]) for r in instances_rows}) if instances_rows else [],
        "n_tiers_seen": sorted({int(r["n_tiers"]) for r in instances_rows}) if instances_rows else [],
        "tiered_rate_grid": list(cx.TIERED_RATE_GRID),
        "flip_grid": "d/len(K_ref) for d in 1..len(K_ref), per instance",
        "reps_per_grid_point": cx.N_REPS_DEFAULT,
        "bin_width": cx.BIN_WIDTH,
        "min_usable_n": cx.MIN_USABLE_N,
    }
    tmp_dir = out_dir / "_xarm_manifest_tmp"
    write_manifest(tmp_dir, seed=0, grid=grid, extra={"gate": cx.GATE_NAME, "design": "matched-instance cross-arm"})
    shutil.move(str(tmp_dir / "manifest.json"), str(out_dir / "xarm_manifest.json"))
    shutil.rmtree(tmp_dir, ignore_errors=True)

    return 0


# --- time-one --------------------------------------------------------------------


def cmd_time_one(args: argparse.Namespace) -> int:
    sep = args.separation
    if sep is None:
        seps = cx.achievable_separations(args.component_size)
        sep = seps[len(seps) // 2]

    t0 = time.perf_counter()
    inst = None
    reason = ""
    for seed in range(args.seed_budget):
        inst, reason = cx.build_matched_instance(args.component_size, sep, seed, args.n_tiers)
        if inst is not None:
            break
    t_build = time.perf_counter() - t0
    if inst is None:
        print(json.dumps({"error": "no instance admitted", "last_reason": reason, "seed_budget": args.seed_budget}))
        return 1

    t1 = time.perf_counter()
    n_k = inst["n_k"]
    n_tiered_samples = len(cx.TIERED_RATE_GRID) * args.reps
    n_flip_samples = n_k * args.reps
    for rate in cx.TIERED_RATE_GRID:
        for rep in range(args.reps):
            cx.sample_tiered_state(inst, rate, rep)
    t_tiered = time.perf_counter() - t1

    t2 = time.perf_counter()
    for d in range(1, n_k + 1):
        for rep in range(args.reps):
            cx.sample_flip_state(inst, d, rep)
    t_flip = time.perf_counter() - t2

    t3 = time.perf_counter()
    rv = cx.compute_r_val(inst["cpdag"], inst["g0"], inst["x"], inst["y"], inst["z_star"], time_limit_s=args.time_limit_s)
    t_rval = time.perf_counter() - t3

    total = time.perf_counter() - t0
    result = {
        "instance_id": inst["instance_id"],
        "n_k": n_k,
        "n_dir_g0": inst["n_dir_g0"],
        "n_tiered_samples": n_tiered_samples,
        "n_flip_samples": n_flip_samples,
        "t_build_s": t_build,
        "t_tiered_sampling_s": t_tiered,
        "t_flip_sampling_s": t_flip,
        "t_r_val_s": t_rval,
        "r_status": rv.get("r_status"),
        "r_val": rv.get("r_val"),
        "total_s": total,
        "per_sample_ms_tiered": 1000 * t_tiered / n_tiered_samples if n_tiered_samples else "",
        "per_sample_ms_flip": 1000 * t_flip / n_flip_samples if n_flip_samples else "",
    }
    print(json.dumps(result, indent=2, default=str))
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if args.cmd == "sweep":
        return cmd_sweep(args)
    if args.cmd == "finalize":
        return cmd_finalize(args)
    if args.cmd == "time-one":
        return cmd_time_one(args)
    raise ValueError(f"unknown cmd {args.cmd!r}")


if __name__ == "__main__":
    raise SystemExit(main())
