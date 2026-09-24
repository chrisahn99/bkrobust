"""Addendum to paired_synth_v2.py, fixing two issues the coordinator flagged
in review of the v2 output:

(A) Within-separation-strata mis-specification. paired_synth_v2.py's
``task2_separation`` stratified by (component_size, coverage, base_wrongness)
for flip / (component_size, n_tiers) for tiered -- i.e. it let separation
vary *within* each "stratum", the opposite of "holding separation fixed".
This module redefines the stratum as
``(paired_table_stratum, component_size, separation)`` -- the generator's own
admission cell (~40 instances each, drawn at one exact target separation) --
so separation truly is pinned within a stratum, and the only thing that
varies is the random seed draw. Only cells where r_val itself varies
contribute a defined tau_b (a constant column has no defined Kendall tau);
this is reported, not hidden.

(B) AUC support. Documents exactly how ``survival.auc_frac`` computes the
endpoint (a fixed 10-point *fractional* grid over n_k, each point filled by
the *nearest depth with a defined S* -- never a trapezoid, never a raw-depth
average -- see survival.py:559-583) and which grid points/depths tend to
have undefined S. It then builds, per paired-table stratum, a common-support
depth set (the intersection of {d : S(d) defined} across every instance in
the stratum, or -- when that intersection is empty -- the longest low-to-high
prefix of depths with >=90% per-depth definedness, dropping instances not
fully defined on that prefix) and recomputes AUC_frac restricted to search
only that common-support depth set (same nearest-neighbour algorithm,
candidate pool narrowed). tau_b(r_val, AUC_common) and the paired deltas vs.
separation/neg_phi_1/n_k/shd_truth are then recomputed on that common-support
population.

Writes (all NEW, under results/axis_robustness_p6_paired_v2/):
  - separation_within_strata_fixed.json  (task A)
  - auc_grid_diagnostics.json            (task B, part 1)
  - auc_common_support_paired.csv        (task B, part 2: paired comparisons)
  - auc_common_support_marginal.json     (task B, part 2: marginal tau_b + support sizes)

Usage::

    PYTHONPATH=src .venv/bin/python experiments/paired_synth_v2_addendum.py
"""

from __future__ import annotations

import csv
import json
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from bkrobust.robustness.paired_resample import (  # noqa: E402
    _kendall_tau_b,
    _to_float_or_none,
    paired_comparison,
)
from bkrobust.robustness.survival import FRAC_GRID  # noqa: E402

from paired_synth_v2 import (  # noqa: E402
    ENDPOINT,
    PAIRED_COLUMNS,
    degeneracy_note,
    load_instances,
    marginal_tau_bootstrap,
    stratified_tau_bootstrap,
    stratum_of,
)

IN_INSTANCES = ROOT / "results" / "axis_robustness_p6" / "survival_instances.csv"
IN_CURVES = ROOT / "results" / "axis_robustness_p6" / "survival_curves.csv"
OUT_DIR = ROOT / "results" / "axis_robustness_p6_paired_v2"

N_BOOT_COMMON_SUPPORT = 2000  # slower per-resample loop (curve rebuild); see FINDINGS note
SEED = 0
COMMON_SUPPORT_BASELINES = ["separation", "neg_phi_1", "n_k", "shd_truth"]


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in fieldnames})


# ---------------------------------------------------------------------------
# Task A: within-separation-strata, correctly holding separation fixed
# ---------------------------------------------------------------------------

def fixed_separation_key(row: dict[str, Any]) -> tuple:
    """(paired-table stratum, component_size, separation) -- the generator's
    own admission cell. Separation is pinned; only the seed draw varies."""
    return (stratum_of(row), row["component_size"], row["separation"])


def task_a(rows: list[dict[str, Any]]) -> dict[str, Any]:
    other_preds = ["n_k", "k_g0", "shd_truth", "neg_phi_1"]

    def group_stats(subset_rows):
        groups: dict[tuple, list] = defaultdict(list)
        for r in subset_rows:
            groups[fixed_separation_key(r)].append(r)
        n_groups = len(groups)
        varying = {k: v for k, v in groups.items()
                   if len(set(r["r_val"] for r in v)) > 1}
        n_instances_varying = sum(len(v) for v in varying.values())
        return n_groups, len(varying), n_instances_varying

    out: dict[str, Any] = {"by_arm": {}, "by_paired_table_stratum": {}}
    for arm in ("flip", "tiered"):
        arm_rows = [r for r in rows if r["arm"] == arm]
        n_groups, n_varying, n_inst_varying = group_stats(arm_rows)
        entry = {
            "n_instances_total": len(arm_rows),
            "n_groups_total": n_groups,
            "n_groups_with_rval_variation": n_varying,
            "n_instances_in_varying_groups": n_inst_varying,
            "r_val": stratified_tau_bootstrap(
                arm_rows, "r_val", ENDPOINT[arm], fixed_separation_key,
                n_boot=10000, seed=SEED),
        }
        for p in other_preds:
            entry[p] = stratified_tau_bootstrap(
                arm_rows, p, ENDPOINT[arm], fixed_separation_key,
                n_boot=10000, seed=SEED)
        out["by_arm"][arm] = entry

    strata: dict[str, list] = defaultdict(list)
    for r in rows:
        strata[stratum_of(r)].append(r)
    for s in sorted(strata):
        srows = strata[s]
        arm = srows[0]["arm"]
        n_groups, n_varying, n_inst_varying = group_stats(srows)
        key_fn_local = lambda r: (r["component_size"], r["separation"])
        entry = {
            "n_instances_total": len(srows),
            "n_groups_total": n_groups,
            "n_groups_with_rval_variation": n_varying,
            "n_instances_in_varying_groups": n_inst_varying,
            "r_val": stratified_tau_bootstrap(
                srows, "r_val", ENDPOINT[arm], key_fn_local, n_boot=10000, seed=SEED),
        }
        for p in other_preds:
            entry[p] = stratified_tau_bootstrap(
                srows, p, ENDPOINT[arm], key_fn_local, n_boot=10000, seed=SEED)
        out["by_paired_table_stratum"][s] = entry
    return out


# ---------------------------------------------------------------------------
# Task B: AUC support -- diagnostics + common-support recompute
# ---------------------------------------------------------------------------

def load_curves() -> dict[str, dict[int, dict[str, Any]]]:
    """{instance_id: {d: {"S": float|None, "n_eval": int}}}"""
    curves: dict[str, dict[int, dict[str, Any]]] = defaultdict(dict)
    with IN_CURVES.open() as fh:
        for row in csv.DictReader(fh):
            d = int(row["d"])
            s = row["S"]
            curves[row["instance_id"]][d] = {
                "S": None if s == "" else float(s),
                "n_eval": int(row["n_eval"]),
            }
    return curves


def grid_diagnostics(rows: list[dict[str, Any]], curves: dict[str, dict]) -> dict[str, Any]:
    """Which FRAC_GRID positions tend to hit an undefined raw depth, and
    which raw depths tend to be undefined, per arm."""
    n_k_by_inst = {r["instance_id"]: int(r["n_k"]) for r in rows}
    arm_by_inst = {r["instance_id"]: r["arm"] for r in rows}

    # (a) per raw-depth undefined rate, per arm (depth relabelled by rank
    # order 1..len since raw d is not comparable across instances of
    # different n_k -- see run_survival_p6.py: flip's d is 1..n_k, tiered's
    # d is the *measured* claim delta, irregular by construction).
    rank_stats: dict[tuple[str, int], list[int]] = defaultdict(list)  # (arm,rank) -> [0/1 undefined]
    for inst, curve in curves.items():
        arm = arm_by_inst.get(inst)
        if arm is None:
            continue
        for rank, d in enumerate(sorted(curve), start=1):
            rank_stats[(arm, rank)].append(0 if curve[d]["S"] is not None else 1)

    by_rank = []
    for (arm, rank), vals in sorted(rank_stats.items()):
        if len(vals) < 5:
            continue
        by_rank.append({"arm": arm, "depth_rank": rank, "n": len(vals),
                         "frac_undefined": sum(vals) / len(vals)})

    # (b) per FRAC_GRID position, whether the EXACT target depth
    # (round(frac*n_k)) is itself defined -- i.e. how often the nearest-
    # neighbour fallback has to reach elsewhere.
    frac_stats: dict[tuple[str, float], list[int]] = defaultdict(list)
    for inst, curve in curves.items():
        arm = arm_by_inst.get(inst)
        n_k = n_k_by_inst.get(inst)
        if arm is None or not n_k:
            continue
        defined_depths = sorted(d for d, v in curve.items() if v["S"] is not None)
        for frac in FRAC_GRID:
            target = frac * n_k
            # auc_frac picks the nearest depth with a defined S to `target`;
            # a nonzero distance here means the exact target position had no
            # defined S and the algorithm had to fall back to a neighbour.
            if defined_depths:
                nearest = min(defined_depths, key=lambda d: (abs(d - target), d))
                dist = abs(nearest - target)
            else:
                dist = None
            frac_stats[(arm, frac)].append(1 if (dist is None or dist > 0) else 0)

    by_frac = []
    for (arm, frac), vals in sorted(frac_stats.items()):
        by_frac.append({"arm": arm, "frac": frac, "n": len(vals),
                         "frac_needed_fallback": sum(vals) / len(vals)})

    return {"by_depth_rank": by_rank, "by_frac_grid_position": by_frac,
            "note": "depth_rank = 1st,2nd,... smallest measured depth for that "
                    "instance (raw d not comparable across instances -- flip's d "
                    "is 1..n_k, tiered's d is a measured, irregular claim delta); "
                    "frac_needed_fallback = fraction of instances where the exact "
                    "target depth round(frac*n_k) was NOT itself a defined S, so "
                    "auc_frac's nearest-neighbour search had to move to another depth."}


def common_support_depths(stratum_rows: list[dict[str, Any]],
                           curves: dict[str, dict]) -> tuple[list[int], list[str], str]:
    """Returns (common_depths, kept_instance_ids, method_note)."""
    inst_ids = [r["instance_id"] for r in stratum_rows if r["instance_id"] in curves]
    defined_sets = {iid: set(d for d, v in curves[iid].items() if v["S"] is not None)
                     for iid in inst_ids}
    if not defined_sets:
        return [], [], "no curves found"
    intersection = set.intersection(*defined_sets.values())
    if intersection:
        kept = inst_ids
        return sorted(intersection), kept, "exact intersection across all instances"

    # Fallback: longest low-to-high prefix of depths with >=90% per-depth
    # definedness (checked depth by depth, independent of prefix so far);
    # instances not fully defined on the resulting prefix are dropped.
    all_depths = sorted(set().union(*defined_sets.values()) | set().union(
        *(set(curves[iid].keys()) for iid in inst_ids)))
    n = len(inst_ids)
    prefix: list[int] = []
    for d in all_depths:
        n_defined_here = sum(1 for iid in inst_ids if d in curves[iid] and curves[iid][d]["S"] is not None)
        frac_defined = n_defined_here / n if n else 0.0
        if frac_defined >= 0.90:
            prefix.append(d)
        else:
            break
    kept = [iid for iid in inst_ids if all(d in defined_sets[iid] for d in prefix)]
    return prefix, kept, f"90%-coverage prefix (len={len(prefix)}), {len(kept)}/{n} instances retained"


def auc_on_support(curve: dict[int, dict[str, Any]], n_k: int, support: list[int]) -> float | None:
    """auc_frac's own nearest-neighbour algorithm, candidate pool narrowed to
    `support` (a subset of depths, all defined by construction of the caller
    -- see common_support_depths)."""
    defined = [d for d in support if d in curve and curve[d]["S"] is not None]
    if not defined or n_k <= 0:
        return None
    vals = []
    for frac in FRAC_GRID:
        target = frac * n_k
        nearest = min(defined, key=lambda d: (abs(d - target), d))
        vals.append(curve[nearest]["S"])
    return sum(vals) / len(vals)


def task_b_common_support(rows: list[dict[str, Any]], curves: dict[str, dict]) -> dict[str, Any]:
    strata: dict[str, list] = defaultdict(list)
    for r in rows:
        strata[stratum_of(r)].append(r)

    marginal_out: dict[str, Any] = {}
    paired_rows: list[dict[str, Any]] = []
    rows_by_instance = {r["instance_id"]: r for r in rows}

    for s in sorted(strata):
        srows = strata[s]
        arm = srows[0]["arm"]
        support, kept_ids, note = common_support_depths(srows, curves)
        kept_rows = []
        for iid in kept_ids:
            r = dict(rows_by_instance[iid])
            n_k = int(r["n_k"])
            auc = auc_on_support(curves[iid], n_k, support)
            r["AUC_common"] = "" if auc is None else repr(auc)
            kept_rows.append(r)
        marginal_out[s] = {
            "arm": arm, "stratum_n": len(srows), "n_common_support_depths": len(support),
            "common_support_depths": support, "n_kept": len(kept_rows), "note": note,
        }
        res = marginal_tau_bootstrap(kept_rows, "r_val", "AUC_common", n_boot=N_BOOT_COMMON_SUPPORT, seed=SEED)
        marginal_out[s]["marginal_tau_r_val"] = res

        for b in COMMON_SUPPORT_BASELINES:
            pres = paired_comparison(kept_rows, "r_val", b, "AUC_common", "instance_id",
                                      n_boot=N_BOOT_COMMON_SUPPORT, seed=SEED)
            paired_rows.append({
                "stratum": s, "arm": arm, "predictor_a": "r_val", "predictor_b": b,
                "endpoint": "AUC_common", "status": pres.get("status", "ok"),
                "stratum_n": len(srows), "n": pres.get("n"), "n_excluded": pres.get("n_excluded"),
                "n_clusters": pres.get("n_clusters"), "tau_a_point": pres.get("tau_a_point"),
                "tau_b_point": pres.get("tau_b_point"), "delta_point": pres.get("delta_point"),
                "ci_lo_2p5": pres.get("ci_lo_2p5"), "ci_hi_97p5": pres.get("ci_hi_97p5"),
                "frac_delta_gt_0": pres.get("frac_delta_gt_0"), "n_boot": N_BOOT_COMMON_SUPPORT,
                "seed": SEED, "degeneracy_note": degeneracy_note(kept_rows, b) if kept_rows else "",
            })
    return {"marginal": marginal_out, "paired": paired_rows}


def main() -> None:
    t0 = time.perf_counter()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = load_instances()
    print(f"[addendum] loaded {len(rows)} instances", flush=True)

    print("[addendum] Task A: within-separation-strata (separation held fixed)", flush=True)
    task_a_out = task_a(rows)
    (OUT_DIR / "separation_within_strata_fixed.json").write_text(json.dumps(task_a_out, indent=1, default=str))

    print("[addendum] Task B: loading curves, grid diagnostics", flush=True)
    curves = load_curves()
    diag = grid_diagnostics(rows, curves)
    (OUT_DIR / "auc_grid_diagnostics.json").write_text(json.dumps(diag, indent=1, default=str))

    print(f"[addendum] Task B: common-support recompute (n_boot={N_BOOT_COMMON_SUPPORT})", flush=True)
    cs = task_b_common_support(rows, curves)
    (OUT_DIR / "auc_common_support_marginal.json").write_text(
        json.dumps(cs["marginal"], indent=1, default=str))
    write_csv(OUT_DIR / "auc_common_support_paired.csv", PAIRED_COLUMNS, cs["paired"])

    print(f"[addendum] done in {time.perf_counter() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
