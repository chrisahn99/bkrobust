"""When r_val and separation DISAGREE, which one ranks survival correctly?

Context: on the committed 4,104-instance synthetic sweep
(``results/axis_robustness_p6/survival_instances.csv``), separation (the
generator's own design variable) ranks survival about as well as r_val
overall, and r_val != separation on 55.7% of instances (see
``results/axis_robustness_p6_paired_v2/FINDINGS.md`` section 2, computed by
``experiments/paired_synth_v2.py``). That leaves open the question this
script answers: on the instances/pairs where the two predictors disagree,
which one actually tracks survival?

Four analyses, per paired-table stratum (flip cov x bw, tiered nt2/3/4) and
pooled per arm, all against the ``AUC_frac`` endpoint, all using the exact
instance-clustered bootstrap machinery already committed in
``src/bkrobust/robustness/paired_resample.py``
(:func:`paired_comparison`) or ``experiments/paired_synth_v2.py``
(:func:`stratified_tau_bootstrap`), n_boot=10000/seed=0 unless noted:

A. Discordant-subset paired tau_b: restrict to instances with
   ``r_val != separation`` (exact, same string-equality test
   ``paired_synth_v2.separation_neq_facts`` already uses for the 55.7%
   figure); ``paired_comparison(rows, "r_val", "separation", "AUC_frac",
   "instance_id")`` on that subset -- the paired delta IS the "who wins
   among the disagreers" statistic at the instance level.

B. Pair-level test (does not condition on a selected instance subset): for
   every stratum (and pooled per arm), over all instance pairs (i, j) whose
   r_val-order and separation-order of {i, j} DISAGREE (reported both as
   "strictly opposite" -- both predictors give a strict, opposite verdict --
   and "any disagreement" -- includes a tie in one predictor while the other
   is strict), the fraction of such pairs where the survival-AUC order
   agrees with r_val's order vs. with separation's order (a third "neither"
   bucket catches an AUC tie or a match with neither -- see docstring of
   :func:`pairwise_stats`). This is the cleanest "who wins when they
   disagree" statistic: it does not throw away information the instance-
   level equality test discards (e.g. a pair where the two INSTANCES have
   r_val_i == r_val_j but separation_i != separation_j is a genuine
   disagreement about the ORDER of that pair, even though neither instance
   individually satisfies "r_val != separation"). CI by resampling
   instances (cluster) with replacement, 2,000 times (a full O(n^2) pair
   matrix rebuilt per resample; see runtime note in the module docstring
   below), same seed=0.

C. Distribution of r_val - separation (on the "measured" separation-status
   population), by arm: is r_val usually larger or smaller, and by how much.

D. Partial association, the symmetric counterpart to the already-committed
   tau_b(r_val, AUC | separation) (``results/axis_robustness_p6_paired_v2/
   separation_within_strata_fixed.json``, section (A) of FINDINGS.md there):
   tau_b(separation, AUC | r_val), stratifying on the EXACT r_val value
   within each paired-table design cell (the generator's own admission cell
   swapped role for role), pair-weighted, instance-cluster bootstrap CI --
   using the identical ``stratified_tau_bootstrap`` helper, only with
   predictor and stratify-key swapped.

Runtime note (Task B): a full n x n pairwise-comparison matrix is rebuilt
on every bootstrap resample (no O(n log n) shortcut exists for the joint
3-way statistic used here -- agreement of AUC's order with each of two
OTHER predictors' orders, restricted to pairs where those two disagree).
For the pooled flip arm (n=2,730) this is a ~3.7M-pair upper triangle,
rebuilt 2,000 times; still tractable in vectorized numpy (single-digit
minutes), see ``manifest.json:task_b_seconds``.

Outputs (all NEW, under
``results/axis_robustness_p6_paired_v2/discordant/``):
  - discordant_subset_paired.csv    (task A)
  - pairwise_disagreement.csv       (task B)
  - diff_distribution.json          (task C)
  - rval_conditioned_separation.json (task D)
  - manifest.json
  - FINDINGS.md (written separately, not by this script)

Usage::

    PYTHONPATH=src .venv/bin/python experiments/discordant_separation_v2.py
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

from bkrobust.robustness.paired_resample import (  # noqa: E402
    _to_float_or_none,
    paired_comparison,
)

from paired_synth_v2 import (  # noqa: E402
    ENDPOINT,
    PAIRED_COLUMNS,
    load_instances,
    stratified_tau_bootstrap,
    stratum_of,
)

IN_INSTANCES = ROOT / "results" / "axis_robustness_p6" / "survival_instances.csv"
OUT_DIR = ROOT / "results" / "axis_robustness_p6_paired_v2" / "discordant"

N_BOOT_A = 10000
N_BOOT_B = 2000
N_BOOT_D = 10000
SEED = 0


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in fieldnames})


# ---------------------------------------------------------------------------
# Task A: discordant-subset (r_val != separation) paired tau_b
# ---------------------------------------------------------------------------

def is_discordant_row(r: dict[str, Any]) -> bool:
    """Same exact-string-inequality test paired_synth_v2.separation_neq_facts
    uses for the committed 55.7% figure: r_val and separation are both
    present (separation_status == 'measured') and their raw CSV values are
    not textually identical (a zero-decimal identity is exact by
    construction of the generator, so string equality is the right test --
    matching the existing committed figure exactly rather than introducing a
    new float-tolerance rule)."""
    if r.get("separation_status") != "measured":
        return False
    return r["r_val"] != r["separation"]


def task_a(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    discordant = [r for r in rows if is_discordant_row(r)]
    out_rows = []

    strata: dict[str, list[dict[str, Any]]] = {}
    for r in discordant:
        strata.setdefault(stratum_of(r), []).append(r)

    for stratum in sorted(strata):
        srows = strata[stratum]
        arm = srows[0]["arm"]
        endpoint = ENDPOINT[arm]
        res = paired_comparison(srows, "r_val", "separation", endpoint, "instance_id",
                                 n_boot=N_BOOT_A, seed=SEED)
        out_rows.append(_paired_row(stratum, arm, endpoint, srows, res))

    for arm in ("flip", "tiered"):
        arm_rows = [r for r in discordant if r["arm"] == arm]
        endpoint = ENDPOINT[arm]
        res = paired_comparison(arm_rows, "r_val", "separation", endpoint, "instance_id",
                                 n_boot=N_BOOT_A, seed=SEED)
        out_rows.append(_paired_row(f"{arm}_POOLED", arm, endpoint, arm_rows, res))

    return out_rows


def _paired_row(stratum: str, arm: str, endpoint: str, srows: list, res: dict) -> dict[str, Any]:
    return {
        "stratum": stratum, "arm": arm, "predictor_a": "r_val", "predictor_b": "separation",
        "endpoint": endpoint, "status": res.get("status", "ok"),
        "stratum_n": len(srows), "n": res.get("n"), "n_excluded": res.get("n_excluded"),
        "n_clusters": res.get("n_clusters"), "tau_a_point": res.get("tau_a_point"),
        "tau_b_point": res.get("tau_b_point"), "delta_point": res.get("delta_point"),
        "ci_lo_2p5": res.get("ci_lo_2p5"), "ci_hi_97p5": res.get("ci_hi_97p5"),
        "frac_delta_gt_0": res.get("frac_delta_gt_0"), "n_boot": N_BOOT_A, "seed": SEED,
        "degeneracy_note": (
            "r_val outranks separation among discordant instances" if (
                res.get("delta_point") is not None and res.get("ci_lo_2p5") is not None
                and res["ci_lo_2p5"] > 0)
            else "separation outranks r_val among discordant instances" if (
                res.get("delta_point") is not None and res.get("ci_hi_97p5") is not None
                and res["ci_hi_97p5"] < 0)
            else ""
        ),
    }


# ---------------------------------------------------------------------------
# Task B: pair-level "who wins when they order a pair differently" test
# ---------------------------------------------------------------------------

def pairwise_stats(r: np.ndarray, s: np.ndarray, a: np.ndarray,
                    iu: np.ndarray, ju: np.ndarray) -> dict[str, dict]:
    """For every pair (iu[k], ju[k]) -- the upper triangle, i<j, no self-pairs
    -- compares the sign of (r_i-r_j), (s_i-s_j), (a_i-a_j). Returns
    {"strict_opposite": {...}, "any_disagreement": {...}}, each with n_pairs,
    n_agree_r, n_agree_sep, n_neither and their fractions of n_pairs.

    "strict_opposite": r_val and separation ORDER the pair with strict,
    opposite verdicts (sr, ss) = (+1,-1) or (-1,+1) -- neither predictor is
    tied on this pair, and they disagree about which of i/j is larger.

    "any_disagreement": sr != ss, which additionally admits a pair where one
    predictor is tied (sign 0) and the other is not.

    A pair "agrees with r_val" (resp. separation) iff the AUC sign on that
    pair equals r_val's (resp. separation's) sign AND that predictor's sign
    is itself non-zero (a tied predictor cannot "win" a pair). Because the
    two predictors disagree by construction of either mask (sr != ss), a
    pair can agree with at most one of them; anything else (an AUC tie, or
    an AUC sign matching neither disagreeing predictor) falls into
    "neither".
    """
    dr = r[iu] - r[ju]
    ds = s[iu] - s[ju]
    da = a[iu] - a[ju]
    sr = np.sign(dr)
    ss = np.sign(ds)
    sa = np.sign(da)

    disagree_any = sr != ss
    disagree_strict = disagree_any & (sr != 0) & (ss != 0)

    def frac_stats(mask: np.ndarray) -> dict[str, Any]:
        n_pairs = int(mask.sum())
        if n_pairs == 0:
            return {"n_pairs": 0, "n_agree_r": 0, "n_agree_sep": 0, "n_neither": 0,
                     "frac_agree_r": None, "frac_agree_sep": None, "frac_neither": None}
        sub_sr, sub_ss, sub_sa = sr[mask], ss[mask], sa[mask]
        agree_r = (sub_sr != 0) & (sub_sa == sub_sr)
        agree_sep = (sub_ss != 0) & (sub_sa == sub_ss)
        n_r = int(agree_r.sum())
        n_s = int(agree_sep.sum())
        n_neither = n_pairs - n_r - n_s
        return {
            "n_pairs": n_pairs, "n_agree_r": n_r, "n_agree_sep": n_s, "n_neither": n_neither,
            "frac_agree_r": n_r / n_pairs, "frac_agree_sep": n_s / n_pairs,
            "frac_neither": n_neither / n_pairs,
        }

    return {"strict_opposite": frac_stats(disagree_strict), "any_disagreement": frac_stats(disagree_any)}


def task_b_one_group(rows: list[dict[str, Any]], endpoint: str,
                      n_boot: int = N_BOOT_B, seed: int = SEED) -> dict[str, Any]:
    cleaned = []
    n_excluded = 0
    for row in rows:
        rv = _to_float_or_none(row.get("r_val"))
        sv = _to_float_or_none(row.get("separation")) if row.get("separation_status") == "measured" else None
        av = _to_float_or_none(row.get(endpoint))
        if rv is None or sv is None or av is None:
            n_excluded += 1
            continue
        cleaned.append((rv, sv, av))
    out: dict[str, Any] = {"stratum_n": len(rows), "n": len(cleaned), "n_excluded": n_excluded,
                            "n_boot": n_boot, "seed": seed, "status": "ok"}
    if len(cleaned) < 2:
        out["status"] = "insufficient_n"
        return out

    r = np.array([c[0] for c in cleaned])
    s = np.array([c[1] for c in cleaned])
    a = np.array([c[2] for c in cleaned])
    n = len(cleaned)
    iu, ju = np.triu_indices(n, k=1)

    point = pairwise_stats(r, s, a, iu, ju)
    out["point"] = point

    rng = np.random.default_rng(seed)
    boots = {
        "strict_opposite": {"frac_agree_r": [], "frac_agree_sep": []},
        "any_disagreement": {"frac_agree_r": [], "frac_agree_sep": []},
    }
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        rb, sb, ab = r[idx], s[idx], a[idx]
        stats = pairwise_stats(rb, sb, ab, iu, ju)
        for key in ("strict_opposite", "any_disagreement"):
            fr = stats[key]["frac_agree_r"]
            fs = stats[key]["frac_agree_sep"]
            if fr is not None:
                boots[key]["frac_agree_r"].append(fr)
            if fs is not None:
                boots[key]["frac_agree_sep"].append(fs)

    ci: dict[str, Any] = {}
    for key in ("strict_opposite", "any_disagreement"):
        ci[key] = {}
        for stat_name in ("frac_agree_r", "frac_agree_sep"):
            vals = np.array(boots[key][stat_name])
            if len(vals) >= 100:
                lo, hi = np.percentile(vals, [2.5, 97.5])
                ci[key][stat_name] = {"ci_lo_2p5": float(lo), "ci_hi_97p5": float(hi),
                                       "n_boot_valid": int(len(vals))}
            else:
                ci[key][stat_name] = {"ci_lo_2p5": None, "ci_hi_97p5": None,
                                       "n_boot_valid": int(len(vals))}
    out["ci"] = ci
    return out


def task_b(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out_rows = []
    strata: dict[str, list[dict[str, Any]]] = {}
    for r in rows:
        strata.setdefault(stratum_of(r), []).append(r)

    def flatten(stratum: str, arm: str, res: dict[str, Any]) -> dict[str, Any]:
        row = {"stratum": stratum, "arm": arm, "stratum_n": res.get("stratum_n"),
               "n": res.get("n"), "n_excluded": res.get("n_excluded"),
               "status": res.get("status"), "n_boot": res.get("n_boot"), "seed": res.get("seed")}
        if res.get("status") != "ok":
            return row
        for key in ("strict_opposite", "any_disagreement"):
            p = res["point"][key]
            c = res["ci"][key]
            prefix = "strict" if key == "strict_opposite" else "any"
            row[f"{prefix}_n_pairs"] = p["n_pairs"]
            row[f"{prefix}_frac_agree_r"] = p["frac_agree_r"]
            row[f"{prefix}_frac_agree_r_ci_lo"] = c["frac_agree_r"]["ci_lo_2p5"]
            row[f"{prefix}_frac_agree_r_ci_hi"] = c["frac_agree_r"]["ci_hi_97p5"]
            row[f"{prefix}_frac_agree_sep"] = p["frac_agree_sep"]
            row[f"{prefix}_frac_agree_sep_ci_lo"] = c["frac_agree_sep"]["ci_lo_2p5"]
            row[f"{prefix}_frac_agree_sep_ci_hi"] = c["frac_agree_sep"]["ci_hi_97p5"]
            row[f"{prefix}_frac_neither"] = p["frac_neither"]
        return row

    for stratum in sorted(strata):
        srows = strata[stratum]
        arm = srows[0]["arm"]
        endpoint = ENDPOINT[arm]
        res = task_b_one_group(srows, endpoint)
        out_rows.append(flatten(stratum, arm, res))

    for arm in ("flip", "tiered"):
        arm_rows = [r for r in rows if r["arm"] == arm]
        endpoint = ENDPOINT[arm]
        res = task_b_one_group(arm_rows, endpoint)
        out_rows.append(flatten(f"{arm}_POOLED", arm, res))

    return out_rows


# ---------------------------------------------------------------------------
# Task C: distribution of r_val - separation, by arm
# ---------------------------------------------------------------------------

def task_c(rows: list[dict[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for arm in ("flip", "tiered", "all"):
        arm_rows = [r for r in rows if (arm == "all" or r["arm"] == arm)
                    and r.get("separation_status") == "measured"]
        diffs = []
        n_gt, n_lt, n_eq = 0, 0, 0
        for r in arm_rows:
            rv = _to_float_or_none(r["r_val"])
            sv = _to_float_or_none(r["separation"])
            if rv is None or sv is None:
                continue
            d = rv - sv
            diffs.append(d)
            if d > 0:
                n_gt += 1
            elif d < 0:
                n_lt += 1
            else:
                n_eq += 1
        n = len(diffs)
        entry: dict[str, Any] = {"n": n, "n_r_gt_sep": n_gt, "n_r_lt_sep": n_lt, "n_r_eq_sep": n_eq,
                                  "frac_r_gt_sep": n_gt / n if n else None,
                                  "frac_r_lt_sep": n_lt / n if n else None,
                                  "frac_r_eq_sep": n_eq / n if n else None}
        if n:
            arr = np.array(diffs)
            entry.update({
                "mean": float(np.mean(arr)), "median": float(np.median(arr)),
                "std": float(np.std(arr, ddof=1)) if n > 1 else 0.0,
                "q25": float(np.percentile(arr, 25)), "q75": float(np.percentile(arr, 75)),
                "min": float(np.min(arr)), "max": float(np.max(arr)),
            })
        out[arm] = entry
    return out


# ---------------------------------------------------------------------------
# Task D: tau_b(separation, AUC | r_val) -- symmetric to the already-
# committed tau_b(r_val, AUC | separation) (separation_within_strata_fixed.json)
# ---------------------------------------------------------------------------

def fixed_r_val_key(row: dict[str, Any]) -> tuple:
    """(paired-table stratum, component_size, r_val) -- the mirror image of
    paired_synth_v2_addendum.fixed_separation_key, with predictor and
    stratify-key roles swapped: here r_val (not separation) is pinned
    exactly, and separation is free to vary within the cell."""
    return (stratum_of(row), row["component_size"], row["r_val"])


def task_d(rows: list[dict[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {"by_arm": {}, "by_paired_table_stratum": {}}
    for arm in ("flip", "tiered"):
        arm_rows = [r for r in rows if r["arm"] == arm]
        groups: dict[tuple, list] = {}
        for r in arm_rows:
            groups.setdefault(fixed_r_val_key(r), []).append(r)
        n_groups = len(groups)
        varying = {k: v for k, v in groups.items() if len(set(r["separation"] for r in v)) > 1}
        n_inst_varying = sum(len(v) for v in varying.values())
        out["by_arm"][arm] = {
            "n_instances_total": len(arm_rows),
            "n_groups_total": n_groups,
            "n_groups_with_separation_variation": len(varying),
            "n_instances_in_varying_groups": n_inst_varying,
            "separation_given_r_val": stratified_tau_bootstrap(
                arm_rows, "separation", ENDPOINT[arm], fixed_r_val_key,
                n_boot=N_BOOT_D, seed=SEED),
        }

    strata: dict[str, list] = {}
    for r in rows:
        strata.setdefault(stratum_of(r), []).append(r)
    for s in sorted(strata):
        srows = strata[s]
        arm = srows[0]["arm"]
        key_fn_local = lambda r: (r["component_size"], r["r_val"])
        groups: dict[tuple, list] = {}
        for r in srows:
            groups.setdefault(key_fn_local(r), []).append(r)
        varying = {k: v for k, v in groups.items() if len(set(r["separation"] for r in v)) > 1}
        n_inst_varying = sum(len(v) for v in varying.values())
        out["by_paired_table_stratum"][s] = {
            "n_instances_total": len(srows),
            "n_groups_total": len(groups),
            "n_groups_with_separation_variation": len(varying),
            "n_instances_in_varying_groups": n_inst_varying,
            "separation_given_r_val": stratified_tau_bootstrap(
                srows, "separation", ENDPOINT[arm], key_fn_local, n_boot=N_BOOT_D, seed=SEED),
        }
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    t0 = time.perf_counter()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows = load_instances()
    print(f"[discordant] loaded {len(rows)} instances", flush=True)

    n_discordant = sum(1 for r in rows if is_discordant_row(r))
    print(f"[discordant] {n_discordant} discordant (r_val != separation) instances "
          f"({n_discordant / len(rows):.1%})", flush=True)

    print("[discordant] task A: discordant-subset paired tau_b", flush=True)
    ta0 = time.perf_counter()
    task_a_rows = task_a(rows)
    write_csv(OUT_DIR / "discordant_subset_paired.csv", PAIRED_COLUMNS, task_a_rows)
    task_a_seconds = time.perf_counter() - ta0

    print("[discordant] task B: pair-level disagreement test", flush=True)
    tb0 = time.perf_counter()
    task_b_rows = task_b(rows)
    b_fields = ["stratum", "arm", "stratum_n", "n", "n_excluded", "status", "n_boot", "seed",
                "strict_n_pairs", "strict_frac_agree_r", "strict_frac_agree_r_ci_lo",
                "strict_frac_agree_r_ci_hi", "strict_frac_agree_sep", "strict_frac_agree_sep_ci_lo",
                "strict_frac_agree_sep_ci_hi", "strict_frac_neither",
                "any_n_pairs", "any_frac_agree_r", "any_frac_agree_r_ci_lo", "any_frac_agree_r_ci_hi",
                "any_frac_agree_sep", "any_frac_agree_sep_ci_lo", "any_frac_agree_sep_ci_hi",
                "any_frac_neither"]
    write_csv(OUT_DIR / "pairwise_disagreement.csv", b_fields, task_b_rows)
    task_b_seconds = time.perf_counter() - tb0

    print("[discordant] task C: r_val - separation distribution", flush=True)
    task_c_out = task_c(rows)
    (OUT_DIR / "diff_distribution.json").write_text(json.dumps(task_c_out, indent=1, default=str))

    print("[discordant] task D: tau_b(separation, AUC | r_val)", flush=True)
    td0 = time.perf_counter()
    task_d_out = task_d(rows)
    (OUT_DIR / "rval_conditioned_separation.json").write_text(json.dumps(task_d_out, indent=1, default=str))
    task_d_seconds = time.perf_counter() - td0

    manifest = {
        "inputs": {
            "survival_instances_csv": {
                "path": str(IN_INSTANCES.relative_to(ROOT)),
                "sha256": sha256_file(IN_INSTANCES), "n_rows": len(rows),
            },
        },
        "n_instances": len(rows),
        "n_discordant_instances": n_discordant,
        "frac_discordant_instances": n_discordant / len(rows),
        "n_boot_task_a": N_BOOT_A, "n_boot_task_b": N_BOOT_B, "n_boot_task_d": N_BOOT_D,
        "seed": SEED,
        "task_a_seconds": task_a_seconds,
        "task_b_seconds": task_b_seconds,
        "task_d_seconds": task_d_seconds,
        "elapsed_seconds": time.perf_counter() - t0,
    }
    (OUT_DIR / "manifest.json").write_text(json.dumps(manifest, indent=1, default=str))
    print(f"[discordant] done in {manifest['elapsed_seconds']:.1f}s", flush=True)


if __name__ == "__main__":
    main()
