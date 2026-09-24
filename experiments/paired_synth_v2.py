"""Follow-up analysis on top of the committed synthetic sweep and paired
bootstrap (reviewer Tier 0/1 fixes: phi_1 direction, separation baseline
disclosure, current-sweep separation facts, AUC support, and a check for a
"consistent draw left closure orientations unchanged" signal).

Everything this script writes is NEW, under
``results/axis_robustness_p6_paired_v2/``. It never modifies anything under
``results/axis_robustness_p6/`` or ``results/axis_robustness_p6_paired/``
(both read-only inputs here).

Task 1 -- direction fix
------------------------
``phi_1`` (src/bkrobust/robustness/survival_p6.py) is a FRAGILITY score:
larger phi_1 means *more* single-claim fragility, so a good radius predictor
should be *negatively* associated with it. The committed paired table
(``results/axis_robustness_p6_paired/paired_diff_synth.csv``) computed
``tau_b(r_val) - tau_b(phi_1)`` directly, comparing r_val's positive
association with survival against phi_1's own (uninverted) sign -- the wrong
direction per the task brief. This script adds ``neg_phi_1 = -phi_1`` and
reruns ``paired_comparison`` (src/bkrobust/robustness/paired_resample.py, the
exact same function, same n_boot=10000, seed=0, same instance-clustered
design) with ``neg_phi_1`` in place of ``phi_1``, holding every other
baseline (shd_truth, n_k, k_g0, separation, r_claim) and every setting fixed.
Before trusting the rerun, :func:`verify_unchanged_baselines` diffs its own
output for those five untouched baselines against the committed CSV,
row for row, and reports any discrepancy rather than assuming reproduction.

Task 2 -- separation facts on the current 4,104-instance sweep
-----------------------------------------------------------------
The paper's Appendix G still quotes stale numbers from an old 1,978-instance
sweep ("radius differs from separation on 56%", tau_b=0.55/0.24 within
separation strata). No script in this repo recomputes those figures on the
current 4,104-instance sweep, and no script under ``src/bkrobust/robustness``
or ``experiments`` was found (by name or by grep for "1978"/"56%"/"separation
strata"/"stratified tau") that implements the *old* sweep's within-strata
combination method -- the old dataset and its analysis script are gone from
this repo. Per the task brief's explicit fallback, this script defines the
within-separation-strata tau_b as the pair-count-weighted average of
per-stratum Kendall tau_b(r_val, AUC_frac), with an instance-level cluster
bootstrap CI recomputed by re-forming the same weighted average on each
resample (see :func:`stratified_tau_bootstrap`). A "separation stratum" here
is (component_size, coverage, base_wrongness) for the flip arm and
(component_size, n_tiers) for the tiered arm -- the design cell that "holds
the design fixed" per the paper's own wording (\"strata of separation, which
hold the design fixed\"): everything about the generator's target is pinned
except the swept separation level (min/median/max achievable), so tau within
a stratum asks whether r_val still orders survival *beyond* what the
component-size/coverage/tier design alone would predict.

Task 3 -- AUC support (reviewer 4.4)
--------------------------------------
Uses ``n_depths_evaluated`` / ``n_depths_all_contradictory`` from
``survival_instances.csv`` (already computed by ``run_survival_p6.py``) to
report, per arm/stratum, the share of instances with at least one grid point
where every draw was contradictory (undefined survival at that depth), and
the distribution of the number of *defined* grid points
(``n_depths_evaluated - n_depths_all_contradictory``). It then recomputes
tau_b(r_val, AUC_frac) and the paired differences vs the five baselines,
restricted to instances with FULL support (``n_depths_all_contradictory ==
0``), and reports whether the sign/significance of the committed conclusions
changes. Accepted-draw counts per intensity (median/min across instances) per
arm come straight from ``survival_curves.csv``'s own ``n_eval`` column (the
count of non-contradictory, evaluable draws at each depth) -- no need to
touch the 500MB per-draw file for this part.

Task 4 -- closure-orientation-change conditioning
----------------------------------------------------
``survival_samples.csv.gz`` (the per-draw file survival_curves.csv is built
from) DOES carry a "did this consistent draw leave the closure orientations
unchanged" signal: ``closure_inert`` (1 iff ``n_closure_orientations_changed
== 0``), logged only for accepted/consistent draws
(``corruption_accepted == 1``; see survival_p6.py's own docstring). This
script streams that 32.8M-row file once (via a single ``awk`` pass, far
faster than a Python csv reader for a flat groupby-count) to compute, per
(instance_id, arm, d): the count of consistent draws, how many left closure
orientations unchanged (inert) vs changed, and how many of the *changed* ones
survived. From that it reports p0 = P(inert | consistent) per arm/intensity,
survival conditional on a changed state, and tau_b(r_val, AUC of
changed-only survival) using the identical AUC_frac construction
(``survival.auc_frac``, same FRAC_GRID) applied to the changed-only curve.

Usage::

    PYTHONPATH=src .venv/bin/python experiments/paired_synth_v2.py
"""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from bkrobust.robustness.paired_resample import (  # noqa: E402
    _clean_matched,
    _kendall_tau_b,
    _to_float_or_none,
    paired_comparison,
)
from bkrobust.robustness.survival import FRAC_GRID  # noqa: E402

IN_INSTANCES = ROOT / "results" / "axis_robustness_p6" / "survival_instances.csv"
IN_CURVES = ROOT / "results" / "axis_robustness_p6" / "survival_curves.csv"
IN_SAMPLES_GZ = ROOT / "results" / "axis_robustness_p6" / "survival_samples.csv.gz"
COMMITTED_PAIRED = ROOT / "results" / "axis_robustness_p6_paired" / "paired_diff_synth.csv"

OUT_DIR = ROOT / "results" / "axis_robustness_p6_paired_v2"

ENDPOINT = {"flip": "AUC_frac", "tiered": "AUC_frac"}
BASELINES_OLD = ["shd_truth", "n_k", "k_g0", "separation", "phi_1", "r_claim"]
BASELINES_NEW = ["shd_truth", "n_k", "k_g0", "separation", "neg_phi_1", "r_claim"]
UNCHANGED_BASELINES = ["shd_truth", "n_k", "k_g0", "separation", "r_claim"]
ALL_MARGINAL_PREDICTORS = ["r_val", "shd_truth", "n_k", "k_g0", "separation", "neg_phi_1", "r_claim"]

N_BOOT = 10000
SEED = 0

PAIRED_COLUMNS = [
    "stratum", "arm", "predictor_a", "predictor_b", "endpoint", "status",
    "stratum_n", "n", "n_excluded", "n_clusters",
    "tau_a_point", "tau_b_point", "delta_point", "ci_lo_2p5", "ci_hi_97p5",
    "frac_delta_gt_0", "n_boot", "seed", "degeneracy_note",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_instances() -> list[dict[str, str]]:
    with IN_INSTANCES.open() as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        phi1 = _to_float_or_none(r.get("phi_1"))
        r["neg_phi_1"] = "" if phi1 is None else repr(-phi1)
    return rows


def stratum_of(row: dict[str, str]) -> str:
    if row["arm"] == "flip":
        return f"flip_cov{row['coverage']}_bw{row['base_wrongness']}"
    return f"tiered_nt{row['n_tiers']}"


def degeneracy_note(srows: list[dict[str, str]], baseline: str) -> str:
    """Identical to paired_synth_p6.degeneracy_note, extended with neg_phi_1
    (which inherits phi_1's own lack of a special-cased note: the committed
    script never special-cased phi_1 either -- both are always defined here,
    see the phi_1-blank-count check in this script's manifest)."""
    def col(name: str) -> list[str]:
        return [r[name] for r in srows if r[name] not in ("", "None")]

    if baseline == "k_g0":
        if col("k_g0") and col("k_g0") == col("n_k"):
            return "k_g0 == n_k on every row of this stratum; a relabelling, not an independent baseline."
        return "k_g0 is definitionally shd_cpdag, but differs from n_k here."
    if baseline == "separation":
        pairs = [(r["r_val"], r["separation"]) for r in srows
                 if r["r_val"] not in ("", "None") and r["separation"] not in ("", "None")]
        if pairs and all(a == b for a, b in pairs):
            return ("r_val == separation exactly on every row: the generator draws instances at "
                    "target separations, so a zero delta is an identity, not a null result.")
        return ""
    if baseline == "r_claim":
        pairs = [(r["r_val"], r["r_claim"]) for r in srows
                 if r["r_val"] not in ("", "None") and r["r_claim"] not in ("", "None")]
        n_cens = sum(1 for r in srows if r.get("r_claim_censored") == "1")
        note = f"{n_cens} of {len(srows)} rows right-censored and dropped, not imputed."
        if pairs and all(a == b for a, b in pairs):
            note = "r_val == r_claim on every uncensored row; a zero delta is an identity. " + note
        return note
    return ""


def run_paired(rows: list[dict[str, Any]], baselines: list[str],
                filter_fn=None) -> list[dict[str, Any]]:
    strata: dict[str, list[dict[str, Any]]] = {}
    use_rows = [r for r in rows if (filter_fn is None or filter_fn(r))]
    for r in use_rows:
        strata.setdefault(stratum_of(r), []).append(r)
    out_rows = []
    for stratum in sorted(strata):
        srows = strata[stratum]
        arm = srows[0]["arm"]
        endpoint = ENDPOINT[arm]
        for b in baselines:
            res = paired_comparison(srows, "r_val", b, endpoint, "instance_id",
                                     n_boot=N_BOOT, seed=SEED)
            out_rows.append({
                "stratum": stratum, "arm": arm, "predictor_a": "r_val", "predictor_b": b,
                "endpoint": endpoint, "status": res.get("status", "ok"),
                "stratum_n": len(srows), "n": res.get("n"), "n_excluded": res.get("n_excluded"),
                "n_clusters": res.get("n_clusters"), "tau_a_point": res.get("tau_a_point"),
                "tau_b_point": res.get("tau_b_point"), "delta_point": res.get("delta_point"),
                "ci_lo_2p5": res.get("ci_lo_2p5"), "ci_hi_97p5": res.get("ci_hi_97p5"),
                "frac_delta_gt_0": res.get("frac_delta_gt_0"), "n_boot": N_BOOT, "seed": SEED,
                "degeneracy_note": degeneracy_note(srows, b),
            })
    return out_rows


def verify_unchanged_baselines(rerun_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Diffs the rerun's rows for UNCHANGED_BASELINES against the committed CSV."""
    with COMMITTED_PAIRED.open() as fh:
        committed = list(csv.DictReader(fh))
    committed_by_key = {(r["stratum"], r["predictor_b"]): r for r in committed}
    discrepancies = []
    numeric_cols = ["tau_a_point", "tau_b_point", "delta_point", "ci_lo_2p5",
                     "ci_hi_97p5", "frac_delta_gt_0"]
    for r in rerun_rows:
        if r["predictor_b"] not in UNCHANGED_BASELINES:
            continue
        key = (r["stratum"], r["predictor_b"])
        c = committed_by_key.get(key)
        if c is None:
            discrepancies.append({"stratum": key[0], "baseline": key[1], "issue": "missing_in_committed"})
            continue
        row_diff = {}
        for col in numeric_cols:
            rv, cv = r.get(col), c.get(col)
            rv_f = None if rv in (None, "") else float(rv)
            cv_f = None if cv in (None, "") else float(cv)
            if rv_f is None and cv_f is None:
                continue
            if rv_f is None or cv_f is None or abs(rv_f - cv_f) > 1e-9:
                row_diff[col] = {"rerun": rv_f, "committed": cv_f}
        if r["status"] != c["status"]:
            row_diff["status"] = {"rerun": r["status"], "committed": c["status"]}
        if row_diff:
            discrepancies.append({"stratum": key[0], "baseline": key[1], "diff": row_diff})
    return discrepancies


# ---------------------------------------------------------------------------
# Task 1: marginal tau_b per predictor per stratum, instance bootstrap CI
# ---------------------------------------------------------------------------

def marginal_tau_bootstrap(rows: list[dict[str, Any]], predictor: str, endpoint: str,
                            n_boot: int = N_BOOT, seed: int = SEED) -> dict[str, Any]:
    cleaned = []
    n_excluded = 0
    for r in rows:
        pv = _to_float_or_none(r.get(predictor))
        ev = _to_float_or_none(r.get(endpoint))
        if pv is None or ev is None:
            n_excluded += 1
            continue
        cleaned.append((pv, ev))
    out = {"predictor": predictor, "endpoint": endpoint, "stratum_n": len(rows),
           "n": len(cleaned), "n_excluded": n_excluded,
           "tau_b": None, "ci_lo_2p5": None, "ci_hi_97p5": None,
           "n_boot": n_boot, "seed": seed, "status": "ok"}
    if len(cleaned) < 2:
        out["status"] = "insufficient_n"
        return out
    x = np.array([c[0] for c in cleaned])
    y = np.array([c[1] for c in cleaned])
    tau = _kendall_tau_b(x, y)
    out["tau_b"] = tau
    if tau is None:
        out["status"] = "undefined (predictor or endpoint constant)"
        return out
    rng = np.random.default_rng(seed)
    n = len(cleaned)
    boots = np.full(n_boot, np.nan)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        t = _kendall_tau_b(x[idx], y[idx])
        if t is not None:
            boots[b] = t
    valid = boots[~np.isnan(boots)]
    if len(valid) >= 100:
        lo, hi = np.percentile(valid, [2.5, 97.5])
        out["ci_lo_2p5"], out["ci_hi_97p5"] = float(lo), float(hi)
    return out


def task1_marginals(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    strata: dict[str, list[dict[str, Any]]] = {}
    for r in rows:
        strata.setdefault(stratum_of(r), []).append(r)
    out = []
    for stratum in sorted(strata):
        srows = strata[stratum]
        arm = srows[0]["arm"]
        endpoint = ENDPOINT[arm]
        for pred in ALL_MARGINAL_PREDICTORS:
            res = marginal_tau_bootstrap(srows, pred, endpoint)
            res["stratum"] = stratum
            res["arm"] = arm
            out.append(res)
    return out


# ---------------------------------------------------------------------------
# Task 2: separation facts on the current sweep
# ---------------------------------------------------------------------------

def separation_key(row: dict[str, Any]) -> tuple:
    if row["arm"] == "flip":
        return ("flip", row["component_size"], row["coverage"], row["base_wrongness"])
    return ("tiered", row["component_size"], row["n_tiers"])


def separation_neq_facts(rows: list[dict[str, Any]]) -> dict[str, Any]:
    measured = [r for r in rows if r.get("separation_status") == "measured"]
    undefined = [r for r in rows if r.get("separation_status") != "measured"]

    def frac_neq(subset):
        n = len(subset)
        if n == 0:
            return {"n": 0, "n_neq": 0, "frac_neq": None}
        neq = sum(1 for r in subset if r["r_val"] != r["separation"])
        return {"n": n, "n_neq": neq, "frac_neq": neq / n}

    out = {
        "overall": frac_neq(measured),
        "n_undefined_separation": len(undefined),
        "by_arm": {}, "by_stratum": {},
    }
    for arm in ("flip", "tiered"):
        out["by_arm"][arm] = frac_neq([r for r in measured if r["arm"] == arm])
    strata: dict[str, list] = {}
    for r in measured:
        strata.setdefault(stratum_of(r), []).append(r)
    for s in sorted(strata):
        out["by_stratum"][s] = frac_neq(strata[s])
    return out


def stratified_tau_bootstrap(rows: list[dict[str, Any]], predictor: str, endpoint: str,
                              key_fn, n_boot: int = N_BOOT, seed: int = SEED) -> dict[str, Any]:
    """Weighted average of per-separation-stratum tau_b(predictor, endpoint),
    weighted by comparable-pair count n*(n-1)/2 per stratum, with an
    instance-cluster bootstrap CI that recomputes the same weighted average
    on each resample (drawing len(rows) instances with replacement from the
    full pooled set, then re-forming strata and re-weighting)."""
    def clean_by_stratum(rr):
        strata: dict[tuple, list[tuple[float, float]]] = {}
        n_excluded = 0
        for r in rr:
            pv = _to_float_or_none(r.get(predictor))
            ev = _to_float_or_none(r.get(endpoint))
            if pv is None or ev is None:
                n_excluded += 1
                continue
            strata.setdefault(key_fn(r), []).append((pv, ev))
        return strata, n_excluded

    def weighted_tau(strata: dict[tuple, list[tuple[float, float]]]) -> float | None:
        num, den = 0.0, 0.0
        for k, pairs in strata.items():
            if len(pairs) < 2:
                continue
            x = np.array([p[0] for p in pairs])
            y = np.array([p[1] for p in pairs])
            t = _kendall_tau_b(x, y)
            if t is None:
                continue
            n = len(pairs)
            w = n * (n - 1) / 2.0
            num += w * t
            den += w
        if den == 0:
            return None
        return num / den

    strata0, n_excluded = clean_by_stratum(rows)
    point = weighted_tau(strata0)
    n_strata_used = sum(1 for k, v in strata0.items() if len(v) >= 2)
    out = {
        "predictor": predictor, "endpoint": endpoint, "n_rows": len(rows),
        "n_excluded": n_excluded, "n_strata_total": len(strata0),
        "n_strata_used": n_strata_used, "tau_b_stratified": point,
        "ci_lo_2p5": None, "ci_hi_97p5": None, "n_boot": n_boot, "seed": seed,
        "status": "ok" if point is not None else "undefined",
    }
    if point is None:
        return out

    # Instance-level cluster bootstrap: resample instances (each its own
    # cluster) from the full matched pool, re-form separation strata, and
    # recompute the same weighted average.
    cleaned_rows = [r for r in rows if _to_float_or_none(r.get(predictor)) is not None
                     and _to_float_or_none(r.get(endpoint)) is not None]
    n = len(cleaned_rows)
    rng = np.random.default_rng(seed)
    boots = np.full(n_boot, np.nan)
    keys = [key_fn(r) for r in cleaned_rows]  # plain list of tuples; a numpy
    # object array of same-length tuples silently degrades to a 2D array,
    # whose rows are unhashable -- keep this a Python list of tuple keys.
    xs = np.array([_to_float_or_none(r.get(predictor)) for r in cleaned_rows])
    ys = np.array([_to_float_or_none(r.get(endpoint)) for r in cleaned_rows])
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        rx, ry = xs[idx], ys[idx]
        by_k: dict[Any, list[int]] = defaultdict(list)
        for i, orig_i in enumerate(idx):
            by_k[keys[orig_i]].append(i)
        num, den = 0.0, 0.0
        for k, ilist in by_k.items():
            if len(ilist) < 2:
                continue
            xi = rx[ilist]
            yi = ry[ilist]
            t = _kendall_tau_b(xi, yi)
            if t is None:
                continue
            w = len(ilist) * (len(ilist) - 1) / 2.0
            num += w * t
            den += w
        if den > 0:
            boots[b] = num / den
    valid = boots[~np.isnan(boots)]
    if len(valid) >= 100:
        lo, hi = np.percentile(valid, [2.5, 97.5])
        out["ci_lo_2p5"], out["ci_hi_97p5"] = float(lo), float(hi)
    out["n_boot_valid"] = int(len(valid))
    return out


def task2_separation(rows: list[dict[str, Any]]) -> dict[str, Any]:
    neq = separation_neq_facts(rows)
    flip_rows = [r for r in rows if r["arm"] == "flip"]
    tiered_rows = [r for r in rows if r["arm"] == "tiered"]

    within_strata = {}
    other_preds = ["n_k", "k_g0", "shd_truth", "neg_phi_1"]
    for arm_name, arm_rows in (("tiered", tiered_rows), ("flip", flip_rows)):
        within_strata[arm_name] = {"r_val": stratified_tau_bootstrap(
            arm_rows, "r_val", "AUC_frac", separation_key)}
        for p in other_preds:
            within_strata[arm_name][p] = stratified_tau_bootstrap(
                arm_rows, p, "AUC_frac", separation_key)
    return {"neq_facts": neq, "within_separation_strata_tau": within_strata}


# ---------------------------------------------------------------------------
# Task 3: AUC support
# ---------------------------------------------------------------------------

def task3_support(rows: list[dict[str, Any]], curve_rows: list[dict[str, str]]) -> dict[str, Any]:
    def n_defined(r):
        ev = int(r["n_depths_evaluated"])
        contra = int(r["n_depths_all_contradictory"])
        return ev - contra

    out: dict[str, Any] = {"by_arm": {}, "by_stratum": {}}
    for arm in ("flip", "tiered"):
        arm_rows = [r for r in rows if r["arm"] == arm]
        n_full = sum(1 for r in arm_rows if int(r["n_depths_all_contradictory"]) == 0)
        n_partial = len(arm_rows) - n_full
        defined_counts = [n_defined(r) for r in arm_rows]
        out["by_arm"][arm] = {
            "n": len(arm_rows), "n_full_support": n_full,
            "n_some_undefined_gridpoint": n_partial,
            "frac_some_undefined": n_partial / len(arm_rows) if arm_rows else None,
            "defined_gridpoints_mean": float(np.mean(defined_counts)) if defined_counts else None,
            "defined_gridpoints_median": float(np.median(defined_counts)) if defined_counts else None,
            "defined_gridpoints_min": int(np.min(defined_counts)) if defined_counts else None,
            "defined_gridpoints_max": int(np.max(defined_counts)) if defined_counts else None,
        }
    strata: dict[str, list] = {}
    for r in rows:
        strata.setdefault(stratum_of(r), []).append(r)
    for s in sorted(strata):
        srows = strata[s]
        n_full = sum(1 for r in srows if int(r["n_depths_all_contradictory"]) == 0)
        out["by_stratum"][s] = {
            "n": len(srows), "n_full_support": n_full,
            "frac_full_support": n_full / len(srows) if srows else None,
        }

    # Full-support subset: rerun tau_b(r_val, AUC_frac) and the paired
    # differences vs baselines, restricted to n_depths_all_contradictory == 0.
    full_rows = [r for r in rows if int(r["n_depths_all_contradictory"]) == 0]
    out["n_full_support_total"] = len(full_rows)
    out["full_support_paired"] = run_paired(full_rows, BASELINES_NEW)
    out["full_support_marginal_r_val"] = []
    strata_full: dict[str, list] = {}
    for r in full_rows:
        strata_full.setdefault(stratum_of(r), []).append(r)
    for s in sorted(strata_full):
        srows = strata_full[s]
        res = marginal_tau_bootstrap(srows, "r_val", ENDPOINT[srows[0]["arm"]])
        res["stratum"] = s
        out["full_support_marginal_r_val"].append(res)

    # Accepted (consistent) draws per intensity (median/min across
    # instances) per arm, straight from survival_curves.csv's n_eval.
    by_arm_d: dict[tuple, list[int]] = defaultdict(list)
    for r in curve_rows:
        by_arm_d[(r["arm"], int(r["d"]))].append(int(r["n_eval"]))
    accepted_per_intensity = []
    for (arm, d), vals in sorted(by_arm_d.items()):
        accepted_per_intensity.append({
            "arm": arm, "d": d, "n_instances": len(vals),
            "median_n_eval": float(np.median(vals)), "min_n_eval": int(np.min(vals)),
            "max_n_eval": int(np.max(vals)),
        })
    out["accepted_draws_per_intensity"] = accepted_per_intensity
    return out


# ---------------------------------------------------------------------------
# Task 4: closure-orientation-change conditioning
# ---------------------------------------------------------------------------

def stream_changed_state_agg() -> Path:
    """One awk pass over survival_samples.csv.gz aggregating, per
    (instance_id, arm, d): n_consistent, n_inert, n_changed, n_survived_changed.
    Column positions are the file's fixed header order (verified against the
    committed header before running); this avoids a ~500MB-file Python csv
    parse, which would be far slower than a single streaming awk pass.
    """
    out_path = OUT_DIR / "_changed_state_agg.tsv"
    with IN_SAMPLES_GZ.open("rb"):
        pass  # existence check
    # Verify header order matches what the fixed column indices below assume.
    header = subprocess.run(
        ["bash", "-c", f"gzip -dc {IN_SAMPLES_GZ} | head -1"],
        capture_output=True, text=True, check=True,
    ).stdout.strip().split(",")
    expected = ["instance_id", "arm", "gate", "coverage", "base_wrongness", "n_tiers", "d",
                "rep", "flip_rate", "corruption_rate", "seed", "status", "survived",
                "symdiff_proxy_not_distance", "n_claims_attempted", "n_claims_reversed",
                "corruption_rejected", "corruption_accepted", "n_closure_orientations_changed",
                "closure_inert"]
    if header != expected:
        raise RuntimeError(f"survival_samples.csv.gz header changed: {header} != {expected}")

    awk_prog = (
        "NR==1{next}"
        "{"
        "  ca=$18;"
        "  if (ca==1) {"
        "    key=$1 SUBSEP $2 SUBSEP $7;"
        "    n_cons[key]++;"
        "    if ($20==1) { n_inert[key]++ }"
        "    else {"
        "      n_changed[key]++;"
        "      if ($13==\"True\") n_surv_changed[key]++"
        "    }"
        "  }"
        "}"
        "END{ for (k in n_cons) print k, n_cons[k]+0, n_inert[k]+0, n_changed[k]+0, n_surv_changed[k]+0 }"
    )
    cmd = f"gzip -dc {IN_SAMPLES_GZ} | awk -F',' '{awk_prog}' OFS='\\t' > {out_path}"
    t0 = time.perf_counter()
    subprocess.run(["bash", "-c", cmd], check=True)
    return out_path, time.perf_counter() - t0


def auc_frac_from_S(S_by_d: dict[int, float | None], n_k: int) -> float | None:
    defined = sorted(d for d, s in S_by_d.items() if s is not None)
    if not defined or n_k <= 0:
        return None
    vals = []
    for frac in FRAC_GRID:
        target = frac * n_k
        nearest = min(defined, key=lambda d: (abs(d - target), d))
        vals.append(S_by_d[nearest])
    return sum(vals) / len(vals)


def task4_changed_state(rows: list[dict[str, Any]]) -> dict[str, Any]:
    if not IN_SAMPLES_GZ.is_file():
        return {"status": "no_per_draw_file",
                "note": "survival_samples.csv.gz not found; skipped, no corruption sweep rerun."}
    agg_path, elapsed = stream_changed_state_agg()
    # key -> (n_cons, n_inert, n_changed, n_surv_changed)
    per_instance_d: dict[tuple, tuple[int, int, int, int]] = {}
    arm_d_totals: dict[tuple, list[int]] = defaultdict(lambda: [0, 0, 0, 0])
    with agg_path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 5:
                continue
            inst, arm, d, n_cons, n_inert, n_changed, n_surv = None, None, None, None, None, None, None
            # key was joined with SUBSEP (default \034) inside awk's key var, not \t
            key_field = parts[0]
            sub = "\x1c"
            kparts = key_field.split(sub)
            if len(kparts) != 3:
                continue
            inst, arm, d = kparts[0], kparts[1], int(kparts[2])
            n_cons, n_inert, n_changed, n_surv = (int(parts[1]), int(parts[2]),
                                                    int(parts[3]), int(parts[4]))
            per_instance_d[(inst, arm, d)] = (n_cons, n_inert, n_changed, n_surv)
            tot = arm_d_totals[(arm, d)]
            tot[0] += n_cons; tot[1] += n_inert; tot[2] += n_changed; tot[3] += n_surv

    p0_by_arm_intensity = []
    surv_given_changed_by_arm_intensity = []
    for (arm, d), (n_cons, n_inert, n_changed, n_surv) in sorted(arm_d_totals.items()):
        p0_by_arm_intensity.append({
            "arm": arm, "d": d, "n_consistent": n_cons,
            "p0_inert_given_consistent": (n_inert / n_cons) if n_cons else None,
        })
        surv_given_changed_by_arm_intensity.append({
            "arm": arm, "d": d, "n_changed": n_changed,
            "survival_given_changed": (n_surv / n_changed) if n_changed else None,
        })

    # Build the changed-only survival curve per instance and its AUC_frac.
    n_k_by_instance = {r["instance_id"]: int(r["n_k"]) for r in rows}
    arm_by_instance = {r["instance_id"]: r["arm"] for r in rows}
    r_val_by_instance = {r["instance_id"]: _to_float_or_none(r["r_val"]) for r in rows}
    by_instance: dict[str, dict[int, tuple[int, int]]] = defaultdict(dict)
    for (inst, arm, d), (n_cons, n_inert, n_changed, n_surv) in per_instance_d.items():
        by_instance[inst][d] = (n_changed, n_surv)

    changed_auc_rows = []
    for inst, dmap in by_instance.items():
        n_k = n_k_by_instance.get(inst)
        if not n_k:
            continue
        S_by_d = {}
        for d, (n_changed, n_surv) in dmap.items():
            S_by_d[d] = (n_surv / n_changed) if n_changed > 0 else None
        auc = auc_frac_from_S(S_by_d, n_k)
        changed_auc_rows.append({
            "instance_id": inst, "arm": arm_by_instance.get(inst),
            "r_val": r_val_by_instance.get(inst), "AUC_frac_changed_only": auc,
        })

    # tau_b(r_val, AUC_frac_changed_only), per arm and overall.
    tau_changed = {}
    for arm in ("flip", "tiered", "all"):
        sub = [r for r in changed_auc_rows
               if (arm == "all" or r["arm"] == arm) and r["r_val"] is not None
               and r["AUC_frac_changed_only"] is not None]
        if len(sub) < 2:
            tau_changed[arm] = {"n": len(sub), "tau_b": None}
            continue
        x = np.array([r["r_val"] for r in sub])
        y = np.array([r["AUC_frac_changed_only"] for r in sub])
        tau_changed[arm] = {"n": len(sub), "tau_b": _kendall_tau_b(x, y)}

    return {
        "status": "ok",
        "note": "closure_inert (1 iff n_closure_orientations_changed==0) found in "
                "survival_samples.csv.gz, logged only for corruption_accepted==1 rows.",
        "awk_pass_seconds": elapsed,
        "p0_inert_given_consistent_by_arm_intensity": p0_by_arm_intensity,
        "survival_given_changed_by_arm_intensity": surv_given_changed_by_arm_intensity,
        "tau_b_rval_vs_AUC_changed_only": tau_changed,
        "n_instances_with_changed_auc": len(changed_auc_rows),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in fieldnames})


def main() -> None:
    t0 = time.perf_counter()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows = load_instances()
    print(f"[v2] loaded {len(rows)} instances", flush=True)
    n_phi1_blank = sum(1 for r in rows if r["phi_1"] == "")
    print(f"[v2] phi_1 blank count: {n_phi1_blank}", flush=True)

    with IN_CURVES.open() as fh:
        curve_rows = list(csv.DictReader(fh))

    # ---- Task 1 ----
    print("[v2] task1: paired rerun (new baselines incl. neg_phi_1)", flush=True)
    paired_v2 = run_paired(rows, BASELINES_NEW)
    write_csv(OUT_DIR / "paired_diff_synth_v2.csv", PAIRED_COLUMNS, paired_v2)

    print("[v2] task1: verifying unchanged baselines reproduce committed run", flush=True)
    discrepancies = verify_unchanged_baselines(paired_v2)
    (OUT_DIR / "reproduction_check.json").write_text(json.dumps(discrepancies, indent=1))
    print(f"[v2] task1: {len(discrepancies)} discrepancies vs committed run", flush=True)

    print("[v2] task1: marginal tau_b per predictor per stratum", flush=True)
    marginals = task1_marginals(rows)
    write_csv(OUT_DIR / "marginal_tau_by_stratum.csv",
              ["stratum", "arm", "predictor", "endpoint", "stratum_n", "n", "n_excluded",
               "tau_b", "ci_lo_2p5", "ci_hi_97p5", "n_boot", "seed", "status"],
              marginals)

    # r_claim n used, per stratum (for reporting censoring handling).
    r_claim_n = {}
    strata_map: dict[str, list] = {}
    for r in rows:
        strata_map.setdefault(stratum_of(r), []).append(r)
    for s, srows in strata_map.items():
        n_blank = sum(1 for r in srows if r["r_claim"] == "")
        r_claim_n[s] = {"stratum_n": len(srows), "n_r_claim_defined": len(srows) - n_blank,
                         "n_r_claim_censored_blank": n_blank}

    # ---- Task 2 ----
    print("[v2] task2: separation facts", flush=True)
    task2 = task2_separation(rows)
    (OUT_DIR / "separation_facts.json").write_text(json.dumps(task2, indent=1, default=str))

    # ---- Task 3 ----
    print("[v2] task3: AUC support", flush=True)
    task3 = task3_support(rows, curve_rows)
    (OUT_DIR / "auc_support.json").write_text(json.dumps(task3, indent=1, default=str))
    write_csv(OUT_DIR / "full_support_paired.csv", PAIRED_COLUMNS, task3["full_support_paired"])

    # ---- Task 4 ----
    print("[v2] task4: closure-orientation-change conditioning (streaming 500MB file)", flush=True)
    task4 = task4_changed_state(rows)
    (OUT_DIR / "changed_state_analysis.json").write_text(json.dumps(task4, indent=1, default=str))

    # ---- Manifest ----
    manifest = {
        "inputs": {
            "survival_instances_csv": {
                "path": str(IN_INSTANCES.relative_to(ROOT)),
                "sha256": sha256_file(IN_INSTANCES), "n_rows": len(rows),
            },
            "survival_curves_csv": {
                "path": str(IN_CURVES.relative_to(ROOT)),
                "sha256": sha256_file(IN_CURVES), "n_rows": len(curve_rows),
            },
            "committed_paired_diff_synth_csv": {
                "path": str(COMMITTED_PAIRED.relative_to(ROOT)),
                "sha256": sha256_file(COMMITTED_PAIRED),
            },
            "survival_samples_csv_gz": (
                {"path": str(IN_SAMPLES_GZ.relative_to(ROOT)), "sha256": sha256_file(IN_SAMPLES_GZ)}
                if IN_SAMPLES_GZ.is_file() else None
            ),
        },
        "n_boot": N_BOOT, "seed": SEED,
        "baselines_new": BASELINES_NEW, "baselines_old": BASELINES_OLD,
        "n_reproduction_discrepancies": len(discrepancies),
        "r_claim_n_used_by_stratum": r_claim_n,
        "elapsed_seconds": time.perf_counter() - t0,
    }
    (OUT_DIR / "manifest.json").write_text(json.dumps(manifest, indent=1, default=str))
    print(f"[v2] done in {manifest['elapsed_seconds']:.1f}s", flush=True)


if __name__ == "__main__":
    main()
