"""Corrected-cohort Table 1 and Table 11 (App. tab:balanced-full):
old cohort (results/axis_robustness_llm/analysis_units.csv) + the 94
real-naming, GAC-valid "empty Z*" units that commit_z_star's ``if not z:``
bug excluded (results/axis_robustness_llm_v2/empty_valid_units.csv, produced
by experiments/llm_empty_valid_survival_v2.py).

Step 1 (validate the code path): reproduces the paper's OLD-cohort Table 1
numbers exactly, using the identical machinery as
``experiments/llm_paired_v2.py`` (``real_analyse.tau_for_stratum``,
``paired_resample.paired_comparison``, and the same within-state stratified
tau_b + network-cluster bootstrap helpers), same n_boot=10000, seed=0.

Step 2: builds the corrected 2,516-unit (+94) cohort and recomputes every
Table 1 cell (pooled / balanced-8of8 / within-state, for r_val, k_g0, |K|,
SHD, separation, k_accuracy), the pooled paired leads of r_val over the other
predictors, the within-state matched-separation comparison, and Table 11's
">= k/8 scorable" pooled thresholds.

Separation for the 94 empty-Z* units: NOT undefined. Confirmed by reading
``bkrobust.benchmarks.measure.separation`` (used to populate the
`separation`/`separation_status` columns joined into every unit row) --
it measures graph distance, inside X's undirected CPDAG component, to the
nearest member of Z, where Z here is the row's *committed* separation input
(joined upstream from the per-(network,x,y) frame keyed off the analyst's
knowledge K/G0, not off the LLM run's committed Z*). ``commit_z_star``'s
``if not z:`` bug only ever affected which branch of survival-endpoint
scoring an LLM unit went through -- it never touched the separately-computed
`separation` column, which `empty_valid_survival_v2.py` carries through
unchanged from `build_state`. The empty_valid_units.csv data confirms this
directly: all 94 real-naming rows have `separation_status == "measured"`
with values in {1, 2, 6}. So separation IS defined for these units and is
included in both cohorts on the same footing as every other predictor.

Usage::

    PYTHONPATH=src .venv/bin/python experiments/table1_corrected_cohort.py
"""

from __future__ import annotations

import csv
import json
import math
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(REPO_ROOT / "experiments"))

from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.robustness.paired_resample import (  # noqa: E402
    DEFAULT_N_BOOT,
    DEFAULT_SEED,
    paired_comparison,
)
from bkrobust.robustness.real_analyse import tau_for_stratum  # noqa: E402

import llm_paired_v2 as base  # noqa: E402  (reuse the validated helpers)

OLD_UNITS_CSV = REPO_ROOT / "results" / "axis_robustness_llm" / "analysis_units.csv"
NEW_UNITS_CSV = REPO_ROOT / "results" / "axis_robustness_llm_v2" / "empty_valid_units.csv"
OUT_DIR = REPO_ROOT / "results" / "axis_robustness_llm_v2"

ENDPOINT = "AUC_frac"
N_BOOT = DEFAULT_N_BOOT  # 10_000
SEED = DEFAULT_SEED  # 0

PREDICTORS = ["radius", "k_g0", "n_k", "shd_truth", "separation", "k_accuracy"]
PREDICTOR_LABEL = {
    "radius": "r_val", "k_g0": "|K_G0|", "n_k": "|K|", "shd_truth": "SHD",
    "separation": "separation s", "k_accuracy": "claim accuracy k_acc",
}

PAPER_POOLED_OLD = {
    "radius": (0.319, 0.122, 0.485),
    "k_g0": (0.214, 0.079, 0.339),
    "n_k": (0.185, 0.071, 0.327),
    "shd_truth": (0.054, -0.077, 0.200),
}
PAPER_BALANCED_OLD = {
    "radius": (0.499, 0.240, 0.624),
    "k_g0": (0.380, -0.074, 0.562),
    "n_k": (0.333, -0.125, 0.484),
    "shd_truth": (0.285, -0.167, 0.456),
}
PAPER_WITHIN_OLD = (0.419, 0.053, 0.679)
PAPER_PAIRED_OLD = {
    "shd_truth": (0.27, 0.04, 0.46),
    "k_g0": (0.11, -0.09, 0.30),
    "n_k": (0.13, -0.08, 0.32),
}
PAPER_SEP_SUBSET_OLD = {
    "radius": (0.68, 0.46, 0.77),
    "separation": (0.49, -0.14, 0.60),
    "paired_lead": (0.19, 0.07, 0.90),
}

PAPER_CORRECTED = {
    "pooled_radius": (0.32, 0.13, 0.49),
    "pooled_paired": {"k_g0": (0.14, -0.06, 0.32), "n_k": (0.16, -0.07, 0.35),
                       "shd_truth": (0.30, 0.06, 0.49)},
    "within_radius": (0.43, 0.07, 0.69),
}


# ---------------------------------------------------------------------------
# loading / cohort assembly
# ---------------------------------------------------------------------------

NUMERIC_COLS = [
    "radius", "shd_truth", "n_k", "k_g0", "AUC_frac", "AUC_frac_usable",
    "separation", "k_accuracy",
]


def _to_float(v: Any) -> float | None:
    if v is None:
        return None
    s = str(v).strip()
    if s == "":
        return None
    try:
        f = float(s)
    except ValueError:
        return None
    if math.isnan(f):
        return None
    return f


def load_units(path: Path) -> list[dict[str, Any]]:
    rows = []
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            row = dict(r)
            for c in NUMERIC_COLS:
                row[c] = _to_float(row.get(c))
            rows.append(row)
    return rows


def build_cohorts() -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    old_units = load_units(OLD_UNITS_CSV)
    new_units_all = load_units(NEW_UNITS_CSV)
    # Only the 94 real-naming units belong in the paper's real-name cohort;
    # the 28 scrambled-naming empty_valid rows belong to the scrambled control
    # (out of scope for Table 1 / Table 11, which are real-naming only).
    new_units_real = [u for u in new_units_all if u.get("naming") == "real"]
    assert len(new_units_real) == 94, f"expected 94 real-naming empty_valid units, got {len(new_units_real)}"

    # `assumes` wording normalisation: the committed panel's rows carry the
    # older epistemic-status label "Conjecture 2 (hence Anti-Exchange Case B,
    # verified not proved)"; empty_valid_units.csv (built against current
    # bkrobust.hybrid.HybridResult) carries the current label "Conjecture 2
    # (proved: Anti-Exchange Case B, THEOREMS.md section 4)". Same theorem
    # citation, same dispatch legs, same oracle, same exact=True -- only a
    # documentation-string vintage mismatch. This exact normalisation is
    # already established and justified by
    # experiments/llm_empty_valid_sensitivity_v2.py (_COMMITTED_ASSUMES /
    # _CURRENT_ASSUMES), which pools these same 94 rows the same way so that
    # real_analyse.tau_for_stratum's within-stratum `assumes` uniformity
    # guard does not fire spuriously. Reused verbatim here.
    _COMMITTED_ASSUMES = "Conjecture 2 (hence Anti-Exchange Case B, verified not proved)"
    _CURRENT_ASSUMES = "Conjecture 2 (proved: Anti-Exchange Case B, THEOREMS.md section 4)"
    n_normalised = 0
    for u in new_units_real:
        if u.get("assumes") == _CURRENT_ASSUMES:
            u["assumes"] = _COMMITTED_ASSUMES
            n_normalised += 1
    print(f"normalised 'assumes' wording on {n_normalised}/{len(new_units_real)} new units "
          f"(vintage label only; same theorem citation)")

    combined = old_units + new_units_real
    return old_units, new_units_real, combined


# ---------------------------------------------------------------------------
# Table 1 (pooled / balanced-8of8 / within-state) for one cohort
# ---------------------------------------------------------------------------


def table1_for_cohort(units: list[dict[str, Any]]) -> dict[str, Any]:
    pooled = base.real_naming(units)
    ok_by_triple = base.ok_conditions_by_triple(units)
    balanced = base.subset_at_least(units, ok_by_triple, 8)

    out: dict[str, Any] = {"pooled": {}, "balanced_8of8": {}, "within_state": {}}
    for label, rows in (("pooled", pooled), ("balanced_8of8", balanced)):
        for pred in PREDICTORS:
            out[label][pred] = tau_for_stratum(rows, pred, ENDPOINT, n_boot=N_BOOT, seed=SEED)

    # within-state
    by_state: dict[tuple, list[dict[str, Any]]] = defaultdict(list)
    for u in pooled:
        rv, ev = u.get("radius"), u.get(ENDPOINT)
        if rv is None or ev is None or rv == UNREACHED:
            continue
        by_state[(u["network"], u["condition"])].append(u)
    varying_states = {
        k: rows for k, rows in by_state.items() if len({r["radius"] for r in rows}) > 1
    }
    res_radius = base._stratified_tau_b(varying_states, "radius", ENDPOINT)
    lo, hi, nb = base._within_state_bootstrap(varying_states, "radius", ENDPOINT, n_boot=N_BOOT, seed=SEED)
    res_radius["ci_lo_2p5"], res_radius["ci_hi_97p5"], res_radius["n_boot_nan"] = lo, hi, nb
    res_radius["n_units"] = sum(len(r) for r in varying_states.values())
    res_radius["n_networks"] = len({k[0] for k in varying_states})
    res_radius["n_states"] = len(varying_states)
    out["within_state"]["radius"] = res_radius

    # within-state k_accuracy (claim accuracy varies with condition/naming;
    # compute over the same states as radius for a like-for-like column)
    for pred in ("k_g0", "n_k", "shd_truth", "k_accuracy"):
        sub_states: dict[tuple, list[dict[str, Any]]] = {}
        for k, rows in varying_states.items():
            s = [r for r in rows if r.get(pred) is not None]
            if len(s) >= 2 and len({r[pred] for r in s}) > 1:
                sub_states[k] = s
        r = base._stratified_tau_b(sub_states, pred, ENDPOINT)
        lo2, hi2, nb2 = base._within_state_bootstrap(sub_states, pred, ENDPOINT, n_boot=N_BOOT, seed=SEED)
        r["ci_lo_2p5"], r["ci_hi_97p5"], r["n_boot_nan"] = lo2, hi2, nb2
        r["n_units"] = sum(len(v) for v in sub_states.values())
        r["n_networks"] = len({k[0] for k in sub_states})
        r["n_states"] = len(sub_states)
        out["within_state"][pred] = r

    # within-state separation, matched to states where it's defined (same
    # rule the paper used: subset rows within each radius-varying state)
    sep_states: dict[tuple, list[dict[str, Any]]] = {}
    for k, rows in varying_states.items():
        sub = [r for r in rows if r.get("separation") is not None]
        if len(sub) >= 2:
            sep_states[k] = sub
    res_sep = base._stratified_tau_b(sep_states, "separation", ENDPOINT)
    lo_s, hi_s, nb_s = base._within_state_bootstrap(sep_states, "separation", ENDPOINT, n_boot=N_BOOT, seed=SEED)
    res_sep["ci_lo_2p5"], res_sep["ci_hi_97p5"], res_sep["n_boot_nan"] = lo_s, hi_s, nb_s
    res_sep["n_units"] = sum(len(r) for r in sep_states.values())
    res_sep["n_networks"] = len({k[0] for k in sep_states})
    res_sep["n_states"] = len(sep_states)
    out["within_state"]["separation"] = res_sep

    # matched radius-on-separation-defined-subset (paper's App. matched
    # comparison: restrict to the SAME rows separation is defined on --
    # sep_states already holds exactly those rows, grouped the same way --
    # not merely the same state keys with all their (possibly
    # separation-undefined) rows reinstated.
    matched_radius_states = sep_states
    res_radius_matched = base._stratified_tau_b(matched_radius_states, "radius", ENDPOINT)
    lo_m, hi_m, nb_m = base._within_state_bootstrap(matched_radius_states, "radius", ENDPOINT, n_boot=N_BOOT, seed=SEED)
    res_radius_matched["ci_lo_2p5"], res_radius_matched["ci_hi_97p5"], res_radius_matched["n_boot_nan"] = lo_m, hi_m, nb_m
    res_radius_matched["n_units"] = sum(len(r) for r in matched_radius_states.values())
    res_radius_matched["n_networks"] = len({k[0] for k in matched_radius_states})
    res_radius_matched["n_states"] = len(matched_radius_states)
    out["within_state_matched_separation"] = {
        "radius": res_radius_matched, "separation": res_sep,
    }
    # paired within-state bootstrap radius - separation on the matched states
    paired_matched = base._paired_within_state_bootstrap(
        matched_radius_states, "radius", "separation", ENDPOINT, n_boot=N_BOOT, seed=SEED
    ) if hasattr(base, "_paired_within_state_bootstrap") else None
    out["within_state_matched_separation"]["paired_lead_available"] = paired_matched is not None
    if paired_matched is not None:
        out["within_state_matched_separation"]["paired_lead"] = paired_matched

    # pooled separation-defined subset counts (paper reports n/networks where s defined)
    pooled_sep_defined = [u for u in pooled if u.get("separation") is not None]
    out["pooled_separation_defined_n"] = len(pooled_sep_defined)
    out["pooled_separation_defined_n_networks"] = len({u["network"] for u in pooled_sep_defined})
    balanced_sep_defined = [u for u in balanced if u.get("separation") is not None]
    out["balanced_separation_defined_n"] = len(balanced_sep_defined)
    out["balanced_separation_defined_n_networks"] = len({u["network"] for u in balanced_sep_defined})

    out["pooled_n"] = len(pooled)
    out["pooled_n_networks"] = len({u["network"] for u in pooled})
    out["balanced_n"] = len(balanced)
    out["balanced_n_networks"] = len({u["network"] for u in balanced})

    return out


def paired_leads_pooled(units: list[dict[str, Any]]) -> dict[str, Any]:
    pooled = base.real_naming(units)
    out = {}
    for other in ("shd_truth", "k_g0", "n_k"):
        out[other] = paired_comparison(pooled, "radius", other, ENDPOINT, "network",
                                        n_boot=N_BOOT, seed=SEED)
    return out


def table11(units: list[dict[str, Any]]) -> list[dict[str, Any]]:
    ok_by_triple = base.ok_conditions_by_triple(units)
    rows = []
    for k in (8, 7, 6, 5, 4, 3, 2, 1):
        subset = base.subset_at_least(units, ok_by_triple, k)
        n_triples = len({(u["network"], u["x"], u["y"]) for u in subset})
        n_networks = len({u["network"] for u in subset})
        row: dict[str, Any] = {"threshold_k_of_8": k, "n_units": len(subset),
                                "n_triples": n_triples, "n_networks": n_networks}
        for pred in PREDICTORS:
            row[pred] = tau_for_stratum(subset, pred, ENDPOINT, n_boot=N_BOOT, seed=SEED)
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Step 1: reproduce OLD cohort exactly
# ---------------------------------------------------------------------------


def reproduce_old(old_units: list[dict[str, Any]]) -> dict[str, Any]:
    t1 = table1_for_cohort(old_units)
    discrepancies = []
    for label, paper in (("pooled", PAPER_POOLED_OLD), ("balanced_8of8", PAPER_BALANCED_OLD)):
        for pred, (pt, lo, hi) in paper.items():
            r = t1[label][pred]
            ok = (r["tau_b"] is not None and abs(r["tau_b"] - pt) < 0.005
                  and abs((r["ci_lo_2p5"] or float("nan")) - lo) < 0.01
                  and abs((r["ci_hi_97p5"] or float("nan")) - hi) < 0.01)
            if not ok:
                discrepancies.append({
                    "label": label, "predictor": pred, "paper": [pt, lo, hi],
                    "reproduced": [r["tau_b"], r["ci_lo_2p5"], r["ci_hi_97p5"]],
                })
    # within-state
    r = t1["within_state"]["radius"]
    pt, lo, hi = PAPER_WITHIN_OLD
    ok = (r["tau_b"] is not None and abs(r["tau_b"] - pt) < 0.005
          and abs((r["ci_lo_2p5"] or float("nan")) - lo) < 0.02
          and abs((r["ci_hi_97p5"] or float("nan")) - hi) < 0.02)
    if not ok:
        discrepancies.append({
            "label": "within_state", "predictor": "radius", "paper": [pt, lo, hi],
            "reproduced": [r["tau_b"], r["ci_lo_2p5"], r["ci_hi_97p5"]],
        })
    # paired leads
    pl = paired_leads_pooled(old_units)
    for other, (d, lo, hi) in PAPER_PAIRED_OLD.items():
        rr = pl[other]
        ok = (rr["delta_point"] is not None and abs(rr["delta_point"] - d) < 0.01
              and abs((rr["ci_lo_2p5"] or float("nan")) - lo) < 0.02
              and abs((rr["ci_hi_97p5"] or float("nan")) - hi) < 0.02)
        if not ok:
            discrepancies.append({
                "label": "paired_pooled", "predictor": other, "paper": [d, lo, hi],
                "reproduced": [rr["delta_point"], rr["ci_lo_2p5"], rr["ci_hi_97p5"]],
            })
    return {"table1": t1, "paired_leads": pl, "discrepancies": discrepancies}


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def _fmt(tau: float | None, lo: float | None, hi: float | None) -> str:
    if tau is None:
        return "n/a"
    return f"{tau:+.3f} [{lo:.3f}, {hi:.3f}]" if lo is not None and hi is not None else f"{tau:+.3f} [n/a]"


def _json_default(o: Any) -> Any:
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, set):
        return sorted(o)
    return str(o)


def main() -> int:
    t0 = time.perf_counter()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    old_units, new_units_94, combined = build_cohorts()
    print(f"old cohort: {len(old_units)} rows loaded from {OLD_UNITS_CSV}")
    print(f"94 empty-valid real-naming units loaded from {NEW_UNITS_CSV}")
    print(f"combined cohort: {len(combined)} rows")

    print("== Step 1: reproduce OLD cohort Table 1 ==")
    step1 = reproduce_old(old_units)
    print(f"  discrepancies: {len(step1['discrepancies'])}")
    for d in step1["discrepancies"]:
        print(f"    MISMATCH {d}")
    if step1["discrepancies"]:
        print("STOPPING: could not reproduce paper's old-cohort numbers.")
        with (OUT_DIR / "table1_corrected.json").open("w", encoding="utf-8") as fh:
            json.dump({"status": "FAILED_REPRODUCTION", "step1": step1}, fh, indent=2, default=_json_default)
        return 1

    print("== Step 2: corrected cohort (old + 94) Table 1 / paired leads / Table 11 ==")
    t1_new = table1_for_cohort(combined)
    pl_new = paired_leads_pooled(combined)
    t11_new = table11(combined)

    # sanity vs paper's reported corrected numbers
    corrected_checks = []
    r = t1_new["pooled"]["radius"]
    pt, lo, hi = PAPER_CORRECTED["pooled_radius"]
    corrected_checks.append({
        "predictor": "pooled radius", "paper": [pt, lo, hi],
        "reproduced": [r["tau_b"], r["ci_lo_2p5"], r["ci_hi_97p5"]], "n": r["n"],
    })
    for other, (d, lo, hi) in PAPER_CORRECTED["pooled_paired"].items():
        rr = pl_new[other]
        corrected_checks.append({
            "predictor": f"paired radius-{other}", "paper": [d, lo, hi],
            "reproduced": [rr["delta_point"], rr["ci_lo_2p5"], rr["ci_hi_97p5"]],
        })
    rw = t1_new["within_state"]["radius"]
    pt, lo, hi = PAPER_CORRECTED["within_radius"]
    corrected_checks.append({
        "predictor": "within-state radius", "paper": [pt, lo, hi],
        "reproduced": [rw["tau_b"], rw["ci_lo_2p5"], rw["ci_hi_97p5"]],
        "n_states": rw.get("n_states"),
    })

    elapsed = round(time.perf_counter() - t0, 1)
    print(f"total elapsed: {elapsed}s")

    out = {
        "old_cohort": {
            "table1": step1["table1"], "paired_leads": step1["paired_leads"],
            "reproduction_discrepancies": step1["discrepancies"],
        },
        "corrected_cohort": {
            "n_old": len(old_units), "n_added": len(new_units_94), "n_total": len(combined),
            "table1": t1_new, "paired_leads": pl_new, "table11": t11_new,
            "checks_vs_paper_reported_corrected_numbers": corrected_checks,
        },
        "n_boot": N_BOOT, "seed": SEED, "elapsed_seconds": elapsed,
    }
    with (OUT_DIR / "table1_corrected.json").open("w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, default=_json_default)

    # markdown rendering
    lines = []
    lines.append("# Corrected-cohort Table 1 and Table 11\n")
    lines.append(f"Old cohort reproduced exactly: discrepancies={len(step1['discrepancies'])}\n")
    lines.append("## Table 1\n")
    lines.append("| predictor | pooled (old) | pooled (corrected) | balanced (old) | balanced (corrected) | within-state (old) | within-state (corrected) |")
    lines.append("|---|---|---|---|---|---|---|")
    for pred in PREDICTORS:
        old_p = step1["table1"]["pooled"].get(pred, {})
        new_p = t1_new["pooled"].get(pred, {})
        old_b = step1["table1"]["balanced_8of8"].get(pred, {})
        new_b = t1_new["balanced_8of8"].get(pred, {})
        old_w = step1["table1"]["within_state"].get(pred, {})
        new_w = t1_new["within_state"].get(pred, {})
        lines.append(
            f"| {PREDICTOR_LABEL[pred]} "
            f"| {_fmt(old_p.get('tau_b'), old_p.get('ci_lo_2p5'), old_p.get('ci_hi_97p5'))} "
            f"| {_fmt(new_p.get('tau_b'), new_p.get('ci_lo_2p5'), new_p.get('ci_hi_97p5'))} "
            f"| {_fmt(old_b.get('tau_b'), old_b.get('ci_lo_2p5'), old_b.get('ci_hi_97p5'))} "
            f"| {_fmt(new_b.get('tau_b'), new_b.get('ci_lo_2p5'), new_b.get('ci_hi_97p5'))} "
            f"| {_fmt(old_w.get('tau_b'), old_w.get('ci_lo_2p5'), old_w.get('ci_hi_97p5'))} "
            f"| {_fmt(new_w.get('tau_b'), new_w.get('ci_lo_2p5'), new_w.get('ci_hi_97p5'))} |"
        )
    lines.append("")
    lines.append(f"pooled n (old): {step1['table1']['pooled_n']} / {step1['table1']['pooled_n_networks']} networks")
    lines.append(f"pooled n (corrected): {t1_new['pooled_n']} / {t1_new['pooled_n_networks']} networks")
    lines.append(f"balanced n (old): {step1['table1']['balanced_n']} / {step1['table1']['balanced_n_networks']} networks")
    lines.append(f"balanced n (corrected): {t1_new['balanced_n']} / {t1_new['balanced_n_networks']} networks")
    lines.append(f"pooled separation-defined n (old): {step1['table1']['pooled_separation_defined_n']} / {step1['table1']['pooled_separation_defined_n_networks']} networks")
    lines.append(f"pooled separation-defined n (corrected): {t1_new['pooled_separation_defined_n']} / {t1_new['pooled_separation_defined_n_networks']} networks")
    lines.append("")
    lines.append("## Pooled paired leads of r_val (network-cluster bootstrap)\n")
    lines.append("| baseline | old delta | corrected delta |")
    lines.append("|---|---|---|")
    for other in ("shd_truth", "k_g0", "n_k"):
        o = step1["paired_leads"][other]
        n = pl_new[other]
        lines.append(f"| {other} | {_fmt(o['delta_point'], o['ci_lo_2p5'], o['ci_hi_97p5'])} | {_fmt(n['delta_point'], n['ci_lo_2p5'], n['ci_hi_97p5'])} |")
    lines.append("")
    lines.append("## Within-state matched-separation comparison (r_val vs s, on units where s is defined)\n")
    ms_old = step1["table1"]["within_state_matched_separation"]
    ms = t1_new["within_state_matched_separation"]
    lines.append("| | old | corrected |")
    lines.append("|---|---|---|")
    lines.append(
        f"| radius (n_units, n_states) "
        f"| {_fmt(ms_old['radius']['tau_b'], ms_old['radius']['ci_lo_2p5'], ms_old['radius']['ci_hi_97p5'])} (n={ms_old['radius']['n_units']}, states={ms_old['radius']['n_states']}) "
        f"| {_fmt(ms['radius']['tau_b'], ms['radius']['ci_lo_2p5'], ms['radius']['ci_hi_97p5'])} (n={ms['radius']['n_units']}, states={ms['radius']['n_states']}) |"
    )
    lines.append(
        f"| separation (n_units, n_states) "
        f"| {_fmt(ms_old['separation']['tau_b'], ms_old['separation']['ci_lo_2p5'], ms_old['separation']['ci_hi_97p5'])} (n={ms_old['separation']['n_units']}, states={ms_old['separation']['n_states']}) "
        f"| {_fmt(ms['separation']['tau_b'], ms['separation']['ci_lo_2p5'], ms['separation']['ci_hi_97p5'])} (n={ms['separation']['n_units']}, states={ms['separation']['n_states']}) |"
    )
    if ms.get("paired_lead_available"):
        pld = ms["paired_lead"]
        lines.append(f"paired lead (corrected): {_fmt(pld.get('delta_point'), pld.get('ci_lo_2p5'), pld.get('ci_hi_97p5'))}")
    else:
        lines.append("(paired within-state lead delta not computed by this script -- "
                      "no `_paired_within_state_bootstrap` helper existed in llm_paired_v2.py; "
                      "see individual tau_b's above for the r_val vs s comparison instead.)")
    lines.append("")
    lines.append("## Table 11 (corrected cohort, >= k/8 scorable, pooled)\n")
    lines.append("| k/8 | n_units | n_networks | " + " | ".join(PREDICTOR_LABEL[p] for p in PREDICTORS) + " |")
    lines.append("|---" * (3 + len(PREDICTORS)) + "|")
    for row in t11_new:
        cells = " | ".join(_fmt(row[p]["tau_b"], row[p]["ci_lo_2p5"], row[p]["ci_hi_97p5"]) for p in PREDICTORS)
        lines.append(f"| {row['threshold_k_of_8']} | {row['n_units']} | {row['n_networks']} | {cells} |")
    lines.append("")
    lines.append("## Checks vs paper-reported corrected numbers\n")
    for c in corrected_checks:
        lines.append(f"- {c}")

    with (OUT_DIR / "table1_corrected.md").open("w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    print("wrote results/axis_robustness_llm_v2/table1_corrected.json and .md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
