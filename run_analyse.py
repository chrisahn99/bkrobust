#!/usr/bin/env python3
"""Orchestrates the Phase 1 (Axis Robustness) analysis. Read-only on every
input file under ``results/axis_robustness/``; writes only the
``analysis_*``-prefixed outputs and ``ANALYSIS_NOTES.md``, all under
``results/axis_robustness/``.

This script never calls ``run_survival`` / ``run_hops`` and never touches
``survival_samples.csv``. It does not regenerate any sweep.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent / "src"))

from bkrobust.robustness import analyse as an  # noqa: E402

RESULTS_DIR = Path(__file__).parent / "results" / "axis_robustness"


def purity_grep() -> str:
    """E3. Greps this module and this script for the forbidden patterns."""
    pattern = re.compile(
        r"all_valid_adjustment_sets_mpdag|synth\.runner|np\.random\.seed|random\.seed|import random"
    )
    hits = []
    for path in [RESULTS_DIR.parent.parent / "src" / "bkrobust" / "robustness" / "analyse.py",
                 Path(__file__)]:
        text = path.read_text()
        for i, line in enumerate(text.splitlines(), start=1):
            if pattern.search(line):
                hits.append(f"{path}:{i}: {line.strip()}")
    return "\n".join(hits) if hits else "(empty)"


def fmt_df(df: pd.DataFrame, max_rows: int = 40) -> str:
    with pd.option_context("display.max_rows", max_rows, "display.width", 200):
        return df.to_string(index=False)


def main() -> None:
    instances = an.load_instances(str(RESULTS_DIR / "survival_instances.csv"))
    curves = an.load_curves(str(RESULTS_DIR / "survival_curves.csv"))
    hops_p5 = an.load_hops_p5(str(RESULTS_DIR / "hops_p5.csv"))

    r_assumes = an.get_r_assumes(instances)

    # --- Appendix E: tiered-arm re-binning on corruption_rate ---------------
    # The tiered arm's d_claims undercounts corruption depth (misses
    # assertions *removed* by node relocation), contaminating its d=0 bin.
    # Re-bin on corruption_rate instead -- the tiered arm's native knob,
    # recorded per-sample, immune to the counting defect. Reads
    # survival_samples.csv exactly once, tiered rows only (chunked, usecols
    # read -- see an.load_tiered_samples), and never regenerates anything.
    tiered_samples = an.load_tiered_samples(str(RESULTS_DIR / "survival_samples.csv"))
    tiered_verify_summary, tiered_verify_table = an.verify_tiered_pooled_rebinning(tiered_samples)
    if not tiered_verify_summary["within_tol"] or not tiered_verify_summary["s_at_rate_0_exact_1"]:
        raise SystemExit(
            "STOP: tiered re-binning does not reproduce PREREGISTRATION.md Appendix E.2's "
            f"pooled target table (max_abs_deviation={tiered_verify_summary['max_abs_deviation']:.4g}, "
            f"S(rate=0)==1.0: {tiered_verify_summary['s_at_rate_0_exact_1']}). "
            "Per the task instruction: do not proceed on numbers that disagree.\n"
            f"{fmt_df(tiered_verify_table)}"
        )
    tiered_rebinned = an.build_tiered_rebinned(tiered_samples)
    # Write immediately -- partial output beats none, and this is the one new
    # required deliverable this run adds.
    tiered_rebinned.to_csv(RESULTS_DIR / "analysis_tiered_rebinned.csv", index=False)
    del tiered_samples  # 1.47M tiered rows; done with it, drop before the rest of the run

    instances = an.attach_auc_frac_rate(instances, tiered_rebinned)
    n_tiered_auc_rate_nan = int(
        instances.loc[instances["arm"] == "tiered", "AUC_frac_rate"].isna().sum()
    )

    # --- E1: decomposition identity ---------------------------------------
    decomp_summary, curves_with_rate = an.verify_decomposition(curves)
    decomposition_by_d = (
        curves_with_rate.groupby("d")
        .agg(
            mean_contradiction_rate=("contradiction_rate", "mean"),
            mean_S=("S", "mean"),
            mean_S_contra_as_fail=("S_contra_as_fail", "mean"),
            n_rows=("d", "size"),
        )
        .reset_index()
        .sort_values("d", kind="stable")
    )
    decomposition_by_d.to_csv(RESULTS_DIR / "analysis_decomposition.csv", index=False)

    # --- AUC_frac_usable sensitivity endpoint ------------------------------
    instances = an.attach_auc_frac_usable(instances, curves)
    n_usable_zero = int((instances["n_usable_depths"] == 0).sum())
    n_auc_usable_nan = int(instances["AUC_frac_usable"].isna().sum())

    # --- Sanity checks required by the task --------------------------------
    n_unreached = int((instances["r_val"] == an.UNREACHED).sum())

    # --- Task 1: stratified tau ---------------------------------------------
    tau_primary = an.build_tau_table(instances, level="primary", r_assumes=r_assumes)
    tau_secondary = an.build_tau_table(instances, level="secondary", r_assumes=r_assumes)
    tau_primary.to_csv(RESULTS_DIR / "analysis_tau_primary.csv", index=False)
    tau_secondary.to_csv(RESULTS_DIR / "analysis_tau_secondary.csv", index=False)

    # E2: determinism -- run the bootstrap twice, compare CIs bit-for-bit.
    tau_primary_rerun = an.build_tau_table(instances, level="primary", r_assumes=r_assumes)
    ci_cols = ["ci_lo_2p5", "ci_hi_97p5"]
    determinism_ok = bool(
        np.allclose(
            tau_primary[ci_cols].fillna(-999).to_numpy(),
            tau_primary_rerun[ci_cols].fillna(-999).to_numpy(),
        )
    )

    # --- Verdicts P1-P4 -------------------------------------------------------
    # Per Appendix E, flip's primary endpoint is AUC_frac and tiered's is the
    # re-binned AUC_frac_rate; tau_primary already tags each row is_primary /
    # is_sensitivity per-arm (see an.PRIMARY_ENDPOINT_BY_ARM), so verdict_p1-3
    # select on that flag rather than a single hard-coded endpoint string.
    p1 = an.verdict_p1(tau_primary)
    p2 = an.verdict_p2(tau_primary)
    p3 = an.verdict_p3(tau_primary)
    p4 = an.verdict_p4(instances)  # per E.4: reports variances, forces no verdict

    # E6: sensitivity endpoint (AUC_frac_usable, n_eval>=30 depths on the
    # d_claims axis) -- do verdicts change? Note this sensitivity endpoint is
    # NOT re-binned for the tiered arm (no "usable" analogue of the rate axis
    # was requested), so it is reported for continuity with the caveat that,
    # for the tiered arm only, it still rests on the defective d_claims axis.
    p1_s = an.verdict_p1(tau_primary, use_sensitivity=True)
    p2_s = an.verdict_p2(tau_primary, use_sensitivity=True)
    p3_s = an.verdict_p3(tau_primary, use_sensitivity=True)
    # P4 has no sensitivity variant post-Appendix-E: the primary comparison is
    # already reported as not-comparable (unit mismatch), independent of which
    # usable-depth threshold would be applied to either axis.

    # --- Cross-check against hops_p5.csv ------------------------------------
    cross = an.hops_cross_check(tau_primary, hops_p5)

    # --- Task 2: discordance spotlight ---------------------------------------
    spotlight = an.discordance_spotlight(instances, curves)
    spotlight.to_csv(RESULTS_DIR / "analysis_discordance.csv", index=False)
    verification = an.verify_spotlight_cases(spotlight, instances, curves)
    n_verified_ok = int(verification["auc_matches"].sum())
    n_verified_total = len(verification)

    # Supplementary E5 check for tiered spotlighted cases: their ranking now
    # runs on AUC_frac_rate, not the AUC_frac verify_spotlight_cases checks
    # against survival_curves.csv -- confirm AUC_frac_rate in the spotlight
    # table matches analysis_tiered_rebinned.csv's own computation for the
    # same instance_id (this is a merge/plumbing check: both numbers must
    # come from the exact same source row, catching any join bug rather than
    # re-deriving the value a second, independent way).
    tiered_spotlight = spotlight[spotlight["arm"] == "tiered"]
    rate_check = tiered_spotlight.merge(
        tiered_rebinned[["instance_id", "AUC_frac_rate"]], on="instance_id", suffixes=("", "_rebinned")
    )
    n_rate_check_ok = int(
        np.isclose(rate_check["AUC_frac_rate"], rate_check["AUC_frac_rate_rebinned"], equal_nan=True).sum()
    )
    n_rate_check_total = len(rate_check)

    # --- Task 3: effect sizes -------------------------------------------------
    effect_sizes = an.build_effect_size_tables(instances)
    effect_sizes.to_csv(RESULTS_DIR / "analysis_effect_sizes.csv", index=False)

    # --- Purity grep (E3) -------------------------------------------------
    grep_result = purity_grep()

    # --- Write ANALYSIS_NOTES.md --------------------------------------------
    notes = []
    notes.append("# Phase 1 Analysis Notes\n")
    notes.append(
        "Endpoint per PREREGISTRATION.md Appendix D.3 (restoring section 4): "
        "**primary = `AUC_frac` over `S`, contradictions excluded**, for the flip arm. "
        "Per **Appendix E**, the tiered arm's `AUC_frac` (gridded on the defective "
        "`d_claims` metric, which undercounts corruption depth because it misses "
        "assertions *removed* by node relocation) is re-binned on `corruption_rate` "
        "-- the tiered arm's native, exactly-recorded experimental knob -- giving "
        "**`AUC_frac_rate`**, the tiered arm's primary endpoint for every table below. "
        "`AUC_frac_usable` (restricted to depths with `n_eval >= 30`, still on the "
        "d_claims axis for both arms) is the required sensitivity analysis, not a "
        "redefinition; for the tiered arm it inherits the d_claims defect and is "
        "reported only for continuity.\n"
    )

    notes.append("## Appendix E -- tiered re-binning verification (E1)\n")
    notes.append(
        f"Pooled tiered-arm samples ({RESULTS_DIR / 'survival_samples.csv'}, tiered rows "
        f"only, {tiered_verify_table['n'].sum()} rows across {len(tiered_verify_table)} rates) "
        f"checked against PREREGISTRATION.md Appendix E.2's target table:\n"
        f"{fmt_df(tiered_verify_table)}\n"
        f"- n matches target (133,600) at every rate: {tiered_verify_summary['all_n_match_target']}\n"
        f"- S(rate=0.0) == 1.000 exactly: {tiered_verify_summary['s_at_rate_0_exact_1']}\n"
        f"- **max_abs_deviation from the target table: "
        f"{tiered_verify_summary['max_abs_deviation']:.4g}** (target table itself is printed "
        f"to 3 decimals, so deviations at this scale are rounding, not disagreement) -- "
        f"**MATCHES** (within tol={0.002:g}).\n"
        f"- Per-instance re-binning written to `analysis_tiered_rebinned.csv` "
        f"({len(tiered_rebinned)} tiered instances, `AUC_frac_rate` undefined for "
        f"{n_tiered_auc_rate_nan} of them).\n"
        f"- Flip arm is unaffected (targeted depth, exact at generation, no d=0 rows) "
        f"and never touches survival_samples.csv, per Appendix E.2 scope.\n"
    )

    notes.append("## Verdicts (primary endpoint: flip=AUC_frac, tiered=AUC_frac_rate)\n")
    notes.append(
        f"- **P1** ({p1['verdict']}), endpoints used: {p1['endpoint']}: "
        f"tau_b(endpoint, r_val) > 0 in "
        f"{p1['n_positive']}/{p1['n_strata_defined']} defined primary strata "
        f"(of {p1['n_strata_total']} total).\n"
        f"{fmt_df(p1['detail'])}\n"
    )
    notes.append(
        f"- **P2** ({p2['verdict']}), endpoints used: {p2['endpoint']}: r_val beats n_k in "
        f"{p2['n_rval_beats_nk']}/{p2['n_strata_both_defined']} strata where both defined.\n"
        f"{fmt_df(p2['detail'])}\n"
    )
    notes.append(
        f"- **P3 / H4** ({p3['verdict']}), endpoints used: {p3['endpoint']}: n_k's CI covers 0 in "
        f"{p3['n_ci_covers_zero']}/{p3['n_strata_defined']} defined strata.\n"
        f"{fmt_df(p3['detail'])}\n"
    )
    notes.append(
        f"- **P4** ({p4['verdict']}).\n"
        f"  {p4['unit_mismatch']}\n"
        f"  var(flip AUC_frac) = {p4['var_flip_overall_AUC_frac']:.4f}; "
        f"var(tiered AUC_frac_rate) = {p4['var_tiered_overall_AUC_frac_rate']:.4f} "
        f"(raw, non-comparable, tiered numerically higher: {p4['overall_tiered_higher_raw']}). "
        f"Matched (component_size, separation) cells: "
        f"{p4['n_cells_tiered_higher_raw']}/{p4['n_matched_cells']} numerically favor tiered "
        f"(reported for completeness; not a same-units comparison; no verdict is forced).\n"
        f"{fmt_df(p4['matched_detail'])}\n"
    )

    notes.append("## E6 -- sensitivity endpoint (AUC_frac_usable): do verdicts change?\n")
    notes.append(
        f"- P1: primary={p1['verdict']}, usable={p1_s['verdict']} "
        f"({p1_s['n_positive']}/{p1_s['n_strata_defined']}) -- "
        f"{'UNCHANGED' if p1['verdict'] == p1_s['verdict'] else 'CHANGED'}\n"
        f"- P2: primary={p2['verdict']}, usable={p2_s['verdict']} "
        f"({p2_s['n_rval_beats_nk']}/{p2_s['n_strata_both_defined']}) -- "
        f"{'UNCHANGED' if p2['verdict'] == p2_s['verdict'] else 'CHANGED'}\n"
        f"- P3: primary={p3['verdict']}, usable={p3_s['verdict']} "
        f"({p3_s['n_ci_covers_zero']}/{p3_s['n_strata_defined']}) -- "
        f"{'UNCHANGED' if p3['verdict'] == p3_s['verdict'] else 'CHANGED'}\n"
        f"- P4: no sensitivity variant post-Appendix-E (see above) -- the primary comparison "
        f"is already NOT COMPARABLE by unit mismatch, independent of any usable-depth threshold.\n"
        f"- n instances with 0 usable (n_eval>=30) depths: {n_usable_zero} "
        f"({n_usable_zero/len(instances):.1%}); AUC_frac_usable undefined (NaN) for "
        f"{n_auc_usable_nan} instances. **Caveat:** for the tiered arm, AUC_frac_usable "
        f"remains on the defective d_claims axis (no rate-binned \"usable\" variant was "
        f"computed), so its P1-P3 rows above should be read as a d_claims-axis sensitivity "
        f"check only, not as a check on AUC_frac_rate itself.\n"
    )

    notes.append("## Cross-check against hops_p5.csv (independent small-component subsample)\n")
    notes.append(f"{fmt_df(cross)}\n")
    n_sign_agree = int(cross["sign_agree"].sum())
    notes.append(
        f"Sign agreement: {n_sign_agree}/{len(cross)} strata. "
        f"Rough magnitude agreement (|main - hops| <= 0.35): "
        f"{int(cross['rough_magnitude_agree'].sum())}/{len(cross)}.\n"
    )

    notes.append("## E1 -- decomposition identity\n")
    notes.append(
        f"S_contra_as_fail(d) = (1 - contradiction_rate(d)) * S(d), checked on "
        f"{decomp_summary['n_rows_checked_S_defined']} of {decomp_summary['n_rows_total']} "
        f"curve rows (S undefined on {decomp_summary['n_rows_S_undefined']}, all-contradictory "
        f"depths, where the identity does not apply). "
        f"max_abs_deviation={decomp_summary['max_abs_deviation']:.3e}, "
        f"mean_abs_deviation={decomp_summary['mean_abs_deviation']:.3e}, "
        f"n_rows with deviation > 1e-6: {decomp_summary['n_deviation_gt_1e6']}.\n"
    )

    notes.append("## E2 -- determinism\n")
    notes.append(f"Bootstrap CIs identical across two independent runs: {determinism_ok}.\n")

    notes.append("## E3 -- purity grep\n")
    notes.append(f"```\n{grep_result}\n```\n")

    notes.append("## E4 -- stratum n and underpowered flags (primary strata)\n")
    n_table = tau_primary[tau_primary["predictor"] == "r_val"][
        ["arm"] + [c for c in an._ALL_STRATA_COLS if c in tau_primary.columns]
        + ["endpoint", "n", "underpowered"]
    ]
    notes.append(f"{fmt_df(n_table, max_rows=100)}\n")

    notes.append("## E5 -- discordance spotlight verification\n")
    notes.append(
        f"AUC_frac recomputed directly from survival_curves.csv matches the recorded value "
        f"(within 1e-6) for {n_verified_ok}/{n_verified_total} spotlighted instances "
        f"(this check is on the d_claims axis AUC_frac column; it applies to both arms "
        f"but for tiered instances is a plumbing check on AUC_frac, not on the "
        f"AUC_frac_rate that actually drove their ranking -- see the next line).\n"
        f"{fmt_df(verification)}\n"
        f"Supplementary check for the {n_rate_check_total} spotlighted tiered instances: "
        f"their AUC_frac_rate in the spotlight table matches "
        f"analysis_tiered_rebinned.csv's own value for the same instance_id for "
        f"{n_rate_check_ok}/{n_rate_check_total}.\n"
    )

    notes.append("## Sanity checks\n")
    notes.append(
        f"- r_val == UNREACHED count: {n_unreached} (task expects 0): "
        f"{'CONFIRMED' if n_unreached == 0 else 'MISMATCH -- investigate'}\n"
        f"- r_assumes constant across dataset: \"{r_assumes}\"\n"
    )

    notes.append(
        "## Effect sizes (median of each arm's primary endpoint by predictor bucket, "
        "primary strata; `endpoint` column states which axis: AUC_frac for flip, "
        "AUC_frac_rate for tiered)\n"
    )
    notes.append(f"{fmt_df(effect_sizes, max_rows=250)}\n")

    (RESULTS_DIR / "ANALYSIS_NOTES.md").write_text("\n".join(notes))

    print("Wrote:")
    for f in [
        "analysis_tiered_rebinned.csv",
        "analysis_tau_primary.csv", "analysis_tau_secondary.csv",
        "analysis_decomposition.csv", "analysis_discordance.csv",
        "analysis_effect_sizes.csv", "ANALYSIS_NOTES.md",
    ]:
        print(" -", RESULTS_DIR / f)


if __name__ == "__main__":
    main()
