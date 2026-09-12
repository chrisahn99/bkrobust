"""Review round 1, weakness 4, on the anchor table's own endpoints.

`build_tau_comparisons.PREDICTORS` is ("r_val", "shd_truth", "n_k"). None of
those is a cheap, truth-free structural quantity, so the published comparison
never asked whether the radius beats one. This re-runs the anchor builder with
`separation` added, using the anchor's own per-arm endpoints (AUC_frac_usable
for flip, AUC_rate_usable for tiered), so the new column is row-comparable with
the published ones.

Guard: the committed r_val / shd_truth / n_k rows must come back unchanged.
"""
from __future__ import annotations

import pathlib
import sys

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from bkrobust.robustness import build_tau_comparisons as btc  # noqa: E402

RESULTS = ROOT / "results" / "axis_robustness"
OUT = RESULTS / "review1_anchor_with_separation.csv"
CMP = ["tau_b", "ci_lo_2p5", "ci_hi_97p5", "n"]
KEY = ["stratum_label", "predictor", "endpoint"]


def main() -> None:
    committed = pd.read_csv(RESULTS / "table_tau_source.csv")

    original = list(btc.PREDICTORS)
    try:
        btc.PREDICTORS = original + ["separation"]
        annotated, _doc = btc.build(str(RESULTS))
    finally:
        btc.PREDICTORS = original

    a = committed.set_index(KEY)[CMP].sort_index()
    b = annotated.set_index(KEY)[CMP].sort_index()
    shared = a.index.intersection(b.index)
    if len(shared) < len(a):
        raise SystemExit(f"ABORT: {len(a) - len(shared)} committed rows vanished")
    delta = (a.loc[shared] - b.loc[shared]).abs().max().max()
    if delta > 1e-9:
        raise SystemExit(f"ABORT: committed rows moved by {delta}")
    print(f"guard OK: all {len(shared)} committed anchor rows reproduced to {delta:.2e}")

    annotated.to_csv(OUT, index=False)
    print(f"wrote {OUT}\n")

    main_rows = annotated[annotated["is_supplementary"] == False]  # noqa: E712
    w = main_rows.pivot_table(index="stratum_label", columns="predictor", values="tau_b")
    lo = main_rows.pivot_table(index="stratum_label", columns="predictor", values="ci_lo_2p5")
    hi = main_rows.pivot_table(index="stratum_label", columns="predictor", values="ci_hi_97p5")
    wins = 0
    print(f"{'stratum':38s} {'r_hop':>22s}   {'separation':>22s}   verdict")
    for s in w.index:
        rv, sp = w.loc[s, "r_val"], w.loc[s, "separation"]
        f_rv = f"{rv:+.3f} [{lo.loc[s,'r_val']:+.2f},{hi.loc[s,'r_val']:+.2f}]"
        f_sp = f"{sp:+.3f} [{lo.loc[s,'separation']:+.2f},{hi.loc[s,'separation']:+.2f}]"
        better = abs(rv) >= abs(sp)
        wins += better
        print(f"{s:38s} {f_rv:>22s}   {f_sp:>22s}   {'radius' if better else 'SEPARATION'}")
    print(f"\n|r_hop| >= |separation| in {wins}/{len(w.index)} strata")


if __name__ == "__main__":
    main()
