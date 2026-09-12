"""Review round 1, weakness 4: does the radius add anything over cheap structural
baselines an analyst can compute without the truth?

`analyse.PREDICTORS` never contained `separation`, `component_size` or `z_size`.
The published comparison therefore pits the radius only against `shd_truth` (which
needs the true DAG and so is not a baseline a practitioner could ever use) and
`n_k` (pre-registered as inert). This script re-runs the committed tau pipeline
with the structural predictors added, on the committed CSVs and with no
re-simulation.

Guard: it first reproduces the committed `r_val` rows exactly. If that check
fails the script aborts rather than report new numbers from a pipeline that no
longer matches the published one.
"""
from __future__ import annotations

import pathlib
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from bkrobust.robustness import analyse as an  # noqa: E402

RESULTS = pathlib.Path(__file__).resolve().parent.parent / "results" / "axis_robustness"
OUT = RESULTS / "review1_structural_tau.csv"

#: Cheap, truth-free structural quantities. `separation` is the hop distance from
#: the treatment to the nearest member of the adjustment set; the project's own
#: appendix law says it nearly determines the radius on a designed family.
EXTRA = ["separation", "component_size", "z_size"]

KEY = ["arm", "coverage", "base_wrongness", "n_tiers", "predictor", "endpoint"]
CMP = ["tau_b", "ci_lo_2p5", "ci_hi_97p5", "n"]


def main() -> None:
    instances = pd.read_csv(RESULTS / "survival_instances.csv")
    curves = pd.read_csv(RESULTS / "survival_curves.csv")
    instances = an.attach_auc_frac_usable(instances, curves)
    rebinned = pd.read_csv(RESULTS / "analysis_tiered_rebinned.csv")
    rate_col = [c for c in rebinned.columns if c.startswith("AUC_frac_rate")]
    if rate_col:
        instances = instances.merge(
            rebinned[["instance_id"] + rate_col[:1]], on="instance_id", how="left"
        )
    r_assumes = "Conjecture 2 (hence Anti-Exchange Case B, verified not proved)"

    committed = pd.read_csv(RESULTS / "analysis_tau_primary.csv")

    original = list(an.PREDICTORS)
    try:
        an.PREDICTORS = original + [p for p in EXTRA if p not in original]
        table = an.build_tau_table(instances, level="primary", r_assumes=r_assumes)
    finally:
        an.PREDICTORS = original

    # --- guard: the committed r_val rows must come back unchanged -------------
    def keyed(df: pd.DataFrame) -> pd.DataFrame:
        d = df[df["predictor"] == "r_val"].copy()
        cols = [c for c in KEY if c in d.columns]
        return d.set_index(cols)[CMP].sort_index()

    a, b = keyed(committed), keyed(table)
    shared = a.index.intersection(b.index)
    if len(shared) == 0:
        raise SystemExit("ABORT: no shared r_val rows; pipeline key mismatch")
    delta = (a.loc[shared] - b.loc[shared]).abs().max().max()
    if delta > 1e-9:
        raise SystemExit(f"ABORT: r_val rows moved by {delta}; pipeline no longer matches")
    print(f"guard OK: {len(shared)} committed r_val rows reproduced to {delta:.2e}")

    table.to_csv(OUT, index=False)
    print(f"wrote {OUT}")

    # --- the comparison the review asked for ----------------------------------
    main_rows = table[
        (table.get("is_primary", True)) & (table["endpoint"].str.contains("usable"))
    ]
    if main_rows.empty:
        main_rows = table[table["endpoint"].str.contains("usable")]
    wide = main_rows.pivot_table(
        index=[c for c in ["arm", "coverage", "base_wrongness", "n_tiers"] if c in main_rows],
        columns="predictor",
        values="tau_b",
        dropna=False,
    )
    cols = [c for c in ["r_val"] + EXTRA + ["shd_truth", "n_k"] if c in wide.columns]
    print("\ntau_b by stratum, conservative endpoint:\n")
    print(wide[cols].round(3).to_string())

    beats = {}
    for c in EXTRA:
        if c in wide.columns:
            ok = (wide["r_val"].abs() >= wide[c].abs()) | wide[c].isna()
            beats[c] = f"{int(ok.sum())}/{int(ok.notna().sum())} strata where |r_val| >= |{c}|"
    print("\n" + "\n".join(f"  {k}: {v}" for k, v in beats.items()))


if __name__ == "__main__":
    main()
