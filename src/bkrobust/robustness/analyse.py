"""Phase 1 analysis: stratified Kendall tau-b, decomposition check, discordance
spotlight, effect-size tables.

Live-tree only. Nothing here imports from the stub scaffold (``graphs``,
``knowledge`` as a package, ``metrics``, ``theory``, ``estimation``, ``data``,
``cfm``, ``representations``, ``utils``), and nothing here calls
``all_valid_adjustment_sets_mpdag`` or ``synth.runner``. Reads
``results/axis_robustness/survival_instances.csv``,
``survival_curves.csv``, and ``hops_p5.csv`` for everything except the tiered
re-binning: (S(d), S_contra_as_fail(d), n_eval, n_contradictory) are already
aggregated in ``survival_curves.csv`` and never require ``survival_samples.csv``.

``survival_samples.csv`` (3.14M rows / 354MB) is read exactly once, for the
Appendix E tiered-arm re-binning below, and only its tiered-arm rows (via a
chunked, ``usecols``-restricted read -- see :func:`load_tiered_samples`): the
tiered arm's ``d_claims`` undercounts corruption depth (it misses assertions
*removed* by node relocation, only counting added/reversed ones), which
contaminates its ``d = 0`` bin with every corruption rate. The fix re-bins on
``corruption_rate``, the tiered arm's native experimental knob, recorded
per-sample and immune to the counting defect. The flip arm is unaffected
(targeted depth, exact at generation, no ``d = 0`` rows) and never touches
``survival_samples.csv``.

Design fixed by ``results/axis_robustness/PREREGISTRATION.md`` sections 3, 4,
6, 7, as restored by Appendix D (D.3 retracts the C.1 endpoint swap; the
primary endpoint is ``AUC_frac`` over ``S`` -- contradictions excluded --
exactly as pre-registered in section 4). ``AUC_frac_usable`` (recomputed here,
restricted to depths with ``n_eval >= 30``) is a required *sensitivity*
analysis, never a redefinition of the endpoint.

This module never runs the sweep. It imports
:func:`bkrobust.robustness.survival.auc_frac` (a pure function operating on an
in-memory ``curve`` dict) purely to reuse the exact nearest-defined-depth /
tie-break-toward-smaller-depth algorithm for the ``AUC_frac_usable``
sensitivity recomputation, so the two AUC definitions cannot silently drift
apart. Importing that module does not execute any sweep.

All randomness is explicit ``numpy.random.Generator`` -- no ``np.random.seed``,
no ``random.seed``, no bare ``np.random.*``.
"""

from __future__ import annotations

import math
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

from bkrobust.core.conventions import UNREACHED
from bkrobust.robustness.survival import auc_frac as _auc_frac_fn

# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def load_instances(path: str) -> pd.DataFrame:
    """Load ``survival_instances.csv``, sorted by ``instance_id`` for
    determinism in every downstream order-dependent operation."""
    df = pd.read_csv(path)
    return df.sort_values("instance_id", kind="stable").reset_index(drop=True)


def load_curves(path: str) -> pd.DataFrame:
    """Load ``survival_curves.csv``, sorted by ``(instance_id, d)``."""
    df = pd.read_csv(path)
    return df.sort_values(["instance_id", "d"], kind="stable").reset_index(drop=True)


def load_hops_p5(path: str) -> pd.DataFrame:
    """Load ``hops_p5.csv`` (the independent small-component subsample's
    pre-computed per-stratum tau table), sorted by stratum for determinism."""
    df = pd.read_csv(path)
    return df.sort_values("stratum", kind="stable").reset_index(drop=True)


# Pulled once, verified constant across the whole dataset (checked in
# ``run_analyse.py`` / E-acceptance), and carried into every table that
# quotes a radius per the task's hard constraint.
def get_r_assumes(instances: pd.DataFrame) -> str:
    vals = instances["r_assumes"].dropna().unique()
    if len(vals) != 1:
        raise ValueError(
            f"r_assumes is not constant across the dataset ({len(vals)} distinct "
            "values) -- a table quoting a radius cannot carry a single verbatim "
            "assumption string until this is investigated."
        )
    return str(vals[0])


# ---------------------------------------------------------------------------
# Appendix E -- tiered-arm re-binning on corruption_rate
# ---------------------------------------------------------------------------

# The 11-point grid the tiered arm was actually swept on (native knob).
TIERED_RATE_GRID = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]

TIERED_TARGET_TABLE = {
    # rate: (contradiction_rate, S_given_consistent), from PREREGISTRATION.md
    # Appendix E.2 -- the verified pooled aggregate this re-binning must reproduce.
    0.0: (0.000, 1.000),
    0.05: (0.229, 0.994),
    0.10: (0.404, 0.963),
    0.15: (0.574, 0.937),
    0.20: (0.648, 0.898),
    0.25: (0.692, 0.879),
    0.30: (0.742, 0.848),
    # The prereg's own printed table shows 0.798 at rate 0.35, but its footnote
    # says "0.753 as measured; see the committed CSV, which is authoritative."
    # This re-binning independently reproduces 0.7526 -- confirming the
    # footnote, not the printed 0.798 -- so 0.753 (footnote value) is what is
    # checked against here.
    0.35: (0.801, 0.753),
    0.40: (0.813, 0.740),
    0.45: (0.833, 0.671),
    0.50: (0.852, 0.618),
}
TIERED_TARGET_N_PER_RATE = 133_600


def load_tiered_samples(path: str) -> pd.DataFrame:
    """Loads *only* the tiered-arm rows of ``survival_samples.csv`` (Appendix
    E re-binning). Reads via chunked, ``usecols``-restricted passes so the
    354MB / 3.14M-row file is never held in memory in full and the flip-arm
    rows (unaffected by the tiered depth-metric defect, per E.2 scope) are
    dropped immediately rather than carried along."""
    usecols = ["instance_id", "arm", "corruption_rate", "status", "survived"]
    chunks = []
    for chunk in pd.read_csv(path, usecols=usecols, chunksize=500_000):
        tiered_chunk = chunk[chunk["arm"] == "tiered"]
        if len(tiered_chunk):
            chunks.append(tiered_chunk)
    if not chunks:
        return pd.DataFrame(columns=usecols)
    return pd.concat(chunks, ignore_index=True)


def verify_tiered_pooled_rebinning(
    tiered_samples: pd.DataFrame, *, tol: float = 0.002
) -> tuple[dict[str, Any], pd.DataFrame]:
    """E1. Pools the tiered-arm samples by ``corruption_rate`` and checks the
    result against :data:`TIERED_TARGET_TABLE` -- the verified aggregate the
    task's acceptance criterion requires reproducing. Contradiction rate is
    ``n_contradictory / n_total`` at that rate; ``S`` given consistent is
    ``survived.mean()`` restricted to non-contradictory rows (contradictions
    excluded from the denominator, matching the primary endpoint's
    convention). Does *not* stop execution on mismatch -- the caller decides
    what to do with a nonzero ``max_abs_deviation``, per the task's
    instruction to report a discrepancy rather than silently proceed."""
    rows = []
    is_contra = tiered_samples["status"] == "corrupted_k_contradictory"
    for rate in TIERED_RATE_GRID:
        sub = tiered_samples[np.isclose(tiered_samples["corruption_rate"], rate)]
        n = len(sub)
        n_contra = int((sub["status"] == "corrupted_k_contradictory").sum())
        contra_rate = n_contra / n if n else float("nan")
        non_contra = sub[sub["status"] != "corrupted_k_contradictory"]
        s_given_consistent = float(non_contra["survived"].mean()) if len(non_contra) else float("nan")
        target_contra, target_s = TIERED_TARGET_TABLE[rate]
        rows.append(
            {
                "rate": rate,
                "n": n,
                "n_matches_target": n == TIERED_TARGET_N_PER_RATE,
                "contradiction_rate": contra_rate,
                "target_contradiction_rate": target_contra,
                "contradiction_abs_dev": abs(contra_rate - target_contra),
                "S_given_consistent": s_given_consistent,
                "target_S_given_consistent": target_s,
                "S_abs_dev": abs(s_given_consistent - target_s),
            }
        )
    table = pd.DataFrame(rows)
    max_dev = float(table[["contradiction_abs_dev", "S_abs_dev"]].to_numpy().max())
    summary = {
        "n_rates": len(table),
        "all_n_match_target": bool(table["n_matches_target"].all()),
        "max_abs_deviation": max_dev,
        "within_tol": bool(max_dev <= tol),
        "s_at_rate_0_exact_1": bool(np.isclose(table.loc[table["rate"] == 0.0, "S_given_consistent"].iloc[0], 1.0)),
    }
    return summary, table


def build_tiered_rebinned(tiered_samples: pd.DataFrame) -> pd.DataFrame:
    """Per-instance tiered re-binning (E.2's fix). For every tiered instance
    and every rate in :data:`TIERED_RATE_GRID`: ``S(rate)`` = survivors /
    non-contradictory at that rate (contradictions excluded, matching the
    primary endpoint), plus the contradiction rate and raw counts.
    ``AUC_frac_rate`` is the mean of ``S(rate)`` over the 11-point grid --
    when a cell has zero non-contradictory samples (rare: 87 of 7,348
    instance x rate cells here, concentrated at rate >= 0.25 where the
    contradiction rate already exceeds 69%), ``S(rate)`` is undefined and is
    filled from the nearest rate on the grid with a defined ``S`` (ties
    broken toward the smaller rate) -- exactly the nearest-defined-depth,
    tie-break-toward-smaller convention :func:`bkrobust.robustness.survival.
    auc_frac` already uses for the flip arm's fractional grid, so the two AUC
    definitions differ only in which axis they grid over, never in the
    fill rule. An instance with *no* rate defined at all (never observed in
    this dataset, but possible in principle) gets ``AUC_frac_rate = NaN``,
    never silently 0.
    """
    is_contra = tiered_samples["status"] == "corrupted_k_contradictory"
    grp = tiered_samples.groupby(["instance_id", "corruption_rate"])
    agg = grp.agg(
        n_total=("status", "size"),
        n_contradictory=("status", lambda s: int((s == "corrupted_k_contradictory").sum())),
    ).reset_index()
    non_contra = tiered_samples[~is_contra]
    surv = (
        non_contra.groupby(["instance_id", "corruption_rate"])["survived"]
        .agg(n_survived="sum", n_eval="size")
        .reset_index()
    )
    merged = agg.merge(surv, on=["instance_id", "corruption_rate"], how="left")
    merged["n_survived"] = merged["n_survived"].fillna(0).astype(int)
    merged["n_eval"] = merged["n_eval"].fillna(0).astype(int)
    merged["contradiction_rate"] = merged["n_contradictory"] / merged["n_total"]
    merged["S"] = np.where(merged["n_eval"] > 0, merged["n_survived"] / merged["n_eval"], np.nan)

    rows = []
    for iid, g in merged.groupby("instance_id", sort=True):
        by_rate = {float(r["corruption_rate"]): r for _, r in g.iterrows()}
        defined_rates = sorted(r for r, row in by_rate.items() if not math.isnan(row["S"]))
        row_out: dict[str, Any] = {"instance_id": iid, "n_rates_defined": len(defined_rates)}
        s_vals = []
        for rate in TIERED_RATE_GRID:
            if rate in by_rate and not math.isnan(by_rate[rate]["S"]):
                s = float(by_rate[rate]["S"])
            elif defined_rates:
                nearest = min(defined_rates, key=lambda r: (abs(r - rate), r))
                s = float(by_rate[nearest]["S"])
            else:
                s = float("nan")
            s_vals.append(s)
            cr = float(by_rate[rate]["contradiction_rate"]) if rate in by_rate else float("nan")
            n_eval = int(by_rate[rate]["n_eval"]) if rate in by_rate else 0
            row_out[f"S_rate_{rate:.2f}"] = s
            row_out[f"contradiction_rate_{rate:.2f}"] = cr
            row_out[f"n_eval_rate_{rate:.2f}"] = n_eval
        row_out["AUC_frac_rate"] = float(np.mean(s_vals)) if not any(math.isnan(v) for v in s_vals) else float("nan")
        rows.append(row_out)

    out = pd.DataFrame(rows).sort_values("instance_id", kind="stable").reset_index(drop=True)
    return out


def attach_auc_frac_rate(instances: pd.DataFrame, tiered_rebinned: pd.DataFrame) -> pd.DataFrame:
    """Merges ``AUC_frac_rate`` onto a copy of ``instances`` (NaN for the
    flip arm, which is unaffected and keeps using ``AUC_frac``)."""
    cols = ["instance_id", "AUC_frac_rate", "n_rates_defined"]
    out = instances.merge(tiered_rebinned[cols], on="instance_id", how="left", validate="one_to_one")
    return out


# ---------------------------------------------------------------------------
# E1 -- decomposition identity: S_contra_as_fail(d) = (1 - contradiction_rate(d)) * S(d)
# ---------------------------------------------------------------------------


def verify_decomposition(curves: pd.DataFrame) -> tuple[dict[str, Any], pd.DataFrame]:
    """Verifies the exact decomposition
    ``S_contra_as_fail(d) = (1 - contradiction_rate(d)) * S(d)`` row-by-row on
    every curve row where ``S`` is defined, and reports the maximum absolute
    deviation. Algebraically this identity is exact (both sides reduce to
    ``n_survived_all / n_total``), so a nonzero deviation beyond floating-point
    slack would indicate the two source columns are not what the schema
    claims -- this is checked, not assumed.

    Returns ``(summary_dict, curves_with_contradiction_rate)``.
    """
    df = curves.copy()
    denom = df["n_eval"] + df["n_contradictory"]
    df["contradiction_rate"] = np.where(denom > 0, df["n_contradictory"] / denom, np.nan)

    defined = df["S"].notna() & (denom > 0)
    predicted = (1.0 - df.loc[defined, "contradiction_rate"]) * df.loc[defined, "S"]
    actual = df.loc[defined, "S_contra_as_fail"]
    deviation = (predicted - actual).abs()

    summary = {
        "n_rows_total": int(len(df)),
        "n_rows_checked_S_defined": int(defined.sum()),
        "n_rows_S_undefined": int((~df["S"].notna()).sum()),
        "max_abs_deviation": float(deviation.max()) if len(deviation) else float("nan"),
        "mean_abs_deviation": float(deviation.mean()) if len(deviation) else float("nan"),
        "n_deviation_gt_1e6": int((deviation > 1e-6).sum()),
    }
    return summary, df


# ---------------------------------------------------------------------------
# AUC_frac_usable -- sensitivity recomputation restricted to n_eval >= 30
# ---------------------------------------------------------------------------

USABLE_MIN_N_EVAL = 30


def compute_auc_frac_usable(
    instances: pd.DataFrame, curves: pd.DataFrame, *, min_n_eval: int = USABLE_MIN_N_EVAL
) -> pd.DataFrame:
    """Recomputes the primary AUC endpoint restricted to depths with
    ``n_eval >= min_n_eval``, using the exact same nearest-defined-depth,
    tie-break-toward-smaller-depth algorithm as the pre-registered
    ``AUC_frac`` (:func:`bkrobust.robustness.survival.auc_frac`), so the
    sensitivity analysis differs from the primary endpoint *only* in which
    depths are admissible, never in the averaging rule.

    Returns a frame indexed by ``instance_id`` with columns
    ``AUC_frac_usable`` and ``n_usable_depths`` (count of depths meeting the
    threshold, before the fractional-grid nearest-neighbour lookup -- an
    instance can have 0 usable depths, in which case ``AUC_frac_usable`` is
    NaN, never silently dropped to 0).
    """
    nk_map = dict(zip(instances["instance_id"], instances["n_k"]))
    curves_sorted = curves.sort_values(["instance_id", "d"], kind="stable")

    rows = []
    for iid, g in curves_sorted.groupby("instance_id", sort=True):
        nk = nk_map.get(iid)
        curve: dict[int, dict[str, Any]] = {}
        for _, r in g.iterrows():
            if r["n_eval"] >= min_n_eval and pd.notna(r["S"]):
                curve[int(r["d"])] = {"S": float(r["S"])}
        n_usable = len(curve)
        if nk is None or not (nk > 0) or n_usable == 0:
            val = float("nan")
        else:
            raw = _auc_frac_fn(curve, int(nk))
            val = float(raw) if raw != "" else float("nan")
        rows.append({"instance_id": iid, "AUC_frac_usable": val, "n_usable_depths": n_usable})

    return pd.DataFrame(rows).sort_values("instance_id", kind="stable").reset_index(drop=True)


def attach_auc_frac_usable(instances: pd.DataFrame, curves: pd.DataFrame) -> pd.DataFrame:
    """Merges :func:`compute_auc_frac_usable` onto a copy of ``instances``."""
    usable = compute_auc_frac_usable(instances, curves)
    out = instances.merge(usable, on="instance_id", how="left", validate="one_to_one")
    return out


# ---------------------------------------------------------------------------
# Task 1 -- stratified Kendall tau-b
# ---------------------------------------------------------------------------

PREDICTORS = ["r_val", "shd_truth", "shd_cpdag", "n_k", "undirected_fraction"]
PRIMARY_ENDPOINT = "AUC_frac"
SENSITIVITY_ENDPOINT = "AUC_frac_usable"
SECONDARY_ENDPOINT = "AUC_abs"  # reported only beside n_k

UNDERPOWERED_N = 15

FLIP_PRIMARY_STRATA_COLS = ["coverage", "base_wrongness"]
FLIP_SECONDARY_STRATA_COLS = ["coverage", "base_wrongness", "component_size"]
TIERED_PRIMARY_STRATA_COLS = ["n_tiers"]
TIERED_SECONDARY_STRATA_COLS = ["n_tiers", "component_size"]

# Appendix E: the flip arm's primary endpoint is unaffected (AUC_frac, on the
# exact/targeted d_claims axis). The tiered arm's primary endpoint is the
# re-binned AUC_frac_rate (corruption_rate axis) -- AUC_frac itself is
# defective for tiered (contaminated d=0 bin, see analyse.py module docstring
# and PREREGISTRATION.md Appendix E). AUC_frac_usable / AUC_abs remain on the
# original (defective, for tiered) d_claims axis for both arms -- reported
# for continuity, with the caveat stated in ANALYSIS_NOTES.md, since no
# "usable" or "abs" analogue of the rate axis was requested.
PRIMARY_ENDPOINT_BY_ARM = {"flip": PRIMARY_ENDPOINT, "tiered": "AUC_frac_rate"}


def _predictor_endpoint_pairs(
    primary_endpoint: str = PRIMARY_ENDPOINT, sensitivity_endpoint: str = SENSITIVITY_ENDPOINT
) -> list[tuple[str, str]]:
    pairs = [(p, e) for p in PREDICTORS for e in (primary_endpoint, sensitivity_endpoint)]
    pairs.append(("n_k", SECONDARY_ENDPOINT))  # confound-visibility pairing only
    return pairs


def _clean_for_correlation(df: pd.DataFrame, predictor: str, endpoint: str) -> pd.DataFrame:
    sub = df[[predictor, endpoint]].copy()
    if predictor == "r_val":
        sub = sub[sub["r_val"] != UNREACHED]
    sub = sub.dropna(subset=[predictor, endpoint])
    return sub


def _bootstrap_tau_ci(
    x: np.ndarray, y: np.ndarray, *, n_boot: int = 10_000, seed: int = 0
) -> tuple[float, float, int]:
    """10,000 resamples over instances (rows), seed 0, explicit Generator.
    Returns (lo2.5, hi97.5, n_nan_resamples)."""
    rng = np.random.default_rng(seed)
    n = len(x)
    taus = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        t, _ = kendalltau(x[idx], y[idx], variant="b", nan_policy="propagate")
        taus[i] = t
    n_nan = int(np.isnan(taus).sum())
    valid = taus[~np.isnan(taus)]
    if len(valid) == 0:
        return float("nan"), float("nan"), n_nan
    lo, hi = np.percentile(valid, [2.5, 97.5])
    return float(lo), float(hi), n_nan


def tau_for_stratum(
    df: pd.DataFrame, predictor: str, endpoint: str, *, r_assumes: str,
    n_boot: int = 10_000, seed: int = 0,
) -> dict[str, Any]:
    """One stratum x one predictor x one endpoint."""
    total_n = len(df)
    n_excluded_unreached = int((df["r_val"] == UNREACHED).sum()) if predictor == "r_val" else 0
    sub = _clean_for_correlation(df, predictor, endpoint)
    n = len(sub)
    n_excluded_nan = max(total_n - n_excluded_unreached - n, 0)

    row: dict[str, Any] = {
        "predictor": predictor,
        "endpoint": endpoint,
        "stratum_n_instances": total_n,
        "n": n,
        "n_excluded_unreached": n_excluded_unreached,
        "n_excluded_nan": n_excluded_nan,
        "predictor_constant": False,
        "tau_b": float("nan"),
        "p_value": float("nan"),
        "ci_lo_2p5": float("nan"),
        "ci_hi_97p5": float("nan"),
        "n_boot_nan": 0,
        "underpowered": n < UNDERPOWERED_N,
        "status": "ok",
        "r_assumes": r_assumes if predictor == "r_val" else "",
    }

    if n < 2:
        row["status"] = "insufficient_n"
        return row

    x = sub[predictor].to_numpy(dtype=float)
    y = sub[endpoint].to_numpy(dtype=float)

    if np.unique(x).size <= 1:
        row["predictor_constant"] = True
        row["status"] = "undefined (predictor constant)"
        return row
    if np.unique(y).size <= 1:
        row["status"] = "undefined (endpoint constant)"
        return row

    tau, pval = kendalltau(x, y, variant="b", nan_policy="propagate")
    row["tau_b"] = float(tau)
    row["p_value"] = float(pval)

    lo, hi, n_boot_nan = _bootstrap_tau_ci(x, y, n_boot=n_boot, seed=seed)
    row["ci_lo_2p5"] = lo
    row["ci_hi_97p5"] = hi
    row["n_boot_nan"] = n_boot_nan
    return row


def _stratified_tau_arm(
    arm_df: pd.DataFrame, arm_name: str, strata_cols: list[str], *, r_assumes: str,
    n_boot: int = 10_000, seed: int = 0, primary_endpoint: str = PRIMARY_ENDPOINT,
    sensitivity_endpoint: str = SENSITIVITY_ENDPOINT,
) -> pd.DataFrame:
    rows = []
    strata_keys = (
        arm_df[strata_cols]
        .drop_duplicates()
        .sort_values(strata_cols, kind="stable")
        .to_records(index=False)
    )
    for key in strata_keys:
        key = tuple(key)
        mask = np.ones(len(arm_df), dtype=bool)
        for col, val in zip(strata_cols, key):
            if isinstance(val, float) and math.isnan(val):
                mask &= arm_df[col].isna().to_numpy()
            else:
                mask &= (arm_df[col] == val).to_numpy()
        stratum_df = arm_df[mask]
        stratum_label = {col: val for col, val in zip(strata_cols, key)}
        for predictor, endpoint in _predictor_endpoint_pairs(primary_endpoint, sensitivity_endpoint):
            row = tau_for_stratum(
                stratum_df, predictor, endpoint, r_assumes=r_assumes, n_boot=n_boot, seed=seed
            )
            row["arm"] = arm_name
            row["is_primary"] = endpoint == primary_endpoint
            row["is_sensitivity"] = endpoint == sensitivity_endpoint
            row["is_secondary"] = endpoint == SECONDARY_ENDPOINT
            row.update(stratum_label)
            rows.append(row)
    return pd.DataFrame(rows)


_ALL_STRATA_COLS = sorted(
    set(FLIP_PRIMARY_STRATA_COLS)
    | set(FLIP_SECONDARY_STRATA_COLS)
    | set(TIERED_PRIMARY_STRATA_COLS)
    | set(TIERED_SECONDARY_STRATA_COLS)
)


def _order_tau_columns(df: pd.DataFrame) -> pd.DataFrame:
    front = ["arm"] + [c for c in _ALL_STRATA_COLS if c in df.columns]
    rest = [
        "predictor", "endpoint", "is_primary", "is_sensitivity", "is_secondary",
        "stratum_n_instances", "n", "n_excluded_unreached",
        "n_excluded_nan", "predictor_constant", "underpowered", "status", "tau_b",
        "p_value", "ci_lo_2p5", "ci_hi_97p5", "n_boot_nan", "r_assumes",
    ]
    cols = front + rest
    cols = [c for c in cols if c in df.columns]
    sort_cols = front + ["predictor", "endpoint"]
    return df[cols].sort_values(sort_cols, kind="stable").reset_index(drop=True)


def build_tau_table(
    instances: pd.DataFrame, *, level: str, r_assumes: str, n_boot: int = 10_000, seed: int = 0
) -> pd.DataFrame:
    """``level`` in {"primary", "secondary"}. Primary strata: flip arm by
    (coverage, base_wrongness) [6], tiered arm by n_tiers [3] -- 9 total.
    Secondary strata: each of those split further by component_size.

    Per-arm primary endpoint (:data:`PRIMARY_ENDPOINT_BY_ARM`, Appendix E):
    flip uses ``AUC_frac`` (unaffected), tiered uses the re-binned
    ``AUC_frac_rate`` (``instances`` must already carry that column, e.g. via
    :func:`attach_auc_frac_rate`). Rows are tagged ``is_primary`` /
    ``is_sensitivity`` / ``is_secondary`` so downstream verdicts select on
    that flag rather than on a single hard-coded endpoint string, since the
    two arms' primary endpoint is no longer the same column name.
    """
    if level == "primary":
        flip_cols, tiered_cols = FLIP_PRIMARY_STRATA_COLS, TIERED_PRIMARY_STRATA_COLS
    elif level == "secondary":
        flip_cols, tiered_cols = FLIP_SECONDARY_STRATA_COLS, TIERED_SECONDARY_STRATA_COLS
    else:
        raise ValueError(level)

    flip = instances[instances["arm"] == "flip"]
    tiered = instances[instances["arm"] == "tiered"]
    frames = [
        _stratified_tau_arm(
            flip, "flip", flip_cols, r_assumes=r_assumes, n_boot=n_boot, seed=seed,
            primary_endpoint=PRIMARY_ENDPOINT_BY_ARM["flip"],
        ),
        _stratified_tau_arm(
            tiered, "tiered", tiered_cols, r_assumes=r_assumes, n_boot=n_boot, seed=seed,
            primary_endpoint=PRIMARY_ENDPOINT_BY_ARM["tiered"],
        ),
    ]
    return _order_tau_columns(pd.concat(frames, ignore_index=True, sort=False))


# ---------------------------------------------------------------------------
# Verdicts P1-P4
# ---------------------------------------------------------------------------


def _primary_strata_key_cols(tau_table: pd.DataFrame) -> list[str]:
    return [
        c for c in tau_table.columns
        if c not in (
            "predictor", "endpoint", "is_primary", "is_sensitivity", "is_secondary",
            "stratum_n_instances", "n", "n_excluded_unreached",
            "n_excluded_nan", "predictor_constant", "underpowered", "status", "tau_b",
            "p_value", "ci_lo_2p5", "ci_hi_97p5", "n_boot_nan", "r_assumes",
        )
    ]


def _endpoint_mask(tau_table: pd.DataFrame, *, use_sensitivity: bool) -> pd.Series:
    """Selects rows on the ``is_primary`` / ``is_sensitivity`` flag rather
    than a literal endpoint-string match, because per Appendix E the two arms'
    primary endpoint is no longer the same column (flip: ``AUC_frac``;
    tiered: ``AUC_frac_rate``) -- see :data:`PRIMARY_ENDPOINT_BY_ARM`."""
    col = "is_sensitivity" if use_sensitivity else "is_primary"
    return tau_table[col] if col in tau_table.columns else (tau_table["endpoint"] == PRIMARY_ENDPOINT)


def _endpoint_label(tau_table: pd.DataFrame, mask: pd.Series) -> str:
    vals = sorted(tau_table.loc[mask, "endpoint"].unique()) if mask.any() else []
    return " / ".join(vals) if vals else "(none)"


def verdict_p1(tau_table: pd.DataFrame, *, use_sensitivity: bool = False) -> dict[str, Any]:
    """P1: tau_b(AUC_frac, r_val) > 0 in a majority of (defined) strata.
    Selects each arm's own primary (or sensitivity) endpoint via the
    ``is_primary`` / ``is_sensitivity`` flag -- flip and tiered use different
    endpoint columns post-Appendix-E (see :data:`PRIMARY_ENDPOINT_BY_ARM`)."""
    mask = _endpoint_mask(tau_table, use_sensitivity=use_sensitivity)
    all_rval = tau_table[(tau_table["predictor"] == "r_val") & mask]
    sub = all_rval[~all_rval["tau_b"].isna()]
    n_defined = len(sub)
    n_pos = int((sub["tau_b"] > 0).sum())
    supported = n_defined > 0 and n_pos > n_defined / 2
    verdict = "SUPPORTED" if supported else ("UNDERPOWERED" if n_defined == 0 else "NOT SUPPORTED")
    return {
        "prediction": "P1",
        "endpoint": _endpoint_label(tau_table, mask),
        "statement": "tau_b(AUC_frac, r_val) > 0 in a majority of strata",
        "n_strata_total": len(all_rval),
        "n_strata_defined": n_defined,
        "n_positive": n_pos,
        "verdict": verdict,
        "detail": sub[_primary_strata_key_cols(tau_table) + ["endpoint", "tau_b", "ci_lo_2p5", "ci_hi_97p5", "n", "underpowered"]],
    }


def verdict_p2(tau_table: pd.DataFrame, *, use_sensitivity: bool = False) -> dict[str, Any]:
    """P2: tau_b(AUC_frac, r_val) > tau_b(AUC_frac, n_k) in a majority of
    strata where both are defined (same per-arm endpoint selection as P1)."""
    mask = _endpoint_mask(tau_table, use_sensitivity=use_sensitivity)
    key_cols = _primary_strata_key_cols(tau_table)
    rv = tau_table[(tau_table["predictor"] == "r_val") & mask]
    nk = tau_table[(tau_table["predictor"] == "n_k") & mask]
    merged = rv.merge(nk, on=key_cols, suffixes=("_rval", "_nk"))
    both_defined = merged[(~merged["tau_b_rval"].isna()) & (~merged["tau_b_nk"].isna())]
    n_defined = len(both_defined)
    n_rval_wins = int((both_defined["tau_b_rval"] > both_defined["tau_b_nk"]).sum())
    supported = n_defined > 0 and n_rval_wins > n_defined / 2
    verdict = "SUPPORTED" if supported else ("UNDERPOWERED" if n_defined == 0 else "NOT SUPPORTED")
    return {
        "prediction": "P2",
        "endpoint": _endpoint_label(tau_table, mask),
        "statement": "tau_b(AUC_frac, r_val) > tau_b(AUC_frac, n_k) where both defined",
        "n_strata_both_defined": n_defined,
        "n_rval_beats_nk": n_rval_wins,
        "verdict": verdict,
        "detail": both_defined[key_cols + ["endpoint_rval", "tau_b_rval", "tau_b_nk"]],
    }


def verdict_p3(tau_table: pd.DataFrame, *, use_sensitivity: bool = False) -> dict[str, Any]:
    """P3 / H4: tau_b(AUC_frac, n_k) CI covers 0 in a majority of strata
    (same per-arm endpoint selection as P1)."""
    mask = _endpoint_mask(tau_table, use_sensitivity=use_sensitivity)
    all_nk = tau_table[(tau_table["predictor"] == "n_k") & mask]
    sub = all_nk[~all_nk["tau_b"].isna()]
    n_defined = len(sub)
    covers0 = (sub["ci_lo_2p5"] <= 0) & (sub["ci_hi_97p5"] >= 0)
    n_covers = int(covers0.sum())
    supported = n_defined > 0 and n_covers > n_defined / 2
    verdict = "SUPPORTED" if supported else ("UNDERPOWERED" if n_defined == 0 else "NOT SUPPORTED")
    return {
        "prediction": "P3",
        "endpoint": _endpoint_label(tau_table, mask),
        "statement": "tau_b(AUC_frac, n_k) CI covers 0 in a majority of strata (H4)",
        "n_strata_defined": n_defined,
        "n_ci_covers_zero": n_covers,
        "verdict": verdict,
        "detail": sub[_primary_strata_key_cols(tau_table) + ["endpoint", "tau_b", "ci_lo_2p5", "ci_hi_97p5", "n"]],
    }


def verdict_p4(instances: pd.DataFrame) -> dict[str, Any]:
    """P4, per PREREGISTRATION.md Appendix E.4: the tiered arm's corrected
    primary endpoint (``AUC_frac_rate``, gridded on corruption_rate = fraction
    of *nodes* relocated) and the flip arm's (``AUC_frac``, gridded on
    d/n_k = fraction of *claims* reversed) are in different units. This
    function therefore reports both arms' across-instance variance --
    overall and on matched (component_size, separation) cells -- but does
    **not** force a SUPPORTED/NOT SUPPORTED verdict on the comparison: per
    E.4, "not comparable as pre-registered" is the expected and reported
    answer, not a null rounded into a pass."""
    flip = instances[instances["arm"] == "flip"].dropna(subset=["AUC_frac"])
    tiered = instances[instances["arm"] == "tiered"].dropna(subset=["AUC_frac_rate"])
    var_flip_overall = float(flip["AUC_frac"].var(ddof=1))
    var_tiered_overall = float(tiered["AUC_frac_rate"].var(ddof=1))

    flip_cells = set(map(tuple, flip[["component_size", "separation"]].drop_duplicates().to_numpy()))
    tiered_cells = set(map(tuple, tiered[["component_size", "separation"]].drop_duplicates().to_numpy()))
    matched = sorted(flip_cells & tiered_cells)

    rows = []
    for cs, sep in matched:
        f = flip[(flip["component_size"] == cs) & (flip["separation"] == sep)]["AUC_frac"]
        t = tiered[(tiered["component_size"] == cs) & (tiered["separation"] == sep)]["AUC_frac_rate"]
        if len(f) >= 2 and len(t) >= 2:
            rows.append(
                {
                    "component_size": cs, "separation": sep,
                    "n_flip": len(f), "n_tiered": len(t),
                    "var_flip_AUC_frac": float(f.var(ddof=1)),
                    "var_tiered_AUC_frac_rate": float(t.var(ddof=1)),
                    "tiered_higher": bool(t.var(ddof=1) > f.var(ddof=1)),
                }
            )
    matched_df = pd.DataFrame(rows)
    n_matched = len(matched_df)
    n_tiered_higher = int(matched_df["tiered_higher"].sum()) if n_matched else 0

    return {
        "prediction": "P4",
        "endpoint": "flip: AUC_frac (d/n_k axis) vs tiered: AUC_frac_rate (corruption_rate axis)",
        "statement": "across-instance variance of AUC_frac higher for tiered than flip",
        "var_flip_overall_AUC_frac": var_flip_overall,
        "var_tiered_overall_AUC_frac_rate": var_tiered_overall,
        "overall_tiered_higher_raw": bool(var_tiered_overall > var_flip_overall),
        "n_matched_cells": n_matched,
        "n_cells_tiered_higher_raw": n_tiered_higher,
        "unit_mismatch": (
            "flip AUC_frac is gridded on d/n_k (fraction of claims reversed); "
            "tiered AUC_frac_rate is gridded on corruption_rate (fraction of nodes "
            "relocated). Different axes, different units -- per Appendix E.4 the "
            "raw variance numbers above are reported for completeness only and are "
            "NOT a same-units comparison."
        ),
        "verdict": "NOT COMPARABLE AS PRE-REGISTERED (unit mismatch, Appendix E.4)",
        "matched_detail": matched_df,
    }


# ---------------------------------------------------------------------------
# Cross-check against hops_p5.csv (P5 subsample)
# ---------------------------------------------------------------------------


def hops_cross_check(tau_primary: pd.DataFrame, hops_p5: pd.DataFrame) -> pd.DataFrame:
    """Compares this dataset's flip-arm, coverage x base_wrongness tau_b(AUC_frac,
    r_val) against ``hops_p5.csv``'s precomputed tau(AUC_claims, r_val) on the
    independent small-component (component_size 3-7) subsample, for the same
    6 strata. Reports sign agreement and whether the magnitude is in rough
    range (main tau within +/-0.35 of the hops tau, generously wide because
    the subsample and the main population have disjoint component sizes)."""
    ours = tau_primary[
        (tau_primary["arm"] == "flip")
        & (tau_primary["predictor"] == "r_val")
        & (tau_primary["endpoint"] == PRIMARY_ENDPOINT)
    ][["coverage", "base_wrongness", "n", "tau_b", "ci_lo_2p5", "ci_hi_97p5"]].copy()
    ours = ours.rename(columns={"n": "n_main", "tau_b": "tau_main"})

    hops = hops_p5[hops_p5["stratum"] != "ALL_POOLED"][
        ["coverage", "base_wrongness", "n_instances", "tau_claims"]
    ].copy()
    hops = hops.rename(columns={"n_instances": "n_hops", "tau_claims": "tau_hops"})

    merged = ours.merge(hops, on=["coverage", "base_wrongness"], how="outer")
    merged["sign_agree"] = np.sign(merged["tau_main"]) == np.sign(merged["tau_hops"])
    merged["abs_diff"] = (merged["tau_main"] - merged["tau_hops"]).abs()
    merged["rough_magnitude_agree"] = merged["abs_diff"] <= 0.35
    return merged.sort_values(["coverage", "base_wrongness"], kind="stable").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Task 2 -- discordance spotlight
# ---------------------------------------------------------------------------

DISCORDANCE_COLS = [
    "discordance_type", "discordance_score", "instance_id", "arm", "component_size",
    "separation", "coverage", "base_wrongness", "n_tiers",
    "r_val", "r_status", "r_assumes", "shd_truth", "shd_cpdag", "n_k",
    "undirected_fraction", "AUC_frac", "AUC_frac_rate", "AUC_frac_usable", "AUC_abs",
    "baseline_pct", "auc_pct", "S_series_d_val", "S_contra_as_fail_series_d_val",
]


def _pct_rank(s: pd.Series) -> pd.Series:
    return s.rank(method="average", pct=True)


def discordance_spotlight(
    instances: pd.DataFrame, curves: pd.DataFrame, *, top_n: int = 15
) -> pd.DataFrame:
    """Task 2. Ranks instances by how much shd_truth/n_k (naive baselines)
    disagree with each arm's own primary AUC endpoint (the actual survival
    outcome), separately within each arm (flip and tiered are different
    instance populations per PREREGISTRATION.md 5.1, and per Appendix E the
    tiered arm's primary endpoint is the re-binned ``AUC_frac_rate``, not the
    defective ``AUC_frac`` -- see :data:`PRIMARY_ENDPOINT_BY_ARM`), then pools
    the most extreme cases.

    discordance score = pct_rank(shd_truth) averaged with pct_rank(n_k)
    ("baseline_pct", low = looks close to truth / few claims), combined with
    pct_rank(arm's primary AUC endpoint) ("auc_pct", low = collapses fast).

    low_low: baseline_pct + auc_pct small (looks fine by naive baselines,
        fails fast) -- the case where naive baselines mislead by false
        reassurance.
    high_high: baseline_pct + auc_pct large (looks bad by naive baselines,
        survives deep corruption) -- the case where naive baselines mislead
        by false alarm.
    """
    df = instances.copy()
    df["baseline_pct"] = np.nan
    df["auc_pct"] = np.nan
    for arm_name, g in df.groupby("arm"):
        idx = g.index
        shd_pct = _pct_rank(g["shd_truth"])
        nk_pct = _pct_rank(g["n_k"])
        df.loc[idx, "baseline_pct"] = (shd_pct + nk_pct) / 2.0
        arm_endpoint = PRIMARY_ENDPOINT_BY_ARM.get(arm_name, PRIMARY_ENDPOINT)
        df.loc[idx, "auc_pct"] = _pct_rank(g[arm_endpoint])

    df["combo"] = df["baseline_pct"] + df["auc_pct"]

    low_low = df.sort_values(["combo", "instance_id"], kind="stable").head(top_n).copy()
    low_low["discordance_type"] = "low_low (looks near truth, collapses fast)"
    low_low["discordance_score"] = low_low["combo"]

    high_high = df.sort_values(["combo", "instance_id"], ascending=[False, True], kind="stable").head(top_n).copy()
    high_high["discordance_type"] = "high_high (looks far from truth, survives deep corruption)"
    high_high["discordance_score"] = high_high["combo"]

    spotlight = pd.concat([low_low, high_high], ignore_index=True)

    series_by_id: dict[str, str] = {}
    series_contra_by_id: dict[str, str] = {}
    curves_sorted = curves.sort_values(["instance_id", "d"], kind="stable")
    relevant = curves_sorted[curves_sorted["instance_id"].isin(spotlight["instance_id"])]
    for inst_id, g in relevant.groupby("instance_id"):
        g = g.sort_values("d", kind="stable")
        series_by_id[inst_id] = ";".join(
            f"{int(d)}:{s:.4f}" if not math.isnan(s) else f"{int(d)}:nan"
            for d, s in zip(g["d"], g["S"])
        )
        series_contra_by_id[inst_id] = ";".join(
            f"{int(d)}:{s:.4f}" for d, s in zip(g["d"], g["S_contra_as_fail"])
        )

    spotlight["S_series_d_val"] = spotlight["instance_id"].map(series_by_id)
    spotlight["S_contra_as_fail_series_d_val"] = spotlight["instance_id"].map(series_contra_by_id)

    cols = [c for c in DISCORDANCE_COLS if c in spotlight.columns]
    return spotlight[cols].sort_values(["discordance_type", "discordance_score"], kind="stable").reset_index(drop=True)


def verify_spotlight_cases(
    spotlight: pd.DataFrame, instances: pd.DataFrame, curves: pd.DataFrame
) -> pd.DataFrame:
    """E5. For every spotlighted instance: (a) independently recomputes
    AUC_frac straight from survival_curves.csv via the live
    :func:`bkrobust.robustness.survival.auc_frac` and checks it reproduces
    the value already in survival_instances.csv (data-integrity check, not a
    narrative heuristic) and (b) reports the first- and last-evaluated S so a
    human can eyeball whether the case actually looks like what its
    discordance_type claims."""
    nk_map = dict(zip(instances["instance_id"], instances["n_k"]))
    rows = []
    for _, r in spotlight.iterrows():
        iid = r["instance_id"]
        g = curves[curves["instance_id"] == iid].sort_values("d", kind="stable")
        curve = {
            int(row["d"]): {"S": float(row["S"])}
            for _, row in g.iterrows() if pd.notna(row["S"])
        }
        nk = nk_map.get(iid)
        recomputed = _auc_frac_fn(curve, int(nk)) if nk else ""
        recomputed = float(recomputed) if recomputed != "" else float("nan")
        recorded = float(r["AUC_frac"])
        matches = (not math.isnan(recomputed)) and abs(recomputed - recorded) < 1e-6

        s_vals = g["S"].dropna()
        first_s = float(s_vals.iloc[0]) if len(s_vals) else float("nan")
        last_s = float(s_vals.iloc[-1]) if len(s_vals) else float("nan")
        rows.append(
            {
                "instance_id": iid,
                "discordance_type": r["discordance_type"],
                "recorded_AUC_frac": recorded,
                "recomputed_AUC_frac": recomputed,
                "auc_matches": bool(matches),
                "first_S": first_s,
                "last_S": last_s,
                "n_curve_rows": len(g),
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Task 3 -- effect sizes (median AUC_frac by predictor bucket)
# ---------------------------------------------------------------------------

EFFECT_SIZE_PREDICTORS = ["r_val", "shd_truth", "n_k"]


def _integer_bucket(v: float, cap: int = 5) -> str:
    if pd.isna(v):
        return "NA"
    v = int(v)
    if v >= cap:
        return f"{cap}+"
    return str(v)


def effect_size_table(
    instances: pd.DataFrame, arm: str, strata_cols: list[str], predictor: str,
    *, endpoint: str = PRIMARY_ENDPOINT,
) -> pd.DataFrame:
    """Median (and mean, n) of ``endpoint`` by ``predictor`` bucket, within
    each stratum of the given arm's PRIMARY stratification."""
    df = instances[instances["arm"] == arm].copy()
    if predictor == "r_val":
        df = df[df["r_val"] != UNREACHED]
    df = df.dropna(subset=[predictor, endpoint])
    df["bucket"] = df[predictor].map(_integer_bucket)

    rows = []
    strata_keys = df[strata_cols].drop_duplicates().sort_values(strata_cols, kind="stable").to_records(index=False)
    for key in strata_keys:
        key = tuple(key)
        mask = np.ones(len(df), dtype=bool)
        for col, val in zip(strata_cols, key):
            if isinstance(val, float) and math.isnan(val):
                mask &= df[col].isna().to_numpy()
            else:
                mask &= (df[col] == val).to_numpy()
        stratum_df = df[mask]
        for bucket, bg in stratum_df.groupby("bucket", observed=True):
            if len(bg) == 0:
                continue
            row = {"arm": arm}
            row.update({col: val for col, val in zip(strata_cols, key)})
            row.update(
                {
                    "predictor": predictor,
                    "bucket": str(bucket),
                    "n": len(bg),
                    f"median_{endpoint}": float(bg[endpoint].median()),
                    f"mean_{endpoint}": float(bg[endpoint].mean()),
                }
            )
            rows.append(row)
    out = pd.DataFrame(rows)
    if len(out) == 0:
        return out
    sort_cols = ["arm"] + strata_cols + ["predictor", "bucket"]
    return out.sort_values(sort_cols, kind="stable").reset_index(drop=True)


def build_effect_size_tables(instances: pd.DataFrame) -> pd.DataFrame:
    """All of Task 3: for each primary stratum (flip: coverage x
    base_wrongness; tiered: n_tiers) and each of r_val / shd_truth / n_k,
    median/mean of the arm's own primary AUC endpoint by bucket. Per
    Appendix E the tiered arm's primary endpoint is the re-binned
    ``AUC_frac_rate`` (:data:`PRIMARY_ENDPOINT_BY_ARM`), not the defective
    ``AUC_frac``; the response column is normalized to ``median_endpoint_value``
    / ``mean_endpoint_value`` with an ``endpoint`` column stating which raw
    column was used per arm, so flip and tiered rows sit in the same columns
    despite using different response axes."""
    frames = []
    for predictor in EFFECT_SIZE_PREDICTORS:
        for arm, strata_cols in (
            ("flip", FLIP_PRIMARY_STRATA_COLS), ("tiered", TIERED_PRIMARY_STRATA_COLS)
        ):
            arm_endpoint = PRIMARY_ENDPOINT_BY_ARM[arm]
            tbl = effect_size_table(instances, arm, strata_cols, predictor, endpoint=arm_endpoint)
            if len(tbl):
                tbl = tbl.rename(
                    columns={
                        f"median_{arm_endpoint}": "median_endpoint_value",
                        f"mean_{arm_endpoint}": "mean_endpoint_value",
                    }
                )
                tbl["endpoint"] = arm_endpoint
            frames.append(tbl)
    frames = [f for f in frames if len(f)]
    combined = pd.concat(frames, ignore_index=True, sort=False)
    front = (
        ["arm"] + [c for c in _ALL_STRATA_COLS if c in combined.columns]
        + ["predictor", "bucket", "endpoint"]
    )
    rest = [c for c in combined.columns if c not in front]
    return combined[front + rest]
