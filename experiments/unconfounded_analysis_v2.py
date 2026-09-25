"""Analysis for the v2 "unconfounded" tiered-arm sweep
(``results/axis_robustness_unconfounded_v2/``).

Read-only against ``survival_instances.csv`` in that directory (one row per
admitted instance; ``AUC_frac`` already computed by the driver). Never
regenerates data, never touches the sweep's own files.

For each stratum (``n_tiers in {2, 3, 4}``, plus pooled across all three):

1. **Marginal tau_b + 95% CI** (instance bootstrap, 10,000 resamples,
   ``np.random.default_rng(0)``, matching
   ``build_tau_comparisons.bootstrap_tau_ci``'s exact scheme) of ``r_val``,
   ``separation`` (restricted to rows with ``separation_status ==
   "measured"``; ``n`` reported alongside so the narrower support is never
   silently implied), ``n_k``, ``k_g0``, ``shd_truth`` against ``AUC_frac``.
2. **Paired delta** (instance-clustered bootstrap, 10,000 resamples) of
   ``r_val`` against each of the same baselines plus ``r_claim``, on rows
   matched for both predictors -- reusing
   :func:`bkrobust.robustness.paired_resample.paired_comparison` verbatim
   (never reimplemented). The cluster column is ``instance_id`` itself: this
   sweep's ``survival_instances.csv`` already carries exactly one row per
   instance, so "instance-clustered" and "row-level" bootstrap coincide here
   (unlike the real corpus, which clusters several units per network).
3. The share of instances with ``r_val != separation`` (denominator: rows
   with ``separation_status == "measured"`` AND a defined, non-sentinel
   ``r_val``) and the share with ``separation`` undefined (denominator: every
   instance in the stratum).
4. The distribution of ``r_val`` (count, mean, std, min, quartiles, max, and
   a value-count histogram, ``UNREACHED`` reported as its own bucket, never
   averaged in).
5. tau_b(r_val, AUC_frac) restricted to the separation-undefined subset
   specifically (the rows the task brief calls out as the ones where r_val is
   the only defined structural predictor).

Writes ``marginal_tau.csv``, ``paired_delta.csv``,
``r_val_distribution.csv``, and prints a human-readable summary to stdout;
the caller is expected to fold key numbers into ``FINDINGS.md`` by hand (this
script does not write prose).

Usage::

    PYTHONPATH=src .venv/bin/python experiments/unconfounded_analysis_v2.py \\
        --results-dir results/axis_robustness_unconfounded_v2
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

from bkrobust.core.conventions import UNREACHED
from bkrobust.robustness.paired_resample import DEFAULT_N_BOOT, paired_comparison

ENDPOINT = "AUC_frac"
MARGINAL_PREDICTORS = ["r_val", "separation", "n_k", "k_g0", "shd_truth"]
PAIRED_BASELINES = ["separation", "n_k", "k_g0", "shd_truth", "r_claim"]
N_BOOT = 10_000
BOOT_SEED = 0

#: Sentinels per predictor: r_val's UNREACHED (-1, a defined non-failure
#: outcome, never a number to correlate) is excluded from every tau/paired
#: computation, matching build_tau_comparisons.py's own convention.
SENTINELS: dict[str, tuple[Any, ...]] = {"r_val": (UNREACHED,)}
CENSORING_COLUMNS: dict[str, str] = {"r_claim": "r_claim_censored"}


def _to_float(v: Any) -> float | None:
    if v is None:
        return None
    if isinstance(v, str):
        s = v.strip()
        if s == "":
            return None
        try:
            v = float(s)
        except ValueError:
            return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    if np.isnan(f):
        return None
    return f


def bootstrap_tau_ci(
    x: np.ndarray, y: np.ndarray, *, n_boot: int = N_BOOT, seed: int = BOOT_SEED
) -> tuple[float | None, float | None, float | None, int]:
    """Instance bootstrap tau_b CI. Returns (tau_point, lo2.5, hi97.5, n_nan)."""
    tau_point, _ = kendalltau(x, y, variant="b", nan_policy="propagate")
    if np.isnan(tau_point):
        return None, None, None, n_boot
    rng = np.random.default_rng(seed)
    n = len(x)
    taus = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        t, _ = kendalltau(x[idx], y[idx], variant="b", nan_policy="propagate")
        taus[i] = t
    n_nan = int(np.isnan(taus).sum())
    valid = taus[~np.isnan(taus)]
    if len(valid) < 100:
        return float(tau_point), None, None, n_nan
    lo, hi = np.percentile(valid, [2.5, 97.5])
    return float(tau_point), float(lo), float(hi), n_nan


def marginal_row(
    df: pd.DataFrame, predictor: str, stratum: str, *, sentinels: dict[str, tuple[Any, ...]]
) -> dict[str, Any]:
    stratum_n = len(df)
    pred = df[predictor].map(_to_float)
    end = df[ENDPOINT].map(_to_float)
    sent = tuple(_to_float(s) for s in sentinels.get(predictor, ()))
    mask = pred.notna() & end.notna() & ~pred.isin(sent)
    x = pred[mask].to_numpy(dtype=float)
    y = end[mask].to_numpy(dtype=float)
    n = len(x)
    if n < 2:
        return {
            "stratum": stratum, "predictor": predictor, "endpoint": ENDPOINT,
            "stratum_n": stratum_n, "n": n, "tau_b": None, "ci_lo_2p5": None,
            "ci_hi_97p5": None, "n_boot": N_BOOT, "n_boot_nan": None, "status": "insufficient_n",
        }
    tau, lo, hi, n_nan = bootstrap_tau_ci(x, y)
    status = "ok" if tau is not None else "undefined (constant predictor or endpoint)"
    return {
        "stratum": stratum, "predictor": predictor, "endpoint": ENDPOINT,
        "stratum_n": stratum_n, "n": n, "tau_b": tau, "ci_lo_2p5": lo,
        "ci_hi_97p5": hi, "n_boot": N_BOOT, "n_boot_nan": n_nan, "status": status,
    }


def paired_rows(df: pd.DataFrame, stratum: str) -> list[dict[str, Any]]:
    rows = df.to_dict(orient="records")
    out = []
    for baseline in PAIRED_BASELINES:
        res = paired_comparison(
            rows, "r_val", baseline, ENDPOINT, "instance_id",
            n_boot=DEFAULT_N_BOOT, seed=BOOT_SEED,
            sentinels=SENTINELS, censoring_columns=CENSORING_COLUMNS,
        )
        res["stratum"] = stratum
        res.pop("loo_by_cluster", None)
        out.append(res)
    return out


def r_val_distribution(df: pd.DataFrame, stratum: str) -> dict[str, Any]:
    vals = df["r_val"].map(_to_float)
    n_total = len(df)
    n_unreached = int((vals == float(UNREACHED)).sum())
    finite = vals[(vals.notna()) & (vals != float(UNREACHED))]
    n_missing_or_timeout = n_total - n_unreached - len(finite)
    hist: dict[str, int] = {}
    for v in finite:
        key = str(int(v)) if float(v).is_integer() else str(v)
        hist[key] = hist.get(key, 0) + 1
    out = {
        "stratum": stratum, "n_total": n_total, "n_unreached": n_unreached,
        "n_missing_or_timeout": n_missing_or_timeout, "n_finite": int(len(finite)),
        "mean": float(finite.mean()) if len(finite) else None,
        "std": float(finite.std()) if len(finite) > 1 else None,
        "min": float(finite.min()) if len(finite) else None,
        "q25": float(finite.quantile(0.25)) if len(finite) else None,
        "median": float(finite.median()) if len(finite) else None,
        "q75": float(finite.quantile(0.75)) if len(finite) else None,
        "max": float(finite.max()) if len(finite) else None,
        "histogram": json.dumps(dict(sorted(hist.items(), key=lambda kv: float(kv[0]))))
        if hist else "{}",
    }
    return out


def run(results_dir: str) -> dict[str, Any]:
    root = Path(results_dir)
    instances = pd.read_csv(root / "survival_instances.csv")
    instances["n_tiers"] = instances["n_tiers"].astype(int)

    strata: list[tuple[str, pd.DataFrame]] = [
        (f"n_tiers={nt}", g) for nt, g in instances.groupby("n_tiers", sort=True)
    ]
    strata.append(("pooled", instances))

    marginal_records: list[dict[str, Any]] = []
    paired_records: list[dict[str, Any]] = []
    dist_records: list[dict[str, Any]] = []
    coverage_notes: list[dict[str, Any]] = []
    undefined_tau_records: list[dict[str, Any]] = []

    for stratum_name, g in strata:
        for predictor in MARGINAL_PREDICTORS:
            marginal_records.append(marginal_row(g, predictor, stratum_name, sentinels=SENTINELS))
        paired_records.extend(paired_rows(g, stratum_name))
        dist_records.append(r_val_distribution(g, stratum_name))

        measured = g[g["separation_status"] == "measured"].copy()
        n_stratum = len(g)
        n_measured = len(measured)
        n_undefined = n_stratum - n_measured
        rv = measured["r_val"].map(_to_float)
        sep = measured["separation"].map(_to_float)
        mask = rv.notna() & sep.notna() & (rv != float(UNREACHED))
        n_comparable = int(mask.sum())
        n_diff = int((rv[mask] != sep[mask]).sum())
        coverage_notes.append(
            {
                "stratum": stratum_name,
                "n_instances": n_stratum,
                "n_separation_measured": n_measured,
                "n_separation_undefined": n_undefined,
                "frac_separation_undefined": n_undefined / n_stratum if n_stratum else None,
                "n_comparable_r_val_separation": n_comparable,
                "n_r_val_ne_separation": n_diff,
                "frac_r_val_ne_separation": n_diff / n_comparable if n_comparable else None,
            }
        )

        undef_df = g[g["separation_status"] != "measured"]
        row = marginal_row(undef_df, "r_val", stratum_name, sentinels=SENTINELS)
        row["subset"] = "separation_undefined"
        undefined_tau_records.append(row)

    marginal_df = pd.DataFrame(marginal_records)
    paired_df = pd.DataFrame(paired_records)
    dist_df = pd.DataFrame(dist_records)
    coverage_df = pd.DataFrame(coverage_notes)
    undefined_tau_df = pd.DataFrame(undefined_tau_records)

    marginal_df.to_csv(root / "marginal_tau.csv", index=False)
    paired_df.to_csv(root / "paired_delta.csv", index=False)
    dist_df.to_csv(root / "r_val_distribution.csv", index=False)
    coverage_df.to_csv(root / "separation_coverage.csv", index=False)
    undefined_tau_df.to_csv(root / "tau_r_val_separation_undefined.csv", index=False)

    print("=== marginal tau_b (instance bootstrap 95% CI) ===")
    print(marginal_df.to_string(index=False))
    print("\n=== paired delta: r_val - baseline (instance-clustered bootstrap) ===")
    cols = [
        "stratum", "predictor_a", "predictor_b", "n", "delta_point",
        "ci_lo_2p5", "ci_hi_97p5", "frac_delta_gt_0", "loo_verdict_flips", "status",
    ]
    print(paired_df[cols].to_string(index=False))
    print("\n=== r_val vs separation coverage ===")
    print(coverage_df.to_string(index=False))
    print("\n=== r_val distribution ===")
    print(dist_df.drop(columns=["histogram"]).to_string(index=False))
    print("\n=== tau_b(r_val) on separation-undefined subset ===")
    print(undefined_tau_df.to_string(index=False))

    return {
        "marginal": marginal_df.to_dict(orient="records"),
        "paired": paired_df.to_dict(orient="records"),
        "coverage": coverage_df.to_dict(orient="records"),
        "distribution": dist_df.to_dict(orient="records"),
        "undefined_subset_tau": undefined_tau_df.to_dict(orient="records"),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", default="results/axis_robustness_unconfounded_v2")
    args = parser.parse_args(argv)
    run(args.results_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
