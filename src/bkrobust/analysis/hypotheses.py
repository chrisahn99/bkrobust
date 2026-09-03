"""The pre-registered statistics for H1-H6, computed from committed result rows.

Each function implements exactly the statistic named in
``results/synth/preregistration.md`` and returns both the number and the
evidence needed to judge it. Nothing here decides whether a hypothesis
"passed" -- verdicts are the orchestrator's, written in the report, so that a
threshold cannot be quietly moved after seeing the data.

Two rules the pre-registration fixed, enforced here:

* ``UNREACHED`` (-1) is never averaged or compared as a number. Rows carrying it
  are counted separately and excluded from distributions.
* Distributions are reported as counts, not means.
"""

from __future__ import annotations

from typing import Any

import pandas as pd

from bkrobust.core.conventions import UNREACHED


def _finite(series: pd.Series) -> pd.Series:
    """Drop the UNREACHED sentinel rather than letting it act as a number."""
    return series[series != UNREACHED]


def h1_saturation(df: pd.DataFrame, by: str | None = None) -> dict[str, Any]:
    """H1: the full distribution of ``r_val``, and how much mass sits at 1.

    Args:
        df: Instance rows.
        by: Optional grouping column, e.g. ``"generator"``.

    Returns:
        Counts per radius, the fraction at 1, and the UNREACHED count. Never a
        bare mean -- the pre-registration asks for the distribution.
    """
    out: dict[str, Any] = {}
    groups = [(None, df)] if by is None else list(df.groupby(by))
    for key, g in groups:
        fin = _finite(g["r_val"])
        hist = fin.value_counts().sort_index().to_dict()
        out[str(key)] = {
            "n_rows": len(g),
            "n_finite": len(fin),
            "n_unreached": int((g["r_val"] == UNREACHED).sum()),
            "hist": {int(k): int(v) for k, v in hist.items()},
            "frac_r_eq_1": float((fin == 1).mean()) if len(fin) else float("nan"),
            "frac_r_ge_2": float((fin >= 2).mean()) if len(fin) else float("nan"),
        }
    return out


def h2_frontier(df: pd.DataFrame, by: str | None = None) -> dict[str, Any]:
    """H2: how often some valid adjustment set is strictly more robust than ``O``.

    Expects columns ``r_val`` (for ``O``) and ``best_r_val`` (max over valid
    sets). Rows where either is UNREACHED are reported separately rather than
    compared numerically.
    """
    out: dict[str, Any] = {}
    groups = [(None, df)] if by is None else list(df.groupby(by))
    for key, g in groups:
        both = g[(g["r_val"] != UNREACHED) & (g["best_r_val"] != UNREACHED)]
        gap = both["best_r_val"] - both["r_val"]
        out[str(key)] = {
            "n_comparable": len(both),
            "n_strictly_better": int((gap > 0).sum()),
            "frac_strictly_better": float((gap > 0).mean()) if len(both) else float("nan"),
            "gap_hist": {int(k): int(v) for k, v in gap.value_counts().sort_index().items()},
            "n_excluded_unreached": int(len(g) - len(both)),
        }
    return out


def h3_coupling(df: pd.DataFrame, coupling_col: str = "param_coupling") -> dict[str, Any]:
    """H3: the frontier statistic as a function of the coupling parameter.

    The key designed test. If the fraction with a strictly more robust set does
    NOT rise as coupling falls, the prior report's structural explanation of the
    flat frontier is wrong or incomplete -- which is a more important finding
    than a confirmation.
    """
    if coupling_col not in df.columns:
        return {"error": f"missing column {coupling_col}"}
    out: dict[str, Any] = {}
    for c, g in df.groupby(coupling_col):
        both = g[(g["r_val"] != UNREACHED) & (g["best_r_val"] != UNREACHED)]
        gap = both["best_r_val"] - both["r_val"]
        out[str(c)] = {
            "n": len(both),
            "frac_strictly_better": float((gap > 0).mean()) if len(both) else float("nan"),
            "mean_gap": float(gap.mean()) if len(both) else float("nan"),
        }
    return out


def h4_naive(df: pd.DataFrame) -> dict[str, Any]:
    """H4: naive K-count radius against the model radius.

    Tests whether ``naive <= model`` is universal, which the prior single
    example suggested but could not establish.
    """
    if "naive_radius" not in df.columns:
        return {"error": "missing column naive_radius"}
    both = df[(df["r_val"] != UNREACHED) & (df["naive_radius"] != UNREACHED)]
    diff = both["r_val"] - both["naive_radius"]
    return {
        "n_comparable": len(both),
        "n_naive_le_model": int((diff >= 0).sum()),
        "n_naive_gt_model": int((diff < 0).sum()),
        "frac_naive_gt_model": float((diff < 0).mean()) if len(both) else float("nan"),
        "diff_hist": {int(k): int(v) for k, v in diff.value_counts().sort_index().items()},
    }


def h5_regimes(df: pd.DataFrame, eps_col: str = "r_eps_eps_medium") -> dict[str, Any]:
    """H5: how often the middle regime occurs -- ``r_val`` low but ``r_eps`` high.

    That regime is fragile identification with a stable estimate, and it decides
    whether ``r_eps`` deserves to carry part of the contribution.
    """
    if eps_col not in df.columns:
        return {"error": f"missing column {eps_col}"}
    both = df[(df["r_val"] != UNREACHED) & (df[eps_col] != UNREACHED)]
    gap = both[eps_col] - both["r_val"]
    return {
        "n_comparable": len(both),
        "n_middle_regime": int((gap > 0).sum()),
        "frac_middle_regime": float((gap > 0).mean()) if len(both) else float("nan"),
        "gap_hist": {int(k): int(v) for k, v in gap.value_counts().sort_index().items()},
    }


def h6_calibration(df: pd.DataFrame) -> dict[str, Any]:
    """H6: coverage and conservativeness, and which descriptors predict them.

    Coverage below 1.0 is a **bug**, not a finding: the caller must halt and
    investigate rather than report it.
    """
    out: dict[str, Any] = {}
    if "coverage" in df.columns:
        cov = df["coverage"].dropna()
        out["coverage_min"] = float(cov.min()) if len(cov) else float("nan")
        out["coverage_all_one"] = bool((cov >= 1.0 - 1e-12).all()) if len(cov) else None
        out["n_coverage_below_one"] = int((cov < 1.0 - 1e-12).sum()) if len(cov) else 0
    if "conservativeness" in df.columns:
        con = df["conservativeness"].dropna()
        out["conservativeness"] = {
            "n": len(con),
            "mean": float(con.mean()) if len(con) else float("nan"),
            "min": float(con.min()) if len(con) else float("nan"),
            "max": float(con.max()) if len(con) else float("nan"),
        }
        desc_cols = [c for c in df.columns if c.startswith("desc_")]
        corrs = {}
        for c in desc_cols:
            sub = df[[c, "conservativeness"]].dropna()
            if len(sub) > 3 and sub[c].nunique() > 1:
                corrs[c] = float(sub[c].corr(sub["conservativeness"], method="spearman"))
        out["conservativeness_spearman_by_descriptor"] = dict(
            sorted(corrs.items(), key=lambda kv: -abs(kv[1]))
        )
    return out
