"""Paired network-level resampling ([Point 6]: "evaluate differences using
paired network-level resampling rather than inferring from separate
confidence intervals").

The committed analysis (``real_analyse.py``) reports one bootstrap CI per
predictor and leaves the reader to compare two predictors' CIs by eye. That
is exactly the error the reviewer names: two overlapping (or non-overlapping)
marginal intervals say nothing rigorous about which predictor wins on the
*same* data, because the two intervals are not computed from the same
resamples and do not share the sampling noise that a fair comparison must
cancel out.

This module fixes that. Its core is a **paired cluster bootstrap**: on each
resample of the clustering unit (here, networks), it recomputes **both**
predictors' Kendall tau-b against the same endpoint on the **same resampled
rows**, and takes the *within-resample* difference
``delta = tau_b(predictor_a) - tau_b(predictor_b)``. The distribution of
``delta`` across resamples is the paired statement "A outranks B" -- its
point estimate, its percentile CI, and the fraction of resamples with
``delta > 0`` -- which two separate CIs cannot support no matter how they are
squinted at.

Generic by design
------------------
:func:`paired_cluster_bootstrap`, :func:`leave_one_cluster_out_delta` and
:func:`paired_comparison` take **column names**, not this project's specific
schema: predictor columns, an endpoint column, a cluster column, all as
strings naming keys on plain ``list[dict]`` rows. Nothing here imports
``real_analyse`` or reads anything under ``results/axis_robustness_real/`` at
module scope. That is deliberate: the task brief calls for the *same*
function to run later against a **synthetic** N=1000 dataset at
``results/axis_robustness_p6/`` (a directory this module never assumes
exists, and does not import from at load time). The real-corpus glue lives
entirely in :func:`run_on_real_corpus` and ``__main__``, which import
``real_analyse`` and ``real_baselines_p6`` lazily, inside the function body,
so importing this module never requires either the real-corpus outputs or
the (not-yet-existing) synthetic ones to be present.

Matched-sample discipline
--------------------------
A paired comparison is only honest if both predictors are evaluated on
*the same rows*. A row missing predictor A, predictor B, the endpoint, or
carrying a sentinel value (e.g. ``UNREACHED`` on ``radius``) for either is
therefore excluded from **both** tau computations in that resample -- never
imputed, and the exclusion count is reported (``n_excluded``) alongside the
matched ``n``. This is stricter than ``real_analyse.tau_for_stratum``, which
cleans each predictor independently; here that would let the two predictors
be scored on different populations within the same "paired" delta, which is
the mistake this whole module exists to avoid.

Leave-one-network-out
----------------------
:func:`leave_one_cluster_out_delta` recomputes the point-estimate delta with
each cluster's rows removed in turn, exactly mirroring
``real_analyse._leave_one_network_out`` but on the *difference*, not a single
tau. Appendix H of the committed work found only 1 of 9 strata survive
dropping any single network, and ``paths`` alone carries large leverage; the
same caution applies here, so every paired delta this module reports also
carries a LOO verdict, and :data:`paired_comparison`'s ``loo_verdict_flips``
is ``True`` whenever the sign disagrees across a LOO removal or the CI
excludes zero while the LOO range still straddles it.

    PYTHONPATH=src /usr/bin/python3 src/bkrobust/robustness/paired_resample.py
"""

from __future__ import annotations

import argparse
import json
import math
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable

import numpy as np
from scipy.stats import kendalltau

from bkrobust.core.conventions import UNREACHED

DEFAULT_N_BOOT = 10_000
DEFAULT_SEED = 0

#: Sentinel values that must never enter a tau computation, keyed by column
#: name. The default only covers ``radius`` (``UNREACHED = -1``, never a
#: number -- see ``bkrobust/core/conventions.py``), matching
#: ``real_analyse.tau_for_stratum``'s own exclusion. Callers may pass a wider
#: mapping (e.g. for a synthetic dataset with its own sentinel columns).
DEFAULT_SENTINELS: dict[str, tuple[Any, ...]] = {"radius": (UNREACHED,)}


# ---------------------------------------------------------------------------
# Generic core: predictor/endpoint/cluster columns in, paired bootstrap out
# ---------------------------------------------------------------------------


def _kendall_tau_b(x: np.ndarray, y: np.ndarray) -> float | None:
    """Kendall tau-b, or ``None`` if undefined (constant x, constant y, or
    fewer than 2 points). Never raises, never returns ``nan`` -- ``None`` is
    the sole undefined sentinel this module ever produces for a tau value.
    """
    if len(x) < 2 or np.unique(x).size <= 1 or np.unique(y).size <= 1:
        return None
    t, _ = kendalltau(x, y, variant="b", nan_policy="propagate")
    if isinstance(t, float) and math.isnan(t):
        return None
    return float(t)


def _to_float_or_none(v: Any) -> float | None:
    """Coerces a raw cell value to ``float``, or ``None`` if it is missing.

    Handles both in-memory rows (``None``, real ``float``/``int``, ``nan``)
    and CSV-round-tripped rows (every cell a ``str``). A CSV writer emits a
    missing/undefined numeric cell as the **empty string**, never as the
    literal text ``"None"`` -- ``csv.DictWriter`` and ``real_analyse.write_csv``
    both do this (a ``None`` value becomes ``""`` on write). Treating ``""``
    as anything other than missing here is exactly the bug this function
    exists to close: it previously reached ``float('')`` and raised.

    Any string that still fails to parse (e.g. a stray status word landing in
    a numeric column) is *also* treated as missing rather than raising --
    this function is a boundary between "whatever a CSV or a dict handed us"
    and "a clean float for tau", and it is never the right place to crash.

    Args:
        v: Raw cell value.

    Returns:
        A finite ``float``, or ``None`` if ``v`` is missing/blank/unparsable.
    """
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
    if math.isnan(f):
        return None
    return f


def _is_truthy_flag(v: Any) -> bool:
    """Whether a censoring-indicator cell reads as ``True``.

    CSV-robust: ``""``, ``"0"``, ``"false"``/``"False"``, ``None``, ``0``,
    ``0.0`` and ``False`` are all falsy; ``"1"``, ``"true"``/``"True"``, ``1``,
    ``1.0`` and ``True`` are all truthy. Anything else falls back to Python
    truthiness.

    Args:
        v: Raw flag cell value.

    Returns:
        The boolean reading of ``v``.
    """
    if v is None:
        return False
    if isinstance(v, str):
        return v.strip().lower() in ("1", "true", "t", "yes", "y")
    return bool(v)


def _clean_matched(
    rows: list[dict[str, Any]],
    predictor_a: str,
    predictor_b: str,
    endpoint_col: str,
    sentinels: dict[str, tuple[Any, ...]],
    censoring_columns: dict[str, str] | None = None,
) -> tuple[list[dict[str, Any]], int, int]:
    """Rows with both predictors and the endpoint defined, non-sentinel and
    (if a censoring column is configured for that predictor) non-censored.

    Every value -- ``predictor_a``, ``predictor_b`` and ``endpoint_col`` alike
    -- is passed through :func:`_to_float_or_none`, so a CSV's empty-string
    encoding of "missing" is handled identically to an in-memory ``None``:
    the row is dropped from **both** tau computations, never imputed. This is
    deliberately not special-cased to any one predictor's name (e.g.
    ``r_claim``): whatever column carries a blank cell is treated the same
    way, which is what makes this function safe to point at ``separation``
    (whose real referent is ``separation_status``) or a legitimately-absent
    ``phi_1`` without a new bug per column.

    A censoring flag is a *different* reason for exclusion than "missing":
    it means the true value exists but is right-censored, not unmeasured.
    ``censoring_columns`` lets a caller name that flag column per predictor
    (e.g. ``{"r_claim": "r_claim_censored"}``) so the two reasons are counted
    separately rather than both landing in one undifferentiated total.

    Args:
        rows: Candidate rows.
        predictor_a: First predictor column name.
        predictor_b: Second predictor column name.
        endpoint_col: Endpoint column name.
        sentinels: ``{column: (sentinel_values...)}``. Sentinel values are
            themselves coerced through :func:`_to_float_or_none` before
            comparison, so a sentinel defined as the ``int`` ``-1`` still
            matches a CSV cell that reads ``"-1"``.
        censoring_columns: Optional ``{predictor_column: flag_column}``. A row
            whose flag column reads truthy (:func:`_is_truthy_flag`) for
            either ``predictor_a`` or ``predictor_b`` is excluded and counted
            under ``n_excluded_censored`` rather than ``n_excluded_missing``,
            regardless of whether the predictor cell itself is also blank.

    Returns:
        ``(matched_rows, n_excluded_missing, n_excluded_censored)``. Rows in
        ``matched_rows`` carry ``predictor_a``, ``predictor_b`` and
        ``endpoint_col`` already coerced to ``float`` (a shallow copy of the
        input row, so the caller's original rows are never mutated).
    """
    censoring_columns = censoring_columns or {}
    cens_col_a = censoring_columns.get(predictor_a)
    cens_col_b = censoring_columns.get(predictor_b)
    sa = tuple(_to_float_or_none(s) for s in sentinels.get(predictor_a, ()))
    sb = tuple(_to_float_or_none(s) for s in sentinels.get(predictor_b, ()))

    cleaned = []
    n_excluded_missing = 0
    n_excluded_censored = 0
    for r in rows:
        if (cens_col_a and _is_truthy_flag(r.get(cens_col_a))) or (
            cens_col_b and _is_truthy_flag(r.get(cens_col_b))
        ):
            n_excluded_censored += 1
            continue
        pa = _to_float_or_none(r.get(predictor_a))
        pb = _to_float_or_none(r.get(predictor_b))
        ev = _to_float_or_none(r.get(endpoint_col))
        if pa is None or pb is None or ev is None:
            n_excluded_missing += 1
            continue
        if pa in sa or pb in sb:
            n_excluded_missing += 1
            continue
        row_copy = dict(r)
        row_copy[predictor_a] = pa
        row_copy[predictor_b] = pb
        row_copy[endpoint_col] = ev
        cleaned.append(row_copy)
    return cleaned, n_excluded_missing, n_excluded_censored


def paired_cluster_bootstrap(
    rows: list[dict[str, Any]],
    predictor_a: str,
    predictor_b: str,
    endpoint_col: str,
    cluster_col: str,
    *,
    n_boot: int = DEFAULT_N_BOOT,
    seed: int = DEFAULT_SEED,
    sentinels: dict[str, tuple[Any, ...]] | None = None,
    censoring_columns: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Paired cluster bootstrap for ``tau_b(predictor_a) - tau_b(predictor_b)``.

    Each resample draws ``len(distinct clusters)`` clusters **with
    replacement** and takes *all* rows of each drawn cluster (a cluster drawn
    twice contributes its rows twice) -- the same resampling scheme
    ``real_analyse._cluster_bootstrap_tau`` uses for a single predictor, here
    applied to both predictors on the *same* draw so their difference is
    paired.

    Args:
        rows: Rows carrying ``predictor_a``, ``predictor_b``, ``endpoint_col``
            and ``cluster_col`` keys.
        predictor_a: First predictor's column name (e.g. ``"radius"``).
        predictor_b: Second (baseline) predictor's column name.
        endpoint_col: Endpoint column name (e.g. ``"AUC_frac_usable"``).
        cluster_col: The resampling unit's column name (e.g. ``"network"``).
        n_boot: Number of resamples. The task brief requires at least 10,000;
            the default matches it and ``real_analyse.py``'s own ``--n-boot``.
        seed: Seed for ``numpy.random.default_rng``. Recorded on the output
            for reproducibility.
        sentinels: See :data:`DEFAULT_SENTINELS`; defaults to it.
        censoring_columns: Optional ``{predictor_column: flag_column}``; see
            :func:`_clean_matched`. A row flagged censored on either
            predictor is excluded and counted under ``n_excluded_censored``,
            reported separately from ``n_excluded_missing`` -- the task's
            requirement that censoring be "a denominator, not silently done."

    Returns:
        A dict: ``predictor_a``, ``predictor_b``, ``endpoint``, ``cluster_col``,
        ``stratum_n`` (rows offered), ``n`` (matched, used), ``n_excluded``
        (total dropped, for backward-compatible summaries), ``n_excluded_missing``,
        ``n_excluded_censored``, ``n_clusters``, ``tau_a_point``, ``tau_b_point``,
        ``delta_point`` (point estimate on the full matched sample, not the
        bootstrap mean), ``ci_lo_2p5``, ``ci_hi_97p5``, ``frac_delta_gt_0``,
        ``n_boot``, ``n_boot_valid``, ``n_boot_nan``, ``seed``, ``status``.
    """
    sentinels = sentinels if sentinels is not None else DEFAULT_SENTINELS
    stratum_n = len(rows)
    cleaned, n_excluded_missing, n_excluded_censored = _clean_matched(
        rows, predictor_a, predictor_b, endpoint_col, sentinels, censoring_columns
    )
    n_excluded = n_excluded_missing + n_excluded_censored

    out: dict[str, Any] = {
        "predictor_a": predictor_a, "predictor_b": predictor_b, "endpoint": endpoint_col,
        "cluster_col": cluster_col,
        "stratum_n": stratum_n, "n": len(cleaned), "n_excluded": n_excluded,
        "n_excluded_missing": n_excluded_missing, "n_excluded_censored": n_excluded_censored,
        "n_clusters": len(set(r[cluster_col] for r in cleaned)) if cleaned else 0,
        "tau_a_point": None, "tau_b_point": None, "delta_point": None,
        "ci_lo_2p5": None, "ci_hi_97p5": None, "frac_delta_gt_0": None,
        "n_boot": n_boot, "n_boot_valid": 0, "n_boot_nan": 0, "seed": seed,
        "status": "ok",
    }
    if len(cleaned) < 2:
        out["status"] = "insufficient_n"
        return out

    xa = np.array([float(r[predictor_a]) for r in cleaned])
    xb = np.array([float(r[predictor_b]) for r in cleaned])
    y = np.array([float(r[endpoint_col]) for r in cleaned])
    clusters = [r[cluster_col] for r in cleaned]

    tau_a_full = _kendall_tau_b(xa, y)
    tau_b_full = _kendall_tau_b(xb, y)
    out["tau_a_point"], out["tau_b_point"] = tau_a_full, tau_b_full
    if tau_a_full is None or tau_b_full is None:
        out["status"] = "undefined (predictor or endpoint constant on matched sample)"
        return out
    out["delta_point"] = tau_a_full - tau_b_full

    uniq_clusters = sorted(set(clusters))
    cluster_to_idx: dict[Any, list[int]] = defaultdict(list)
    for i, c in enumerate(clusters):
        cluster_to_idx[c].append(i)
    cluster_arr = np.array(uniq_clusters, dtype=object)
    n_clusters = len(uniq_clusters)

    rng = np.random.default_rng(seed)
    deltas = np.full(n_boot, np.nan, dtype=float)
    for b in range(n_boot):
        drawn = cluster_arr[rng.integers(0, n_clusters, size=n_clusters)]
        idx: list[int] = []
        for c in drawn:
            idx.extend(cluster_to_idx[c])
        idx_arr = np.array(idx)
        ta = _kendall_tau_b(xa[idx_arr], y[idx_arr])
        tb = _kendall_tau_b(xb[idx_arr], y[idx_arr])
        if ta is None or tb is None:
            continue
        deltas[b] = ta - tb

    valid = deltas[~np.isnan(deltas)]
    out["n_boot_valid"] = int(len(valid))
    out["n_boot_nan"] = int(n_boot - len(valid))
    if len(valid) < 100:
        out["status"] = "bootstrap_degenerate"
        return out

    lo, hi = np.percentile(valid, [2.5, 97.5])
    out["ci_lo_2p5"], out["ci_hi_97p5"] = float(lo), float(hi)
    out["frac_delta_gt_0"] = float(np.mean(valid > 0))
    return out


def leave_one_cluster_out_delta(
    rows: list[dict[str, Any]],
    predictor_a: str,
    predictor_b: str,
    endpoint_col: str,
    cluster_col: str,
    *,
    sentinels: dict[str, tuple[Any, ...]] | None = None,
    censoring_columns: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Point-estimate ``delta`` with each cluster's rows removed in turn.

    Mirrors ``real_analyse._leave_one_network_out``, applied to the paired
    difference rather than a single tau.

    Args:
        rows: Candidate rows.
        predictor_a: First predictor column name.
        predictor_b: Second predictor column name.
        endpoint_col: Endpoint column name.
        cluster_col: Cluster column name.
        sentinels: See :data:`DEFAULT_SENTINELS`.
        censoring_columns: See :func:`_clean_matched`.

    Returns:
        ``{"n", "n_excluded", "n_excluded_missing", "n_excluded_censored",
        "loo": {cluster: delta}, "delta_min", "delta_max", "cluster_at_min",
        "cluster_at_max", "status"}``. ``loo`` only carries clusters whose
        removal still leaves both predictors and the endpoint non-constant.
    """
    sentinels = sentinels if sentinels is not None else DEFAULT_SENTINELS
    cleaned, n_excluded_missing, n_excluded_censored = _clean_matched(
        rows, predictor_a, predictor_b, endpoint_col, sentinels, censoring_columns
    )
    out: dict[str, Any] = {
        "n": len(cleaned), "n_excluded": n_excluded_missing + n_excluded_censored,
        "n_excluded_missing": n_excluded_missing, "n_excluded_censored": n_excluded_censored,
        "loo": {},
        "delta_min": None, "delta_max": None, "cluster_at_min": None, "cluster_at_max": None,
        "status": "ok",
    }
    if len(cleaned) < 2:
        out["status"] = "insufficient_n"
        return out

    xa = np.array([float(r[predictor_a]) for r in cleaned])
    xb = np.array([float(r[predictor_b]) for r in cleaned])
    y = np.array([float(r[endpoint_col]) for r in cleaned])
    clusters = np.array([r[cluster_col] for r in cleaned], dtype=object)

    loo: dict[Any, float] = {}
    for c in sorted(set(clusters.tolist())):
        keep = clusters != c
        if keep.sum() < 2:
            continue
        ta = _kendall_tau_b(xa[keep], y[keep])
        tb = _kendall_tau_b(xb[keep], y[keep])
        if ta is None or tb is None:
            continue
        loo[c] = ta - tb
    out["loo"] = loo
    if loo:
        c_min = min(loo, key=lambda k: loo[k])
        c_max = max(loo, key=lambda k: loo[k])
        out["delta_min"], out["cluster_at_min"] = loo[c_min], c_min
        out["delta_max"], out["cluster_at_max"] = loo[c_max], c_max
    return out


def _verdict_flips(
    delta_point: float | None, ci_lo: float | None, ci_hi: float | None, loo: dict[Any, float]
) -> bool | None:
    """Whether the paired-delta verdict is sensitive to a single cluster.

    ``True`` iff either (a) some LOO removal flips the sign of ``delta_point``
    relative to the full-sample sign, or (b) the full-sample CI excludes zero
    while the LOO range still straddles zero. Same two-part rule
    ``real_analyse.tau_for_stratum`` uses for a single tau, applied here to
    the paired difference.

    Args:
        delta_point: Full-sample point estimate.
        ci_lo: Full-sample CI lower bound (or ``None``).
        ci_hi: Full-sample CI upper bound (or ``None``).
        loo: ``{cluster: delta}`` from :func:`leave_one_cluster_out_delta`.

    Returns:
        ``None`` if there is nothing to check (no LOO values or no point
        estimate), else the boolean verdict.
    """
    if not loo or delta_point is None:
        return None
    full_sign = 0.0 if delta_point == 0 else math.copysign(1.0, delta_point)
    opposite_sign = any(
        v != 0 and full_sign != 0 and math.copysign(1.0, v) != full_sign for v in loo.values()
    )
    ci_excludes_zero = ci_lo is not None and (ci_lo > 0 or ci_hi < 0)
    vmin, vmax = min(loo.values()), max(loo.values())
    flips_by_ci = ci_excludes_zero and (vmin <= 0 <= vmax)
    return bool(opposite_sign or flips_by_ci)


def paired_comparison(
    rows: list[dict[str, Any]],
    predictor_a: str,
    predictor_b: str,
    endpoint_col: str,
    cluster_col: str,
    *,
    n_boot: int = DEFAULT_N_BOOT,
    seed: int = DEFAULT_SEED,
    sentinels: dict[str, tuple[Any, ...]] | None = None,
    censoring_columns: dict[str, str] | None = None,
) -> dict[str, Any]:
    """The full paired statement: bootstrap delta + CI + P(delta>0), plus LOO.

    This is the primary entry point: predictor columns, cluster column and
    endpoint column in; one dict out, ready to become one row of a paired
    comparison table.

    Args:
        rows: Rows carrying ``predictor_a``, ``predictor_b``, ``endpoint_col``
            and ``cluster_col``.
        predictor_a: First predictor's column name.
        predictor_b: Second (baseline) predictor's column name.
        endpoint_col: Endpoint column name.
        cluster_col: Resampling unit's column name.
        n_boot: Bootstrap resamples (>= 10,000 recommended; see task brief).
        seed: Seed, recorded for reproducibility.
        sentinels: See :data:`DEFAULT_SENTINELS`.
        censoring_columns: Optional ``{predictor_column: flag_column}``. A
            predictor whose values can be right-censored (e.g. a depth-capped
            search that stops without finding a break) should name its own
            boolean flag column here so censored rows are dropped from the
            matched sample **and** counted under ``n_excluded_censored``,
            distinct from genuinely missing values. See :func:`_clean_matched`.

    Returns:
        The merged dict of :func:`paired_cluster_bootstrap` plus
        ``loo_delta_min``, ``loo_delta_max``, ``loo_cluster_at_min``,
        ``loo_cluster_at_max``, ``loo_n_clusters_with_delta``,
        ``loo_verdict_flips``.
    """
    boot = paired_cluster_bootstrap(
        rows, predictor_a, predictor_b, endpoint_col, cluster_col,
        n_boot=n_boot, seed=seed, sentinels=sentinels, censoring_columns=censoring_columns,
    )
    loo = leave_one_cluster_out_delta(
        rows, predictor_a, predictor_b, endpoint_col, cluster_col,
        sentinels=sentinels, censoring_columns=censoring_columns,
    )
    boot["loo_delta_min"] = loo["delta_min"]
    boot["loo_delta_max"] = loo["delta_max"]
    boot["loo_cluster_at_min"] = loo["cluster_at_min"]
    boot["loo_cluster_at_max"] = loo["cluster_at_max"]
    boot["loo_n_clusters_with_delta"] = len(loo["loo"])
    boot["loo_verdict_flips"] = _verdict_flips(
        boot["delta_point"], boot["ci_lo_2p5"], boot["ci_hi_97p5"], loo["loo"]
    )
    boot["loo_by_cluster"] = loo["loo"]  # full detail; not a CSV column, kept for JSON output
    return boot


# ---------------------------------------------------------------------------
# Real-corpus glue (imports real_analyse / real_baselines_p6 lazily)
# ---------------------------------------------------------------------------

REAL_FRAME_DIR = "results/axis_robustness_real"
REAL_BASELINES_DIR = "results/axis_robustness_real_p6"

#: r_val (radius) against each baseline, per the task brief.
REAL_BASELINE_PREDICTORS: tuple[str, ...] = (
    "shd_truth", "n_k", "k_g0", "separation", "phi_1", "r_claim",
)

#: phi_1/r_claim are only defined on flip-arm units at base_wrongness == 0.0
#: (see real_baselines_p6.build_extended_units's docstring for why: that is
#: the only case where the unit's operative K is the frame's own K). Rows
#: outside that scope carry phi_1 = r_claim = None by construction and are
#: reported as such, not silently substituted.
PHI1_RCLAIM_SCOPED_PREDICTORS = {"phi_1", "r_claim"}

PAIRED_COLUMNS: list[str] = [
    "stratum", "arm", "predictor_a", "predictor_b", "endpoint", "status",
    "stratum_n", "n", "n_excluded", "n_clusters",
    "tau_a_point", "tau_b_point", "delta_point",
    "ci_lo_2p5", "ci_hi_97p5", "frac_delta_gt_0",
    "n_boot", "n_boot_valid", "n_boot_nan", "seed",
    "loo_delta_min", "loo_delta_max", "loo_cluster_at_min", "loo_cluster_at_max",
    "loo_n_clusters_with_delta", "loo_verdict_flips",
    "scope_note",
]


def run_on_real_corpus(
    *,
    frame_dir: str = REAL_FRAME_DIR,
    baselines_dir: str = REAL_BASELINES_DIR,
    n_boot: int = DEFAULT_N_BOOT,
    seed: int = DEFAULT_SEED,
    log: Callable[[str], None] = print,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Runs :func:`paired_comparison` for ``radius`` against every baseline
    in :data:`REAL_BASELINE_PREDICTORS`, over every stratum of the real
    survival sweep.

    Requires ``results/axis_robustness_real_p6/baselines_per_instance.csv``
    to already exist (written by ``real_baselines_p6.py``'s ``__main__``);
    this function does not recompute ``phi_1``/``r_claim`` itself, since that
    is a many-minute combinatorial search better run once and cached.

    Args:
        frame_dir: The frozen real-survival results directory (read-only;
            never written here).
        baselines_dir: Where ``baselines_per_instance.csv`` lives.
        n_boot: Bootstrap resamples per (stratum, baseline) cell.
        seed: Seed for every bootstrap and LOO call.
        log: Verbose progress callback.

    Returns:
        ``(paired_rows, run_info)``.
    """
    import csv as _csv

    # Lazy imports: this keeps the module's top level free of any dependency
    # on the real-corpus results existing, per the module docstring.
    from bkrobust.robustness.real_analyse import group_units_by_stratum
    from bkrobust.robustness.real_baselines_p6 import build_extended_units
    from bkrobust.robustness.real_analyse import load_all_shards, verify_frame

    frame_path = Path(frame_dir)
    baselines_path = Path(baselines_dir)
    baselines_csv = baselines_path / "baselines_per_instance.csv"
    if not baselines_csv.is_file():
        raise FileNotFoundError(
            f"{baselines_csv} does not exist -- run "
            "`PYTHONPATH=src /usr/bin/python3 src/bkrobust/robustness/real_baselines_p6.py` first."
        )

    frame_rows, frame_info = verify_frame(frame_path)
    log(f"frame verified: {len(frame_rows)} rows, sha256={frame_info['frame_sha256'][:16]}...")

    baseline_by_row_id: dict[str, dict[str, Any]] = {}
    with baselines_csv.open(encoding="utf-8") as fh:
        for row in _csv.DictReader(fh):
            row_id = row["row_id"]
            phi_1 = float(row["phi_1"]) if row["phi_1"] not in ("", None) else None
            r_claim = float(row["r_claim"]) if row["r_claim"] not in ("", None) else None
            baseline_by_row_id[row_id] = {"phi_1": phi_1, "r_claim": r_claim}
    log(f"loaded {len(baseline_by_row_id)} baseline rows from {baselines_csv}")

    loaded_shards, skipped_shards = load_all_shards(frame_path)
    log(f"loaded {len(loaded_shards)} shards ({len(skipped_shards)} skipped)")
    units = build_extended_units(loaded_shards, baseline_by_row_id)
    stratum_groups = group_units_by_stratum(units)
    log(f"{len(stratum_groups)} strata")

    from bkrobust.robustness.real_analyse import ENDPOINTS_BY_ARM

    paired_rows: list[dict[str, Any]] = []
    t0 = time.perf_counter()
    for stratum in sorted(stratum_groups):
        s_units = stratum_groups[stratum]
        arm = s_units[0]["arm"]
        bw = s_units[0].get("base_wrongness")
        endpoint = ENDPOINTS_BY_ARM[arm][1]  # the *_usable endpoint (Sec. 5.1 / Appendix E)
        for baseline in REAL_BASELINE_PREDICTORS:
            scope_note = ""
            if baseline in PHI1_RCLAIM_SCOPED_PREDICTORS and not (arm == "flip" and bw == 0.0):
                scope_note = (
                    f"{baseline} is only defined where the unit's operative K is the frame's "
                    "own K (flip arm, base_wrongness=0.0); this stratum's units carry "
                    f"{baseline}=None throughout, so this row is expected to read "
                    "insufficient_n / n=0, not a genuine null finding."
                )
            r = paired_comparison(
                s_units, "radius", baseline, endpoint, "network",
                n_boot=n_boot, seed=seed,
            )
            r["stratum"], r["arm"], r["scope_note"] = stratum, arm, scope_note
            paired_rows.append(r)
        log(f"[stratum done] {stratum} ({arm}, n_units={len(s_units)}) "
            f"{round(time.perf_counter() - t0, 1)}s elapsed")

    run_info = {
        "frame_dir": frame_dir, "baselines_dir": baselines_dir,
        "n_boot": n_boot, "seed": seed,
        "frame_sha256": frame_info["frame_sha256"],
        "n_strata": len(stratum_groups),
        "strata": sorted(stratum_groups),
        "predictor_a": "radius",
        "baseline_predictors": list(REAL_BASELINE_PREDICTORS),
        "phi1_rclaim_scoped_to": "flip arm, base_wrongness == 0.0",
        "n_shards_read": len(loaded_shards),
        "n_shards_skipped": len(skipped_shards),
        "total_seconds": round(time.perf_counter() - t0, 1),
    }
    return paired_rows, run_info


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    """Builds the CLI parser."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--frame-dir", default=REAL_FRAME_DIR)
    p.add_argument("--baselines-dir", default=REAL_BASELINES_DIR)
    p.add_argument("--out-dir", default=REAL_BASELINES_DIR)
    p.add_argument("--n-boot", type=int, default=DEFAULT_N_BOOT)
    return p


def main(argv: list[str] | None = None) -> int:
    """Entry point: runs the paired real-corpus comparison and writes
    ``paired_diff_real.csv`` and ``paired_diff_manifest.json`` under
    ``--out-dir``.
    """
    from bkrobust.robustness.real_analyse import write_csv, write_json

    args = build_arg_parser().parse_args(argv)
    out_dir = Path(args.out_dir)
    seed = DEFAULT_SEED

    try:
        paired_rows, run_info = run_on_real_corpus(
            frame_dir=args.frame_dir, baselines_dir=args.baselines_dir,
            n_boot=args.n_boot, seed=seed, log=print,
        )
    except FileNotFoundError as exc:
        print(f"REFUSING TO RUN: {exc}")
        return 1

    write_csv(out_dir / "paired_diff_real.csv", PAIRED_COLUMNS, paired_rows)
    print(f"wrote paired_diff_real.csv: {len(paired_rows)} rows")

    # Full LOO detail (per-cluster deltas), for transparency beyond min/max.
    loo_detail = [
        {
            "stratum": r["stratum"], "predictor_a": r["predictor_a"], "predictor_b": r["predictor_b"],
            "loo_by_cluster": r.get("loo_by_cluster", {}),
        }
        for r in paired_rows
    ]
    write_json(out_dir / "paired_diff_loo_detail.json", loo_detail)
    print("wrote paired_diff_loo_detail.json")

    write_json(out_dir / "paired_diff_manifest.json", run_info)
    print("wrote paired_diff_manifest.json")

    n_positive_beats = sum(
        1 for r in paired_rows
        if r["status"] == "ok" and r["ci_lo_2p5"] is not None and r["ci_lo_2p5"] > 0
    )
    n_negative_beats = sum(
        1 for r in paired_rows
        if r["status"] == "ok" and r["ci_hi_97p5"] is not None and r["ci_hi_97p5"] < 0
    )
    n_loo_flips = sum(1 for r in paired_rows if r["loo_verdict_flips"] is True)
    print(json.dumps({
        "n_rows": len(paired_rows),
        "n_radius_significantly_outranks_baseline": n_positive_beats,
        "n_radius_significantly_underranks_baseline": n_negative_beats,
        "n_loo_verdict_flips": n_loo_flips,
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
