"""Null-stratum isolation (Axis Robustness, post-session-7 follow-up).

Session 7's flip arm, coverage = 0.5, base_wrongness = 0.25 stratum returned
tau_b(AUC_frac, r_val) = +0.031, CI [-0.114, +0.169], n = 220 -- the one flip
stratum whose CI did not exclude zero. This module isolates whether that null
is genuine (GAC validity collapses so fast under this much base wrongness
that r_val stops discriminating) or an artifact of measurement noise
(AUC_frac at N=200 draws/depth is attenuated toward zero because this cell's
contradiction rate is unusually high, so few evaluable samples survive per
depth -- Appendix F of PREREGISTRATION.md already found
tau(r_val, usable depths) = -0.100 here, i.e. high-r_val instances rest on
*less* data in this cell specifically).

Design: hold the *instances* fixed -- same ComponentSpec grid, same seed
sequence, same base-wrongness corruption, same fast_gate, same Z* -- and vary
only the number of draws per depth, from 200 (session 7) to 1000. See
``results/axis_robustness/NULLCELL_RUN_NOTES.md`` for the full writeup.

This module only *imports* :mod:`bkrobust.robustness.survival` (read-only
reuse of its construction/scoring/AUC functions, per the task's process
rules) and never calls, imports, or otherwise runs
:mod:`bkrobust.robustness.run_survival` or any other ``run_*`` driver in this
package -- those sweeps are frozen, committed, and must not be re-executed.

Gate discipline, RNG discipline: identical to :mod:`bkrobust.robustness.
survival` -- this module adds no gate and no RNG source of its own. Every
stochastic draw here is made *inside* ``survival.py``'s already-audited
functions (``build_flip_instance``, ``sample_flip_state``), seeded via
``survival.derived_seed`` -> ``np.random.default_rng``. This module never
seeds a process-global RNG and never touches Python's stdlib random module.
"""

from __future__ import annotations

from typing import Any, Iterable

import numpy as np
from scipy.stats import kendalltau

from bkrobust.core.conventions import UNREACHED
from bkrobust.robustness import survival as sv

# --- the null cell, fixed to session 7's manifest.json grid -------------------

#: Component sizes swept in Phase 1 (manifest.json: component_sizes).
COMPONENT_SIZES: tuple[int, ...] = (6, 8, 10, 12)

#: The one flip-arm stratum under investigation.
COVERAGE = 0.5
BASE_WRONGNESS = 0.25

#: Session 7's per-cell instance quota and seed search budget (manifest.json:
#: instances_per_cell=20, seed_budget=150). Reusing these exactly is what
#: makes the regenerated instance set provably identical to session 7's.
INSTANCES_PER_CELL = 20
SEED_BUDGET = 150

#: Draw counts compared in this run.
N_DRAWS_CHECK = 200  # session 7's N -- reproduction check (1)
N_DRAWS_FULL = 1000  # this run's N -- checks (2) and (3)

#: Sensitivity threshold for "usable" depths, matching analyse.py's
#: AUC_frac_usable (n_eval >= 30).
USABLE_MIN_N_EVAL = 30

#: Bootstrap discipline (PREREGISTRATION.md section 7 / task spec): 10,000
#: resamples over instances, seed 0, 2.5/97.5 percentiles.
N_BOOT = 10_000
BOOT_SEED = 0


# --- instance population -------------------------------------------------------


def cells() -> list[tuple[int, int]]:
    """The (component_size, separation) cells that make up this stratum.

    Same enumeration order ``run_survival.py`` used: for each component size,
    the sorted set of ``separations_for`` (min, median, max -- deduplicated).
    """
    out: list[tuple[int, int]] = []
    for c in sorted(COMPONENT_SIZES):
        for s in sorted(set(sv.separations_for(c))):
            out.append((c, s))
    return out


def build_cell_instances(
    component_size: int, separation: int, *, exclusions: dict[str, int]
) -> list[dict[str, Any]]:
    """Fill one (component_size, separation) cell up to quota.

    Identical acceptance loop to ``run_survival.run_flip_cell``: walk
    ``seed in range(SEED_BUDGET)``, keep the first ``INSTANCES_PER_CELL``
    accepted via :func:`survival.build_flip_instance`. Rejections are
    tallied into ``exclusions`` by reason, never silently dropped.
    """
    accepted: list[dict[str, Any]] = []
    for seed in range(SEED_BUDGET):
        if len(accepted) >= INSTANCES_PER_CELL:
            break
        inst, reason = sv.build_flip_instance(component_size, separation, COVERAGE, BASE_WRONGNESS, seed)
        if inst is None:
            exclusions[reason] = exclusions.get(reason, 0) + 1
            continue
        accepted.append(inst)
    return accepted


def build_all_instances() -> tuple[list[dict[str, Any]], dict[str, int]]:
    """All accepted instances across every cell in this stratum, in order."""
    exclusions: dict[str, int] = {}
    instances: list[dict[str, Any]] = []
    for c, s in cells():
        instances.extend(build_cell_instances(c, s, exclusions=exclusions))
    return instances, exclusions


def check_s0(instances: Iterable[dict[str, Any]]) -> list[str]:
    """K2: S(0) = 1.0 for every instance (Z* is valid at its own G0 by
    construction -- ``score_committed_graph`` already checked this at build
    time, this is an independent re-check on the built instances)."""
    failures = []
    for inst in instances:
        ok = sv.is_gac_valid_mpdag(inst["g0"], inst["x"], inst["y"], inst["z_star"])
        if not ok:
            failures.append(inst["instance_id"])
    return failures


# --- sampling --------------------------------------------------------------------


def sample_instance(inst: dict[str, Any], n_draws: int = N_DRAWS_FULL) -> list[dict[str, Any]]:
    """Draw ``n_draws`` flip-arm samples at every depth ``d in 1..n_k``.

    Reuses :func:`survival.sample_flip_state` verbatim, so ``rep in
    range(N_DRAWS_CHECK)`` reproduces session 7's exact draws bit-for-bit
    (same ``derived_seed(instance_id, "flip", d, rep)`` derivation) and
    ``rep in range(N_DRAWS_CHECK, n_draws)`` extends them with new,
    independently seeded draws.
    """
    rows: list[dict[str, Any]] = []
    n_k = inst["n_k"]
    for d in range(1, n_k + 1):
        for rep in range(n_draws):
            rows.append(sv.sample_flip_state(inst, d, rep))
    return rows


# --- curves, AUC, and its analytic SE -------------------------------------------


def curve_from_rows(rows: list[dict[str, Any]], *, rep_limit: int | None = None) -> dict[int, dict[str, Any]]:
    """Build a per-depth curve (:func:`survival.build_curve`) from a subset of
    ``rows``, restricted to ``rep < rep_limit`` when given. Used to build the
    N=200 curve as the literal first-200-reps prefix of the N=1000 draws.
    """
    if rep_limit is not None:
        rows = [r for r in rows if r["rep"] < rep_limit]
    return sv.build_curve(rows)


def contradiction_rate_by_depth(curve: dict[int, dict[str, Any]]) -> dict[int, float]:
    """``n_contradictory / (n_eval + n_contradictory)`` per depth."""
    out = {}
    for d, row in curve.items():
        total = row["n_eval"] + row["n_contradictory"]
        out[d] = (row["n_contradictory"] / total) if total else float("nan")
    return out


def _auc_frac_weights(defined_depths: list[int], n_k: int) -> dict[int, float]:
    """Replica of :func:`survival.auc_frac`'s fractional-grid weighting.

    ``auc_frac`` maps each of the 10 ``FRAC_GRID`` points to the nearest
    *defined* depth (ties toward the smaller depth) and averages those 10
    values. This returns the resulting ``{depth: weight}`` (weights sum to
    1.0) so the same weights can be used to analytically propagate per-depth
    binomial variance into a standard error for AUC_frac -- it must track
    ``auc_frac``'s own nearest-depth rule exactly, or the SE would not
    correspond to the number actually reported.
    """
    weights: dict[int, float] = {}
    n_grid = len(sv.FRAC_GRID)
    for frac in sv.FRAC_GRID:
        target = frac * n_k
        nearest = min(defined_depths, key=lambda d: (abs(d - target), d))
        weights[nearest] = weights.get(nearest, 0.0) + 1.0 / n_grid
    return weights


def auc_frac_se(curve: dict[int, dict[str, Any]], n_k: int) -> float:
    """Analytic standard error of ``AUC_frac`` from per-depth binomial
    sampling variance, propagated through :func:`_auc_frac_weights`.

    Depths are drawn with independent derived seeds (different ``d`` ->
    different ``derived_seed`` -> independent ``Generator`` streams), so no
    covariance term between depths is needed; ``S(d)`` is a Bernoulli
    proportion over ``n_eval(d)`` non-contradictory draws, variance
    ``p(1-p)/n_eval``.
    """
    defined_depths = sorted(d for d, row in curve.items() if row["S"] != "")
    if not defined_depths or n_k <= 0:
        return float("nan")
    weights = _auc_frac_weights(defined_depths, n_k)
    var = 0.0
    for d, w in weights.items():
        p = curve[d]["S"]
        n_eval = curve[d]["n_eval"]
        if n_eval <= 0:
            continue
        var += (w ** 2) * (p * (1 - p) / n_eval)
    return float(var ** 0.5)


def auc_frac_usable(curve: dict[int, dict[str, Any]], n_k: int, *, min_n_eval: int = USABLE_MIN_N_EVAL):
    """``AUC_frac`` restricted to depths with ``n_eval >= min_n_eval`` --
    same algorithm as ``analyse.py``'s ``compute_auc_frac_usable``, reusing
    :func:`survival.auc_frac` for the averaging rule so the sensitivity
    endpoint cannot drift from the primary one.
    """
    restricted = {d: row for d, row in curve.items() if row["S"] != "" and row["n_eval"] >= min_n_eval}
    if not restricted:
        return ""
    return sv.auc_frac(restricted, n_k)


def n_usable_depths(curve: dict[int, dict[str, Any]], *, min_n_eval: int = USABLE_MIN_N_EVAL) -> int:
    return sum(1 for row in curve.values() if row["S"] != "" and row["n_eval"] >= min_n_eval)


# --- tau_b with bootstrap CI (mirrors analyse.py's tau_for_stratum exactly) ----


def bootstrap_tau_ci(x: np.ndarray, y: np.ndarray, *, n_boot: int = N_BOOT, seed: int = BOOT_SEED):
    """10,000 resamples over instances, explicit Generator, seed 0, 2.5/97.5
    percentiles -- identical procedure to ``analyse.py::_bootstrap_tau_ci``.
    """
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


def tau_b_stratum(r_val: list[float], endpoint: list[float], *, n_boot: int = N_BOOT, seed: int = BOOT_SEED) -> dict[str, Any]:
    """tau_b(endpoint, r_val) with bootstrap CI, r_val==UNREACHED and NaN
    endpoints excluded -- identical exclusion + statistic to
    ``analyse.py::tau_for_stratum`` for ``predictor="r_val"``.
    """
    r = np.asarray(r_val, dtype=float)
    e = np.asarray(endpoint, dtype=float)
    total_n = len(r)
    n_excluded_unreached = int((r == UNREACHED).sum())
    mask = (r != UNREACHED) & ~np.isnan(r) & ~np.isnan(e)
    r, e = r[mask], e[mask]
    n = len(r)
    n_excluded_nan = max(total_n - n_excluded_unreached - n, 0)

    row: dict[str, Any] = {
        "n_total": total_n,
        "n": n,
        "n_excluded_unreached": n_excluded_unreached,
        "n_excluded_nan": n_excluded_nan,
        "tau_b": float("nan"),
        "p_value": float("nan"),
        "ci_lo_2p5": float("nan"),
        "ci_hi_97p5": float("nan"),
        "n_boot_nan": 0,
    }
    if n < 2 or np.unique(r).size <= 1 or np.unique(e).size <= 1:
        return row
    tau, pval = kendalltau(r, e, variant="b", nan_policy="propagate")
    lo, hi, n_boot_nan = bootstrap_tau_ci(r, e, n_boot=n_boot, seed=seed)
    row.update({"tau_b": float(tau), "p_value": float(pval), "ci_lo_2p5": lo, "ci_hi_97p5": hi, "n_boot_nan": n_boot_nan})
    return row
