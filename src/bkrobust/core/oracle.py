"""The validity / optimality / bias oracle. Frozen interface.

``Z`` is fixed once from the analyst's graph ``G0`` and then evaluated in every
visited graph, each treated as a hypothetical truth. Validity in an MPDAG is a
**for-all over its represented DAGs**: ``Z`` is valid in ``G`` iff it is a valid
adjustment set in every DAG of ``[G]``. An MPDAG with no extensions is treated
as invalid rather than vacuously valid -- a graph representing no model cannot
certify anything.

Determinism notes, each of which was a bug once in this repository:

* extensions are consumed in the sorted order produced by
  :func:`bkrobust.demo.meek.enumerate_dag_extensions`, never in set order;
* SEM coefficient draws come only from the passed generator;
* floating-point accumulation is over sorted, deterministically ordered lists.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from bkrobust.demo.evaluate import (
    adjusted_estimand,
    asymptotic_variance,
    is_valid_adjustment_set_dag,
    optimal_adjustment_set_dag,
    random_sem,
)
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import enumerate_dag_extensions

#: Memoised ``[G]``. Keyed on the MPDAG, which hashes on a canonical key, so the
#: cache never changes iteration order or results -- only cost.
_EXT_CACHE: dict[MPDAG, tuple[MPDAG, ...]] = {}


def extensions(g: MPDAG) -> tuple[MPDAG, ...]:
    """Return ``[G]``, the consistent DAG extensions of ``g``, memoised.

    Args:
        g: An MPDAG.

    Returns:
        The extensions in the deterministic sorted order the enumerator emits.
    """
    hit = _EXT_CACHE.get(g)
    if hit is None:
        hit = tuple(enumerate_dag_extensions(g))
        _EXT_CACHE[g] = hit
    return hit


def clear_cache() -> None:
    """Drop the extension cache. Only for tests and memory management."""
    _EXT_CACHE.clear()


def is_valid(z: frozenset[str], g: MPDAG, x: str, y: str) -> bool:
    """Whether ``z`` is a valid adjustment set for ``(x, y)`` in every DAG of ``[g]``.

    Args:
        z: The adjustment set, fixed once from ``G0``.
        g: The graph being treated as a hypothetical truth.
        x: Treatment.
        y: Outcome.

    Returns:
        ``True`` iff ``z`` is valid in every extension. ``False`` if ``[g]`` is
        empty.
    """
    exts = extensions(g)
    if not exts:
        return False
    return all(is_valid_adjustment_set_dag(d, x, y, z) for d in exts)


def is_optimal(z: frozenset[str], g: MPDAG, x: str, y: str) -> bool:
    """Whether ``z`` equals the optimal adjustment set in every DAG of ``[g]``.

    Returns:
        ``True`` iff every extension's optimal set equals ``z``. ``False`` if
        ``[g]`` is empty or the extensions disagree.
    """
    exts = extensions(g)
    if not exts:
        return False
    return all(frozenset(optimal_adjustment_set_dag(d, x, y)) == z for d in exts)


@dataclass(frozen=True)
class BiasStats:
    """Bias of the ``z``-adjusted estimator at one graph, over sampled SEMs.

    Attributes:
        mean_abs_bias: Mean ``|estimand - true effect|`` over extensions and
            coefficient draws. **This is what radii are keyed off**: it is
            stable across runs.
        max_abs_bias: Maximum over the same product. A **sampled proxy for a
            worst case, not a supremum** -- coefficients are drawn, not
            optimised over -- and measurably unstable. Diagnostic only; never
            quote it as a bound.
        mean_abs_effect: Mean ``|true total effect|``, for putting bias on a
            scale.
        n_evaluations: Extensions times draws, the sample size behind the above.
    """

    mean_abs_bias: float
    max_abs_bias: float
    mean_abs_effect: float
    n_evaluations: int


def bias_stats(
    z: frozenset[str],
    g: MPDAG,
    x: str,
    y: str,
    rng: np.random.Generator,
    n_draws: int,
) -> BiasStats:
    """Sampled bias of adjusting for ``z`` at ``g``.

    Instantiates a linear-Gaussian SEM on **every** DAG of ``[g]`` crossed with
    ``n_draws`` coefficient draws, and reports the mean and max absolute bias.

    Args:
        z: The adjustment set.
        g: The graph treated as a hypothetical truth.
        x: Treatment.
        y: Outcome.
        rng: The sole source of randomness.
        n_draws: Coefficient draws per extension.

    Returns:
        A :class:`BiasStats`. All-zero with ``n_evaluations == 0`` when ``[g]``
        is empty.
    """
    exts = extensions(g)
    if not exts:
        return BiasStats(0.0, 0.0, 0.0, 0)
    biases: list[float] = []
    effects: list[float] = []
    for d in exts:
        for _ in range(n_draws):
            sem = random_sem(d, rng)
            tau = sem.true_total_effect(x, y)
            est = adjusted_estimand(sem, x, y, z)
            biases.append(abs(est - tau))
            effects.append(abs(tau))
    return BiasStats(
        mean_abs_bias=float(np.mean(biases)),
        max_abs_bias=float(np.max(biases)),
        mean_abs_effect=float(np.mean(effects)),
        n_evaluations=len(biases),
    )


def mean_asymptotic_variance(
    z: frozenset[str],
    dag: MPDAG,
    x: str,
    y: str,
    rng: np.random.Generator,
    n_draws: int,
) -> float:
    """Mean n-free asymptotic variance of the ``z``-adjusted estimator on ``dag``.

    Used for the efficiency axis of the robustness frontier.
    """
    return float(
        np.mean([asymptotic_variance(random_sem(dag, rng), x, y, z) for _ in range(n_draws)])
    )
