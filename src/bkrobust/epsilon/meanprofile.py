r"""The mean bias profile: ``mu(d)``, its radius, its bounds and its estimators.

``beta_up(d)`` -- the worst case over the retraction shell -- is what
``r_eps`` thresholds, and it is exact, cheap and proved. It is also, on this
corpus, almost perfectly flat above ``r_val``: it jumps to ``B(Chat)`` at the
validity boundary and stays there, so ``r_eps`` reproduces ``r_val`` and the
epsilon axis carries no information (``docs/RESULT_EPSILON_ABOVE_ONE.md``).

This module computes a companion quantity that does vary with depth:

.. math::

    \mu(d) \;=\; \mathbb{E}\big[\,B(\mathrm{Meek}(\hat C,\;K_{G_0}\setminus S))\,\big],
    \qquad S \sim \mathrm{Uniform}\big(\{S \subseteq K_{G_0} : |S| = d\}\big)

in words: **if exactly ``d`` of the orientations the analyst asserted are wrong,
and nothing says which ``d``, this is the expected bias.** The expectation is
over subsets, not over distinct closures: several subsets can Meek-close to the
same state and each one is a distinct way of being wrong, so each is counted.

What this is and is not
-----------------------

``mu`` is **not** an estimate of the mean of ``B`` over the sphere of the full
space at distance ``d``, and must not be reported as one. Those two differ, in
both directions, and a radius read off ``mu`` disagrees with one read off the
sphere mean on 29% of instances -- see ``docs/RESULT_MEAN_SOUNDNESS_GATE.md``.
The domain here is deliberately the retraction up-set: it is the set of ways the
*analyst's own claims* can be wrong, which is the question a practitioner asks.
States elsewhere in the space assert orientations the analyst never made and are
not ways for the analyst to be mistaken.

Consequences of that choice, all of which must travel with any reported number:

* ``mu`` is an **average-case** quantity under a **uniform prior** over which
  ``d`` claims are wrong. That prior is an assumption, and it is not implied by
  anything else in the framework.
* ``mu`` does **not** bound the worst case. ``beta_up`` does, and remains the
  quantity any guarantee must be stated against. ``mu(d) <= beta_up(d)`` always.
* Monotonicity of ``mu`` in ``d`` is **observed, not proved**
  (``docs/MEAN_REPORTING.md``). Nothing here assumes it; the radius helper
  reports whether the profile it was handed is monotone so a caller can say so.

Two facts make this affordable
------------------------------

* **Below ``r_val`` the profile is exactly zero, for free.** Theorem A gives
  ``beta(d) = 0`` for ``d < r_val``; a maximum of zero over non-negative values
  forces every state in those shells to zero, hence ``mu(d) = 0``. So the walk
  starts at ``r_val``, never below, and this is exact rather than a heuristic.
* **``mu`` is bracketed by cheap statistics.** With ``p(d)`` the fraction of the
  shell carrying non-zero bias,

  .. math:: p(d)\cdot \mathrm{min}^{+}(d) \;\le\; \mu(d) \;\le\; p(d)\cdot \beta_{\mathrm{up}}(d)

  and the upper bound is *exact* whenever every contaminated state in the shell
  carries the same bias, which is the common case. See :func:`mean_bounds`.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Iterator, Sequence
from dataclasses import dataclass

import numpy as np

from bkrobust.core.conventions import UNREACHED
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.epsilon.bias import BiasContext, bias_at, knowledge_of

Node = str
Edge = tuple[str, str]


@dataclass
class ShellStat:
    """The bias distribution over one retraction shell.

    Attributes:
        d: Retraction depth.
        n_subsets_total: ``C(|K_G0|, d)``, the size of the population.
        n_evaluated: Subsets actually drawn. Equals ``n_subsets_total`` when
            exhaustive.
        n_states: Distinct Meek closures among the drawn subsets.
        mean: ``mu(d)`` in the context's reporting units -- the subset-weighted
            mean, which is the estimand. An estimate when ``not exhaustive``.
        se: Standard error of ``mean``. Exactly ``0.0`` when exhaustive: the
            population has been enumerated, so there is nothing to be uncertain
            about. Finite-population corrected otherwise.
        beta_up: ``max B`` over the drawn subsets. A lower bound on the true
            shell maximum when sampled.
        min_nonzero: Smallest non-zero ``B`` seen, or ``0.0`` if none.
        frac_nonzero: Fraction of drawn subsets whose state carries ``B > 0``.
        exhaustive: Whether the whole shell was enumerated.
        seconds: Wall time.
    """

    d: int
    n_subsets_total: int
    n_evaluated: int
    n_states: int
    mean: float
    se: float
    beta_up: float
    min_nonzero: float
    frac_nonzero: float
    exhaustive: bool
    seconds: float = 0.0

    def ci(self, z: float = 1.96) -> tuple[float, float]:
        """A normal confidence interval for ``mu(d)``, clipped at zero.

        Args:
            z: Normal quantile; 1.96 is the usual 95% two-sided value.

        Returns:
            ``(lo, hi)``. Degenerate to ``(mean, mean)`` when exhaustive.
        """
        if self.exhaustive or self.se == 0.0:
            return (self.mean, self.mean)
        return (max(0.0, self.mean - z * self.se), self.mean + z * self.se)


def _subsets(n: int, d: int) -> Iterator[tuple[int, ...]]:
    """All ``d``-subsets of ``range(n)``."""
    return itertools.combinations(range(n), d)


def _sample_subsets(
    n: int, d: int, m: int, rng: np.random.Generator
) -> list[tuple[int, ...]]:
    """``m`` distinct ``d``-subsets of ``range(n)``, drawn uniformly.

    Rejection-samples distinct subsets, which is sampling **without**
    replacement from the population of subsets -- the scheme the finite
    population correction in :func:`shell_stat` assumes.

    Args:
        n: Population size ``|K_G0|``.
        d: Subset size.
        m: How many to draw.
        rng: Source of randomness.

    Returns:
        ``m`` distinct subsets, or every subset when ``m`` reaches the total.
    """
    total = math.comb(n, d)
    if m >= total:
        return list(_subsets(n, d))
    seen: set[tuple[int, ...]] = set()
    # Rejection sampling is efficient while m is a modest fraction of the
    # population; above that, enumerate and permute instead of spinning on
    # collisions.
    if m > total // 2:
        allsub = list(_subsets(n, d))
        idx = rng.choice(len(allsub), size=m, replace=False)
        return [allsub[i] for i in idx]
    while len(seen) < m:
        pick = tuple(sorted(rng.choice(n, size=d, replace=False).tolist()))
        seen.add(pick)
    return sorted(seen)


def shell_stat(
    ctx: BiasContext,
    cpdag: MPDAG,
    g0: MPDAG,
    d: int,
    *,
    scale: float = 1.0,
    sample: int | None = None,
    rng: np.random.Generator | None = None,
    cache: dict[str, float] | None = None,
    method: str = "semilocal",
) -> ShellStat:
    """The bias distribution over retraction shell ``d``, exhaustive or sampled.

    Args:
        ctx: The per-instance bias context.
        cpdag: The CPDAG.
        g0: The analyst's state.
        d: Retraction depth.
        scale: Divide every bias by this. Pass ``abs(theta_z)`` for relative
            units, matching what ``certify`` reports.
        sample: Draw this many subsets instead of enumerating. ``None`` or a
            value at or above ``C(|K_G0|, d)`` enumerates exhaustively.
        rng: Randomness for sampling; required when ``sample`` bites.
        cache: Edge-string to bias cache, shared across shells of one instance.
            ``B`` depends only on the closed state and many subsets collapse
            onto the same one, so this is where most of the saving is.
        method: Passed to :func:`~bkrobust.epsilon.bias.bias_at`.

    Returns:
        A :class:`ShellStat`.

    Raises:
        ValueError: If ``d`` is out of range, or sampling is requested with no
            generator.
    """
    import time

    t0 = time.perf_counter()
    k0 = knowledge_of(cpdag, g0)
    n = len(k0)
    if d < 0 or d > n:
        raise ValueError(f"d must be in 0..{n}, got {d}")
    total = math.comb(n, d)
    if cache is None:
        cache = {}

    exhaustive = sample is None or sample >= total
    if exhaustive:
        picks: Sequence[tuple[int, ...]] = list(_subsets(n, d))
    else:
        if rng is None:
            raise ValueError("sampling requires rng=")
        picks = _sample_subsets(n, d, int(sample), rng)

    values: list[float] = []
    states: set[str] = set()
    for drop in picks:
        dropped = set(drop)
        keep = [k0[i] for i in range(n) if i not in dropped]
        h = apply_orientations(cpdag, keep)
        if h is None:  # unreachable: a subset of a consistent set is consistent
            continue
        es = h.edge_string()
        if es not in cache:
            cache[es] = bias_at(ctx, h, method=method).worst / (scale or 1.0)
        values.append(cache[es])
        states.add(es)

    arr = np.asarray(values, dtype=float)
    mean = float(arr.mean()) if arr.size else 0.0
    if exhaustive or arr.size < 2:
        se = 0.0
    else:
        # Finite population correction: the population of subsets is C(n, d),
        # not infinite, and at a large sampling fraction the naive s/sqrt(m)
        # overstates the uncertainty badly.
        var = float(arr.var(ddof=1))
        se = math.sqrt(max(0.0, var / arr.size) * max(0.0, 1.0 - arr.size / total))
    nonzero = arr[arr > 0]
    return ShellStat(
        d=d,
        n_subsets_total=total,
        n_evaluated=int(arr.size),
        n_states=len(states),
        mean=mean,
        se=se,
        beta_up=float(arr.max()) if arr.size else 0.0,
        min_nonzero=float(nonzero.min()) if nonzero.size else 0.0,
        frac_nonzero=float(nonzero.size / arr.size) if arr.size else 0.0,
        exhaustive=exhaustive,
        seconds=time.perf_counter() - t0,
    )


def mean_bounds(stat: ShellStat) -> tuple[float, float]:
    """Bracket ``mu(d)`` from the shell's contamination and extremes.

    With ``p`` the contaminated fraction, ``m+`` the smallest non-zero bias and
    ``beta`` the shell maximum, every contaminated state contributes at least
    ``m+`` and at most ``beta``, and the clean ones contribute zero:

    ``p * m+  <=  mu  <=  p * beta``

    The upper bound is attained exactly when every contaminated state in the
    shell carries the same bias, which is the common case on this corpus (81%
    of shells), making it a tight bound rather than a loose one.

    Args:
        stat: The shell's statistics.

    Returns:
        ``(lower, upper)``.
    """
    return (stat.frac_nonzero * stat.min_nonzero, stat.frac_nonzero * stat.beta_up)


def r_mean(
    profile: Sequence[ShellStat], eps: float, *, use: str = "mean"
) -> int:
    """First depth whose ``mu`` exceeds ``eps``.

    Args:
        profile: Shells in increasing ``d``. Depths below ``r_val`` may be
            omitted: they are identically zero and can never cross a
            non-negative threshold.
        eps: Threshold, in the same units as the profile.
        use: ``"mean"``, or ``"ci_lo"`` / ``"ci_hi"`` to threshold the
            confidence bound instead -- ``"ci_lo"`` gives the conservative
            radius (crosses latest), ``"ci_hi"`` the optimistic one.

    Returns:
        The depth, or :data:`~bkrobust.core.conventions.UNREACHED`.

    Raises:
        ValueError: On an unknown ``use``.
    """
    for s in profile:
        if use == "mean":
            v = s.mean
        elif use == "ci_lo":
            v = s.ci()[0]
        elif use == "ci_hi":
            v = s.ci()[1]
        else:
            raise ValueError(f"unknown use {use!r}")
        if v > eps:
            return s.d
    return UNREACHED


def is_monotone(profile: Sequence[ShellStat], tol: float = 1e-12) -> bool:
    """Whether ``mu`` is non-decreasing over the profile as given.

    Monotonicity is not guaranteed by any result in this codebase, so callers
    that want to describe a crossing as a radius should check it and say so.

    Args:
        profile: Shells in increasing ``d``.
        tol: Slack for floating-point noise.

    Returns:
        ``True`` if no shell's mean falls below its predecessor's.
    """
    vals = [s.mean for s in profile]
    return all(vals[i] >= vals[i - 1] - tol for i in range(1, len(vals)))
