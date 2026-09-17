"""``r_eps`` from data rather than from a population covariance.

``r_val`` is a purely graphical quantity: it depends on the CPDAG, the asserted
knowledge and the query, and on no numbers at all. ``r_eps`` is not. It reads a
covariance matrix, and in practice that matrix is estimated, so the radius
inherits sampling uncertainty that ``r_val`` does not have. Pretending otherwise
would be the one mistake this module exists to prevent: a plug-in ``r_eps``
reported as though it were exact is a confidence statement with the confidence
removed.

**Why the uncertainty has a direction.** ``r_eps`` is the first depth at which
the worst-case bias crosses ``eps``. Under-estimating the bias pushes the
crossing outwards, which *overstates* robustness -- the dangerous direction, and
the same one the paper flags for the upward search. Over-estimating the bias
pulls the crossing inwards and understates robustness, which is safe. So a
conservative radius is obtained by thresholding an **upper** bound on the bias,
equivalently by taking a **lower** quantile of the radius's own sampling
distribution. That is what :func:`conservative_radius` returns.

**Why the whole staircase is resampled, not each state separately.** A per-state
confidence band would need a multiplicity correction across the shell, and the
correction would have to respect the monotone structure or it would produce a
non-monotone band and hence an ill-defined crossing. Resampling the *radius*
sidesteps both problems: each bootstrap replicate computes a complete, internally
consistent staircase from one resampled covariance, so every replicate's radius
is a legitimate radius, and the quantile of legitimate radii is a legitimate
bound. It costs ``R`` profile evaluations, which is affordable because a profile
is cheap once shells below ``r_val`` are skipped.

Determinism: every draw comes from the passed generator; the global NumPy RNG is
never touched, and bootstrap indices are drawn in a fixed order.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from bkrobust.core.conventions import UNREACHED
from bkrobust.demo.evaluate import LinearSEM
from bkrobust.demo.graph import MPDAG
from bkrobust.epsilon.bias import make_context_from_covariance
from bkrobust.epsilon.profile import DEFAULT_SHELL_CAP, epsilon_grid

Node = str


def sample_data(sem: LinearSEM, n: int, rng: np.random.Generator) -> np.ndarray:
    """Draw ``n`` i.i.d. observations from a linear-Gaussian SEM.

    Drawn from the SEM's exact population covariance rather than by forward
    simulation through the structural equations. The two are distributionally
    identical for a linear-Gaussian model, and going through the covariance
    keeps this function insensitive to the node ordering of the structural
    equations, which forward simulation is not.

    Args:
        sem: The SEM.
        n: Number of observations.
        rng: Sole source of randomness.

    Returns:
        An ``(n, p)`` array, columns in ``sem.dag.nodes`` order.

    Raises:
        ValueError: If ``n`` is not positive.
    """
    if n <= 0:
        raise ValueError(f"n must be positive, got {n}")
    sigma = sem.covariance()
    mean = np.zeros(sigma.shape[0])
    return rng.multivariate_normal(mean, sigma, size=n, method="cholesky")


def sample_covariance(data: np.ndarray) -> np.ndarray:
    """The unbiased sample covariance of ``data``.

    Args:
        data: An ``(n, p)`` array.

    Returns:
        The ``(p, p)`` sample covariance.

    Raises:
        ValueError: If fewer than two observations are supplied, where the
            unbiased estimate is undefined.
    """
    if data.shape[0] < 2:
        raise ValueError("at least two observations are needed for a covariance")
    return np.cov(data, rowvar=False, ddof=1)


@dataclass(frozen=True)
class RadiusDistribution:
    """The bootstrap sampling distribution of ``r_eps`` at one threshold.

    Attributes:
        eps: The threshold, in effect units.
        point: The plug-in radius from the observed data. A point estimate, not
            a guarantee.
        replicates: One radius per bootstrap replicate, including
            :data:`~bkrobust.core.conventions.UNREACHED` entries.
        n_unreached: How many replicates found no crossing. Reported separately
            because ``UNREACHED`` is a status and must never enter a quantile
            as a number -- it is treated as ``+infinity`` for ordering, which is
            what it means, and that is handled in :meth:`conservative`.
        n_replicates: Total replicates.
    """

    eps: float
    point: int
    replicates: tuple[int, ...]
    n_unreached: int
    n_replicates: int

    def conservative(self, alpha: float = 0.05) -> int:
        """The ``alpha``-quantile radius: the conservative, reportable number.

        A radius that is too small understates robustness, which is the safe
        error; so the lower quantile is the guarantee and the point estimate is
        not.

        Args:
            alpha: Tail probability. ``0.05`` gives a radius that the bootstrap
                distribution exceeds about 95% of the time.

        Returns:
            The radius at the ``alpha`` quantile, ordering ``UNREACHED`` above
            every finite radius (it means "no crossing anywhere", which is the
            largest possible answer, not a small one). Returns ``UNREACHED``
            only when at least ``1 - alpha`` of the replicates found no
            crossing at all.

        Raises:
            ValueError: If ``alpha`` is not in ``(0, 1)``.
        """
        if not 0.0 < alpha < 1.0:
            raise ValueError(f"alpha must be in (0,1), got {alpha}")
        big = max((r for r in self.replicates if r != UNREACHED), default=0) + 1
        ordered = sorted(big if r == UNREACHED else r for r in self.replicates)
        if not ordered:
            return UNREACHED
        idx = int(np.floor(alpha * (len(ordered) - 1)))
        value = ordered[idx]
        return UNREACHED if value == big else value


def bootstrap_radius(
    cpdag: MPDAG,
    g0: MPDAG,
    x: Node,
    y: Node,
    z: frozenset[Node],
    data: np.ndarray,
    nodes: Sequence[Node],
    epsilons: Sequence[float],
    *,
    r_val: int = 0,
    n_boot: int = 200,
    rng: np.random.Generator | None = None,
    method: str = "semilocal",
    cap: int = DEFAULT_SHELL_CAP,
) -> dict[float, RadiusDistribution]:
    """Bootstrap the sampling distribution of ``r_eps`` at each threshold.

    Each replicate resamples the rows of ``data`` with replacement, recomputes
    the covariance, and runs a complete staircase from it. Because the whole
    staircase is recomputed, every replicate's radius is internally consistent
    with its own covariance -- which a per-state confidence band would not be.

    ``r_val`` is **not** bootstrapped: it is a graphical quantity with no
    dependence on the data, so it is the same in every replicate and is passed
    through to skip the certified shells, exactly as in the population case.

    Args:
        cpdag: The estimated CPDAG.
        g0: The analyst's state.
        x: Treatment.
        y: Outcome.
        z: The adjustment set.
        data: An ``(n, p)`` observation array.
        nodes: Column names of ``data``.
        epsilons: Thresholds, in effect units. Note that a *relative* threshold
            would move between replicates, since the reported estimate does;
            convert to effect units once, on the observed data, before calling
            this, so that every replicate is judged against the same bar.
        r_val: The validity radius; shells below it are skipped.
        n_boot: Number of replicates.
        rng: Sole source of randomness; defaults to a fixed-seed generator so a
            forgotten seed cannot silently make a run irreproducible.
        method: Passed through.
        cap: Passed through.

    Returns:
        One :class:`RadiusDistribution` per threshold.
    """
    generator = np.random.default_rng(0) if rng is None else rng
    nodes = tuple(nodes)
    n = data.shape[0]

    observed = make_context_from_covariance(sample_covariance(data), nodes, cpdag, x, y, z)
    point_radii, _ = epsilon_grid(observed, g0, list(epsilons), r_val=r_val, method=method, cap=cap)

    draws: dict[float, list[int]] = {float(e): [] for e in epsilons}
    for _ in range(n_boot):
        idx = generator.integers(0, n, size=n)
        try:
            ctx = make_context_from_covariance(sample_covariance(data[idx]), nodes, cpdag, x, y, z)
        except np.linalg.LinAlgError:
            # A degenerate resample is a property of the resample, not of the
            # instance; dropping it silently would bias the quantile, so it is
            # recorded as a non-crossing, which is the conservative reading.
            for e in epsilons:
                draws[float(e)].append(UNREACHED)
            continue
        radii, _ = epsilon_grid(ctx, g0, list(epsilons), r_val=r_val, method=method, cap=cap)
        for e in epsilons:
            draws[float(e)].append(radii[e])

    return {
        float(e): RadiusDistribution(
            eps=float(e),
            point=point_radii[e],
            replicates=tuple(draws[float(e)]),
            n_unreached=sum(1 for r in draws[float(e)] if r == UNREACHED),
            n_replicates=n_boot,
        )
        for e in epsilons
    }
