"""Asymptotic variance of an adjustment set, and the optimality test.

For a linear-Gaussian SCM the asymptotic variance of the adjusted least-squares
estimator of the total effect of ``X`` on ``Y`` given a valid adjustment set
``Z`` has a closed form in terms of residual variances. Comparing that quantity
across valid sets is what makes ``O*`` "optimal", and it is what
``delta_opt`` -- the radius at which ``O*`` stops being the minimiser -- is
defined against.

The distinction from validity matters and is easy to blur: an estimator built
on a valid-but-suboptimal set is still consistent. It converges to the right
number, just more slowly. Between ``delta_opt`` and ``delta_valid`` the analyst
pays in variance and nothing else; past ``delta_valid`` they pay in bias, and
no sample size rescues them. See :mod:`bkrobust.theory.efficiency_gap`.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    from bkrobust.graphs.mpdag import MPDAG, Node


def asymptotic_variance(
    graph: MPDAG,
    treatment: Node,
    outcome: Node,
    z: Iterable[Node],
    scm: Any,
) -> float:
    """Asymptotic variance of the ``Z``-adjusted total-effect estimator.

    Args:
        graph: The graph the estimator's adjustment set was read off. Used to
            check validity of ``z`` before the variance is meaningful.
        treatment: Treatment node.
        outcome: Outcome node.
        z: The adjustment set.
        scm: The structural causal model supplying edge coefficients and noise
            variances -- a :class:`bkrobust.data.synthetic.LinearGaussianSCM` or
            anything exposing the same covariance interface.

    Returns:
        The asymptotic variance, scaled so that it is comparable across
        adjustment sets on the same SCM (i.e. ``n * Var(tau_hat)`` in the limit).

    Raises:
        ValueError: If ``z`` is not a valid adjustment set for ``(treatment,
            outcome)`` in ``graph``, in which case the estimator is biased and
            its variance is not the quantity of interest.
    """
    raise NotImplementedError


def is_optimal_adjustment_set(
    graph: MPDAG,
    treatment: Node,
    outcome: Node,
    z: Iterable[Node],
    scm: Any,
    tol: float = 1e-9,
) -> bool:
    """Whether ``z`` attains the minimum asymptotic variance among valid sets.

    Args:
        graph: The graph defining which sets are valid.
        treatment: Treatment node.
        outcome: Outcome node.
        z: The candidate set.
        scm: The SCM supplying the variance.
        tol: Absolute tolerance on the variance comparison; ties within ``tol``
            count as optimal, since several sets can be genuinely tied.

    Returns:
        ``True`` if ``z`` is valid and no valid set has strictly smaller
        asymptotic variance.
    """
    raise NotImplementedError


def variance_ratio(
    graph: MPDAG,
    treatment: Node,
    outcome: Node,
    z_a: Iterable[Node],
    z_b: Iterable[Node],
    scm: Any,
) -> float:
    """Ratio ``asVar(z_a) / asVar(z_b)``.

    The headline efficiency number for the band between the two radii: with
    ``z_a`` the set chosen under false knowledge and ``z_b`` the set chosen
    under true knowledge, this is the multiplicative sample-size penalty the
    analyst pays for being wrong but not yet biased.

    Raises:
        ValueError: If either set is invalid, or ``asVar(z_b)`` is zero.
    """
    raise NotImplementedError


def minimum_variance_adjustment_set(
    graph: MPDAG,
    treatment: Node,
    outcome: Node,
    scm: Any,
) -> tuple[set[Node], float]:
    """Return the variance-minimising valid adjustment set, by enumeration.

    The brute-force reference for
    :func:`bkrobust.graphs.adjustment.optimal_adjustment_set`. On a graph where
    the theory applies these must agree, and ``tests/test_adjustment.py``
    asserts it; where they disagree, the graph is outside the theory's scope
    and that is worth knowing.

    Returns:
        The minimising set and its asymptotic variance.
    """
    raise NotImplementedError


def efficiency_loss(
    graph: MPDAG,
    treatment: Node,
    outcome: Node,
    z: Iterable[Node],
    scm: Any,
) -> float:
    """Excess asymptotic variance of ``z`` over the optimum, as a fraction.

    Returns:
        ``asVar(z) / asVar(O*) - 1``; zero when ``z`` is optimal.
    """
    raise NotImplementedError
