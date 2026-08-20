"""The two breakdown radii, delta_valid and delta_opt.

Both are minima over a ball of consistent-but-false knowledge sets. Writing
``K_delta(C, G)`` for the knowledge sets that are consistent with the CPDAG
``C``, false with respect to the true DAG ``G``, and at distance ``delta`` from
the truth, and ``O*(K)`` for the optimal adjustment set of the MPDAG that ``K``
induces:

    delta_valid = min { delta : exists K in K_delta with O*(K) invalid in G }
    delta_opt   = min { delta : exists K in K_delta with O*(K) not optimal in G }

Both are worst-case over the ball -- the smallest perturbation at which
*something* can go wrong, not the typical one. That is the right definition for
a guarantee ("below this radius you are safe regardless of which wrong claims
the expert made") and the wrong one for a forecast, so the sweeps also report
the empirical distribution over the ball alongside these minima.

**Two computational paths.**

*Exact.* Enumerate ``K_delta`` for increasing ``delta`` and test each. Correct
by construction, and combinatorial: the ball grows like the number of
``delta``-subsets of orientable edges, times the DAG extensions needed to
evaluate optimality. Tractable to roughly a dozen nodes, gated by the
``radius.exact_max_nodes`` config knob.

*Heuristic.* Search the ball greedily, ordering candidate perturbations by
their effect on ``O*`` -- perturbations that touch the parents of the causal
nodes first, since only those can move ``O*`` at all. The result is an **upper
bound**: the search finds a breaking knowledge set at some radius, but cannot
certify that no smaller one exists. Every reported heuristic radius must be
labelled as an upper bound, and :func:`radius_report` carries the ``exact`` flag
that says which path produced it. Reporting a heuristic radius as if it were
exact would overstate the guarantee, which is the one error this module must
not make.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    from bkrobust.graphs.mpdag import MPDAG, Node
    from bkrobust.knowledge.base import BackgroundKnowledge


#: Returned when no perturbation within the search budget breaks the property.
#: Distinct from a large finite radius: it means "not found", not "infinite".
UNREACHED: int = -1


@dataclass(frozen=True)
class RadiusReport:
    """A computed radius together with everything needed to interpret it.

    Attributes:
        radius: The radius, or :data:`UNREACHED` if no breaking perturbation was
            found within budget.
        exact: Whether the exact enumeration produced this. When ``False`` the
            value is an **upper bound** on the true radius.
        witness: A knowledge set that achieves the breakdown at ``radius``, kept
            so the result can be inspected and reproduced. ``None`` when
            ``radius`` is :data:`UNREACHED`.
        n_evaluated: How many knowledge sets were tested.
        metric: Distance metric the radius is measured in.
        budget_exhausted: Whether the search stopped on budget rather than on
            finding a witness -- the flag that separates "safe up to here" from
            "we ran out of time".
    """

    radius: int
    exact: bool
    witness: BackgroundKnowledge | None
    n_evaluated: int
    metric: str
    budget_exhausted: bool


def delta_valid(
    cpdag: MPDAG,
    true_graph: MPDAG,
    treatment: Node,
    outcome: Node,
    distance: str = "orientation_flips",
    *,
    method: str = "exact",
    max_delta: int | None = None,
    max_evaluations: int = 100_000,
) -> int:
    """Smallest perturbation radius at which ``O*`` stops being a valid adjustment set.

    Args:
        cpdag: The CPDAG discovery produced.
        true_graph: The ground-truth DAG; must be a consistent extension of
            ``cpdag``.
        treatment: Treatment node.
        outcome: Outcome node.
        distance: Metric the radius is measured in.
        method: ``"exact"`` or ``"search"``. ``"search"`` returns an upper bound.
        max_delta: Stop searching above this radius; ``None`` uses the number of
            undirected edges in ``cpdag``.
        max_evaluations: Budget on knowledge sets tested.

    Returns:
        The radius, or :data:`UNREACHED` if none was found within budget.

    Raises:
        ValueError: If ``true_graph`` is not a consistent extension of
            ``cpdag``, or if no valid adjustment set exists for ``(treatment,
            outcome)`` even under the truth -- there is then nothing to break.
    """
    raise NotImplementedError


def delta_opt(
    cpdag: MPDAG,
    true_graph: MPDAG,
    treatment: Node,
    outcome: Node,
    scm: Any,
    distance: str = "orientation_flips",
    *,
    method: str = "exact",
    max_delta: int | None = None,
    max_evaluations: int = 100_000,
) -> int:
    """Smallest perturbation radius at which ``O*`` stops being optimal.

    Unlike :func:`delta_valid`, this needs an SCM: optimality is a statement
    about asymptotic variance, and variance depends on the edge coefficients and
    noise scales, not on the graph alone. Two SCMs on the same graph can have
    different ``delta_opt``, and the sweeps therefore report it as a
    distribution over sampled SCMs rather than as a single number per graph.

    Args:
        cpdag: The CPDAG discovery produced.
        true_graph: The ground-truth DAG.
        treatment: Treatment node.
        outcome: Outcome node.
        scm: The structural causal model supplying the variances.
        distance: Metric the radius is measured in.
        method: ``"exact"`` or ``"search"``. ``"search"`` returns an upper bound.
        max_delta: Stop searching above this radius.
        max_evaluations: Budget on knowledge sets tested.

    Returns:
        The radius, or :data:`UNREACHED` if none was found within budget.

    Raises:
        ValueError: As for :func:`delta_valid`.
    """
    raise NotImplementedError


def radius_report(
    cpdag: MPDAG,
    true_graph: MPDAG,
    treatment: Node,
    outcome: Node,
    scm: Any | None = None,
    distance: str = "orientation_flips",
    *,
    method: str = "exact",
    max_delta: int | None = None,
    max_evaluations: int = 100_000,
) -> tuple[RadiusReport, RadiusReport | None]:
    """Compute both radii in one pass, sharing the enumeration.

    The two searches walk the same ball, so computing them together roughly
    halves the work relative to two independent calls.

    Args:
        cpdag: The CPDAG discovery produced.
        true_graph: The ground-truth DAG.
        treatment: Treatment node.
        outcome: Outcome node.
        scm: SCM for the optimality test; when ``None``, only ``delta_valid`` is
            computed and the second element of the return is ``None``.
        distance: Metric the radii are measured in.
        method: ``"exact"`` or ``"search"``.
        max_delta: Stop searching above this radius.
        max_evaluations: Shared budget.

    Returns:
        ``(valid_report, opt_report)``.
    """
    raise NotImplementedError


def verify_ordering(valid: RadiusReport, opt: RadiusReport) -> bool:
    """Check the conjectured ordering ``delta_opt <= delta_valid``.

    Args:
        valid: Report from :func:`delta_valid`.
        opt: Report from :func:`delta_opt`.

    Returns:
        ``True`` if the ordering holds, or if either radius is
        :data:`UNREACHED`, in which case the comparison is vacuous.

    Note:
        A ``False`` here on an **exact** pair of reports is a counterexample to
        the conjecture and is the single most valuable output this codebase can
        produce. The sweeps must persist the full witness for any such case
        rather than aggregating it into a rate; see
        ``configs/experiment/synthetic_sweep.yaml``. On a heuristic pair it
        means nothing, since both values are upper bounds.
    """
    raise NotImplementedError


def radius_distribution(
    cpdag: MPDAG,
    true_graph: MPDAG,
    treatment: Node,
    outcome: Node,
    delta: int,
    scm: Any | None = None,
    distance: str = "orientation_flips",
    *,
    limit: int | None = None,
) -> dict[str, float]:
    """Fraction of the radius-``delta`` ball that breaks each property.

    The typical-case companion to the worst-case radii: the radii say the ball
    at ``delta`` *contains* a breaking knowledge set, this says how much of it
    breaks. A radius of 2 where 1% of the ball breaks is a very different
    practical situation from one where 80% does, and only the second is a reason
    to distrust knowledge-informed discovery in general.

    Args:
        cpdag: The CPDAG discovery produced.
        true_graph: The ground-truth DAG.
        treatment: Treatment node.
        outcome: Outcome node.
        delta: The radius to characterise.
        scm: SCM for the optimality fraction; ``None`` omits it.
        distance: Metric the radius is measured in.
        limit: Cap the ball enumeration; the fractions are then estimates over a
            deterministic prefix and must be reported as such.

    Returns:
        Keys ``"invalid_fraction"``, ``"suboptimal_fraction"``,
        ``"benign_fraction"``, ``"n_evaluated"``.
    """
    raise NotImplementedError
