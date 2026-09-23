"""The knowledge-sensitivity certificate -- the practitioner-facing deliverable.

Everything else in :mod:`bkrobust.theory` needs the true DAG. A practitioner
does not have one; if they did, they would not be running discovery. So the
theory as stated is unusable by the people whose problem it describes, and the
certificate is the bridge.

The move is the one sensitivity analysis always makes: stop asking "is my
knowledge true" -- unanswerable -- and ask "how wrong would it have to be
before my conclusion changed". That question is answerable from the CPDAG and
the asserted knowledge alone, because both radii are minima over a ball defined
by the *observed* CPDAG, and nothing in the definition requires knowing which
member of the ball is true.

What a certificate reports:

* the radii, computed against the analyst's own knowledge as the reference
  point rather than against an unknown truth;
* the worst-case bias if the knowledge is wrong by that much, in units of the
  estimate's own standard error;
* which specific constraints are load-bearing -- the ones whose removal changes
  ``O*``. Usually a small subset, and usually not the ones the analyst was most
  worried about, which is what makes the list worth printing;
* a verdict in plain language.

A certificate is not a guarantee that the knowledge is right. It says how much
of the conclusion rests on it, which is the honest thing to say and, if the
analyst is going to publish an effect estimate, the thing they should have to
say.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    from bkrobust.graphs.mpdag import MPDAG, Node
    from bkrobust.knowledge.base import BackgroundKnowledge


@dataclass(frozen=True)
class Certificate:
    """How much of a conclusion rests on background knowledge that was never verified.

    Attributes:
        treatment: Treatment node.
        outcome: Outcome node.
        adjustment_set: ``O*`` under the knowledge as supplied -- the set the
            analyst would actually adjust on.
        delta_opt: Radius at which optimality is lost, relative to the supplied
            knowledge. :data:`bkrobust.theory.radius.UNREACHED` if not found.
        delta_valid: Radius at which validity is lost.
        exact: Whether the radii are exact or upper bounds from the heuristic
            search. Printed in the report, because an upper bound that reads as
            exact would overstate the safety margin.
        worst_case_bias: Largest bias achievable by knowledge one step outside
            the safe radius, in units of the estimate's standard error;
            ``None`` if no estimate was supplied.
        critical_constraints: The load-bearing constraints -- those whose
            removal changes ``adjustment_set``. These are what to re-examine.
        n_constraints: Total constraints supplied, for context on how many of
            them turned out to matter.
        cascade_amplification: Forced orientations per imposed constraint. A
            high value means the supplied knowledge reaches much further than it
            appears to, which is itself a reason for caution.
        verdict: Machine-readable summary -- ``"robust"``, ``"fragile"``,
            ``"broken"``, or ``"undetermined"``.
        metric: Distance metric the radii are measured in.
        notes: Caveats attached during computation -- budget exhaustion,
            heuristic fallbacks, assumptions on sensitivity parameters.
    """

    treatment: Node
    outcome: Node
    adjustment_set: frozenset[Node]
    delta_opt: int
    delta_valid: int
    exact: bool
    worst_case_bias: float | None
    critical_constraints: tuple[tuple[str, tuple[Node, Node]], ...]
    n_constraints: int
    cascade_amplification: float
    verdict: str
    metric: str
    notes: tuple[str, ...] = field(default_factory=tuple)

    def __str__(self) -> str:
        """Render a report an analyst can paste into an appendix.

        Sections: what was assumed, what it bought, how far it can be wrong
        before the conclusion moves, which constraints are doing the work, and
        the verdict. Plain text, fixed width, no colour -- it goes into papers
        and logs, not terminals.

        Heuristic radii must be marked as upper bounds in the rendered text,
        not only in :attr:`exact`.
        """
        raise NotImplementedError

    def to_dict(self) -> dict[str, Any]:
        """JSON-serialisable form, for the run manifest."""
        raise NotImplementedError


def certify(
    graph: MPDAG,
    bk: BackgroundKnowledge,
    treatment: Node,
    outcome: Node,
    scm: Any | None = None,
    *,
    distance: str = "orientation_flips",
    method: str = "exact",
    r2_treatment: float | None = None,
    r2_outcome: float | None = None,
    se_tau: float | None = None,
) -> Certificate:
    """Certify how sensitive a conclusion is to the background knowledge behind it.

    Args:
        graph: The CPDAG discovery produced. May also be an already-closed
            MPDAG, in which case ``bk`` should be the knowledge that produced it.
        bk: The knowledge the analyst supplied.
        treatment: Treatment node.
        outcome: Outcome node.
        scm: Optional fitted SCM. Without it the optimality radius cannot be
            computed and ``delta_opt`` is :data:`bkrobust.theory.radius.UNREACHED`.
        distance: Metric the radii are measured in.
        method: ``"exact"`` or ``"search"``.
        r2_treatment: Assumed treatment-side partial :math:`R^2` for the bias
            bound.
        r2_outcome: Assumed outcome-side partial :math:`R^2`.
        se_tau: Standard error of the estimate, so the bias can be reported in
            standard-error units.

    Returns:
        A :class:`Certificate`.

    Raises:
        ValueError: If ``bk`` is inconsistent with ``graph`` -- the analyst's
            own tooling would already have rejected it -- or if the effect is
            not identifiable by adjustment even under the supplied knowledge.
    """
    raise NotImplementedError


def critical_constraints(
    graph: MPDAG,
    bk: BackgroundKnowledge,
    treatment: Node,
    outcome: Node,
) -> list[tuple[str, tuple[Node, Node]]]:
    """Return the constraints whose removal changes ``O*``.

    Leave-one-out over the supplied constraints. Cheap, and the most directly
    actionable line in the certificate: it turns "please re-examine your
    assumptions" into a specific short list.

    Args:
        graph: The CPDAG the knowledge is imposed on.
        bk: The supplied knowledge.
        treatment: Treatment node.
        outcome: Outcome node.

    Returns:
        Load-bearing constraints as ``(kind, (a, b))``, in a deterministic order.
    """
    raise NotImplementedError


def verdict_from_radii(
    delta_opt: int,
    delta_valid: int,
    *,
    fragile_below: int = 2,
) -> str:
    """Map the two radii to a plain-language verdict.

    Args:
        delta_opt: The optimality radius.
        delta_valid: The validity radius.
        fragile_below: Radii strictly below this count as fragile. Config-driven
            rather than fixed, because what counts as fragile depends on how
            much the analyst trusts the elicitation -- a threshold of 2 says a
            single wrong claim is tolerable and two are not.

    Returns:
        ``"robust"``, ``"fragile"``, ``"broken"`` (``delta_valid == 0``: the
        knowledge as supplied already invalidates ``O*``), or
        ``"undetermined"`` (a radius was not found within budget).
    """
    raise NotImplementedError


def certificate_table(certificates: Iterable[Certificate]) -> Any:
    """Collect certificates into a DataFrame, one row each.

    For the real-data experiment, where a certificate is produced per bootstrap
    replicate of the discovery step and the spread across replicates is itself
    part of the answer.

    Returns:
        A ``pandas.DataFrame``.
    """
    raise NotImplementedError
