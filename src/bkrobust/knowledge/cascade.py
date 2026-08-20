"""Meek-closure amplification: k imposed constraints, m forced orientations.

Background knowledge does not stay where it is put. Impose one orientation and
Meek's rules propagate it: R1 forbids new v-structures, R2 forbids cycles, and
the consequences travel along the undirected components of the CPDAG until they
run out of edges to orient. A single expert assertion can determine a dozen
edges nowhere near the pair it named.

For true knowledge that is a feature -- it is the whole reason to elicit
knowledge at all. For false knowledge it is the mechanism of harm, and it is
why the two breakdown radii can be small even on large graphs: the analyst
supplies one wrong claim, the closure turns it into many wrong orientations,
and one of those lands on the adjustment set.

The quantities this module measures:

* **amplification factor** ``m / k`` -- forced orientations per imposed
  constraint;
* **cascade reach** -- graph distance from the imposed constraint to the
  farthest edge it forces, i.e. how far from the claim the damage can appear;
* **hitting probability** -- the chance a cascade reaches ``O*`` at all, which
  is what converts amplification into estimation error.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    import numpy as np

    from bkrobust.graphs.mpdag import MPDAG, Edge, Node
    from bkrobust.knowledge.base import BackgroundKnowledge


@dataclass(frozen=True)
class CascadeReport:
    """What imposing one body of knowledge on one CPDAG actually did.

    Attributes:
        n_imposed: Constraints asserted directly (``k``).
        n_forced: Additional edges oriented by the closure (``m``).
        amplification: ``n_forced / n_imposed``; ``0.0`` when ``n_imposed`` is 0.
        forced_edges: The cascade-oriented edges themselves.
        reach: Maximum graph distance from an imposed constraint to a forced
            edge; ``0`` when nothing cascaded.
        hits_optimal_set: Whether any forced edge changed ``O*``. ``None`` when
            no treatment/outcome pair was supplied.
        rule_counts: How many orientations each of R1-R4 contributed, keyed
            ``"R1"``-``"R4"``. Diagnoses *which* rule does the damage, which R4
            in particular is expected to, since it only fires once background
            knowledge is present.
    """

    n_imposed: int
    n_forced: int
    amplification: float
    forced_edges: frozenset[Edge]
    reach: int
    hits_optimal_set: bool | None
    rule_counts: dict[str, int]


def cascade_report(
    cpdag: MPDAG,
    bk: BackgroundKnowledge,
    *,
    treatment: Node | None = None,
    outcome: Node | None = None,
) -> CascadeReport:
    """Measure the amplification of imposing ``bk`` on ``cpdag``.

    Args:
        cpdag: The CPDAG the knowledge is imposed on.
        bk: The knowledge to impose.
        treatment: Treatment node; supply with ``outcome`` to populate
            ``hits_optimal_set``.
        outcome: Outcome node.

    Returns:
        A :class:`CascadeReport`.

    Raises:
        ValueError: If ``bk`` is inconsistent with ``cpdag``.
    """
    raise NotImplementedError


def amplification_factor(cpdag: MPDAG, bk: BackgroundKnowledge) -> float:
    """Forced orientations per imposed constraint.

    Returns:
        ``m / k``, or ``0.0`` when ``bk`` asserts nothing.
    """
    raise NotImplementedError


def cascade_reach(cpdag: MPDAG, bk: BackgroundKnowledge) -> int:
    """Graph distance from the imposed constraints to the farthest forced edge.

    Answers how far from the stated claim a wrong orientation can surface, and
    hence how implausible it is that an analyst would notice.
    """
    raise NotImplementedError


def hits_optimal_set(
    cpdag: MPDAG,
    bk: BackgroundKnowledge,
    treatment: Node,
    outcome: Node,
) -> bool:
    """Whether imposing ``bk`` changes ``O*`` relative to the unaided CPDAG.

    The bridge from structure to estimation: a cascade that misses ``O*``
    costs nothing measurable, however large it is.
    """
    raise NotImplementedError


def expected_amplification(
    cpdag: MPDAG,
    k: int,
    rng: np.random.Generator,
    n_draws: int = 100,
) -> tuple[float, float]:
    """Monte-Carlo estimate of amplification for ``k`` random consistent constraints.

    Characterises a CPDAG's *susceptibility* independently of any particular
    knowledge set, so graph families can be compared -- the expectation is that
    scale-free graphs amplify harder than Erdos-Renyi at matched density,
    because a hub orientation cascades through the whole hub neighbourhood.

    Args:
        cpdag: The CPDAG to characterise.
        k: Number of constraints imposed per draw.
        rng: Seeded generator.
        n_draws: Number of Monte-Carlo draws.

    Returns:
        ``(mean, standard_error)`` of the amplification factor.
    """
    raise NotImplementedError
