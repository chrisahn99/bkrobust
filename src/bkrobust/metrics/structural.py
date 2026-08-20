"""Structural metrics: validity flips, optimality loss, cascade size.

These sit between the perturbation and the estimation error and explain the
link. A perturbation that changes many orientations but never touches ``O*``
produces no estimation error at all, and without structural metrics that case is
indistinguishable in the results from a perturbation that did nothing.

The three headline rates, all measured over a ball of sampled knowledge sets at
a fixed radius:

``validity_flip_rate``
    Fraction whose ``O*`` is invalid in the true graph. Rises from zero at
    ``delta_valid``.
``optimality_loss_rate``
    Fraction whose ``O*`` is valid but suboptimal. Rises from zero at
    ``delta_opt``, and should rise first if the conjectured ordering holds.
``benign_rate``
    Fraction where ``O*`` is unchanged. Falls as the radius grows; how fast it
    falls is how likely wrong knowledge is to matter at all.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    from bkrobust.graphs.mpdag import MPDAG, Node
    from bkrobust.knowledge.base import BackgroundKnowledge


def validity_flip_rate(
    cpdag: MPDAG,
    true_graph: MPDAG,
    treatment: Node,
    outcome: Node,
    knowledge_sets: Sequence[BackgroundKnowledge],
) -> float:
    """Fraction of knowledge sets whose ``O*`` is invalid in ``true_graph``.

    Raises:
        ValueError: If ``knowledge_sets`` is empty.
    """
    raise NotImplementedError


def optimality_loss_rate(
    cpdag: MPDAG,
    true_graph: MPDAG,
    treatment: Node,
    outcome: Node,
    knowledge_sets: Sequence[BackgroundKnowledge],
    scm: Any,
) -> float:
    """Fraction whose ``O*`` is valid but not optimal.

    Needs an SCM: optimality is a variance statement, not a graph statement.

    Raises:
        ValueError: If ``knowledge_sets`` is empty.
    """
    raise NotImplementedError


def benign_rate(
    cpdag: MPDAG,
    true_graph: MPDAG,
    treatment: Node,
    outcome: Node,
    knowledge_sets: Sequence[BackgroundKnowledge],
) -> float:
    """Fraction whose ``O*`` is unchanged by the perturbation.

    Wrong knowledge that costs nothing. Not a footnote -- if this rate stays
    high across the radii, the practical message of the paper is narrower than
    the theory suggests, and that needs saying.
    """
    raise NotImplementedError


def adjustment_set_distance(z1: set[Node], z2: set[Node]) -> int:
    """Symmetric difference size between two adjustment sets."""
    raise NotImplementedError


def cascade_size(cpdag: MPDAG, bk: BackgroundKnowledge) -> int:
    """Number of orientations forced beyond those directly imposed.

    Thin wrapper over :mod:`bkrobust.knowledge.cascade`, here so that every
    metric an experiment records comes from one namespace.
    """
    raise NotImplementedError


def severity_distribution(
    cpdag: MPDAG,
    true_graph: MPDAG,
    treatment: Node,
    outcome: Node,
    knowledge_sets: Sequence[BackgroundKnowledge],
    scm: Any | None = None,
) -> dict[str, float]:
    """Fraction of the ball falling in each severity class.

    Returns:
        Fractions keyed by the values of
        :class:`~bkrobust.knowledge.taxonomy.Severity`, summing to 1.
    """
    raise NotImplementedError


def structural_summary(
    cpdag: MPDAG,
    true_graph: MPDAG,
    treatment: Node,
    outcome: Node,
    knowledge_sets: Sequence[BackgroundKnowledge],
    scm: Any | None = None,
) -> dict[str, float]:
    """Every structural metric at once, for one experimental cell.

    Returns:
        Keys ``"validity_flip_rate"``, ``"optimality_loss_rate"``,
        ``"benign_rate"``, ``"mean_cascade_size"``,
        ``"mean_adjustment_set_distance"``, ``"n"``.
    """
    raise NotImplementedError
