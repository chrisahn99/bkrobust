"""Distances between graphs, and the metric that the radius delta is measured in.

The breakdown radii are only meaningful relative to a stated distance, because
"how far off is this knowledge" has several defensible answers and they do not
agree:

``shd``
    Structural Hamming distance. The default in the discovery literature;
    counts every edge disagreement including additions and deletions.
``orientation_flips``
    Counts only edges present in both graphs but oriented differently. This is
    the natural unit for background knowledge: knowledge does not add or remove
    edges from a CPDAG, it orients them, so the flip count is the number of
    orientation claims the knowledge gets wrong.
``poset``
    Model-oriented distance on the constraint poset: two knowledge sets are
    close when they entail similar orientation closures. This is the one that
    respects the Meek cascade -- two knowledge sets one flip apart in
    ``orientation_flips`` can be far apart here if that flip cascades -- and it
    is the metric that the conditional representation-learning component is
    trying to learn a geometry for.

Experiments must record which metric produced a reported delta;
:func:`get_distance` exists so the choice comes from config and never from a
default buried in a call site.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from bkrobust.graphs.mpdag import MPDAG
    from bkrobust.knowledge.base import BackgroundKnowledge


def structural_hamming_distance(g1: MPDAG, g2: MPDAG) -> int:
    """Structural Hamming distance between two graphs on the same node set.

    Counts, over all node pairs: edges present in one graph and absent in the
    other, and edges present in both with different orientation states.

    Raises:
        ValueError: If the two graphs have different node sets.
    """
    raise NotImplementedError


def orientation_flip_count(g1: MPDAG, g2: MPDAG) -> int:
    """Number of edges oriented one way in ``g1`` and the other way in ``g2``.

    Edges undirected in either graph, and edges absent from either, do not
    count. This is the canonical unit of perturbation radius for this project.

    Raises:
        ValueError: If the two graphs have different node sets.
    """
    raise NotImplementedError


def poset_distance(
    bk1: BackgroundKnowledge,
    bk2: BackgroundKnowledge,
    cpdag: MPDAG,
) -> float:
    """Model-oriented distance between two knowledge sets on a common CPDAG.

    Compares the *closures* the two knowledge sets induce rather than their
    literal constraints, so that constraints which entail one another are near
    and constraints whose cascades diverge are far. Defined as the normalised
    symmetric difference of the induced orientation sets.

    Args:
        bk1: First knowledge set.
        bk2: Second knowledge set.
        cpdag: The CPDAG both are imposed on; the distance is relative to it
            because the same constraint cascades differently on different CPDAGs.

    Returns:
        A value in ``[0, 1]``; zero iff the two induce the same MPDAG.

    Raises:
        ValueError: If either knowledge set is inconsistent with ``cpdag``,
            in which case it induces no MPDAG and the distance is undefined.
    """
    raise NotImplementedError


def knowledge_distance(
    bk1: BackgroundKnowledge,
    bk2: BackgroundKnowledge,
    cpdag: MPDAG,
    metric: str = "orientation_flips",
) -> float:
    """Distance between two knowledge sets under the named metric.

    Args:
        bk1: First knowledge set.
        bk2: Second knowledge set.
        cpdag: The CPDAG both are imposed on.
        metric: ``"shd"``, ``"orientation_flips"`` or ``"poset"``.

    Returns:
        The distance. Integral for the first two metrics, in ``[0, 1]`` for
        ``"poset"``.

    Raises:
        KeyError: If ``metric`` is not a registered metric name.
    """
    raise NotImplementedError


def get_distance(metric: str) -> Callable[..., float]:
    """Look up a distance function by config name.

    Args:
        metric: ``"shd"``, ``"orientation_flips"`` or ``"poset"``.

    Returns:
        The corresponding callable.

    Raises:
        KeyError: If ``metric`` is unknown; the message lists the valid names.
    """
    raise NotImplementedError
