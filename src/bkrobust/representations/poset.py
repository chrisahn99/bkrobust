"""Build the constraint poset from a body of background knowledge.

Constraints are partially ordered by entailment on a fixed CPDAG: ``c1 <= c2``
when imposing ``c1`` forces ``c2`` through Meek's rules. The order depends on
the CPDAG, not on the constraints alone -- the same pair of assertions can be
comparable on one graph and incomparable on another, because entailment runs
through the closure.

That dependence is what makes the poset informative and what makes it expensive:
building it requires a closure computation per pair, so :func:`build_poset`
scales quadratically in the number of constraints and the ``max_elements`` cap
in the config exists to keep it tractable.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    import networkx as nx
    import numpy as np

    from bkrobust.graphs.mpdag import MPDAG, Node
    from bkrobust.knowledge.base import BackgroundKnowledge


@dataclass(frozen=True)
class ConstraintPoset:
    """A partial order over constraints, relative to one CPDAG.

    Attributes:
        elements: The constraints, as ``(kind, (a, b))`` pairs, in a fixed order
            that indexes the relation matrix.
        relation: Boolean matrix; ``relation[i, j]`` is ``True`` when element
            ``i`` entails element ``j``.
        cpdag_signature: Hash of the CPDAG the order was built against, so a
            poset can never be silently reused on a different graph.
        cover_relation: The transitive reduction -- the Hasse diagram edges --
            which is what the embedding losses train on.
    """

    elements: tuple[tuple[str, tuple[Node, Node]], ...]
    relation: np.ndarray
    cpdag_signature: str
    cover_relation: frozenset[tuple[int, int]]

    def hasse_diagram(self) -> nx.DiGraph:
        """Return the Hasse diagram as a NetworkX graph, for plotting and probes."""
        raise NotImplementedError

    def rank(self, index: int) -> int:
        """Longest chain below element ``index``; a cheap scalar summary of depth."""
        raise NotImplementedError

    def __len__(self) -> int:
        """Number of elements."""
        raise NotImplementedError


def build_poset(
    cpdag: MPDAG,
    constraints: Sequence[tuple[str, tuple[Node, Node]]] | None = None,
    *,
    include_cascade: bool = True,
    max_elements: int = 5000,
) -> ConstraintPoset:
    """Construct the entailment poset over constraints on ``cpdag``.

    Args:
        cpdag: The CPDAG entailment is computed relative to.
        constraints: The constraints to order; ``None`` uses every constraint
            that is consistent with ``cpdag``, which is the full orientable pool.
        include_cascade: Include the orientations Meek's rules force as poset
            elements in their own right. Without them the order records what was
            asserted; with them it records what was entailed, and the second is
            the object the hypothesis is actually about.
        max_elements: Cap on poset size.

    Returns:
        A :class:`ConstraintPoset`.

    Raises:
        ValueError: If the constraint set exceeds ``max_elements``.
    """
    raise NotImplementedError


def poset_from_knowledge(
    cpdag: MPDAG,
    bk: BackgroundKnowledge,
    *,
    include_cascade: bool = True,
) -> ConstraintPoset:
    """Build the poset induced by one knowledge set's own constraints.

    The per-sample object the embeddings consume: one poset per sampled
    knowledge set, labelled with that set's radius and truth value.
    """
    raise NotImplementedError


def entails(
    cpdag: MPDAG,
    c1: tuple[str, tuple[Node, Node]],
    c2: tuple[str, tuple[Node, Node]],
) -> bool:
    """Whether imposing ``c1`` on ``cpdag`` forces ``c2``.

    The primitive the whole poset is built from.

    Raises:
        ValueError: If either constraint is inconsistent with ``cpdag``.
    """
    raise NotImplementedError
