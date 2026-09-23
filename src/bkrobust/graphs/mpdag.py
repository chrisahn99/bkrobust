"""MPDAG representation: adjacency and orientation storage.

A *maximally oriented PDAG* (MPDAG) is a partially directed graph that is the
Meek closure of a CPDAG plus background knowledge. It has directed edges
(``a -> b``, orientation determined by data or by knowledge) and undirected
edges (``a - b``, orientation not determined either way). A DAG and a CPDAG are
both MPDAGs -- a DAG has no undirected edges, a CPDAG is the closure of the
empty knowledge set -- so this one container serves for all three, and callers
distinguish them by predicate rather than by type.

Conventions used throughout the package:

* Nodes are hashable labels (:data:`Node`); string labels are canonical.
* A directed edge is the ordered pair ``(tail, head)`` meaning ``tail -> head``.
* An undirected edge is stored once, under a canonical ordering of its
  endpoints, so ``a - b`` and ``b - a`` are the same object.
* No graph in this package has both ``a -> b`` and ``a - b``: an edge is
  either oriented or not.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from typing import TYPE_CHECKING, TypeAlias

if TYPE_CHECKING:  # pragma: no cover
    import networkx as nx

#: A node label. Strings are canonical; any hashable is accepted.
Node: TypeAlias = str

#: An ordered pair of nodes. Read as ``tail -> head`` for directed edges and as
#: a canonically ordered endpoint pair for undirected edges.
Edge: TypeAlias = tuple[Node, Node]


class MPDAG:
    """A partially directed acyclic graph with directed and undirected edges.

    Args:
        nodes: The vertex set. Order is retained for deterministic iteration
            and for the row/column order of any matrix encoding.
        directed: Directed edges as ``(tail, head)`` pairs.
        undirected: Undirected edges as unordered endpoint pairs; stored
            canonically so that insertion order does not matter.

    Raises:
        ValueError: If an edge references a node not in ``nodes``, if the same
            pair appears both directed and undirected, or if the directed
            subgraph already contains a cycle.
    """

    def __init__(
        self,
        nodes: Iterable[Node],
        directed: Iterable[Edge] = (),
        undirected: Iterable[Edge] = (),
    ) -> None:
        raise NotImplementedError

    # --- construction ------------------------------------------------------

    @classmethod
    def from_cpdag(cls, cpdag: MPDAG | nx.DiGraph) -> MPDAG:
        """Build an MPDAG from a CPDAG.

        Accepts either an existing :class:`MPDAG` (copied) or a NetworkX
        directed graph in the common CPDAG encoding, where an undirected edge
        is represented by the pair of arcs ``a -> b`` and ``b -> a``.

        Args:
            cpdag: The CPDAG to convert.

        Returns:
            An MPDAG with the same skeleton, mutually-directed arc pairs folded
            into single undirected edges.

        Raises:
            ValueError: If ``cpdag`` is not a valid CPDAG encoding.
        """
        raise NotImplementedError

    @classmethod
    def from_dag(cls, dag: MPDAG | nx.DiGraph) -> MPDAG:
        """Build an MPDAG with no undirected edges from a DAG.

        Args:
            dag: An acyclic directed graph.

        Returns:
            The DAG viewed as a fully oriented MPDAG.

        Raises:
            ValueError: If ``dag`` contains a cycle.
        """
        raise NotImplementedError

    def copy(self) -> MPDAG:
        """Return an independent copy; mutating it must not touch ``self``."""
        raise NotImplementedError

    # --- mutation ----------------------------------------------------------

    def orient(self, tail: Node, head: Node) -> None:
        """Orient the edge between ``tail`` and ``head`` as ``tail -> head``.

        Orients in place and does **not** apply Meek closure; call
        :func:`bkrobust.graphs.meek.meek_closure` for that.

        Args:
            tail: Source endpoint.
            head: Target endpoint.

        Raises:
            KeyError: If there is no edge between the two nodes.
            ValueError: If the edge is already oriented the other way, or if
                the orientation would create a directed cycle.
        """
        raise NotImplementedError

    def unorient(self, a: Node, b: Node) -> None:
        """Replace the directed edge between ``a`` and ``b`` with an undirected one.

        Used by the perturbation samplers when constructing a knowledge set
        relative to a partially oriented graph.

        Raises:
            KeyError: If there is no directed edge between the two nodes.
        """
        raise NotImplementedError

    # --- predicates --------------------------------------------------------

    def is_consistent(self) -> bool:
        """Whether this graph is a well-formed MPDAG.

        Checks that the directed subgraph is acyclic and that the graph
        contains no partially directed cycle -- both conditions must hold for
        the graph to admit at least one consistent DAG extension.

        Returns:
            ``True`` if the graph is a well-formed MPDAG.
        """
        raise NotImplementedError

    def is_dag(self) -> bool:
        """Whether every edge is oriented."""
        raise NotImplementedError

    def has_edge(self, a: Node, b: Node) -> bool:
        """Whether ``a`` and ``b`` are adjacent, in either orientation state."""
        raise NotImplementedError

    def is_directed_edge(self, tail: Node, head: Node) -> bool:
        """Whether the graph contains exactly the directed edge ``tail -> head``."""
        raise NotImplementedError

    def is_undirected_edge(self, a: Node, b: Node) -> bool:
        """Whether the edge between ``a`` and ``b`` exists and is unoriented."""
        raise NotImplementedError

    # --- accessors ---------------------------------------------------------

    @property
    def nodes(self) -> tuple[Node, ...]:
        """The vertex set in a fixed, deterministic order."""
        raise NotImplementedError

    @property
    def directed_edges(self) -> frozenset[Edge]:
        """Directed edges as ``(tail, head)`` pairs."""
        raise NotImplementedError

    @property
    def undirected_edges(self) -> frozenset[Edge]:
        """Undirected edges, each stored once under the canonical ordering."""
        raise NotImplementedError

    def parents(self, node: Node) -> frozenset[Node]:
        """Nodes with a directed edge into ``node``."""
        raise NotImplementedError

    def children(self, node: Node) -> frozenset[Node]:
        """Nodes with a directed edge out of ``node``."""
        raise NotImplementedError

    def neighbors(self, node: Node) -> frozenset[Node]:
        """Nodes joined to ``node`` by an *undirected* edge."""
        raise NotImplementedError

    def adjacent(self, node: Node) -> frozenset[Node]:
        """Every node adjacent to ``node``, regardless of orientation."""
        raise NotImplementedError

    def ancestors(self, nodes: Node | Iterable[Node]) -> frozenset[Node]:
        """Nodes reachable by a directed path into ``nodes`` (exclusive of ``nodes``)."""
        raise NotImplementedError

    def descendants(self, nodes: Node | Iterable[Node]) -> frozenset[Node]:
        """Nodes reachable by a directed path from ``nodes`` (exclusive of ``nodes``)."""
        raise NotImplementedError

    def possible_descendants(self, nodes: Node | Iterable[Node]) -> frozenset[Node]:
        """Nodes reachable by a *possibly* directed path from ``nodes``.

        A possibly directed path may traverse undirected edges in the forward
        direction. Needed by the adjustment criterion, where the forbidden set
        is defined in terms of possible descendants rather than descendants.
        """
        raise NotImplementedError

    # --- extensions --------------------------------------------------------

    def to_dag(self) -> MPDAG:
        """Return one consistent DAG extension of this MPDAG.

        Any extension will do; the choice must be deterministic given the graph
        so that experiments are reproducible.

        Raises:
            ValueError: If no consistent extension exists.
        """
        raise NotImplementedError

    def enumerate_dags(self, limit: int | None = None) -> Iterator[MPDAG]:
        """Yield the consistent DAG extensions of this MPDAG.

        The count is exponential in the number of undirected edges; ``limit``
        caps it. Used by the exact path of :mod:`bkrobust.theory.radius`, which
        needs to quantify over the whole Markov equivalence class.

        Args:
            limit: Stop after yielding this many extensions; ``None`` yields all.

        Yields:
            Fully oriented MPDAGs.
        """
        raise NotImplementedError

    # --- interop -----------------------------------------------------------

    def to_networkx(self) -> nx.DiGraph:
        """Export in the standard CPDAG/MPDAG NetworkX encoding.

        Undirected edges become mutual arc pairs (``a -> b`` and ``b -> a``),
        which is the convention used by the common discovery libraries.

        Returns:
            A ``networkx.DiGraph``.
        """
        raise NotImplementedError

    def to_adjacency_matrix(self) -> object:
        """Export as an integer adjacency matrix over :attr:`nodes` order.

        Encoding: ``M[i, j] == 1`` for ``i -> j``; ``M[i, j] == M[j, i] == 1``
        for an undirected edge; ``0`` otherwise.

        Returns:
            A ``numpy.ndarray`` of shape ``(n_nodes, n_nodes)`` and integer dtype.
        """
        raise NotImplementedError

    # --- dunders -----------------------------------------------------------

    def __eq__(self, other: object) -> bool:
        raise NotImplementedError

    def __hash__(self) -> int:
        raise NotImplementedError

    def __repr__(self) -> str:
        raise NotImplementedError
