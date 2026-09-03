"""The shared graph container for the demonstration.

Deliberately minimal: a container plus accessors, with the invariants stated and
checked. All algorithms -- Meek's rules, validity, enumeration, evaluation --
live in sibling modules and operate on this type.

Conventions, fixed here so every module agrees:

* Nodes are strings.
* A directed edge is the ordered pair ``(tail, head)`` meaning ``tail -> head``.
* An undirected edge is stored once, under a canonical (sorted) ordering of its
  endpoints, so ``a - b`` and ``b - a`` are the same object.
* No pair is both directed and undirected.
* The container is immutable in practice: mutators return new graphs.

Written to run on Python 3.9 as well as the repository's 3.11 target -- the
demonstration had to execute in the available interpreter. No runtime-only 3.10+
syntax is used.
"""

from __future__ import annotations

from collections.abc import Iterable

Node = str
Edge = tuple[str, str]


def canon(a: Node, b: Node) -> Edge:
    """Canonical ordering of an undirected edge's endpoints."""
    return (a, b) if a <= b else (b, a)


class MPDAG:
    """A partially directed graph: directed edges plus undirected edges.

    A DAG (no undirected edges) and a CPDAG are both representable, so one
    container serves for all three and callers distinguish by predicate.

    Args:
        nodes: The vertex set. Iteration order is sorted for determinism.
        directed: ``(tail, head)`` pairs.
        undirected: Endpoint pairs; stored canonically.

    Raises:
        ValueError: On an unknown node, a self-loop, or a pair that is both
            directed and undirected.
    """

    __slots__ = ("_directed", "_hash", "_nodes", "_undirected")

    def __init__(
        self,
        nodes: Iterable[Node],
        directed: Iterable[Edge] = (),
        undirected: Iterable[Edge] = (),
    ) -> None:
        self._nodes: tuple[Node, ...] = tuple(sorted(set(nodes)))
        node_set = set(self._nodes)

        d: set[Edge] = set()
        for a, b in directed:
            if a not in node_set or b not in node_set:
                raise ValueError(f"directed edge references unknown node: {(a, b)}")
            if a == b:
                raise ValueError(f"self-loop: {a}")
            d.add((a, b))

        u: set[Edge] = set()
        for a, b in undirected:
            if a not in node_set or b not in node_set:
                raise ValueError(f"undirected edge references unknown node: {(a, b)}")
            if a == b:
                raise ValueError(f"self-loop: {a}")
            u.add(canon(a, b))

        for a, b in d:
            if canon(a, b) in u:
                raise ValueError(f"pair is both directed and undirected: {(a, b)}")
        for a, b in d:
            if (b, a) in d:
                raise ValueError(f"mutually directed pair (encode as undirected): {(a, b)}")

        self._directed: frozenset[Edge] = frozenset(d)
        self._undirected: frozenset[Edge] = frozenset(u)
        self._hash: int | None = None

    # --- accessors ---------------------------------------------------------

    @property
    def nodes(self) -> tuple[Node, ...]:
        """Vertices, sorted."""
        return self._nodes

    @property
    def directed_edges(self) -> frozenset[Edge]:
        """Directed edges as ``(tail, head)``."""
        return self._directed

    @property
    def undirected_edges(self) -> frozenset[Edge]:
        """Undirected edges, canonically ordered, each stored once."""
        return self._undirected

    def has_edge(self, a: Node, b: Node) -> bool:
        """Whether ``a`` and ``b`` are adjacent in either orientation state."""
        return (
            (a, b) in self._directed or (b, a) in self._directed or canon(a, b) in self._undirected
        )

    def is_directed_edge(self, a: Node, b: Node) -> bool:
        """Whether the graph contains exactly ``a -> b``."""
        return (a, b) in self._directed

    def is_undirected_edge(self, a: Node, b: Node) -> bool:
        """Whether ``a`` and ``b`` are joined by an unoriented edge."""
        return canon(a, b) in self._undirected

    def parents(self, node: Node) -> set[Node]:
        """Nodes with a directed edge into ``node``."""
        return {a for a, b in self._directed if b == node}

    def children(self, node: Node) -> set[Node]:
        """Nodes with a directed edge out of ``node``."""
        return {b for a, b in self._directed if a == node}

    def neighbors(self, node: Node) -> set[Node]:
        """Nodes joined to ``node`` by an *undirected* edge."""
        out: set[Node] = set()
        for a, b in self._undirected:
            if a == node:
                out.add(b)
            elif b == node:
                out.add(a)
        return out

    def adjacent(self, node: Node) -> set[Node]:
        """Every node adjacent to ``node``, regardless of orientation."""
        return self.parents(node) | self.children(node) | self.neighbors(node)

    def skeleton(self) -> frozenset[Edge]:
        """Adjacencies with orientation forgotten, canonically ordered."""
        return frozenset([canon(a, b) for a, b in self._directed] + list(self._undirected))

    # --- derived graphs ----------------------------------------------------

    def oriented(self, tail: Node, head: Node) -> MPDAG:
        """Return a copy with the ``tail -- head`` edge oriented ``tail -> head``.

        Raises:
            KeyError: If no undirected edge joins the two nodes.
        """
        key = canon(tail, head)
        if key not in self._undirected:
            raise KeyError(f"no undirected edge between {tail} and {head}")
        return MPDAG(
            self._nodes,
            set(self._directed) | {(tail, head)},
            set(self._undirected) - {key},
        )

    def unoriented(self, a: Node, b: Node) -> MPDAG:
        """Return a copy with the directed edge between ``a`` and ``b`` made undirected.

        Raises:
            KeyError: If no directed edge joins the two nodes.
        """
        if (a, b) in self._directed:
            drop = (a, b)
        elif (b, a) in self._directed:
            drop = (b, a)
        else:
            raise KeyError(f"no directed edge between {a} and {b}")
        return MPDAG(
            self._nodes,
            set(self._directed) - {drop},
            set(self._undirected) | {canon(a, b)},
        )

    def copy(self) -> MPDAG:
        """An independent copy."""
        return MPDAG(self._nodes, self._directed, self._undirected)

    # --- predicates --------------------------------------------------------

    def is_dag(self) -> bool:
        """Whether every edge is oriented (and the graph is acyclic)."""
        return not self._undirected and self.is_acyclic()

    def is_acyclic(self) -> bool:
        """Whether the *directed* subgraph contains no directed cycle."""
        colour: dict[Node, int] = {n: 0 for n in self._nodes}
        adj: dict[Node, list[Node]] = {n: [] for n in self._nodes}
        for a, b in self._directed:
            adj[a].append(b)

        def visit(start: Node) -> bool:
            stack: list[tuple[Node, int]] = [(start, 0)]
            while stack:
                node, state = stack.pop()
                if state == 0:
                    if colour[node] == 2:
                        continue
                    if colour[node] == 1:
                        return False
                    colour[node] = 1
                    stack.append((node, 1))
                    for nxt in adj[node]:
                        if colour[nxt] == 1:
                            return False
                        if colour[nxt] == 0:
                            stack.append((nxt, 0))
                else:
                    colour[node] = 2
            return True

        for n in self._nodes:
            if colour[n] == 0 and not visit(n):
                return False
        return True

    # --- reachability ------------------------------------------------------

    def descendants(self, node: Node) -> set[Node]:
        """Nodes reachable from ``node`` by a directed path (exclusive)."""
        seen: set[Node] = set()
        stack = list(self.children(node))
        while stack:
            cur = stack.pop()
            if cur in seen:
                continue
            seen.add(cur)
            stack.extend(self.children(cur))
        return seen

    def ancestors(self, node: Node) -> set[Node]:
        """Nodes with a directed path into ``node`` (exclusive)."""
        seen: set[Node] = set()
        stack = list(self.parents(node))
        while stack:
            cur = stack.pop()
            if cur in seen:
                continue
            seen.add(cur)
            stack.extend(self.parents(cur))
        return seen

    # --- identity ----------------------------------------------------------

    def key(self) -> tuple[tuple[Node, ...], frozenset[Edge], frozenset[Edge]]:
        """Hashable canonical identity, used as a dict key across modules."""
        return (self._nodes, self._directed, self._undirected)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, MPDAG):
            return NotImplemented
        return self.key() == other.key()

    def __hash__(self) -> int:
        if self._hash is None:
            self._hash = hash(self.key())
        return self._hash

    def __repr__(self) -> str:
        d = ", ".join(f"{a}->{b}" for a, b in sorted(self._directed))
        u = ", ".join(f"{a}-{b}" for a, b in sorted(self._undirected))
        return f"MPDAG(directed=[{d}], undirected=[{u}])"

    def edge_string(self) -> str:
        """Compact one-line rendering, used in tables and figure annotations."""
        parts = [f"{a}->{b}" for a, b in sorted(self._directed)]
        parts += [f"{a}-{b}" for a, b in sorted(self._undirected)]
        return " ".join(parts)


def induced_subgraph(graph: MPDAG, nodes: Iterable[Node]) -> MPDAG:
    """The subgraph induced on ``nodes``, keeping orientation states."""
    keep = set(nodes)
    return MPDAG(
        keep,
        [(a, b) for a, b in graph.directed_edges if a in keep and b in keep],
        [(a, b) for a, b in graph.undirected_edges if a in keep and b in keep],
    )


def undirected_components(graph: MPDAG) -> list[frozenset[Node]]:
    """Connected components of the undirected subgraph, singletons omitted.

    Sorted for determinism. These are the *chordal components* whose orientations
    the perturbation space ranges over.
    """
    parent: dict[Node, Node] = {n: n for n in graph.nodes}

    def find(x: Node) -> Node:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in graph.undirected_edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    groups: dict[Node, set[Node]] = {}
    for a, b in graph.undirected_edges:
        groups.setdefault(find(a), set()).update((a, b))
    return sorted(
        (frozenset(g) for g in groups.values()),
        key=lambda s: sorted(s),
    )
