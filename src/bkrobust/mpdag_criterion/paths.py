"""Possibly-causal paths, definite status, and blocking -- the path layer.

Why this layer exists separately from :mod:`bkrobust.mpdag_criterion.criterion`:
the graphical adjustment criterion for an MPDAG is stated in terms of three
notions that are each independently error-prone, and each of which is a
statement about a *path*, not about the graph as a whole:

* whether a path is **possibly causal** (no edge on it points back towards the
  treatment) -- this is what makes "descendant" ambiguous in a partially
  oriented graph and is the basis of ``possible_descendants``;
* whether a node on a path is of **definite status** -- a collider or a
  *definite* non-collider. A node whose role as collider/non-collider differs
  between DAG extensions is of indefinite status, and paths through such nodes
  are deliberately *not* required to be blocked: the criterion's completeness
  proof shows a separate definite-status path always witnesses any real
  failure;
* whether a definite-status path is **blocked** by a set ``Z``.

Splitting them out means each can be unit-tested against hand-worked examples
rather than only through the end-to-end predicate.

Two implementations of the same notions
---------------------------------------
Every notion above is *defined* on a path, and the obvious way to decide a
"for every path" statement is to enumerate the paths. That is what
:func:`simple_paths`, :func:`possibly_causal_paths`,
:func:`unshielded_possibly_causal_paths` and ``_unshielded_reach_by_paths`` do,
and it is
**exponential in the number of vertices**: a dense graph on ``n`` nodes has
``Theta((n-2)!)`` simple paths between any two of them. Measured on this
repository's own dense Erdos-Renyi instances, the enumerating decision took 48
seconds for a single query at ``n = 12`` and blew a 60-second cap at ``n = 14``.

So the module carries a second implementation of the two notions the *decision*
needs, phrased as **reachability over edge states** rather than enumeration --
exactly the move that turns "is every path blocked?" into Bayes-ball, and the
same move :func:`bkrobust.demo.evaluate.is_dseparated` makes for d-separation:

* :func:`unshielded_reachable` replaces the path walk behind
  :func:`possible_descendants`;
* :func:`open_definite_status_non_causal_path` replaces the path walk behind
  condition (c) of the criterion.

A state is an *edge traversal* ``(u, v)`` -- "we arrived at ``v`` from ``u``" --
optionally carrying one bit of history. The identity of ``u`` has to be in the
state, not just the kind of edge traversed: both the unshielded-triple test
(``u`` non-adjacent to the next node) and the definite-non-collider test are
statements about ``u``, ``v`` and the next node together. There are ``O(E)``
states and ``O(V*E)`` transitions, each visited once, so the searches are
polynomial and allocate nothing per path.

The enumerating functions are **kept**, not deleted: they are the reference the
state searches are differentially tested against (see
``tests/criterion/test_criterion.py``), and they are the only way to *name* the
paths for a report. Nothing on the decision path calls them; each one says so in
its own docstring.

Determinism: every iteration over a graph accessor (which return plain ``set``s)
goes through ``sorted``, so path enumeration and search order never depend on
``PYTHONHASHSEED``.

Written to run on Python 3.9 as well as the repository's 3.11 target, matching
the rest of the package: builtin generics and ``X | Y`` unions appear only in
annotations, and ``zip`` is never called with ``strict=`` (3.10+).
"""

from __future__ import annotations

from collections import deque
from collections.abc import Iterable

from bkrobust.demo.graph import MPDAG

Node = str
Path = tuple[Node, ...]

#: How the edge between two nodes is oriented *relative to the direction of
#: travel* ``u -> v``. ``"out"`` is ``u -> v`` (an arrow leaving ``u``),
#: ``"in"`` is ``v -> u`` (an arrow pointing back the way we came), ``"und"``
#: is ``u - v``. A path is non-causal exactly when some step on it is ``"in"``.
_OUT = "out"
_IN = "in"
_UND = "und"


class _Index:
    """Adjacency tables for one MPDAG, built once and reused.

    :class:`~bkrobust.demo.graph.MPDAG` stores edges as flat frozensets, so
    ``parents``/``children``/``neighbors``/``adjacent`` each cost ``O(E)`` per
    call. The state searches call them inside their inner loop, which would turn
    an ``O(V*E)`` search into ``O(V*E^2)``. Materialising the tables once makes
    every accessor ``O(1)`` (or ``O(deg)`` to iterate).

    Attributes:
        adj: Node -> sorted tuple of adjacent nodes. Sorted, so every search
            visits neighbours in a fixed order.
        adj_set: Node -> set of adjacent nodes, for ``O(1)`` adjacency tests.
        children: Node -> set of heads of directed edges out of it.
        neighbors: Node -> set of undirected neighbours.
        possde: Memo for :func:`possible_descendants` of a single node on this
            graph, filled lazily.
    """

    __slots__ = ("adj", "adj_set", "children", "neighbors", "possde")

    def __init__(self, g: MPDAG) -> None:
        adj_set: dict[Node, set[Node]] = {n: set() for n in g.nodes}
        children: dict[Node, set[Node]] = {n: set() for n in g.nodes}
        neighbors: dict[Node, set[Node]] = {n: set() for n in g.nodes}
        for a, b in g.directed_edges:
            children[a].add(b)
            adj_set[a].add(b)
            adj_set[b].add(a)
        for a, b in g.undirected_edges:
            neighbors[a].add(b)
            neighbors[b].add(a)
            adj_set[a].add(b)
            adj_set[b].add(a)
        self.adj_set = adj_set
        self.children = children
        self.neighbors = neighbors
        self.adj: dict[Node, tuple[Node, ...]] = {n: tuple(sorted(adj_set[n])) for n in g.nodes}
        self.possde: dict[Node, set[Node]] = {}

    def step_kind(self, u: Node, v: Node) -> str:
        """How the ``u -- v`` edge is oriented relative to travelling ``u`` to ``v``.

        Args:
            u: The node the step leaves.
            v: The node the step enters. Assumed adjacent to ``u``.

        Returns:
            :data:`_OUT`, :data:`_IN` or :data:`_UND`.
        """
        if v in self.children[u]:
            return _OUT
        if u in self.children[v]:
            return _IN
        return _UND


#: Memoised adjacency tables, keyed on the MPDAG (which hashes on a canonical
#: key, so the cache changes cost and nothing else). Capped, because a
#: perturbation search walks through very many distinct graphs and an unbounded
#: table would be a slow memory leak.
_INDEX_CACHE: dict[MPDAG, _Index] = {}

#: Entries retained before :func:`_index` drops the whole table. A plain
#: high-water mark rather than an LRU: the searches that matter re-query the
#: same graph many times in a row, so recency is not worth tracking.
_INDEX_CACHE_LIMIT = 4096


def clear_index_cache() -> None:
    """Drop the memoised adjacency tables. Only for tests and memory management."""
    _INDEX_CACHE.clear()


def _index(g: MPDAG) -> _Index:
    """Return ``g``'s adjacency tables, building them on first use.

    Args:
        g: The MPDAG.

    Returns:
        The memoised :class:`_Index`.
    """
    idx = _INDEX_CACHE.get(g)
    if idx is None:
        if len(_INDEX_CACHE) >= _INDEX_CACHE_LIMIT:
            _INDEX_CACHE.clear()
        idx = _Index(g)
        _INDEX_CACHE[g] = idx
    return idx


def _as_node_set(value: Node | Iterable[Node]) -> set[Node]:
    """Normalise a single node or an iterable of nodes to a plain set.

    Args:
        value: A node label, or any iterable of node labels.

    Returns:
        A plain ``set`` of node labels.
    """
    if isinstance(value, str):
        return {value}
    return set(value)


def is_possibly_causal_step(g: MPDAG, a: Node, b: Node) -> bool:
    """Whether the edge between ``a`` and ``b`` may be traversed forward.

    A step ``a -> b`` or ``a - b`` is compatible with the path being causal from
    ``a``'s side; a step ``a <- b`` is not.

    Args:
        g: The MPDAG.
        a: The node the step leaves.
        b: The node the step enters.

    Returns:
        True iff the edge is ``a -> b`` or ``a - b``.
    """
    return g.is_directed_edge(a, b) or g.is_undirected_edge(a, b)


def is_possibly_causal(g: MPDAG, path: Path) -> bool:
    """Whether ``path`` is a possibly causal path (from ``path[0]``).

    A path ``<v0, ..., vk>`` is possibly causal iff no edge on it points back
    towards ``v0``, i.e. no ``vi <- vi+1`` occurs. Equivalently, every edge is
    either ``vi -> vi+1`` or ``vi - vi+1``. A path that is not possibly causal
    is called **non-causal**.

    Args:
        g: The MPDAG.
        path: A sequence of distinct, consecutively adjacent nodes. Adjacency is
            assumed, not re-checked.

    Returns:
        True iff no edge on ``path`` points backwards.
    """
    for i in range(len(path) - 1):
        if g.is_directed_edge(path[i + 1], path[i]):
            return False
    return True


def is_non_causal(g: MPDAG, path: Path) -> bool:
    """Whether ``path`` is non-causal, i.e. not possibly causal.

    Args:
        g: The MPDAG.
        path: A sequence of distinct, consecutively adjacent nodes.

    Returns:
        True iff some edge on ``path`` points back towards ``path[0]``.
    """
    return not is_possibly_causal(g, path)


def simple_paths(g: MPDAG, x: Node, y: Node) -> list[Path]:
    """Every simple path from ``x`` to ``y`` in the skeleton of ``g``.

    .. warning::
       **Exponential in the number of vertices, and not used by the decision.**
       A dense graph on ``n`` nodes has ``Theta((n-2)!)`` simple paths between
       any two of them. Condition (c) of the criterion is decided by
       :func:`open_definite_status_non_causal_path`, which is a polynomial
       state search; this function is kept as its slow reference and for
       reports that need to name the paths. Do not call it at scale.

    Orientation is ignored during the walk: what a path *is* does not depend on
    arrowheads, only on adjacency. Since ``x`` is a single node here, every such
    path is automatically *proper* (a path from a set ``X`` is proper when only
    its first node lies in ``X``), so the criterion's "proper" qualifier is
    vacuous throughout this package.

    Args:
        g: The MPDAG.
        x: The start node.
        y: The end node.

    Returns:
        The paths as node tuples, in a deterministic (depth-first, sorted-
        neighbour) order. Empty if ``x == y``.
    """
    if x == y:
        return []
    out: list[Path] = []
    current: list[Node] = [x]
    visited: set[Node] = {x}

    def walk() -> None:
        node = current[-1]
        if node == y:
            out.append(tuple(current))
            return
        for nxt in sorted(g.adjacent(node)):
            if nxt in visited:
                continue
            visited.add(nxt)
            current.append(nxt)
            walk()
            current.pop()
            visited.discard(nxt)

    walk()
    return out


def possibly_causal_paths(g: MPDAG, x: Node, y: Node) -> list[Path]:
    """Every simple possibly causal path from ``x`` to ``y``.

    .. warning::
       **Exponential in the number of vertices, and not used by the decision.**
       See :func:`simple_paths`. Kept as a reference and a reporting aid.

    Enumerated directly (the walk only ever takes ``->`` or ``-`` steps) rather
    than by filtering :func:`simple_paths`, so the search space stays small on
    graphs where most paths are non-causal.

    Args:
        g: The MPDAG.
        x: The start node.
        y: The end node.

    Returns:
        The possibly causal paths as node tuples, deterministically ordered.
    """
    if x == y:
        return []
    out: list[Path] = []
    current: list[Node] = [x]
    visited: set[Node] = {x}

    def walk() -> None:
        node = current[-1]
        if node == y:
            out.append(tuple(current))
            return
        for nxt in sorted(g.adjacent(node)):
            if nxt in visited or not is_possibly_causal_step(g, node, nxt):
                continue
            visited.add(nxt)
            current.append(nxt)
            walk()
            current.pop()
            visited.discard(nxt)

    walk()
    return out


def is_unshielded(g: MPDAG, path: Path) -> bool:
    """Whether ``path`` is unshielded: no two nodes two apart on it are adjacent.

    Args:
        g: The MPDAG.
        path: A sequence of distinct, consecutively adjacent nodes.

    Returns:
        True iff ``path[i-1]`` and ``path[i+1]`` are non-adjacent for every
        interior index ``i``.
    """
    for i in range(1, len(path) - 1):
        if g.has_edge(path[i - 1], path[i + 1]):
            return False
    return True


def unshielded_possibly_causal_paths(g: MPDAG, x: Node, y: Node) -> list[Path]:
    """Every simple possibly causal path from ``x`` to ``y`` that is also unshielded.

    .. warning::
       **Exponential in the number of vertices, and not used by the decision.**
       :func:`possible_descendants` -- the only consumer that mattered -- now
       goes through :func:`unshielded_reachable` instead. Kept as that
       function's slow reference.

    Enumerated with the shielding test applied as a pruning rule during the
    walk, so shielded prefixes are never extended.

    Args:
        g: The MPDAG.
        x: The start node.
        y: The end node.

    Returns:
        The unshielded possibly causal paths, deterministically ordered.
    """
    if x == y:
        return []
    out: list[Path] = []
    current: list[Node] = [x]
    visited: set[Node] = {x}

    def walk() -> None:
        node = current[-1]
        if node == y:
            out.append(tuple(current))
            return
        for nxt in sorted(g.adjacent(node)):
            if nxt in visited or not is_possibly_causal_step(g, node, nxt):
                continue
            if len(current) >= 2 and g.has_edge(current[-2], nxt):
                continue
            visited.add(nxt)
            current.append(nxt)
            walk()
            current.pop()
            visited.discard(nxt)

    walk()
    return out


def possible_descendants(g: MPDAG, x: Node | Iterable[Node]) -> set[Node]:
    r"""Nodes that are descendants of ``x`` in *some* DAG extension, ``x`` included.

    ``possde(x, G)``. Computed as the nodes reachable from ``x`` by an
    **unshielded** possibly causal path. The set is inclusive: ``x`` (or all of
    ``x``) is always in the result, matching the convention under which ``forb``
    is built.

    Why unshielded, and not just "possibly causal"
    ----------------------------------------------
    The textbook shortcut -- reach by any ``->``/``-`` walk -- is **wrong for
    MPDAGs carrying background knowledge**, and getting this wrong was the only
    substantive bug found while validating this package. Minimal witness, an
    element of this repository's own perturbation space::

        V0 -> V1,  V0 - V2,  V1 - V2      (a triangle with one edge known)

    ``<V1, V2, V0>`` is a possibly causal path (no arrowhead points back at
    ``V1``), so the naive rule puts ``V0`` in ``possde(V1)``. But all three DAG
    extensions of this graph keep ``V0 -> V1``, so ``V0`` is a descendant of
    ``V1`` in *none* of them: orienting the path forward would close the cycle
    ``V0 -> V1 -> V2 -> V0``. The path is shielded (``V1`` and ``V0`` are
    adjacent), and dropping shielded paths repairs it exactly.

    Correctness, both directions:

    * ``possde \subseteq`` unshielded-reachable: take a directed path
      ``x -> ... -> v`` in an extension ``D``. Whenever two of its nodes two
      apart are adjacent, that edge must run forwards in ``D`` (backwards would
      close a cycle), so it can be short-circuited; iterating leaves an
      unshielded directed path in ``D``, which is an unshielded possibly causal
      path in ``G`` (same skeleton).
    * the converse is the standard MPDAG lemma that an unshielded possibly
      causal path is causal in some extension.

    Verified equal to ``union over D in [G] of de(x, D)`` on every graph of the
    n=3 and n=4 sweep scope (1,588 graphs) and on a 25,268-graph n=5 sample,
    with zero mismatches; the naive walk rule mismatched on 570 and 7,401 of
    those respectively.

    Cost: polynomial. Decided by :func:`unshielded_reachable`, one edge-state
    search per element of ``x``, i.e. ``O(V*E)`` per source. Results for single
    nodes are memoised on the graph's index.

    Args:
        g: The MPDAG.
        x: A node, or a set of nodes, to start from.

    Returns:
        The inclusive set of possible descendants.
    """
    idx = _index(g)
    out: set[Node] = set()
    for src in sorted(_as_node_set(x)):
        cached = idx.possde.get(src)
        if cached is None:
            cached = unshielded_reachable(g, src)
            idx.possde[src] = cached
        out |= cached
    return out


def unshielded_reachable(g: MPDAG, src: Node) -> set[Node]:
    """Nodes reachable from ``src`` by an unshielded possibly causal path, inclusive.

    The polynomial replacement for ``_unshielded_reach_by_paths``, and the
    engine behind :func:`possible_descendants`.

    How the exponent goes away
    --------------------------
    The path formulation asks for a *sequence*; the answer only needs the set of
    endpoints. The three constraints on an extension of a walk -- the step must
    be ``->`` or ``-``, the step must not return to the previous node, and the
    previous node must not be adjacent to the next -- all refer to at most the
    last **two** nodes. So the walk's whole relevant history is the traversed
    edge ``(u, v)``, and a search over those ``O(E)`` states, each expanded once
    over its ``O(deg)`` continuations, visits every reachable endpoint in
    ``O(V*E)`` total.

    Walks versus paths
    ------------------
    Dropping the sequence drops the record of which nodes were already visited,
    so the search ranges over unshielded possibly causal *walks*, a superset of
    the paths. On a Meek-closed MPDAG the two have the same endpoint set, and
    the reason is worth stating because it is what licenses the whole move:

    * On such a graph an unshielded step ``u -> v`` can never be followed by
      ``v - w``. That triple is unshielded by construction, so Meek's rule 1
      would have fired and oriented ``v -> w`` already. Hence along any
      unshielded possibly causal walk the undirected steps all precede the
      directed ones.
    * If a walk repeats a node, take a repeat with the shortest gap; the segment
      between the two occurrences is a closed walk, all of whose steps are
      ``->`` or ``-`` forwards, whose nodes are distinct, and which is
      unshielded including at the wrap-around. By the previous point it cannot
      mix step kinds -- going round the loop, a directed step would eventually
      have to be followed by an undirected one. All-directed is a directed cycle,
      which an MPDAG does not have. All-undirected is a cycle of length >= 3
      (length 3 is impossible: an unshielded triple has non-adjacent ends) with
      no chord at distance 2, i.e. a chordless cycle in a chain component, which
      an MPDAG does not have either.

    That argument is an argument, not a proof-checked theorem, so it is backed
    by measurement: this function is asserted equal to the path enumeration on
    every graph in the exhaustive ``n <= 4`` scope and on sampled ``n = 5``
    graphs (``tests/criterion/test_criterion.py``).

    Args:
        g: The MPDAG.
        src: The node to start from.

    Returns:
        ``{src}`` together with every node such a path ends at.
    """
    idx = _index(g)
    if src not in idx.adj:
        return {src}
    out: set[Node] = {src}
    seen: set[tuple[Node, Node]] = set()
    stack: list[tuple[Node, Node]] = []
    children, neighbors, adj, adj_set = idx.children, idx.neighbors, idx.adj, idx.adj_set

    for nxt in adj[src]:
        if nxt in children[src] or nxt in neighbors[src]:
            out.add(nxt)
            seen.add((src, nxt))
            stack.append((src, nxt))

    while stack:
        prev, node = stack.pop()
        prev_adj = adj_set[prev]
        for nxt in adj[node]:
            if nxt == prev or nxt in prev_adj:
                continue
            if nxt not in children[node] and nxt not in neighbors[node]:
                continue
            state = (node, nxt)
            if state in seen:
                continue
            seen.add(state)
            out.add(nxt)
            stack.append(state)
    return out


def _unshielded_reach_by_paths(g: MPDAG, src: Node) -> set[Node]:
    """Endpoints of unshielded possibly causal paths from ``src``, by enumeration.

    .. warning::
       **Exponential in the number of vertices, and not used by the decision.**
       The slow reference :func:`unshielded_reachable` is tested against; kept
       for exactly that purpose.

    Args:
        g: The MPDAG.
        src: The node to start from.

    Returns:
        ``{src}`` together with every node such a path ends at.
    """
    out: set[Node] = {src}
    current: list[Node] = [src]
    visited: set[Node] = {src}

    def walk() -> None:
        node = current[-1]
        for nxt in sorted(g.adjacent(node)):
            if nxt in visited or not is_possibly_causal_step(g, node, nxt):
                continue
            if len(current) >= 2 and g.has_edge(current[-2], nxt):
                continue
            out.add(nxt)
            visited.add(nxt)
            current.append(nxt)
            walk()
            current.pop()
            visited.discard(nxt)

    walk()
    return out


def is_collider(g: MPDAG, path: Path, i: int) -> bool:
    """Whether ``path[i]`` is a collider on ``path``.

    A collider is an interior node entered by an arrowhead from both sides:
    ``path[i-1] -> path[i] <- path[i+1]``. Endpoints are never colliders.

    Args:
        g: The MPDAG.
        path: A sequence of distinct, consecutively adjacent nodes.
        i: An index into ``path``.

    Returns:
        True iff ``0 < i < len(path) - 1`` and both incident path edges are
        directed into ``path[i]``.
    """
    if i <= 0 or i >= len(path) - 1:
        return False
    return g.is_directed_edge(path[i - 1], path[i]) and g.is_directed_edge(path[i + 1], path[i])


def is_definite_non_collider(g: MPDAG, path: Path, i: int) -> bool:
    """Whether ``path[i]`` is a *definite* non-collider on ``path``.

    Two ways an interior node ``v = path[i]`` can be certified a non-collider
    without knowing which DAG extension is the truth:

    1. **An edge out of ``v`` on the path.** If ``v -> path[i-1]`` or
       ``v -> path[i+1]``, then in every extension at least one incident path
       edge has its tail at ``v``, so ``v`` cannot be a collider.
    2. **The unshielded undirected pattern.** If ``path[i-1] - v - path[i+1]``
       (both incident path edges undirected) *and* ``path[i-1]`` and
       ``path[i+1]`` are non-adjacent, then orienting both edges into ``v``
       would create a new unshielded collider, which no consistent extension of
       an MPDAG may introduce. So ``v`` is a non-collider in every extension.

    Everything else is of **indefinite** status. The mixed pattern
    ``path[i-1] -> v - path[i+1]`` with ``path[i-1]``, ``path[i+1]``
    non-adjacent does not arise in a Meek-closed graph (Meek's rule 1 would
    already have oriented ``v -> path[i+1]``, landing in case 1); when it does
    arise -- because a caller passed a graph that is not Meek-closed -- this
    function correctly reports "not definite", which is the conservative side.

    Args:
        g: The MPDAG.
        path: A sequence of distinct, consecutively adjacent nodes.
        i: An index into ``path``.

    Returns:
        True iff ``0 < i < len(path) - 1`` and one of the two patterns holds.
    """
    if i <= 0 or i >= len(path) - 1:
        return False
    left, mid, right = path[i - 1], path[i], path[i + 1]
    if g.is_directed_edge(mid, left) or g.is_directed_edge(mid, right):
        return True
    return (
        g.is_undirected_edge(left, mid)
        and g.is_undirected_edge(mid, right)
        and not g.has_edge(left, right)
    )


def is_definite_status_node(g: MPDAG, path: Path, i: int) -> bool:
    """Whether ``path[i]`` is of definite status on ``path``.

    Endpoints are of definite status by convention (they have no role to be
    ambiguous about). Interior nodes qualify iff they are a collider or a
    definite non-collider.

    Args:
        g: The MPDAG.
        path: A sequence of distinct, consecutively adjacent nodes.
        i: An index into ``path``.

    Returns:
        True iff ``path[i]`` is an endpoint, a collider, or a definite
        non-collider.
    """
    if i <= 0 or i >= len(path) - 1:
        return True
    return is_collider(g, path, i) or is_definite_non_collider(g, path, i)


def is_definite_status_path(g: MPDAG, path: Path) -> bool:
    """Whether every node on ``path`` is of definite status.

    Args:
        g: The MPDAG.
        path: A sequence of distinct, consecutively adjacent nodes.

    Returns:
        True iff no node on ``path`` has indefinite collider status.
    """
    return all(is_definite_status_node(g, path, i) for i in range(len(path)))


def is_blocked(g: MPDAG, path: Path, z: Iterable[Node]) -> bool:
    """Whether the definite-status ``path`` is blocked by ``z`` in every extension.

    ``path`` is assumed to be of definite status (see
    :func:`is_definite_status_path`); on such a path every interior node is a
    collider or a definite non-collider *in every DAG extension alike*, so
    "blocked" is well defined at the MPDAG level. The path is blocked iff some
    interior node blocks it:

    * a definite non-collider that lies in ``z``; or
    * a collider ``b`` with ``possde(b, G) & z`` empty.

    The collider clause uses **possible** descendants, not the descendants
    reachable by directed edges alone. If ``b - d`` is undirected and ``d`` is
    in ``z``, then the extension that orients ``b -> d`` makes ``d`` a genuine
    descendant of the collider and re-opens the path; a predicate that must hold
    in *every* extension therefore may not count such a ``b`` as blocking. Note
    ``possde`` is inclusive, so ``b`` itself being in ``z`` already opens the
    collider, as it should.

    Args:
        g: The MPDAG.
        path: A definite-status path.
        z: The conditioning set.

    Returns:
        True iff ``path`` is blocked by ``z`` in every DAG extension of ``g``.
    """
    z_set = set(z)
    for i in range(1, len(path) - 1):
        node = path[i]
        if is_collider(g, path, i):
            if not (possible_descendants(g, node) & z_set):
                return True
        elif node in z_set:
            return True
    return False


# --------------------------------------------------------------------------------
# Condition (c) as a state search
# --------------------------------------------------------------------------------


def open_definite_status_non_causal_path(
    g: MPDAG, x: Node, y: Node, z: Iterable[Node]
) -> Path | None:
    r"""A proper definite-status non-causal path from ``x`` to ``y`` left open by ``z``.

    This is condition (c) of the adjustment criterion, decided in polynomial
    time. It is what :func:`bkrobust.mpdag_criterion.criterion.why_invalid`
    calls; the enumerating
    :func:`~bkrobust.mpdag_criterion.criterion.open_non_causal_path` is kept
    only as its reference.

    Why a state search decides a "for every path" question
    -----------------------------------------------------
    Condition (c) reads as a quantifier over paths, but every clause in it is
    **local**: whether ``v`` is a collider, whether it is a *definite*
    non-collider, and whether it blocks, are all determined by the triple
    ``(u, v, w)`` of the previous, current and next node -- plus, for the
    unshielded-undirected pattern, whether ``u`` and ``w`` are adjacent. The one
    non-local clause, "the path is non-causal", is a single bit: has any step so
    far pointed backwards? So a walk's entire relevant history is
    ``(previous node, current node, seen-a-backward-step)``, and the question
    "does an open one reach ``y``?" is reachability over those states. This is
    the same reduction Bayes-ball makes for d-separation, and the same one
    :func:`bkrobust.demo.evaluate.is_dseparated` relies on.

    Concretely, writing a step ``u`` to ``v`` as ``out`` (``u -> v``), ``in``
    (``v -> u``) or ``und`` (``u - v``), the interior node ``v`` between an
    incoming step of kind ``i`` and an outgoing step of kind ``o`` is:

    * a **collider** iff ``i == out`` and ``o == in``; it leaves the path open
      only when ``possde(v, G)`` meets ``z``;
    * a **definite non-collider** iff ``i == in`` (the edge ``v -> u`` is an
      arrow out of ``v``), or ``o == out``, or both steps are ``und`` and ``u``
      is non-adjacent to ``w``; it leaves the path open only when ``v`` is not
      in ``z``;
    * of **indefinite status** otherwise, and the walk is dropped -- the two
      dropped combinations are ``out``/``und`` and ``und``/``in``.

    The path is non-causal iff at least one step is ``in``, which is the bit
    carried in the state, and is required before ``y`` is accepted.

    Search over walks, answer about paths
    -------------------------------------
    As in :func:`unshielded_reachable`, states forget which nodes have already
    been used, so the search ranges over *walks*. Two constraints that a simple
    path satisfies for free are re-imposed explicitly, because they are cheap
    and they tighten the relaxation: a step never returns to the node it just
    came from, and no walk re-enters ``x`` or continues through ``y``. The
    residual gap is measured rather than assumed -- the differential test in
    ``tests/criterion/test_criterion.py`` asserts this function agrees with the
    enumerating decision on every ``(graph, x, y, z)`` in the exhaustive
    ``n <= 4`` scope, and that the witness it returns is always a genuine simple
    path satisfying :func:`is_definite_status_path`, :func:`is_non_causal` and
    ``not`` :func:`is_blocked`.

    Cost: ``O(E)`` states times ``O(deg)`` continuations, so ``O(V*E)``, plus at
    most one ``O(V*E)`` :func:`possible_descendants` call per distinct collider
    node (memoised on the graph). Breadth-first, so the witness returned is a
    shortest one, and neighbours are visited in sorted order, so it is
    deterministic.

    Args:
        g: The MPDAG, assumed Meek-closed.
        x: The treatment node.
        y: The outcome node.
        z: The conditioning set.

    Returns:
        A witness path, or ``None`` if every proper definite-status non-causal
        path from ``x`` to ``y`` is blocked by ``z``.
    """
    idx = _index(g)
    if x == y or x not in idx.adj or y not in idx.adj:
        return None
    z_set = set(z)
    adj, adj_set, children = idx.adj, idx.adj_set, idx.children
    kind = idx.step_kind

    def collider_is_open(node: Node) -> bool:
        """Whether a collider at ``node`` is opened by ``z`` (memoised via the index)."""
        cached = idx.possde.get(node)
        if cached is None:
            cached = unshielded_reachable(g, node)
            idx.possde[node] = cached
        return bool(cached & z_set)

    # State: (previous node, current node, has a backward step been taken?).
    # ``parent`` doubles as the visited set and as the back-pointers the witness
    # is reconstructed from.
    parent: dict[tuple[Node, Node, bool], tuple[Node, Node, bool] | None] = {}
    queue: deque = deque()
    for nxt in adj[x]:
        backward = x in children[nxt]
        state = (x, nxt, backward)
        parent[state] = None
        if nxt == y:
            if backward:
                return _rebuild(parent, state)
            continue
        queue.append(state)

    while queue:
        state = queue.popleft()
        prev, node, backward = state
        in_kind = kind(prev, node)
        prev_adj = adj_set[prev]
        node_in_z = node in z_set
        for nxt in adj[node]:
            if nxt == prev or nxt == x:
                continue
            if nxt in children[node]:
                out_kind = _OUT
            elif node in children[nxt]:
                out_kind = _IN
            else:
                out_kind = _UND
            if in_kind == _OUT and out_kind == _IN:
                if not collider_is_open(node):
                    continue
            else:
                definite = (
                    in_kind == _IN
                    or out_kind == _OUT
                    or (in_kind == _UND and out_kind == _UND and nxt not in prev_adj)
                )
                if not definite or node_in_z:
                    continue
            nxt_backward = backward or out_kind == _IN
            nxt_state = (node, nxt, nxt_backward)
            if nxt_state in parent:
                continue
            parent[nxt_state] = state
            if nxt == y:
                if nxt_backward:
                    return _rebuild(parent, nxt_state)
                continue
            queue.append(nxt_state)
    return None


def _rebuild(
    parent: dict[tuple[Node, Node, bool], tuple[Node, Node, bool] | None],
    state: tuple[Node, Node, bool],
) -> Path:
    """Turn a chain of search states back into the node sequence it traversed.

    Args:
        parent: Back-pointers written by
            :func:`open_definite_status_non_causal_path`.
        state: The accepting state.

    Returns:
        The nodes of the witness, from the start node to the end node.
    """
    nodes: list[Node] = [state[1]]
    cur: tuple[Node, Node, bool] | None = state
    while cur is not None:
        nodes.append(cur[0])
        cur = parent[cur]
    nodes.reverse()
    return tuple(nodes)


def has_open_definite_status_non_causal_path(g: MPDAG, x: Node, y: Node, z: Iterable[Node]) -> bool:
    """Whether condition (c) of the criterion *fails* for ``z``.

    Args:
        g: The MPDAG, assumed Meek-closed.
        x: The treatment node.
        y: The outcome node.
        z: The conditioning set.

    Returns:
        True iff some proper definite-status non-causal path from ``x`` to ``y``
        is left open by ``z``.
    """
    return open_definite_status_non_causal_path(g, x, y, z) is not None
