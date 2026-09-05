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

Determinism: every iteration over a graph accessor (which return plain ``set``s)
goes through ``sorted``, so path enumeration order never depends on
``PYTHONHASHSEED``.

Written to run on Python 3.9 as well as the repository's 3.11 target, matching
the rest of the package: builtin generics and ``X | Y`` unions appear only in
annotations, and ``zip`` is never called with ``strict=`` (3.10+).
"""

from __future__ import annotations

from collections.abc import Iterable

from bkrobust.demo.graph import MPDAG

Node = str
Path = tuple[Node, ...]


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

    Args:
        g: The MPDAG.
        x: A node, or a set of nodes, to start from.

    Returns:
        The inclusive set of possible descendants.
    """
    out: set[Node] = set()
    for src in sorted(_as_node_set(x)):
        out |= _unshielded_reach(g, src)
    return out


def _unshielded_reach(g: MPDAG, src: Node) -> set[Node]:
    """Nodes reachable from ``src`` by an unshielded possibly causal path, inclusive.

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
