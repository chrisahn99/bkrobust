"""Meek's rules, MPDAG validity, and DAG-extension enumeration.

Everything here is a pure function over :class:`~bkrobust.demo.graph.MPDAG` --
no mutation, no hidden state. The four functions ``meek_rule_1``..``meek_rule_4``
each compute the orientations a single Meek rule forces in one pass; the rest of
the module builds on them: closing a partial graph to a fixpoint
(:func:`meek_closure`), imposing a batch of background-knowledge orientations
(:func:`apply_orientations`), checking Meek-closedness and chordality, and
enumerating -- or validating -- the DAGs consistent with a given MPDAG.

Meek's rules are confluent: applying them in any legal order, one at a time or
in batches, drives a PDAG to the same fixpoint (or the same FAIL). The internal
:func:`_meek_closure_ordered` exposes the choice of order as a parameter so that
property is directly testable rather than merely assumed.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from itertools import combinations, product

from bkrobust.demo.graph import MPDAG, Edge, Node, undirected_components

# --- individual rules --------------------------------------------------------


def meek_rule_1(g: MPDAG) -> set[Edge]:
    """Orientations forced by Meek's rule 1 in a single pass.

    Pattern: ``a -> b``, ``b - c``, and ``a``, ``c`` not adjacent. Leaving
    ``b - c`` unoriented as ``c -> b`` would create a new, unshielded
    v-structure ``a -> b <- c`` that the background knowledge did not license,
    so ``b - c`` must be ``b -> c``.

    Args:
        g: The current (possibly partial) graph. Not mutated.

    Returns:
        The set of ``(tail, head)`` orientations this rule forces, each for an
        edge that is currently undirected in ``g``.
    """
    forced: set[Edge] = set()
    for a, b in g.directed_edges:
        for c in g.neighbors(b):
            if not g.has_edge(a, c):
                forced.add((b, c))
    return forced


def meek_rule_2(g: MPDAG) -> set[Edge]:
    """Orientations forced by Meek's rule 2 in a single pass.

    Pattern: ``a -> c -> b`` and ``a - b``. Leaving ``a - b`` as ``b -> a``
    would close the directed path into a cycle, so ``a - b`` must be
    ``a -> b``.

    Args:
        g: The current (possibly partial) graph. Not mutated.

    Returns:
        The set of ``(tail, head)`` orientations this rule forces, each for an
        edge that is currently undirected in ``g``.
    """
    forced: set[Edge] = set()
    for u, v in g.undirected_edges:
        if g.children(u) & g.parents(v):
            forced.add((u, v))
        if g.children(v) & g.parents(u):
            forced.add((v, u))
    return forced


def meek_rule_3(g: MPDAG) -> set[Edge]:
    """Orientations forced by Meek's rule 3 in a single pass.

    Pattern: ``a - b``, ``a - c``, ``a - d``, ``c -> b``, ``d -> b``, and
    ``c``, ``d`` not adjacent. Leaving ``a - b`` as ``b -> a`` would force
    ``a - c`` and ``a - d`` to also orient into ``a`` (by rule 1, since
    ``c``/``d`` and ``b`` are non-adjacent-safe already-directed neighbours),
    ultimately compelling a new unshielded v-structure at ``a``; ``a - b``
    must be ``a -> b``.

    Args:
        g: The current (possibly partial) graph. Not mutated.

    Returns:
        The set of ``(tail, head)`` orientations this rule forces, each for an
        edge that is currently undirected in ``g``.
    """
    forced: set[Edge] = set()
    for a in g.nodes:
        nbrs = g.neighbors(a)
        for b in nbrs:
            candidates = sorted(nbrs & g.parents(b))
            for c, d in combinations(candidates, 2):
                if not g.has_edge(c, d):
                    forced.add((a, b))
                    break
    return forced


def meek_rule_4(g: MPDAG) -> set[Edge]:
    """Orientations forced by Meek's rule 4 in a single pass.

    Pattern: ``a - b``, ``a - c``, ``c -> d``, ``d -> b``, and ``c``, ``b``
    not adjacent. Leaving ``a - b`` as ``b -> a`` would, together with
    ``a - c`` and the directed path ``c -> d -> b``, compel a new unshielded
    v-structure; ``a - b`` must be ``a -> b``.

    Args:
        g: The current (possibly partial) graph. Not mutated.

    Returns:
        The set of ``(tail, head)`` orientations this rule forces, each for an
        edge that is currently undirected in ``g``.
    """
    forced: set[Edge] = set()
    for a in g.nodes:
        nbrs = g.neighbors(a)
        for b in nbrs:
            for c in nbrs:
                if c == b or g.has_edge(c, b):
                    continue
                if g.children(c) & g.parents(b):
                    forced.add((a, b))
    return forced


def _all_forced(g: MPDAG) -> set[Edge]:
    """Union of everything R1-R4 force in one pass over ``g``."""
    return meek_rule_1(g) | meek_rule_2(g) | meek_rule_3(g) | meek_rule_4(g)


# --- closure -------------------------------------------------------------


def _meek_closure_ordered(g: MPDAG, key: Callable[[Edge], object]) -> MPDAG | None:
    """Close ``g`` under Meek's rules, applying one forced edge at a time.

    At each round every rule is re-evaluated against the current graph, the
    full set of currently-forced orientations is collected, and -- if none
    conflict -- exactly one of them (the one sorted first by ``key``) is
    applied before the next round. Recomputing from scratch each round, rather
    than applying a whole round's orientations at once, keeps the procedure
    obviously correct at the cost of some redundant work, which is immaterial
    at the sizes this demonstration targets.

    Because each round strictly converts one undirected edge to directed, the
    loop terminates in at most ``len(g.undirected_edges)`` rounds.

    This is the order-parameterised engine that lets order-independence (a
    consequence of Meek's rules being confluent) be tested directly: different
    ``key`` functions impose different processing orders, and every legal
    order must reach the same fixpoint or the same FAIL.

    Args:
        g: The starting graph. Not mutated.
        key: A total order over ``Edge`` used to pick, among the orientations
            currently forced, which one to apply next.

    Returns:
        The closed graph, or ``None`` on FAIL (the rules force a directed
        cycle, or force the same edge in both directions -- including the
        case where ``g`` itself is already cyclic).
    """
    if not g.is_acyclic():
        return None

    current = g
    for _ in range(len(g.undirected_edges) + 1):
        forced = _all_forced(current)
        if not forced:
            return current
        if any((b, a) in forced for a, b in forced):
            return None
        tail, head = sorted(forced, key=key)[0]
        if current.is_directed_edge(head, tail):
            return None
        try:
            nxt = current.oriented(tail, head)
        except KeyError:
            return None
        if not nxt.is_acyclic():
            return None
        current = nxt
    # Unreachable in practice: each round consumes one undirected edge, so the
    # loop above always returns within its bound. Kept as a defensive FAIL
    # rather than silently returning a non-fixpoint graph.
    return None


def meek_closure(g: MPDAG) -> MPDAG | None:
    """Apply Meek's rules R1-R4 to a fixpoint.

    Args:
        g: The starting (possibly partial) graph. Not mutated.

    Returns:
        The Meek-closed graph, or ``None`` on FAIL: the rules force a directed
        cycle, or force the same edge in both directions (this includes the
        case where ``g`` already contains a directed cycle before any rule is
        applied).
    """
    return _meek_closure_ordered(g, key=lambda e: e)


def apply_orientations(g: MPDAG, orientations: Iterable[Edge]) -> MPDAG | None:
    """Impose a batch of background-knowledge orientations, closing after each.

    Orientations are processed one at a time in a deterministic (sorted)
    order, independent of any ordering in ``orientations`` itself, so the
    verdict never depends on iteration order of the input.

    Args:
        g: The starting (possibly partial) graph. Not mutated.
        orientations: ``(tail, head)`` pairs to impose, in addition to
            whatever ``g`` already entails.

    Returns:
        The resulting Meek-closed graph, or ``None`` on FAIL: an orientation
        names a pair that is not adjacent in the current graph, an
        orientation contradicts an already-directed edge, or the subsequent
        closure FAILs.
    """
    current = g
    for tail, head in sorted(set(orientations)):
        if not current.has_edge(tail, head):
            return None
        if current.is_directed_edge(head, tail):
            return None
        if not current.is_directed_edge(tail, head):
            current = current.oriented(tail, head)
        closed = meek_closure(current)
        if closed is None:
            return None
        current = closed
    return current


def is_meek_closed(g: MPDAG) -> bool:
    """Whether no Meek rule fires against ``g``.

    Args:
        g: The graph to check. Not mutated.

    Returns:
        True iff R1-R4 collectively force nothing new.
    """
    return not _all_forced(g)


# --- chordality --------------------------------------------------------


def _is_chordal_undirected(nodes: list[Node], edges: set[Edge]) -> bool:
    """Whether the undirected graph ``(nodes, edges)`` is chordal.

    Implemented via Maximum Cardinality Search (Tarjan & Yannakakis 1984):
    visit vertices one at a time, numbering the first visit ``1`` and the
    last ``n``, always choosing next an unvisited vertex adjacent to the most
    already-numbered vertices (ties broken by node label, for determinism).
    The resulting numbering is a perfect elimination ordering, read in
    *increasing* number, iff the graph is chordal, which is checked directly:
    for every vertex, its neighbours numbered before it must form a clique.

    Args:
        nodes: The vertex set.
        edges: Undirected edges as unordered pairs (either orientation).

    Returns:
        True iff the graph is chordal (vacuously true for 0-3 vertices, since
        no chordless cycle of length >= 4 can exist).
    """
    n = len(nodes)
    if n <= 3:
        return True

    adj = {v: set() for v in nodes}
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)

    weight = {v: 0 for v in nodes}
    numbered = {}
    remaining = set(nodes)
    for i in range(1, n + 1):
        v = max(remaining, key=lambda x: (weight[x], x))
        numbered[v] = i
        remaining.discard(v)
        for w in adj[v]:
            if w in remaining:
                weight[w] += 1

    for v in nodes:
        earlier = [w for w in adj[v] if numbered[w] < numbered[v]]
        for x, y in combinations(earlier, 2):
            if y not in adj[x]:
                return False
    return True


def is_chordal_components(g: MPDAG) -> bool:
    """Whether every undirected connected component of ``g`` is chordal.

    Directed edges play no role here: chordality is a property of the
    undirected skeleton within each undirected-connected component.

    Args:
        g: The graph to check. Not mutated.

    Returns:
        True iff every undirected component induces a chordal undirected
        graph (vacuously true if there are no undirected edges at all).
    """
    for comp in undirected_components(g):
        comp_nodes = sorted(comp)
        comp_edges = {e for e in g.undirected_edges if e[0] in comp and e[1] in comp}
        if not _is_chordal_undirected(comp_nodes, comp_edges):
            return False
    return True


# --- DAG extensions ------------------------------------------------------


def is_consistent_extension(dag: MPDAG, g: MPDAG) -> bool:
    """Whether ``dag`` is a consistent DAG extension of ``g``.

    ``dag`` qualifies iff it: shares ``g``'s node set and skeleton (every
    edge of ``g`` appears in ``dag``, oriented one way or the other, and no
    edge is added or removed); is acyclic; keeps every already-directed edge
    of ``g``; and introduces no new v-structure, i.e. every collider
    ``a -> b <- c`` with ``a``, ``c`` non-adjacent in ``dag`` was already a
    v-structure in ``g`` (both ``a -> b`` and ``c -> b`` already directed
    there).

    Args:
        dag: A candidate fully oriented graph.
        g: The (possibly partial) graph it should extend.

    Returns:
        True iff ``dag`` is a consistent extension of ``g`` per the above.
    """
    if dag.nodes != g.nodes:
        return False
    if dag.skeleton() != g.skeleton():
        return False
    if not dag.is_dag():
        return False
    if not (g.directed_edges <= dag.directed_edges):
        return False

    for b in dag.nodes:
        parents = sorted(dag.parents(b))
        for a, c in combinations(parents, 2):
            if not dag.has_edge(a, c):
                if not (g.is_directed_edge(a, b) and g.is_directed_edge(c, b)):
                    return False
    return True


def enumerate_dag_extensions(g: MPDAG, limit: int | None = None) -> list[MPDAG]:
    """Every DAG consistent with ``g`` per :func:`is_consistent_extension`.

    Brute forces all ``2**k`` orientations of ``g``'s ``k`` undirected edges
    (``g``'s directed edges are always kept as is), keeping those that are
    acyclic and introduce no new v-structure. Brute force is deliberate: it is
    the clearest possible implementation and this demonstration only ever
    targets graphs small enough (``k`` up to about a dozen undirected edges)
    for it to be cheap.

    Args:
        g: The (possibly partial) graph to extend. Not mutated.
        limit: If given, stop once this many valid extensions have been
            found (a safety valve for larger graphs, not a request for the
            lexicographically-first extensions). The returned list is always
            sorted, regardless of whether ``limit`` was hit.

    Returns:
        The consistent extensions, sorted deterministically by
        ``MPDAG.edge_string()``.
    """
    undirected = sorted(g.undirected_edges)
    k = len(undirected)
    results: list[MPDAG] = []

    for bits in product((False, True), repeat=k):
        # zip()'s strict= kwarg needs Python 3.10+; this module targets 3.9, and
        # the lengths are equal by construction (`bits` is drawn from
        # `product(..., repeat=k)` with `k == len(undirected)`).
        assigned = [
            (b, a) if flip else (a, b)
            for (a, b), flip in zip(undirected, bits)  # noqa: B905
        ]
        directed_all = set(g.directed_edges) | set(assigned)
        try:
            candidate = MPDAG(g.nodes, directed_all, [])
        except ValueError:
            continue
        if not is_consistent_extension(candidate, g):
            continue
        results.append(candidate)
        if limit is not None and len(results) >= limit:
            break

    results.sort(key=lambda d: d.edge_string())
    return results


def is_valid_mpdag(g: MPDAG) -> bool:
    """Whether ``g`` is a valid MPDAG.

    A graph qualifies iff: its directed part is acyclic; it is Meek-closed
    (nothing left for R1-R4 to force); every undirected connected component
    is chordal; and it admits at least one consistent DAG extension.

    Args:
        g: The graph to check. Not mutated.

    Returns:
        True iff all four conditions hold.
    """
    return (
        g.is_acyclic()
        and is_meek_closed(g)
        and is_chordal_components(g)
        and len(enumerate_dag_extensions(g)) > 0
    )
