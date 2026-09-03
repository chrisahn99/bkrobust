r"""The perturbation space of MPDAGs refining a CPDAG.

Given a (Meek-closed) CPDAG :math:`\hat C`, this module brute-force
constructs the space :math:`\mathfrak G_{\hat C}` of every valid MPDAG
reachable by orienting some subset of :math:`\hat C`'s undirected edges,
the model-inclusion partial order over that space (:math:`G_1 \preceq G_2`
iff :math:`[G_1] \subseteq [G_2]`, where :math:`[G]` is the set of DAGs
:func:`~bkrobust.demo.meek.enumerate_dag_extensions` finds consistent with
``G``), its covering relation, and the metric induced by BFS distance on the
resulting neighbour graph.

Everything here is brute force by design: the covering relation is computed
by enumerating and comparing represented-DAG sets directly, never via a
conjectured shortcut (such as "covering pairs differ by exactly one
orientation", which is false in general once Meek propagation is taken into
account -- see :func:`covering_pairs`).

A note on :func:`enumerate_space`'s filter (read this before changing it):
``meek.is_valid_mpdag(G)`` alone checks only that ``G`` is *internally*
self-consistent -- it says nothing about whether ``G`` is actually a
legitimate refinement of :math:`\hat C`. Concretely, on the unshielded-triple
skeleton ``A-B-C`` (no ``A-C`` edge), the fully oriented ``A->B<-C`` is
acyclic, trivially Meek-closed (there are no undirected edges left for R1-4
to act on), and trivially admits itself as its own sole "consistent
extension" -- so ``is_valid_mpdag(A->B<-C)`` is ``True`` in isolation. But if
:math:`\hat C` is the fully undirected chain ``A-B-C``, that collider belongs
to a *different* Markov equivalence class: orienting both edges into ``B`` is
exactly the new v-structure Meek's rule 1 exists to forbid, and
``is_consistent_extension(A->B<-C, \hat C)`` is correctly ``False``. Filtering
:func:`enumerate_space` by ``is_valid_mpdag`` alone therefore admits this
foreign-class graph and breaks the "``cpdag`` is the unique maximum"
invariant (its extension set is *not* a subset of :math:`\hat C`'s). This is
verified directly in ``tests/demo/test_space.py`` via the same construction.
The fix, also brute force, is the second half of the filter: every one of
``G``'s own consistent extensions must *also* be a consistent extension of
:math:`\hat C`, i.e. genuine model inclusion :math:`[G] \subseteq [\hat C]`,
checked by direct enumeration and comparison, not assumed.
"""

from __future__ import annotations

from collections import deque
from itertools import product

from bkrobust.demo.graph import MPDAG, Edge, Node
from bkrobust.demo.meek import (
    enumerate_dag_extensions,
    is_consistent_extension,
    is_valid_mpdag,
)

# --- space enumeration -------------------------------------------------------


def enumerate_space(cpdag: MPDAG) -> list[MPDAG]:
    """Enumerate every valid MPDAG refinement of ``cpdag``.

    Brute forces all ``3**k`` assignments of {leave undirected, orient one
    way, orient the other} to ``cpdag``'s ``k`` undirected edges (``cpdag``'s
    own directed edges are always kept fixed), keeping a candidate ``G`` iff:

    1. ``meek.is_valid_mpdag(G)``: ``G`` is acyclic, Meek-closed, chordal in
       every undirected component, and admits at least one of its own
       consistent DAG extensions.
    2. Every one of ``G``'s own consistent DAG extensions is *also* a
       consistent extension of ``cpdag`` -- genuine model inclusion
       ``[G] subseteq [cpdag]``, checked by direct enumeration. See the
       module docstring for why condition 1 alone is insufficient (it can
       admit graphs from a different Markov equivalence class than
       ``cpdag``'s).

    Args:
        cpdag: The background-knowledge MPDAG the space is anchored at. Must
            itself be a valid, Meek-closed MPDAG (not mutated).

    Returns:
        The space, deduplicated (each of the ``3**k`` raw assignments
        determines a distinct graph by construction, so no collisions can
        occur) and sorted by ``(len(directed_edges), edge_string())``.
    """
    undirected = sorted(cpdag.undirected_edges)
    k = len(undirected)
    base_directed = set(cpdag.directed_edges)

    space: list[MPDAG] = []
    for state in product((0, 1, 2), repeat=k):
        chosen: list[Edge] = []
        remaining: list[Edge] = []
        for (a, b), s in zip(undirected, state):  # noqa: B905
            if s == 1:
                chosen.append((a, b))
            elif s == 2:
                chosen.append((b, a))
            else:
                remaining.append((a, b))
        try:
            candidate = MPDAG(cpdag.nodes, base_directed | set(chosen), remaining)
        except ValueError:
            # A chosen orientation contradicted another (shouldn't happen,
            # since each edge is assigned independently, but stay defensive).
            continue
        if not is_valid_mpdag(candidate):
            continue
        extensions = enumerate_dag_extensions(candidate)
        if not all(is_consistent_extension(d, cpdag) for d in extensions):
            continue
        space.append(candidate)

    space.sort(key=lambda g: (len(g.directed_edges), g.edge_string()))
    return space


def represented_dags(space: list[MPDAG]) -> dict[MPDAG, frozenset[MPDAG]]:
    """``[G]`` for every ``G`` in ``space``, computed once and cached.

    Args:
        space: The enumerated space, as returned by :func:`enumerate_space`.

    Returns:
        A mapping from each ``G`` to the frozenset of its consistent DAG
        extensions (:func:`~bkrobust.demo.meek.enumerate_dag_extensions`).
        Every other function in this module that needs ``[G]`` takes this
        mapping as an argument rather than recomputing it.
    """
    return {g: frozenset(enumerate_dag_extensions(g)) for g in space}


# --- model inclusion and covering --------------------------------------------


def covering_pairs(
    space: list[MPDAG], reps: dict[MPDAG, frozenset[MPDAG]]
) -> set[tuple[MPDAG, MPDAG]]:
    """Covering pairs of the model-inclusion order, found by brute force.

    ``lower`` is covered by ``upper`` iff ``reps[lower]`` is a *strict*
    subset of ``reps[upper]`` (model inclusion, ``lower`` the more-refined,
    smaller model) and no third element of ``space`` has an extension set
    strictly between the two.

    No shortcut characterisation of covering (e.g. "differs by exactly one
    edge orientation") is used. That characterisation is false in general:
    orienting one edge can force further orientations via Meek propagation,
    so two MPDAGs that cover each other in the model-inclusion order can
    legitimately differ in several edges at once. Covering here is computed
    exactly as defined -- a subset check plus an exhaustive nothing-strictly-
    between search over every other element of the space.

    Args:
        space: The enumerated space.
        reps: ``[G]`` for each ``G`` in ``space``, from
            :func:`represented_dags`.

    Returns:
        The set of ``(lower, upper)`` pairs with ``lower`` covered by
        ``upper``.
    """
    covers: set[tuple[MPDAG, MPDAG]] = set()
    for lower in space:
        lower_ext = reps[lower]
        for upper in space:
            if lower is upper:
                continue
            upper_ext = reps[upper]
            if not (lower_ext < upper_ext):
                continue
            between = False
            for h in space:
                if h is lower or h is upper:
                    continue
                h_ext = reps[h]
                if lower_ext < h_ext < upper_ext:
                    between = True
                    break
            if not between:
                covers.add((lower, upper))
    return covers


# --- neighbour graph and distances -------------------------------------------


def neighbour_graph(
    space: list[MPDAG], covers: set[tuple[MPDAG, MPDAG]]
) -> dict[MPDAG, set[MPDAG]]:
    """The undirected adjacency graph of the covering relation.

    Args:
        space: The enumerated space (defines the full vertex set, including
            any vertex with no covering edges).
        covers: Covering pairs, from :func:`covering_pairs`.

    Returns:
        A symmetric adjacency mapping: ``g in nbrs[h]`` iff ``h in
        nbrs[g]`` iff ``(g, h)`` or ``(h, g)`` is in ``covers``.
    """
    nbrs: dict[MPDAG, set[MPDAG]] = {g: set() for g in space}
    for lower, upper in covers:
        nbrs[lower].add(upper)
        nbrs[upper].add(lower)
    return nbrs


def bfs_distances(nbrs: dict[MPDAG, set[MPDAG]], source: MPDAG) -> dict[MPDAG, int]:
    """BFS hop distances from ``source`` on the neighbour graph.

    Args:
        nbrs: Symmetric adjacency, from :func:`neighbour_graph`.
        source: The starting vertex.

    Returns:
        A mapping from every vertex reachable from ``source`` (``source``
        itself included, at distance 0) to its hop count. Unreachable
        vertices are simply absent from the mapping -- never recorded as
        infinity.
    """
    if source not in nbrs:
        return {}
    dist: dict[MPDAG, int] = {source: 0}
    queue: deque[MPDAG] = deque([source])
    while queue:
        cur = queue.popleft()
        for nxt in nbrs[cur]:
            if nxt not in dist:
                dist[nxt] = dist[cur] + 1
                queue.append(nxt)
    return dist


def all_pairs_distances(nbrs: dict[MPDAG, set[MPDAG]]) -> dict[tuple[MPDAG, MPDAG], int]:
    """BFS distance for every ordered pair reachable from one another.

    Args:
        nbrs: Symmetric adjacency, from :func:`neighbour_graph`.

    Returns:
        A mapping ``{(a, b): dist}`` covering ``(a, a): 0`` for every vertex
        and both ``(a, b)`` and ``(b, a)`` (necessarily equal, since the
        neighbour graph is undirected) whenever ``b`` is reachable from
        ``a``. Unreachable pairs are absent, matching
        :func:`bfs_distances`'s convention.
    """
    result: dict[tuple[MPDAG, MPDAG], int] = {}
    for a in nbrs:
        for b, d in bfs_distances(nbrs, a).items():
            result[(a, b)] = d
    return result


# --- metric verification ------------------------------------------------------


def check_metric_axioms(
    dists: dict[tuple[MPDAG, MPDAG], int], space: list[MPDAG]
) -> dict[str, object]:
    """Verify BFS distance is a genuine metric on ``space``, checking all triples.

    Because :func:`bfs_distances` omits unreachable pairs rather than
    recording infinity, a pairwise distance can be undefined; a triple only
    counts toward ``n_triples_checked`` when all three of its pairwise
    distances are defined (undefined distances cannot be compared).

    Args:
        dists: All-pairs distances, from :func:`all_pairs_distances`.
        space: The enumerated space.

    Returns:
        A dict with:

        * ``identity_of_indiscernibles`` (bool): ``d(x, x) == 0`` for every
          ``x`` in ``space``, and ``d(x, y) == 0`` implies ``x == y`` for
          every defined pair.
        * ``symmetry`` (bool): ``d(x, y) == d(y, x)`` for every defined pair.
        * ``triangle_inequality`` (bool): ``d(x, z) <= d(x, y) + d(y, z)``
          for every triple with all three distances defined.
        * ``n_triples_checked`` (int): how many of the ``len(space) ** 3``
          ordered triples (with repetition) had all three distances defined
          and were actually evaluated for the triangle inequality.
        * ``violations`` (list[str]): human-readable descriptions of every
          failure found, empty if none.
    """

    def d(x: MPDAG, y: MPDAG) -> int | None:
        if x == y:
            return 0
        return dists.get((x, y))

    violations: list[str] = []
    identity_ok = True
    symmetry_ok = True
    triangle_ok = True

    for x in space:
        if d(x, x) != 0:
            identity_ok = False
            violations.append(f"identity: d({x.edge_string()!r}, self) != 0")

    for (a, b), dab in dists.items():
        if a != b and dab == 0:
            identity_ok = False
            violations.append(
                f"identity: d({a.edge_string()!r}, {b.edge_string()!r}) == 0 for distinct a, b"
            )

    for (a, b), dab in dists.items():
        dba = dists.get((b, a))
        if dba is None or dba != dab:
            symmetry_ok = False
            violations.append(
                f"symmetry: d({a.edge_string()!r},{b.edge_string()!r})={dab} "
                f"!= d({b.edge_string()!r},{a.edge_string()!r})={dba}"
            )

    n_triples_checked = 0
    for x in space:
        for y in space:
            dxy = d(x, y)
            if dxy is None:
                continue
            for z in space:
                dyz = d(y, z)
                dxz = d(x, z)
                if dyz is None or dxz is None:
                    continue
                n_triples_checked += 1
                if dxz > dxy + dyz:
                    triangle_ok = False
                    violations.append(
                        f"triangle: d({x.edge_string()!r},{z.edge_string()!r})={dxz} > "
                        f"d(x,y)={dxy} + d(y,z)={dyz} "
                        f"(y={y.edge_string()!r})"
                    )

    return {
        "identity_of_indiscernibles": identity_ok,
        "symmetry": symmetry_ok,
        "triangle_inequality": triangle_ok,
        "n_triples_checked": n_triples_checked,
        "violations": violations,
    }


# --- atomic moves --------------------------------------------------------


def _edge_state(g: MPDAG, a: Node, b: Node) -> Edge | None:
    """The oriented direction of the ``a``/``b`` edge in ``g``, or ``None`` if undirected."""
    if g.is_directed_edge(a, b):
        return (a, b)
    if g.is_directed_edge(b, a):
        return (b, a)
    return None


def atomic_moves(g_from: MPDAG, g_to: MPDAG) -> list[str]:
    """Human-readable atomic move sequence transforming ``g_from`` into ``g_to``.

    ``g_from`` and ``g_to`` need not be adjacent in the neighbour graph. Each
    edge of their shared skeleton is in exactly one of four situations,
    contributing zero, one, or two moves:

    * unchanged orientation state -- no move;
    * undirected in ``g_from``, directed in ``g_to`` -- one assertion,
      ``"orient A->B"``;
    * directed in ``g_from``, undirected in ``g_to`` -- one retraction,
      ``"un-orient A-B"``;
    * directed in ``g_from``, directed the *opposite* way in ``g_to`` (a
      flip) -- exactly **two** moves, the retraction followed by the new
      assertion. A flip is not a primitive move: an already-oriented edge
      must first be un-oriented before it can be asserted the other way, so
      it always costs two atomic moves, never one.

    Args:
        g_from: The starting graph.
        g_to: The target graph. Must share ``g_from``'s skeleton.

    Returns:
        The move sequence, one string per atomic step, in a deterministic
        (skeleton-sorted) order.

    Raises:
        ValueError: If ``g_from`` and ``g_to`` do not share a skeleton.
    """
    if g_from.skeleton() != g_to.skeleton():
        raise ValueError("atomic_moves requires g_from and g_to to share a skeleton")

    moves: list[str] = []
    for a, b in sorted(g_from.skeleton()):
        state_from = _edge_state(g_from, a, b)
        state_to = _edge_state(g_to, a, b)
        if state_from == state_to:
            continue
        if state_from is not None:
            moves.append(f"un-orient {a}-{b}")
        if state_to is not None:
            tail, head = state_to
            moves.append(f"orient {tail}->{head}")
    return moves
