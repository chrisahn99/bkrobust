"""Accelerated exact radius computation, and the bounds that certify it.

Shell-by-shell BFS requires the whole space, which costs ``3^m`` in the number
of undirected edges. Two structural facts let us do better:

* **Failure is upward-closed** (a theorem; see :mod:`bkrobust.search.conjectures`).
  Once a graph fails, everything above it fails, so the search never needs to
  expand above a known failure.
* **Conjecture 2** -- empirically supported, not proved -- says the nearest
  failure lies in the *up-set* of ``G0``. If so the search is over retraction
  subsets, roughly ``2^k`` in the number of knowledge-oriented edges, and it
  never needs the rest of the space.

Method labels are carried on every result so no number is silently attributed to
an unverified shortcut:

``bfs_exact``
    Ground truth: full space, full BFS. Always correct.
``up_bfs_exact``
    BFS restricted to upward covers, still using the enumerated space. Exact
    **iff** Conjecture 2 holds; differentially tested against ``bfs_exact``.
``local_up_exact``
    As above but generating upper covers locally, so the space is never
    enumerated. This is where the real speedup lives. Exact iff Conjecture 2
    holds AND local cover generation is complete -- both differentially tested.
"""

from __future__ import annotations

import itertools
from collections import deque
from collections.abc import Callable
from dataclasses import dataclass, field

from bkrobust.core.conventions import UNREACHED
from bkrobust.demo.graph import MPDAG, canon
from bkrobust.demo.meek import enumerate_dag_extensions, meek_closure


@dataclass
class SearchStats:
    """Machine-independent cost counters, reported alongside wall-clock.

    Attributes:
        elements_visited: Graphs dequeued.
        closures: Meek closures computed.
        validity_checks: Calls to the validity oracle.
        extensions_enumerated: DAG-extension enumerations performed.
    """

    elements_visited: int = 0
    closures: int = 0
    validity_checks: int = 0
    extensions_enumerated: int = 0

    def as_dict(self) -> dict[str, int]:
        """Flatten for a results row."""
        return {
            "elements_visited": self.elements_visited,
            "closures": self.closures,
            "validity_checks": self.validity_checks,
            "extensions_enumerated": self.extensions_enumerated,
        }


@dataclass
class RadiusResult:
    """A radius together with how it was obtained.

    Attributes:
        radius: Under the normative convention, or ``UNREACHED``.
        witness: A nearest failing graph, as an edge string.
        method: Which method produced it.
        exact: Whether the value is exact or a labelled bound.
        stats: Cost counters.
    """

    radius: int
    witness: str | None
    method: str
    exact: bool
    stats: SearchStats = field(default_factory=SearchStats)


def retractable_edges(cpdag: MPDAG, g: MPDAG) -> list[tuple[str, str]]:
    """Edges of ``g`` that knowledge oriented, and which may therefore be retracted.

    An edge is retractable iff it is directed in ``g`` and undirected in the
    CPDAG. Data-compelled edges are frozen and never appear here.

    Returns:
        ``(tail, head)`` pairs in a deterministic sorted order.
    """
    return sorted((a, b) for (a, b) in g.directed_edges if canon(a, b) in cpdag.undirected_edges)


def local_up_covers(cpdag: MPDAG, g: MPDAG, stats: SearchStats | None = None) -> list[MPDAG]:
    """Upper covers of ``g``, generated locally without enumerating the space.

    Candidates come from un-orienting one knowledge-oriented edge and re-closing
    under Meek. Every graph strictly above ``g`` is reachable by retractions, so
    this candidate set contains every cover; the minimal candidates are then
    kept, which is what makes them covers rather than merely upper bounds.

    Args:
        cpdag: The CPDAG, which fixes which edges are frozen.
        g: The graph to step up from.
        stats: Optional counters to increment.

    Returns:
        The upper covers, deterministically ordered.
    """
    cands: dict[str, MPDAG] = {}
    for a, b in retractable_edges(cpdag, g):
        h = meek_closure(g.unoriented(a, b))
        if stats is not None:
            stats.closures += 1
        if h is None or h == g:
            continue
        cands[h.edge_string()] = h
    graphs = [cands[k] for k in sorted(cands)]

    # Keep only the minimal candidates: h is a cover unless another candidate
    # sits strictly between g and h.
    reps = {}
    for h in graphs:
        reps[h] = frozenset(enumerate_dag_extensions(h))
        if stats is not None:
            stats.extensions_enumerated += 1
    covers = []
    for h in graphs:
        if not any(other is not h and reps[other] < reps[h] for other in graphs):
            covers.append(h)
    return covers


def radius_local_up(
    cpdag: MPDAG,
    g0: MPDAG,
    fails: Callable[[MPDAG], bool],
    *,
    max_depth: int | None = None,
    stats: SearchStats | None = None,
) -> RadiusResult:
    """Exact radius by upward BFS with locally generated covers -- no space enumeration.

    Exploits upward-closure by never expanding above a graph already known to
    fail: everything above it fails too, so it can contribute no nearer witness.

    Args:
        cpdag: The CPDAG.
        g0: The analyst's graph.
        fails: Predicate that is ``True`` where the property fails.
        max_depth: Stop after this many retractions. On exhaustion of the budget
            without a failure the result is the **anytime** statement
            ``radius >= max_depth + 1``, returned with ``exact=False``.
        stats: Optional counters.

    Returns:
        A :class:`RadiusResult` labelled ``local_up_exact``.
    """
    st = stats if stats is not None else SearchStats()
    if fails(g0):
        st.validity_checks += 1
        return RadiusResult(0, g0.edge_string(), "local_up_exact", True, st)
    st.validity_checks += 1

    seen = {g0.edge_string()}
    frontier = deque([(g0, 0)])
    budget_hit = False
    while frontier:
        cur, d = frontier.popleft()
        st.elements_visited += 1
        if max_depth is not None and d >= max_depth:
            budget_hit = True
            continue
        for h in local_up_covers(cpdag, cur, st):
            key = h.edge_string()
            if key in seen:
                continue
            seen.add(key)
            st.validity_checks += 1
            if fails(h):
                return RadiusResult(d + 1, key, "local_up_exact", True, st)
            # h is clean; only clean graphs are expanded. Anything above a
            # failing graph also fails and cannot be nearer, so it is pruned.
            frontier.append((h, d + 1))
    if budget_hit:
        return RadiusResult(UNREACHED, None, "local_up_budget", False, st)
    return RadiusResult(UNREACHED, None, "local_up_exact", True, st)


def descendant_lower_bound(
    cpdag: MPDAG,
    g0: MPDAG,
    z: frozenset[str],
    x: str,
) -> int:
    """Minimum retractions before any ``z in Z`` can become a descendant of ``X``.

    One of the two ways ``Z`` can stop being a valid adjustment set is for a
    member to enter the forbidden set, i.e. to become a possible descendant of
    the treatment. Realising that needs a directed path ``X -> ... -> z``, so
    every edge on some ``X``-to-``z`` path that currently points the wrong way
    must first be retracted. Shortest such path, over all ``z``, is the cost of
    this failure mode.

    Edge costs on a path from ``X`` towards ``z``: 0 if the edge is undirected or
    already points forward, 1 if it points backward and is retractable,
    unreachable if it points backward and is data-compelled.

    Args:
        cpdag: The CPDAG (fixes which edges are frozen).
        g0: The analyst's graph.
        z: The adjustment set.
        x: Treatment.

    Returns:
        The minimum retraction count for this failure mode, or
        :data:`~bkrobust.core.conventions.UNREACHED` if no ``z`` can be reached.

    Note:
        This is a bound for **one** failure mode. It is a lower bound on the
        radius only in combination with a bound for the other mode (a back-door
        path becoming unblocked); on its own it may exceed the true radius. It is
        exposed for analysis and for the ``L == U`` study, not as a standalone
        certificate.
    """
    import heapq

    retractable = set(retractable_edges(cpdag, g0))
    dist: dict[str, int] = {x: 0}
    heap: list[tuple[int, str]] = [(0, x)]
    best = UNREACHED
    while heap:
        d, node = heapq.heappop(heap)
        if d > dist.get(node, 1 << 30):
            continue
        if node in z and (best == UNREACHED or d < best):
            best = d
        for nxt in sorted(g0.adjacent(node)):
            if g0.is_directed_edge(node, nxt) or g0.is_undirected_edge(node, nxt):
                w = 0
            elif g0.is_directed_edge(nxt, node):
                w = 1 if (nxt, node) in retractable else None
            else:
                w = None
            if w is None:
                continue
            nd = d + w
            if nd < dist.get(nxt, 1 << 30):
                dist[nxt] = nd
                heapq.heappush(heap, (nd, nxt))
    return best


def all_retraction_subsets(cpdag: MPDAG, g0: MPDAG, k: int) -> list[tuple[tuple[str, str], ...]]:
    """All size-``k`` subsets of the retractable edges, deterministically ordered."""
    return list(itertools.combinations(retractable_edges(cpdag, g0), k))


def guided_upper_bound(
    cpdag: MPDAG,
    g0: MPDAG,
    fails: Callable[[MPDAG], bool],
    z: frozenset[str],
    x: str,
    *,
    max_tries: int = 64,
) -> RadiusResult:
    """Construct a failing graph directly and measure its distance: ``r <= U``.

    Rather than searching, this aims at a known failure mode -- making some
    ``z in Z`` a possible descendant of ``X`` -- by retracting the
    backward-pointing edges along a cheapest ``X``-to-``z`` route, then checking
    whether the result actually fails. Retraction sets are tried in increasing
    size, so the first success is the best this construction can offer.

    Args:
        cpdag: The CPDAG.
        g0: The analyst's graph.
        fails: Predicate that is ``True`` where the property fails.
        z: The adjustment set.
        x: Treatment.
        max_tries: Cap on candidate retraction sets examined.

    Returns:
        A :class:`RadiusResult` with ``exact=False``: the value is an **upper
        bound**, not a radius. ``UNREACHED`` means the construction found no
        failure, which does not mean none exists.
    """
    st = SearchStats()
    retr = retractable_edges(cpdag, g0)
    tries = 0
    for k in range(1, len(retr) + 1):
        for subset in itertools.combinations(retr, k):
            tries += 1
            if tries > max_tries:
                return RadiusResult(UNREACHED, None, "guided_upper_budget", False, st)
            g = g0
            ok = True
            for a, b in subset:
                if not g.is_directed_edge(a, b):
                    ok = False
                    break
                g = g.unoriented(a, b)
            if not ok:
                continue
            h = meek_closure(g)
            st.closures += 1
            if h is None:
                continue
            st.validity_checks += 1
            if fails(h):
                return RadiusResult(k, h.edge_string(), "guided_upper", False, st)
    return RadiusResult(UNREACHED, None, "guided_upper", False, st)
