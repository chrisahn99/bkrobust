"""Cover generation without enumerating DAG extensions.

Session 4 profiled ``radius_local_up`` and found it **96.6% enumeration-bound**
— but not where the plan expected. The validity oracle is only **23.0%** of its
time; **73.6%** is spent inside :func:`~bkrobust.search.exact.local_up_covers`,
which calls ``enumerate_dag_extensions`` once per candidate purely to decide
minimality:

.. code-block:: python

    reps[h] = frozenset(enumerate_dag_extensions(h))  # the 73.6%
    if not any(reps[other] < reps[h] for other in graphs):
        covers.append(h)

That comparison asks only whether one candidate is strictly below another in
**model inclusion**. It does not need the extension sets themselves, because on
elements of the space the order is a set comparison on directed edges:

    Lemma O.  For elements ``G``, ``H`` of the space:  [G] ⊆ [H]  ⟺  dir(H) ⊆ dir(G).

*Proof.* (⇐) If ``dir(H) ⊆ dir(G)`` then every DAG consistent with ``G`` is
consistent with the weaker constraint set ``H``, so ``[G] ⊆ [H]``.
(⇒) Let ``a→b ∈ dir(H)``. Every ``D ∈ [H]`` orients it ``a→b``, and
``[G] ⊆ [H]``, so every ``D ∈ [G]`` does too. Elements of the space are
**maximally oriented** — an edge oriented identically across all of ``[G]`` is
directed in ``G`` (this is the defining property established in §1 of
``THEOREMS.md`` and checked by ``space_fixed.is_maximally_oriented``) — so
``a→b ∈ dir(G)``. ∎

Both directions are strict together, since the correspondence between an element
and its extension set is injective. So

    reps[other] < reps[h]   ⟺   dir(h) ⊊ dir(other),

and cover generation becomes pure set arithmetic.

:func:`local_up_covers_fast` is a drop-in replacement returning the same covers
in the same order. ``exact.py`` is left untouched, so prior sessions' results
stay reproducible from the code that produced them, and the two implementations
can be differentially tested against each other.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Callable

from bkrobust.core.conventions import UNREACHED
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import meek_closure
from bkrobust.search.exact import RadiusResult, SearchStats, retractable_edges


def local_up_covers_fast(cpdag: MPDAG, g: MPDAG, stats: SearchStats | None = None) -> list[MPDAG]:
    """Upper covers of ``g``, with minimality decided by set comparison.

    Identical in output to :func:`bkrobust.search.exact.local_up_covers`,
    including ordering, but performs no DAG-extension enumeration. See the module
    docstring for Lemma O, which licenses the substitution.

    Args:
        cpdag: The CPDAG, which fixes which edges are frozen.
        g: The graph to step up from.
        stats: Optional counters to increment. ``extensions_enumerated`` is
            deliberately **not** incremented, because none are.

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

    dirs = [frozenset(h.directed_edges) for h in graphs]
    covers = []
    for i, h in enumerate(graphs):
        # h is a cover unless another candidate sits strictly between g and h,
        # i.e. unless some other candidate is strictly MORE oriented than h
        # while still containing h's orientations.
        if not any(j != i and dirs[i] < dirs[j] for j in range(len(graphs))):
            covers.append(h)
    return covers


def radius_local_up_fast(
    cpdag: MPDAG,
    g0: MPDAG,
    fails: Callable[[MPDAG], bool],
    *,
    max_depth: int | None = None,
    stats: SearchStats | None = None,
) -> RadiusResult:
    """Exact radius by upward BFS, using enumeration-free cover generation.

    Semantically identical to :func:`bkrobust.search.exact.radius_local_up`,
    including its exactness caveat: it is exact **iff Conjecture 2 holds**, and
    its error is one-sided — it can only ever return radii that are too large,
    never too small.

    Args:
        cpdag: The CPDAG.
        g0: The analyst's graph.
        fails: Predicate that is ``True`` where the property fails.
        max_depth: Stop after this many retractions; on exhaustion the result is
            the anytime statement ``radius >= max_depth + 1``, with
            ``exact=False``.
        stats: Optional counters.

    Returns:
        A :class:`~bkrobust.search.exact.RadiusResult` labelled
        ``local_up_fast_exact``.
    """
    st = stats if stats is not None else SearchStats()
    st.validity_checks += 1
    if fails(g0):
        return RadiusResult(0, g0.edge_string(), "local_up_fast_exact", True, st)

    seen = {g0.edge_string()}
    frontier = deque([(g0, 0)])
    budget_hit = False
    while frontier:
        cur, d = frontier.popleft()
        st.elements_visited += 1
        if max_depth is not None and d >= max_depth:
            budget_hit = True
            continue
        for h in local_up_covers_fast(cpdag, cur, st):
            key = h.edge_string()
            if key in seen:
                continue
            seen.add(key)
            st.validity_checks += 1
            if fails(h):
                return RadiusResult(d + 1, key, "local_up_fast_exact", True, st)
            frontier.append((h, d + 1))
    if budget_hit:
        return RadiusResult(UNREACHED, None, "local_up_fast_budget", False, st)
    return RadiusResult(UNREACHED, None, "local_up_fast_exact", True, st)
