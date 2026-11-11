"""The upward search of the submitted hybrid, with a closure budget and a wall clock.

``bkrobust.search.exact_fast.radius_local_up_fast`` has no budget other than depth. On graphs of
700 to 1041 nodes one cover step can cost hundreds of Meek closures, so an unbudgeted depth-3
search may not return. This is the same breadth-first search over the same covers, in the same
order (``local_up_covers_fast`` from the submitted tree), so where it finds a failure it finds the
one the hybrid's search leg finds. It stops when the closure count or the clock runs out, and
then certifies only what it checked.

What a stop certifies: BFS pops states in non-decreasing depth, and a state is checked when it is
generated, so when the search stops before expanding a state at depth d, every state at depth
<= d has been checked and none fails: r >= d + 1.
"""

from __future__ import annotations

import time
from collections import deque
from dataclasses import dataclass, field

from bkrobust.core.conventions import UNREACHED
from bkrobust.demo.graph import MPDAG
from bkrobust.search.exact import SearchStats
from bkrobust.search.exact_fast import local_up_covers_fast


@dataclass
class BudgetedResult:
    radius: int  # exact radius, or UNREACHED
    exact: bool
    stop: str  # found | degenerate | exhausted | depth_budget | closure_budget | time_budget
    lower_bound: int | None  # r >= lower_bound when not exact
    witness: str | None = None
    depth_reached: int = 0
    seconds: float = 0.0
    stats: SearchStats = field(default_factory=SearchStats)


def budgeted_up_search(
    cpdag: MPDAG,
    g0: MPDAG,
    fails,
    *,
    max_depth: int,
    closure_budget: int,
    time_limit_s: float,
) -> BudgetedResult:
    t0 = time.perf_counter()
    st = SearchStats()
    st.validity_checks += 1
    if fails(g0):
        return BudgetedResult(0, True, "degenerate", None, g0.edge_string(), 0,
                              time.perf_counter() - t0, st)
    seen = {g0.edge_string()}
    frontier = deque([(g0, 0)])
    depth_hit = False
    deepest = 0
    while frontier:
        cur, d = frontier.popleft()
        deepest = max(deepest, d)
        if st.closures >= closure_budget or time.perf_counter() - t0 > time_limit_s:
            stop = "closure_budget" if st.closures >= closure_budget else "time_budget"
            return BudgetedResult(UNREACHED, False, stop, d + 1, None, d,
                                  time.perf_counter() - t0, st)
        st.elements_visited += 1
        if d >= max_depth:
            # mirrors radius_local_up_fast: only an unseen cover leaves anything unexplored
            if not depth_hit and any(
                h.edge_string() not in seen for h in local_up_covers_fast(cpdag, cur, st)
            ):
                depth_hit = True
            continue
        for h in local_up_covers_fast(cpdag, cur, st):
            key = h.edge_string()
            if key in seen:
                continue
            seen.add(key)
            st.validity_checks += 1
            if fails(h):
                return BudgetedResult(d + 1, True, "found", None, key, d + 1,
                                      time.perf_counter() - t0, st)
            frontier.append((h, d + 1))
    if depth_hit:
        return BudgetedResult(UNREACHED, False, "depth_budget", max_depth + 1, None, deepest,
                              time.perf_counter() - t0, st)
    return BudgetedResult(UNREACHED, True, "exhausted", None, None, deepest,
                          time.perf_counter() - t0, st)
