"""Budgeted upward BFS, vendored from E1 of the same rebuttal.

Source: experiments/rebuttal/E1_cpdag_check_on_lower_bound_rows/search.py
commit 5421cabc23b28af8f7c6e4124a698d66a76e0a08
(the source file is untracked at that commit; copied unchanged below the next line on 2026-11-06).

The traversal is the submitted ``radius_local_up_fast`` step for step, over the submitted
``local_up_covers_fast``; it adds the closure budget and wall-clock limit that the submitted
function lacks and returns the sound lower bound ``lb`` when it stops early. E10 uses it on
networks above 60 nodes, where the submitted hybrid would enter the SAT ladder, and as the
independent cross-check of every ladder answer.
"""

from __future__ import annotations

import time
from collections import deque
from collections.abc import Callable

from bkrobust.core.conventions import UNREACHED
from bkrobust.demo.graph import MPDAG
from bkrobust.search.exact import SearchStats
from bkrobust.search.exact_fast import local_up_covers_fast


def bfs_radius(
    cpdag: MPDAG,
    g0: MPDAG,
    fails: Callable[[MPDAG], bool],
    *,
    max_depth: int | None,
    max_closures: int | None,
    time_limit_s: float | None,
) -> dict:
    """Radius by budgeted upward BFS.

    Returns a dict with ``radius`` (int, or ``UNREACHED``), ``exact``, ``stop``
    (``found``, ``exhausted``, ``depth``, ``closures``, ``time``), ``lb`` and
    ``lb_release_formula`` (None when exact), the witness edge string, the
    counters and the seconds spent.
    """
    t0 = time.perf_counter()
    deadline = None if time_limit_s is None else t0 + time_limit_s
    st = SearchStats()

    def out(radius, exact, stop, lb=None, lb_rel=None, witness=None):
        return {
            "radius": radius,
            "exact": exact,
            "stop": stop,
            "lb": lb,
            "lb_release_formula": lb_rel,
            "witness": witness,
            "max_depth": max_depth,
            "max_closures": max_closures,
            "time_limit_s": time_limit_s,
            "seconds": round(time.perf_counter() - t0, 4),
            **st.as_dict(),
        }

    st.validity_checks += 1
    if fails(g0):
        return out(0, True, "found", witness=g0.edge_string())

    seen = {g0.edge_string()}
    frontier: deque[tuple[MPDAG, int]] = deque([(g0, 0)])
    max_popped = 0
    min_skipped: int | None = None  # smallest depth of an element left unexpanded
    stop = None
    while frontier:
        cur, d = frontier.popleft()
        st.elements_visited += 1
        max_popped = max(max_popped, d)
        if max_depth is not None and d >= max_depth:
            # Skipped, but the walk continues: shallower elements may remain.
            min_skipped = d if min_skipped is None else min(min_skipped, d)
            stop = stop or "depth"
            continue
        out_of_closures = max_closures is not None and st.closures >= max_closures
        out_of_time = deadline is not None and time.perf_counter() >= deadline
        if out_of_closures or out_of_time:
            # Nothing more is expanded: cur and the rest of the queue stay unexpanded.
            # FIFO order means cur has the smallest depth among them.
            min_skipped = d if min_skipped is None else min(min_skipped, d)
            # The release formula would have popped every remaining element too.
            max_popped = max([max_popped, *(dd for _, dd in frontier)])
            stop = "closures" if out_of_closures else "time"
            break
        for h in local_up_covers_fast(cpdag, cur, st):
            key = h.edge_string()
            if key in seen:
                continue
            seen.add(key)
            st.validity_checks += 1
            if fails(h):
                return out(d + 1, True, "found", witness=key)
            frontier.append((h, d + 1))

    if min_skipped is None:
        return out(UNREACHED, True, "exhausted")
    return out(UNREACHED, False, stop, lb=min_skipped + 1, lb_rel=max_popped + 1)
