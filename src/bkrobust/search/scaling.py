"""Scaling of the space-free search against full-space BFS, at larger graphs.

Two accountings are kept apart throughout, because they disagree and both are
real:

``single query``
    One radius query against a fresh CPDAG. BFS must build the space first;
    ``radius_local_up`` never builds one. This is the regime a practitioner is
    in when they certify one analysis.
``amortised``
    Many queries against one CPDAG, so BFS pays for the space once. Here the
    space-free method can be, and is, **slower**.

The prediction under test is that the two methods scale against *different*
governing parameters: BFS against ``3^k`` in the undirected-edge count ``k``,
and ``radius_local_up`` against the size of the **up-set** of ``G0``, which is
roughly ``2^(#knowledge-oriented edges)``. Both are measured directly rather
than assumed, and machine-independent counters are reported alongside
wall-clock, which does not survive a change of machine.
"""

from __future__ import annotations

import itertools
import time
import tracemalloc
from dataclasses import asdict, dataclass
from typing import Any

from bkrobust.core.oracle import clear_cache, is_valid
from bkrobust.core.spacelib import distances_from, radius
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.search.conjectures import up_distances
from bkrobust.search.exact import SearchStats, radius_local_up, retractable_edges
from bkrobust.search.space_fixed import build_space_correct


@dataclass
class ScalingRow:
    """One measured radius query, both methods."""

    cpdag: str
    n_nodes: int
    k_undirected: int
    space_size: int
    upset_size: int
    n_retractable: int
    build_seconds: float
    bfs_query_seconds: float
    fast_seconds: float
    bfs_single_query_seconds: float
    speedup_single_query: float
    speedup_amortised: float
    crossover_queries: float
    fast_elements_visited: int
    fast_closures: int
    fast_validity_checks: int
    fast_peak_kib: float
    r_bfs: int
    r_fast: int
    agree: bool


def measure_instance(
    cpdag: MPDAG,
    max_targets: int = 4,
) -> list[ScalingRow]:
    """Measure both methods on one CPDAG, over a few target pairs.

    Returns:
        One :class:`ScalingRow` per measured query.
    """
    clear_cache()
    t0 = time.perf_counter()
    space = build_space_correct(cpdag)
    build_s = time.perf_counter() - t0

    rows: list[ScalingRow] = []
    nodes = list(cpdag.nodes)
    done = 0
    for g0 in space.elements:
        if done >= max_targets:
            break
        dists = distances_from(space, g0)
        ups = up_distances(space, g0)
        for x, y in itertools.permutations(nodes, 2):
            if done >= max_targets:
                break
            o = optimal_adjustment_set_mpdag(g0, x, y)
            if o is None:
                continue
            z = frozenset(o)
            if not z or not is_valid(z, g0, x, y):
                continue

            def fails(g: MPDAG, _z: frozenset[str] = z, _x: str = x, _y: str = y) -> bool:
                return not is_valid(_z, g, _x, _y)

            t0 = time.perf_counter()
            r_bfs, _ = radius(space, dists, fails)
            bfs_q = time.perf_counter() - t0

            st = SearchStats()
            tracemalloc.start()
            t0 = time.perf_counter()
            res = radius_local_up(cpdag, g0, fails, stats=st)
            fast_s = time.perf_counter() - t0
            _cur, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()

            bfs_single = build_s + bfs_q
            # queries at which amortised BFS overtakes: build/(fast - bfs_query)
            denom = fast_s - bfs_q
            crossover = (build_s / denom) if denom > 0 else float("inf")
            rows.append(
                ScalingRow(
                    cpdag=cpdag.edge_string(),
                    n_nodes=len(nodes),
                    k_undirected=len(cpdag.undirected_edges),
                    space_size=len(space),
                    upset_size=len(ups),
                    n_retractable=len(retractable_edges(cpdag, g0)),
                    build_seconds=build_s,
                    bfs_query_seconds=bfs_q,
                    fast_seconds=fast_s,
                    bfs_single_query_seconds=bfs_single,
                    speedup_single_query=bfs_single / max(fast_s, 1e-9),
                    speedup_amortised=bfs_q / max(fast_s, 1e-9),
                    crossover_queries=crossover,
                    fast_elements_visited=st.elements_visited,
                    fast_closures=st.closures,
                    fast_validity_checks=st.validity_checks,
                    fast_peak_kib=peak / 1024.0,
                    r_bfs=r_bfs,
                    r_fast=res.radius,
                    agree=bool(res.radius == r_bfs),
                )
            )
            done += 1
    return rows


def row_dicts(rows: list[ScalingRow]) -> list[dict[str, Any]]:
    """Flatten rows for CSV writing."""
    return [asdict(r) for r in rows]
