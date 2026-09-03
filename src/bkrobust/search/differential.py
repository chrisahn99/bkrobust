"""Differential testing: every accelerated method must reproduce exact BFS radii.

Any disagreement is a bug, not a tuning parameter. The harness records a minimal
counterexample when one appears so it can be turned into a regression test.

Cost is reported both as wall-clock and as machine-independent counters
(elements visited, Meek closures, validity checks), because wall-clock alone
does not survive a change of machine and is not a usable claim in a paper.
"""

from __future__ import annotations

import itertools
import time
from dataclasses import dataclass, field
from typing import Any

from bkrobust.core.oracle import is_valid
from bkrobust.core.spacelib import build_space, distances_from, radius
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.search.conjecture_study import all_cpdags
from bkrobust.search.exact import SearchStats, radius_local_up


@dataclass
class DiffTotals:
    """Totals for one differential sweep."""

    n_cases: int = 0
    n_agree: int = 0
    n_disagree: int = 0
    bfs_seconds: float = 0.0
    fast_seconds: float = 0.0
    bfs_space_elements: int = 0
    fast_elements_visited: int = 0
    fast_closures: int = 0
    fast_validity_checks: int = 0
    disagreements: list[dict[str, Any]] = field(default_factory=list)


def differential_sweep(
    n: int,
    *,
    max_undirected: int = 6,
    max_cases: int | None = None,
    per_cpdag_targets: int | None = None,
) -> tuple[DiffTotals, list[dict[str, Any]]]:
    """Race ``radius_local_up`` against full-space BFS on every feasible instance.

    Args:
        n: Node count.
        max_undirected: Skip CPDAGs above this undirected-edge count.
        max_cases: Stop after this many comparisons.
        per_cpdag_targets: Cap on ``(G0, X, Y)`` combinations per CPDAG.

    Returns:
        ``(totals, rows)`` where each row records both radii and both costs.
    """
    totals = DiffTotals()
    rows: list[dict[str, Any]] = []

    for cpdag in all_cpdags(n):
        if not cpdag.undirected_edges or len(cpdag.undirected_edges) > max_undirected:
            continue
        space = build_space(cpdag)
        totals.bfs_space_elements += len(space)
        nodes = list(cpdag.nodes)
        combos = 0
        for g0 in space.elements:
            dists = distances_from(space, g0)
            for x, y in itertools.permutations(nodes, 2):
                o = optimal_adjustment_set_mpdag(g0, x, y)
                if o is None:
                    continue
                z = frozenset(o)
                if not is_valid(z, g0, x, y):
                    continue

                def fails(g: MPDAG, _z: frozenset[str] = z, _x: str = x, _y: str = y) -> bool:
                    return not is_valid(_z, g, _x, _y)

                t0 = time.perf_counter()
                r_bfs, w_bfs = radius(space, dists, fails)
                t1 = time.perf_counter()
                st = SearchStats()
                res = radius_local_up(cpdag, g0, fails, stats=st)
                t2 = time.perf_counter()

                totals.n_cases += 1
                combos += 1
                totals.bfs_seconds += t1 - t0
                totals.fast_seconds += t2 - t1
                totals.fast_elements_visited += st.elements_visited
                totals.fast_closures += st.closures
                totals.fast_validity_checks += st.validity_checks

                agree = res.radius == r_bfs
                if agree:
                    totals.n_agree += 1
                else:
                    totals.n_disagree += 1
                    totals.disagreements.append(
                        {
                            "n": n,
                            "cpdag": cpdag.edge_string(),
                            "g0": g0.edge_string(),
                            "x": x,
                            "y": y,
                            "z": sorted(z),
                            "r_bfs": r_bfs,
                            "r_fast": res.radius,
                            "witness_bfs": None if w_bfs is None else w_bfs.edge_string(),
                            "witness_fast": res.witness,
                            "space_size": len(space),
                        }
                    )
                rows.append(
                    {
                        "n": n,
                        "cpdag": cpdag.edge_string(),
                        "n_undirected": len(cpdag.undirected_edges),
                        "space_size": len(space),
                        "g0": g0.edge_string(),
                        "x": x,
                        "y": y,
                        "z": "|".join(sorted(z)),
                        "r_bfs": r_bfs,
                        "r_fast": res.radius,
                        "agree": agree,
                        "bfs_seconds": t1 - t0,
                        "fast_seconds": t2 - t1,
                        "fast_elements_visited": st.elements_visited,
                        "fast_closures": st.closures,
                        "fast_validity_checks": st.validity_checks,
                    }
                )
                if max_cases is not None and totals.n_cases >= max_cases:
                    return totals, rows
                if per_cpdag_targets is not None and combos >= per_cpdag_targets:
                    break
            if per_cpdag_targets is not None and combos >= per_cpdag_targets:
                break
    return totals, rows
