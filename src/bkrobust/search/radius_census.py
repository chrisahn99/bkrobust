"""An exhaustive census of radii over ALL CPDAGs on n nodes.

Complements the sampled synthetic ensembles of Axis A. Where an ensemble answers
"what happens for graphs drawn like this", a census answers "what happens for
*every* graph this small" -- no generator bias, no sampling error. It is only
possible because the number of Markov equivalence classes is small for
``n <= 5`` (11, 185, 8782).

Directly addresses:

* **H1 (saturation)** -- the full distribution of ``r_val``, in particular how
  much mass sits at 1.
* **H2 (frontier)** -- how often some valid adjustment set is strictly more
  robust than the optimal one, ``max_Z r_val(Z) > r_val(O)``.

The census uses the exact BFS radius throughout (``method = "bfs_exact"``), so
no claim here rests on an accelerated method.
"""

from __future__ import annotations

import itertools
from collections import Counter
from dataclasses import dataclass, field
from typing import Any

from bkrobust.core.conventions import UNREACHED
from bkrobust.core.oracle import is_valid
from bkrobust.core.spacelib import build_space, distances_from, radius
from bkrobust.demo.evaluate import (
    all_valid_adjustment_sets_mpdag,
    optimal_adjustment_set_mpdag,
)
from bkrobust.demo.graph import MPDAG
from bkrobust.search.conjecture_study import all_cpdags


@dataclass
class CensusTotals:
    """Aggregates for one census sweep."""

    n_nodes: int = 0
    n_cpdags_total: int = 0
    n_cpdags_used: int = 0
    n_instances: int = 0
    r_val_hist: dict[int, int] = field(default_factory=dict)
    r_val_unreached: int = 0
    n_frontier_checked: int = 0
    n_frontier_strict: int = 0
    frontier_gap_hist: dict[int, int] = field(default_factory=dict)
    skipped_no_optimal: int = 0
    skipped_z_invalid_at_g0: int = 0
    skipped_empty_z: int = 0

    def as_dict(self) -> dict[str, Any]:
        """JSON-friendly view with histogram keys as strings."""
        d = dict(self.__dict__)
        d["r_val_hist"] = {str(k): v for k, v in sorted(self.r_val_hist.items())}
        d["frontier_gap_hist"] = {str(k): v for k, v in sorted(self.frontier_gap_hist.items())}
        return d


def census(
    n: int,
    *,
    max_undirected: int = 6,
    frontier_max_sets: int = 12,
    exclude_empty_z: bool = True,
) -> tuple[CensusTotals, list[dict[str, Any]]]:
    """Compute the exact radius for every (CPDAG, G0, X, Y) on ``n`` nodes.

    Args:
        n: Node count.
        max_undirected: Skip CPDAGs above this undirected-edge count, since
            space construction is ``3^k``. Skips are counted.
        frontier_max_sets: Cap on valid adjustment sets examined per instance
            for the H2 frontier statistic.
        exclude_empty_z: Skip instances whose optimal adjustment set is empty.
            An empty ``Z`` cannot become invalid by acquiring a descendant of
            ``X``, so including them would bias the radius distribution upward
            for reasons unrelated to the phenomenon of interest. Counted
            separately.

    Returns:
        ``(totals, rows)``, one row per instance.
    """
    totals = CensusTotals(n_nodes=n)
    rows: list[dict[str, Any]] = []
    cpdags = all_cpdags(n)
    totals.n_cpdags_total = len(cpdags)

    for cpdag in cpdags:
        if not cpdag.undirected_edges or len(cpdag.undirected_edges) > max_undirected:
            continue
        space = build_space(cpdag)
        totals.n_cpdags_used += 1
        nodes = list(cpdag.nodes)
        for g0 in space.elements:
            dists = distances_from(space, g0)
            for x, y in itertools.permutations(nodes, 2):
                o = optimal_adjustment_set_mpdag(g0, x, y)
                if o is None:
                    totals.skipped_no_optimal += 1
                    continue
                z = frozenset(o)
                if exclude_empty_z and not z:
                    totals.skipped_empty_z += 1
                    continue
                if not is_valid(z, g0, x, y):
                    totals.skipped_z_invalid_at_g0 += 1
                    continue

                def fails(g: MPDAG, _z: frozenset[str] = z, _x: str = x, _y: str = y) -> bool:
                    return not is_valid(_z, g, _x, _y)

                r_opt_set, _ = radius(space, dists, fails)
                totals.n_instances += 1
                if r_opt_set == UNREACHED:
                    totals.r_val_unreached += 1
                else:
                    totals.r_val_hist[r_opt_set] = totals.r_val_hist.get(r_opt_set, 0) + 1

                # H2: is any valid adjustment set strictly more robust than O?
                best_r, best_z = r_opt_set, z
                cands = all_valid_adjustment_sets_mpdag(g0, x, y)[:frontier_max_sets]
                for cand in cands:
                    zc = frozenset(cand)
                    if zc == z or (exclude_empty_z and not zc):
                        continue
                    if not is_valid(zc, g0, x, y):
                        continue

                    def failsc(g: MPDAG, _z: frozenset[str] = zc, _x: str = x, _y: str = y) -> bool:
                        return not is_valid(_z, g, _x, _y)

                    rc, _ = radius(space, dists, failsc)
                    if rc == UNREACHED:
                        best_r, best_z = UNREACHED, zc
                        break
                    if best_r != UNREACHED and rc > best_r:
                        best_r, best_z = rc, zc
                totals.n_frontier_checked += 1
                if best_r == UNREACHED and r_opt_set != UNREACHED:
                    strict, gap = True, 999
                elif best_r != UNREACHED and r_opt_set != UNREACHED:
                    gap = best_r - r_opt_set
                    strict = gap > 0
                else:
                    strict, gap = False, 0
                if strict:
                    totals.n_frontier_strict += 1
                totals.frontier_gap_hist[gap] = totals.frontier_gap_hist.get(gap, 0) + 1

                rows.append(
                    {
                        "n": n,
                        "cpdag": cpdag.edge_string(),
                        "n_undirected": len(cpdag.undirected_edges),
                        "space_size": len(space),
                        "g0": g0.edge_string(),
                        "x": x,
                        "y": y,
                        "z_opt": "|".join(sorted(z)),
                        "r_val_opt": r_opt_set,
                        "best_r_val": best_r,
                        "best_z": "|".join(sorted(best_z)),
                        "frontier_gap": gap,
                        "n_valid_sets_examined": len(cands),
                        "method": "bfs_exact",
                    }
                )
    return totals, rows


def summarise(totals: CensusTotals) -> str:
    """Human-readable summary of a census, for the log and the report."""
    hist = Counter(totals.r_val_hist)
    total = sum(hist.values())
    lines = [
        f"n={totals.n_nodes}: {totals.n_instances} instances "
        f"from {totals.n_cpdags_used}/{totals.n_cpdags_total} CPDAGs",
        f"  r_val distribution (n={total}):",
    ]
    for r in sorted(hist):
        lines.append(f"    r={r}: {hist[r]:6d}  ({100 * hist[r] / max(total, 1):5.1f}%)")
    if totals.r_val_unreached:
        lines.append(f"    UNREACHED: {totals.r_val_unreached}")
    lines.append(
        f"  H2 frontier: {totals.n_frontier_strict}/{totals.n_frontier_checked} "
        f"instances have a strictly more robust valid set "
        f"({100 * totals.n_frontier_strict / max(totals.n_frontier_checked, 1):.1f}%)"
    )
    return "\n".join(lines)
