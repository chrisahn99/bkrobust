"""The efficiency-robustness frontier: does any valid adjustment set beat O*?

This module exists because the first attempt at the H2 statistic was wrong, in a
way worth recording. That version compared ``r_val`` across adjustment sets but
**excluded every comparison in which either radius was UNREACHED**, on the
grounds that UNREACHED is not a number. It is not a number -- but it is not
missing data either. A set that never fails anywhere in the fully enumerated
space is *strictly more robust* than one that fails at distance 1, and that is
the largest frontier gap there is. Discarding those comparisons discarded
precisely the cases where the frontier is non-flat, and produced a confident
"0 out of 134,140" that was an artefact of the analysis.

The correct ordering is

    UNREACHED  >  ...  >  3  >  2  >  1

i.e. UNREACHED is treated as "more robust than any finite radius", never as -1
and never as absent.

Two further choices, both material and both stated rather than buried:

* The **empty adjustment set is excluded** from the comparison. An empty set
  cannot acquire a descendant of the treatment, so it is trivially robust for a
  reason unrelated to the phenomenon, and including it inflates the non-flat
  rate (23.05% with it, 10.96% without, at n=4). Instances where the empty set
  is valid are in any case gated out of the ensembles as unconfounded.
* **No cap** on the number of valid adjustment sets examined. The earlier census
  capped at 10 sorted by size, which can hide a more robust larger set.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from typing import Any

from bkrobust.core.conventions import UNREACHED
from bkrobust.core.oracle import is_valid
from bkrobust.core.spacelib import build_space, distances_from, radius
from bkrobust.demo.evaluate import all_valid_adjustment_sets_mpdag, optimal_adjustment_set_mpdag
from bkrobust.search.conjecture_study import all_cpdags


def more_robust(a: int, b: int) -> bool:
    """Whether radius ``a`` is strictly more robust than radius ``b``.

    UNREACHED beats every finite radius; otherwise larger is more robust.
    """
    if a == UNREACHED:
        return b != UNREACHED
    if b == UNREACHED:
        return False
    return a > b


@dataclass
class FrontierTotals:
    """Aggregates for one frontier sweep."""

    n_nodes: int = 0
    n_instances: int = 0
    n_nonflat: int = 0
    n_optimal_is_strictly_beaten: int = 0
    n_optimal_unreached: int = 0
    gap_kinds: dict[str, int] = field(default_factory=dict)
    examples: list[dict[str, Any]] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        """JSON-friendly view."""
        return {
            "n_nodes": self.n_nodes,
            "n_instances": self.n_instances,
            "n_nonflat": self.n_nonflat,
            "frac_nonflat": (self.n_nonflat / self.n_instances) if self.n_instances else None,
            "n_optimal_is_strictly_beaten": self.n_optimal_is_strictly_beaten,
            "frac_optimal_beaten": (self.n_optimal_is_strictly_beaten / self.n_instances)
            if self.n_instances
            else None,
            "n_optimal_unreached": self.n_optimal_unreached,
            "gap_kinds": dict(sorted(self.gap_kinds.items())),
            "examples": self.examples[:5],
        }


def frontier_sweep(
    n: int,
    *,
    max_undirected: int = 5,
    include_empty_z: bool = False,
) -> FrontierTotals:
    """Exhaustively compare the radii of every valid adjustment set, per instance.

    Args:
        n: Node count.
        max_undirected: Skip CPDAGs above this undirected-edge count.
        include_empty_z: Include the empty adjustment set. Default ``False``;
            see the module docstring for why it distorts the statistic.

    Returns:
        A :class:`FrontierTotals`.
    """
    totals = FrontierTotals(n_nodes=n)
    for cpdag in all_cpdags(n):
        if not cpdag.undirected_edges or len(cpdag.undirected_edges) > max_undirected:
            continue
        space = build_space(cpdag)
        for g0 in space.elements:
            dists = distances_from(space, g0)
            for x, y in itertools.permutations(list(cpdag.nodes), 2):
                cands = [frozenset(z) for z in all_valid_adjustment_sets_mpdag(g0, x, y)]
                zs = [z for z in cands if (include_empty_z or z) and is_valid(z, g0, x, y)]
                if len(zs) < 2:
                    continue
                radii = []
                for z in zs:
                    r, _ = radius(
                        space,
                        dists,
                        lambda g, _z=z, _x=x, _y=y: not is_valid(_z, g, _x, _y),
                    )
                    radii.append(r)
                totals.n_instances += 1
                if len(set(radii)) > 1:
                    totals.n_nonflat += 1

                opt = optimal_adjustment_set_mpdag(g0, x, y)
                if opt is None:
                    continue
                zopt = frozenset(opt)
                if zopt not in zs:
                    continue
                r_opt = radii[zs.index(zopt)]
                if r_opt == UNREACHED:
                    totals.n_optimal_unreached += 1
                beaten = any(more_robust(r, r_opt) for r in radii)
                if beaten:
                    totals.n_optimal_is_strictly_beaten += 1
                    kind = (
                        "unreached_beats_finite"
                        if r_opt != UNREACHED and any(r == UNREACHED for r in radii)
                        else "finite_gap"
                    )
                    totals.gap_kinds[kind] = totals.gap_kinds.get(kind, 0) + 1
                    if len(totals.examples) < 5:
                        totals.examples.append(
                            {
                                "cpdag": cpdag.edge_string(),
                                "g0": g0.edge_string(),
                                "x": x,
                                "y": y,
                                "z_optimal": sorted(zopt),
                                "r_optimal": r_opt,
                                "sets": [
                                    (sorted(z), r)
                                    for z, r in zip(zs, radii)  # noqa: B905
                                ][:6],
                            }
                        )
    return totals
