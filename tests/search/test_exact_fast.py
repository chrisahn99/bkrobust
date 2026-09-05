"""Regression tests for enumeration-free cover generation (Lemma O).

These pin the substitution that removes 70.4% of ``local_up``'s measured cost.
A wrong minimality test does not crash: it silently returns the wrong covers,
which shifts radii. So the tests compare against the frozen implementation
directly, including ordering, and against brute-force BFS.
"""

from __future__ import annotations

import itertools

import pytest

from bkrobust.core.oracle import extensions, is_valid
from bkrobust.core.spacelib import distances_from, radius
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag
from bkrobust.search.conjecture_study import all_cpdags
from bkrobust.search.exact import local_up_covers, radius_local_up
from bkrobust.search.exact_fast import local_up_covers_fast, radius_local_up_fast
from bkrobust.search.space_fixed import build_space_correct


@pytest.mark.parametrize("n", [3, 4])
def test_lemma_o_model_inclusion_is_reverse_edge_containment(n: int) -> None:
    """[G] ⊆ [H] iff dir(H) ⊆ dir(G), on every ordered pair of space elements."""
    for cpdag in list(all_cpdags(n))[:20]:
        if not cpdag.undirected_edges:
            continue
        space = build_space_correct(cpdag)
        reps = {g: frozenset(extensions(g)) for g in space.elements}
        for g, h in itertools.permutations(space.elements, 2):
            assert (reps[g] <= reps[h]) == (h.directed_edges <= g.directed_edges), (
                cpdag.edge_string(),
                g.edge_string(),
                h.edge_string(),
            )


@pytest.mark.parametrize("n", [3, 4])
def test_fast_covers_identical_to_frozen(n: int) -> None:
    """Same covers, same order, with no DAG-extension enumeration."""
    for cpdag in list(all_cpdags(n))[:25]:
        if not cpdag.undirected_edges:
            continue
        for g in build_space_correct(cpdag).elements:
            slow = [h.edge_string() for h in local_up_covers(cpdag, g)]
            fast = [h.edge_string() for h in local_up_covers_fast(cpdag, g)]
            assert slow == fast, (cpdag.edge_string(), g.edge_string())


def test_fast_covers_enumerate_nothing() -> None:
    """The counter that should stay at zero, stays at zero.

    This is the point of the change, so it is asserted rather than assumed.
    """
    from bkrobust.search.exact import SearchStats

    cpdag = next(c for c in all_cpdags(4) if len(c.undirected_edges) >= 3)
    stats = SearchStats()
    for g in build_space_correct(cpdag).elements:
        local_up_covers_fast(cpdag, g, stats)
    assert stats.extensions_enumerated == 0
    assert stats.closures > 0


def test_fast_radius_matches_frozen_and_bfs() -> None:
    """Radii are unchanged against both the frozen search and brute force."""
    checked = 0
    for cpdag in all_cpdags(4):
        if not cpdag.undirected_edges:
            continue
        space = build_space_correct(cpdag)
        for g0 in space.elements:
            dists = distances_from(space, g0)
            for x, y in itertools.permutations(sorted(cpdag.nodes), 2):
                o = optimal_adjustment_set_mpdag(g0, x, y)
                if o is None:
                    continue
                z = frozenset(o)
                if not z or not is_valid(z, g0, x, y):
                    continue
                fails = lambda g, _z=z, _x=x, _y=y: not is_valid(_z, g, _x, _y)  # noqa: E731
                r_bfs, _ = radius(space, dists, fails)
                assert radius_local_up(cpdag, g0, fails).radius == r_bfs
                assert radius_local_up_fast(cpdag, g0, fails).radius == r_bfs
                checked += 1
                if checked >= 60:
                    return
    assert checked > 0
