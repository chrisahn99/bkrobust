"""Regression tests for :mod:`bkrobust.hybrid`'s GAC default.

``breakdown_radius`` switched from back-door to the generalised adjustment
criterion (GAC) by default. ``tests/hybrid/test_hybrid.py`` covers the
back-door predicate (pinned there with ``criterion="backdoor"``); this module
covers the new default end to end, including the SAT-ladder leg (forced with
``search_budget=0``, since on this small scope the bounded search alone
almost always answers first and the ladder would otherwise go untested).
"""

from __future__ import annotations

import itertools

from bkrobust.core.conventions import UNREACHED
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag
from bkrobust.hybrid import breakdown_radius
from bkrobust.search.conjecture_study import all_cpdags
from bkrobust.search.space_fixed import build_space_correct


def _brute_gac_radius(cpdag, g0, x, y, z) -> int:
    k0 = set(g0.directed_edges)
    best = None
    for g in (cpdag, *build_space_correct(cpdag).elements):
        if is_gac_valid_mpdag(g, x, y, z):
            continue
        d = len(k0 ^ set(g.directed_edges))
        if best is None or d < best:
            best = d
    return UNREACHED if best is None else best


def test_hybrid_default_is_gac() -> None:
    """No ``criterion`` argument means GAC, labelled ``gac_criterion``."""
    from bkrobust.demo.graph import MPDAG

    cpdag = MPDAG(["W", "X", "Y"], [], [("W", "X"), ("X", "Y")])
    g0 = MPDAG(["W", "X", "Y"], [("W", "X"), ("X", "Y")], [])
    z = frozenset({"W"})
    got = breakdown_radius(cpdag, [("W", "X"), ("X", "Y")], "X", "Y", z, g0=g0)
    assert got.oracle == "gac_criterion"
    assert got.radius == 2


def test_hybrid_gac_matches_brute_force_n3() -> None:
    """Default-budget hybrid agrees with brute-force GAC radius, n=3 scope."""
    checked = 0
    for cpdag in all_cpdags(3):
        if not cpdag.undirected_edges:
            continue
        for g0 in build_space_correct(cpdag).elements:
            for x, y in itertools.permutations(sorted(cpdag.nodes), 2):
                o = optimal_adjustment_set_mpdag(g0, x, y)
                if o is None:
                    continue
                z = frozenset(o)
                if not z or not is_gac_valid_mpdag(g0, x, y, z):
                    continue
                want = _brute_gac_radius(cpdag, g0, x, y, z)
                got = breakdown_radius(cpdag, None, x, y, z, g0=g0)
                assert got.radius == want
                assert got.oracle == "gac_criterion"
                checked += 1
    assert checked > 0


def test_hybrid_gac_ladder_leg_matches_brute_force_n3() -> None:
    """Same, but with ``search_budget=0``.

    That forces the E1 ladder to be exercised for every non-degenerate,
    non-``top_state`` instance.
    """
    checked = 0
    ladder_hits = 0
    for cpdag in all_cpdags(3):
        if not cpdag.undirected_edges:
            continue
        for g0 in build_space_correct(cpdag).elements:
            for x, y in itertools.permutations(sorted(cpdag.nodes), 2):
                o = optimal_adjustment_set_mpdag(g0, x, y)
                if o is None:
                    continue
                z = frozenset(o)
                if not z or not is_gac_valid_mpdag(g0, x, y, z):
                    continue
                want = _brute_gac_radius(cpdag, g0, x, y, z)
                got = breakdown_radius(cpdag, None, x, y, z, g0=g0, search_budget=0)
                assert got.radius == want
                if got.method == "e1_ladder":
                    ladder_hits += 1
                checked += 1
    assert checked > 0
    assert ladder_hits > 0
