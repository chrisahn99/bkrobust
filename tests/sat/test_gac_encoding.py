"""Regression tests for the GAC failure mode of the SAT encoding.

``bkrobust.sat.failure.add_failure(..., criterion="gac")`` is a second failure
predicate on the same witness DAG, checked against
:func:`bkrobust.gac.mpdag_level.is_gac_valid_mpdag` (the polynomial-time
graphical decision procedure, itself validated against enumeration
elsewhere). These tests pin the two facts that matter most:

* the reviewer's counterexample, where GAC and back-door genuinely disagree
  (``tests/sat/test_encodings.py`` and ``tests/hybrid/test_hybrid.py`` only
  ever exercise back-door, so nothing else in the suite would catch a
  regression here);
* full agreement between the SAT E1 ladder (forced onto the SAT path, not the
  hybrid's search shortcut) and the brute-force GAC radius, over a slice of
  the exhaustive small-graph scope, for both the optimal adjustment set and
  arbitrary GAC-valid candidate sets.
"""

from __future__ import annotations

import itertools

from bkrobust.core.conventions import UNREACHED
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag
from bkrobust.hybrid import breakdown_radius
from bkrobust.sat.e1 import radius_e1
from bkrobust.search.conjecture_study import all_cpdags
from bkrobust.search.space_fixed import build_space_correct


def _brute_gac_radius(cpdag: MPDAG, g0: MPDAG, x: str, y: str, z: frozenset[str]) -> int:
    """The brute-force GAC radius.

    Min over ``[cp, *build_space_correct(cp).elements]`` of ``|K_G0 Δ K_G|``
    over GAC-invalid ``G`` -- the reference this module's SAT tests are
    measured against.
    """
    k0 = set(g0.directed_edges)
    best = None
    for g in (cpdag, *build_space_correct(cpdag).elements):
        if is_gac_valid_mpdag(g, x, y, z):
            continue
        d = len(k0 ^ set(g.directed_edges))
        if best is None or d < best:
            best = d
    return UNREACHED if best is None else best


def test_reviewer_counterexample_gac_radius_is_two_backdoor_is_one() -> None:
    """``W - X - Y``, ``G0: W -> X -> Y``, ``Z = {W}``.

    ``W`` is a descendant of ``X`` but not on (or below) the causal route to
    ``Y``, so back-door bars it (radius 1: flip ``W -> X``, i.e. the *only*
    other retraction, to fail) while GAC does not (radius 2: also break the
    ``X -> Y`` orientation before ``W`` starts blocking a back-door path).
    This is the minimal witness from the module docstrings of
    ``bkrobust.gac.dag_level`` and ``bkrobust.mpdag_criterion``.
    """
    cpdag = MPDAG(["W", "X", "Y"], [], [("W", "X"), ("X", "Y")])
    g0 = MPDAG(["W", "X", "Y"], [("W", "X"), ("X", "Y")], [])
    z = frozenset({"W"})

    assert is_gac_valid_mpdag(g0, "X", "Y", z)

    r_gac = radius_e1(cpdag, g0, "X", "Y", z, criterion="gac")
    r_bd = radius_e1(cpdag, g0, "X", "Y", z, criterion="backdoor")
    assert r_gac.radius == 2
    assert r_bd.radius == 1

    h_gac = breakdown_radius(cpdag, [("W", "X"), ("X", "Y")], "X", "Y", z)
    h_bd = breakdown_radius(cpdag, [("W", "X"), ("X", "Y")], "X", "Y", z, criterion="backdoor")
    assert h_gac.radius == 2
    assert h_gac.oracle == "gac_criterion"
    assert h_bd.radius == 1
    assert h_bd.oracle == "mpdag_criterion"


def test_sat_e1_gac_matches_brute_force_optimal_set_n3() -> None:
    """E1's GAC radius equals the brute-force GAC radius.

    Forced onto the SAT path, over every CPDAG on 3 nodes with an undirected
    edge, every element, every ordered ``(x, y)``, for ``Z = O(G0)``.
    """
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
                got = radius_e1(cpdag, g0, x, y, z, criterion="gac", max_k=None)
                assert got.radius == want, (cpdag.edge_string(), g0.edge_string(), x, y, sorted(z))
                checked += 1
    assert checked > 0


def test_sat_e1_gac_matches_brute_force_arbitrary_sets_n3() -> None:
    """Same as above but for arbitrary GAC-valid ``Z``.

    Size <= 2, not just the optimal set -- the scope where GAC and back-door
    actually diverge.
    """
    checked = 0
    for cpdag in all_cpdags(3):
        if not cpdag.undirected_edges:
            continue
        for g0 in build_space_correct(cpdag).elements:
            for x, y in itertools.permutations(sorted(cpdag.nodes), 2):
                others = sorted(set(cpdag.nodes) - {x, y})
                for size in (0, 1, 2):
                    for combo in itertools.combinations(others, size):
                        z = frozenset(combo)
                        if not is_gac_valid_mpdag(g0, x, y, z):
                            continue
                        want = _brute_gac_radius(cpdag, g0, x, y, z)
                        got = radius_e1(cpdag, g0, x, y, z, criterion="gac", max_k=None)
                        assert got.radius == want, (
                            cpdag.edge_string(),
                            g0.edge_string(),
                            x,
                            y,
                            sorted(z),
                        )
                        checked += 1
    assert checked > 0
