"""Tests for :mod:`bkrobust.hybrid`, the single entry point.

This is the artefact the project hands to someone else, so the tests cover not
only that the radius is right but that the *labelling* is right: a caller who
trusts ``method``/``oracle``/``assumes`` and gets them wrong is worse off than
one who had to choose a method themselves.
"""

from __future__ import annotations

import itertools

from bkrobust.core.conventions import UNREACHED
from bkrobust.core.oracle import is_valid
from bkrobust.core.spacelib import distances_from, radius
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.scenario import OUTCOME, TREATMENT, scenarios, true_dag
from bkrobust.hybrid import breakdown_radius, describe
from bkrobust.search.conjecture_study import all_cpdags
from bkrobust.search.space_fixed import build_space_correct


def test_worked_example_radii() -> None:
    """The published 3 / 3 / 2 on the session-1 worked example."""
    cpdag = dag_to_cpdag(true_dag())
    want = {"A": 3, "B": 3, "C": 2}
    for name, spec in scenarios().items():
        g0 = apply_orientations(cpdag, spec["knowledge"])
        z = frozenset(optimal_adjustment_set_mpdag(g0, TREATMENT, OUTCOME))
        got = breakdown_radius(cpdag, spec["knowledge"], TREATMENT, OUTCOME, z)
        assert got.radius == want[name], name
        assert got.exact


def test_agrees_with_brute_force_under_both_oracles() -> None:
    """Radius unchanged whether validity comes from the criterion or enumeration."""
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
                want, _ = radius(space, dists, fails)
                for flag in (True, False):
                    got = breakdown_radius(cpdag, None, x, y, z, g0=g0, use_criterion=flag)
                    assert got.radius == want
                    assert got.oracle == ("mpdag_criterion" if flag else "enumeration")
                checked += 1
                if checked >= 40:
                    return
    assert checked > 0


def test_result_carries_its_assumption() -> None:
    """The Conjecture 2 dependency travels with every answer."""
    cpdag = dag_to_cpdag(true_dag())
    spec = scenarios()["A"]
    g0 = apply_orientations(cpdag, spec["knowledge"])
    z = frozenset(optimal_adjustment_set_mpdag(g0, TREATMENT, OUTCOME))
    got = breakdown_radius(cpdag, None, TREATMENT, OUTCOME, z, g0=g0)
    assert "Conjecture 2" in got.assumes
    assert "Anti-Exchange" in got.assumes


def test_budget_does_not_change_the_answer() -> None:
    """Dispatch is a performance choice, never a semantic one.

    With a budget of 0 the search cannot find anything and every instance falls
    through to the ladder; the radii must be identical either way.
    """
    cpdag = dag_to_cpdag(true_dag())
    for spec in scenarios().values():
        g0 = apply_orientations(cpdag, spec["knowledge"])
        z = frozenset(optimal_adjustment_set_mpdag(g0, TREATMENT, OUTCOME))
        wide = breakdown_radius(cpdag, None, TREATMENT, OUTCOME, z, g0=g0, search_budget=9)
        none = breakdown_radius(cpdag, None, TREATMENT, OUTCOME, z, g0=g0, search_budget=0)
        assert wide.radius == none.radius
        assert none.method == "e1_ladder"
        assert wide.method == "local_up_fast"


def test_describe_never_reports_a_sentinel_as_a_number() -> None:
    """UNREACHED must read as 'no failure exists', not as a radius of -1."""
    from bkrobust.hybrid import HybridResult

    text = describe(
        HybridResult(radius=UNREACHED, method="local_up_fast", oracle="mpdag_criterion")
    )
    assert "-1" not in text
    assert "robust" in text
    inexact = describe(
        HybridResult(radius=UNREACHED, method="e1_ladder", oracle="mpdag_criterion", exact=False)
    )
    assert "not" in inexact.lower()


def test_inconsistent_knowledge_is_rejected() -> None:
    """Knowledge that admits no G0 raises rather than silently returning a radius."""
    import pytest

    cpdag = dag_to_cpdag(true_dag())
    und = sorted(cpdag.undirected_edges)
    bad = [(a, b) for a, b in und] + [(b, a) for a, b in und]
    with pytest.raises(ValueError, match="inconsistent"):
        breakdown_radius(cpdag, bad, TREATMENT, OUTCOME, frozenset())
