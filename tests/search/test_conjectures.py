"""Tests for the two structural conjectures and the machinery that checks them."""

from __future__ import annotations

import pytest

from bkrobust.core.oracle import is_valid
from bkrobust.core.spacelib import build_space, distances_from
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.scenario import OUTCOME, TREATMENT, scenarios, true_dag
from bkrobust.search.conjectures import (
    check_conjecture1,
    check_conjecture2,
    check_no_empty_extensions,
    up_distances,
    up_neighbours,
)


@pytest.fixture(scope="module")
def worked() -> tuple:
    cpdag = dag_to_cpdag(true_dag())
    return cpdag, build_space(cpdag)


def _setup(worked: tuple, label: str) -> tuple:
    cpdag, space = worked
    g0 = apply_orientations(cpdag, scenarios()[label]["knowledge"])
    z = frozenset(optimal_adjustment_set_mpdag(g0, TREATMENT, OUTCOME) or ())
    return space, g0, z


# --- Conjecture 1 -----------------------------------------------------------


@pytest.mark.parametrize("label", ["A", "B", "C"])
def test_conjecture1_holds_on_worked_example(worked, label):
    """Upward-closure of failure. This is a theorem; a violation means a bug."""
    space, _g0, z = _setup(worked, label)
    res = check_conjecture1(space, lambda g: not is_valid(z, g, TREATMENT, OUTCOME))
    assert res.holds, res.violations[:3]
    assert res.n_pairs_checked > 0
    assert res.n_failing > 0, "a vacuous check would pass trivially"


def test_space_has_no_empty_extension_elements(worked):
    """The precondition that makes Conjecture 1 hold under this is_valid convention."""
    _cpdag, space = worked
    assert check_no_empty_extensions(space)


def test_empty_extension_graph_is_the_documented_caveat():
    """A graph representing no DAG is reported invalid, which is why the caveat exists.

    A chordless 4-cycle admits no consistent extension. Under the for-all
    convention it "fails" vacuously, so upward-closure would break for it -- but
    such graphs never enter the space, which is what makes the theorem apply.
    """
    cycle = MPDAG(
        ["A", "B", "C", "D"],
        undirected=[("A", "B"), ("B", "C"), ("C", "D"), ("A", "D")],
    )
    assert enumerate_dag_extensions(cycle) == []
    assert not is_valid(frozenset(), cycle, "A", "C")


# --- Conjecture 2 -----------------------------------------------------------


@pytest.mark.parametrize("label", ["A", "B", "C"])
def test_conjecture2_holds_on_worked_example(worked, label):
    """The nearest failure is reachable by retractions alone."""
    space, g0, z = _setup(worked, label)
    dists = distances_from(space, g0)
    res = check_conjecture2(space, g0, lambda g: not is_valid(z, g, TREATMENT, OUTCOME), dists)
    assert res.holds
    assert res.r_full >= 1


def test_up_distance_is_never_below_full_distance(worked):
    """Restricting to upward moves can only lengthen a path, never shorten it."""
    space, g0, _z = _setup(worked, "B")
    full = distances_from(space, g0)
    ups = up_distances(space, g0)
    for g, d in ups.items():
        assert d >= full[g]


def test_up_neighbours_are_strictly_above(worked):
    """Every upward neighbour represents strictly more DAGs -- a retracted claim."""
    space, g0, _z = _setup(worked, "B")
    for h in up_neighbours(space, g0):
        assert space.reps[g0] < space.reps[h]


def test_up_set_of_cpdag_is_only_itself(worked):
    """The CPDAG is the maximum: nothing can be retracted from it."""
    cpdag, space = worked
    assert up_neighbours(space, cpdag) == []
    assert up_distances(space, cpdag) == {cpdag: 0}


# --- study machinery --------------------------------------------------------


def test_dag_and_cpdag_counts_match_published_values():
    """Independent validation: labelled DAG counts and MEC counts are known.

    OEIS A003024 gives 25 and 543 labelled DAGs on 3 and 4 nodes; the numbers of
    Markov equivalence classes are 11 and 185. Matching both is a strong check on
    the enumerator and on dag_to_cpdag.
    """
    from bkrobust.search.conjecture_study import all_cpdags, all_dags

    assert len(all_dags(3)) == 25
    assert len(all_dags(4)) == 543
    assert len(all_cpdags(3)) == 11
    assert len(all_cpdags(4)) == 185
