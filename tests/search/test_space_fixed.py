"""Task 0: the corrected space, and why the old chordality filter was wrong.

These tests pin the resolution. If they fail, the space definition has moved and
every radius in the session moves with it.
"""

from __future__ import annotations

import itertools

import pytest

from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import (
    apply_orientations,
    enumerate_dag_extensions,
    is_meek_closed,
    is_valid_mpdag,
)
from bkrobust.demo.scenario import true_dag
from bkrobust.demo.space import enumerate_space
from bkrobust.search.conjecture_study import all_cpdags
from bkrobust.search.space_fixed import (
    build_space_correct,
    enumerate_space_correct,
    is_knowledge_state,
    is_maximally_oriented,
    knowledge_of,
)

#: The minimal witness found in Task 0: K4 minus one edge, plus ONE assertion.
WITNESS_CPDAG = MPDAG(
    ["V0", "V1", "V2", "V3"],
    undirected=[("V0", "V1"), ("V0", "V2"), ("V0", "V3"), ("V1", "V2"), ("V1", "V3")],
)
WITNESS_KNOWLEDGE = [("V0", "V1")]


def test_minimal_witness_is_excluded_by_the_old_filter():
    """A single assertion produces a state the old enumerate_space drops."""
    g = apply_orientations(WITNESS_CPDAG, WITNESS_KNOWLEDGE)
    assert g is not None
    assert g.edge_string() not in {h.edge_string() for h in enumerate_space(WITNESS_CPDAG)}


def test_the_excluded_witness_is_a_legitimate_mpdag():
    """It is Meek-closed, represents DAGs, and is MAXIMALLY ORIENTED.

    Maximal orientation is the defining property of an MPDAG: every undirected
    edge must be genuinely undecided across the represented DAGs. The excluded
    graph satisfies it, so excluding it was wrong.
    """
    g = apply_orientations(WITNESS_CPDAG, WITNESS_KNOWLEDGE)
    assert is_meek_closed(g)
    assert len(enumerate_dag_extensions(g)) > 0
    assert is_maximally_oriented(g) is True
    assert not is_valid_mpdag(g), "the old predicate is the thing that rejects it"


def test_corrected_space_contains_the_witness():
    g = apply_orientations(WITNESS_CPDAG, WITNESS_KNOWLEDGE)
    assert g.edge_string() in {h.edge_string() for h in enumerate_space_correct(WITNESS_CPDAG)}


def test_fixpoint_predicate_characterises_reachability():
    """Every element of the corrected space satisfies the membership predicate."""
    n_checked = 0
    for n in (3, 4):
        for cpdag in all_cpdags(n):
            if not cpdag.undirected_edges:
                continue
            for g in enumerate_space_correct(cpdag):
                assert is_knowledge_state(cpdag, g)
                n_checked += 1
    assert n_checked > 1000, "the sweep must actually exercise something"


def test_corrected_space_is_a_superset_of_the_old_one():
    """The correction only ever ADDS elements; it never removes a valid one.

    Direction matters: if the corrected space dropped an element the old one had,
    the fix would itself be losing knowledge states.
    """
    for n in (3, 4):
        for cpdag in all_cpdags(n):
            if not cpdag.undirected_edges:
                continue
            old = {g.edge_string() for g in enumerate_space(cpdag)}
            new = {g.edge_string() for g in enumerate_space_correct(cpdag)}
            assert old <= new, cpdag.edge_string()


def test_every_corrected_element_is_reachable_by_explicit_knowledge():
    """Reachability is constructive: each element has a K that produces it."""
    cpdag = WITNESS_CPDAG
    for g in enumerate_space_correct(cpdag):
        k = knowledge_of(cpdag, g)
        assert apply_orientations(cpdag, k) == g


def test_worked_example_space_is_unchanged():
    """The prior report's example has no excluded states, so its radii stand."""
    cpdag = dag_to_cpdag(true_dag())
    assert len(enumerate_space_correct(cpdag)) == len(enumerate_space(cpdag)) == 48


def test_corrected_space_matches_brute_force_reachability():
    """Independent re-derivation: 3^k closures, deduplicated, must match."""
    for cpdag in all_cpdags(4):
        if not cpdag.undirected_edges or len(cpdag.undirected_edges) > 4:
            continue
        und = sorted(cpdag.undirected_edges)
        brute = set()
        for assign in itertools.product((0, 1, 2), repeat=len(und)):
            ors = []
            for (a, b), v in zip(und, assign):  # noqa: B905
                if v == 1:
                    ors.append((a, b))
                elif v == 2:
                    ors.append((b, a))
            h = apply_orientations(cpdag, ors)
            if h is not None:
                brute.add(h.edge_string())
        assert brute == {g.edge_string() for g in enumerate_space_correct(cpdag)}


def test_build_space_correct_has_consistent_order():
    """The covering relation is still computed from model inclusion."""
    sp = build_space_correct(WITNESS_CPDAG)
    assert len(sp) == len(enumerate_space_correct(WITNESS_CPDAG))
    for lo, hi in sp.covers:
        assert sp.reps[lo] < sp.reps[hi]


@pytest.mark.parametrize("k", [5, 6])
def test_exclusion_rate_is_nonzero_at_higher_density(k):
    """The defect is not hypothetical: at k>=5 it affects real CPDAGs."""
    pool = [c for c in all_cpdags(5) if len(c.undirected_edges) == k][:6]
    affected = 0
    for cpdag in pool:
        old = {g.edge_string() for g in enumerate_space(cpdag)}
        new = {g.edge_string() for g in enumerate_space_correct(cpdag)}
        if new - old:
            affected += 1
    assert affected > 0, f"expected the defect to bite at k={k}"
