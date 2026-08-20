"""Adjustment sets, and O* against worked examples.

The `nontrivial_optimal` fixture carries the case that matters: O* and the
canonical back-door set differ there, so a test that only ever sees graphs where
they coincide would pass against an implementation that returns the back-door
set and calls it optimal.
"""

from __future__ import annotations

import pytest

from bkrobust.graphs import adjustment


def test_empty_set_valid_for_chain(mpdag, chain):
    """A -> B -> C is unconfounded, so the empty set is valid for A on C."""
    graph = mpdag(chain)
    with pytest.raises(NotImplementedError):
        assert adjustment.is_valid_adjustment_set(graph, "A", "C", set())


def test_mediator_is_forbidden(mpdag, chain):
    """Adjusting for B blocks the very path being estimated, so it is invalid."""
    graph = mpdag(chain)
    with pytest.raises(NotImplementedError):
        assert not adjustment.is_valid_adjustment_set(graph, "A", "C", {"B"})


def test_common_cause_must_be_adjusted(mpdag, fork):
    """On the fork, the empty set is invalid and {B} is valid."""
    graph = mpdag(fork)
    with pytest.raises(NotImplementedError):
        assert not adjustment.is_valid_adjustment_set(graph, "A", "C", set())
        assert adjustment.is_valid_adjustment_set(graph, "A", "C", {"B"})


def test_optimal_set_on_nontrivial_example(mpdag, nontrivial_optimal):
    """O* is {Z1, Z2} where the canonical back-door set is {Z1}.

    Z2 touches only the outcome. It is not needed for identification and it is
    needed for efficiency, which is the distinction the whole optimality theory
    turns on.
    """
    graph = mpdag(nontrivial_optimal)
    with pytest.raises(NotImplementedError):
        result = adjustment.optimal_adjustment_set(graph, "X", "Y")
        assert result == nontrivial_optimal["optimal_adjustment_set"]


def test_canonical_set_differs_from_optimal(mpdag, nontrivial_optimal):
    """Both are valid; they are not the same set."""
    graph = mpdag(nontrivial_optimal)
    with pytest.raises(NotImplementedError):
        canonical = adjustment.canonical_adjustment_set(graph, "X", "Y")
        optimal = adjustment.optimal_adjustment_set(graph, "X", "Y")
        assert canonical == nontrivial_optimal["canonical_adjustment_set"]
        assert canonical != optimal


def test_forbidden_set_contains_mediators_and_treatment(mpdag, nontrivial_optimal):
    """The forbidden set is the possible descendants of the causal nodes, plus X."""
    graph = mpdag(nontrivial_optimal)
    with pytest.raises(NotImplementedError):
        assert adjustment.forbidden_set(graph, "X", "Y") == nontrivial_optimal["forbidden_set"]


def test_optimal_set_is_valid(mpdag, canonical_graphs):
    """O* must be valid on every graph where it is defined.

    An optimal set that is not valid is not an optimisation, it is a bug -- and
    the failure would look like a small variance improvement, which is exactly
    how it would get missed.
    """
    for spec in canonical_graphs.values():
        graph = mpdag(spec)
        with pytest.raises(NotImplementedError):
            z = adjustment.optimal_adjustment_set(graph, spec["treatment"], spec["outcome"])
            assert adjustment.is_valid_adjustment_set(graph, spec["treatment"], spec["outcome"], z)


def test_all_valid_sets_includes_optimal(mpdag, nontrivial_optimal):
    """Enumeration must contain O*, or one of the two is wrong."""
    graph = mpdag(nontrivial_optimal)
    with pytest.raises(NotImplementedError):
        sets = adjustment.all_valid_adjustment_sets(graph, "X", "Y")
        assert nontrivial_optimal["optimal_adjustment_set"] in sets


@pytest.mark.slow
def test_optimal_matches_brute_force_minimum(mpdag, nontrivial_optimal):
    """O* must equal the enumerated variance minimiser under a sampled SCM.

    The strongest available check on the O* construction: it compares a graphical
    shortcut against the quantity that shortcut is supposed to compute.
    """
    from bkrobust.graphs import optimality

    graph = mpdag(nontrivial_optimal)
    with pytest.raises(NotImplementedError):
        best, _ = optimality.minimum_variance_adjustment_set(graph, "X", "Y", scm=None)
        assert best == adjustment.optimal_adjustment_set(graph, "X", "Y")
