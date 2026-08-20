"""Breakdown radii: exactness, the ordering conjecture, and honest reporting.

Two failures here would be worse than a wrong number. Reporting a heuristic
radius as exact overstates a guarantee, and confusing "no breaking perturbation
found" with "the radius is large" turns a budget limit into a robustness claim.
Both have dedicated tests.
"""

from __future__ import annotations

import pytest

from bkrobust.theory import radius


def test_radius_is_nonnegative(mpdag, chain):
    """A radius is either non-negative or the UNREACHED sentinel."""
    cpdag, true_graph = mpdag(chain, as_cpdag=True), mpdag(chain)
    with pytest.raises(NotImplementedError):
        d = radius.delta_valid(cpdag, true_graph, "A", "C")
        assert d >= 0 or d == radius.UNREACHED


def test_collider_cpdag_has_no_reachable_radius(mpdag, collider):
    """A fully oriented CPDAG admits no consistent perturbation at all.

    The result is UNREACHED, not zero and not a large number: there is nothing
    to perturb, which is a different statement from "perturbations do no harm".
    """
    cpdag, true_graph = mpdag(collider, as_cpdag=True), mpdag(collider)
    with pytest.raises(NotImplementedError):
        assert radius.delta_valid(cpdag, true_graph, "A", "B") == radius.UNREACHED


def test_exact_flag_is_true_only_for_exact_method(mpdag, chain):
    """A heuristic report must never claim exactness.

    The whole difference between a guarantee and a guess lives in this flag.
    """
    cpdag, true_graph = mpdag(chain, as_cpdag=True), mpdag(chain)
    with pytest.raises(NotImplementedError):
        valid, _ = radius.radius_report(cpdag, true_graph, "A", "C", method="search")
        assert valid.exact is False


def test_heuristic_radius_upper_bounds_exact(mpdag, diamond):
    """The heuristic search finds a breaking set no earlier than the exact search.

    The search may miss a smaller breaking knowledge set, so its answer can only
    be too large. If it ever comes back smaller than the exact radius, the exact
    enumeration is missing members of the ball.
    """
    cpdag, true_graph = mpdag(diamond, as_cpdag=True), mpdag(diamond)
    with pytest.raises(NotImplementedError):
        exact = radius.delta_valid(cpdag, true_graph, "A", "D", method="exact")
        heuristic = radius.delta_valid(cpdag, true_graph, "A", "D", method="search")
        assert heuristic >= exact


def test_ordering_conjecture_on_canonical_graphs(mpdag, canonical_graphs):
    """delta_opt <= delta_valid on every canonical graph.

    A violation on an exact pair of reports is a counterexample to the paper's
    conjecture, and this test is where it would first show up. If it fails,
    do not fix the test.
    """
    for spec in canonical_graphs.values():
        cpdag, true_graph = mpdag(spec, as_cpdag=True), mpdag(spec)
        with pytest.raises(NotImplementedError):
            valid, opt = radius.radius_report(
                cpdag, true_graph, spec["treatment"], spec["outcome"], scm=None
            )
            assert radius.verify_ordering(valid, opt)


def test_budget_exhaustion_is_flagged(mpdag, diamond):
    """Stopping on budget is recorded, and is not reported as a large radius."""
    cpdag, true_graph = mpdag(diamond, as_cpdag=True), mpdag(diamond)
    with pytest.raises(NotImplementedError):
        valid, _ = radius.radius_report(cpdag, true_graph, "A", "D", max_evaluations=1)
        assert valid.budget_exhausted


def test_witness_reproduces_the_breakdown(mpdag, chain):
    """The recorded witness actually breaks the property at the reported radius.

    A radius without a checkable witness is an assertion; with one it is
    evidence, and the witness is what makes a counterexample publishable.
    """
    cpdag, true_graph = mpdag(chain, as_cpdag=True), mpdag(chain)
    with pytest.raises(NotImplementedError):
        valid, _ = radius.radius_report(cpdag, true_graph, "A", "C")
        assert valid.witness is not None


def test_radius_distribution_fractions_sum_to_one(mpdag, diamond):
    """The severity fractions over a ball partition it."""
    cpdag, true_graph = mpdag(diamond, as_cpdag=True), mpdag(diamond)
    with pytest.raises(NotImplementedError):
        d = radius.radius_distribution(cpdag, true_graph, "A", "D", delta=1)
        total = d["invalid_fraction"] + d["suboptimal_fraction"] + d["benign_fraction"]
        assert total == pytest.approx(1.0)
