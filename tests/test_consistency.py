"""Consistency is not truth.

The tests that matter most here are the ones asserting the two predicates come
apart. If a change ever makes `is_consistent` and `is_true` agree everywhere,
the project's premise has been implemented away.
"""

from __future__ import annotations

import pytest

from bkrobust.graphs import consistency
from bkrobust.knowledge.base import BackgroundKnowledge


def test_empty_knowledge_is_consistent(mpdag, chain):
    """Asserting nothing is always consistent with anything."""
    cpdag = mpdag(chain, as_cpdag=True)
    with pytest.raises(NotImplementedError):
        consistency.is_consistent(cpdag, BackgroundKnowledge.empty())


def test_true_knowledge_is_consistent(mpdag, chain):
    """Knowledge read off the true DAG is always consistent with its CPDAG."""
    cpdag = mpdag(chain, as_cpdag=True)
    with pytest.raises(NotImplementedError):
        bk = BackgroundKnowledge.from_dag(mpdag(chain))
        consistency.is_consistent(cpdag, bk)


def test_reversed_chain_is_consistent_but_false(mpdag, chain):
    """Asserting C -> B -> A on the chain's CPDAG is consistent and false.

    The canonical example, and the one to reach for when explaining the project.
    A - B - C is equally compatible with A -> B -> C and C -> B -> A; the data
    cannot separate them. So an expert asserting the reverse passes every check
    a practitioner runs, and is wrong about both edges.
    """
    cpdag = mpdag(chain, as_cpdag=True)
    true_graph = mpdag(chain)
    with pytest.raises(NotImplementedError):
        bk = BackgroundKnowledge().with_edge("C", "B").with_edge("B", "A")
        assert consistency.is_consistent(cpdag, bk)
        assert not consistency.is_true(true_graph, bk)
        assert consistency.is_consistent_but_false(cpdag, true_graph, bk)


def test_contradicting_a_v_structure_is_inconsistent(mpdag, collider):
    """Asserting B -> A against the collider's oriented A -> B is inconsistent.

    The boundary of the study: this knowledge is caught by the existing check,
    so it is not the failure mode under investigation.
    """
    cpdag = mpdag(collider, as_cpdag=True)
    with pytest.raises(NotImplementedError):
        bk = BackgroundKnowledge().with_edge("B", "A")
        assert not consistency.is_consistent(cpdag, bk)


def test_check_consistency_reports_the_failing_constraint(mpdag, collider):
    """A failure report names which constraint broke and why."""
    cpdag = mpdag(collider, as_cpdag=True)
    with pytest.raises(NotImplementedError):
        bk = BackgroundKnowledge().with_edge("B", "A")
        report = consistency.check_consistency(cpdag, bk)
        assert report.failing_constraint == ("B", "A")


def test_false_constraints_counts_only_violations(mpdag, chain):
    """Only the constraints the true DAG violates are returned."""
    true_graph = mpdag(chain)
    with pytest.raises(NotImplementedError):
        bk = BackgroundKnowledge().with_edge("C", "B").with_edge("A", "B")
        assert consistency.false_constraints(true_graph, bk) == {("C", "B")}


def test_is_consistent_but_false_rejects_mismatched_true_graph(mpdag, chain, collider):
    """A true graph outside the CPDAG's equivalence class is a caller error."""
    cpdag = mpdag(chain, as_cpdag=True)
    wrong_truth = mpdag(collider)
    with pytest.raises((NotImplementedError, ValueError)):
        consistency.is_consistent_but_false(cpdag, wrong_truth, BackgroundKnowledge.empty())
