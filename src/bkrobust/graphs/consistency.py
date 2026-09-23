"""Consistency of background knowledge with a CPDAG.

This module is deliberately narrow, and its narrowness is the point of the
paper. *Consistency* asks one question: does there exist a DAG in the Markov
equivalence class described by the CPDAG that satisfies every constraint in the
knowledge set? Equivalently, does Meek's Algorithm 1 terminate without FAIL?

That is the entire check background knowledge receives in current practice.
It is a statement about the knowledge and the *observed* CPDAG. It is not a
statement about the true DAG, and knowledge can be consistent with the CPDAG
while contradicting the true DAG on every constraint it asserts.

:func:`is_consistent` is what a practitioner runs. :func:`is_true` is what they
cannot run, because it needs the true DAG; it exists here because the
simulation studies do have the true DAG and the whole design rests on holding
these two predicates apart.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from bkrobust.graphs.mpdag import MPDAG, Node
    from bkrobust.knowledge.base import BackgroundKnowledge


@dataclass(frozen=True)
class ConsistencyReport:
    """Outcome of a consistency check, with enough detail to act on a failure.

    Attributes:
        consistent: Whether Meek's Algorithm 1 terminated without FAIL.
        failing_constraint: The first constraint that could not be imposed, or
            ``None`` when ``consistent``.
        reason: Machine-readable failure category -- one of ``"cycle"``,
            ``"conflicting_orientation"``, ``"contradicts_cpdag"``,
            ``"unknown_node"``, or ``None`` when ``consistent``.
        n_constraints_imposed: How many constraints were imposed before the
            failure, so callers can report how deep the contradiction lay.
    """

    consistent: bool
    failing_constraint: tuple[Node, Node] | None
    reason: str | None
    n_constraints_imposed: int


def is_consistent(cpdag: MPDAG, bk: BackgroundKnowledge) -> bool:
    """Whether ``bk`` is consistent with ``cpdag`` in Perkovic's sense.

    True exactly when at least one DAG in the equivalence class represented by
    ``cpdag`` satisfies every constraint of ``bk``; operationally, when
    :func:`bkrobust.graphs.meek.apply_background_knowledge` does not FAIL.

    Args:
        cpdag: The CPDAG output by discovery.
        bk: The knowledge to check.

    Returns:
        ``True`` if consistent.

    Note:
        Says nothing about truth. See :func:`is_true`.
    """
    raise NotImplementedError


def check_consistency(cpdag: MPDAG, bk: BackgroundKnowledge) -> ConsistencyReport:
    """Run the consistency check and report why it failed if it did.

    Args:
        cpdag: The CPDAG output by discovery.
        bk: The knowledge to check.

    Returns:
        A :class:`ConsistencyReport`.
    """
    raise NotImplementedError


def is_true(true_graph: MPDAG, bk: BackgroundKnowledge) -> bool:
    """Whether every constraint in ``bk`` holds in ``true_graph``.

    The oracle predicate. Only available in simulation, where the true DAG is
    known; it is precisely the check a practitioner cannot perform, which is
    why consistent-but-false knowledge survives in the wild.

    Args:
        true_graph: The ground-truth DAG.
        bk: The knowledge to check.

    Returns:
        ``True`` if ``bk`` is true of ``true_graph``.
    """
    raise NotImplementedError


def false_constraints(
    true_graph: MPDAG,
    bk: BackgroundKnowledge,
) -> set[tuple[Node, Node]]:
    """Return the constraints of ``bk`` that ``true_graph`` violates.

    The size of this set is the natural "how wrong is this knowledge" counter,
    distinct from the perturbation radius delta, which counts *induced*
    orientation changes. The two differ whenever the Meek cascade amplifies a
    single false assertion, and quantifying that gap is the amplification study.

    Args:
        true_graph: The ground-truth DAG.
        bk: The knowledge to check.

    Returns:
        The violated constraints as node pairs.
    """
    raise NotImplementedError


def is_consistent_but_false(
    cpdag: MPDAG,
    true_graph: MPDAG,
    bk: BackgroundKnowledge,
) -> bool:
    """Whether ``bk`` passes the consistency check yet contradicts the truth.

    The defining predicate of this project: ``is_consistent(cpdag, bk)`` and
    ``not is_true(true_graph, bk)``. Every sample the perturbation module emits
    must satisfy it, and the tests in ``tests/test_perturb.py`` assert exactly
    that.

    Args:
        cpdag: The CPDAG output by discovery.
        true_graph: The ground-truth DAG, which must be a member of the
            equivalence class ``cpdag`` represents.
        bk: The knowledge to check.

    Returns:
        ``True`` if the knowledge is consistent but false.

    Raises:
        ValueError: If ``true_graph`` is not a consistent extension of ``cpdag``.
    """
    raise NotImplementedError
