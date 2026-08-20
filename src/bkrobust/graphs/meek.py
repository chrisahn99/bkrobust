"""Meek's orientation rules, the closure, and Algorithm 1 with FAIL detection.

Meek's four rules propagate orientations that are forced by the ones already
present, given that the result must stay acyclic and must not introduce a new
v-structure:

* **R1** ``a -> b - c``  with ``a`` and ``c`` non-adjacent  =>  ``b -> c``
  (otherwise ``a -> b <- c`` is a new v-structure).
* **R2** ``a -> c -> b`` and ``a - b``  =>  ``a -> b`` (otherwise a cycle).
* **R3** ``a - b``, ``a - c``, ``a - d``, ``c -> b``, ``d -> b``, ``c`` and
  ``d`` non-adjacent  =>  ``a -> b``.
* **R4** ``a - b``, ``a - c``, ``c -> d``, ``d -> b``, with ``c`` and ``b``
  non-adjacent  =>  ``a -> b``. Only bites in the presence of background
  knowledge, which is exactly this project's setting.

Meek's Algorithm 1 imposes a body of background knowledge one constraint at a
time, taking the closure after each, and reports **FAIL** if any constraint
cannot be imposed -- because it contradicts an already-forced orientation, or
because it would create a cycle. Termination without FAIL is precisely what
current practice means by "the knowledge is consistent", and it is the only
check the knowledge ever receives. Note what it does not say: nothing about
whether the knowledge is *true*. That gap is the subject of this package.

Implementation note: the b-LOAD Meek/MPDAG closure is the intended starting
point. Do not vendor an implementation from a third-party discovery library;
see ``docs/ANONYMITY_CHECKLIST.md`` and the licence notes in ``pyproject.toml``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from bkrobust.graphs.mpdag import MPDAG

if TYPE_CHECKING:  # pragma: no cover
    from bkrobust.graphs.mpdag import Edge, Node
    from bkrobust.knowledge.base import BackgroundKnowledge


def meek_rule_1(graph: MPDAG) -> set[Edge]:
    """Return the orientations forced by R1 in a single pass.

    Args:
        graph: The graph to inspect. Not mutated.

    Returns:
        Directed edges ``(tail, head)`` that R1 forces and that are not already
        present.
    """
    raise NotImplementedError


def meek_rule_2(graph: MPDAG) -> set[Edge]:
    """Return the orientations forced by R2 in a single pass. Does not mutate ``graph``."""
    raise NotImplementedError


def meek_rule_3(graph: MPDAG) -> set[Edge]:
    """Return the orientations forced by R3 in a single pass. Does not mutate ``graph``."""
    raise NotImplementedError


def meek_rule_4(graph: MPDAG) -> set[Edge]:
    """Return the orientations forced by R4 in a single pass. Does not mutate ``graph``."""
    raise NotImplementedError


def meek_closure(graph: MPDAG) -> MPDAG | None:
    """Apply R1-R4 to a fixpoint.

    Repeatedly applies the four rules until no further edge is oriented. The
    result is the unique maximally oriented PDAG consistent with the input.

    Args:
        graph: The graph to close. Not mutated; a new graph is returned.

    Returns:
        The closed MPDAG, or ``None`` on **FAIL** -- that is, if the rules force
        a directed cycle or force an edge in both directions. ``None`` is the
        signal that the input was not consistent to begin with.
    """
    raise NotImplementedError


def apply_background_knowledge(
    cpdag: MPDAG,
    bk: BackgroundKnowledge,
) -> MPDAG | None:
    """Impose background knowledge on a CPDAG -- Meek's Algorithm 1.

    Imposes each constraint in ``bk`` in turn, taking the Meek closure after
    each, and aborts the moment a constraint cannot be imposed. Constraints are
    processed in a deterministic order so that the FAIL/no-FAIL verdict does not
    depend on iteration order of the underlying containers.

    Args:
        cpdag: The CPDAG output by discovery. Not mutated.
        bk: The knowledge to impose: required edges, forbidden edges, tiers and
            ancestral constraints.

    Returns:
        The resulting MPDAG, or ``None`` on FAIL.

    Note:
        A non-``None`` return means only that ``bk`` is *consistent* with
        ``cpdag``. It carries no information about whether ``bk`` is true. Use
        :func:`bkrobust.knowledge.perturb.sample_consistent_but_false` to
        construct knowledge that this function accepts and that is nevertheless
        wrong.
    """
    raise NotImplementedError


def forced_orientations(
    cpdag: MPDAG,
    bk: BackgroundKnowledge,
) -> set[Edge]:
    """Return the edges oriented by cascade rather than directly imposed.

    Partitions the orientations gained by
    :func:`apply_background_knowledge` into those that ``bk`` states outright
    and those that Meek's rules then force. Only the latter are returned.

    This is the measurement the amplification study rests on: ``k`` imposed
    constraints produce ``m >= k`` oriented edges, and the ratio ``m / k`` is
    how far the consequences of a single false assertion travel. See
    :mod:`bkrobust.knowledge.cascade`.

    Args:
        cpdag: The CPDAG the knowledge is imposed on.
        bk: The knowledge to impose.

    Returns:
        Directed edges present in the closure, undirected in ``cpdag``, and not
        directly asserted by ``bk``. Empty if imposing ``bk`` FAILs.
    """
    raise NotImplementedError


def orientable_edges(cpdag: MPDAG) -> set[Edge]:
    """Return the undirected edges of ``cpdag`` that could be oriented either way.

    The candidate pool the perturbation samplers draw from: asserting an
    orientation contrary to a *directed* CPDAG edge is inconsistent, not
    consistent-but-false, and so is out of scope for this project.

    Args:
        cpdag: The CPDAG to inspect.

    Returns:
        One ``(a, b)`` pair per undirected edge, in canonical endpoint order.
    """
    raise NotImplementedError


def is_meek_closed(graph: MPDAG) -> bool:
    """Whether no Meek rule fires on ``graph``.

    A cheap invariant check for tests and for asserting that a graph handed
    between modules really is an MPDAG and not a half-processed PDAG.
    """
    raise NotImplementedError


def orientation_conflict(cpdag: MPDAG, bk: BackgroundKnowledge) -> tuple[Node, Node] | None:
    """Return the first constraint that cannot be imposed, or ``None``.

    Diagnostic companion to :func:`apply_background_knowledge`: when that
    function returns ``None``, this says which constraint broke it. Used to
    give the rejection sampler in :mod:`bkrobust.knowledge.perturb` an
    actionable reason for each rejection rather than a bare count.

    Args:
        cpdag: The CPDAG the knowledge is imposed on.
        bk: The knowledge to impose.

    Returns:
        The offending node pair, or ``None`` if ``bk`` is consistent.
    """
    raise NotImplementedError
