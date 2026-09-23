"""Sampling background knowledge that is consistent but false.

This is the workhorse of every experiment in the project. Everything else
measures what happens downstream of a sample drawn here, so the sampling
contract needs to be exact.

**The contract.** :func:`sample_consistent_but_false` returns a
:class:`~bkrobust.knowledge.base.BackgroundKnowledge` satisfying all four of:

1. **Consistent.** ``apply_background_knowledge(cpdag, bk)`` does not FAIL.
   This is what the practitioner's check would report, and passing it is what
   makes the sample interesting rather than merely broken.
2. **False.** At least one constraint of ``bk`` is violated by ``true_graph``.
3. **At radius delta.** The distance between the MPDAG that ``bk`` induces and
   the one the true knowledge induces equals ``delta`` under the requested
   metric -- or is at most ``delta`` when ``exact_radius`` is off, which the
   tier kind needs because moving one node rarely lands on an exact flip count.
4. **Of the requested kind.** Every false constraint belongs to ``kind``, so a
   sweep over kinds is a sweep over one variable at a time.

**Why rejection sampling.** There is no useful direct parameterisation of the
consistent-but-false set: consistency is a property of the Meek closure, which
is a global fixpoint, so whether a candidate qualifies is only known after
running it. The samplers therefore propose from a structured candidate pool --
the orientable edges, the tier assignments -- and reject. Two consequences the
callers must respect:

* Rejection rates are informative, not noise. A configuration where almost
  every draw is rejected is one where the CPDAG already pins down nearly
  everything, and that fact belongs in the results. It is recorded in
  ``metadata["n_rejections"]`` on every returned sample.
* The samplers can fail. At large ``delta`` on a well-determined CPDAG the
  target set is genuinely empty, and after ``max_rejections`` the sampler
  raises :class:`PerturbationExhaustedError` rather than looping or silently
  returning something off-radius.
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    import numpy as np

    from bkrobust.graphs.mpdag import MPDAG, Node
    from bkrobust.knowledge.base import BackgroundKnowledge
    from bkrobust.knowledge.taxonomy import MisspecificationKind


class PerturbationExhaustedError(RuntimeError):
    """Raised when no consistent-but-false knowledge at the requested radius was found.

    Carries the rejection count and the last rejection reason so a caller can
    tell "this configuration has an empty target set" apart from "the budget
    was too small", which are different findings and want different write-ups.
    """


def sample_consistent_but_false(
    true_graph: MPDAG,
    cpdag: MPDAG,
    delta: int,
    kind: MisspecificationKind | str,
    rng: np.random.Generator,
    *,
    treatment: Node | None = None,
    outcome: Node | None = None,
    distance: str = "orientation_flips",
    exact_radius: bool = True,
    max_rejections: int = 10_000,
    scope: dict[str, object] | None = None,
) -> BackgroundKnowledge:
    """Draw one knowledge set that is consistent with ``cpdag`` and false of ``true_graph``.

    Satisfies the four-part contract documented at module level. Proposals come
    from a pool determined by ``kind`` -- orientable edges for ``ORIENTATION``,
    undetermined node pairs for ``ANCESTRAL``, tier reassignments for ``TIER``
    -- and are rejected until one is both consistent and false at radius
    ``delta``.

    Args:
        true_graph: The ground-truth DAG. Must be a consistent extension of
            ``cpdag``.
        cpdag: The CPDAG the knowledge will be imposed on.
        delta: The perturbation radius, in units of ``distance``. ``delta=0``
            is a special case: it returns *true* knowledge, the control arm, and
            requirement 2 is waived.
        kind: Which misspecification kind to sample.
        rng: Seeded NumPy generator. The sole source of randomness -- no module
            may touch global RNG state.
        treatment: Treatment node, when ``scope`` restricts sampling relative to
            the target pair.
        outcome: Outcome node, likewise.
        distance: Metric that ``delta`` is measured in; see
            :mod:`bkrobust.graphs.distances`.
        exact_radius: Require distance exactly ``delta`` rather than at most.
        max_rejections: Give up after this many rejected proposals.
        scope: Sampling restrictions from the perturbation config -- the
            ``scope`` block of ``configs/perturbation/*.yaml``.

    Returns:
        The sampled knowledge, with ``source="sampled"`` and ``metadata``
        carrying ``kind``, ``delta``, ``realised_delta`` and ``n_rejections``.

    Raises:
        PerturbationExhaustedError: If no qualifying sample was found within
            ``max_rejections``.
        ValueError: If ``true_graph`` is not a consistent extension of
            ``cpdag``, or ``delta`` is negative, or ``kind`` is unknown.
    """
    raise NotImplementedError


def sample_consistent_and_true(
    true_graph: MPDAG,
    cpdag: MPDAG,
    n_constraints: int,
    kind: MisspecificationKind | str,
    rng: np.random.Generator,
    *,
    scope: dict[str, object] | None = None,
) -> BackgroundKnowledge:
    """Draw knowledge of the same shape and size as a perturbation, but true.

    The control arm. Matching on constraint count matters: knowledge shrinks the
    equivalence class whether or not it is correct, and some of that shrinkage
    changes ``O*`` on its own. Without a size-matched true arm, an experiment
    cannot attribute a change to falsity rather than to informativeness.

    Args:
        true_graph: The ground-truth DAG.
        cpdag: The CPDAG the knowledge will be imposed on.
        n_constraints: How many constraints to assert.
        kind: Which constraint kind to draw.
        rng: Seeded generator.
        scope: Sampling restrictions, as in
            :func:`sample_consistent_but_false`.

    Returns:
        Knowledge that is true of ``true_graph`` and consistent with ``cpdag``.

    Raises:
        PerturbationExhaustedError: If ``cpdag`` does not admit ``n_constraints``
            distinct true constraints of that kind.
    """
    raise NotImplementedError


def enumerate_consistent_but_false(
    true_graph: MPDAG,
    cpdag: MPDAG,
    delta: int,
    kind: MisspecificationKind | str,
    *,
    distance: str = "orientation_flips",
    limit: int | None = None,
) -> Iterator[BackgroundKnowledge]:
    """Yield every consistent-but-false knowledge set at radius ``delta``.

    The exhaustive counterpart to :func:`sample_consistent_but_false`, used by
    the exact path of :mod:`bkrobust.theory.radius`, where a radius is defined
    as a minimum over the whole ball and sampling would only ever give an upper
    bound on it.

    Combinatorial in ``delta`` and in the number of undirected edges; guard
    calls with ``limit`` and with the ``radius.exact_max_nodes`` config knob.

    Args:
        true_graph: The ground-truth DAG.
        cpdag: The CPDAG the knowledge will be imposed on.
        delta: The exact radius to enumerate at.
        kind: Which misspecification kind to enumerate.
        distance: Metric that ``delta`` is measured in.
        limit: Stop after this many; ``None`` enumerates all.

    Yields:
        Knowledge sets in a deterministic order.
    """
    raise NotImplementedError


def perturbation_pool(
    cpdag: MPDAG,
    kind: MisspecificationKind | str,
    *,
    treatment: Node | None = None,
    outcome: Node | None = None,
    scope: dict[str, object] | None = None,
) -> Sequence[tuple[Node, Node]]:
    """Return the candidate constraints a sampler of this kind may propose from.

    Exposed separately because the pool size determines whether a target set can
    be non-empty at a given ``delta``, and reporting "the pool had 3 members, so
    delta=5 was unreachable" is more useful than reporting a rejection count.

    Args:
        cpdag: The CPDAG the knowledge will be imposed on.
        kind: Which misspecification kind the pool is for.
        treatment: Treatment node, for target-relative scoping.
        outcome: Outcome node, likewise.
        scope: Restrictions from the perturbation config.

    Returns:
        Candidate node pairs in a deterministic order.
    """
    raise NotImplementedError
