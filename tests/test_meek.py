"""Meek's rules against examples whose answers are known by hand.

Every assertion here should be verifiable without running the code. That is what
makes these tests worth having: an implementation checked only against its own
behaviour on random graphs can be confidently, consistently wrong.
"""

from __future__ import annotations

import pytest

from bkrobust.graphs import meek


def test_r1_orients_away_from_collider(mpdag, chain):
    """R1: given A -> B - C with A, C non-adjacent, B -> C is forced.

    Otherwise A -> B <- C would be a v-structure the data did not show.
    """
    graph = mpdag(chain, as_cpdag=True)
    graph.orient("A", "B")
    with pytest.raises(NotImplementedError):
        meek.meek_rule_1(graph)


def test_r2_prevents_cycle(mpdag, diamond):
    """R2: given A -> C -> B and A - B, A -> B is forced, since B -> A would cycle."""
    graph = mpdag(diamond, as_cpdag=True)
    with pytest.raises(NotImplementedError):
        meek.meek_rule_2(graph)


def test_r3_and_r4_fire_on_their_configurations(mpdag, diamond):
    """R3 and R4 fire only on their own configurations and are otherwise silent.

    R4 in particular does not fire on a bare CPDAG -- it needs background
    knowledge to have oriented something first, which is exactly this project's
    setting and why R4 cannot be skipped here as it sometimes is elsewhere.
    """
    graph = mpdag(diamond, as_cpdag=True)
    with pytest.raises(NotImplementedError):
        meek.meek_rule_3(graph)
    with pytest.raises(NotImplementedError):
        meek.meek_rule_4(graph)


def test_closure_is_idempotent(mpdag, diamond):
    """Closing an already-closed graph changes nothing."""
    graph = mpdag(diamond, as_cpdag=True)
    with pytest.raises(NotImplementedError):
        meek.meek_closure(graph)


def test_closure_of_collider_cpdag_is_itself(mpdag, collider):
    """The collider's CPDAG is fully oriented, so the closure is a fixpoint already."""
    graph = mpdag(collider, as_cpdag=True)
    with pytest.raises(NotImplementedError):
        meek.meek_closure(graph)


def test_closure_returns_none_on_cycle(mpdag, chain):
    """A knowledge set forcing a directed cycle must yield FAIL, i.e. None.

    Not an exception: FAIL is an expected outcome that the rejection sampler
    hits constantly, and making it exceptional would put a try/except in the
    sampler's inner loop.
    """
    graph = mpdag(chain, as_cpdag=True)
    with pytest.raises(NotImplementedError):
        meek.meek_closure(graph)


def test_apply_background_knowledge_is_order_independent(mpdag, diamond):
    """The FAIL verdict must not depend on the order constraints are imposed in.

    Algorithm 1 imposes constraints one at a time, so a naive implementation can
    return different answers for the same knowledge set depending on set
    iteration order. That would make consistency non-deterministic and every
    rejection-sampled result irreproducible.
    """
    from bkrobust.knowledge.base import BackgroundKnowledge

    graph = mpdag(diamond, as_cpdag=True)
    with pytest.raises(NotImplementedError):
        bk = BackgroundKnowledge.empty()
        meek.apply_background_knowledge(graph, bk)


def test_forced_orientations_excludes_imposed(mpdag, diamond):
    """Directly asserted edges are not counted as cascade-forced.

    The amplification factor is meaningless if the numerator includes the
    denominator.
    """
    from bkrobust.knowledge.base import BackgroundKnowledge

    graph = mpdag(diamond, as_cpdag=True)
    with pytest.raises(NotImplementedError):
        bk = BackgroundKnowledge.empty()
        meek.forced_orientations(graph, bk)


def test_orientable_edges_excludes_directed_cpdag_edges(mpdag, collider):
    """Edges already oriented in the CPDAG are not perturbable.

    Contradicting a v-structure is inconsistent, not consistent-but-false, so it
    is out of scope and must never enter the sampler's candidate pool.
    """
    graph = mpdag(collider, as_cpdag=True)
    with pytest.raises(NotImplementedError):
        meek.orientable_edges(graph)
