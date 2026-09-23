"""Distance metrics on graphs and knowledge.

Distances define what a radius means, so an error here silently rescales every
reported radius rather than producing a visible failure.
"""

from __future__ import annotations

import pytest

from bkrobust.graphs import distances


def test_distance_to_self_is_zero(mpdag, canonical_graphs):
    """Every metric is zero between a graph and itself."""
    for spec in canonical_graphs.values():
        graph = mpdag(spec)
        with pytest.raises(NotImplementedError):
            assert distances.structural_hamming_distance(graph, graph) == 0
            assert distances.orientation_flip_count(graph, graph) == 0


def test_shd_is_symmetric(mpdag, chain, fork):
    """SHD does not depend on argument order."""
    g1, g2 = mpdag(chain), mpdag(fork)
    with pytest.raises(NotImplementedError):
        assert distances.structural_hamming_distance(g1, g2) == (
            distances.structural_hamming_distance(g2, g1)
        )


def test_chain_and_fork_differ_by_one_flip(mpdag, chain, fork):
    """A -> B -> C and B -> A, B -> C differ in the orientation of A - B only.

    One flip, one SHD. They are Markov equivalent, which is why background
    knowledge is the only thing that could separate them -- and why getting that
    knowledge wrong is undetectable from the data.
    """
    g1, g2 = mpdag(chain), mpdag(fork)
    with pytest.raises(NotImplementedError):
        assert distances.orientation_flip_count(g1, g2) == 1
        assert distances.structural_hamming_distance(g1, g2) == 1


def test_orientation_flips_ignore_undirected_edges(mpdag, chain):
    """An undirected edge is not oriented "differently" from a directed one.

    SHD counts that disagreement; the flip count does not. The two metrics are
    meant to differ here, and a radius reported without naming its metric is
    ambiguous by exactly this amount.
    """
    dag = mpdag(chain)
    cpdag = mpdag(chain, as_cpdag=True)
    with pytest.raises(NotImplementedError):
        assert distances.orientation_flip_count(dag, cpdag) == 0
        assert distances.structural_hamming_distance(dag, cpdag) > 0


def test_poset_distance_is_bounded(mpdag, chain):
    """The poset distance lies in [0, 1] and is zero for identical closures."""
    from bkrobust.knowledge.base import BackgroundKnowledge

    cpdag = mpdag(chain, as_cpdag=True)
    with pytest.raises(NotImplementedError):
        bk = BackgroundKnowledge.empty()
        d = distances.poset_distance(bk, bk, cpdag)
        assert d == 0.0


def test_distances_reject_mismatched_node_sets(mpdag, chain, diamond):
    """Graphs on different node sets are not comparable."""
    g1, g2 = mpdag(chain), mpdag(diamond)
    with pytest.raises((NotImplementedError, ValueError)):
        distances.structural_hamming_distance(g1, g2)


def test_get_distance_rejects_unknown_metric():
    """An unknown metric name fails loudly rather than falling back to a default."""
    with pytest.raises((NotImplementedError, KeyError)):
        distances.get_distance("not_a_metric")
