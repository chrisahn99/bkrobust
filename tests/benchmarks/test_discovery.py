"""A learned CPDAG round-trips through the MPDAG class and is Meek-closed.

The estimated panel converts causal-learn's endpoint matrix into the
repository's graph and then runs the same instruments the ledger ran on the
oracle CPDAG. Two things must hold for that to mean anything: the conversion
must be lossless on a conflict-free graph, and the graph the instruments see
must be Meek-closed, since every certificate is computed by re-closing it.
"""

from __future__ import annotations

import numpy as np
import pytest

from bkrobust.benchmarks.discovery import (
    from_endpoint_matrix,
    partition_labels,
    sample_ancestral,
    structure_scores,
    to_endpoint_matrix,
)
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import is_meek_closed, meek_closure


def _collider_chain() -> MPDAG:
    """``A -> B -> C <- D`` with ``C -> E`` and ``E -> F``: one v-structure, one open tail."""
    nodes = ("A", "B", "C", "D", "E", "F")
    edges = (("A", "B"), ("B", "C"), ("D", "C"), ("C", "E"), ("E", "F"))
    return MPDAG(nodes, directed=edges)


def test_hand_built_matrix_converts_and_counts_conflicts():
    nodes = ["A", "B", "C", "D"]
    m = np.zeros((4, 4), dtype=int)
    m[0, 1], m[1, 0] = -1, 1  # A -> B
    m[1, 2], m[2, 1] = -1, -1  # B -- C
    m[2, 3], m[3, 2] = 1, 1  # C <-> D, a conflict
    g, counts = from_endpoint_matrix(m, nodes)
    assert g.is_directed_edge("A", "B")
    assert g.is_undirected_edge("B", "C")
    assert g.is_undirected_edge("C", "D")
    assert counts == {"directed": 1, "undirected": 1, "conflict": 1, "other": 0}


def test_conflict_free_graph_round_trips_exactly():
    cpdag = dag_to_cpdag(_collider_chain())
    back, counts = from_endpoint_matrix(to_endpoint_matrix(cpdag), cpdag.nodes)
    assert back == cpdag
    assert counts["conflict"] == 0
    assert is_meek_closed(back)


def test_partition_labels_give_the_oracle_an_index_of_one():
    from sklearn.metrics import adjusted_rand_score

    cpdag = dag_to_cpdag(_collider_chain())
    labels = partition_labels(cpdag)
    assert len(labels) == len(cpdag.nodes)
    assert adjusted_rand_score(labels, partition_labels(cpdag)) == 1.0
    # A -- B is the only chain component: A and B share a label, nobody else does
    idx = {v: i for i, v in enumerate(cpdag.nodes)}
    assert labels[idx["A"]] == labels[idx["B"]]
    assert len(set(labels)) == len(cpdag.nodes) - 1


def test_learned_cpdag_is_meek_closed_and_round_trips():
    pytest.importorskip("causallearn")
    from bkrobust.benchmarks.discovery import learn_graph

    dag = _collider_chain()
    weights = {e: 1.0 for e in dag.directed_edges}
    noise = {v: 1.0 for v in dag.nodes}
    data = sample_ancestral(dag, weights, noise, 3000, np.random.default_rng(0))
    matrix = learn_graph(data, "pc", alpha=0.01)
    learned, counts = from_endpoint_matrix(matrix, dag.nodes)
    # lossless on what the algorithm returned
    back, _ = from_endpoint_matrix(to_endpoint_matrix(learned), dag.nodes)
    assert back == learned
    # the graph the instruments see is Meek-closed, whether or not it arrived so
    scores, used = structure_scores(learned, dag, dag_to_cpdag(dag))
    assert scores["closure_ok"] == 1
    assert is_meek_closed(used)
    assert meek_closure(used) == used
    # at n = 3000 on six variables PC recovers the skeleton and the collider
    assert used.skeleton() == dag.skeleton()
    assert used.is_directed_edge("B", "C") and used.is_directed_edge("D", "C")
    assert scores["skeleton_precision"] == 1.0 and scores["skeleton_recall"] == 1.0
    assert counts["conflict"] == 0
