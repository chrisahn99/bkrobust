"""The ancestral reduction rule, against graphs whose answers are known by hand.

Soundness is the requirement: an emitted orientation must hold in every DAG of
the class in which the claimed ancestor really is one. Completeness is not, and
one test records a graph where the rule stays silent on an orientation that is
in fact entailed, so that the incompleteness is a documented choice rather than
a surprise.
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import enumerate_dag_extensions, is_meek_closed
from bkrobust.knowledge.ancestral import (
    NOT_REDUCIBLE,
    REDUCED,
    first_edges,
    last_edges,
    reduce_ancestral_claim,
)


def _class_with_ancestor(cpdag: MPDAG, a: str, b: str) -> list[MPDAG]:
    """Every DAG of the class in which ``a`` is an ancestor of ``b``."""
    return [d for d in enumerate_dag_extensions(cpdag) if a in d.ancestors(b)]


def _assert_sound(cpdag: MPDAG, a: str, b: str) -> None:
    red = reduce_ancestral_claim(cpdag, a, b)
    dags = _class_with_ancestor(cpdag, a, b)
    if not dags:
        # the claim contradicts the class; the rule must emit nothing
        assert red.orientations == ()
        return
    for tail, head in red.orientations:
        assert all(d.is_directed_edge(tail, head) for d in dags), (a, b, tail, head)


def test_rule_fires_at_both_ends_on_the_chain():
    """A - C - B: "A causes B" leaves only A -> C -> B in the class."""
    cpdag = MPDAG("ABC", undirected=[("A", "C"), ("C", "B")])
    assert is_meek_closed(cpdag)
    red = reduce_ancestral_claim(cpdag, "A", "B")
    assert red.status == REDUCED
    assert red.reason == "unique_first_and_last_edge"
    assert set(red.orientations) == {("A", "C"), ("C", "B")}
    dags = _class_with_ancestor(cpdag, "A", "B")
    assert len(dags) == 1
    _assert_sound(cpdag, "A", "B")


def test_rule_fires_on_the_longer_chain_and_meek_closes_the_middle():
    """A - C - D - B: the two end edges fire; C -> D is left to the closure."""
    cpdag = MPDAG("ABCD", undirected=[("A", "C"), ("C", "D"), ("D", "B")])
    red = reduce_ancestral_claim(cpdag, "A", "B")
    assert set(red.orientations) == {("A", "C"), ("D", "B")}
    _assert_sound(cpdag, "A", "B")


def test_rule_must_not_fire_when_the_paths_diverge():
    """Two triangles sharing the edge C - D, A on one side and B on the other.

    Possibly directed paths from A to B start with A - C or A - D and end with
    C - B or D - B, so neither end is unique and the rule emits nothing. That
    is right for the first edge: the class contains DAGs with A -> C and with
    C -> A in which A is still an ancestor of B. It is silent on the last
    edges even though both C -> B and D -> B are in fact entailed, since a DAG
    with B -> C would need C -> A and A -> D -> B, a cycle. The rule is sound
    and not complete, and this is the graph that shows the difference.
    """
    dag = MPDAG("ABCD", directed=[("A", "C"), ("A", "D"), ("C", "D"), ("C", "B"), ("D", "B")])
    cpdag = dag_to_cpdag(dag)
    assert not cpdag.directed_edges  # no v-structure, everything open
    red = reduce_ancestral_claim(cpdag, "A", "B")
    assert red.status == NOT_REDUCIBLE
    assert red.reason == "end_edges_diverge"
    assert red.orientations == ()
    assert red.n_first_edges == 2 and red.n_last_edges == 2
    dags = _class_with_ancestor(cpdag, "A", "B")
    assert any(d.is_directed_edge("A", "C") for d in dags)
    assert any(d.is_directed_edge("C", "A") for d in dags)
    # the entailment the rule leaves on the table, on record
    assert all(d.is_directed_edge("C", "B") and d.is_directed_edge("D", "B") for d in dags)


def test_no_possibly_directed_path_is_a_contradiction_not_a_guess():
    """E -> A <- C compels both edges into A, so nothing leaves A: "A causes B" has no path."""
    cpdag = MPDAG("ABCE", directed=[("E", "A"), ("C", "A")], undirected=[("C", "B")])
    assert is_meek_closed(cpdag)
    red = reduce_ancestral_claim(cpdag, "A", "B")
    assert red.status == NOT_REDUCIBLE
    assert red.reason == "no_possibly_directed_path"
    assert red.orientations == ()
    assert _class_with_ancestor(cpdag, "A", "B") == []


def test_unique_first_edge_fires_when_the_last_edge_is_already_compelled():
    """Same graph, "B causes A": the only path is B - C -> A, so B -> C and nothing else."""
    cpdag = MPDAG("ABCE", directed=[("E", "A"), ("C", "A")], undirected=[("C", "B")])
    red = reduce_ancestral_claim(cpdag, "B", "A")
    assert red.status == REDUCED
    assert red.reason == "unique_first_edge"
    assert red.orientations == (("B", "C"),)
    _assert_sound(cpdag, "B", "A")


def test_end_edge_sets_are_read_off_reachability():
    cpdag = MPDAG("ABCD", undirected=[("A", "C"), ("C", "D"), ("D", "B")])
    assert first_edges(cpdag, "A", "B") == {"C"}
    assert last_edges(cpdag, "A", "B") == {"D"}
    assert first_edges(cpdag, "B", "A") == {"D"}
    assert last_edges(cpdag, "B", "A") == {"C"}


def test_rejects_degenerate_claims():
    cpdag = MPDAG("AB", undirected=[("A", "B")])
    with pytest.raises(ValueError):
        reduce_ancestral_claim(cpdag, "A", "A")
    with pytest.raises(ValueError):
        reduce_ancestral_claim(cpdag, "A", "Z")


def test_sound_on_random_small_classes():
    """Brute force: on random CPDAGs every emitted orientation holds wherever the ancestry does."""
    rng = np.random.default_rng(20260912)
    n_reduced = 0
    for _ in range(40):
        n = int(rng.integers(4, 7))
        nodes = [f"v{i}" for i in range(n)]
        order = list(rng.permutation(n))
        directed = [
            (nodes[order[i]], nodes[order[j]])
            for i, j in itertools.combinations(range(n), 2)
            if rng.random() < 0.45
        ]
        cpdag = dag_to_cpdag(MPDAG(nodes, directed=directed))
        if len(cpdag.undirected_edges) > 10:
            continue
        for a, b in itertools.permutations(nodes, 2):
            if cpdag.has_edge(a, b):
                continue
            _assert_sound(cpdag, a, b)
            if reduce_ancestral_claim(cpdag, a, b).status == REDUCED:
                n_reduced += 1
    assert n_reduced > 0  # the sweep exercised the firing branch, not only the silent one
