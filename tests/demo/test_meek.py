"""Tests for :mod:`bkrobust.demo.meek`.

Hand-checkable rule instances first (each rule's docstring in ``meek.py``
states the pattern it implements; the instances here are chosen to be
verifiable by inspection), then the structural properties -- order
independence, idempotence, FAIL detection, extension enumeration, chordality,
and end-to-end validity -- that a demo relying on this module needs to trust.
"""

from __future__ import annotations

from collections.abc import Callable
from itertools import combinations

import networkx as nx
import numpy as np
import pytest

from bkrobust.demo.graph import MPDAG, Edge, Node
from bkrobust.demo.meek import (
    _meek_closure_ordered,
    apply_orientations,
    enumerate_dag_extensions,
    is_chordal_components,
    is_consistent_extension,
    is_meek_closed,
    is_valid_mpdag,
    meek_closure,
    meek_rule_1,
    meek_rule_2,
    meek_rule_3,
    meek_rule_4,
)

# --- individual rules ------------------------------------------------------


def test_rule1_fires_when_nonadjacent():
    """R1: A -> B, B - C, A and C non-adjacent => B -> C is forced."""
    g = MPDAG(["A", "B", "C"], directed=[("A", "B")], undirected=[("B", "C")])
    assert meek_rule_1(g) == {("B", "C")}


def test_rule1_silent_when_adjacent():
    """R1 does not fire when A and C are adjacent: no v-structure is at risk."""
    g = MPDAG(
        ["A", "B", "C"],
        directed=[("A", "B")],
        undirected=[("B", "C"), ("A", "C")],
    )
    assert meek_rule_1(g) == set()


def test_rule2_fires():
    """R2: A -> C -> B, A - B => A -> B is forced (else B -> A would cycle)."""
    g = MPDAG(
        ["A", "B", "C"],
        directed=[("A", "C"), ("C", "B")],
        undirected=[("A", "B")],
    )
    assert meek_rule_2(g) == {("A", "B")}


def test_rule2_silent_without_directed_path():
    """R2 needs a directed a->c->b path; a lone undirected edge forces nothing."""
    g = MPDAG(["A", "B", "C"], directed=[("A", "C")], undirected=[("A", "B")])
    assert meek_rule_2(g) == set()


def test_rule3_fires_on_minimal_instance():
    """R3: A-B, A-C, A-D, C->B, D->B, C and D non-adjacent => A->B is forced.

    Instance: nodes A, B, C, D; undirected A-B, A-C, A-D; directed C->B, D->B;
    C and D not adjacent to each other. Orienting A-B as B->A would, via R1
    applied to C->B/A-... chains, ultimately compel an unshielded v-structure
    at A -- so A-B must be A->B.
    """
    g = MPDAG(
        ["A", "B", "C", "D"],
        directed=[("C", "B"), ("D", "B")],
        undirected=[("A", "B"), ("A", "C"), ("A", "D")],
    )
    assert ("A", "B") in meek_rule_3(g)


def test_rule3_silent_near_miss_when_c_d_adjacent():
    """Same instance as above but C and D are adjacent: R3 must not fire."""
    g = MPDAG(
        ["A", "B", "C", "D"],
        directed=[("C", "B"), ("D", "B")],
        undirected=[("A", "B"), ("A", "C"), ("A", "D"), ("C", "D")],
    )
    assert ("A", "B") not in meek_rule_3(g)


def test_rule4_fires_on_minimal_instance():
    """R4: A-B, A-C, C->D, D->B, C and B non-adjacent => A->B is forced.

    Instance: nodes A, B, C, D; undirected A-B, A-C; directed C->D, D->B; C
    and B not adjacent. Orienting A-B as B->A, together with A-C and the
    directed path C->D->B, would compel a new unshielded v-structure -- so
    A-B must be A->B.
    """
    g = MPDAG(
        ["A", "B", "C", "D"],
        directed=[("C", "D"), ("D", "B")],
        undirected=[("A", "B"), ("A", "C")],
    )
    assert ("A", "B") in meek_rule_4(g)


def test_rule4_silent_near_miss_when_c_b_adjacent():
    """Same instance as above but C and B are adjacent: R4 must not fire."""
    g = MPDAG(
        ["A", "B", "C", "D"],
        directed=[("C", "D"), ("D", "B")],
        undirected=[("A", "B"), ("A", "C"), ("B", "C")],
    )
    assert ("A", "B") not in meek_rule_4(g)


# --- closure: order independence and idempotence ---------------------------


def _random_pdag(rng: np.random.Generator, n: int, p_edge: float, p_directed: float) -> MPDAG:
    """A random small MPDAG respecting every container invariant.

    Each unordered node pair independently becomes a directed edge (in a
    random direction), an undirected edge, or no edge -- so it can never
    violate the "not both directed and undirected" or "no mutually directed
    pair" invariants. The resulting directed part may or may not be acyclic;
    :func:`meek_closure` is expected to handle both (returning ``None`` in the
    cyclic case), and order-independence is meaningful either way.
    """
    nodes = [str(i) for i in range(n)]
    directed: list[Edge] = []
    undirected: list[Edge] = []
    for i, j in combinations(range(n), 2):
        r = rng.random()
        if r >= p_edge:
            continue
        if rng.random() < p_directed:
            if rng.random() < 0.5:
                directed.append((nodes[i], nodes[j]))
            else:
                directed.append((nodes[j], nodes[i]))
        else:
            undirected.append((nodes[i], nodes[j]))
    return MPDAG(nodes, directed, undirected)


def _random_key(rng: np.random.Generator, nodes: list[Node]) -> Callable[[Edge], object]:
    """A random total order over all ordered node pairs, as a sort key."""
    pairs = [(a, b) for a in nodes for b in nodes if a != b]
    perm = rng.permutation(len(pairs))
    rank: dict[Edge, int] = {pairs[i]: int(perm[i]) for i in range(len(pairs))}
    return lambda e: rank[e]


def test_closure_is_order_independent():
    """meek_closure agrees across 10 random application orders, on 200 PDAGs.

    Confluence of Meek's rules means the fixpoint (or the FAIL) reached does
    not depend on the order forced orientations are applied in. This drives
    that order through ``_meek_closure_ordered`` directly via random sort
    keys, rather than relying on Python's incidental set-iteration order,
    which would make the test vacuous (a single fixed iteration order could
    pass by accident).
    """
    master = np.random.default_rng(0)
    for _ in range(200):
        n = int(master.integers(3, 6))
        p_edge = float(master.uniform(0.2, 0.7))
        p_directed = float(master.uniform(0.2, 0.8))
        g = _random_pdag(master, n, p_edge, p_directed)

        reference = _meek_closure_ordered(g, key=lambda e: e)
        for _ in range(10):
            key = _random_key(master, list(g.nodes))
            result = _meek_closure_ordered(g, key=key)
            assert result == reference


def test_closure_idempotent():
    """meek_closure(meek_closure(g)) == meek_closure(g) for several graphs."""
    graphs = [
        MPDAG(["A", "B", "C"], directed=[("A", "B")], undirected=[("B", "C")]),
        MPDAG(
            ["A", "B", "C", "D"],
            directed=[("C", "B"), ("D", "B")],
            undirected=[("A", "B"), ("A", "C"), ("A", "D")],
        ),
        MPDAG(["A", "B", "C"], directed=[], undirected=[("A", "B"), ("B", "C")]),
    ]
    master = np.random.default_rng(1)
    for _ in range(50):
        n = int(master.integers(3, 6))
        p_edge = float(master.uniform(0.2, 0.6))
        p_directed = float(master.uniform(0.2, 0.8))
        graphs.append(_random_pdag(master, n, p_edge, p_directed))

    for g in graphs:
        once = meek_closure(g)
        if once is None:
            continue
        twice = meek_closure(once)
        assert twice == once


def test_meek_closure_does_not_mutate_input():
    """meek_closure must not mutate its argument."""
    g = MPDAG(["A", "B", "C"], directed=[("A", "B")], undirected=[("B", "C")])
    before_directed = g.directed_edges
    before_undirected = g.undirected_edges
    meek_closure(g)
    assert g.directed_edges == before_directed
    assert g.undirected_edges == before_undirected


# --- FAIL detection ----------------------------------------------------


def test_fail_when_rules_force_a_cycle():
    """A case where the rules force the same edge in both directions.

    Nodes A, B, C, D. Directed A->B, B->C (a path) and D->C. Undirected A-C.
    D and A are not adjacent.

    R2 sees A->B->C with undirected A-C and forces A->C.
    R1 sees D->C, undirected C-A, and D/A non-adjacent, and forces C->A.
    Both fire in the same pass: A-C is forced in both directions at once,
    which is exactly the FAIL condition, and is caught without the directed
    part of the input graph itself containing a cycle (A->B, B->C, D->C is
    acyclic on its own).
    """
    g = MPDAG(
        ["A", "B", "C", "D"],
        directed=[("A", "B"), ("B", "C"), ("D", "C")],
        undirected=[("A", "C")],
    )
    assert g.is_acyclic()
    assert ("A", "C") in meek_rule_2(g)
    assert ("C", "A") in meek_rule_1(g)
    assert meek_closure(g) is None


def test_fail_when_input_already_cyclic():
    """A directed part that is already cyclic FAILs immediately."""
    g = MPDAG(["A", "B", "C"], directed=[("A", "B"), ("B", "C"), ("C", "A")])
    assert meek_closure(g) is None


def test_apply_orientations_contradiction_fails():
    """Imposing an orientation opposite an existing directed edge FAILs."""
    g = MPDAG(["A", "B"], directed=[("A", "B")])
    assert apply_orientations(g, [("B", "A")]) is None


def test_apply_orientations_nonadjacent_pair_fails():
    """Imposing an orientation between a non-adjacent pair FAILs."""
    g = MPDAG(["A", "B", "C"], directed=[("A", "B")])
    assert apply_orientations(g, [("A", "C")]) is None


def test_apply_orientations_deterministic_regardless_of_input_order():
    """apply_orientations processes in sorted order, so input order is moot."""
    g = MPDAG(["A", "B", "C", "D"], undirected=[("A", "B"), ("C", "D")])
    forward = apply_orientations(g, [("A", "B"), ("C", "D")])
    backward = apply_orientations(g, [("C", "D"), ("A", "B")])
    assert forward is not None
    assert forward == backward


# --- DAG extension enumeration -------------------------------------------


def test_enumerate_three_node_chain_excludes_collider():
    """The 3-node chain CPDAG A-B-C (no v-structure) has exactly 3 extensions.

    A-B-C fully undirected, with A and C non-adjacent, is the CPDAG for the
    Markov equivalence class of chain/fork DAGs on this skeleton. Its three
    members are A->B->C, C->B->A, and A<-B->C; the fourth orientation,
    A->B<-C, is the collider and must be excluded since it is not a
    v-structure of the original (fully undirected) graph.
    """
    g = MPDAG(["A", "B", "C"], undirected=[("A", "B"), ("B", "C")])
    extensions = enumerate_dag_extensions(g)
    strings = {d.edge_string() for d in extensions}
    assert strings == {"A->B B->C", "B->A C->B", "B->A B->C"}
    for d in extensions:
        assert not (d.is_directed_edge("A", "B") and d.is_directed_edge("C", "B"))


def test_extensions_are_dags_and_consistent():
    """Every returned extension is a DAG and passes is_consistent_extension.

    Uses a richer graph: an existing v-structure X->W<-Y plus a hanging
    undirected edge W-Z, so the check exercises both "keep existing directed
    edges" and "no new v-structure" simultaneously. Of the two orientations
    of W-Z, only W->Z avoids introducing a new collider at W (Z->W would
    make X->W<-Z and Y->W<-Z new, unlicensed v-structures, since X/Z and
    Y/Z are non-adjacent and neither was directed in the original graph) --
    so exactly one extension is expected.
    """
    g = MPDAG(
        ["W", "X", "Y", "Z"],
        directed=[("X", "W"), ("Y", "W")],
        undirected=[("W", "Z")],
    )
    extensions = enumerate_dag_extensions(g)
    assert len(extensions) == 1
    d = extensions[0]
    assert d.is_dag()
    assert is_consistent_extension(d, g)
    assert d.is_directed_edge("X", "W")
    assert d.is_directed_edge("Y", "W")
    assert d.is_directed_edge("W", "Z")


def test_chordless_4cycle_has_no_extensions():
    """A chordless 4-cycle A-B-C-D-A admits zero consistent DAG extensions.

    Any acyclic orientation of a chordless cycle necessarily creates at
    least one collider (a node with two non-adjacent in-neighbours), by a
    standard parity argument: the orientation cannot be consistently
    "clockwise" (that would be a directed cycle) so it must reverse
    direction an even number of times going around, producing at least one
    local sink with two incoming arrows. Since none of these edges is
    directed in the original graph, every such collider is a new,
    unlicensed v-structure -- so no orientation survives. This is also why
    the chordality gate in :func:`is_valid_mpdag` matters: this same graph
    is Meek-closed but not chordal, and correspondingly has no extensions.
    """
    g = MPDAG(
        ["A", "B", "C", "D"],
        undirected=[("A", "B"), ("B", "C"), ("C", "D"), ("D", "A")],
    )
    assert enumerate_dag_extensions(g) == []


def test_extensions_of_chordal_component_are_dags_and_consistent():
    """A chordal skeleton (4-cycle plus one chord) has multiple extensions.

    Every extension returned must be a DAG, pass is_consistent_extension,
    and -- since the source graph has no directed edges at all -- have no
    two non-adjacent parents at any node.
    """
    g = MPDAG(
        ["A", "B", "C", "D"],
        undirected=[("A", "B"), ("B", "C"), ("C", "D"), ("D", "A"), ("A", "C")],
    )
    extensions = enumerate_dag_extensions(g)
    assert len(extensions) > 0
    for d in extensions:
        assert d.is_dag()
        assert is_consistent_extension(d, g)
        for b in d.nodes:
            parents = sorted(d.parents(b))
            for a, c in combinations(parents, 2):
                assert d.has_edge(a, c)


# --- chordality --------------------------------------------------------


def test_chordal_cross_check_against_networkx():
    """is_chordal_components agrees with networkx.is_chordal on 100 graphs."""
    rng = np.random.default_rng(2)
    for _ in range(100):
        n = int(rng.integers(3, 9))
        p = float(rng.uniform(0.1, 0.6))
        seed = int(rng.integers(0, 1_000_000))
        nx_graph = nx.gnp_random_graph(n, p, seed=seed)
        nodes = [str(v) for v in nx_graph.nodes]
        edges = [(str(a), str(b)) for a, b in nx_graph.edges]
        mp = MPDAG(nodes, undirected=edges)
        assert is_chordal_components(mp) == nx.is_chordal(nx_graph)


def test_is_chordal_components_true_for_chordal_and_false_for_4cycle():
    """A triangle-plus-chord is chordal; a chordless 4-cycle is not."""
    chordal = MPDAG(
        ["A", "B", "C", "D"],
        undirected=[("A", "B"), ("B", "C"), ("C", "D"), ("D", "A"), ("A", "C")],
    )
    four_cycle = MPDAG(
        ["A", "B", "C", "D"],
        undirected=[("A", "B"), ("B", "C"), ("C", "D"), ("D", "A")],
    )
    assert is_chordal_components(chordal) is True
    assert is_chordal_components(four_cycle) is False


# --- end-to-end validity -------------------------------------------------


def test_cpdag_is_meek_closed_and_valid():
    """A genuine CPDAG (v-structure plus a disjoint undetermined edge) is closed.

    X->W<-Y (v-structure, X and Y non-adjacent) has no undirected edges
    touching W, so no rule can fire from it. A disjoint undirected edge P-Q
    (representing the 2-node MEC {P->Q, Q->P}) similarly has no directed
    edges anywhere near it. No rule premise is satisfiable, so this graph is
    vacuously Meek-closed, chordal in its only (2-node) undirected component,
    acyclic, and admits DAG extensions -- i.e. it is a valid MPDAG.
    """
    g = MPDAG(
        ["W", "X", "Y", "P", "Q"],
        directed=[("X", "W"), ("Y", "W")],
        undirected=[("P", "Q")],
    )
    assert is_meek_closed(g)
    assert is_valid_mpdag(g)


def test_non_closed_graph_is_not_valid_mpdag():
    """A graph on which a rule still fires is not Meek-closed, hence invalid."""
    g = MPDAG(["A", "B", "C"], directed=[("A", "B")], undirected=[("B", "C")])
    assert not is_meek_closed(g)
    assert not is_valid_mpdag(g)


def test_cyclic_directed_part_is_not_valid_mpdag():
    """A graph whose directed part contains a cycle is not a valid MPDAG."""
    g = MPDAG(["A", "B", "C"], directed=[("A", "B"), ("B", "C"), ("C", "A")])
    assert not is_valid_mpdag(g)


def test_non_chordal_component_is_not_valid_mpdag():
    """A chordless 4-cycle undirected component fails the chordality gate."""
    g = MPDAG(
        ["A", "B", "C", "D"],
        undirected=[("A", "B"), ("B", "C"), ("C", "D"), ("D", "A")],
    )
    assert is_meek_closed(g)
    assert not is_chordal_components(g)
    assert not is_valid_mpdag(g)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
