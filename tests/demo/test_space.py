"""Tests for :mod:`bkrobust.demo.space`.

Two worked examples anchor everything:

* ``CHAIN`` -- the fully undirected 3-node chain CPDAG ``A-B-C`` (no ``A-C``
  edge, no v-structure). Its Markov equivalence class has exactly the 3 DAGs
  without a collider at ``B``.
* ``CLIQUE`` -- the fully undirected 4-clique CPDAG on ``{W, X, Y, Z}``. Since
  every pair of nodes is adjacent, no acyclic orientation can ever create an
  unshielded collider, so its DAGs are exactly its ``4! = 24`` linear orders.

A judgement call, made explicit and tested directly (see
``test_self_validity_alone_is_unsound_for_the_chain`` and the
:mod:`bkrobust.demo.space` module docstring): the spec's formal definition of
the space as "same skeleton, contains every directed edge of the CPDAG, and
is a valid MPDAG" is, read completely literally with
``meek.is_valid_mpdag`` as the *only* filter, unsound for skeletons with an
unshielded triple. On ``CHAIN``, the fully oriented collider ``A->B<-C`` is
self-consistent in isolation (``is_valid_mpdag`` doesn't compare against the
original CPDAG at all) but belongs to a different Markov equivalence class
than the fully-undirected chain -- it is exactly the new v-structure Meek's
rule 1 forbids. Concretely this means the naive "9 raw states filtered only
by is_valid_mpdag" reading of the space has 7 elements for the chain, one of
which ([A->B<-C]) is not a subset of [CHAIN]. Since "the CPDAG is the unique
maximum" is a required invariant, ``enumerate_space`` (and every test below,
including the required 3**k cross-check) additionally requires genuine model
inclusion of the candidate's own extensions in the CPDAG's -- computed by
brute-force enumeration and comparison, never assumed. Under that (necessary)
reading the chain space has 6 elements, not the illustrative 4 one gets by
assuming every single-edge orientation must immediately resolve the whole
2-edge chain; it doesn't, because Meek's rule 1 only forces propagation for
an edge oriented *into* the shared node B, not out of it (orienting B->A
forces nothing about B-C, and vice versa) -- both intermediate MPDAGs are
demonstrably genuine, distinct, valid, non-maximal members of the space.
"""

from __future__ import annotations

from itertools import combinations, product

import pytest

from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import enumerate_dag_extensions, is_consistent_extension, is_valid_mpdag
from bkrobust.demo.space import (
    all_pairs_distances,
    atomic_moves,
    bfs_distances,
    check_metric_axioms,
    covering_pairs,
    enumerate_space,
    neighbour_graph,
    represented_dags,
)

# --- fixtures ----------------------------------------------------------------


def _chain_cpdag() -> MPDAG:
    return MPDAG(["A", "B", "C"], directed=[], undirected=[("A", "B"), ("B", "C")])


def _clique_cpdag() -> MPDAG:
    nodes = ["W", "X", "Y", "Z"]
    return MPDAG(nodes, directed=[], undirected=list(combinations(nodes, 2)))


@pytest.fixture(scope="module")
def chain_pipeline():
    """The full pipeline computed once for the chain CPDAG."""
    cpdag = _chain_cpdag()
    space = enumerate_space(cpdag)
    reps = represented_dags(space)
    covers = covering_pairs(space, reps)
    nbrs = neighbour_graph(space, covers)
    dists = all_pairs_distances(nbrs)
    axioms = check_metric_axioms(dists, space)
    return {
        "cpdag": cpdag,
        "space": space,
        "reps": reps,
        "covers": covers,
        "nbrs": nbrs,
        "dists": dists,
        "axioms": axioms,
    }


@pytest.fixture(scope="module")
def clique_pipeline():
    """The full pipeline computed once for the 4-clique CPDAG."""
    cpdag = _clique_cpdag()
    space = enumerate_space(cpdag)
    reps = represented_dags(space)
    covers = covering_pairs(space, reps)
    nbrs = neighbour_graph(space, covers)
    dists = all_pairs_distances(nbrs)
    axioms = check_metric_axioms(dists, space)
    return {
        "cpdag": cpdag,
        "space": space,
        "reps": reps,
        "covers": covers,
        "nbrs": nbrs,
        "dists": dists,
        "axioms": axioms,
    }


# --- enumerate_space -----------------------------------------------------


def test_chain_space_size_and_members(chain_pipeline):
    """CHAIN's space has exactly 6 elements.

    The CPDAG, 2 non-maximal single-edge refinements, and the 3 DAGs of the MEC.

    See the module docstring for why this is 6 and not the naively-expected
    4: orienting B-C alone or A-B alone (away from the shared node B) is a
    legitimate, distinct, valid, non-maximal member of the space, since
    Meek's rule 1 only forces further propagation when an edge is oriented
    *into* B, not out of it.
    """
    space = chain_pipeline["space"]
    edge_strings = {g.edge_string() for g in space}
    assert edge_strings == {
        "A-B B-C",  # the CPDAG itself
        "B->A B-C",  # A-B oriented away from B only
        "B->C A-B",  # B-C oriented away from B only
        "A->B B->C",  # chain A->B->C
        "B->A B->C",  # fork A<-B->C
        "B->A C->B",  # chain A<-B<-C
    }
    assert len(space) == 6


def test_clique_maximal_dag_count(clique_pipeline):
    """CLIQUE's maximal elements (fully oriented DAGs) number 4! = 24.

    A clique's v-structure-free acyclic orientations are exactly its linear
    orders, since every pair of nodes is adjacent (no unshielded triple can
    ever exist), and there are 4! linear orders on 4 nodes.
    """
    maximal = [g for g in clique_pipeline["space"] if len(g.undirected_edges) == 0]
    assert len(maximal) == 24
    # Every maximal element really is a full DAG (no undirected edges) and
    # acyclic (guaranteed by MPDAG's own constructor, re-asserted here).
    for g in maximal:
        assert g.is_dag()


def test_self_validity_alone_is_unsound_for_the_chain():
    """Direct evidence for the judgement call.

    On CHAIN, is_valid_mpdag alone (no comparison to the CPDAG) admits the
    collider A->B<-C, which belongs to a different Markov equivalence class
    than the fully-undirected chain.
    """
    collider = MPDAG(["A", "B", "C"], directed=[("A", "B"), ("C", "B")], undirected=[])
    assert is_valid_mpdag(collider)  # self-consistent in isolation
    cpdag = _chain_cpdag()
    # ... yet it is not a consistent extension of the chain CPDAG: orienting
    # both A-B and B-C into B is precisely the new v-structure forbidden by
    # background-knowledge closure starting from the fully undirected chain.
    assert not is_consistent_extension(collider, cpdag)


# --- represented_dags / model inclusion -----------------------------------


def test_cpdag_is_unique_maximum(chain_pipeline, clique_pipeline):
    """[cpdag] is a superset of [G] for every G in the space, on both examples."""
    for pipeline in (chain_pipeline, clique_pipeline):
        cpdag = pipeline["cpdag"]
        reps = pipeline["reps"]
        cpdag_ext = reps[cpdag]
        for g in pipeline["space"]:
            assert reps[g] <= cpdag_ext, f"{g.edge_string()} not subset of cpdag's extensions"
        # And it is the *unique* maximum: no other element has an equal or
        # larger extension set.
        for g in pipeline["space"]:
            if g == cpdag:
                continue
            assert reps[g] < cpdag_ext


def test_every_dag_in_space_is_minimal(chain_pipeline, clique_pipeline):
    """Every fully oriented DAG in the space has a singleton [G] (it represents only itself)."""
    for pipeline in (chain_pipeline, clique_pipeline):
        reps = pipeline["reps"]
        for g in pipeline["space"]:
            if len(g.undirected_edges) == 0:
                assert reps[g] == frozenset({g}), f"{g.edge_string()} is not minimal"


# --- covering pairs, re-verified independently ----------------------------


def _independent_covering_check(space, reps, lower, upper) -> bool:
    """Re-derive, from scratch, whether (lower, upper) is a genuine covering pair.

    Deliberately does not call bkrobust.demo.space.covering_pairs -- this
    re-verifies its output by direct, independent search over the space.
    """
    lower_ext = reps[lower]
    upper_ext = reps[upper]
    if not (lower_ext < upper_ext):
        return False
    for h in space:
        if h is lower or h is upper:
            continue
        if lower_ext < reps[h] < upper_ext:
            return False
    return True


def test_covering_pairs_reverified_independently(chain_pipeline, clique_pipeline):
    """Every reported covering pair is re-verified by an independent search.

    Each pair differs by >= 1 orientation and has nothing strictly between it,
    re-checked by a fresh search (not by calling covering_pairs' own
    internals) -- and nothing that *should* be reported as covering is
    missing.
    """
    for pipeline in (chain_pipeline, clique_pipeline):
        space = pipeline["space"]
        reps = pipeline["reps"]
        covers = pipeline["covers"]

        for lower, upper in covers:
            assert lower != upper
            assert lower.edge_string() != upper.edge_string()
            assert _independent_covering_check(space, reps, lower, upper)

        # Completeness: every pair the independent search calls a cover is
        # actually present in covers (checked on the smaller chain example
        # to keep this O(n^3) sanity pass cheap).
        if len(space) <= 10:
            for lower in space:
                for upper in space:
                    if lower is upper:
                        continue
                    if _independent_covering_check(space, reps, lower, upper):
                        assert (lower, upper) in covers


def test_covering_pairs_all_strict_subset(chain_pipeline, clique_pipeline):
    """Every covering pair is a strict model-inclusion pair.

    ``lower`` always has strictly fewer represented DAGs than ``upper``.
    """
    for pipeline in (chain_pipeline, clique_pipeline):
        reps = pipeline["reps"]
        for lower, upper in pipeline["covers"]:
            assert reps[lower] < reps[upper]


# --- neighbour graph / distances ------------------------------------------


def test_bfs_distance_from_cpdag_to_itself_is_zero(chain_pipeline, clique_pipeline):
    for pipeline in (chain_pipeline, clique_pipeline):
        nbrs = pipeline["nbrs"]
        cpdag = pipeline["cpdag"]
        dist = bfs_distances(nbrs, cpdag)
        assert dist[cpdag] == 0


def test_bfs_distances_omit_unreachable_not_infinity():
    """A vertex with no neighbours reports distance 0 to itself and nothing else."""
    isolated = MPDAG(["P", "Q"], directed=[], undirected=[("P", "Q")])
    other = MPDAG(["P", "Q"], directed=[("P", "Q")], undirected=[])
    nbrs = {isolated: set(), other: set()}
    dist = bfs_distances(nbrs, isolated)
    assert dist == {isolated: 0}
    assert other not in dist


# --- metric axioms -----------------------------------------------------


def test_metric_axioms_hold_over_all_triples(chain_pipeline, clique_pipeline):
    """Metric axioms hold for both examples; report exactly how many triples were checked."""
    for pipeline in (chain_pipeline, clique_pipeline):
        axioms = pipeline["axioms"]
        space = pipeline["space"]
        assert axioms["identity_of_indiscernibles"] is True
        assert axioms["symmetry"] is True
        assert axioms["triangle_inequality"] is True
        assert axioms["violations"] == []
        assert isinstance(axioms["n_triples_checked"], int)
        # The neighbour graph is connected in both examples (every pair of
        # MPDAGs in the space is mutually reachable via covering steps), so
        # literally every one of the n**3 ordered triples has all three
        # pairwise distances defined and gets checked.
        assert axioms["n_triples_checked"] == len(space) ** 3
        assert axioms["n_triples_checked"] > 0


def test_metric_axioms_report_smaller_count_when_disconnected():
    """n_triples_checked drops when part of the space is unreachable.

    Rather than silently treating a missing distance as some default value.
    """
    a = MPDAG(["A", "B"], directed=[], undirected=[("A", "B")])
    b = MPDAG(["A", "B"], directed=[("A", "B")], undirected=[])
    c = MPDAG(["C", "D"], directed=[], undirected=[("C", "D")])
    space = [a, b, c]
    nbrs = {a: {b}, b: {a}, c: set()}
    dists = all_pairs_distances(nbrs)
    axioms = check_metric_axioms(dists, space)
    # c is isolated: only triples drawn from {a, b} (plus c alone, trivially
    # d(c,c)=0) have every pairwise distance defined.
    assert axioms["n_triples_checked"] == 2 * 2 * 2 + 1  # {a,b}^3 plus (c,c,c)
    assert axioms["identity_of_indiscernibles"] is True
    assert axioms["symmetry"] is True
    assert axioms["triangle_inequality"] is True


# --- atomic_moves --------------------------------------------------------


def test_atomic_moves_single_assertion_is_one_move():
    g_from = MPDAG(["A", "B"], directed=[], undirected=[("A", "B")])
    g_to = MPDAG(["A", "B"], directed=[("A", "B")], undirected=[])
    moves = atomic_moves(g_from, g_to)
    assert moves == ["orient A->B"]


def test_atomic_moves_single_retraction_is_one_move():
    g_from = MPDAG(["A", "B"], directed=[("A", "B")], undirected=[])
    g_to = MPDAG(["A", "B"], directed=[], undirected=[("A", "B")])
    moves = atomic_moves(g_from, g_to)
    assert moves == ["un-orient A-B"]


def test_atomic_moves_flip_is_exactly_two_moves():
    """A flip is not a primitive move.

    A->B to B->A costs a retraction plus a fresh assertion: exactly 2 atomic
    moves, never 1.
    """
    g_from = MPDAG(["A", "B"], directed=[("A", "B")], undirected=[])
    g_to = MPDAG(["A", "B"], directed=[("B", "A")], undirected=[])
    moves = atomic_moves(g_from, g_to)
    assert moves == ["un-orient A-B", "orient B->A"]
    assert len(moves) == 2


def test_atomic_moves_no_change_is_zero_moves():
    g = MPDAG(["A", "B"], directed=[("A", "B")], undirected=[])
    assert atomic_moves(g, g) == []


def test_atomic_moves_mixed_case_covers_chain(chain_pipeline):
    """Check a mix of assertions and flips on the chain space.

    Moving from the CPDAG to the chain-forward DAG is two independent
    assertions (no flip needed); moving between the two opposite-chain DAGs
    is two flips (four moves).
    """
    space_by_edges = {g.edge_string(): g for g in chain_pipeline["space"]}
    cpdag = space_by_edges["A-B B-C"]
    chain_fwd = space_by_edges["A->B B->C"]
    chain_bwd = space_by_edges["B->A C->B"]

    moves_from_cpdag = atomic_moves(cpdag, chain_fwd)
    assert sorted(moves_from_cpdag) == sorted(["orient A->B", "orient B->C"])

    moves_between_dags = atomic_moves(chain_fwd, chain_bwd)
    # Both edges flip: 2 moves each = 4 total.
    assert len(moves_between_dags) == 4
    assert moves_between_dags.count("un-orient A-B") == 1
    assert moves_between_dags.count("un-orient B-C") == 1
    assert "orient B->A" in moves_between_dags
    assert "orient C->B" in moves_between_dags


def test_atomic_moves_requires_matching_skeleton():
    g1 = MPDAG(["A", "B"], directed=[], undirected=[("A", "B")])
    g2 = MPDAG(["A", "B", "C"], directed=[], undirected=[("A", "B"), ("B", "C")])
    with pytest.raises(ValueError):
        atomic_moves(g1, g2)


# --- cross-check against an independent brute-force enumeration ----------


def _independent_enumerate_space(cpdag: MPDAG) -> set[MPDAG]:
    """A from-scratch re-implementation of the 3**k brute-force enumeration.

    Written independently of bkrobust.demo.space.enumerate_space, used only
    to cross-check it in this test.

    As required by the brief: iterates all 3**k assignments of {leave
    undirected, orient one way, orient the other} to the k undirected edges,
    filtering by is_valid_mpdag. As demonstrated by
    test_self_validity_alone_is_unsound_for_the_chain, is_valid_mpdag alone
    is not sufficient on skeletons with an unshielded triple (it admits
    graphs from a different Markov equivalence class than cpdag's), so this
    independent enumeration also requires genuine model inclusion of each
    candidate's own DAG extensions in cpdag's -- the same necessary
    correction enumerate_space itself makes, verified here via a separately
    written loop.
    """
    undirected = sorted(cpdag.undirected_edges)
    k = len(undirected)
    base_directed = set(cpdag.directed_edges)
    found: set[MPDAG] = set()

    for state in product((0, 1, 2), repeat=k):
        directed = set(base_directed)
        remaining = []
        for (a, b), s in zip(undirected, state):  # noqa: B905
            if s == 1:
                directed.add((a, b))
            elif s == 2:
                directed.add((b, a))
            else:
                remaining.append((a, b))
        try:
            candidate = MPDAG(cpdag.nodes, directed, remaining)
        except ValueError:
            continue
        if not is_valid_mpdag(candidate):
            continue
        if all(is_consistent_extension(d, cpdag) for d in enumerate_dag_extensions(candidate)):
            found.add(candidate)
    return found


def test_enumerate_space_matches_independent_brute_force_chain():
    cpdag = _chain_cpdag()
    assert set(enumerate_space(cpdag)) == _independent_enumerate_space(cpdag)


def test_enumerate_space_matches_independent_brute_force_clique(clique_pipeline):
    cpdag = clique_pipeline["cpdag"]
    assert set(clique_pipeline["space"]) == _independent_enumerate_space(cpdag)


def test_naive_is_valid_mpdag_only_filter_disagrees_on_the_chain():
    """Documents, as a test, the exact size of the gap.

    Literally filtering the chain's 3**2 raw states by is_valid_mpdag alone
    (no model-inclusion check against the CPDAG) yields 7 states, one more
    than the true space of 6 -- the extra one is the collider A->B<-C.
    """
    cpdag = _chain_cpdag()
    undirected = sorted(cpdag.undirected_edges)
    naive: set[MPDAG] = set()
    for state in product((0, 1, 2), repeat=len(undirected)):
        directed = set()
        remaining = []
        for (a, b), s in zip(undirected, state):  # noqa: B905
            if s == 1:
                directed.add((a, b))
            elif s == 2:
                directed.add((b, a))
            else:
                remaining.append((a, b))
        candidate = MPDAG(cpdag.nodes, directed, remaining)
        if is_valid_mpdag(candidate):
            naive.add(candidate)

    assert len(naive) == 7
    correct = set(enumerate_space(cpdag))
    assert len(correct) == 6
    extra = naive - correct
    assert {g.edge_string() for g in extra} == {"A->B C->B"}
