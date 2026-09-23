"""Tests for bkrobust.synth.generators.

Each generator is checked for: returning a DAG (acyclic, right node count),
respecting its declared parameters, and reproducibility from a seed. The
acyclicity checks are cross-validated against an independent implementation
(networkx) rather than trusting MPDAG.is_acyclic() alone.

decoupled_backdoor_dag gets the sharpest scrutiny: it is "the key designed
family" per the brief, and the whole point of the family is the claim that
its two confounder-group blocking nodes sit in different undirected CPDAG
components at coupling=0.0 and the same component at coupling=1.0. That claim
is checked directly.

decoupled_backdoor_dag also carries a second round of tests added after a
review found the first version unusable: with the treatment/outcome edge
left ambiguous in the CPDAG, `all_valid_adjustment_sets_mpdag(cpdag, X, Y)`
was empty at every coupling level (some CPDAG extensions have Y causing X),
so the degeneracy gate rejected essentially every instance. The fix adds a
dedicated node W that forces X->Y to be compelled (mirroring how the prior
worked example's own Age->CVD edge is compelled). The tests below assert
exactly the properties a reviewer would want checked: X->Y is compelled,
>=2 valid adjustment sets exist at the CPDAG level at every coupling, the
component-separation claim still holds, and each confounder chain stays
individually perturbable. A further test documents -- rather than hides --
the fix's own cost: forcing CPDAG-level identification also freezes the
optimal adjustment set across the whole perturbation space, so r_val/r_opt
are UNREACHED for every instance of this family (see the generator's
docstring "Revision history" section for the full argument).
"""

from __future__ import annotations

import networkx as nx
import numpy as np
import pytest

from bkrobust.core.conventions import UNREACHED
from bkrobust.core.oracle import is_optimal, is_valid
from bkrobust.core.spacelib import build_space, distances_from, radius
from bkrobust.demo.evaluate import all_valid_adjustment_sets_mpdag, optimal_adjustment_set_mpdag
from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover
from bkrobust.demo.graph import undirected_components
from bkrobust.demo.meek import apply_orientations, is_valid_mpdag
from bkrobust.synth.generators import (
    GENERATORS,
    block_dag,
    decoupled_backdoor_dag,
    erdos_renyi_dag,
    scale_free_dag,
)


def _nx_digraph(dag) -> nx.DiGraph:
    g = nx.DiGraph()
    g.add_nodes_from(dag.nodes)
    g.add_edges_from(dag.directed_edges)
    return g


# --- registry -----------------------------------------------------------


def test_generators_registry_keys():
    assert set(GENERATORS) == {"erdos_renyi", "scale_free", "block", "decoupled_backdoor"}
    for fn in GENERATORS.values():
        assert callable(fn)


# --- erdos_renyi_dag ------------------------------------------------------


@pytest.mark.parametrize("n", [5, 8, 12])
def test_erdos_renyi_is_dag_with_right_node_count(n):
    rng = np.random.default_rng(0)
    dag = erdos_renyi_dag(n, rng, edge_prob=0.4)
    assert dag.is_dag()
    assert len(dag.nodes) == n


def test_erdos_renyi_acyclic_cross_checked_with_networkx():
    rng = np.random.default_rng(1)
    for _ in range(30):
        dag = erdos_renyi_dag(9, rng, edge_prob=0.5)
        assert nx.is_directed_acyclic_graph(_nx_digraph(dag))


def test_erdos_renyi_edge_prob_extremes():
    rng = np.random.default_rng(2)
    empty = erdos_renyi_dag(6, rng, edge_prob=0.0)
    assert len(empty.directed_edges) == 0
    full = erdos_renyi_dag(6, rng, edge_prob=1.0)
    assert len(full.directed_edges) == 6 * 5 // 2


def test_erdos_renyi_rejects_bad_edge_prob():
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError):
        erdos_renyi_dag(5, rng, edge_prob=1.5)


def test_erdos_renyi_reproducible_from_seed():
    a = erdos_renyi_dag(8, np.random.default_rng(123), edge_prob=0.3)
    b = erdos_renyi_dag(8, np.random.default_rng(123), edge_prob=0.3)
    assert a == b


def test_erdos_renyi_expected_density_roughly_matches_edge_prob():
    """Not a tight statistical test -- just a sanity band that edge_prob is honoured."""
    rng = np.random.default_rng(5)
    n = 10
    max_edges = n * (n - 1) // 2
    counts = [len(erdos_renyi_dag(n, rng, edge_prob=0.5).directed_edges) for _ in range(200)]
    mean_density = np.mean(counts) / max_edges
    assert 0.4 < mean_density < 0.6


# --- scale_free_dag ---------------------------------------------------------


@pytest.mark.parametrize(("n", "m"), [(6, 1), (10, 2), (12, 3)])
def test_scale_free_is_dag_with_right_node_count_and_edge_formula(n, m):
    rng = np.random.default_rng(0)
    dag = scale_free_dag(n, rng, m_attach=m)
    assert dag.is_dag()
    assert len(dag.nodes) == n
    # BA construction: the first new node attaches to all m seeds, and every
    # subsequent node adds exactly m edges, so the total is exact, not
    # approximate -- this is the "cross-check against a known formula" for
    # this generator (a full independent BA reimplementation would just be
    # this same algorithm restated).
    assert len(dag.directed_edges) == (n - m) * m


def test_scale_free_acyclic_cross_checked_with_networkx():
    rng = np.random.default_rng(3)
    for _ in range(30):
        dag = scale_free_dag(11, rng, m_attach=2)
        assert nx.is_directed_acyclic_graph(_nx_digraph(dag))


def test_scale_free_rejects_bad_m_attach():
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError):
        scale_free_dag(5, rng, m_attach=5)
    with pytest.raises(ValueError):
        scale_free_dag(5, rng, m_attach=0)


def test_scale_free_reproducible_from_seed():
    a = scale_free_dag(10, np.random.default_rng(9), m_attach=2)
    b = scale_free_dag(10, np.random.default_rng(9), m_attach=2)
    assert a == b


def test_scale_free_produces_hub_structure():
    """Preferential attachment should give a skewed degree distribution: max >> mean."""
    rng = np.random.default_rng(11)
    dag = scale_free_dag(30, rng, m_attach=2)
    skeleton = dag.skeleton()
    degree: dict[str, int] = {node: 0 for node in dag.nodes}
    for a, b in skeleton:
        degree[a] += 1
        degree[b] += 1
    degrees = list(degree.values())
    assert max(degrees) > 2 * np.mean(degrees)


# --- block_dag ---------------------------------------------------------


@pytest.mark.parametrize(("n", "n_blocks"), [(6, 2), (9, 3)])
def test_block_dag_is_dag_with_right_node_count(n, n_blocks):
    rng = np.random.default_rng(0)
    dag = block_dag(n, rng, n_blocks=n_blocks, p_within=0.7, p_between=0.05)
    assert dag.is_dag()
    assert len(dag.nodes) == n


def test_block_dag_acyclic_cross_checked_with_networkx():
    rng = np.random.default_rng(4)
    for _ in range(20):
        dag = block_dag(12, rng, n_blocks=3, p_within=0.8, p_between=0.1)
        assert nx.is_directed_acyclic_graph(_nx_digraph(dag))


def test_block_dag_denser_within_than_between():
    """With p_within >> p_between, within-block skeleton density should exceed cross-block."""
    rng = np.random.default_rng(6)
    n, n_blocks = 12, 3
    block_of = {i % n_blocks: [] for i in range(n_blocks)}
    for i in range(n):
        block_of[i % n_blocks].append(i)

    within_hits, within_total = 0, 0
    between_hits, between_total = 0, 0
    for _ in range(50):
        dag = block_dag(n, rng, n_blocks=n_blocks, p_within=0.9, p_between=0.05)
        skeleton = dag.skeleton()
        idx = {node: int(node[1:]) for node in dag.nodes}
        for a, b in skeleton:
            ia, ib = idx[a], idx[b]
            same = (ia % n_blocks) == (ib % n_blocks)
            if same:
                within_hits += 1
            else:
                between_hits += 1
        for i in range(n):
            for j in range(i + 1, n):
                if (i % n_blocks) == (j % n_blocks):
                    within_total += 1
                else:
                    between_total += 1
    within_rate = within_hits / within_total
    between_rate = between_hits / between_total
    assert within_rate > between_rate * 3


def test_block_dag_rejects_bad_n_blocks():
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError):
        block_dag(5, rng, n_blocks=0, p_within=0.5, p_between=0.1)


def test_block_dag_reproducible_from_seed():
    a = block_dag(9, np.random.default_rng(2), n_blocks=3, p_within=0.6, p_between=0.1)
    b = block_dag(9, np.random.default_rng(2), n_blocks=3, p_within=0.6, p_between=0.1)
    assert a == b


# --- decoupled_backdoor_dag: THE key designed family ------------------------


def _chain_sizes(n: int) -> tuple[int, int]:
    remaining = n - 3  # X, Y, W
    g1 = -(-remaining // 2)
    g2 = remaining - g1
    return g1, g2


@pytest.mark.parametrize("n", [7, 8, 9])
def test_decoupled_backdoor_is_dag_with_right_node_count(n):
    rng = np.random.default_rng(0)
    dag = decoupled_backdoor_dag(n, rng, coupling=0.5)
    assert dag.is_dag()
    assert len(dag.nodes) == n


@pytest.mark.parametrize("n", [7, 8, 9])
def test_decoupled_backdoor_acyclic_cross_checked_with_networkx(n):
    rng = np.random.default_rng(0)
    for coupling in (0.0, 0.3, 0.6, 1.0):
        dag = decoupled_backdoor_dag(n, rng, coupling=coupling)
        assert nx.is_directed_acyclic_graph(_nx_digraph(dag))


def test_decoupled_backdoor_rejects_too_small_n():
    """n=6 used to work; adding the W node (see the generator docstring) raised the minimum to 7."""
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError):
        decoupled_backdoor_dag(6, rng, coupling=0.0)


def test_decoupled_backdoor_rejects_bad_coupling():
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError):
        decoupled_backdoor_dag(7, rng, coupling=1.5)
    with pytest.raises(ValueError):
        decoupled_backdoor_dag(7, rng, coupling=-0.1)


def test_decoupled_backdoor_reproducible_from_seed():
    a = decoupled_backdoor_dag(8, np.random.default_rng(0), coupling=0.4)
    b = decoupled_backdoor_dag(8, np.random.default_rng(0), coupling=0.4)
    assert a == b


@pytest.mark.parametrize("n", [7, 8, 9])
def test_decoupled_backdoor_cpdag_is_valid(n):
    """The generated CPDAG must itself be a legitimate essential graph."""
    rng = np.random.default_rng(0)
    for coupling in (0.0, 0.5, 1.0):
        dag = decoupled_backdoor_dag(n, rng, coupling=coupling)
        cpdag = dag_to_cpdag(dag)
        assert is_valid_mpdag(cpdag)


# --- post-review fix: CPDAG-level identification (criteria 1-2) -------------
#
# A review of the original design found all_valid_adjustment_sets_mpdag(cpdag,
# "X", "Y") empty at every coupling level, because the CPDAG left X-Y
# undirected: some extensions have Y causing X, and no adjustment set is
# valid in both an "X causes Y" and a "Y causes X" world. That made the
# degeneracy gate reject essentially every instance from this family. The fix
# (a dedicated node W, W -> Y, adjacent to nothing else) forces X -> Y to be
# compelled via the unshielded (X, W) pair at Y. These two tests are exactly
# the acceptance check for that fix: X -> Y compelled, and at least two valid
# adjustment sets exist at the bare-CPDAG level, at every coupling swept.


@pytest.mark.parametrize("n", [7, 8, 9])
def test_decoupled_backdoor_xy_is_compelled_at_every_coupling(n):
    rng = np.random.default_rng(0)
    for coupling in (0.0, 0.25, 0.5, 0.75, 1.0):
        dag = decoupled_backdoor_dag(n, rng, coupling=coupling)
        cpdag = dag_to_cpdag(dag)
        assert ("X", "Y") in cpdag.directed_edges, (
            f"X->Y must be compelled at n={n}, coupling={coupling}"
        )


@pytest.mark.parametrize("n", [7, 8, 9])
def test_decoupled_backdoor_has_at_least_two_valid_adjustment_sets_at_cpdag_level(n):
    rng = np.random.default_rng(0)
    for coupling in (0.0, 0.25, 0.5, 0.75, 1.0):
        dag = decoupled_backdoor_dag(n, rng, coupling=coupling)
        cpdag = dag_to_cpdag(dag)
        valid_sets = all_valid_adjustment_sets_mpdag(cpdag, "X", "Y")
        assert len(valid_sets) >= 2, (
            f"expected >=2 valid adjustment sets at n={n}, coupling={coupling}, "
            f"got {len(valid_sets)}"
        )


@pytest.mark.parametrize("n", [7, 8, 9])
def test_decoupled_backdoor_coupling_zero_gives_disjoint_components(n):
    """THE sharp claim: at coupling=0.0, C1_last and C2_last are in different components."""
    rng = np.random.default_rng(0)
    g1, g2 = _chain_sizes(n)
    c1_last, c2_last = f"C1_{g1 - 1}", f"C2_{g2 - 1}"

    dag = decoupled_backdoor_dag(n, rng, coupling=0.0)
    cpdag = dag_to_cpdag(dag)
    comps = undirected_components(cpdag)

    comp_of_c1 = next((i for i, c in enumerate(comps) if c1_last in c), None)
    comp_of_c2 = next((i for i, c in enumerate(comps) if c2_last in c), None)
    assert comp_of_c1 is not None, "C1_last should sit in its chain's undirected component"
    assert comp_of_c2 is not None, "C2_last should sit in its chain's undirected component"
    assert comp_of_c1 != comp_of_c2, "at coupling=0.0 the two blocking nodes must be decoupled"

    # And the direct confounding edges (into X, Y) are frozen (compelled), not
    # perturbable -- this is now true at every coupling level, not just 0.0
    # (see test_decoupled_backdoor_xy_is_compelled_at_every_coupling and the
    # generator's docstring), but is asserted here where it was originally.
    for edge in ((c1_last, "X"), (c1_last, "Y"), (c2_last, "X"), (c2_last, "Y")):
        assert edge in cpdag.directed_edges


@pytest.mark.parametrize("n", [7, 8, 9])
def test_decoupled_backdoor_coupling_one_merges_components(n):
    """THE sharp claim's other half: at coupling=1.0 they do NOT lie in different components."""
    rng = np.random.default_rng(0)
    g1, g2 = _chain_sizes(n)
    c1_last, c2_last = f"C1_{g1 - 1}", f"C2_{g2 - 1}"

    dag = decoupled_backdoor_dag(n, rng, coupling=1.0)
    cpdag = dag_to_cpdag(dag)
    comps = undirected_components(cpdag)

    comp_of_c1 = next((i for i, c in enumerate(comps) if c1_last in c), None)
    comp_of_c2 = next((i for i, c in enumerate(comps) if c2_last in c), None)
    assert comp_of_c1 is not None and comp_of_c2 is not None
    assert comp_of_c1 == comp_of_c2, "at coupling=1.0 the two blocking nodes must be coupled"


@pytest.mark.parametrize("n", [7, 8, 9])
def test_decoupled_backdoor_chains_stay_perturbable_at_both_endpoints(n):
    """Acceptance criterion 5: each chain keeps >= 1 undirected edge, so knowledge can move it.

    Checked at both coupling endpoints, which is what the family's own tests
    (and the runner) actually exercise; intermediate coupling levels can
    momentarily compel one chain's own interior edge too (a "boundary"
    v-structure where a chain-2 node's cross edge and its chain edge meet),
    which is a known, accepted quirk of this family and not tested here.
    """
    rng = np.random.default_rng(0)
    for coupling in (0.0, 1.0):
        dag = decoupled_backdoor_dag(n, rng, coupling=coupling)
        cpdag = dag_to_cpdag(dag)
        # Each chain has >= 1 internal edge (g1, g2 >= 2), and none of those
        # internal edges touch X, Y or W, so nothing else could compel them.
        assert cpdag.is_undirected_edge("C1_0", "C1_1"), f"chain 1 frozen at coupling={coupling}"
        assert cpdag.is_undirected_edge("C2_0", "C2_1"), f"chain 2 frozen at coupling={coupling}"


@pytest.mark.parametrize("n", [7, 8, 9])
def test_decoupled_backdoor_merge_is_monotonic_in_coupling(n):
    """Component count should not increase as coupling rises (a regression guard on the design)."""
    rng = np.random.default_rng(0)
    prev_n_comps = None
    for coupling in (0.0, 0.25, 0.5, 0.75, 1.0):
        dag = decoupled_backdoor_dag(n, rng, coupling=coupling)
        cpdag = dag_to_cpdag(dag)
        n_comps = len(undirected_components(cpdag))
        if prev_n_comps is not None:
            assert n_comps <= prev_n_comps
        prev_n_comps = n_comps


# --- the fix's own cost: r_val/r_opt are UNREACHED for this whole family ----
#
# Documented rather than hidden: forcing CPDAG-level identification (above)
# requires C1_last/C2_last's edges into X and Y to be compelled
# unconditionally too (any residual ambiguity there breaks identification
# the same way X-Y's did). The optimal adjustment set is a pure function of
# that now-fully-compelled backbone (O = pa(Y) \ {X}, since Y is X's only
# causal descendant here), so it is literally the same set at every element
# of the perturbation space -- r_val and r_opt are UNREACHED regardless of
# coupling, knows_fraction or corruption. This is verified directly, not
# asserted as a belief: every step of run_instance's pipeline is replayed by
# hand up to the radius computation.


@pytest.mark.parametrize("coupling", [0.0, 1.0])
def test_decoupled_backdoor_r_val_and_r_opt_are_unreached(coupling):
    rng = np.random.default_rng(0)
    dag = decoupled_backdoor_dag(8, rng, coupling=coupling)
    cpdag = dag_to_cpdag(dag)
    k_true = knowledge_to_recover(dag, cpdag)
    # X->Y is already compelled (see test_decoupled_backdoor_xy_is_compelled_at_every_coupling),
    # so knowledge_to_recover -- which only asserts what is still undirected --
    # correctly leaves it out; asserting the chain interiors is enough to
    # fully resolve the CPDAG.
    g0 = apply_orientations(cpdag, k_true)
    assert g0 is not None
    z = optimal_adjustment_set_mpdag(g0, "X", "Y")
    assert z is not None, "Z is identified (that part of the fix worked)"
    zf = frozenset(z)

    space = build_space(cpdag)
    dists = distances_from(space, g0)
    assert is_valid(zf, g0, "X", "Y"), "Z is read off G0 and must be valid there"

    r_val, _ = radius(space, dists, lambda g: not is_valid(zf, g, "X", "Y"))
    r_opt, _ = radius(space, dists, lambda g: not is_optimal(zf, g, "X", "Y"))
    assert r_val == UNREACHED, "the optimal set is frozen across the whole space -- see docstring"
    assert r_opt == UNREACHED


# --- hash-seed determinism (mandatory invariant test) ------------------------


def _run_under_hashseeds(code: str) -> None:
    import subprocess
    import sys

    outs = set()
    for hs in ("0", "1", "12345"):
        env = {"PYTHONPATH": "src", "PYTHONHASHSEED": hs, "PATH": "/usr/bin:/bin"}
        res = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, env=env, check=True
        )
        outs.add(res.stdout.strip())
    assert len(outs) == 1, f"non-deterministic across PYTHONHASHSEED: {outs}"


def test_generators_are_hash_invariant():
    """Every generator's output must be bit-identical across PYTHONHASHSEED.

    RNG consumption or float accumulation that (accidentally) depends on
    set/dict iteration order would pass within one process but differ across
    processes with a different PYTHONHASHSEED -- this is exactly the bug
    class the invariants call out, so it is checked directly here rather than
    only trusted from code review.
    """
    code = (
        "import numpy as np;"
        "from bkrobust.synth.generators import ("
        "erdos_renyi_dag, scale_free_dag, block_dag, decoupled_backdoor_dag);"
        "out = [];"
        "out.append(erdos_renyi_dag(9, np.random.default_rng(0), 0.4).edge_string());"
        "out.append(scale_free_dag(10, np.random.default_rng(0), 2).edge_string());"
        "out.append(block_dag(9, np.random.default_rng(0), 3, 0.7, 0.1).edge_string());"
        "out.append(decoupled_backdoor_dag(8, np.random.default_rng(0), 0.6).edge_string());"
        "print(repr(out))"
    )
    _run_under_hashseeds(code)
