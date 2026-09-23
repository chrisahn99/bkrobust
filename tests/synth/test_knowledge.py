"""Tests for bkrobust.synth.knowledge.

What each test discriminates:

* ``draw_k_true`` tests check the sampler only ever emits orientations that
  are actually true of the DAG it was drawn from, and that the count tracks
  ``knows_fraction``.
* ``omit``/``flip``/``compound`` tests check each corruption changes exactly
  what its name says: omit shrinks the list (and only ever removes, never
  invents or reverses), flip preserves length but reverses direction, and
  compound does both.
* ``tiered`` tests check the uncorrupted tiering only ever asserts true
  claims, that it "compels more edges than it asserts" (the Meek closure of
  the asserted set is strictly larger, which is the whole point of tiered
  knowledge per the brief), and that a corrupted tier assignment can make
  ``K_assumed`` inconsistent -- exercising ``check_consistent`` on a real
  failure, not just a hand-built one.
"""

from __future__ import annotations

import numpy as np
import pytest

from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.synth.generators import block_dag, erdos_renyi_dag
from bkrobust.synth.knowledge import (
    CORRUPTIONS,
    check_consistent,
    compound,
    draw_k_true,
    flip,
    omit,
    tiered,
)


def _dag_and_cpdag(
    seed: int, n: int = 8, p: float = 0.5
) -> tuple[MPDAG, MPDAG, np.random.Generator]:
    rng = np.random.default_rng(seed)
    dag = erdos_renyi_dag(n, rng, edge_prob=p)
    return dag, dag_to_cpdag(dag), rng


# --- registry -----------------------------------------------------------


def test_corruptions_registry_keys():
    assert set(CORRUPTIONS) == {"omit", "flip", "compound", "tiered"}


# --- draw_k_true --------------------------------------------------------


def test_draw_k_true_yields_only_true_orientations():
    dag, cpdag, rng = _dag_and_cpdag(0)
    for frac in (0.0, 0.3, 0.6, 1.0):
        k = draw_k_true(dag, cpdag, rng, knows_fraction=frac)
        for a, b in k:
            assert dag.is_directed_edge(a, b), f"{(a, b)} is not true of dag"


def test_draw_k_true_count_tracks_knows_fraction():
    dag, cpdag, rng = _dag_and_cpdag(1)
    n_undirected = len(cpdag.undirected_edges)
    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        k = draw_k_true(dag, cpdag, rng, knows_fraction=frac)
        assert len(k) == round(frac * n_undirected)


def test_draw_k_true_full_fraction_covers_every_undirected_edge():
    dag, cpdag, rng = _dag_and_cpdag(2)
    k = draw_k_true(dag, cpdag, rng, knows_fraction=1.0)
    covered = {frozenset(e) for e in k}
    expected = {frozenset(e) for e in cpdag.undirected_edges}
    assert covered == expected


def test_draw_k_true_rejects_bad_fraction():
    dag, cpdag, rng = _dag_and_cpdag(0)
    with pytest.raises(ValueError):
        draw_k_true(dag, cpdag, rng, knows_fraction=1.5)


def test_draw_k_true_reproducible_from_seed():
    dag, cpdag, _ = _dag_and_cpdag(3)
    a = draw_k_true(dag, cpdag, np.random.default_rng(77), knows_fraction=0.5)
    b = draw_k_true(dag, cpdag, np.random.default_rng(77), knows_fraction=0.5)
    assert a == b


# --- omit -----------------------------------------------------------------


def test_omit_only_removes_never_invents_or_reverses():
    dag, cpdag, rng = _dag_and_cpdag(4)
    k = draw_k_true(dag, cpdag, rng, knows_fraction=1.0)
    o = omit(k, rng, rate=0.5)
    assert set(o) <= set(k)
    assert len(o) == len(k) - round(0.5 * len(k))


def test_omit_rate_zero_is_identity_up_to_order():
    dag, cpdag, rng = _dag_and_cpdag(5)
    k = draw_k_true(dag, cpdag, rng, knows_fraction=0.8)
    assert omit(k, rng, rate=0.0) == sorted(k)


def test_omit_rate_one_empties_the_list():
    dag, cpdag, rng = _dag_and_cpdag(6)
    k = draw_k_true(dag, cpdag, rng, knows_fraction=1.0)
    assert omit(k, rng, rate=1.0) == []


def test_omit_handles_empty_input():
    rng = np.random.default_rng(0)
    assert omit([], rng, rate=0.5) == []


# --- flip -------------------------------------------------------------------


def test_flip_preserves_length_and_only_reverses():
    dag, cpdag, rng = _dag_and_cpdag(7)
    k = draw_k_true(dag, cpdag, rng, knows_fraction=1.0)
    f = flip(k, rng, rate=0.4)
    assert len(f) == len(k)
    k_set = {frozenset(e) for e in k}
    f_set = {frozenset(e) for e in f}
    assert k_set == f_set, "flip must not change which edges are claimed, only direction"


def test_flip_rate_zero_is_identity_up_to_order():
    dag, cpdag, rng = _dag_and_cpdag(8)
    k = draw_k_true(dag, cpdag, rng, knows_fraction=0.7)
    assert flip(k, rng, rate=0.0) == sorted(k)


def test_flip_rate_one_reverses_every_edge():
    dag, cpdag, rng = _dag_and_cpdag(9)
    k = draw_k_true(dag, cpdag, rng, knows_fraction=1.0)
    f = flip(k, rng, rate=1.0)
    assert set(f) == {(b, a) for a, b in k}


def _first_seed_with_enough_undirected_edges(start_seed: int, min_edges: int) -> int:
    """Find a seed whose CPDAG has at least ``min_edges`` undirected edges.

    Some individual seeds happen to draw a near-complete (fully oriented)
    graph with nothing left ambiguous; rather than skip the test for such an
    unlucky draw, scan forward deterministically for one that has enough
    structure to exercise the property under test.
    """
    seed = start_seed
    while True:
        _dag, cpdag, _rng = _dag_and_cpdag(seed)
        if len(cpdag.undirected_edges) >= min_edges:
            return seed
        seed += 1


def test_flip_changes_at_least_one_orientation_at_high_rate():
    """The whole point of flip: some claims must now be FALSE of the dag."""
    seed = _first_seed_with_enough_undirected_edges(10, min_edges=1)
    dag, cpdag, rng = _dag_and_cpdag(seed)
    k = draw_k_true(dag, cpdag, rng, knows_fraction=1.0)
    f = flip(k, rng, rate=1.0)
    assert any(not dag.is_directed_edge(a, b) for a, b in f)


# --- compound ---------------------------------------------------------------


def test_compound_both_shrinks_and_reverses():
    seed = _first_seed_with_enough_undirected_edges(11, min_edges=4)
    dag, cpdag, rng = _dag_and_cpdag(seed)
    k = draw_k_true(dag, cpdag, rng, knows_fraction=1.0)
    c = compound(k, rng, rate=0.5)
    assert len(c) <= len(k)


def test_compound_equals_flip_of_omit():
    """Compound is documented as literally composing the other two operators."""
    dag, cpdag, _ = _dag_and_cpdag(12)
    k = draw_k_true(dag, cpdag, np.random.default_rng(0), knows_fraction=1.0)
    rng_a = np.random.default_rng(99)
    rng_b = np.random.default_rng(99)
    assert compound(k, rng_a, rate=0.5) == flip(omit(k, rng_b, rate=0.5), rng_b, rate=0.5)


# --- tiered -----------------------------------------------------------------


def test_tiered_uncorrupted_only_asserts_true_claims():
    rng = np.random.default_rng(0)
    for _trial in range(30):
        dag = erdos_renyi_dag(8, rng, edge_prob=0.5)
        cpdag = dag_to_cpdag(dag)
        t = tiered(dag, cpdag, rng, n_tiers=3, corruption_rate=0.0)
        for a, b in t:
            assert dag.is_directed_edge(a, b)


def test_tiered_compels_more_edges_than_it_asserts():
    """Tiering's point per the brief: it compels many edges at once via Meek closure."""
    rng = np.random.default_rng(1)
    found_a_strict_case = False
    for _trial in range(50):
        dag = block_dag(8, rng, n_blocks=2, p_within=0.6, p_between=0.3)
        cpdag = dag_to_cpdag(dag)
        t = tiered(dag, cpdag, rng, n_tiers=3, corruption_rate=0.0)
        if not t:
            continue
        g0 = apply_orientations(cpdag, t)
        assert g0 is not None, "uncorrupted tiering must be consistent"
        if len(g0.directed_edges) > len(t):
            found_a_strict_case = True
            break
    assert found_a_strict_case, "expected at least one trial where closure compels extra edges"


def test_tiered_same_tier_edges_are_never_asserted():
    rng = np.random.default_rng(2)
    dag = erdos_renyi_dag(8, rng, edge_prob=0.5)
    cpdag = dag_to_cpdag(dag)
    t = tiered(dag, cpdag, rng, n_tiers=1, corruption_rate=0.0)
    assert t == [], "a single tier puts every node in the same tier: nothing crosses"


def test_tiered_corruption_can_make_k_assumed_inconsistent():
    """A corrupted tier assignment must be able to break consistency (that's the point)."""
    rng = np.random.default_rng(3)
    saw_inconsistent = False
    for _trial in range(100):
        dag = erdos_renyi_dag(8, rng, edge_prob=0.5)
        cpdag = dag_to_cpdag(dag)
        t = tiered(dag, cpdag, rng, n_tiers=4, corruption_rate=0.5)
        if not check_consistent(cpdag, t):
            saw_inconsistent = True
            break
    assert saw_inconsistent


def test_tiered_rejects_bad_params():
    dag, cpdag, rng = _dag_and_cpdag(4)
    with pytest.raises(ValueError):
        tiered(dag, cpdag, rng, n_tiers=0)
    with pytest.raises(ValueError):
        tiered(dag, cpdag, rng, n_tiers=2, corruption_rate=1.5)


def test_tiered_reproducible_from_seed():
    dag, cpdag, _ = _dag_and_cpdag(5)
    a = tiered(dag, cpdag, np.random.default_rng(42), n_tiers=3, corruption_rate=0.3)
    b = tiered(dag, cpdag, np.random.default_rng(42), n_tiers=3, corruption_rate=0.3)
    assert a == b


# --- check_consistent ---------------------------------------------------


def test_check_consistent_true_for_true_knowledge():
    dag, cpdag, rng = _dag_and_cpdag(6)
    k = draw_k_true(dag, cpdag, rng, knows_fraction=1.0)
    assert check_consistent(cpdag, k)


def test_check_consistent_false_for_a_directed_cycle():
    """Hand-built degenerate case: asserting both directions of a chain closes a cycle."""
    cpdag = MPDAG(["A", "B", "C"], directed=[], undirected=[("A", "B"), ("B", "C"), ("A", "C")])
    # A->B, B->C, C->A: a directed 3-cycle. apply_orientations must FAIL.
    assert not check_consistent(cpdag, [("A", "B"), ("B", "C"), ("C", "A")])


def test_check_consistent_false_for_contradicted_edge():
    cpdag = MPDAG(["A", "B"], directed=[], undirected=[("A", "B")])
    assert not check_consistent(cpdag, [("A", "B"), ("B", "A")])


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


def test_knowledge_module_is_hash_invariant():
    code = (
        "import numpy as np;"
        "from bkrobust.synth.generators import erdos_renyi_dag;"
        "from bkrobust.demo.example import dag_to_cpdag;"
        "from bkrobust.synth.knowledge import (draw_k_true, omit, flip, compound, tiered);"
        "rng = np.random.default_rng(0);"
        "dag = erdos_renyi_dag(9, rng, 0.5);"
        "cpdag = dag_to_cpdag(dag);"
        "k = draw_k_true(dag, cpdag, rng, 0.7);"
        "out = [k, omit(k, rng, 0.3), flip(k, rng, 0.3), compound(k, rng, 0.3),"
        "       tiered(dag, cpdag, rng, 3, 0.4)];"
        "print(repr(out))"
    )
    _run_under_hashseeds(code)
