"""Tests for the running example and the three analyst knowledge states.

These guard design properties the report's narrative depends on. If one fails,
the narrative is wrong, not just the code.
"""

from __future__ import annotations

from bkrobust.demo.evaluate import is_valid_adjustment_set_mpdag, optimal_adjustment_set_mpdag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations, is_valid_mpdag
from bkrobust.demo.scenario import (
    K_TRUE,
    OUTCOME,
    PIVOT,
    TREATMENT,
    scenarios,
    true_dag,
)


def _cpdag() -> MPDAG:
    return dag_to_cpdag(true_dag())


def test_truth_is_a_dag():
    assert true_dag().is_dag()


def test_cpdag_has_five_undirected_edges():
    """The design gate: a 4-6 edge chordal component. Five here."""
    assert len(_cpdag().undirected_edges) == 5


def test_cpdag_is_valid_mpdag():
    assert is_valid_mpdag(_cpdag())


def test_frozen_edges_are_never_perturbed():
    """Edges compelled in the CPDAG stay directed in every scenario."""
    cpdag = _cpdag()
    for s in scenarios().values():
        g = apply_orientations(cpdag, s["knowledge"])
        assert g is not None
        assert cpdag.directed_edges <= g.directed_edges


def test_k_true_recovers_the_truth():
    """The four clinical claims imply the ground-truth DAG exactly."""
    assert apply_orientations(_cpdag(), K_TRUE) == true_dag()


def test_k_true_is_redundant_two_claims_suffice():
    """Only Smoke->BMI and CRP->Statin are needed; Meek forces the other three.

    This redundancy is the cascade the report points at, so it is pinned here.
    """
    minimal = [("Smoke", "BMI"), PIVOT]
    assert apply_orientations(_cpdag(), minimal) == true_dag()


def test_all_scenarios_are_consistent_with_the_cpdag():
    """Every scenario passes the check an analyst would actually run."""
    for s in scenarios().values():
        assert apply_orientations(_cpdag(), s["knowledge"]) is not None


def test_scenarios_are_distinct():
    """A, B and C must give three different G0.

    Regression guard: the first draft of C omitted BMI->Chol, which Meek
    re-derives, making C identical to B.
    """
    cpdag = _cpdag()
    gs = {k: apply_orientations(cpdag, s["knowledge"]) for k, s in scenarios().items()}
    assert gs["A"] != gs["B"]
    assert gs["B"] != gs["C"]
    assert gs["A"] != gs["C"]


def test_scenario_b_is_false_but_consistent():
    """B asserts Statin->CRP; the truth says CRP->Statin. Consistent, and wrong."""
    truth = true_dag()
    b = scenarios()["B"]["knowledge"]
    assert apply_orientations(_cpdag(), b) is not None
    assert (PIVOT[1], PIVOT[0]) in b
    assert truth.is_directed_edge(*PIVOT)


def test_optimal_set_is_age_smoke_in_truth():
    """Hand-derived: cn={BP,CVD}, pa(cn)={Statin,Age,BP,Smoke}, so O*={Age,Smoke}."""
    assert optimal_adjustment_set_mpdag(true_dag(), TREATMENT, OUTCOME) == {"Age", "Smoke"}


def test_optimal_set_is_valid_in_truth():
    assert is_valid_adjustment_set_mpdag(true_dag(), TREATMENT, OUTCOME, {"Age", "Smoke"})


def test_empty_set_is_confounded():
    """There is genuine confounding, or the example would be vacuous."""
    assert not is_valid_adjustment_set_mpdag(true_dag(), TREATMENT, OUTCOME, set())


def test_pipeline_is_bit_reproducible_across_hash_seeds() -> None:
    """The whole per-element evaluation must be identical across PYTHONHASHSEED.

    Regression test for two separate set-iteration-order bugs, both of which
    made results stable within a process and silently different between runs:

    * ``random_sem`` drew coefficients while iterating a frozenset of edges, so
      the RNG stream itself differed (a large effect -- it changed which radius
      was reported);
    * ``LinearSEM.true_total_effect`` summed path products in set order, and
      floating-point addition is not associative (a last-bit effect).
    """
    import subprocess
    import sys

    code = (
        "from bkrobust.demo.pipeline import run_scenario;"
        "r = run_scenario('A', n_draws=25);"
        "print([(x['edge_string'], repr(x['mean_abs_bias'])) for x in r['rows']]);"
        "print(r['radii'])"
    )
    outs = set()
    for hashseed in ("0", "1", "12345"):
        env = {"PYTHONPATH": "src", "PYTHONHASHSEED": hashseed, "PATH": "/usr/bin:/bin"}
        res = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, env=env, check=True
        )
        outs.add(res.stdout.strip())
    assert len(outs) == 1, f"pipeline is not bit-reproducible: {len(outs)} distinct results"
