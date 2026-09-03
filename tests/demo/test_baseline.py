"""Tests for :mod:`bkrobust.demo.baseline` -- naive-distance baseline and calibration.

Two things are guarded here: that the naive Hamming metric is actually a
metric (so comparing it to the model-oriented BFS distance is a fair
comparison, not an apples-to-oranges one), and that the calibration machinery
that checks the certificate against a hypothetical truth is itself correct --
in particular that coverage really is 1.0 (a bug otherwise, per the module's
own framing) and that the UNREACHED sentinel (-1) is never silently averaged
into an aggregate as if it were a real, small distance.
"""

from __future__ import annotations

import itertools
import random

import pandas as pd
import pytest

from bkrobust.demo.baseline import (
    _agreement_from_frame,
    calibration_summary,
    compare_radii,
    disagreement_examples,
    inconsistent_fraction,
    k_sets_for_space,
    naive_distance,
    naive_radius,
    naive_vs_model_frame,
)
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.pipeline import UNREACHED, run_scenario
from bkrobust.demo.scenario import PIVOT, scenarios, true_dag

LABELS = ("A", "B", "C")


def _cpdag() -> MPDAG:
    return dag_to_cpdag(true_dag())


# --- naive_distance is a metric -----------------------------------------------


def _sample_k_sets(pool_size: int = 20, max_edge_pairs: int = 3, seed: int = 7) -> list:
    """A varied sample of small edge-assertion sets drawn from the CPDAG's own edges."""
    cpdag = _cpdag()
    undirected = sorted(cpdag.undirected_edges)
    universe = []
    for a, b in undirected:
        universe.append((a, b))
        universe.append((b, a))

    rng = random.Random(seed)
    samples = []
    for _ in range(pool_size):
        size = rng.randint(0, max_edge_pairs)
        samples.append(frozenset(rng.sample(universe, size)))
    # Always include the empty set and the full K_TRUE-derived sets for scenarios.
    samples.append(frozenset())
    for spec in scenarios().values():
        samples.append(frozenset(spec["knowledge"]))
    return samples


def test_naive_distance_identity():
    for k in _sample_k_sets():
        assert naive_distance(k, k) == 0


def test_naive_distance_symmetry():
    samples = _sample_k_sets()
    for k1, k2 in itertools.islice(itertools.product(samples, samples), 400):
        assert naive_distance(k1, k2) == naive_distance(k2, k1)


def test_naive_distance_zero_implies_equal():
    samples = _sample_k_sets()
    for k1, k2 in itertools.product(samples, samples):
        if naive_distance(k1, k2) == 0:
            assert set(k1) == set(k2)


def test_naive_distance_triangle_inequality():
    samples = _sample_k_sets()
    rng = random.Random(11)
    triples = [tuple(rng.sample(samples, 3)) for _ in range(300)]
    for k1, k2, k3 in triples:
        d13 = naive_distance(k1, k3)
        d12 = naive_distance(k1, k2)
        d23 = naive_distance(k2, k3)
        assert d13 <= d12 + d23


def test_flip_has_naive_distance_two():
    """A flip costs 2 in the naive metric, mirroring the model side's two-move cost.

    ``atomic_moves`` (bkrobust.demo.space) charges a flip as a retraction plus
    a fresh assertion -- two atomic moves, never one, because an already
    oriented edge must be un-oriented before it can be asserted the other way.
    The naive Hamming metric on K agrees for exactly the same reason: flipping
    one assertion removes one element and adds a different one, so the
    symmetric difference has size exactly 2.
    """
    k_true = scenarios()["B"]["knowledge"]
    # PIVOT was already flipped from K_true to build scenario B's knowledge; flip
    # it back to get a concrete before/after pair.
    flipped_pivot = (PIVOT[1], PIVOT[0])
    assert flipped_pivot in k_true
    original = [e for e in k_true if e != flipped_pivot] + [PIVOT]
    assert naive_distance(k_true, original) == 2


# --- inconsistent_fraction ----------------------------------------------------


@pytest.mark.parametrize("label", LABELS)
def test_inconsistent_fraction_counts_and_range(label):
    cpdag = _cpdag()
    k_assumed = scenarios()[label]["knowledge"]
    result = inconsistent_fraction(cpdag, k_assumed)

    total_from_counts = sum(c["total"] for c in result["counts"].values())
    inconsistent_from_counts = sum(c["inconsistent"] for c in result["counts"].values())

    assert total_from_counts == result["total_perturbations"]
    assert inconsistent_from_counts == result["total_inconsistent"]
    assert 0 <= result["fraction_inconsistent"] <= 1
    assert result["total_perturbations"] > 0
    # Sanity: recompute "remove" count independently.
    assert result["counts"]["remove"]["total"] == len(k_assumed)


def test_inconsistent_fraction_matches_manual_check_on_scenario_a():
    """Cross-check inconsistent_fraction's bookkeeping against a hand-rolled enumeration."""
    cpdag = _cpdag()
    k_assumed = scenarios()["A"]["knowledge"]
    result = inconsistent_fraction(cpdag, k_assumed)

    k_set = set(k_assumed)
    manual_bad = 0
    manual_total = 0
    for e in k_set:
        manual_total += 1
        if apply_orientations(cpdag, sorted(k_set - {e})) is None:
            manual_bad += 1
    for tail, head in k_set:
        manual_total += 1
        flipped = (k_set - {(tail, head)}) | {(head, tail)}
        if apply_orientations(cpdag, sorted(flipped)) is None:
            manual_bad += 1
    from bkrobust.demo.graph import canon

    asserted_canon = {canon(a, b) for a, b in k_set}
    for a, b in cpdag.undirected_edges:
        if (a, b) in asserted_canon:
            continue
        for cand in ((a, b), (b, a)):
            manual_total += 1
            if apply_orientations(cpdag, sorted(k_set | {cand})) is None:
                manual_bad += 1

    assert result["total_perturbations"] == manual_total
    assert result["total_inconsistent"] == manual_bad


# --- naive_vs_model_frame / k_sets_for_space -----------------------------------


@pytest.mark.parametrize("label", LABELS)
def test_naive_vs_model_frame_shape(label):
    result = run_scenario(label, n_draws=5)
    frame = naive_vs_model_frame(result)
    assert list(frame.columns) == [
        "edge_string",
        "model_distance",
        "naive_distance",
        "z_valid",
        "z_is_optimal",
        "mean_abs_bias",
    ]
    assert len(frame) == len(result["space"])
    assert set(frame["edge_string"]) == {g.edge_string() for g in result["space"]}


@pytest.mark.parametrize("label", LABELS)
def test_k_sets_for_space_g0_is_zero(label):
    """G0 itself is reachable from k_assumed at naive distance 0 (assert exactly K)."""
    result = run_scenario(label, n_draws=5)
    k_sets = k_sets_for_space(result["cpdag"], result["spec"]["knowledge"], result["space"])
    assert k_sets[result["g0"]] == 0


@pytest.mark.parametrize("label", LABELS)
def test_k_sets_for_space_no_negative_other_than_unreached(label):
    result = run_scenario(label, n_draws=5)
    k_sets = k_sets_for_space(result["cpdag"], result["spec"]["knowledge"], result["space"])
    for d in k_sets.values():
        assert d == UNREACHED or d >= 0


# --- calibration ---------------------------------------------------------------


@pytest.mark.parametrize("label", LABELS)
def test_coverage_is_one(label):
    """Coverage must be 1.0 by construction: r_val is the first shell with a failure.

    If this fails it is a genuine bug in the radius computation or in the
    calibration check itself -- not a robustness finding -- per the module's
    own documented framing. Do not weaken this test to accommodate a failure;
    investigate instead.
    """
    result = run_scenario(label, n_draws=5)
    summary = calibration_summary(result)
    assert summary["counts"]["n_total"] == len(result["space"])
    assert not summary["coverage_is_bug"]
    if summary["counts"]["inside_valid"] + summary["counts"]["inside_invalid"] > 0:
        assert summary["coverage"] == pytest.approx(1.0)


@pytest.mark.parametrize("label", LABELS)
def test_compare_radii_smoke(label):
    result = run_scenario(label, n_draws=5)
    comparison = compare_radii(result)
    assert comparison["model_radius"] == result["radii"]["r_val"]
    assert isinstance(comparison["radii_agree"], bool)
    assert comparison["n_compared"] >= 0
    assert comparison["n_disagree_by_2_or_more"] >= 0


@pytest.mark.parametrize("label", LABELS)
def test_disagreement_examples_shape(label):
    result = run_scenario(label, n_draws=5)
    examples = disagreement_examples(result, n=3)
    assert len(examples) <= 3
    for ex in examples:
        assert set(ex) == {
            "edge_string",
            "model_distance",
            "naive_distance",
            "difference",
            "atomic_moves",
            "description",
        }
        assert ex["model_distance"] != UNREACHED
        assert ex["naive_distance"] != UNREACHED
        assert isinstance(ex["description"], str) and len(ex["description"]) > 0
    # Sorted most-disagreeing first.
    diffs = [ex["difference"] for ex in examples]
    assert diffs == sorted(diffs, reverse=True)


# --- UNREACHED is never silently treated as a real distance -------------------


def test_agreement_from_frame_excludes_unreached_rows():
    """A row carrying UNREACHED in either distance column must not enter any aggregate.

    Built directly against the frame-shaped helper (rather than a full
    ``run_scenario`` result) so the exclusion can be checked even if neither
    real scenario happens to produce an UNREACHED naive distance: two clean,
    perfectly correlated rows are given a real signal, and a third,
    wildly-inconsistent-looking row is marked UNREACHED in one column. If
    UNREACHED (-1) were folded into the correlation or the disagreement count
    as an ordinary small integer, the results below would differ from the
    "compute only on the two clean rows" ground truth.
    """
    clean = pd.DataFrame(
        {
            "model_distance": [0, 1, 2, 3],
            "naive_distance": [0, 1, 2, 3],
            "z_valid": [True, True, False, False],
        }
    )
    contaminated = pd.concat(
        [
            clean,
            pd.DataFrame(
                {
                    "model_distance": [UNREACHED, 5],
                    "naive_distance": [7, UNREACHED],
                    "z_valid": [False, False],
                }
            ),
        ],
        ignore_index=True,
    )

    clean_stats = _agreement_from_frame(clean)
    contaminated_stats = _agreement_from_frame(contaminated)

    assert contaminated_stats["n_compared"] == clean_stats["n_compared"] == 4
    assert contaminated_stats["pearson"] == pytest.approx(clean_stats["pearson"])
    assert contaminated_stats["spearman"] == pytest.approx(clean_stats["spearman"])
    assert contaminated_stats["n_disagree_by_2_or_more"] == 0
    assert clean_stats["n_disagree_by_2_or_more"] == 0

    # Sanity: if UNREACHED were treated as a real value, the two UNREACHED rows
    # would register a large |model_distance - naive_distance| (e.g. |-1 - 7| = 8)
    # and inflate n_disagree_by_2_or_more, which the assertion above rules out.


def test_naive_radius_ignores_unreached_naive_distance():
    """naive_radius must not report UNREACHED (-1) as if it were the smallest real distance."""
    for label in LABELS:
        result = run_scenario(label, n_draws=5)
        frame = naive_vs_model_frame(result)
        r = naive_radius(result)
        if r == UNREACHED:
            # No invalid element with a defined naive distance exists at all.
            invalid = frame[~frame["z_valid"]]
            assert (invalid["naive_distance"] == UNREACHED).all() or invalid.empty
        else:
            assert r >= 0
            invalid_defined = frame[(~frame["z_valid"]) & (frame["naive_distance"] != UNREACHED)]
            assert r == invalid_defined["naive_distance"].min()
