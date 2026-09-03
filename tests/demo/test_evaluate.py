"""Tests for :mod:`bkrobust.demo.evaluate`.

Each test is chosen to discriminate a specific failure mode: a confounder DAG
that must be adjusted for, a mediator that must *not* be adjusted for, a
collider that must not be conditioned on, an independent cross-check of
d-separation against networkx, a hand-worked optimal-adjustment-set example,
and closed-form SEM identities that are exact by construction (population
moments, not samples).

MPDAG-level tests (`is_valid_adjustment_set_mpdag`, `all_valid_adjustment_sets_mpdag`,
`optimal_adjustment_set_mpdag`) depend on `bkrobust.demo.meek.enumerate_dag_extensions`,
implemented by a concurrently-developed sibling module. If it is not yet
available, those specific tests are skipped rather than failing the whole
suite -- everything else in this file is independent of it.
"""

from __future__ import annotations

import networkx as nx
import numpy as np
import pytest

from bkrobust.demo.evaluate import (
    LinearSEM,
    adjusted_estimand,
    all_valid_adjustment_sets_mpdag,
    asymptotic_variance,
    bias,
    is_dseparated,
    is_valid_adjustment_set_dag,
    is_valid_adjustment_set_mpdag,
    optimal_adjustment_set_dag,
    optimal_adjustment_set_mpdag,
    random_sem,
)
from bkrobust.demo.graph import MPDAG

try:
    from bkrobust.demo.meek import enumerate_dag_extensions  # noqa: F401

    _HAVE_MEEK = True
except ImportError:
    _HAVE_MEEK = False

requires_meek = pytest.mark.skipif(
    not _HAVE_MEEK, reason="bkrobust.demo.meek.enumerate_dag_extensions not yet available"
)


# ----------------------------------------------------------------------------
# Fixture graphs
# ----------------------------------------------------------------------------


def confounder_dag() -> MPDAG:
    """Z -> X, Z -> Y, X -> Y: the classic confounder."""
    return MPDAG(["X", "Y", "Z"], directed=[("Z", "X"), ("Z", "Y"), ("X", "Y")])


def mediator_dag() -> MPDAG:
    """X -> M -> Y."""
    return MPDAG(["X", "M", "Y"], directed=[("X", "M"), ("M", "Y")])


def collider_dag() -> MPDAG:
    """X -> C <- Y, plus X -> Y."""
    return MPDAG(["X", "Y", "C"], directed=[("X", "C"), ("Y", "C"), ("X", "Y")])


def hpm_example_dag() -> MPDAG:
    """Z1 -> X, Z1 -> Y, Z2 -> Y, X -> M -> Y, X -> Y.

    The hand-worked optimal-adjustment-set example from Henckel, Perkovic &
    Maathuis: cn(X, Y) = {M, Y}, pa({M, Y}) = {X, Z1, Z2, M}, so
    O = {Z1, Z2}.
    """
    return MPDAG(
        ["X", "Y", "M", "Z1", "Z2"],
        directed=[
            ("Z1", "X"),
            ("Z1", "Y"),
            ("Z2", "Y"),
            ("X", "M"),
            ("M", "Y"),
            ("X", "Y"),
        ],
    )


# ----------------------------------------------------------------------------
# Back-door criterion
# ----------------------------------------------------------------------------


def test_confounder_empty_set_invalid_full_set_valid() -> None:
    dag = confounder_dag()
    assert not is_valid_adjustment_set_dag(dag, "X", "Y", set())
    assert is_valid_adjustment_set_dag(dag, "X", "Y", {"Z"})


def test_mediator_adjusting_for_mediator_is_invalid() -> None:
    dag = mediator_dag()
    assert not is_valid_adjustment_set_dag(dag, "X", "Y", {"M"})
    assert is_valid_adjustment_set_dag(dag, "X", "Y", set())


def test_collider_conditioning_on_collider_is_invalid() -> None:
    dag = collider_dag()
    # C is a collider on the X -> C <- Y path; it is not on any back-door path
    # from X to Y at all (its only path-role is as a collider), so conditioning
    # on it *opens* a spurious association without blocking anything needed.
    assert not is_valid_adjustment_set_dag(dag, "X", "Y", {"C"})
    assert is_valid_adjustment_set_dag(dag, "X", "Y", set())


def test_descendant_of_x_exclusion_enforced() -> None:
    # D is a descendant of X but not on the X -> Y path at all: condition (a)
    # must reject {D} regardless, distinct from the mediator case above where
    # the descendant is also the thing blocking the effect.
    dag = MPDAG(["X", "Y", "D"], directed=[("X", "Y"), ("X", "D")])
    assert is_valid_adjustment_set_dag(dag, "X", "Y", set())
    assert not is_valid_adjustment_set_dag(dag, "X", "Y", {"D"})


def test_x_or_y_in_z_is_invalid() -> None:
    dag = confounder_dag()
    assert not is_valid_adjustment_set_dag(dag, "X", "Y", {"X"})
    assert not is_valid_adjustment_set_dag(dag, "X", "Y", {"Y"})


# ----------------------------------------------------------------------------
# d-separation cross-check against networkx (independent second implementation)
# ----------------------------------------------------------------------------


def _random_dag(rng: np.random.Generator, n: int, edge_prob: float) -> MPDAG:
    """A random DAG on n nodes: edges only go from lower to higher topological index."""
    nodes = [f"n{i}" for i in range(n)]
    order = list(nodes)
    rng.shuffle(order)
    edges: list[tuple[str, str]] = []
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < edge_prob:
                edges.append((order[i], order[j]))
    return MPDAG(nodes, directed=edges)


def _to_networkx(dag: MPDAG) -> nx.DiGraph:
    g = nx.DiGraph()
    g.add_nodes_from(dag.nodes)
    g.add_edges_from(dag.directed_edges)
    return g


def test_dseparation_matches_networkx_on_random_dags() -> None:
    rng = np.random.default_rng(20260903)
    n_trials = 220
    checked = 0
    for _ in range(n_trials):
        n = int(rng.integers(3, 7))
        edge_prob = float(rng.uniform(0.2, 0.6))
        dag = _random_dag(rng, n, edge_prob)
        nxg = _to_networkx(dag)

        nodes = list(dag.nodes)
        a, b = rng.choice(nodes, size=2, replace=False)
        remaining = [x for x in nodes if x not in (a, b)]
        rng.shuffle(remaining)
        z_size = int(rng.integers(0, len(remaining) + 1))
        z = set(remaining[:z_size])

        expected = nx.algorithms.d_separated(nxg, {a}, {b}, z)
        actual = is_dseparated(dag, a, b, z)
        assert actual == expected, f"mismatch for dag={dag}, a={a}, b={b}, z={z}"
        checked += 1
    assert checked >= 200


# ----------------------------------------------------------------------------
# Optimal adjustment set
# ----------------------------------------------------------------------------


def test_optimal_adjustment_set_hand_worked_example() -> None:
    dag = hpm_example_dag()
    optimal = optimal_adjustment_set_dag(dag, "X", "Y")
    assert optimal == {"Z1", "Z2"}

    # {Z1} alone is a valid (but not optimal, and not equal to O) adjustment set:
    # it blocks the X <- Z1 -> Y back-door path, and Z1 is not a descendant of X.
    assert is_valid_adjustment_set_dag(dag, "X", "Y", {"Z1"})
    assert frozenset({"Z1"}) != frozenset(optimal)


def test_optimal_adjustment_set_no_causal_path_is_empty() -> None:
    # X and Y have no directed path at all; cn(X, Y) is empty, so O is empty.
    dag = MPDAG(["X", "Y", "Z"], directed=[("Z", "X"), ("Z", "Y")])
    assert optimal_adjustment_set_dag(dag, "X", "Y") == set()


# ----------------------------------------------------------------------------
# MPDAG-level lifting (requires bkrobust.demo.meek)
# ----------------------------------------------------------------------------


@requires_meek
def test_mpdag_validity_requires_agreement_across_extensions() -> None:
    # Fully-oriented confounder DAG has exactly one "extension" (itself).
    dag = confounder_dag()
    assert is_valid_adjustment_set_mpdag(dag, "X", "Y", {"Z"})
    assert not is_valid_adjustment_set_mpdag(dag, "X", "Y", set())


@requires_meek
def test_all_valid_adjustment_sets_mpdag_matches_dag_level() -> None:
    dag = confounder_dag()
    valid = all_valid_adjustment_sets_mpdag(dag, "X", "Y")
    assert frozenset({"Z"}) in valid
    assert frozenset() not in valid
    # Every returned set actually satisfies the DAG-level criterion too.
    for z in valid:
        assert is_valid_adjustment_set_dag(dag, "X", "Y", z)
    # Sorted by (size, sorted labels).
    sizes = [len(z) for z in valid]
    assert sizes == sorted(sizes)


@requires_meek
def test_optimal_adjustment_set_mpdag_agrees_on_fully_oriented_dag() -> None:
    dag = hpm_example_dag()
    assert optimal_adjustment_set_mpdag(dag, "X", "Y") == {"Z1", "Z2"}


# ----------------------------------------------------------------------------
# LinearSEM: closed-form correctness
# ----------------------------------------------------------------------------


def test_sem_two_node_hand_checkable_covariance() -> None:
    # X -> Y, weight w, Var(X) = sigma_x^2, Var(Y) = w^2 sigma_x^2 + sigma_y^2,
    # Cov(X, Y) = w * sigma_x^2. Verifies the (I - B)^-1 Omega (I - B)^-T orientation.
    dag = MPDAG(["X", "Y"], directed=[("X", "Y")])
    w = 2.0
    sigma_x2 = 1.5
    sigma_y2 = 0.7
    sem = LinearSEM(dag=dag, weights={("X", "Y"): w}, noise_var={"X": sigma_x2, "Y": sigma_y2})
    sigma = sem.covariance()
    nodes = list(dag.nodes)  # sorted: ["X", "Y"]
    ix, iy = nodes.index("X"), nodes.index("Y")
    assert sigma[ix, ix] == pytest.approx(sigma_x2)
    assert sigma[iy, iy] == pytest.approx(w**2 * sigma_x2 + sigma_y2)
    assert sigma[ix, iy] == pytest.approx(w * sigma_x2)
    assert sigma[iy, ix] == pytest.approx(w * sigma_x2)


def test_sem_simple_edge_true_effect_and_zero_bias() -> None:
    dag = MPDAG(["X", "Y"], directed=[("X", "Y")])
    sem = LinearSEM(dag=dag, weights={("X", "Y"): 2.0}, noise_var={"X": 1.0, "Y": 1.0})
    assert sem.true_total_effect("X", "Y") == pytest.approx(2.0)
    assert adjusted_estimand(sem, "X", "Y", set()) == pytest.approx(2.0)
    assert bias(sem, "X", "Y", set()) == pytest.approx(0.0, abs=1e-10)


def test_sem_confounding_bias_nonzero_unadjusted_zero_adjusted() -> None:
    dag = confounder_dag()
    sem = LinearSEM(
        dag=dag,
        weights={("Z", "X"): 1.0, ("Z", "Y"): 1.0, ("X", "Y"): 2.0},
        noise_var={"Z": 1.0, "X": 1.0, "Y": 1.0},
    )
    assert sem.true_total_effect("X", "Y") == pytest.approx(2.0)
    unadjusted_bias = bias(sem, "X", "Y", set())
    adjusted_bias = bias(sem, "X", "Y", {"Z"})
    assert abs(unadjusted_bias) > 1e-3
    assert adjusted_bias == pytest.approx(0.0, abs=1e-10)


def test_sem_multi_path_true_total_effect() -> None:
    # X -> M -> Y and X -> Y directly: total effect is the sum over both paths.
    dag = mediator_dag()
    dag = MPDAG(["X", "M", "Y"], directed=[("X", "M"), ("M", "Y"), ("X", "Y")])
    sem = LinearSEM(
        dag=dag,
        weights={("X", "M"): 2.0, ("M", "Y"): 3.0, ("X", "Y"): 0.5},
        noise_var={"X": 1.0, "M": 1.0, "Y": 1.0},
    )
    # Direct path X->Y (0.5) plus indirect X->M->Y (2.0 * 3.0 = 6.0) = 6.5.
    assert sem.true_total_effect("X", "Y") == pytest.approx(6.5)


def test_covariance_symmetric_positive_definite() -> None:
    dag = hpm_example_dag()
    rng = np.random.default_rng(7)
    for _ in range(25):
        sem = random_sem(dag, rng)
        sigma = sem.covariance()
        assert np.allclose(sigma, sigma.T)
        # Raises LinAlgError if not positive definite.
        np.linalg.cholesky(sigma)


def test_random_sem_uses_only_passed_generator(monkeypatch: pytest.MonkeyPatch) -> None:
    # random_sem must not touch the global numpy RNG state.
    dag = confounder_dag()
    np.random.seed(12345)
    state_before = np.random.get_state()[1].copy()
    rng = np.random.default_rng(0)
    random_sem(dag, rng)
    state_after = np.random.get_state()[1]
    assert np.array_equal(state_before, state_after)


def test_asymptotic_variance_is_positive_and_finite() -> None:
    dag = confounder_dag()
    rng = np.random.default_rng(3)
    sem = random_sem(dag, rng)
    v_empty = asymptotic_variance(sem, "X", "Y", set())
    v_adjusted = asymptotic_variance(sem, "X", "Y", {"Z"})
    for v in (v_empty, v_adjusted):
        assert np.isfinite(v)
        assert v > 0.0


def test_valid_adjustment_gives_near_zero_mean_bias_invalid_does_not() -> None:
    # A single invalid draw's bias can accidentally be ~0 by coefficient
    # cancellation, so this asserts the *mean absolute* bias over many draws,
    # not any single draw's bias.
    dag = confounder_dag()
    rng = np.random.default_rng(99)
    n_draws = 80
    valid_biases = []
    invalid_biases = []
    for _ in range(n_draws):
        sem = random_sem(dag, rng)
        valid_biases.append(bias(sem, "X", "Y", {"Z"}))
        invalid_biases.append(bias(sem, "X", "Y", set()))

    mean_abs_valid = float(np.mean(np.abs(valid_biases)))
    mean_abs_invalid = float(np.mean(np.abs(invalid_biases)))

    assert mean_abs_valid < 1e-8
    assert mean_abs_invalid > 0.1


def test_random_sem_is_reproducible_across_hash_seeds() -> None:
    """The same seed must give the same SEM regardless of PYTHONHASHSEED.

    Regression test. `random_sem` originally drew coefficients while iterating
    `dag.directed_edges`, a frozenset, so the draw order -- and therefore every
    coefficient -- depended on the per-process string hash seed. Results were
    stable within a process and silently differed between processes.
    """
    import subprocess
    import sys

    code = (
        "import numpy as np;"
        "from bkrobust.demo.scenario import true_dag;"
        "from bkrobust.demo.evaluate import random_sem;"
        "s=random_sem(true_dag(), np.random.default_rng(0));"
        "print(sorted((f'{a}->{b}', round(w, 12)) for (a, b), w in s.weights.items()))"
    )
    outs = set()
    for hashseed in ("0", "1", "12345"):
        env = {"PYTHONPATH": "src", "PYTHONHASHSEED": hashseed, "PATH": "/usr/bin:/bin"}
        out = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, env=env, check=True
        )
        outs.add(out.stdout.strip())
    assert len(outs) == 1, f"random_sem differs across PYTHONHASHSEED: {len(outs)} distinct results"
