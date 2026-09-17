"""The sortability functional against the worked example of Reisach et al. (2021).

Section 3.1 of that paper computes var-sortability 3/4 on ``A -> B``, ``A -> C``,
``B -> C`` with ``Var(A) = 2``, ``Var(B) = 1``, ``Var(C) = 3``. The probe must
reproduce it, count a tie as one half, count one term per path length rather
than per path, and, on a large ancestral sample, agree with the population value
computed from the exact covariance.
"""

from __future__ import annotations

import numpy as np
from experiments.sortability_probe import (
    GaussianNetwork,
    order_alignment,
    population_r2,
    r2_on_all_others,
)


def _adjacency(nodes: list[str], arcs: list[tuple[str, str]]) -> np.ndarray:
    index = {v: i for i, v in enumerate(nodes)}
    adj = np.zeros((len(nodes), len(nodes)), dtype=bool)
    for s, t in arcs:
        adj[index[s], index[t]] = True
    return adj


def test_worked_example_is_three_quarters():
    adj = _adjacency(["A", "B", "C"], [("A", "B"), ("A", "C"), ("B", "C")])
    value, terms = order_alignment(adj, np.array([2.0, 1.0, 3.0]))
    assert terms == 4
    assert value == 0.75


def test_a_tie_counts_one_half():
    adj = _adjacency(["A", "B"], [("A", "B")])
    value, terms = order_alignment(adj, np.array([1.0, 1.0]))
    assert (value, terms) == (0.5, 1)


def test_one_term_per_path_length_not_per_path():
    diamond = _adjacency(["A", "B", "C", "D"], [("A", "B"), ("A", "C"), ("B", "D"), ("C", "D")])
    _, terms = order_alignment(diamond, np.arange(4.0))
    # Four arcs at length one, and (A, D) once at length two although two paths reach it.
    assert terms == 5


def test_no_arcs_is_undefined():
    value, terms = order_alignment(np.zeros((3, 3), dtype=bool), np.arange(3.0))
    assert terms == 0 and np.isnan(value)


def test_sample_matches_population_on_a_fitted_network():
    net = GaussianNetwork(
        name="toy",
        nodes=["A", "B", "C", "D"],
        arcs=[("A", "B"), ("A", "C"), ("B", "D"), ("C", "D")],
        intercept={"A": 0.3, "B": -1.0, "C": 0.0, "D": 2.0},
        coef={"A": {}, "B": {"A": 0.8}, "C": {"A": -0.5}, "D": {"B": 0.6, "C": 1.2}},
        noise_var={"A": 1.0, "B": 0.5, "C": 2.0, "D": 0.25},
    )
    sigma = net.population_covariance()
    x = net.sample(200_000, np.random.default_rng(0))
    assert np.allclose(x.var(axis=0, ddof=1), np.diag(sigma), rtol=0.02)
    assert np.allclose(r2_on_all_others(x), population_r2(sigma), atol=0.01)
    adj = net.adjacency
    assert order_alignment(adj, x.var(axis=0))[0] == order_alignment(adj, np.diag(sigma))[0]
