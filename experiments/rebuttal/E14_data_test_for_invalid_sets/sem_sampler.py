"""Non-Gaussian linear SEM on a DAG, standardised to unit marginal variances.

The submitted ``bkrobust.demo.evaluate`` only has the linear-Gaussian SEM and no sampler. The
plan asked for a sampler there; it lives here instead so that src/ stays as submitted. Weights and
noise variances come from the submitted ``random_sem``; this module adds three things:

* unit marginal variances, by rescaling each node in topological order (a reparametrisation of the
  same SEM on the same DAG). Without it variances grow geometrically along the deep paths of the
  larger networks and the regressions become ill-conditioned;
* non-Gaussian noise (uniform, Laplace or centred exponential per node), because the test's
  covariance must not lean on Gaussian moments;
* the population covariance, so the population coefficient of every adjustment set is exact.
"""

from __future__ import annotations

import numpy as np

from bkrobust.demo.evaluate import LinearSEM, random_sem
from bkrobust.demo.graph import MPDAG, induced_subgraph

FAMILIES = ("uniform", "laplace", "exponential")


def ancestral_closure(dag: MPDAG, nodes) -> set:
    """``nodes`` and all their ancestors: the marginal of a linear SEM on it is the same SEM."""
    out = set(nodes)
    for v in list(nodes):
        out |= dag.ancestors(v)
    return out


def topological_order(dag: MPDAG) -> list:
    """Kahn's algorithm with sorted ties, so the order (and the RNG stream) is deterministic."""
    indeg = {v: len(dag.parents(v)) for v in dag.nodes}
    ready = sorted(v for v, d in indeg.items() if d == 0)
    order = []
    while ready:
        v = ready.pop(0)
        order.append(v)
        for c in sorted(dag.children(v)):
            indeg[c] -= 1
            if indeg[c] == 0:
                ready.append(c)
        ready.sort()
    if len(order) != len(dag.nodes):
        raise ValueError("graph is not acyclic")
    return order


def draw_sem(dag: MPDAG, keep, rng: np.random.Generator) -> dict:
    """A standardised non-Gaussian linear SEM on ``dag`` restricted to the ancestral set ``keep``."""
    sub = induced_subgraph(dag, keep)
    base = random_sem(sub, rng)
    order = topological_order(sub)
    idx = {v: i for i, v in enumerate(order)}
    p = len(order)
    sigma = np.zeros((p, p))
    parents, weights, noise_sd = [], [], np.zeros(p)
    new_w: dict = {}
    new_var: dict = {}
    for j, v in enumerate(order):
        pa = sorted(sub.parents(v))
        pi = np.array([idx[u] for u in pa], dtype=int)
        w = np.array([base.weights[(u, v)] for u in pa], dtype=float)
        omega = base.noise_var[v]
        var = float(w @ sigma[np.ix_(pi, pi)] @ w) + omega if len(pi) else omega
        s = 1.0 / np.sqrt(var)
        w, omega = w * s, omega * s * s
        if len(pi):
            sigma[j, :j] = w @ sigma[pi, :j]
            sigma[:j, j] = sigma[j, :j]
        sigma[j, j] = 1.0
        parents.append(pi)
        weights.append(w)
        noise_sd[j] = np.sqrt(omega)
        for u, wu in zip(pa, w):
            new_w[(u, v)] = float(wu)
        new_var[v] = float(omega)
    family = rng.integers(len(FAMILIES), size=p)
    sem = LinearSEM(dag=sub, weights=new_w, noise_var=new_var)
    return {"order": order, "idx": idx, "parents": parents, "weights": weights,
            "noise_sd": noise_sd, "family": family, "sigma": sigma, "sem": sem}


def check_covariance(model: dict) -> float:
    """Largest gap between the incremental covariance and the submitted closed form."""
    sem = model["sem"]
    closed = sem.covariance()  # in sorted node order
    nodes = list(sem.dag.nodes)
    perm = [model["idx"][v] for v in nodes]
    return float(np.max(np.abs(model["sigma"][np.ix_(perm, perm)] - closed)))


def unit_noise(kind: int, n: int, rng: np.random.Generator) -> np.ndarray:
    """Mean 0, variance 1, non-Gaussian."""
    if FAMILIES[kind] == "uniform":
        return rng.uniform(-np.sqrt(3.0), np.sqrt(3.0), size=n)
    if FAMILIES[kind] == "laplace":
        return rng.laplace(0.0, 1.0 / np.sqrt(2.0), size=n)
    return rng.exponential(1.0, size=n) - 1.0


def sample(model: dict, n: int, rng: np.random.Generator) -> np.ndarray:
    """``n`` draws, columns in ``model['order']``."""
    p = len(model["order"])
    data = np.empty((n, p))
    for j in range(p):
        e = model["noise_sd"][j] * unit_noise(int(model["family"][j]), n, rng)
        pi = model["parents"][j]
        data[:, j] = data[:, pi] @ model["weights"][j] + e if len(pi) else e
    return data


def population_coefficient(model: dict, x, y, z) -> float:
    """Coefficient of ``x`` in the population regression of ``y`` on ``x`` and ``z``."""
    idx, sigma = model["idx"], model["sigma"]
    reg = [idx[x]] + [idx[v] for v in sorted(z)]
    return float(np.linalg.solve(sigma[np.ix_(reg, reg)], sigma[reg, idx[y]])[0])
