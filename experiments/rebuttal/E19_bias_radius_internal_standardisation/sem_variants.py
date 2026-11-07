"""Three parameterisations of one coefficient draw, for E19 (and the same idea as E18).

Starts from the submitted ``bkrobust.demo.evaluate.random_sem`` draw and returns
another ``LinearSEM`` over the same DAG:

* ``naive``    -- the draw unchanged.
* ``posthoc``  -- every variable divided by its marginal sd after the fact.
* ``internal`` -- iSCM (Ormaniec, Sussex, Lorch, Rothfuss, Krause, arXiv:2406.11601):
  in topological order, each variable is standardised before its children see it.

The internal construction follows the pilot constructor ``make_linear_iscm`` in
the pilot tree, rho_breakdown/pilot/code/scm.py
(commit 8486bb3509decf17adf490c078607b225d8df1b0;
last commit touching the file 69656dcceaaf4a75d3b2473f66338d5f1b256382). No code is copied:
that constructor estimates sd_j from a 100,000-row reference sample, while here
sd_j is the exact population value, because in the linear case the frozen iSCM is
again a linear SEM with coefficients w_ij / sd_j and noise variance sigma_j^2 / sd_j^2,
and a Monte-Carlo sd would add noise to a quantity known in closed form.
"""

from __future__ import annotations

import numpy as np

from bkrobust.demo.evaluate import LinearSEM

ARMS = ("naive", "posthoc", "internal")
UNIT_VAR_TOL = 1e-9


def topological_order(dag) -> list:
    # Kahn with sorted ties, so the order (and the float sums) never depend on set order.
    remaining = {v: set(dag.parents(v)) for v in dag.nodes}
    order: list = []
    while remaining:
        ready = sorted(v for v, pa in remaining.items() if not pa)
        if not ready:
            raise ValueError("graph is not acyclic")
        for v in ready:
            order.append(v)
            del remaining[v]
        for pa in remaining.values():
            pa.difference_update(ready)
    return order


def posthoc_standardise(sem: LinearSEM) -> LinearSEM:
    nodes = list(sem.dag.nodes)
    sd = dict(zip(nodes, np.sqrt(np.diag(sem.covariance())), strict=True))
    weights = {(p, v): w * sd[p] / sd[v] for (p, v), w in sem.weights.items()}
    noise = {v: sem.noise_var[v] / sd[v] ** 2 for v in nodes}
    return LinearSEM(dag=sem.dag, weights=weights, noise_var=noise)


def internal_standardise(sem: LinearSEM) -> LinearSEM:
    dag = sem.dag
    nodes = list(dag.nodes)
    idx = {v: i for i, v in enumerate(nodes)}
    cov = np.zeros((len(nodes), len(nodes)))  # covariance of the standardised variables so far
    done: list = []
    weights: dict = {}
    noise: dict = {}
    for v in topological_order(dag):
        pa = sorted(dag.parents(v))
        w = np.array([sem.weights[(p, v)] for p in pa])
        pi = [idx[p] for p in pa]
        raw_var = float(w @ cov[np.ix_(pi, pi)] @ w) + sem.noise_var[v] if pa else sem.noise_var[v]
        sd = np.sqrt(raw_var)
        for p, wp in zip(pa, w, strict=True):
            weights[(p, v)] = float(wp / sd)
        noise[v] = sem.noise_var[v] / raw_var
        a = w / sd
        for k in done:
            c = float(a @ cov[pi, idx[k]]) if pa else 0.0
            cov[idx[v], idx[k]] = cov[idx[k], idx[v]] = c
        cov[idx[v], idx[v]] = 1.0
        done.append(v)
    return LinearSEM(dag=dag, weights=weights, noise_var=noise)


def make_arm(sem: LinearSEM, arm: str) -> LinearSEM:
    if arm == "naive":
        return sem
    if arm == "posthoc":
        return posthoc_standardise(sem)
    if arm == "internal":
        return internal_standardise(sem)
    raise ValueError(f"unknown arm {arm!r}")


def max_unit_var_error(sem: LinearSEM) -> float:
    return float(np.max(np.abs(np.diag(sem.covariance()) - 1.0)))


def _path_counts(dag) -> tuple[list, np.ndarray]:
    # entry (i, j) = number of directed paths i -> ... -> j of length >= 1
    nodes = list(dag.nodes)
    idx = {v: i for i, v in enumerate(nodes)}
    a = np.zeros((len(nodes), len(nodes)))
    for t, h in dag.directed_edges:
        a[idx[t], idx[h]] = 1.0
    total = np.zeros_like(a)
    m = np.eye(len(nodes))
    for _ in range(len(nodes) - 1):
        m = m @ a
        total += m
    return nodes, total


def _sortability(score: np.ndarray, counts: np.ndarray) -> float:
    # Reisach et al. (arXiv:2303.18211): share of directed paths along which the score
    # increases, ties counted one half; population version, read off Sigma.
    tot = inc = 0.0
    tol = 1e-9 * max(1.0, float(np.max(np.abs(score))))
    for i, j in zip(*np.nonzero(counts), strict=True):
        c = counts[i, j]
        tot += c
        d = score[j] - score[i]
        inc += c if d > tol else (0.5 * c if abs(d) <= tol else 0.0)
    return 0.5 if tot == 0 else inc / tot


def sortabilities(sem: LinearSEM) -> dict[str, float]:
    sigma = sem.covariance()
    _, counts = _path_counts(sem.dag)
    var = np.diag(sigma)
    r2 = 1.0 - 1.0 / (var * np.diag(np.linalg.inv(sigma)))
    return {"var_sortability": _sortability(var, counts), "r2_sortability": _sortability(r2, counts)}
