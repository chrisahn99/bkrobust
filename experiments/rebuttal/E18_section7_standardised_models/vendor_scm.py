"""Internal standardisation and sortability, adapted for population covariances.

Source: the pilot tree, rho_breakdown/pilot/code/scm.py
        (make_linear_iscm, var_sortability, r2_sortability)
Source commit: commit 8486bb3509decf17adf490c078607b225d8df1b0 (file clean at that commit)

What changed, and why:
- make_linear_iscm fixed its standardisation constants on a 100,000-row Monte-Carlo
  reference sample. Here the constants are the exact population sds, computed in
  topological order, because Section 7 works on population covariances and a
  sampled constant would add noise that is not part of the model (Ormaniec et al.
  2025 define iSCM by the population moments; the sample was an approximation).
- The weights are not redrawn: the function takes the naive LinearSEM of the
  submitted code and standardises it, so the naive and internal arms of one draw
  share every coefficient and noise variance.
- The frozen iSCM is again a linear SEM with coefficients w_pj / sd_j and noise
  variances s_j^2 / sd_j^2 (scm.py's own docstring), returned as a LinearSEM.
- Sortability is computed from the population covariance: variances from its
  diagonal, R^2 of each node on all others from its inverse. Paths are counted with
  multiplicity, as in scm.py, but with matrix powers instead of a Python triple loop.
"""

from __future__ import annotations

import numpy as np

from bkrobust.demo.evaluate import LinearSEM


def topological_order(dag) -> list[str]:
    """Kahn's algorithm, ties broken by name so the order is reproducible."""
    indeg = {n: len(dag.parents(n)) for n in dag.nodes}
    ready = sorted(n for n, d in indeg.items() if d == 0)
    out: list[str] = []
    while ready:
        n = ready.pop(0)
        out.append(n)
        for c in sorted(dag.children(n)):
            indeg[c] -= 1
            if indeg[c] == 0:
                ready.append(c)
        ready.sort()
    if len(out) != len(dag.nodes):
        raise ValueError("graph is not acyclic")
    return out


def iscm_population(sem: LinearSEM) -> tuple[LinearSEM, dict[str, float]]:
    """The internally standardised version of ``sem`` and the sd used at each node.

    Each node's raw value sum_p w_pj x_p + e_j is divided by its population sd,
    where the parents x_p are already standardised. Intercepts do not affect
    covariances or effects and are dropped, as in scm.py.
    """
    dag = sem.dag
    order = topological_order(dag)
    idx = {n: i for i, n in enumerate(order)}
    p = len(order)
    sigma = np.zeros((p, p))  # covariance of the standardised variables, built row by row
    sd: dict[str, float] = {}
    weights: dict[tuple[str, str], float] = {}
    noise: dict[str, float] = {}
    for j in order:
        pa = sorted(dag.parents(j))
        w = np.array([sem.weights[(q, j)] for q in pa])
        pidx = [idx[q] for q in pa]
        if pa:
            var_raw = float(w @ sigma[np.ix_(pidx, pidx)] @ w) + sem.noise_var[j]
            cov_raw = w @ sigma[pidx, :]  # with every earlier node (later rows are still 0)
        else:
            var_raw = sem.noise_var[j]
            cov_raw = np.zeros(p)
        s = float(np.sqrt(var_raw))
        sd[j] = s
        row = cov_raw / s
        row[idx[j]] = 1.0
        sigma[idx[j], :] = row
        sigma[:, idx[j]] = row
        for q in pa:
            weights[(q, j)] = sem.weights[(q, j)] / s
        noise[j] = sem.noise_var[j] / s ** 2
    return LinearSEM(dag=dag, weights=weights, noise_var=noise), sd


def _path_weighted(adj: np.ndarray, score: np.ndarray) -> float | None:
    """Fraction of directed paths (any length, with multiplicity) along which ``score`` rises.

    Ties count one half, as in scm.py. Returns None when path counts overflow.
    """
    up = (score[None, :] > score[:, None]).astype(float)
    tie = (score[None, :] == score[:, None]).astype(float)
    m = adj.copy()
    inc = tot = 0.0
    for _ in range(adj.shape[0]):
        if not m.any():
            break
        tot += m.sum()
        inc += (m * up).sum() + 0.5 * (m * tie).sum()
        m = m @ adj
    if not np.isfinite(tot) or not np.isfinite(inc):
        return None
    return 0.5 if tot == 0 else float(inc / tot)


def sortability(sigma: np.ndarray, nodes: tuple[str, ...], dag) -> dict[str, float | None]:
    """Population var-sortability and R^2-sortability (Reisach et al. 2021, 2023)."""
    idx = {n: i for i, n in enumerate(nodes)}
    adj = np.zeros((len(nodes), len(nodes)))
    for a, b in dag.directed_edges:
        adj[idx[a], idx[b]] = 1.0
    var = np.diag(sigma).copy()
    # Exact ties matter here: a standardised model has all variances 1 up to rounding,
    # and must score 0.5, not a value decided by the last bit.
    var = np.round(var, 12)
    prec = np.linalg.inv(sigma)
    r2 = np.round(1.0 - 1.0 / (np.diag(prec) * np.diag(sigma)), 12)
    return {"var_sortability": _path_weighted(adj, var),
            "r2_sortability": _path_weighted(adj, r2)}
