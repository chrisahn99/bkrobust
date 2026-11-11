"""Coefficient generators for E20, all returning a bkrobust LinearSEM over the true DAG.

naive  -- the submitted ``bkrobust.demo.evaluate.random_sem`` (unchanged, imported).
iscm   -- that same draw, internally standardised at the population level: each node's
          incoming weights and noise variance are divided by the standard deviation the node
          would have had given its already-standardised parents (Ormaniec et al., arXiv:2406.11601,
          in the infinite reference-sample limit).
uumc   -- Herman, Wahl, Ninad, Runge (CLeaR 2025, arXiv:2503.17037v3), Algorithm 1, in
          population form: for a node with d parents draw a'' ~ N(0, I_d), r = u^(1/d) with
          u ~ U(0, 1), a' = r a''/||a''||, s' = sqrt(1 - r^2), sigma' = sqrt(a'^T R a' + s'^2),
          coefficients a'/sigma', noise sd s'/sigma'; roots get unit noise. Ported from the paper's
          pseudocode, not from the UUMCdata package (not installed).

Sortability (Herman et al. App. B, pair-counted; Reisach et al. 2021/2023, path-length-counted)
is computed on the exact covariance, so it carries no sampling error.
"""
from __future__ import annotations

import numpy as np

from bkrobust.demo.evaluate import LinearSEM, random_sem


def topological_order(dag) -> list[str]:
    """Kahn's algorithm with a sorted frontier, so the draw order never depends on set order."""
    parents = {n: set(dag.parents(n)) for n in dag.nodes}
    done: list[str] = []
    placed: set[str] = set()
    while len(done) < len(parents):
        ready = sorted(n for n in parents if n not in placed and parents[n] <= placed)
        if not ready:
            raise ValueError("graph has a directed cycle")
        done.extend(ready)
        placed.update(ready)
    return done


def _standardised_build(dag, coef_of):
    """Build a unit-variance LinearSEM node by node.

    ``coef_of(node, parents, sigma_pp)`` returns (raw weights over ``parents``, raw noise variance)
    given the covariance of the already-standardised parents; the node is then rescaled to unit
    variance. The covariance is tracked incrementally because each node only needs its parents'.
    """
    nodes = list(dag.nodes)
    idx = {n: i for i, n in enumerate(nodes)}
    sigma = np.zeros((len(nodes), len(nodes)))
    placed: list[int] = []
    weights: dict[tuple[str, str], float] = {}
    noise_var: dict[str, float] = {}
    for v in topological_order(dag):
        pa = sorted(dag.parents(v))
        pi = [idx[p] for p in pa]
        s_pp = sigma[np.ix_(pi, pi)]
        w, nv = coef_of(v, pa, s_pp)
        w = np.asarray(w, dtype=float)
        var = float(w @ s_pp @ w) + nv if pa else nv
        sd = np.sqrt(var)
        w = w / sd
        nv = nv / var
        i = idx[v]
        if pa:
            cov_row = w @ sigma[np.ix_(pi, placed)] if placed else np.zeros(0)
            sigma[i, placed] = cov_row
            sigma[placed, i] = cov_row
        sigma[i, i] = 1.0
        placed.append(i)
        for p, wv in zip(pa, w, strict=True):
            weights[(p, v)] = float(wv)
        noise_var[v] = float(nv)
    return LinearSEM(dag=dag, weights=weights, noise_var=noise_var)


def naive_sem(dag, rng: np.random.Generator) -> LinearSEM:
    return random_sem(dag, rng)


def iscm_from(sem: LinearSEM) -> LinearSEM:
    """Population iSCM of a given raw draw: same signs and ratios, unit marginal variances."""
    def coef_of(v, pa, _s_pp):
        return [sem.weights[(p, v)] for p in pa], sem.noise_var[v]
    return _standardised_build(sem.dag, coef_of)


def uumc_sem(dag, rng: np.random.Generator) -> LinearSEM:
    """UUMC, Algorithm 1. Draws consume ``rng`` in topological order, parents sorted."""
    def coef_of(v, pa, _s_pp):
        d = len(pa)
        if d == 0:
            return [], 1.0
        a2 = rng.standard_normal(d)
        r = rng.uniform() ** (1.0 / d)
        a1 = r * a2 / np.linalg.norm(a2)
        s1 = np.sqrt(max(1.0 - r * r, 0.0))
        # _standardised_build divides by sqrt(a1' R a1 + s1^2) = sigma', as in the paper.
        return a1, s1 * s1
    return _standardised_build(dag, coef_of)


# --- sortability -----------------------------------------------------------


def path_structure(dag) -> tuple[np.ndarray, np.ndarray]:
    """(ancestor pairs, number of distinct path lengths per pair), node order = dag.nodes.

    Float matmuls with re-binarisation, because numpy's boolean matmul does not use BLAS.
    """
    nodes = list(dag.nodes)
    idx = {n: i for i, n in enumerate(nodes)}
    p = len(nodes)
    e = np.zeros((p, p))
    for a, b in dag.directed_edges:
        e[idx[a], idx[b]] = 1.0
    ek = e.copy()
    n_lengths = np.zeros((p, p))
    for _ in range(max(p - 1, 0)):
        if not ek.any():
            break
        n_lengths += ek
        ek = ((ek @ e) > 0).astype(float)
    return n_lengths > 0, n_lengths


def _ordered_share(metric: np.ndarray, weight: np.ndarray, tol: float = 1e-9) -> float:
    """Weighted share of (ancestor j, descendant i) pairs with metric[i] > metric[j]; ties 1/2."""
    tot = weight.sum()
    if tot == 0:
        return 0.5
    # Relative tolerance as in Herman et al.'s M/M.T > 1 + tol, written as a difference so that
    # a metric of 0 (an R2 of an isolated node) cannot divide by zero.
    diff = metric[None, :] - metric[:, None]
    thr = tol * np.maximum(np.abs(metric)[:, None], 1e-12)
    up = diff > thr
    tie = np.abs(diff) <= thr
    return float(((up + 0.5 * tie) * weight).sum() / tot)


def sortabilities(sigma: np.ndarray, anc: np.ndarray, n_lengths: np.ndarray) -> dict[str, float]:
    var = np.diag(sigma).copy()
    prec = np.linalg.inv(sigma)
    r2 = 1.0 - 1.0 / (var * np.diag(prec))
    w_pair = anc.astype(float)
    return {
        "var_sortability_pair": _ordered_share(var, w_pair),
        "r2_sortability_pair": _ordered_share(r2, w_pair),
        "var_sortability_reisach": _ordered_share(var, n_lengths),
        "r2_sortability_reisach": _ordered_share(r2, n_lengths),
    }
