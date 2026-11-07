"""The coefficient-agreement test of Su and Henckel (2022), rank-estimation variant.

Under the null that the candidate graph is right, every set valid in it gives an OLS coefficient
of x that converges to the same total effect. The test stacks those coefficients, takes the
contrasts with the committed set, and builds a Wald statistic from their joint covariance. That
covariance is generally degenerate (the paper's central difficulty), so its rank is estimated:
eigen-directions whose variance is below ``tol`` times the variance of the committed set's
influence function are treated as exactly null. The covariance comes from influence functions
(a sandwich), so it holds for non-Gaussian and heteroscedastic errors.

Written from the paper's description; no reference code was used (see PREREG.md, deviations).
"""

from __future__ import annotations

import numpy as np
from scipy.stats import chi2, norm


def coefficient_and_influence(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[float, np.ndarray]:
    """OLS coefficient of x in y ~ 1 + x + z, and its influence function (length n).

    Frisch-Waugh: with rx, ry the residuals of x and y on (1, z), beta = <rx, ry> / <rx, rx> and
    the influence of unit i is rx_i * e_i / mean(rx^2), e the full-regression residual.
    """
    n = len(x)
    a = np.column_stack([np.ones(n), z]) if z.size else np.ones((n, 1))
    q, _ = np.linalg.qr(a)
    rx = x - q @ (q.T @ x)
    ry = y - q @ (q.T @ y)
    sxx = float(rx @ rx)
    beta = float(rx @ ry) / sxx
    e = ry - beta * rx
    return beta, rx * e / (sxx / n)


def agreement_test(betas: np.ndarray, psi: np.ndarray, tol: float, alpha: float) -> dict:
    """Wald test that all coefficients equal the first (the committed set's).

    Args:
        betas: K coefficients, the committed set first.
        psi: n x K influence functions in the same order.
        tol: relative eigenvalue threshold of the rank estimate.
        alpha: level.
    """
    n, k = psi.shape
    ref = float(np.mean(psi[:, 0] ** 2))
    d = betas[1:] - betas[0]
    if k < 2 or ref <= 0.0:
        return {"stat": 0.0, "df": 0, "p": 1.0, "reject": False}
    pd = psi[:, 1:] - psi[:, [0]]
    pd = pd - pd.mean(axis=0)
    v = pd.T @ pd / n  # asymptotic covariance of sqrt(n) * d
    lam, u = np.linalg.eigh(v)
    keep = lam > tol * ref
    df = int(keep.sum())
    if df == 0:
        return {"stat": 0.0, "df": 0, "p": 1.0, "reject": False}
    proj = u[:, keep].T @ d
    stat = float(n * np.sum(proj ** 2 / lam[keep]))
    p = float(chi2.sf(stat, df))
    return {"stat": stat, "df": df, "p": p, "reject": bool(p < alpha)}


def pairwise_bonferroni(betas: np.ndarray, psi: np.ndarray, tol: float, alpha: float) -> dict:
    """Secondary: one z-test per contrast with the committed set, Bonferroni over usable pairs."""
    n, k = psi.shape
    ref = float(np.mean(psi[:, 0] ** 2))
    ps = []
    for j in range(1, k):
        dj = psi[:, j] - psi[:, 0]
        var = float(np.var(dj))
        if ref > 0.0 and var > tol * ref:
            ps.append(float(2.0 * norm.sf(abs(betas[j] - betas[0]) / np.sqrt(var / n))))
    if not ps:
        return {"m": 0, "min_p": 1.0, "reject": False}
    return {"m": len(ps), "min_p": min(ps), "reject": bool(min(ps) < alpha / len(ps))}
