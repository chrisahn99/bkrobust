"""
The `se` criterion: rho* against SAMPLING ERROR instead of against tau.

The pilot's two criteria both compare the perturbed estimate to `tau`, the TRUE
total effect (analyse.py:57,59). No analyst has tau. This module replaces the
yardstick with the analyst's OWN standard error, which is computable from
(Sigma-hat, n, O0) alone.

    SE_n(est0) = sqrt( sigma0^2 * [Sigma_{S0 S0}^{-1}]_{XX} / n )
    d(rho)     = max_{m in ball(rho)} |est_m - est0|
    rho*_se    = min { rho : d(rho) > z * SE_n }        (|K|+1 == CENSORED)

Nothing here reads `tau`. `tau` is used ONLY by the evaluation code, never by
the instrument -- that separation is the point of the experiment and is
enforced by `assert_no_tau_leak` below.
"""
import numpy as np

Z95 = 1.959963984540054


# ------------------------------------------------------------------ the SE itself
def ols_with_se(Sigma, x, y, Z, n):
    """Population OLS coefficient of x in Y ~ X + Z, plus its analytic SE at
    sample size n.

    Everything is jointly Gaussian here, so Y | X,Z is linear with CONSTANT
    variance for ANY Z (valid or not) and the classical homoskedastic SE is the
    right one. Outside linear-Gaussian a sandwich SE would be required
    (AUDIT.md A4).

    Returns (beta_x, se, sigma2_resid).  se is np.inf for n <= |S| (no df).
    """
    S = [x] + sorted(Z)
    k = len(S)
    Sxx = Sigma[np.ix_(S, S)]
    Sxy = Sigma[np.ix_(S, [y])]
    try:
        Sinv = np.linalg.inv(Sxx)
    except np.linalg.LinAlgError:
        return np.nan, np.inf, np.nan
    beta = Sinv @ Sxy
    sigma2 = float(Sigma[y, y] - (Sxy.T @ beta).item())
    sigma2 = max(sigma2, 0.0)
    if not np.isfinite(n) or n <= k + 1:
        # n = inf  ->  SE = 0 (the population limit); n too small -> undefined
        se = 0.0 if np.isinf(n) else np.inf
    else:
        # df correction: sigma2_hat unbiased uses n-k-1
        se = float(np.sqrt(sigma2 * Sinv[0, 0] / (n - k - 1)))
    return float(beta[0, 0]), se, sigma2


def se_only(Sigma, x, y, Z, n):
    return ols_with_se(Sigma, x, y, Z, n)[1]


# ------------------------------------------------------------------ deviation radius
def deviation_radius(members, est0, max_rho):
    """d(rho) for rho = 1..max_rho, over consistent+amenable members only.

    Non-amenable members carry NO estimate -- identification failed and the
    analyst sees it. They are excluded here and counted separately as the
    `ident` event (AUDIT.md A7); letting them vanish would make the ball look
    quieter than it is.

    Monotone non-decreasing by construction (PREREG.md 5.2) -- gate G4.
    """
    d = np.zeros(max_rho + 1)
    for rho in range(1, max_rho + 1):
        best = d[rho - 1]
        for m in members:
            if m["rho"] != rho or not m.get("consistent") or not m.get("amenable"):
                continue
            best = max(best, abs(m["est"] - est0))
        d[rho] = best
    return d


def ident_radius(members, max_rho):
    """smallest rho at which SOME consistent member is non-amenable."""
    for rho in range(1, max_rho + 1):
        for m in members:
            if m["rho"] == rho and m.get("consistent") and not m.get("amenable"):
                return rho
    return max_rho + 1


# ------------------------------------------------------------------ the criteria
def rho_star_se(d, se, z=Z95, max_rho=None):
    """min rho with d(rho) > z*se ; max_rho+1 == CENSORED."""
    if max_rho is None:
        max_rho = len(d) - 1
    thr = z * se
    for rho in range(1, max_rho + 1):
        if d[rho] > thr:
            return rho
    return max_rho + 1


def rho_star_oracle(members, tau, est0, max_rho, crit):
    """The pilot's ORACLE criteria, kept verbatim for continuity with E1.
    `sign` and `rel50` both read tau -- that is the defect under test."""
    for rho in range(1, max_rho + 1):
        for m in members:
            if m["rho"] != rho or not m.get("consistent"):
                continue
            if crit == "ident":
                if not m.get("amenable"):
                    return rho
                continue
            if not m.get("amenable"):
                continue
            e = m["est"]
            if crit == "sign" and np.sign(e) != np.sign(tau):
                return rho
            if crit == "rel50" and abs(e - tau) > 0.5 * abs(tau):
                return rho
    return max_rho + 1


def ratio_curve(d, se, z=Z95):
    """R(rho, n) = d(rho) / (z*se) -- CONTINUOUS, never censored.

    The reviewer's objection was that a 4-valued top-censored statistic is not
    an instrument. R has no bins and no ceiling. inf when se == 0 and d > 0.
    """
    thr = z * se
    if thr == 0:
        return np.where(d > 0, np.inf, 0.0)
    return d / thr


# ------------------------------------------------------------------ integrity
def assert_no_tau_leak(fn_source):
    """Crude but load-bearing: the analyst-facing path must not mention tau."""
    return "tau" not in fn_source


def check_monotone_n(rhos_by_n, n_grid):
    """Gate G3: rho*_se is non-increasing in n, per SCM, by construction."""
    order = np.argsort(n_grid)
    v = np.asarray(rhos_by_n, dtype=float)[order]
    return bool(np.all(np.diff(v) <= 0))


def check_monotone_rho(d):
    """Gate G4."""
    return bool(np.all(np.diff(np.asarray(d, dtype=float)) >= -1e-12))
