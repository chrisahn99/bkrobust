"""Self-contained PC for linear-Gaussian data. No external causal-discovery dependency.

Two CI oracles, deliberately sharing ONE code path so the finite-sample arm cannot silently
differ from the population arm in anything but the test:
  * `oracle`     -- partial correlation from the TRUE covariance, |r| < tol  =>  independent.
  * `fisher_z`   -- partial correlation from the SAMPLE covariance, Fisher-z at level alpha.

GATE (test_population_recovers_cpdag): with the oracle test, PC must return EXACTLY
dag_to_cpdag(D). If that ever fails, this file is wrong and no number produced with it means
anything. It is run before every experiment, not once.
"""
import sys, os
from itertools import combinations
import numpy as np

sys.path.insert(0, os.environ.get("X2CODE",
    "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/code"))
from graphs import meek_closure, dag_to_cpdag, skeleton  # noqa: E402

_SQRT2 = np.sqrt(2.0)


def _norm_ppf(q):
    """Inverse standard normal CDF via erfinv (avoids a scipy dependency)."""
    from math import erf  # noqa: F401
    # bisection on erf is overkill; use scipy if present, else a rational approximation
    try:
        from scipy.stats import norm
        return float(norm.ppf(q))
    except Exception:
        # Acklam's approximation, accurate to ~1e-9
        a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
             1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
        b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
             6.680131188771972e+01, -1.328068155288572e+01]
        c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
             -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
        d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
             3.754408661907416e+00]
        pl, ph = 0.02425, 1 - 0.02425
        if q < pl:
            x = np.sqrt(-2 * np.log(q))
            return (((((c[0]*x+c[1])*x+c[2])*x+c[3])*x+c[4])*x+c[5]) / ((((d[0]*x+d[1])*x+d[2])*x+d[3])*x+1)
        if q > ph:
            x = np.sqrt(-2 * np.log(1 - q))
            return -(((((c[0]*x+c[1])*x+c[2])*x+c[3])*x+c[4])*x+c[5]) / ((((d[0]*x+d[1])*x+d[2])*x+d[3])*x+1)
        x = q - 0.5
        r = x * x
        return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*x / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)


def partial_corr(S, i, j, cond):
    """Partial correlation of i,j given `cond`, from covariance S."""
    idx = [i, j] + list(cond)
    M = S[np.ix_(idx, idx)]
    try:
        P = np.linalg.inv(M + 1e-12 * np.eye(len(idx)))
    except np.linalg.LinAlgError:
        return 0.0
    den = np.sqrt(P[0, 0] * P[1, 1])
    if den <= 0:
        return 0.0
    r = -P[0, 1] / den
    return float(np.clip(r, -0.999999, 0.999999))


def ci_oracle(S, i, j, cond, n=None, tol=1e-8, zcrit=None):
    return abs(partial_corr(S, i, j, cond)) < tol


def ci_fisher_z(S, i, j, cond, n=None, tol=None, zcrit=None):
    r = partial_corr(S, i, j, cond)
    dof = n - len(cond) - 3
    if dof <= 0:
        return True
    z = 0.5 * np.log((1 + r) / (1 - r))
    return abs(np.sqrt(dof) * z) < zcrit


def pc(S, p, test, n=None, alpha=0.01, tol=1e-8, max_cond=None):
    """PC. Returns (cpdag_amat, n_tests). Adjacency convention matches graphs.py."""
    zcrit = _norm_ppf(1 - alpha / 2.0)
    adj = np.ones((p, p), dtype=bool)
    np.fill_diagonal(adj, False)
    sep = {}
    ntest = 0
    d = 0
    cap = p - 2 if max_cond is None else min(max_cond, p - 2)
    while d <= cap:
        edges = [(i, j) for i in range(p) for j in range(i + 1, p) if adj[i, j]]
        if not edges:
            break
        progressed = False
        for (i, j) in edges:
            if not adj[i, j]:
                continue
            for a, b in ((i, j), (j, i)):
                nb = [k for k in range(p) if adj[a, k] and k != b]
                if len(nb) < d:
                    continue
                progressed = True
                for cond in combinations(sorted(nb), d):
                    ntest += 1
                    if test(S, a, b, cond, n=n, tol=tol, zcrit=zcrit):
                        adj[i, j] = adj[j, i] = False
                        sep[(i, j)] = sep[(j, i)] = set(cond)
                        break
                if not adj[i, j]:
                    break
        if not progressed:
            break
        d += 1

    G = np.zeros((p, p), dtype=np.int8)
    G[adj] = 1
    for k in range(p):
        for i, j in combinations([m for m in range(p) if adj[m, k]], 2):
            if adj[i, j]:
                continue
            if k not in sep.get((i, j), set()):
                G[k, i] = 0
                G[k, j] = 0
    return meek_closure(G), ntest


def test_population_recovers_cpdag(n_graphs=200, ps=(4, 5, 6, 7), degs=(1.5, 2.0, 2.5, 3.0), seed=20260823):
    """THE GATE. Oracle-CI PC must reproduce dag_to_cpdag(D) exactly."""
    from graphs import random_dag
    from scm import make_linear_iscm
    rng = np.random.default_rng(seed)
    bad = 0
    for _ in range(n_graphs):
        p = int(rng.choice(ps)); deg = float(rng.choice(degs))
        D = random_dag(p, deg, rng)
        scm = make_linear_iscm(D, rng)
        A, om = scm["A"], scm["omega"]
        I = np.eye(p)
        M = np.linalg.inv(I - A.T)
        S = M @ np.diag(om) @ M.T
        Chat, _ = pc(S, p, ci_oracle)
        if not np.array_equal(Chat, dag_to_cpdag(D)):
            bad += 1
    return bad, n_graphs


if __name__ == "__main__":
    bad, tot = test_population_recovers_cpdag(int(sys.argv[1]) if len(sys.argv) > 1 else 200)
    print(f"GATE population-PC == dag_to_cpdag : {tot - bad}/{tot} exact"
          f"  -> {'PASS' if bad == 0 else 'FAIL (' + str(bad) + ' mismatches)'}")
    sys.exit(0 if bad == 0 else 1)
