"""
SECONDARY ARM: nonlinear iSCM, finite-sample, nonparametric g-formula.

Everything the linear arm computes exactly is here estimated, so this arm carries
estimation noise. It exists to answer one question the linear arm cannot: does the
O*-relevance / ignorance-interval picture survive when the mechanisms are
nonlinear and the estimator is nonparametric?

The rho=0 row IS the estimator-noise baseline (correct O*, so any deviation from
tau is pure estimation error). Read the rho>=1 spreads against it.

Usage: python run_nonlinear.py [n_scm] [out.json]
"""
import json
import sys
from itertools import combinations
from multiprocessing import Pool

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor

from adjust import (is_valid_adjustment_set, optimal_adjustment_set,
                    possibly_causal_paths, is_amenable)
from graphs import (MeekFail, apply_background_knowledge, dag_agrees_with,
                    dag_to_cpdag, random_dag, undirected_edges)
from run_linear import generic_knowledge, tiered_knowledge, MAX_K
from scm import make_nonlinear_iscm, r2_sortability, sample_nonlinear, var_sortability

N_FIT = 4000
N_MC = 100_000
X1, X0 = 1.0, -1.0   # do(X=+1) vs do(X=-1), on the standardized scale


def gformula(X, x, y, Z, seed):
    """E[Y|do(x=X1)] - E[Y|do(x=X0)] by regression adjustment on (x, Z)."""
    S = [x] + sorted(Z)
    m = HistGradientBoostingRegressor(max_iter=200, random_state=seed)
    m.fit(X[:, S], X[:, y])
    A = X[:, S].copy()
    A[:, 0] = X1
    B = X[:, S].copy()
    B[:, 0] = X0
    return float(m.predict(A).mean() - m.predict(B).mean())


def analyse_one(args):
    seed, p, deg = args
    rng = np.random.default_rng(seed)
    D = random_dag(p, deg, rng)
    C = dag_to_cpdag(D)
    if len(undirected_edges(C)) < 3:
        return None
    pairs = [(x, y) for x in range(p) for y in range(p)
             if x != y and len(possibly_causal_paths(D, x, y)) > 0]
    if not pairs:
        return None
    x, y = pairs[int(rng.integers(len(pairs)))]

    scm = make_nonlinear_iscm(D, rng)
    # true effect: interventional Monte Carlo on the frozen SCM
    r1 = np.random.default_rng(seed + 1)
    r0 = np.random.default_rng(seed + 1)  # common random numbers
    tau = float(sample_nonlinear(scm, N_MC, r1, intervene={x: X1})[:, y].mean()
                - sample_nonlinear(scm, N_MC, r0, intervene={x: X0})[:, y].mean())
    if abs(tau) < 0.05:
        return None

    Xd = sample_nonlinear(scm, N_FIT, np.random.default_rng(seed + 2))
    diag = dict(varsort_iscm=var_sortability(Xd, D), r2sort_iscm=r2_sortability(Xd, D))

    out = dict(seed=seed, p=p, deg=deg, x=x, y=y, tau=tau, diag=diag,
               cpdag_amenable=bool(is_amenable(C, x, y)), arms={})
    cache = {}

    def est(Z):
        k = tuple(sorted(Z))
        if k not in cache:
            cache[k] = gformula(Xd, x, y, Z, seed)
        return cache[k]

    for arm in ("generic", "tiered"):
        K, _ = (generic_knowledge(D, C, rng), None) if arm == "generic" \
            else tiered_knowledge(D, C, rng)
        if len(K) < 3:
            continue
        try:
            G0 = apply_background_knowledge(C, K)
        except MeekFail:
            continue
        O0 = optimal_adjustment_set(G0, x, y)
        if O0 is None:
            continue
        members = []
        for rho in (1, 2):
            if rho > len(K):
                break
            for flip in combinations(range(len(K)), rho):
                Kp = [(v, u) if k in flip else (u, v) for k, (u, v) in enumerate(K)]
                rec = dict(rho=rho)
                try:
                    G = apply_background_knowledge(C, Kp)
                except MeekFail:
                    rec.update(consistent=False)
                    members.append(rec)
                    continue
                O = optimal_adjustment_set(G, x, y)
                rec.update(consistent=True, amenable=O is not None,
                           expels_true_dag=not dag_agrees_with(D, G))
                if O is not None:
                    rec["ostar_changed"] = (O != O0)
                    rec["ostar_valid"] = bool(is_valid_adjustment_set(D, x, y, O))
                    rec["est"] = est(O)
                    rec["bias"] = rec["est"] - tau
                members.append(rec)
        out["arms"][arm] = dict(n_K=len(K), O0=sorted(O0), est0=est(O0), members=members)
    return out if out["arms"] else None


def main():
    n_scm = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    out_path = sys.argv[2] if len(sys.argv) > 2 else "../results/nonlinear_raw.json"
    rng = np.random.default_rng(20260718)
    jobs = []
    for s in range(n_scm * 4):
        p = int(rng.integers(5, 9))
        deg = float(rng.choice([1.5, 2.0, 2.5]))
        jobs.append((s + 500_000, p, deg))
    with Pool(16) as pool:
        res = [r for r in pool.imap_unordered(analyse_one, jobs, chunksize=1) if r]
    res = res[:n_scm]
    json.dump(res, open(out_path, "w"))
    print(f"kept {len(res)} SCMs -> {out_path}")


if __name__ == "__main__":
    main()
