"""
PRIMARY ARM (exact, population-level, zero Monte-Carlo error).

For each random linear iSCM:
  - true DAG D*, CPDAG C, treatment X, outcome Y with a directed path
  - background knowledge K: orientations of undirected edges of C, TRUE by
    construction (hence Meek-consistent)
  - for rho = 1,2,3: enumerate ALL flip-sets of size rho (the expert got rho of
    their claims backwards), keep the Meek-CONSISTENT ones, and for each compute
    O* and the population effect estimate  ->  ignorance interval as f(rho)

The estimate under adjustment set Z is the POPULATION OLS coefficient of X in
Y ~ X + Z, computed in closed form from Sigma. No sampling noise anywhere: every
number below is exact for the stated SCM.

Usage: python run_linear.py [n_scm] [out.json]
"""
import json
import sys
from itertools import combinations
from multiprocessing import Pool

import numpy as np

from adjust import (cov_linear, is_valid_adjustment_set, optimal_adjustment_set,
                    ols_coefficient, possibly_causal_paths, total_effect_linear,
                    is_amenable)
from graphs import (MeekFail, apply_background_knowledge, dag_agrees_with,
                    dag_to_cpdag, directed_edges, random_dag, topological_order,
                    undirected_edges)
from scm import (make_linear_iscm, make_linear_naive, r2_sortability,
                 sample_linear, var_sortability)

MAX_K = 4        # |K| ; ball sizes C(4,rho) = 4, 6, 4
N_DIAG = 5000    # sample size for the sortability diagnostics


def tiered_knowledge(D, C, rng, max_stmt=MAX_K):
    """Bang & Didelez (arXiv:2306.01638) tiered knowledge: partition a topological
    order of D* into contiguous tiers; K = the cross-tier orientations of C's
    undirected edges. Such K needs only Meek's R1."""
    order = topological_order(D)
    p = len(order)
    n_tiers = int(rng.integers(2, min(4, p) + 1))
    cuts = sorted(rng.choice(np.arange(1, p), size=n_tiers - 1, replace=False))
    tier = np.zeros(p, dtype=int)
    t = 0
    for pos, node in enumerate(order):
        while t < len(cuts) and pos >= cuts[t]:
            t += 1
        tier[node] = t
    K = [((i, j) if tier[i] < tier[j] else (j, i))
         for (i, j) in undirected_edges(C) if tier[i] != tier[j]]
    # sanity: tiered K must agree with D*
    K = [(u, v) for (u, v) in K if D[u, v] == 1]
    if len(K) > max_stmt:
        idx = rng.choice(len(K), size=max_stmt, replace=False)
        K = [K[i] for i in idx]
    return K, tier.tolist()


def generic_knowledge(D, C, rng, max_stmt=MAX_K):
    U = undirected_edges(C)
    K_all = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U]
    if len(K_all) > max_stmt:
        idx = rng.choice(len(K_all), size=max_stmt, replace=False)
        K_all = [K_all[i] for i in idx]
    return K_all


def n_extra_orientations(G, C, K):
    """cascade size: edges oriented by the Meek closure beyond the |K| asserted."""
    return len(directed_edges(G)) - len(directed_edges(C)) - len(K)


def analyse_one(args):
    seed, p, deg = args
    rng = np.random.default_rng(seed)
    D = random_dag(p, deg, rng)
    C = dag_to_cpdag(D)
    U = undirected_edges(C)
    if len(U) < 3:
        return None

    pairs = [(x, y) for x in range(p) for y in range(p)
             if x != y and len(possibly_causal_paths(D, x, y)) > 0]
    if not pairs:
        return None
    x, y = pairs[int(rng.integers(len(pairs)))]

    scm = make_linear_iscm(D, rng)
    Sigma = cov_linear(scm["A"], scm["omega"])
    tau = total_effect_linear(scm["A"], x, y)
    if abs(tau) < 1e-3:
        return None

    # ---- sortability diagnostics (Reisach/Tami 2303.18211): iSCM vs naive control
    rng_d = np.random.default_rng(seed + 10 ** 7)
    Xi = sample_linear(scm, N_DIAG, rng_d)
    scm_naive = make_linear_naive(D, np.random.default_rng(seed))
    Xn = sample_linear(scm_naive, N_DIAG, rng_d)
    diag = dict(
        varsort_iscm=var_sortability(Xi, D), r2sort_iscm=r2sortability_safe(Xi, D),
        varsort_naive=var_sortability(Xn, D), r2sort_naive=r2sortability_safe(Xn, D),
    )

    out = dict(seed=seed, p=p, deg=deg, x=x, y=y, tau=tau,
               n_undirected=len(U), n_edges=int(D.sum()),
               cpdag_amenable=bool(is_amenable(C, x, y)), diag=diag, arms={})

    for arm in ("generic", "tiered"):
        if arm == "generic":
            K = generic_knowledge(D, C, rng)
            tier = None
        else:
            K, tier = tiered_knowledge(D, C, rng)
        if len(K) < 3:
            continue
        try:
            G0 = apply_background_knowledge(C, K)
        except MeekFail:
            continue  # cannot happen for true K, but be safe
        O0 = optimal_adjustment_set(G0, x, y)
        if O0 is None:
            continue  # MPDAG not amenable even under TRUE knowledge -> excluded
        est0 = ols_coefficient(Sigma, x, y, O0)

        members = []
        for rho in (1, 2, 3):
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
                rec.update(
                    consistent=True,
                    expels_true_dag=not dag_agrees_with(D, G),
                    cascade=n_extra_orientations(G, C, Kp),
                    amenable=O is not None,
                )
                if O is not None:
                    rec["ostar_changed"] = (O != O0)
                    rec["ostar_valid"] = bool(is_valid_adjustment_set(D, x, y, O))
                    rec["est"] = ols_coefficient(Sigma, x, y, O)
                    rec["bias"] = rec["est"] - tau
                members.append(rec)

        out["arms"][arm] = dict(K=[list(e) for e in K], n_K=len(K),
                                O0=sorted(O0), est0=est0, tier=tier,
                                mpdag_amenable=True, members=members)
    if not out["arms"]:
        return None
    return out


def r2sortability_safe(X, D):
    try:
        return r2_sortability(X, D)
    except np.linalg.LinAlgError:
        return float("nan")


def main():
    n_scm = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    out_path = sys.argv[2] if len(sys.argv) > 2 else "../results/linear_raw.json"
    rng = np.random.default_rng(20260717)
    jobs = []
    for s in range(n_scm * 3):  # oversample; many are rejected
        p = int(rng.integers(5, 9))
        deg = float(rng.choice([1.5, 2.0, 2.5]))
        jobs.append((s, p, deg))
    with Pool(16) as pool:
        res = [r for r in pool.imap_unordered(analyse_one, jobs, chunksize=4) if r]
    res = res[:n_scm]
    with open(out_path, "w") as f:
        json.dump(res, f)
    print(f"kept {len(res)} SCMs -> {out_path}")


if __name__ == "__main__":
    main()
