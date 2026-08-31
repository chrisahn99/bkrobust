"""
ARM 1 -- population Sigma, exact, zero Monte-Carlo.

Does the `se` criterion un-degenerate rho* WITHOUT the |K|=8 ensemble that
exists on 1.5% of the licensed box?

MAX_K = 4 (the pilot's own cap) and the ball is FULL: rho = 1..|K|, so
"censored" means "no combination of expert errors overturns", not "we stopped
at 3" (PREREG.md D5).

Emits per SCM everything needed to recompute rho* at any n, any z, any
criterion, plus the placebo ball and the locality features.

Usage: python run_arm1.py <n_scm> <ensemble> <out.json> [n_workers]
"""
import json
import sys
from itertools import combinations
from multiprocessing import Pool

import numpy as np

from adjust import (cov_linear, is_amenable, is_valid_adjustment_set,
                    optimal_adjustment_set, ols_coefficient,
                    possibly_causal_paths, total_effect_linear)
from graphs import (MeekFail, apply_background_knowledge, dag_agrees_with,
                    dag_to_cpdag, directed_edges, random_dag, skeleton,
                    undirected_edges)
from scm import make_linear_iscm
from se import deviation_radius, ident_radius, ols_with_se, rho_star_oracle

MAX_K = 4
TAU_FLOOR = 1e-3          # the pilot's filter; sensitivity at 1e-6 reported (AUDIT A8)

ENSEMBLES = {
    "original": dict(ps=range(5, 9),   degs=[1.5, 2.0, 2.5],           need_u=3, ks=(4,),       max_rho=4),
    "licensed": dict(ps=range(5, 11),  degs=[1.5, 2.0, 2.5, 4.0, 6.0], need_u=3, ks=(4,),       max_rho=4),
    "large":    dict(ps=range(15, 26), degs=[1.5, 2.0, 2.5],           need_u=3, ks=(4,),       max_rho=4),
    # only to evaluate the registered falsifier as literally written (|K|=8,
    # matched |K|=4 control, MAX_RHO=3 exactly as E1)
    "k8":       dict(ps=range(5, 11),  degs=[1.5, 2.0, 2.5, 4.0, 6.0], need_u=8, ks=(4, 6, 8),  max_rho=3),
}
OVERSAMPLE = {"original": 3, "licensed": 3, "large": 3, "k8": 120}


def hop_dist(C, srcs, p):
    """BFS on the skeleton of C from {x,y}."""
    S = skeleton(C)
    dist = np.full(p, 1 << 20, dtype=int)
    frontier = list(srcs)
    for v in frontier:
        dist[v] = 0
    while frontier:
        nxt = []
        for v in frontier:
            for w in np.flatnonzero(S[v]):
                if dist[w] > dist[v] + 1:
                    dist[w] = dist[v] + 1
                    nxt.append(int(w))
        frontier = nxt
    return dist


def build_arm(C, D, Sigma, x, y, K, max_rho, n_cpdag_dir):
    """Enumerate the ball around K; return the per-arm record or None."""
    k = len(K)
    try:
        G0 = apply_background_knowledge(C, K)
    except MeekFail:
        return None
    O0 = optimal_adjustment_set(G0, x, y)
    if O0 is None:
        return dict(n_K=k, mpdag_amenable=False)
    est0, se_unit0, sig2_0 = ols_with_se(Sigma, x, y, O0, np.inf)
    # se_factor: SE_n = se_factor / sqrt(n - k_reg - 1)
    S0 = [x] + sorted(O0)
    Sinv = np.linalg.inv(Sigma[np.ix_(S0, S0)])
    se_factor = float(np.sqrt(max(sig2_0, 0.0) * Sinv[0, 0]))

    members, members_frozen = [], []
    top = min(max_rho, k)
    for rho in range(1, top + 1):
        for flip in combinations(range(k), rho):
            Kp = [(v, u) if idx in flip else (u, v) for idx, (u, v) in enumerate(K)]
            rec = dict(rho=rho)
            try:
                G = apply_background_knowledge(C, Kp)
            except MeekFail:
                rec.update(consistent=False, amenable=False)
                members.append(rec)
                # placebo P1: the frozen-O0 ball keeps the SAME consistency
                # bookkeeping but never moves the adjustment set
                members_frozen.append(dict(rho=rho, consistent=False, amenable=False))
                continue
            O = optimal_adjustment_set(G, x, y)
            rec.update(consistent=True,
                       expels_true_dag=not dag_agrees_with(D, G),
                       cascade=len(directed_edges(G)) - n_cpdag_dir - len(Kp),
                       amenable=O is not None)
            if O is not None:
                e, _, s2 = ols_with_se(Sigma, x, y, O, np.inf)
                Sm = [x] + sorted(O)
                Sminv = np.linalg.inv(Sigma[np.ix_(Sm, Sm)])
                rec["est"] = e
                rec["se_factor"] = float(np.sqrt(max(s2, 0.0) * Sminv[0, 0]))
                rec["k_reg"] = len(Sm)
                rec["ostar_changed"] = (O != O0)
                rec["ostar_valid"] = bool(is_valid_adjustment_set(D, x, y, O))
            members.append(rec)
            # PLACEBO P1 -- adjustment set frozen at O0, so est == est0 exactly
            members_frozen.append(dict(rho=rho, consistent=True, amenable=True,
                                       est=est0))
    return dict(n_K=k, mpdag_amenable=True, K=[list(e) for e in K],
                O0=sorted(O0), est0=est0, se_factor=se_factor,
                k_reg=len(S0), sigma2=sig2_0, members=members,
                members_frozen=members_frozen)


def analyse_one(args):
    seed, p, deg, ens = args
    cfg = ENSEMBLES[ens]
    rng = np.random.default_rng(seed)
    D = random_dag(p, deg, rng)
    C = dag_to_cpdag(D)
    U = undirected_edges(C)
    if len(U) < cfg["need_u"]:
        return None
    pairs = [(a, b) for a in range(p) for b in range(p)
             if a != b and len(possibly_causal_paths(D, a, b)) > 0]
    if not pairs:
        return None
    x, y = pairs[int(rng.integers(len(pairs)))]
    scm = make_linear_iscm(D, rng)
    Sigma = cov_linear(scm["A"], scm["omega"])
    tau = total_effect_linear(scm["A"], x, y)
    if abs(tau) < 1e-6:                    # keep 1e-6 rows; TAU_FLOOR applied at analysis
        return None

    K_all = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U]
    order = rng.permutation(len(K_all))
    K_all = [K_all[i] for i in order]

    dist = hop_dist(C, [x, y], p)
    n_cpdag_dir = len(directed_edges(C))

    out = dict(seed=int(seed), p=int(p), deg=float(deg), x=int(x), y=int(y),
               tau=float(tau), n_undirected=len(U), n_edges=int(D.sum()),
               cpdag_amenable=bool(is_amenable(C, x, y)),
               max_rho=cfg["max_rho"], arms={})
    for k in cfg["ks"]:
        if len(K_all) < k:
            continue
        K = K_all[:k]
        arm = build_arm(C, D, Sigma, x, y, K, cfg["max_rho"], n_cpdag_dir)
        if arm is None:
            continue
        arm["stmt_hopdist"] = [int(min(dist[u], dist[v])) for (u, v) in K]
        out["arms"][f"K{k}"] = arm
    if not out["arms"]:
        return None
    return out


def main():
    n_scm = int(sys.argv[1]) if len(sys.argv) > 1 else 1200
    ens = sys.argv[2] if len(sys.argv) > 2 else "licensed"
    out_path = sys.argv[3] if len(sys.argv) > 3 else f"../results/arm1_{ens}.json"
    nw = int(sys.argv[4]) if len(sys.argv) > 4 else 12
    cfg = ENSEMBLES[ens]
    ps = list(cfg["ps"])
    rng = np.random.default_rng(20260819)
    over = OVERSAMPLE[ens]
    jobs = [(int(s), int(rng.choice(ps)), float(rng.choice(cfg["degs"])), ens)
            for s in range(n_scm * over)]
    with Pool(nw) as pool:
        res = []
        for r in pool.imap(analyse_one, jobs, chunksize=4):
            if r:
                res.append(r)
                if len(res) >= n_scm:
                    break
        pool.terminate()
    with open(out_path, "w") as f:
        json.dump(res, f)
    print(f"[arm1/{ens}] kept {len(res)} SCMs from <= {len(jobs)} draws -> {out_path}",
          flush=True)


if __name__ == "__main__":
    main()
