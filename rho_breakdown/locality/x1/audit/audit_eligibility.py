"""X1 DESIGN AUDIT -- eligibility / pool-composition measurement ONLY.

Structural properties of the ensemble. Computes NOTHING about the S operator's
outcomes (no catch rate, no damage rate, no O* after perturbation).

Usage: python audit_eligibility.py <n_scm> <ensemble> <out.json> [nw]
"""
import json, sys
from multiprocessing import Pool
import numpy as np

from adjust import (cov_linear, is_amenable, optimal_adjustment_set,
                    possibly_causal_paths, total_effect_linear)
from graphs import (MeekFail, apply_background_knowledge, dag_to_cpdag,
                    directed_edges, is_directed, random_dag, skeleton,
                    undirected_edges)
from scm import make_linear_iscm
from run_arm1 import ENSEMBLES, OVERSAMPLE, hop_dist, MAX_K


def directed_reach(G, p):
    """Reachability using ONLY directed edges of G (a property of C alone)."""
    A = np.zeros((p, p), dtype=bool)
    for i, j in directed_edges(G):
        A[i, j] = True
    R = A.copy()
    for _ in range(p):
        R2 = R | (R @ A)
        if np.array_equal(R2, R):
            break
        R = R2
    return R


def one(args):
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
    if abs(tau) < 1e-6:
        return None
    K_all = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U]
    order = rng.permutation(len(K_all))
    K_all = [K_all[i] for i in order]
    if len(K_all) < 4:
        return None
    K = K_all[:4]
    try:
        G0 = apply_background_knowledge(C, K)
    except MeekFail:
        return None
    O0 = optimal_adjustment_set(G0, x, y)
    if O0 is None:
        return None
    if abs(tau) < 1e-3:
        return None                      # TAU_FLOOR applied at analysis in E1'

    S = skeleton(C)
    dist = hop_dist(C, [x, y], p)
    N = [(a, b) for a in range(p) for b in range(a + 1, p) if S[a, b] == 0]
    hopsN = [int(min(dist[a], dist[b])) for (a, b) in N]
    hopsK = [int(min(dist[u], dist[v])) for (u, v) in K]

    Rdir = directed_reach(C, p)          # directed-only reachability in C
    # ordered non-adjacent pairs (a,b): would a->b close a directed cycle in C?
    n_ord = 2 * len(N)
    n_ord_cyc = sum(1 for (a, b) in N for (u, v) in ((a, b), (b, a)) if Rdir[v, u])

    trueE = {(i, j) for i in range(p) for j in range(p) if D[i, j] == 1}
    loc_pool = [(a, b) for a in (x, y) for b in range(p)
                if b != a and (a, b) not in trueE]
    def label(a, b):
        if D[b, a] == 1:
            return "reversal_of_true_edge"
        if S[a, b] == 1:
            return "reorient_cpdag_edge"   # adjacent in C but not a true edge dir
        return "pure_spurious"
    loc_lab = [label(a, b) for (a, b) in loc_pool]
    loc_query = sum(1 for (a, b) in loc_pool if {a, b} == {x, y})

    # v-structures of C destroyed by adding {a,b}: a,b are endpoints of an
    # unshielded collider in C  -> structural, property of C alone
    from graphs import v_structures
    vs = v_structures(C)
    endpoints = {frozenset((a, c)) for (a, _b, c) in vs}
    n_N_kills_vstruct = sum(1 for (a, b) in N if frozenset((a, b)) in endpoints)

    return dict(seed=int(seed), p=int(p), deg=float(deg),
                n_U=len(U), n_N=len(N), n_A=int(S.sum() // 2),
                xy_nonadj=bool(S[x, y] == 0),
                hopsN=hopsN, hopsK=hopsK,
                n_ord=n_ord, n_ord_cyc=int(n_ord_cyc),
                n_N_kills_vstruct=int(n_N_kills_vstruct),
                n_vstruct_C=len(vs),
                loc_pool=len(loc_pool),
                loc_pure=sum(1 for l in loc_lab if l == "pure_spurious"),
                loc_rev=sum(1 for l in loc_lab if l == "reversal_of_true_edge"),
                loc_reo=sum(1 for l in loc_lab if l == "reorient_cpdag_edge"),
                loc_query=int(loc_query),
                O0=len(O0))


def main():
    n_scm = int(sys.argv[1]); ens = sys.argv[2]; out = sys.argv[3]
    nw = int(sys.argv[4]) if len(sys.argv) > 4 else 12
    cfg = ENSEMBLES[ens]; ps = list(cfg["ps"])
    rng = np.random.default_rng(20260819)
    over = OVERSAMPLE[ens]
    jobs = [(int(s), int(rng.choice(ps)), float(rng.choice(cfg["degs"])), ens)
            for s in range(n_scm * over)]
    res = []
    with Pool(nw) as pool:
        for r in pool.imap(one, jobs, chunksize=4):
            if r:
                res.append(r)
                if len(res) >= n_scm:
                    break
        pool.terminate()
    json.dump(res, open(out, "w"))
    print(f"[audit/{ens}] kept {len(res)} of <= {len(jobs)} draws -> {out}", flush=True)


if __name__ == "__main__":
    main()
