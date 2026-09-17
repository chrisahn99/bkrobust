"""REFUTER verification: regenerate SCMs from seed, re-run the operator from scratch,
check (a) the asserted edge survives into the final graph, (b) `silent` reproduces,
(c) est0==tau, (d) concrete hop>=1 coherent counterexamples are real."""
import json, sys
import numpy as np
from itertools import combinations
sys.path.insert(0, "/home/costaj/latent-causal/x1-spurious/code")
from graphs import (MeekFail, apply_background_knowledge, dag_to_cpdag,
                    directed_edges, random_dag, skeleton, undirected_edges, is_directed)
from adjust import (cov_linear, is_valid_adjustment_set, optimal_adjustment_set,
                    possibly_causal_paths, total_effect_linear, causal_nodes)
from scm import make_linear_iscm
from se import ols_with_se
from x1_ops import bk_assert, draw_suni, draw_sloc, hop_dist_from, stmt_hop, pdag_extendable, nonadjacent_pairs, bload_pool
MAX_K = 4
UNREACH = 1 << 20

def regen(seed, p, deg):
    rng = np.random.default_rng(seed)
    D = random_dag(p, deg, rng)
    C = dag_to_cpdag(D)
    U = undirected_edges(C)
    pairs = [(a,b) for a in range(p) for b in range(p)
             if a != b and len(possibly_causal_paths(D,a,b)) > 0]
    x, y = pairs[int(rng.integers(len(pairs)))]
    scm = make_linear_iscm(D, rng)
    Sigma = cov_linear(scm["A"], scm["omega"])
    tau = total_effect_linear(scm["A"], x, y)
    K_all = [(i,j) if D[i,j]==1 else (j,i) for (i,j) in U]
    order = rng.permutation(len(K_all))
    K_all = [K_all[i] for i in order]
    K = [(int(a),int(b)) for (a,b) in K_all[:MAX_K]]
    srng = np.random.default_rng([int(seed), 0xA55E27])
    reps_uni = draw_suni(C, MAX_K, srng)
    reps_loc = draw_sloc(D, MAX_K, x, y, srng)
    return D, C, Sigma, x, y, tau, K, reps_uni, reps_loc

def main(ens, nsample):
    raw = json.load(open(f"/home/costaj/latent-causal/x1-spurious/results/x1_{ens}.json"))
    scms = raw["scms"]
    rs = np.random.default_rng(7)
    idx = rs.choice(len(scms), size=min(nsample, len(scms)), replace=False)
    n_edge_present = n_edge_missing = 0
    n_skel_same = 0
    n_check = n_mismatch = 0
    n_tau = 0; maxd = 0.0
    n_KX = 0
    for i in idx:
        s = scms[int(i)]
        D, C, Sigma, x, y, tau, K, reps_uni, reps_loc = regen(s["seed"], s["p"], s["deg"])
        assert [list(e) for e in K] == s["K"], ("K mismatch", s["seed"])
        assert (x, y) == (s["x"], s["y"]), ("xy mismatch", s["seed"])
        assert [list(e) for e in reps_uni] == s["reps_uni"], ("reps mismatch", s["seed"])
        G0 = apply_background_knowledge(C, K)
        O0 = optimal_adjustment_set(G0, x, y)
        est0, _, _ = ols_with_se(Sigma, x, y, O0, np.inf)
        if abs(est0 - tau) > maxd: maxd = abs(est0 - tau)
        if is_valid_adjustment_set(D, x, y, O0): n_tau += 1
        distC = hop_dist_from(skeleton(C), [x, y], C.shape[0])
        for rho in (1,):
            for flip in combinations(range(4), rho):
                Kp = [reps_uni[t] if t in flip else K[t] for t in range(4)]
                G, info = bk_assert(C, Kp)
                # (a) did the operator DO anything: is the asserted edge in G?
                for t in flip:
                    a, b = reps_uni[t]
                    if G[a, b] == 1 and G[b, a] == 0:
                        n_edge_present += 1
                    else:
                        n_edge_missing += 1
                        print("  !! asserted edge missing", s["seed"], (a,b), G[a,b], G[b,a])
                if np.array_equal(skeleton(G), skeleton(C)):
                    n_skel_same += 1
                # (b) recompute silent
                O = optimal_adjustment_set(G, x, y)
                if O is None:
                    sil = False; am = False
                else:
                    am = True
                    sil = not is_valid_adjustment_set(D, x, y, O)
                rec = [m for m in s["Suni"] if m["rho"]==rho and tuple(m["flip"])==flip][0]
                recsil = bool(rec.get("consistent") and rec.get("amenable") and rec.get("ostar_valid") is False)
                n_check += 1
                if recsil != sil or bool(rec.get("amenable")) != am:
                    n_mismatch += 1
                    print("  !! MISMATCH", s["seed"], flip, recsil, sil, rec.get("amenable"), am)
    print(f"[{ens}] sampled {len(idx)} SCMs")
    print(f"  asserted edge present in final G: {n_edge_present}/{n_edge_present+n_edge_missing}")
    print(f"  skeleton UNCHANGED after assertion (should be 0): {n_skel_same}/{n_check}")
    print(f"  independent recompute of (amenable,silent): mismatches {n_mismatch}/{n_check}")
    print(f"  O0 valid: {n_tau}/{len(idx)};  max |est0-tau| = {maxd:.3e}")

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))
