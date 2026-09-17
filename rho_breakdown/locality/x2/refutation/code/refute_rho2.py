"""X2 REFUTER, run 2: the rho>=2 "sharp radius" claim.

X2's §4.4 claims: knowledge sets whose CLOSEST misstatement is at dmin>=2 gave
0 O*-changes in 6,726 amenable perturbations, therefore "the radius is sharp".
My re-analysis of the archive shows that stratum is 83% `large`, whose own
dmin=1 rate is 2e-4, and whose within-SCM matched near rate is 1/5,253:
the stratified expectation is ~6 events, and the MATCHED expectation ~0.5.

This run supplies the matched design X2 lacks: on the SAME graph, with the SAME
base K, flip subsets of NEAR statements and subsets of FAR statements and
compare.  Families are chosen so the near rate is LARGE (dense pockets) while
dmin>=2 edges exist (chain of pockets).
"""
import sys, json, time, itertools, collections
import numpy as np
sys.path.insert(0, "/home/costaj/latent-causal/x2-locality/code")
sys.path.insert(0, "/home/costaj/latent-causal/x2-locality-refute/code")
from graphs import (dag_to_cpdag, undirected_edges, apply_background_knowledge,
                    MeekFail, skeleton, random_dag)
from adjust import optimal_adjustment_set, is_valid_adjustment_set, possibly_causal_paths
from refute_search import hop, draw_family, dag_from_skel, INF
from multiprocessing import Pool

MAXK = 10          # cap |K|
MAXRHO = 3
MAXFAR = 7         # cap statements enumerated per stratum (combinatorics)

def one(args):
    seed, fam, qrule = args
    rng = np.random.default_rng(seed)
    S = draw_family(fam, rng); p = S.shape[0]
    D = dag_from_skel(S, rng)
    C = dag_to_cpdag(D); U = undirected_edges(C)
    if len(U) < 3: return None
    pairs = [(a, b) for a in range(p) for b in range(p)
             if a != b and len(possibly_causal_paths(D, a, b)) > 0]
    if not pairs: return None
    if qrule == "argmax_delta":
        best, bd = [], -1
        for (a, b) in pairs:
            d = hop(C, [a], p)[b]
            if d < INF and d > bd: best, bd = [(a, b)], d
            elif d == bd: best.append((a, b))
        x, y = best[int(rng.integers(len(best)))]
    else:
        x, y = pairs[int(rng.integers(len(pairs)))]
    dist = hop(C, [x, y], p)
    K_all = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U]
    dmins = [int(min(dist[u], dist[v])) for (u, v) in U]
    idx = list(rng.permutation(len(K_all)))[:MAXK]
    K = [K_all[i] for i in idx]; dm = [dmins[i] for i in idx]
    try: G0 = apply_background_knowledge(C, K)
    except MeekFail: return None
    O0 = optimal_adjustment_set(G0, x, y)
    if O0 is None: return None                       # base must be identified
    near = [i for i, d in enumerate(dm) if d <= 1]
    far  = [i for i, d in enumerate(dm) if 2 <= d < INF]
    out = collections.Counter(); hits = []
    for label, pool_ in (("near", near[:MAXFAR]), ("far", far[:MAXFAR])):
        if not pool_: continue
        for rho in range(1, min(MAXRHO, len(pool_)) + 1):
            for flip in itertools.combinations(pool_, rho):
                Kp = [(v, u) if i in flip else (u, v) for i, (u, v) in enumerate(K)]
                key = f"{label}/rho{rho}"
                out[key + "/D_all"] += 1
                try: G = apply_background_knowledge(C, Kp)
                except MeekFail: continue
                out[key + "/D_con"] += 1
                O = optimal_adjustment_set(G, x, y)
                if O is None: continue
                out[key + "/D_am"] += 1
                if O != O0:
                    out[key + "/N_set"] += 1
                    if not is_valid_adjustment_set(D, x, y, O):
                        out[key + "/N_bias"] += 1
                    if label == "far" and len(hits) < 30:
                        hits.append(dict(seed=int(seed), fam=fam, p=int(p), x=int(x), y=int(y),
                                         rho=rho, flip=[list(K[i]) for i in flip],
                                         dmins=[dm[i] for i in flip],
                                         dmaxs=[int(max(dist[K[i][0]], dist[K[i][1]])) for i in flip],
                                         K=[list(e) for e in K], K_dmin=dm,
                                         O0=sorted(O0), O1=sorted(O),
                                         valid=bool(is_valid_adjustment_set(D, x, y, O)),
                                         D_edges=[[int(a), int(b)] for a in range(p) for b in range(p)
                                                  if D[a, b] == 1 and D[b, a] == 0]))
    return dict(counts=dict(out), hits=hits,
                has_both=bool(near and far))

def main():
    n = int(sys.argv[1]); fam = sys.argv[2]; qrule = sys.argv[3]
    nw = int(sys.argv[4]); s0 = int(sys.argv[5]); tag = sys.argv[6]
    jobs = [(s0 + i, fam, qrule) for i in range(n)]
    tot = collections.Counter(); hits = []; nk = 0; nboth = 0; t0 = time.time()
    with Pool(nw) as pool:
        for r in pool.imap_unordered(one, jobs, chunksize=8):
            if not r: continue
            nk += 1; nboth += r["has_both"]; tot.update(r["counts"]); hits.extend(r["hits"])
    res = dict(fam=fam, qrule=qrule, n=n, n_keep=nk, n_both=nboth,
               secs=time.time()-t0, counts=dict(tot), n_hits=len(hits), hits=hits[:100])
    json.dump(res, open(f"/home/costaj/latent-causal/x2-locality-refute/results/{tag}.json", "w"))
    print(f"[{tag}] fam={fam} q={qrule} kept {nk}/{n} (both strata {nboth}) {res['secs']:.0f}s")
    for lab in ("near", "far"):
        for rho in (1, 2, 3):
            k = f"{lab}/rho{rho}"
            da = tot.get(k+"/D_am", 0)
            if da: print(f"   {lab:4s} rho={rho}  N_set={tot.get(k+'/N_set',0):7d} "
                         f"N_bias={tot.get(k+'/N_bias',0):7d} / D_am={da:8d} "
                         f"D_con={tot.get(k+'/D_con',0):8d} D_all={tot.get(k+'/D_all',0):8d} "
                         f"= {tot.get(k+'/N_set',0)/da:.5f}")
    print(f"   FAR HITS: {len(hits)}")

if __name__ == "__main__":
    main()
