"""Instrumented search: does cn / poss_de / pa move at all under `add`?
A violation needs cn to shrink (LOSS) or forb to shrink while pa gains (GAIN).
"""
import sys, json, random
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/code")
import numpy as np
import x2lib as X
from adjust import causal_nodes, poss_de, parents_of_set
from graphs import (meek_closure, v_structures, has_directed_cycle,
                    undirected_edges, dag_to_cpdag, random_dag)


def succ(G, ref):
    out = []
    for (u, v) in undirected_edges(G):
        for (a, b) in ((u, v), (v, u)):
            H = G.copy(); H[b, a] = 0
            H = meek_closure(H)
            if has_directed_cycle(H) or v_structures(H) != ref:
                continue
            out.append(((a, b), H))
    return out


def run(p, n_dags, seed, degs):
    rng = np.random.default_rng(seed); prng = random.Random(seed)
    s = dict(D_am=0, cn_shrink=0, cn_same=0, pa_grow=0, forb_shrink=0,
             pa_grow_and_forb_shrink=0, N_set=0, gain_cand=0)
    ex = []
    for _ in range(n_dags):
        D = random_dag(p, prng.choice(degs), rng)
        C = dag_to_cpdag(D); ref = v_structures(C)
        G = C
        for _step in range(6):
            S = succ(G, ref)
            if not S: break
            for x in range(p):
                for y in range(p):
                    if x == y: continue
                    O0, _, a0 = X.ostar_and_paths(G, x, y)
                    if not a0: continue
                    cn0 = causal_nodes(G, x, y)
                    pa0 = parents_of_set(G, cn0); fb0 = poss_de(G, cn0) | {x}
                    for (ab, H) in S:
                        O1, _, a1 = X.ostar_and_paths(H, x, y)
                        if not a1: continue
                        s["D_am"] += 1
                        cn1 = causal_nodes(H, x, y)
                        pa1 = parents_of_set(H, cn1); fb1 = poss_de(H, cn1) | {x}
                        if cn1 != cn0:
                            s["cn_shrink"] += 1
                        else:
                            s["cn_same"] += 1
                        pg = bool(pa1 - pa0); fs = bool(fb0 - fb1)
                        s["pa_grow"] += pg; s["forb_shrink"] += fs
                        if pg and fs: s["pa_grow_and_forb_shrink"] += 1
                        # gain candidate: a node newly a parent of cn AND newly out of forb
                        cand = (pa1 - pa0) & (fb0 - fb1)
                        if cand:
                            s["gain_cand"] += 1
                            if len(ex) < 3:
                                ex.append(dict(x=x, y=y, e=list(ab), cand=sorted(cand),
                                               G0=G.tolist(), H=H.tolist(),
                                               O0=sorted(O0), O1=sorted(O1)))
                        if O1 != O0:
                            s["N_set"] += 1
            G = S[prng.randrange(len(S))][1]
    return s, ex


if __name__ == "__main__":
    p = int(sys.argv[1]); n = int(sys.argv[2]); seed = int(sys.argv[3])
    degs = [float(d) for d in sys.argv[4].split(",")]
    s, ex = run(p, n, seed, degs)
    print(json.dumps(dict(p=p, seed=seed, stats=s, ex=ex), default=float)[:3000])
