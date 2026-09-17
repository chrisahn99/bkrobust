"""Large multiprocess refutation search with full instrumentation.

For every sampled MPDAG G0 (reachable state of a CPDAG), every ordered query
(x,y) with G0 amenable, every undirected edge, both orientations, Meek-closed:
  - guard rejections (cycle / new v-structure)      -> is the guard vacuous?
  - amenability of H                                -> is amenability preserved?
  - cn_H vs cn_0, forb_H vs forb_0, pa vs pa
  - PARTIAL LOSS near-miss: u in pa_0(cn_0)\forb_0 whose set of children in cn
    strictly shrinks (a TOTAL loss would be a counterexample)
  - N_set: the violation itself
"""
import sys, json, random
from multiprocessing import Pool
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/code")
import numpy as np
import x2lib as X
from adjust import causal_nodes, poss_de, parents_of_set
from graphs import (meek_closure, v_structures, has_directed_cycle, is_directed,
                    undirected_edges, dag_to_cpdag, random_dag)

VIOL = []


def job(args):
    p, degs, seed, n_dags, walk_len = args
    rng = np.random.default_rng(seed); prng = random.Random(seed)
    s = dict(n_cpdag=0, n_mpdag=0, D_all=0, N_cycle=0, N_vstruct=0, D_con=0,
             D_am=0, N_ident=0, cn_shrink=0, forb_shrink=0, pa_grow=0,
             partial_loss=0, N_set=0, L1_fail=0)
    viol = []
    for _ in range(n_dags):
        D = random_dag(p, prng.choice(degs), rng)
        C = dag_to_cpdag(D); ref = v_structures(C)
        s["n_cpdag"] += 1
        G = C
        for _step in range(walk_len):
            s["n_mpdag"] += 1
            # ---- L1: parents of a chain component are shared
            for (u, v) in undirected_edges(G):
                for c in range(p):
                    if is_directed(G, c, u) and not is_directed(G, c, v):
                        s["L1_fail"] += 1
            succ = []
            for (u, v) in undirected_edges(G):
                for (a, b) in ((u, v), (v, u)):
                    Hh = G.copy(); Hh[b, a] = 0
                    Hh = meek_closure(Hh)
                    bad = 0
                    if has_directed_cycle(Hh): bad = 1
                    elif v_structures(Hh) != ref: bad = 2
                    succ.append(((a, b), Hh, bad))
            for x in range(p):
                for y in range(p):
                    if x == y: continue
                    O0, _, a0 = X.ostar_and_paths(G, x, y)
                    if not a0: continue
                    cn0 = causal_nodes(G, x, y)
                    pa0 = parents_of_set(G, cn0); fb0 = poss_de(G, cn0) | {x}
                    for (ab, Hh, bad) in succ:
                        s["D_all"] += 1
                        if bad == 1: s["N_cycle"] += 1; continue
                        if bad == 2: s["N_vstruct"] += 1; continue
                        s["D_con"] += 1
                        O1, _, a1 = X.ostar_and_paths(Hh, x, y)
                        if not a1:
                            s["N_ident"] += 1; continue
                        s["D_am"] += 1
                        cn1 = causal_nodes(Hh, x, y)
                        pa1 = parents_of_set(Hh, cn1); fb1 = poss_de(Hh, cn1) | {x}
                        if cn1 != cn0: s["cn_shrink"] += 1
                        if fb1 != fb0: s["forb_shrink"] += 1
                        if pa1 - pa0: s["pa_grow"] += 1
                        for u in (pa0 - fb0):
                            ch0 = {w for w in cn0 if is_directed(G, u, w)}
                            ch1 = {w for w in cn1 if is_directed(Hh, u, w)}
                            if ch1 < ch0:
                                s["partial_loss"] += 1
                        if O1 != O0:
                            s["N_set"] += 1
                            if len(viol) < 5:
                                viol.append(dict(p=p, x=x, y=y, e=list(ab),
                                                 G0=G.tolist(), H=Hh.tolist(),
                                                 O0=sorted(O0), O1=sorted(O1)))
            ok = [t for t in succ if t[2] == 0]
            if not ok: break
            G = ok[prng.randrange(len(ok))][1]
    return s, viol


def main():
    p = int(sys.argv[1]); degs = [float(d) for d in sys.argv[2].split(",")]
    n_jobs = int(sys.argv[3]); n_dags = int(sys.argv[4]); seed0 = int(sys.argv[5])
    walk = int(sys.argv[6]) if len(sys.argv) > 6 else 8
    nw = 12
    tot = {}; viols = []
    with Pool(nw) as pool:
        for s, v in pool.imap_unordered(job, [(p, degs, seed0 + i, n_dags, walk)
                                              for i in range(n_jobs)]):
            for k, val in s.items(): tot[k] = tot.get(k, 0) + val
            viols.extend(v)
    print(json.dumps(dict(p=p, degs=degs, stats=tot, n_viol=len(viols),
                          viol=viols[:2]), default=float)[:3000], flush=True)


if __name__ == "__main__":
    main()
