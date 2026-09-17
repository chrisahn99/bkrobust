"""BOUNDARY TEST: does Conjecture A survive when G0 is NOT a genuine MPDAG but
the graph an analyst actually holds after SPURIOUS background knowledge
(x1's `bk_assert`: assert i->j on a NON-ADJACENT pair, growing the skeleton)?

G0 := meek_closure(C + spurious required edges).  Reference v-structures are
taken to be v_structures(G0) (the analyst has no other reference).
Then: orient one undirected edge of G0, Meek-close, keep if no cycle and no new
v-structure and both graphs amenable.  Compare O*.
"""
import sys, json, random
from multiprocessing import Pool
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/code")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x1/code")
import numpy as np
import x2lib as X
from x1_ops import bk_assert, nonadjacent_pairs, pdag_extendable
from graphs import (meek_closure, v_structures, has_directed_cycle, random_dag,
                    undirected_edges, dag_to_cpdag)
from adjust import causal_nodes, poss_de, parents_of_set


def job(args):
    p, seed, n, nspur = args
    rng = np.random.default_rng(seed); prng = random.Random(seed)
    s = dict(n_G0=0, n_extendable=0, D_all=0, N_cycle=0, N_vstruct=0, D_con=0,
             D_am=0, N_ident=0, forb_shrink=0, cn_shrink=0, N_set=0)
    viol = []
    for _ in range(n):
        D0 = random_dag(p, prng.choice([1.5, 2.0, 2.5, 3.0]), rng)
        C = dag_to_cpdag(D0)
        na = nonadjacent_pairs(C)
        if len(na) < nspur: continue
        pairs = prng.sample(na, nspur)
        K = [(a, b) if prng.random() < .5 else (b, a) for (a, b) in pairs]
        G0, info = bk_assert(C, K)
        if has_directed_cycle(G0): continue
        s["n_G0"] += 1
        ext = pdag_extendable(G0)
        s["n_extendable"] += bool(ext)
        ref = v_structures(G0)
        succ = []
        for (u, v) in undirected_edges(G0):
            for (a, b) in ((u, v), (v, u)):
                H = G0.copy(); H[b, a] = 0; H = meek_closure(H)
                bad = 1 if has_directed_cycle(H) else (2 if v_structures(H) != ref else 0)
                succ.append(((a, b), H, bad))
        for x in range(p):
            for y in range(p):
                if x == y: continue
                O0, _, a0 = X.ostar_and_paths(G0, x, y)
                if not a0: continue
                cn0 = causal_nodes(G0, x, y); fb0 = poss_de(G0, cn0) | {x}
                for (ab, H, bad) in succ:
                    s["D_all"] += 1
                    if bad == 1: s["N_cycle"] += 1; continue
                    if bad == 2: s["N_vstruct"] += 1; continue
                    s["D_con"] += 1
                    O1, _, a1 = X.ostar_and_paths(H, x, y)
                    if not a1: s["N_ident"] += 1; continue
                    s["D_am"] += 1
                    cn1 = causal_nodes(H, x, y); fb1 = poss_de(H, cn1) | {x}
                    if cn1 != cn0: s["cn_shrink"] += 1
                    if fb1 != fb0: s["forb_shrink"] += 1
                    if O1 != O0:
                        s["N_set"] += 1
                        if len(viol) < 3:
                            viol.append(dict(p=p, x=x, y=y, e=list(ab), K=[list(k) for k in K],
                                             C=C.tolist(), G0=G0.tolist(), H=H.tolist(),
                                             O0=sorted(O0), O1=sorted(O1),
                                             G0_extendable=bool(ext)))
    return s, viol


if __name__ == "__main__":
    tasks = [(p, 3000 + 17 * i + p, 400, ns)
             for p in (5, 6, 7, 8) for ns in (1, 2, 3) for i in range(4)]
    tot = {}; viols = []
    with Pool(12) as pool:
        for s, v in pool.imap_unordered(job, tasks):
            for k, val in s.items(): tot[k] = tot.get(k, 0) + val
            viols.extend(v)
    print(json.dumps(dict(stats=tot, n_viol=len(viols), viol=viols[:2]), default=float)[:4000])
