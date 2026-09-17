"""GUARD ABLATION from genuine MPDAG bases.
Which hypothesis of Conjecture A is load-bearing?  Drop one at a time and see
whether O* moves / forb shrinks.
modes: full | nomeek | novstruct | nomeek_novstruct
"""
import sys, json, random
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/code")
import numpy as np
import x2lib as X
from adjust import causal_nodes, poss_de, parents_of_set
from graphs import (meek_closure, v_structures, has_directed_cycle,
                    undirected_edges, dag_to_cpdag, random_dag)


def succ_full(G, ref):
    out = []
    for (u, v) in undirected_edges(G):
        for (a, b) in ((u, v), (v, u)):
            H = G.copy(); H[b, a] = 0
            H = meek_closure(H)
            if has_directed_cycle(H) or v_structures(H) != ref:
                continue
            out.append(((a, b), H))
    return out


def succ_mode(G, ref, mode):
    out = []
    for (u, v) in undirected_edges(G):
        for (a, b) in ((u, v), (v, u)):
            H = G.copy(); H[b, a] = 0
            if "nomeek" not in mode:
                H = meek_closure(H)
            if has_directed_cycle(H):
                continue
            if "novstruct" not in mode and v_structures(H) != ref:
                continue
            out.append(((a, b), H))
    return out


def run(mode, p, n_dags, seed, degs):
    rng = np.random.default_rng(seed); prng = random.Random(seed)
    s = dict(D_am=0, forb_shrink=0, cn_shrink=0, N_set=0, N_gain_node=0, N_loss_node=0)
    ex = []
    for _ in range(n_dags):
        D = random_dag(p, prng.choice(degs), rng)
        C = dag_to_cpdag(D); ref = v_structures(C)
        G = C
        for _step in range(6):
            walk = succ_full(G, ref)          # base states are ALWAYS legit MPDAGs
            S = succ_mode(G, ref, mode)
            for x in range(p):
                for y in range(p):
                    if x == y: continue
                    O0, _, a0 = X.ostar_and_paths(G, x, y)
                    if not a0: continue
                    cn0 = causal_nodes(G, x, y); fb0 = poss_de(G, cn0) | {x}
                    for (ab, H) in S:
                        O1, _, a1 = X.ostar_and_paths(H, x, y)
                        if not a1: continue
                        s["D_am"] += 1
                        cn1 = causal_nodes(H, x, y); fb1 = poss_de(H, cn1) | {x}
                        if fb1 != fb0: s["forb_shrink"] += 1
                        if cn1 != cn0: s["cn_shrink"] += 1
                        if O1 != O0:
                            s["N_set"] += 1
                            s["N_gain_node"] += len(O1 - O0) > 0
                            s["N_loss_node"] += len(O0 - O1) > 0
                            if len(ex) < 2:
                                ex.append(dict(x=x, y=y, e=list(ab), G0=G.tolist(),
                                               H=H.tolist(), O0=sorted(O0), O1=sorted(O1)))
            if not walk: break
            G = walk[prng.randrange(len(walk))][1]
    return s, ex


if __name__ == "__main__":
    p = int(sys.argv[1]); n = int(sys.argv[2]); seed = int(sys.argv[3])
    degs = [float(d) for d in sys.argv[4].split(",")]
    for mode in ("full", "novstruct", "nomeek", "nomeek_novstruct"):
        s, ex = run(mode, p, n, seed, degs)
        print(json.dumps(dict(mode=mode, p=p, stats=s, ex=ex[:1]), default=float)[:2200])
