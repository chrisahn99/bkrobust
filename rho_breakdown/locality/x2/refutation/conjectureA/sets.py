"""Conjecture A with SET-VALUED queries (|X|>=1, |Y|>=1) -- the one regime the
x2 code never touches (ostar_and_paths is single-node X, single-node Y).

Henckel-Perkovic-Maathuis definitions, lifted to sets:
  proper possibly-causal path: from some x in X to some y in Y, possibly causal,
    only its FIRST node in X;
  cn(X,Y)  = nodes on such paths, minus X;
  forb     = poss_de(cn) u X;
  amenable = every such path starts with a directed edge out of X;
  O*       = pa(cn) \ forb.
"""
import sys, json, random
from itertools import combinations
from multiprocessing import Pool
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/code")
import numpy as np
from graphs import (meek_closure, v_structures, has_directed_cycle, random_dag,
                    undirected_edges, dag_to_cpdag, is_directed)
from adjust import possibly_directed_neighbors, poss_de, parents_of_set


def proper_pc_paths(G, Xs, Ys):
    Xs, Ys = set(Xs), set(Ys)
    out = []
    for x0 in Xs:
        def dfs(path, visited):
            v = path[-1]
            if v in Ys and len(path) > 1:
                out.append(list(path))
                return
            for w in possibly_directed_neighbors(G, v):
                if w in visited or w in Xs:
                    continue
                visited.add(w); path.append(w)
                dfs(path, visited)
                path.pop(); visited.remove(w)
        dfs([x0], {x0})
    return out


def ostar_sets(G, Xs, Ys):
    Xs = set(Xs)
    paths = proper_pc_paths(G, Xs, Ys)
    if not paths:
        return None, False
    if not all(is_directed(G, p[0], p[1]) for p in paths):
        return None, False
    cn = set()
    for p in paths:
        cn.update(p[1:])
    cn -= Xs
    fb = (poss_de(G, cn) | Xs) if cn else set(Xs)
    return frozenset(parents_of_set(G, cn) - fb), True


def job(args):
    p, seed, n, kx, ky = args
    rng = np.random.default_rng(seed); prng = random.Random(seed)
    s = dict(n_mpdag=0, D_all=0, N_cycle=0, N_vstruct=0, D_con=0, D_am=0,
             N_ident=0, N_set=0, n_queries=0)
    viol = []
    for _ in range(n):
        D0 = random_dag(p, prng.choice([1.5, 2.0, 2.5, 3.0, 4.0]), rng)
        C = dag_to_cpdag(D0); ref = v_structures(C); G = C
        for _step in range(6):
            s["n_mpdag"] += 1
            succ = []
            for (u, v) in undirected_edges(G):
                for (a, b) in ((u, v), (v, u)):
                    H = G.copy(); H[b, a] = 0; H = meek_closure(H)
                    bad = 1 if has_directed_cycle(H) else (2 if v_structures(H) != ref else 0)
                    succ.append(((a, b), H, bad))
            nodes = list(range(p))
            for Xs in combinations(nodes, kx):
                rest = [v for v in nodes if v not in Xs]
                for Ys in combinations(rest, ky):
                    s["n_queries"] += 1
                    O0, a0 = ostar_sets(G, Xs, Ys)
                    if not a0: continue
                    for (ab, H, bad) in succ:
                        s["D_all"] += 1
                        if bad == 1: s["N_cycle"] += 1; continue
                        if bad == 2: s["N_vstruct"] += 1; continue
                        s["D_con"] += 1
                        O1, a1 = ostar_sets(H, Xs, Ys)
                        if not a1: s["N_ident"] += 1; continue
                        s["D_am"] += 1
                        if O1 != O0:
                            s["N_set"] += 1
                            if len(viol) < 3:
                                viol.append(dict(p=p, X=list(Xs), Y=list(Ys), e=list(ab),
                                                 G0=G.tolist(), H=H.tolist(),
                                                 O0=sorted(O0), O1=sorted(O1)))
            ok = [t for t in succ if t[2] == 0]
            if not ok: break
            G = ok[prng.randrange(len(ok))][1]
    return (kx, ky, p), s, viol


if __name__ == "__main__":
    tasks = []
    sd = 5000
    for (kx, ky) in ((1, 2), (2, 1), (2, 2), (3, 1), (1, 3), (2, 3)):
        for p in (6, 7, 8):
            for i in range(4):
                sd += 1
                tasks.append((p, sd, 120, kx, ky))
    agg = {}; viols = []
    with Pool(12) as pool:
        for key, s, v in pool.imap_unordered(job, tasks):
            k = f"|X|={key[0]},|Y|={key[1]},p={key[2]}"
            a = agg.setdefault(k, {})
            for kk, vv in s.items(): a[kk] = a.get(kk, 0) + vv
            viols.extend(v)
    print(json.dumps(dict(agg=agg, n_viol=len(viols), viol=viols[:2]), default=float)[:4000])
