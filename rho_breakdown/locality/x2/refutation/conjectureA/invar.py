"""Attack the LEMMA that makes Conjecture A trivial:
   G amenable MPDAG  =>  O*(G) == O*(D) for every consistent DAG extension D.
Break this and Conjecture A very likely breaks with it.
Families: ER random, plus structured skeletons (complete, cycle, grid, star,
bipartite, chordal-heavy) whose CPDAGs have big undirected chain components.
"""
import sys, json, random
from multiprocessing import Pool
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/code")
import numpy as np
import x2lib as X
from graphs import (meek_closure, v_structures, has_directed_cycle, random_dag,
                    undirected_edges, dag_to_cpdag, consistent_dag_extensions)
from adjust import optimal_adjustment_set


def dag_from_skeleton(S, rng):
    p = S.shape[0]
    order = rng.permutation(p)
    pos = {int(v): i for i, v in enumerate(order)}
    D = np.zeros((p, p), dtype=np.int8)
    for i in range(p):
        for j in range(p):
            if i < j and S[i, j]:
                if pos[i] < pos[j]: D[i, j] = 1
                else: D[j, i] = 1
    return D


def skel(kind, p, rng):
    S = np.zeros((p, p), dtype=np.int8)
    if kind == "complete":
        S[:] = 1; np.fill_diagonal(S, 0)
    elif kind == "cycle":
        for i in range(p): S[i, (i + 1) % p] = S[(i + 1) % p, i] = 1
    elif kind == "star":
        for i in range(1, p): S[0, i] = S[i, 0] = 1
    elif kind == "grid":
        w = int(np.ceil(np.sqrt(p)))
        for i in range(p):
            if (i % w) + 1 < w and i + 1 < p: S[i, i + 1] = S[i + 1, i] = 1
            if i + w < p: S[i, i + w] = S[i + w, i] = 1
    elif kind == "bipart":
        h = p // 2
        for i in range(h):
            for j in range(h, p):
                if rng.random() < 0.7: S[i, j] = S[j, i] = 1
    elif kind == "twoclique":
        h = p // 2 + 1
        for i in range(h):
            for j in range(i + 1, h): S[i, j] = S[j, i] = 1
        for i in range(h - 1, p):
            for j in range(i + 1, p): S[i, j] = S[j, i] = 1
    return S


def job(args):
    kind, p, seed, n = args
    rng = np.random.default_rng(seed); prng = random.Random(seed)
    s = dict(n_graph=0, n_triples=0, n_mismatch=0, n_amen_q=0, max_undir=0)
    ex = []
    for _ in range(n):
        if kind == "er":
            D0 = random_dag(p, prng.choice([1.5, 2.0, 2.5, 3.0, 4.0, 5.0]), rng)
        else:
            D0 = dag_from_skeleton(skel(kind, p, prng), rng)
        C = dag_to_cpdag(D0); ref = v_structures(C)
        G = C
        for _step in range(6):
            s["n_graph"] += 1
            U = undirected_edges(G)
            s["max_undir"] = max(s["max_undir"], len(U))
            if len(U) <= 14:
                exts = consistent_dag_extensions(G, ref)
                for x in range(p):
                    for y in range(p):
                        if x == y: continue
                        O0, _, a0 = X.ostar_and_paths(G, x, y)
                        if not a0: continue
                        s["n_amen_q"] += 1
                        for D in exts:
                            OD = optimal_adjustment_set(D, x, y)
                            s["n_triples"] += 1
                            if OD != O0:
                                s["n_mismatch"] += 1
                                if len(ex) < 2:
                                    ex.append(dict(kind=kind, p=p, x=x, y=y,
                                                   G=G.tolist(), D=D.tolist(),
                                                   O_G=sorted(O0),
                                                   O_D=sorted(OD) if OD is not None else None))
            nxt = []
            for (u, v) in U:
                for (a, b) in ((u, v), (v, u)):
                    H = G.copy(); H[b, a] = 0; H = meek_closure(H)
                    if has_directed_cycle(H) or v_structures(H) != ref: continue
                    nxt.append(H)
            if not nxt: break
            G = nxt[prng.randrange(len(nxt))]
    return kind, p, s, ex


if __name__ == "__main__":
    tasks = []
    sd = 900
    for kind in ("er", "complete", "cycle", "star", "grid", "bipart", "twoclique"):
        for p in (6, 7, 8, 9, 10):
            for k in range(2):
                sd += 1
                tasks.append((kind, p, sd, 250))
    agg = {}
    exs = []
    with Pool(12) as pool:
        for kind, p, s, ex in pool.imap_unordered(job, tasks):
            key = f"{kind}-p{p}"
            a = agg.setdefault(key, dict(n_graph=0, n_triples=0, n_mismatch=0, n_amen_q=0, max_undir=0))
            for k2 in ("n_graph", "n_triples", "n_mismatch", "n_amen_q"): a[k2] += s[k2]
            a["max_undir"] = max(a["max_undir"], s["max_undir"])
            exs.extend(ex)
    print(json.dumps(dict(agg=agg, n_mismatch_total=sum(v["n_mismatch"] for v in agg.values()),
                          ex=exs[:2]), default=float)[:4000])
