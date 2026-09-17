"""EXHAUSTIVE p=6 scan of CONJECTURE A (the `add` operator only).

Stage 1: enumerate every labelled DAG on 6 nodes, map to CPDAG, dedupe.
Stage 2: for each CPDAG, BFS over EVERY reachable MPDAG; for every ordered
         query with G0 amenable, every undirected edge, both orientations,
         Meek-closed, guards checked -> compare O*(H) vs O*(G0).
"""
import sys, json, pickle, os
from itertools import combinations, product
from multiprocessing import Pool
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/code")
import numpy as np
import x2lib as X
from adjust import causal_nodes, poss_de, parents_of_set
from graphs import (meek_closure, v_structures, has_directed_cycle, is_directed,
                    undirected_edges, dag_to_cpdag)

P = 6
IDX = [(i, j) for i in range(P) for j in range(i + 1, P)]


def acyclic(A):
    R = A.astype(bool).copy()
    B = R.copy()
    for _ in range(P):
        R = R | (R @ B)
    return not np.any(np.diag(R))


def stage1_job(mask):
    """all DAGs whose skeleton is `mask` (bitmask over IDX) -> set of CPDAG bytes"""
    edges = [IDX[k] for k in range(15) if mask >> k & 1]
    out = set()
    for bits in product([0, 1], repeat=len(edges)):
        D = np.zeros((P, P), dtype=np.int8)
        for (i, j), b in zip(edges, bits):
            if b: D[i, j] = 1
            else: D[j, i] = 1
        if not acyclic(D):
            continue
        out.add(dag_to_cpdag(D).tobytes())
    return out


def reachable_mpdags(C):
    ref = v_structures(C)
    seen = {C.tobytes(): C}
    frontier = [C]
    while frontier:
        nxt = []
        for G in frontier:
            for (u, v) in undirected_edges(G):
                for (a, b) in ((u, v), (v, u)):
                    H = G.copy(); H[b, a] = 0
                    H = meek_closure(H)
                    if has_directed_cycle(H) or v_structures(H) != ref:
                        continue
                    kb = H.tobytes()
                    if kb not in seen:
                        seen[kb] = H; nxt.append(H)
        frontier = nxt
    return list(seen.values()), ref


def stage2_job(Cb):
    C = np.frombuffer(Cb, dtype=np.int8).reshape(P, P).copy()
    s = dict(n_mpdag=0, D_all=0, N_cycle=0, N_vstruct=0, D_con=0, D_am=0,
             N_ident=0, cn_shrink=0, forb_shrink=0, pa_grow=0, total_loss=0,
             partial_loss=0, N_set=0)
    viol = []
    mp, ref = reachable_mpdags(C)
    for G in mp:
        s["n_mpdag"] += 1
        U = undirected_edges(G)
        if not U: continue
        succ = []
        for (u, v) in U:
            for (a, b) in ((u, v), (v, u)):
                H = G.copy(); H[b, a] = 0
                H = meek_closure(H)
                bad = 1 if has_directed_cycle(H) else (2 if v_structures(H) != ref else 0)
                succ.append(((a, b), H, bad))
        for x in range(P):
            for y in range(P):
                if x == y: continue
                O0, _, a0 = X.ostar_and_paths(G, x, y)
                if not a0: continue
                cn0 = causal_nodes(G, x, y)
                pa0 = parents_of_set(G, cn0); fb0 = poss_de(G, cn0) | {x}
                for (ab, H, bad) in succ:
                    s["D_all"] += 1
                    if bad == 1: s["N_cycle"] += 1; continue
                    if bad == 2: s["N_vstruct"] += 1; continue
                    s["D_con"] += 1
                    O1, _, a1 = X.ostar_and_paths(H, x, y)
                    if not a1: s["N_ident"] += 1; continue
                    s["D_am"] += 1
                    cn1 = causal_nodes(H, x, y)
                    pa1 = parents_of_set(H, cn1); fb1 = poss_de(H, cn1) | {x}
                    if cn1 != cn0: s["cn_shrink"] += 1
                    if fb1 != fb0: s["forb_shrink"] += 1
                    if pa1 - pa0: s["pa_grow"] += 1
                    for u in (pa0 - fb0):
                        ch0 = {w for w in cn0 if is_directed(G, u, w)}
                        ch1 = {w for w in cn1 if is_directed(H, u, w)}
                        if ch1 < ch0:
                            s["partial_loss"] += 1
                            if not ch1: s["total_loss"] += 1
                    if O1 != O0:
                        s["N_set"] += 1
                        if len(viol) < 3:
                            viol.append(dict(x=x, y=y, e=list(ab), G0=G.tolist(),
                                             H=H.tolist(), O0=sorted(O0), O1=sorted(O1)))
    return s, viol


def main():
    cache = "/private/tmp/claude-501/-Users-josecosta-mugango/bceac0a4-93a7-4757-b3b3-a3d4e6c42d55/scratchpad/refute/cpdags6.pkl"
    if os.path.exists(cache):
        cps = pickle.load(open(cache, "rb"))
    else:
        cps = set()
        with Pool(12) as pool:
            for out in pool.imap_unordered(stage1_job, range(1 << 15), chunksize=64):
                cps |= out
        cps = sorted(cps)
        pickle.dump(cps, open(cache, "wb"))
    print(f"[exh6] {len(cps)} distinct CPDAGs", flush=True)
    tot = {}; viols = []; done = 0
    with Pool(12) as pool:
        for s, v in pool.imap_unordered(stage2_job, cps, chunksize=16):
            for k, val in s.items(): tot[k] = tot.get(k, 0) + val
            viols.extend(v); done += 1
            if done % 20000 == 0:
                print(f"[exh6] {done}/{len(cps)}  N_set={tot.get('N_set',0)} "
                      f"forb_shrink={tot.get('forb_shrink',0)} D_am={tot.get('D_am',0)}", flush=True)
    out = dict(p=P, n_cpdags=len(cps), stats=tot, n_viol=len(viols), viol=viols[:3])
    json.dump(out, open("/private/tmp/claude-501/-Users-josecosta-mugango/bceac0a4-93a7-4757-b3b3-a3d4e6c42d55/scratchpad/refute/exh6_result.json", "w"), indent=1, default=float)
    print(json.dumps({k: v for k, v in out.items() if k != "viol"}, default=float), flush=True)


if __name__ == "__main__":
    main()
