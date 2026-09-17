"""REFUTER for CONJECTURE A.

G0 MPDAG (reachable from a CPDAG by Meek-consistent single-edge orientations),
amenable for (x,y).  Orient one undirected edge e, Meek-close -> H.
Guards: no directed cycle, v_structures(H) == v_structures(G0), H amenable.
Claim: O*(H) == O*(G0).   We hunt for O*(H) != O*(G0).
"""
import sys, json, itertools, random
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/code")
import numpy as np
import x2lib as X
from graphs import (meek_closure, v_structures, has_directed_cycle,
                    undirected_edges, dag_to_cpdag, random_dag, skeleton)


def reachable_step(G, ref):
    """all one-edge Meek-consistent successors of G"""
    out = []
    for (u, v) in undirected_edges(G):
        for (a, b) in ((u, v), (v, u)):
            H = G.copy(); H[b, a] = 0
            H = meek_closure(H)
            if has_directed_cycle(H) or v_structures(H) != ref:
                continue
            out.append(((a, b), H))
    return out


def scan_mpdag(G0, ref, hits, stats, tag):
    """test conjecture A on every (x,y) and every undirected edge of G0"""
    p = G0.shape[0]
    succ = reachable_step(G0, ref)
    if not succ:
        return
    for x in range(p):
        for y in range(p):
            if x == y:
                continue
            O0, np0, amen0 = X.ostar_and_paths(G0, x, y)
            if not amen0:
                continue
            dist = X.hop_dist_inf(G0, [x, y])
            for (ab, H) in succ:
                stats["D_con"] += 1
                O1, np1, amen1 = X.ostar_and_paths(H, x, y)
                if not amen1:
                    continue
                stats["D_am"] += 1
                if O1 != O0:
                    stats["N_set"] += 1
                    dmin, dmax = X.stmt_dist(dist, ab)
                    hits.append(dict(tag=tag, p=int(p), x=x, y=y, e=list(ab),
                                     dmin=float(dmin),
                                     G0=G0.tolist(), H=H.tolist(),
                                     O0=sorted(O0), O1=sorted(O1)))


def random_walk_mpdags(C, rng, n_states):
    """random walk over reachable MPDAGs; return list of visited states"""
    ref = v_structures(C)
    states = [C]
    G = C
    for _ in range(n_states):
        succ = reachable_step(G, ref)
        if not succ:
            break
        G = succ[rng.randrange(len(succ))][1]
        states.append(G)
    return states, ref


def main():
    p = int(sys.argv[1]); n_dags = int(sys.argv[2]); seed = int(sys.argv[3])
    degs = [float(d) for d in sys.argv[4].split(",")]
    rng = np.random.default_rng(seed)
    prng = random.Random(seed)
    hits = []
    stats = dict(D_con=0, D_am=0, N_set=0, n_cpdag=0, n_mpdag=0)
    seen = set()
    for it in range(n_dags):
        deg = degs[prng.randrange(len(degs))]
        D = random_dag(p, deg, rng)
        C = dag_to_cpdag(D)
        kb = C.tobytes()
        if kb in seen:
            continue
        seen.add(kb)
        stats["n_cpdag"] += 1
        states, ref = random_walk_mpdags(C, prng, 8)
        for G0 in states:
            stats["n_mpdag"] += 1
            scan_mpdag(G0, ref, hits, stats, f"rand-p{p}")
        if hits:
            break
    print(json.dumps(dict(p=p, seed=seed, degs=degs, stats=stats,
                          n_hits=len(hits), hits=hits[:3]), default=float))


if __name__ == "__main__":
    main()
