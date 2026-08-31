"""
X2 EXHAUSTIVE SMALL-p EXISTENCE SEARCH -- the proof-substitute M2 asks for.

The question S5 poses is a pure graph question: does there exist an MPDAG G0
(= M(C,K) for SOME CPDAG C and SOME Meek-consistent knowledge set K), a query
(X,Y), and ONE further admissible statement s whose skeleton distance to {X,Y}
is >= 1, such that O*(X,Y) changes?

No SCM, no Sigma and no truth are needed: E_set, E_ident and E_gain are all
functions of (G0, x, y, s) alone.  So at small p the question can be settled by
ENUMERATION rather than by sampling:

  1. enumerate every labelled DAG on p nodes, map to CPDAG, dedupe;
  2. from each CPDAG, BFS over every reachable MPDAG (orient one undirected
     edge, Meek-close, reject a directed cycle or a new v-structure) -- this is
     exactly the set { M(C,K) : K consistent }, i.e. EVERY knowledge state;
  3. for every reachable G0, every ordered query, every undirected edge of G0
     and both of its orientations, classify the trial.

At p <= 5 this is complete.  At p = 6 the DAG count is 3.8M and step 3 is
~10^8 O* evaluations, so p = 6 is a random DAG SAMPLE, deduped by CPDAG, and is
labelled as a sample -- NOT an enumeration.  That deviation from M2's wording
is recorded in RESULTS.md.

Usage: python run_lemma.py <p> <exhaustive|sample> <n_sample> <out.json> [nw]
"""
import json
import sys
from collections import defaultdict
from itertools import product
from multiprocessing import Pool

import numpy as np

import x2lib as X
from graphs import (dag_to_cpdag, has_directed_cycle, is_undirected,
                    meek_closure, random_dag, skeleton, undirected_edges,
                    v_structures)


def all_dags(p):
    idx = [(i, j) for i in range(p) for j in range(i + 1, p)]
    for state in product([0, 1, 2], repeat=len(idx)):
        D = np.zeros((p, p), dtype=np.int8)
        for (i, j), s in zip(idx, state):
            if s == 1:
                D[i, j] = 1
            elif s == 2:
                D[j, i] = 1
        A = D.astype(bool)
        R = A.copy()
        for _ in range(p):
            R = R | (R @ A)
        if np.any(np.diag(R)):
            continue
        yield D


def reachable_mpdags(C, cap=0):
    """Every M(C,K): BFS over single-edge orientations closed under Meek, with
    the same cycle / v-structure guards as apply_background_knowledge.

    `cap` > 0 aborts (returns None) once more than `cap` distinct MPDAGs are
    reachable.  A CPDAG that trips the cap is SKIPPED WHOLE and counted -- never
    partially scanned, because a partial scan would drop states non-uniformly.
    The cap is 0 (off) for the exhaustive p <= 5 runs and is only used at p = 6,
    where the reachable set blows up; the skip count is reported."""
    ref = v_structures(C)
    seen = {C.tobytes(): C}
    frontier = [C]
    while frontier:
        if cap and len(seen) > cap:
            return None
        nxt = []
        for G in frontier:
            for (u, v) in undirected_edges(G):
                for (a, b) in ((u, v), (v, u)):
                    H = G.copy()
                    H[b, a] = 0
                    H = meek_closure(H)
                    if has_directed_cycle(H) or v_structures(H) != ref:
                        continue
                    kb = H.tobytes()
                    if kb not in seen:
                        seen[kb] = H
                        nxt.append(H)
        frontier = nxt
    return list(seen.values())


def scan_cpdag(C, cap=0):
    p = C.shape[0]
    ref = v_structures(C)
    acc = defaultdict(lambda: defaultdict(int))
    hits = []
    n_mpdag = 0
    mp = reachable_mpdags(C, cap)
    if mp is None:
        return {}, [], -1                      # skipped: tripped the MPDAG cap
    for G0 in mp:
        n_mpdag += 1
        U0 = undirected_edges(G0)
        if not U0:
            continue
        for x in range(p):
            for y in range(p):
                if x == y:
                    continue
                O0, np0, amen0 = X.ostar_and_paths(G0, x, y)
                dist = X.hop_dist_inf(G0, [x, y])
                for (u, v) in U0:
                    dmin, dmax = X.stmt_dist(dist, (u, v))
                    bmin = X.dbucket(dmin)
                    key = acc[bmin]
                    sides = {}
                    for (a, b) in ((u, v), (v, u)):
                        H = G0.copy()
                        H[b, a] = 0
                        H = meek_closure(H)
                        key["D_all"] += 1
                        if has_directed_cycle(H) or v_structures(H) != ref:
                            key["N_meek"] += 1
                            continue
                        key["D_con"] += 1
                        O1, np1, amen1 = X.ostar_and_paths(H, x, y)
                        sides[(a, b)] = (O1, amen1)
                        if amen0 and not amen1:
                            key["N_ident"] += 1
                            key["nopath" if np1 == 0 else "unamen"] += 1
                        elif amen1 and not amen0:
                            key["N_gain"] += 1
                            key["gain_paths" if np0 == 0 else "gain_orient"] += 1
                        elif amen0 and amen1:
                            key["D_am"] += 1
                            if O1 != O0:
                                key["N_set"] += 1
                                if dmin >= 1:
                                    hits.append(dict(
                                        kind="add", p=int(p), x=int(x), y=int(y),
                                        s=[int(a), int(b)], dmin=float(dmin),
                                        dmax=float(dmax),
                                        G0=G0.tolist(), C=C.tolist(),
                                        O0=sorted(int(z) for z in O0),
                                        O1=sorted(int(z) for z in O1)))
                    # ---- the REVERSAL operator (E1'/pilot): compare the two
                    # orientations of the SAME edge against each other, i.e.
                    # "the rest of K is fixed and the expert states e backwards".
                    # This is the operator the lemma's `add` sweep does NOT cover
                    # when the common base G0 is itself non-amenable.
                    if len(sides) == 2:
                        (Oa, ama), (Ob, amb) = sides[(u, v)], sides[(v, u)]
                        key["D_con_rev"] += 1
                        if ama and amb:
                            key["D_am_rev"] += 1
                            if Oa != Ob:
                                key["N_set_rev"] += 1
                                if dmin >= 1:
                                    hits.append(dict(
                                        kind="reverse", p=int(p), x=int(x), y=int(y),
                                        s=[int(u), int(v)], dmin=float(dmin),
                                        dmax=float(dmax),
                                        G0=G0.tolist(), C=C.tolist(),
                                        O0=sorted(int(z) for z in Oa),
                                        O1=sorted(int(z) for z in Ob)))
                        elif ama != amb:
                            key["N_ident_rev"] += 1
    return {b: dict(v) for b, v in acc.items()}, hits, n_mpdag


def _job(args):
    Cb, p, cap = args
    C = np.frombuffer(Cb, dtype=np.int8).reshape(p, p).copy()
    return scan_cpdag(C, cap)


def main():
    p = int(sys.argv[1])
    mode = sys.argv[2]
    n_sample = int(sys.argv[3])
    out_path = sys.argv[4]
    nw = int(sys.argv[5]) if len(sys.argv) > 5 else 12
    cap = int(sys.argv[6]) if len(sys.argv) > 6 else 0

    cpdags = {}
    if mode == "exhaustive":
        for D in all_dags(p):
            C = dag_to_cpdag(D)
            cpdags.setdefault(C.tobytes(), C)
    else:
        rng = np.random.default_rng(20260822 + p)
        for _ in range(n_sample):
            deg = float(rng.choice([1.0, 1.5, 2.0, 2.5, 3.0, 4.0]))
            D = random_dag(p, deg, rng)
            C = dag_to_cpdag(D)
            cpdags.setdefault(C.tobytes(), C)
    keys = list(cpdags.keys())
    print(f"[lemma p={p} {mode}] {len(keys)} distinct CPDAGs", flush=True)

    tot = defaultdict(lambda: defaultdict(int))
    hits = []
    n_mpdag = 0
    n_skip = 0
    with Pool(nw) as pool:
        for acc, h, nm in pool.imap_unordered(_job, [(k, p, cap) for k in keys], chunksize=1):
            if nm < 0:
                n_skip += 1
                continue
            n_mpdag += nm
            hits.extend(h)
            for b, d in acc.items():
                for kk, vv in d.items():
                    tot[b][kk] += vv
    out = dict(p=p, mode=mode, mpdag_cap=cap, n_cpdags_scanned=len(keys) - n_skip,
               n_cpdags_skipped=n_skip, n_cpdags=len(keys), n_mpdags=n_mpdag,
               by_dmin={b: dict(v) for b, v in sorted(tot.items())},
               n_hits_far=len(hits), hits=hits[:20])
    with open(out_path, "w") as f:
        json.dump(out, f, indent=1, default=float)
    print(f"[lemma p={p}] {len(keys)-n_skip}/{len(keys)} CPDAGs scanned "
          f"({n_skip} skipped by cap={cap}), {n_mpdag} MPDAGs, "
          f"far E_set hits = {len(hits)} -> {out_path}", flush=True)
    for b, v in sorted(tot.items()):
        print(f"   dmin={b:>4}: " + " ".join(f"{a}={c}" for a, c in sorted(v.items())),
              flush=True)


if __name__ == "__main__":
    main()
