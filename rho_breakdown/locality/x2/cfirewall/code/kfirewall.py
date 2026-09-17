"""
C-FIREWALL with k >= 2 spurious statements -- does the firewall SURVIVE a second
false required edge?

C-FIREWALL as stated is a ONE-statement theorem, and it does NOT self-induct:
after the first assertion H has a skeleton no DAG of the class has, so H is not
an MPDAG of any CPDAG containing D and the theorem cannot be re-applied to it.
Whether the property nevertheless survives is therefore an independent question,
and it is the question a practitioner actually faces -- b-LOAD stamps a whole
background-knowledge matrix, not one edge.

Same census skeleton as cfirewall.py, but the operator is an ORDERED SEQUENCE of
k statements on DISTINCT non-adjacent pairs of skeleton(G0) (order matters: bk_assert
re-stamps the whole applied prefix after each Meek closure).  The three coherence
checks are applied to the FINAL graph only, exactly as a tool would.

Usage: python kfirewall.py <census|sample> <p> <k> ... (see cfirewall.py)
"""
import json
import sys
import time
import importlib.util
from collections import defaultdict
from itertools import permutations

import numpy as np

BASE = "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious"
sys.path.insert(0, BASE + "/x2/code")
sys.path.insert(0, BASE + "/x2/cfirewall/code")
from graphs import (consistent_dag_extensions, dag_to_cpdag, has_directed_cycle,
                    random_dag, skeleton, undirected_edges, v_structures)
import adjust
import x2lib as X
from run_lemma import all_dags, reachable_mpdags
import cfirewall as CF

x1_ops = CF.x1_ops
K = 2


def scan_mpdag_k(G0, ref, C, acc, witnesses, vcache, k, rng=None, nstmt=None):
    p = G0.shape[0]
    NA = x1_ops.nonadjacent_pairs(G0)
    if len(NA) < k:
        return
    Q = []
    for x in range(p):
        for y in range(p):
            if x == y:
                continue
            O0, np0, amen0 = X.ostar_and_paths(G0, x, y)
            if amen0:
                Q.append((x, y, O0))
    if not Q:
        return
    ext = consistent_dag_extensions(G0, ref_vstructs=ref)
    if not ext:
        return

    # ordered sequences of k distinct pairs, each independently oriented
    ordered = [(a, b) for (a0, b0) in NA for (a, b) in ((a0, b0), (b0, a0))]
    seqs = []
    if nstmt is None:
        for seq in permutations(range(len(ordered)), k):
            if len({frozenset(ordered[i]) for i in seq}) < k:
                continue
            seqs.append([ordered[i] for i in seq])
    else:
        for _ in range(nstmt):
            idx = rng.choice(len(NA), size=k, replace=False)
            seq = []
            for t in idx:
                a, b = NA[int(t)]
                seq.append((a, b) if rng.random() < 0.5 else (b, a))
            seqs.append(seq)

    cand = {}
    for seq in seqs:
        H, info = x1_ops.bk_assert(G0, seq)
        acc["stmt_all"] += 1
        if info["conflict"] or has_directed_cycle(H) or not x1_ops.pdag_extendable(H):
            acc["stmt_rejected"] += 1
            continue
        acc["stmt_coherent"] += 1
        for (x, y, O0) in Q:
            O1, np1, amen1 = X.ostar_and_paths(H, x, y)
            acc["trial"] += 1
            if not amen1:
                acc["abort"] += 1
                continue
            acc["report"] += 1
            moved = (O1 != O0)
            acc["report_moved" if moved else "report_same"] += 1
            kk = (x, y, O1)
            e = cand.get(kk)
            if e is None:
                cand[kk] = e = [0, 0, []]
            e[0] += 1
            e[1] += int(moved)
            if len(e[2]) < 3:
                e[2].append([list(map(int, s)) for s in seq])

    for (x, y, O1), (ntr, nmv, seqs_w) in cand.items():
        bad = []
        for D in ext:
            key = (D.tobytes(), x, y, O1)
            v = vcache.get(key)
            if v is None:
                v = adjust.is_valid_adjustment_set(D, x, y, set(O1))
                vcache[key] = v
            if not v:
                bad.append(D)
        acc["checks"] += ntr * len(ext)
        if bad:
            acc["violating_trials"] += ntr
            acc["violating_trials_moved"] += nmv
            acc["violating_keys"] += 1
            acc["violating_keys_allD" if len(bad) == len(ext) else "violating_keys_someD"] += 1
            if len(witnesses) < 60:
                witnesses.append(dict(p=int(p), C=C.tolist(), G0=G0.tolist(),
                                      x=int(x), y=int(y),
                                      Ostar_H=sorted(int(z) for z in O1),
                                      n_class=len(ext), n_bad_D=len(bad),
                                      bad_D=[D.tolist() for D in bad[:4]],
                                      statement_seqs=seqs_w,
                                      n_trials=ntr, n_trials_moved=nmv))
        else:
            acc["clean_trials"] += ntr


def _job(args):
    kb, p, k = args
    C = np.frombuffer(kb, dtype=np.int8).reshape(p, p).copy()
    ref = v_structures(C)
    acc, wit, vc = defaultdict(int), [], {}
    for G0 in reachable_mpdags(C):
        scan_mpdag_k(G0, ref, C, acc, wit, vc, k)
    return dict(acc), wit


def _sblock(args):
    p, k, n, seed, max_u, nstmt = args
    rng = np.random.default_rng(seed)
    acc, wit, vc = defaultdict(int), [], {}
    for _ in range(n):
        deg = float(rng.choice([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0]))
        C = dag_to_cpdag(random_dag(p, deg, rng))
        ref = v_structures(C)
        G0 = CF.random_mpdag(C, ref, rng)
        if len(undirected_edges(G0)) > max_u:
            acc["draw_skipped_bigclass"] += 1
            continue
        acc["draw"] += 1
        scan_mpdag_k(G0, ref, C, acc, wit, vc, k, rng, nstmt)
        if len(vc) > 400000:
            vc.clear()
    return dict(acc), wit


def main():
    mode, p, k = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    if mode == "census":
        out_path, nw = sys.argv[4], int(sys.argv[5])
        cp = {}
        for D in all_dags(p):
            C = dag_to_cpdag(D)
            cp.setdefault(C.tobytes(), C)
        jobs = [(kb, p, k) for kb in cp]
        fn, ck = _job, 4
        n_draw = None
    else:
        n_draw, out_path, nw = int(sys.argv[4]), sys.argv[5], int(sys.argv[6])
        seed = int(sys.argv[7]) if len(sys.argv) > 7 else 424242
        nstmt = int(sys.argv[8]) if len(sys.argv) > 8 else 40
        nblk = 200
        jobs = [(p, k, max(1, n_draw // nblk), seed + 7919 * i, 14, nstmt)
                for i in range(nblk)]
        fn, ck = _sblock, 1
    print(f"[kfirewall p={p} k={k} {mode}] {len(jobs)} units", flush=True)
    tot, wit = defaultdict(int), []
    t0 = time.time()
    from multiprocessing import Pool
    done = 0
    with Pool(nw) as pool:
        for acc, w in pool.imap_unordered(fn, jobs, chunksize=ck):
            for kk, vv in acc.items():
                tot[kk] += vv
            wit.extend(w)
            done += 1
            if done % max(1, len(jobs) // 20) == 0:
                print(f"  {done}/{len(jobs)}  {time.time()-t0:.0f}s trials={tot['trial']} "
                      f"viol={tot['violating_trials']}", flush=True)
    out = dict(conjecture="C-FIREWALL k-statement extension", p=p, k=k, mode=mode,
               exhaustive=(mode == "census"), n_dag_draws=n_draw,
               elapsed_s=round(time.time() - t0, 1),
               counts={kk: int(vv) for kk, vv in sorted(tot.items())},
               n_witnesses=len(wit), witnesses=wit)
    json.dump(out, open(out_path, "w"), indent=1)
    print(json.dumps(out["counts"], indent=1))
    print(f"VIOLATIONS: {tot['violating_trials']} trials / {tot['violating_keys']} keys")


if __name__ == "__main__":
    main()
