"""
C-FIREWALL -- exhaustive census.

CONJECTURE (as stated).  Let G0 be an MPDAG amenable rel. (X,Y), {a,b} NON-ADJACENT
in skeleton(G0), H = bk_assert(G0, a->b) (add the edge, Meek-close, re-stamp).
If H passes the three INTRINSIC coherence checks
    (c1) no orientation conflict, (c2) no directed cycle, (c3) Dor-Tarsi extendable
then EITHER H is not amenable rel. (X,Y) (a visible abort), OR O*(X,Y,H) is a valid
adjustment set for (X,Y) in D -- for EVERY D in [G0], not merely for one sampled D.

This script decides that by CENSUS, not by sampling:

  every labelled DAG on p nodes  ->  CPDAG (deduped)
  every reachable MPDAG G0 = M(C,K)          (run_lemma.reachable_mpdags: exactly
                                              { M(C,K) : K Meek-consistent })
  every ordered query (X,Y), G0 amenable
  every NON-ADJACENT unordered pair {a,b} of skeleton(G0), BOTH orientations
  -> bk_assert -> three coherence checks -> if H amenable, test O*(X,Y,H)
     against EVERY D in [G0] = consistent_dag_extensions(G0).

[G0] is intrinsic to G0: a DAG D lies in [G0] iff skeleton(D)=skeleton(G0),
D contains every directed edge of G0, D is acyclic and v_structures(D) =
v_structures(G0) ( = v_structures(C), asserted at run time ).  That is exactly
graphs.consistent_dag_extensions(G0, ref_vstructs=v_structures(C)).

Usage:
  python cfirewall.py census <p> <out.json> [nworkers]
  python cfirewall.py sample <p> <n_dag_draws> <out.json> [nworkers] [seed]
"""
import json
import sys
import time
import importlib.util
from collections import defaultdict

import numpy as np

BASE = "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious"
sys.path.insert(0, BASE + "/x2/code")

from graphs import (consistent_dag_extensions, dag_to_cpdag, directed_edges,
                    has_directed_cycle, meek_closure, random_dag, skeleton,
                    undirected_edges, v_structures)
import adjust
import x2lib as X
from run_lemma import all_dags, reachable_mpdags

_spec = importlib.util.spec_from_file_location("x1_ops", BASE + "/x1/code/x1_ops.py")
x1_ops = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(x1_ops)
bk_assert = x1_ops.bk_assert
pdag_extendable = x1_ops.pdag_extendable
nonadjacent_pairs = x1_ops.nonadjacent_pairs

MAX_WITNESS = 200


# ------------------------------------------------------------------ one MPDAG
def scan_mpdag(G0, ref, C, acc, witnesses, vcache):
    p = G0.shape[0]

    NA = nonadjacent_pairs(G0)              # unordered non-adjacent pairs of skeleton(G0)
    if not NA:
        return

    # --- amenable ordered queries in G0
    Q = []
    for x in range(p):
        for y in range(p):
            if x == y:
                continue
            O0, np0, amen0 = X.ostar_and_paths(G0, x, y)
            acc["query_all"] += 1
            if amen0:
                Q.append((x, y, O0))
    acc["query_amenable"] += len(Q)
    if not Q:
        return

    # --- [G0]: EVERY DAG in the equivalence class
    ext = consistent_dag_extensions(G0, ref_vstructs=ref)
    acc["ext_total"] += len(ext)
    if not ext:                              # cannot happen for a reachable MPDAG
        acc["ext_empty"] += 1
        return

    acc["ext_hist_%d" % min(len(ext), 24)] += 1
    if len(ext) > acc["ext_max"]:
        acc["ext_max"] = len(ext)

    # --- the statements: both orientations of every non-adjacent pair
    Hs = []
    Hs_rej = []
    for (a0, b0) in NA:
        for (a, b) in ((a0, b0), (b0, a0)):
            acc["stmt_all"] += 1
            H, info = bk_assert(G0, [(a, b)])
            c1 = not info["conflict"]
            c2 = not has_directed_cycle(H)
            c3 = pdag_extendable(H)
            if not (c1 and c2 and c3):
                acc["stmt_rejected"] += 1
                acc["rej_conflict"] += (not c1)
                acc["rej_cycle"] += (not c2)
                acc["rej_notext"] += (not c3)
                Hs_rej.append((a, b, H))
                continue
            acc["stmt_coherent"] += 1
            # sanity: the operator really did grow the skeleton
            if np.array_equal(skeleton(H), skeleton(G0)):
                acc["skel_unchanged_BUG"] += 1
            Hs.append((a, b, H))

    # --- CONTROL ARM: the same predicate on graphs the coherence check REJECTED
    #     (edge stamped anyway).  Establishes that a 0/N treatment result is not
    #     the artefact of a detector that never fires.
    for (a, b, H) in Hs_rej:
        for (x, y, O0) in Q:
            O1, np1, amen1 = X.ostar_and_paths(H, x, y)
            acc["ctrl_trial"] += 1
            if not amen1:
                acc["ctrl_abort"] += 1
                continue
            acc["ctrl_report"] += 1
            if O1 != O0:
                acc["ctrl_moved"] += 1
            nbad = 0
            for D in ext:
                kk = (D.tobytes(), x, y, O1)
                v = vcache.get(kk)
                if v is None:
                    v = adjust.is_valid_adjustment_set(D, x, y, set(O1))
                    vcache[kk] = v
                if not v:
                    nbad += 1
            if nbad:
                acc["ctrl_violating"] += 1
                if O1 != O0:
                    acc["ctrl_violating_moved"] += 1
                if nbad == len(ext):
                    acc["ctrl_violating_allD"] += 1
                else:
                    acc["ctrl_violating_someD"] += 1

    if not Hs:
        return

    # --- (statement x query) trials.  Group by (x,y,O*) so the |[G0]| validity
    #     sweep is paid once per distinct reported set, not once per trial.
    cand = {}
    for (a, b, H) in Hs:
        for (x, y, O0) in Q:
            O1, np1, amen1 = X.ostar_and_paths(H, x, y)
            acc["trial"] += 1
            if not amen1:
                acc["abort"] += 1
                acc["abort_nopath" if np1 == 0 else "abort_unamen"] += 1
                continue
            acc["report"] += 1
            moved = (O1 != O0)
            if moved:
                acc["report_moved"] += 1
            else:
                acc["report_same"] += 1
            if len(ext) > 1:
                acc["report_multiD"] += 1
                if moved:
                    acc["report_multiD_moved"] += 1
            k = (x, y, O1)
            e = cand.get(k)
            if e is None:
                cand[k] = e = [0, 0, []]
            e[0] += 1
            e[1] += int(moved)
            if len(e[2]) < 4:
                e[2].append((a, b, bool(moved)))

    # --- validity of every reported O* in EVERY D of [G0]
    for (x, y, O1), (ntr, nmv, stmts) in cand.items():
        bad = []
        for D in ext:
            kk = (D.tobytes(), x, y, O1)
            v = vcache.get(kk)
            if v is None:
                v = adjust.is_valid_adjustment_set(D, x, y, set(O1))
                vcache[kk] = v
            if not v:
                bad.append(D)
        acc["checks"] += ntr * len(ext)
        if bad:
            acc["violating_trials"] += ntr
            acc["violating_trials_moved"] += nmv
            acc["violating_keys"] += 1
            # THE distinction the brief asks about: invalid in every D of the
            # class, or invalid only for some D (so a one-D sample would miss it)?
            if len(bad) == len(ext):
                acc["violating_keys_allD"] += 1
            else:
                acc["violating_keys_someD"] += 1
            if len(witnesses) < MAX_WITNESS:
                witnesses.append(dict(
                    p=int(p), C=C.tolist(), G0=G0.tolist(),
                    x=int(x), y=int(y),
                    Ostar_H=sorted(int(z) for z in O1),
                    n_class=len(ext), n_bad_D=len(bad),
                    bad_D=[D.tolist() for D in bad[:6]],
                    statements=[[int(a), int(b), mv] for (a, b, mv) in stmts],
                    n_trials=ntr, n_trials_moved=nmv))
        else:
            acc["clean_trials"] += ntr


def scan_cpdag(C):
    ref = v_structures(C)
    acc = defaultdict(int)
    witnesses = []
    vcache = {}
    mp = reachable_mpdags(C)
    acc["mpdag"] = len(mp)
    for G0 in mp:
        assert v_structures(G0) == ref
        scan_mpdag(G0, ref, C, acc, witnesses, vcache)
    return dict(acc), witnesses


def random_mpdag(C, ref, rng):
    """One uniform-ish draw from { M(C,K) : K Meek-consistent }, by the SAME
    random walk `reachable_mpdags` explores exhaustively: orient one undirected
    edge at a time, Meek-close, reject a directed cycle or a new v-structure.
    A walk of random length t >= 0 from C, so G0 = C (t=0) is in the support."""
    G = C.copy()
    t = int(rng.integers(0, len(undirected_edges(C)) + 1))
    for _ in range(t):
        U = undirected_edges(G)
        if not U:
            break
        u, v = U[int(rng.integers(len(U)))]
        a, b = (u, v) if rng.random() < 0.5 else (v, u)
        H = G.copy()
        H[b, a] = 0
        H = meek_closure(H)
        if has_directed_cycle(H) or v_structures(H) != ref:
            continue
        G = H
    return G


def scan_sample_block(args):
    """One worker block of `n` independent (C, G0) draws at size p."""
    p, n, seed, max_u = args
    rng = np.random.default_rng(seed)
    acc = defaultdict(int)
    witnesses = []
    vcache = {}
    for _ in range(n):
        deg = float(rng.choice([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0]))
        D = random_dag(p, deg, rng)
        C = dag_to_cpdag(D)
        ref = v_structures(C)
        G0 = random_mpdag(C, ref, rng)
        acc["draw"] += 1
        if len(undirected_edges(G0)) > max_u:      # [G0] enumeration too big
            acc["draw_skipped_bigclass"] += 1
            continue
        acc["mpdag"] += 1
        scan_mpdag(G0, ref, C, acc, witnesses, vcache)
        if len(vcache) > 400000:
            vcache.clear()
    return dict(acc), witnesses


def _job(args):
    kb, p = args
    C = np.frombuffer(kb, dtype=np.int8).reshape(p, p).copy()
    return scan_cpdag(C)


def main():
    mode = sys.argv[1]
    p = int(sys.argv[2])
    if mode == "census":
        out_path = sys.argv[3]
        nw = int(sys.argv[4]) if len(sys.argv) > 4 else 10
        cpdags = {}
        for D in all_dags(p):
            C = dag_to_cpdag(D)
            cpdags.setdefault(C.tobytes(), C)
        n_draw = None
        jobs = [(k, p) for k in cpdags]
        fn = _job
        unit = "cpdags"
        nunits = len(cpdags)
        print(f"[cfirewall p={p} census] {len(cpdags)} distinct CPDAGs", flush=True)
    else:
        n_draw = int(sys.argv[3])
        out_path = sys.argv[4]
        nw = int(sys.argv[5]) if len(sys.argv) > 5 else 10
        seed = int(sys.argv[6]) if len(sys.argv) > 6 else 20260823
        max_u = int(sys.argv[7]) if len(sys.argv) > 7 else 16
        nblk = 400
        per = max(1, n_draw // nblk)
        jobs = [(p, per, seed + 7919 * i, max_u) for i in range(nblk)]
        fn = scan_sample_block
        unit = "blocks"
        nunits = nblk
        print(f"[cfirewall p={p} SAMPLE] {nblk} blocks x {per} (C,G0) draws "
              f"= {nblk*per}", flush=True)

    tot = defaultdict(int)
    wit = []
    t0 = time.time()
    done = 0
    from multiprocessing import Pool
    with Pool(nw) as pool:
        for acc, w in pool.imap_unordered(fn, jobs, chunksize=1 if mode != "census" else 4):
            for kk, vv in acc.items():
                if kk == "ext_max":
                    tot[kk] = max(tot[kk], vv)
                else:
                    tot[kk] += vv
            wit.extend(w)
            done += 1
            if done % max(1, nunits // 40) == 0:
                print(f"  {done}/{nunits} {unit}  {time.time()-t0:.0f}s  "
                      f"trials={tot['trial']} viol={tot['violating_trials']}", flush=True)

    out = dict(conjecture="C-FIREWALL", p=p, mode=mode, n_dag_draws=n_draw,
               exhaustive=(mode == "census"),
               n_units=nunits, elapsed_s=round(time.time() - t0, 1),
               counts={k: int(v) for k, v in sorted(tot.items())},
               n_witnesses=len(wit), witnesses=wit)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out["counts"], indent=1))
    print(f"VIOLATIONS: {tot['violating_trials']} trials / {tot['violating_keys']} "
          f"(x,y,O*) keys   -> {out_path}", flush=True)


if __name__ == "__main__":
    main()
