"""
C-FIREWALL structural diagnostics -- the PROOF HANDLE, measured.

The synthesiser's handle was: "adding a->b can only GROW cn and forb.  Show that
every parent the addition contributes to cn_H is already in forb_H."

Both halves of that are testable by census, and neither is an assumption we may
make for free: Meek closure in H can orient an edge v_{i+1} -> v_i that was
undirected in G0, which DESTROYS a possibly-causal path and can therefore SHRINK
cn.  This script measures, over the same exhaustive census as cfirewall.py
(coherent statement, G0 amenable, H amenable):

  cn      : cn_H vs cn_G0                     -- superset / subset / equal / incomparable
  forb    : forb_H vs forb_G0                 -- same four-way
  O*      : O*_H vs O*_G0                     -- same four-way
  NEWPA   : pa_H(cn_H) \\ pa_G0(cn_G0), and how many of those nodes are in forb_H
            (the razor: `absorbed` = in forb_H, `escaped` = not, so it lands in O*_H)
  ESCAPE  : trials where O*_H \\ O*_G0 is non-empty -- a node the truth's own class
            never put in the reported set

Usage: python structure.py <p> <out.json> [nworkers]
"""
import json
import sys
import time
import importlib.util
from collections import defaultdict

import numpy as np

BASE = "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious"
sys.path.insert(0, BASE + "/x2/code")
from graphs import (dag_to_cpdag, has_directed_cycle, skeleton, v_structures)
import adjust
import x2lib as X
from run_lemma import all_dags, reachable_mpdags

_spec = importlib.util.spec_from_file_location("x1_ops", BASE + "/x1/code/x1_ops.py")
x1_ops = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(x1_ops)


def rel(A, B):
    """Relation of A to B."""
    if A == B:
        return "equal"
    if A > B:
        return "superset"
    if A < B:
        return "subset"
    return "incomparable"


def parts(G, x, y):
    """(cn, forb, pa(cn), O*, amenable)."""
    paths = X.pcp_capped(G, x, y)
    if not paths:
        return None
    if not all(adjust.is_directed(G, q[0], q[1]) for q in paths):
        return None
    cn = set()
    for q in paths:
        cn.update(q[1:])
    fb = adjust.poss_de(G, cn) | {x}
    pa = adjust.parents_of_set(G, cn)
    return frozenset(cn), frozenset(fb), frozenset(pa), frozenset(pa - fb)


def scan(C):
    acc = defaultdict(int)
    for G0 in reachable_mpdags(C):
        for k, v in scan_one(G0, C).items():
            acc[k] += v
    return dict(acc)


def scan_one(G0, C):
    p = C.shape[0]
    acc = defaultdict(int)
    NA = x1_ops.nonadjacent_pairs(G0)
    if NA:
        Q = []
        for x in range(p):
            for y in range(p):
                if x != y:
                    z = parts(G0, x, y)
                    if z is not None:
                        Q.append((x, y, z))
        for (a0, b0) in (NA if Q else []):
            for (a, b) in ((a0, b0), (b0, a0)):
                H, info = x1_ops.bk_assert(G0, [(a, b)])
                if info["conflict"] or has_directed_cycle(H) or not x1_ops.pdag_extendable(H):
                    continue
                for (x, y, (cn0, fb0, pa0, O0)) in Q:
                    z = parts(H, x, y)
                    if z is None:
                        continue
                    cn1, fb1, pa1, O1 = z
                    acc["n"] += 1
                    if O1 & fb0:
                        acc["F1_VIOLATED"] += 1     # O*_H meets forb_G0
                    acc["cn_" + rel(cn1, cn0)] += 1
                    acc["forb_" + rel(fb1, fb0)] += 1
                    acc["pa_" + rel(pa1, pa0)] += 1
                    acc["ostar_" + rel(O1, O0)] += 1
                    newpa = pa1 - pa0
                    acc["newpa_nodes"] += len(newpa)
                    acc["newpa_absorbed"] += len(newpa & fb1)
                    acc["newpa_escaped"] += len(newpa - fb1)
                    if newpa - fb1:
                        acc["trials_with_escaped_newpa"] += 1
                    esc = O1 - O0
                    acc["ostar_new_nodes"] += len(esc)
                    if esc:
                        acc["trials_ostar_grew"] += 1
                        # is the escaped node forbidden in the TRUTH's own graph?
                        acc["escaped_in_forb_G0"] += len(esc & fb0)
                    if O1 - O0 and O0 - O1:
                        acc["trials_ostar_swapped"] += 1
    return dict(acc)


def _job(args):
    kb, p = args
    return scan(np.frombuffer(kb, dtype=np.int8).reshape(p, p).copy())


def _sblock(args):
    """Sample mode: independent (C, G0) draws, same inner sweep."""
    import cfirewall as CF
    p, n, seed = args
    from graphs import random_dag
    rng = np.random.default_rng(seed)
    tot = defaultdict(int)
    for _ in range(n):
        deg = float(rng.choice([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0]))
        C = dag_to_cpdag(random_dag(p, deg, rng))
        ref = v_structures(C)
        G0 = CF.random_mpdag(C, ref, rng)
        for k, v in scan_one(G0, C).items():
            tot[k] += v
    return dict(tot)


def main():
    p = int(sys.argv[1])
    out_path = sys.argv[2]
    nw = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    ndraw = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    tot = defaultdict(int)
    t0 = time.time()
    from multiprocessing import Pool
    if ndraw:
        jobs = [(p, max(1, ndraw // 200), 5000 + 131 * i) for i in range(200)]
        fn, ck = _sblock, 1
    else:
        cpdags = {}
        for D in all_dags(p):
            C = dag_to_cpdag(D)
            cpdags.setdefault(C.tobytes(), C)
        jobs = [(k, p) for k in cpdags]
        fn, ck = _job, 4
    with Pool(nw) as pool:
        for acc in pool.imap_unordered(fn, jobs, chunksize=ck):
            for k, v in acc.items():
                tot[k] += v
    out = dict(p=p, mode=("sample" if ndraw else "census"), n_draws=ndraw,
               n_units=len(jobs), elapsed_s=round(time.time() - t0, 1),
               counts={k: int(v) for k, v in sorted(tot.items())})
    json.dump(out, open(out_path, "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
