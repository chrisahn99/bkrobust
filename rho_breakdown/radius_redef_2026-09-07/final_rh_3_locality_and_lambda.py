#!/usr/bin/env python3
"""FINAL DEFINITION, part 3.

(a) Does lambda-ordering of the ladder buy anything?  top-1 vs random-1 vs worst-1.
(b) Is the width small where damage is absent?  Stratify by whether any elicited
    claim is incident to X or Y (the paper's own locality measurement).
(c) Cost at the calibrated position, not at the full ladder.
"""
import sys, itertools
import numpy as np
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag, random_sem
from harness import random_dag
from final_rh_1_ladder import Theta, hull, w, merge, ladder, TOL


def make_strat(rng, n=9, p=0.30, kmax=3, want_far=True):
    """A problem whose claims are all FAR from (all NEAR to) the query."""
    dag = random_dag(rng, n, p)
    if not dag.directed_edges:
        return None
    cp = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cp.undirected_edges)
    if not (2 <= len(und) <= 7):
        return None
    nodes = sorted(dag.nodes)
    cand = [(x, y) for x in nodes for y in nodes if x != y]
    rng.shuffle(cand)
    for x, y in cand:
        pool = [e for e in und
                if ((x not in e and y not in e) if want_far else (x in e or y in e))]
        if not pool:
            continue
        idx = rng.permutation(len(pool))[:min(kmax, len(pool))]
        K = []
        for t in idx:
            a, b = pool[t]
            K.append((a, b) if (a, b) in dag.directed_edges else (b, a))
        g0 = apply_orientations(cp, K)
        if g0 is None:
            continue
        if optimal_adjustment_set_mpdag(g0, x, y) is None:
            continue
        sem = random_sem(dag, rng)
        cache = {}
        blanket = hull(Theta(sem, cp, x, y, cache))
        if w(blanket) <= TOL:
            continue
        return dict(dag=dag, cpdag=cp, x=x, y=y, K=[tuple(e) for e in K], g0=g0,
                    sem=sem, blanket=blanket, cache=cache,
                    n_ext_blanket=len(enumerate_dag_extensions(cp)))
    return None


def stratum(seed, N, want_far, q=0.20):
    rng = np.random.default_rng(seed)
    rows, tried = [], 0
    while len(rows) < N and tried < 900000:
        tried += 1
        pr = make_strat(rng, want_far=want_far)
        if pr is None:
            continue
        K = pr["K"]
        flip = [bool(rng.random() < q) for _ in K]
        Kst = [((e[1], e[0]) if f else e) for e, f in zip(K, flip)]
        g0 = apply_orientations(pr["cpdag"], Kst)
        if g0 is None:
            continue
        if optimal_adjustment_set_mpdag(g0, pr["x"], pr["y"]) is None:
            continue
        base, entries, ncl = ladder(pr, Kst, g0, 1)      # singleton tier only
        if not entries:
            continue
        rows.append(dict(pr=pr, base=base, entries=entries, m=len(Kst),
                         b_true=sum(flip), t=pr["sem"].true_total_effect(pr["x"], pr["y"]),
                         rng_pick=int(rng.integers(len(entries)))))
    return rows


def summarise(rows, name):
    if not rows:
        print(f"  [{name}] empty")
        return
    bl = np.mean([w(r["pr"]["blanket"]) for r in rows])
    out = {}
    for tag, sel in (("top-1 (lambda)", lambda r: r["entries"][0]),
                     ("random-1", lambda r: r["entries"][r["rng_pick"]]),
                     ("worst-1 (min lambda)", lambda r: r["entries"][-1]),
                     ("all singletons", None)):
        covs, wids, nus = [], [], []
        for r in rows:
            if sel is None:
                iv = r["base"]
                for e in r["entries"]:
                    iv = merge(iv, e[1])
            else:
                iv = merge(r["base"], sel(r)[1])
            covs.append(iv[0] - 1e-8 <= r["t"] <= iv[1] + 1e-8)
            wids.append(w(iv))
            nus.append(w(iv) / w(r["pr"]["blanket"]))
        out[tag] = (np.mean(covs), np.mean(wids), np.mean(nus))
    print(f"  [{name}] n={len(rows)}  blanket width={bl:.4f}  "
          f"P(b_true>=1)={np.mean([r['b_true']>=1 for r in rows]):.2f}")
    print(f"      base R_0            cov={np.mean([r['base'][0]-1e-8<=r['t']<=r['base'][1]+1e-8 for r in rows]):.3f} "
          f"width={np.mean([w(r['base']) for r in rows]):.4f}")
    for tag, (c, wd, nu) in out.items():
        print(f"      {tag:<22} cov={c:.3f} width={wd:.4f} nu={nu:.3f}")
    # lambda sparsity in this stratum
    lam0 = [sum(1 for e in r["entries"] if e[2] < 1e-9) / len(r["entries"]) for r in rows]
    print(f"      lambda=0 fraction of claims: {np.mean(lam0):.3f}   "
          f"all-zero ranking on {100*np.mean([x>0.999 for x in lam0]):.1f}% of problems")


if __name__ == "__main__":
    print("=== (a)+(b) stratified by locality of the elicited claims, q=0.20 ===")
    summarise(stratum(20260907, 200, want_far=True), "FAR: no claim incident to X or Y")
    summarise(stratum(20260907, 200, want_far=False), "NEAR: every claim incident to X or Y")
