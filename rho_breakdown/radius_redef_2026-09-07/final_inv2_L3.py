#!/usr/bin/env python3
"""FINAL inverse-2, part 3: ABLATION. Does L add anything to P for the width, or is
the whole saving the proximity criterion? Also: does P add anything to L?"""
import sys
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from bkrobust.demo.evaluate import (optimal_adjustment_set_dag, adjusted_estimand,
                                    is_valid_adjustment_set_dag, random_sem)
from bkrobust.demo.meek import enumerate_dag_extensions
import final_inv2_L as F
from final_inv2_L import partition, admissible, pairs_revisions

def run(p_err, n_target, seed):
    F.P_ERR = p_err
    rng = np.random.default_rng(seed); R, tried = [], 0
    while len(R) < n_target and tried < 80000:
        tried += 1
        pr = F.problem(rng)
        if pr is None: continue
        cpdag, K, x, y, z, g0, dag = (pr["cpdag"], pr["K"], pr["x"], pr["y"], pr["z"],
                                      pr["g0"], pr["dag"])
        L = partition(cpdag, K, x, y, z, g0)[0]
        P = [k for k in K if k[0] in (x, y) or k[1] in (x, y)]
        S = sorted(set(L) | set(P)); Lo = sorted(set(L)); Po = sorted(set(P))
        sem = random_sem(dag, rng); cache = {}
        def th(g):
            k = g.edge_string()
            if k not in cache:
                cache[k] = [adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(d, x, y))
                            for d in enumerate_dag_extensions(g)]
            return cache[k]
        point = adjusted_estimand(sem, x, y, z); tau = sem.true_total_effect(x, y)
        def hedge(gs):
            v = [point]
            for g in gs: v += th(g)
            return min(v), max(v)
        sg = lambda C: [g for k in C for _, g in admissible(cpdag, K, k)]
        R.append(dict(tau=tau, point=point, nL=len(L), nP=len(P), nS=len(S), nK=len(K),
                      L2=hedge(sg(Lo) + pairs_revisions(cpdag, K, Lo)),
                      P2=hedge(sg(Po) + pairs_revisions(cpdag, K, Po)),
                      S2=hedge(sg(S)  + pairs_revisions(cpdag, K, S)),
                      K2=hedge(sg(K)  + pairs_revisions(cpdag, K, K)),
                      BL=hedge([cpdag]),
                      zbad=not is_valid_adjustment_set_dag(dag, x, y, z)))
    tau = np.array([r["tau"] for r in R]); pts = np.array([r["point"] for r in R])
    def cw(key):
        lo = np.minimum(np.array([r[key][0] for r in R]), pts)
        hi = np.maximum(np.array([r[key][1] for r in R]), pts)
        half = (hi-lo)/2; mid = (lo+hi)/2; need = np.abs(tau-mid)
        raw = float(np.mean(need <= half + 1e-12))
        c = np.where(half > 1e-12, need/np.where(half > 1e-12, half, 1.0),
                     np.where(need < 1e-9, 0.0, np.inf))
        fin = np.isfinite(c)
        if fin.mean() < 0.95: return raw, None, float("nan"), 1-fin.mean()
        q = float(np.quantile(np.sort(c[fin]), min(0.95/fin.mean(), 1.0)))
        return raw, q, float(np.mean(2*half*q)), 1-fin.mean()
    print(f"\n=== ablation, per-claim error {p_err}, n={len(R)} "
          f"(Z invalid in truth on {100*np.mean([r['zbad'] for r in R]):.1f}%) ===")
    print(f"    mean |L| {np.mean([r['nL'] for r in R]):.2f}  |P| {np.mean([r['nP'] for r in R]):.2f}  "
          f"|LuP| {np.mean([r['nS'] for r in R]):.2f}  |K| {np.mean([r['nK'] for r in R]):.2f}")
    for name, key in (("order-2 on L only", "L2"), ("order-2 on P only (proximity)", "P2"),
                      ("order-2 on L u P  (PROPOSED)", "S2"), ("order-2 on all of K", "K2"),
                      ("blanket over the class", "BL")):
        raw, q, w, fi = cw(key)
        print(f"    {name:32s} raw {raw:.4f}  q {'  n/a' if q is None else f'{q:5.3f}'}"
              f"  width {w:7.4f}  uncoverable {fi:.4f}")
    d = sum(1 for r in R if r["S2"] != r["P2"]); e = sum(1 for r in R if r["S2"] != r["L2"])
    print(f"    L u P differs from P alone on {100*d/len(R):.1f}% of problems; "
          f"from L alone on {100*e/len(R):.1f}%")

for pe, sd in ((0.35, 22), (0.50, 23)):
    run(pe, 250, sd)
