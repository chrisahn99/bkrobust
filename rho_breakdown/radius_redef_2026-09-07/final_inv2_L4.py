#!/usr/bin/env python3
"""FINAL inverse-2, part 4: the fallback variant L-else-P, two seeds, n=400."""
import sys
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
    while len(R) < n_target and tried < 200000:
        tried += 1
        pr = F.problem(rng)
        if pr is None: continue
        cpdag, K, x, y, z, g0, dag = (pr["cpdag"], pr["K"], pr["x"], pr["y"], pr["z"],
                                      pr["g0"], pr["dag"])
        L = partition(cpdag, K, x, y, z, g0)[0]
        P = [k for k in K if k[0] in (x, y) or k[1] in (x, y)]
        sem = random_sem(dag, rng); cache = {}
        def th(g):
            k = g.edge_string()
            if k not in cache:
                cache[k] = [adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(d, x, y))
                            for d in enumerate_dag_extensions(g)]
            return cache[k]
        point = adjusted_estimand(sem, x, y, z); tau = sem.true_total_effect(x, y)
        def H(C):
            C = sorted(set(C)); v = [point]
            for k in C:
                for _, g in admissible(cpdag, K, k): v += th(g)
            for g in pairs_revisions(cpdag, K, C): v += th(g)
            return min(v), max(v)
        R.append(dict(tau=tau, point=point, nL=len(L), nP=len(P), nK=len(K),
                      A=H(L), B=H(L if L else P), C=H(set(L)|set(P)), D=H(P), E=H(K),
                      Fb=(min([point]+th(cpdag)), max([point]+th(cpdag))),
                      zbad=not is_valid_adjustment_set_dag(dag, x, y, z)))
    tau = np.array([r["tau"] for r in R]); pts = np.array([r["point"] for r in R])
    def cw(key):
        lo = np.minimum(np.array([r[key][0] for r in R]), pts)
        hi = np.maximum(np.array([r[key][1] for r in R]), pts)
        half = (hi-lo)/2; mid = (lo+hi)/2; need = np.abs(tau-mid)
        raw = float(np.mean(need <= half+1e-12))
        c = np.where(half > 1e-12, need/np.where(half > 1e-12, half, 1.0),
                     np.where(need < 1e-9, 0.0, np.inf))
        fin = np.isfinite(c)
        if fin.mean() < 0.95: return raw, None, float("nan"), 1-fin.mean()
        q = float(np.quantile(np.sort(c[fin]), min(0.95/fin.mean(), 1.0)))
        return raw, q, float(np.mean(2*half*q)), 1-fin.mean()
    print(f"\n=== p_err {p_err}, seed {seed}, n={len(R)}  (Z invalid in truth "
          f"{100*np.mean([r['zbad'] for r in R]):.1f}%; mean |L| {np.mean([r['nL'] for r in R]):.2f}, "
          f"|P| {np.mean([r['nP'] for r in R]):.2f}, |K| {np.mean([r['nK'] for r in R]):.2f}) ===")
    for name, key in (("order<=2 on L", "A"), ("order<=2 on L, else P if L empty", "B"),
                      ("order<=2 on L u P", "C"), ("order<=2 on P", "D"),
                      ("order<=2 on all of K", "E"), ("blanket over class", "Fb")):
        raw, q, w, fi = cw(key)
        print(f"    {name:34s} raw {raw:.4f}  q {'  n/a' if q is None else f'{q:5.3f}'}"
              f"  width {w:7.4f}  uncoverable {fi:.4f}")

for pe, sd in ((0.50, 31), (0.50, 32), (0.35, 33)):
    run(pe, 400, sd)
