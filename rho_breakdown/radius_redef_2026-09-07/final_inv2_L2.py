#!/usr/bin/env python3
"""FINAL inverse-2, part 2: width table across error regimes + padding direction."""
import sys, itertools
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_dag, adjusted_estimand,
                                    is_valid_adjustment_set_dag, random_sem)
import final_inv2_L as F
from final_inv2_L import partition, admissible, pairs_revisions, rev

def run(p_err, n_target, seed):
    F.P_ERR = p_err
    rng = np.random.default_rng(seed)
    R, tried = [], 0
    up = dn = same = 0
    while len(R) < n_target and tried < 80000:
        tried += 1
        pr = F.problem(rng)
        if pr is None: continue
        cpdag, K, x, y, z, g0, dag = (pr["cpdag"], pr["K"], pr["x"], pr["y"], pr["z"],
                                      pr["g0"], pr["dag"])
        L, Fr, I, D = partition(cpdag, K, x, y, z, g0)
        P = [k for k in K if k[0] in (x, y) or k[1] in (x, y)]
        S = sorted(set(L) | set(P))
        sem = random_sem(dag, rng); cache = {}
        def thetas(g):
            key = g.edge_string()
            if key not in cache:
                cache[key] = [adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(d, x, y))
                              for d in enumerate_dag_extensions(g)]
            return cache[key]
        point = adjusted_estimand(sem, x, y, z); tau = sem.true_total_effect(x, y)
        def hedge(gs):
            v = [point]
            for g in gs: v += thetas(g)
            return min(v), max(v)
        singles = lambda C: [g for k in C for _, g in admissible(cpdag, K, k)]
        R.append(dict(tau=tau, point=point, nL=len(L), nK=len(K),
                      hL=hedge(singles(L)), hK=hedge(singles(K)),
                      hS2=hedge(singles(S) + pairs_revisions(cpdag, K, S)),
                      hK2=hedge(singles(K) + pairs_revisions(cpdag, K, K)),
                      hbl=hedge([cpdag]),
                      zbad=not is_valid_adjustment_set_dag(dag, x, y, z)))
        ent = [tuple(e) for e in sorted(g0.directed_edges) if tuple(e) not in set(K)]
        if ent and apply_orientations(cpdag, list(K)+ent) == g0:
            L2 = partition(cpdag, list(K)+ent, x, y, z, g0)[0]
            up += len(L2) > len(L); dn += len(L2) < len(L); same += len(L2) == len(L)
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
    print(f"\n=== per-claim error rate {p_err}  (n={len(R)}) ===")
    print(f"   Z invalid in the true DAG on {100*np.mean([r['zbad'] for r in R]):.1f}% of problems")
    print(f"   point estimate raw coverage {np.mean(np.abs(tau-pts) < 1e-9):.3f}")
    cL = Counter(r["nL"] for r in R)
    print(f"   |L| profile {dict(sorted(cL.items()))}  max cell {100*max(cL.values())/len(R):.1f}%")
    for name, key in (("L only, order 1", "hL"), ("all claims, order 1", "hK"),
                      ("PROPOSED  L u P, order 1 + order 2", "hS2"),
                      ("all claims, order 1 + order 2", "hK2"),
                      ("blanket over the class", "hbl")):
        raw, q, w, fi = cw(key)
        print(f"   {name:38s} raw {raw:.4f}  q {'  n/a' if q is None else f'{q:5.3f}'}"
              f"  width {w:7.4f}  uncoverable {fi:.4f}")
    print(f"   padding attack on |L|: up {up}, down {dn}, unchanged {same}")

for pe, sd in ((0.20, 11), (0.35, 12), (0.50, 13)):
    run(pe, 250, sd)
