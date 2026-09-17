#!/usr/bin/env python3
"""r_prod audit 4: the load-bearing stratum at n large enough to quote.
Only problems where reversing the claim genuinely invalidates Z in the truth.
Compares F = K (claim flippable) against F = K \\ {false claim} (declared unflippable).
"""
import sys
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import make_problem
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, optimal_adjustment_set_dag,
                                    is_valid_adjustment_set_dag, adjusted_estimand, random_sem)

def hull(cpdag, Kp, F, t, sem, x, y):
    v = []
    for d in enumerate_dag_extensions(cpdag):
        C = frozenset(c for c in Kp if not d.is_directed_edge(c[0], c[1]))
        if len(C) <= t and C <= F:
            v.append(adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(d, x, y)))
    return (min(v), max(v)) if v else None

rng = np.random.default_rng(11071988)
A, B, tried, lb = [], [], 0, 0
while lb < 300 and tried < 400000:
    tried += 1
    pr = make_problem(rng)
    if pr is None: continue
    K, cpdag, dag, x, y = [tuple(e) for e in pr["K"]], pr["cpdag"], pr["dag"], pr["x"], pr["y"]
    if not K or not is_valid_adjustment_set_dag(dag, x, y, pr["z"]): continue
    cstar = K[int(rng.integers(len(K)))]
    Kp = [((c[1], c[0]) if c == cstar else c) for c in K]
    g0 = apply_orientations(cpdag, Kp)
    if g0 is None: continue
    zp = optimal_adjustment_set_mpdag(g0, x, y)
    if zp is None: continue
    if is_valid_adjustment_set_dag(dag, x, y, zp):   # not load-bearing
        continue
    lb += 1
    sem = random_sem(dag, rng); truth = sem.true_total_effect(x, y)
    FULL = frozenset(Kp); cbad = (cstar[1], cstar[0])
    for arm, F, out in (("flippable", FULL, A), ("unflippable", FULL - {cbad}, B)):
        h = hull(cpdag, Kp, F, 1, sem, x, y)
        if h is None: continue
        miss = max(0.0, max(h[0]-truth, truth-h[1]))
        out.append((h[1]-h[0], miss <= 1e-9*max(1, abs(truth)), miss))

print(f"estrato LOAD-BEARING: {lb} problemas de {tried} sorteios\n")
for name, r in (("F = K            (flipavel)", A), ("F = K\\{c*}  (INFLIPAVEL)", B)):
    w = np.mean([a for a,_,_ in r]); cov = np.mean([c for _,c,_ in r])
    ms = [m for _,c,m in r if not c]
    se = np.sqrt(cov*(1-cov)/len(r))
    print(f"{name}  n={len(r)}  cobertura {cov:.3f} +- {se:.3f}   largura media {w:.4f}")
    if ms:
        print(f"{'':30s} erro quando falha: mediana {np.median(ms):.4f}  p90 {np.quantile(ms,.9):.4f}  max {max(ms):.4f}")
        print(f"{'':30s} razao erro/largura: mediana {np.median(ms)/max(w,1e-12):.1f}x")
