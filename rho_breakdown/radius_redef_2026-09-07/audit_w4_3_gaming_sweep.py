#!/usr/bin/env python3
"""r_prod audit 3: (a) coverage when the misdeclared claim is LOAD-BEARING,
(b) the monotone gaming sweep -- width and coverage as a function of how much
of K the protocol declares unflippable.

(a) is the stratum the paper is actually about. Audit 2 draws the false claim
uniformly, so most reversals are benign and the hull survives by accident.
Conditioning on the reversal mattering (Z moves) is the honest read.

(b) r_prod is MONOTONE INCREASING in |K \\ F| and the hedge width is MONOTONE
DECREASING in it, with no data term opposing. That is the shape of a free
parameter, and this measures the slope.
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
    vals = []
    for d in enumerate_dag_extensions(cpdag):
        C = frozenset(c for c in Kp if not d.is_directed_edge(c[0], c[1]))
        if len(C) > t or not (C <= F):
            continue
        vals.append(adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(d, x, y)))
    return (min(vals), max(vals)) if vals else None

def run():
    rng = np.random.default_rng(20260907)
    got = tried = 0
    LB, NB, sweep = [], [], {}
    while got < 700 and tried < 90000:
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
        got += 1
        sem = random_sem(dag, rng); truth = sem.true_total_effect(x, y)
        FULL = frozenset(Kp); cbad = (cstar[1], cstar[0])
        # load-bearing = the reversal actually invalidates the analyst's Z in the truth
        load = not is_valid_adjustment_set_dag(dag, x, y, zp)
        h = hull(cpdag, Kp, FULL - {cbad}, 1, sem, x, y)
        if h:
            miss = max(0.0, max(h[0]-truth, truth-h[1]))
            rec = dict(w=h[1]-h[0], cov=miss <= 1e-9*max(1,abs(truth)), miss=miss)
            (LB if load else NB).append(rec)
        # gaming sweep: declare the first j claims of a fixed order unflippable
        order = sorted(FULL)
        for j in range(len(order)+1):
            F = FULL - set(order[:j])
            hh = hull(cpdag, Kp, F, 1, sem, x, y)
            if hh is None: continue
            miss = max(0.0, max(hh[0]-truth, truth-hh[1]))
            frac = j/len(order)
            b = min(int(frac*4), 3)
            sweep.setdefault(b, []).append((hh[1]-hh[0], miss <= 1e-9*max(1,abs(truth)), miss))
    return got, tried, LB, NB, sweep

if __name__ == "__main__":
    got, tried, LB, NB, sweep = run()
    print(f"problemas: {got} de {tried}\n")
    print("(a) o braco B (a afirmacao falsa foi declarada INFLIPAVEL), estratificado:")
    for name, r in (("LOAD-BEARING (Z invalido na verdade)", LB), ("benigna (Z segue valido)", NB)):
        if not r: continue
        cov = np.mean([d["cov"] for d in r]); w = np.mean([d["w"] for d in r])
        ms = [d["miss"] for d in r if not d["cov"]]
        print(f"   {name:38s} n={len(r):4d}  cobertura {cov:.3f}  largura {w:.4f}")
        if ms:
            print(f"   {'':38s}   erro mediano {np.median(ms):.4f}  p90 {np.quantile(ms,.9):.4f}  max {max(ms):.4f}")
            print(f"   {'':38s}   erro / largura mediano = {np.median(ms)/max(w,1e-12):.1f}x")
    print("\n(b) varredura de gaming: quanto de K o protocolo declara inflipavel ->")
    print("      faixa de |K\\F|/|K|      n     largura media   cobertura")
    lab = {0:"[0.00,0.25)",1:"[0.25,0.50)",2:"[0.50,0.75)",3:"[0.75,1.00]"}
    for b in sorted(sweep):
        v = sweep[b]
        print(f"      {lab[b]:16s} {len(v):5d}   {np.mean([a for a,_,_ in v]):10.4f}    {np.mean([c for _,c,_ in v]):.3f}")
