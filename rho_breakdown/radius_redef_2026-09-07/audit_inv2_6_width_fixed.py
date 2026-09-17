#!/usr/bin/env python3
"""inverse-2 audit, part 6: width at honest 95% coverage, with the estimand range
taken over DAG EXTENSIONS (part 5's O*-on-MPDAG version collapsed to a point)."""
import sys
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import random_dag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, optimal_adjustment_set_dag,
                                    adjusted_estimand, is_valid_adjustment_set_mpdag,
                                    is_valid_adjustment_set_dag, random_sem)
P_ERR = 0.20

def flags(cpdag, Kp, x, y, z, g0):
    lb, den = [], 0
    for k in Kp:
        rest = [e for e in Kp if e != k]
        gr = apply_orientations(cpdag, rest)
        gv = apply_orientations(cpdag, rest + [(k[1], k[0])])
        adm = [g for g in (gr, gv) if g is not None]
        if adm: den += 1
        if gv is None: continue
        if gr is not None and gr == g0: continue
        if any(not is_valid_adjustment_set_mpdag(g, x, y, z) for g in adm):
            lb.append(k)
    return lb, (len(lb)/den if den else 0.0)

def thetas(sem, g, x, y, cache):
    key = g.edge_string()
    if key in cache: return cache[key]
    out = []
    for d in enumerate_dag_extensions(g):
        try: out.append(adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(d, x, y)))
        except Exception: pass
    cache[key] = out
    return out

def problem(rng):
    dag = random_dag(rng, 6, 0.35)
    if not dag.directed_edges: return None
    cpdag = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cpdag.undirected_edges)
    if not (2 <= len(und) <= 5): return None
    nodes = sorted(dag.nodes)
    cand = [(a, b) for a in nodes for b in nodes if a != b]
    rng.shuffle(cand)
    for x, y in cand:
        idx = rng.permutation(len(und))[:min(4, len(und))]
        Kp, wrong = [], []
        for t in idx:
            a, b = und[t]
            tr = (a, b) if (a, b) in dag.directed_edges else (b, a)
            if rng.random() < P_ERR: Kp.append((tr[1], tr[0])); wrong.append((tr[1], tr[0]))
            else: Kp.append(tr)
        g0 = apply_orientations(cpdag, Kp)
        if g0 is None: continue
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None: continue
        return dict(dag=dag, cpdag=cpdag, x=x, y=y, K=Kp, wrong=set(wrong),
                    g0=g0, z=frozenset(z))
    return None

rng = np.random.default_rng(4242)
R = []; tried = 0
while len(R) < 250 and tried < 60000:
    tried += 1
    pr = problem(rng)
    if pr is None: continue
    cpdag, K, x, y, z, g0, dag = (pr["cpdag"], pr["K"], pr["x"], pr["y"], pr["z"],
                                  pr["g0"], pr["dag"])
    sem = random_sem(dag, rng); cache = {}
    tau = sem.true_total_effect(x, y)
    point = adjusted_estimand(sem, x, y, z)
    L, phi = flags(cpdag, K, x, y, z, g0)
    def hedge(claims):
        vals = list(thetas(sem, g0, x, y, cache)) + [point]
        for k in claims:
            rest = [e for e in K if e != k]
            for extra in ([], [(k[1], k[0])]):
                g = apply_orientations(cpdag, rest + extra)
                if g is not None: vals += thetas(sem, g, x, y, cache)
        return min(vals), max(vals)
    sel = hedge(L); alk = hedge(K)
    bl = thetas(sem, cpdag, x, y, cache) or [point]
    R.append(dict(tau=tau, point=point, nL=len(L), nK=len(K), nwrong=len(pr["wrong"]),
                  hit=len(set(L) & pr["wrong"]),
                  zbad=not is_valid_adjustment_set_dag(dag, x, y, z),
                  pt=(point, point), sel=sel, alk=alk, bla=(min(bl), max(bl))))

tau = np.array([r["tau"] for r in R])
def calib(key):
    lo = np.array([r[key][0] for r in R]); hi = np.array([r[key][1] for r in R])
    lo = np.minimum(lo, np.array([r["point"] for r in R]))
    hi = np.maximum(hi, np.array([r["point"] for r in R]))
    mid = (lo+hi)/2; half = (hi-lo)/2; need = np.abs(tau-mid)
    c = np.where(half > 0, need/np.where(half > 0, half, 1), np.where(need == 0, 0.0, np.inf))
    fin = np.isfinite(c)
    if fin.mean() < 0.95:
        return None, None, np.mean(need <= half), 1-fin.mean()
    q = float(np.quantile(np.sort(c[fin]), min(0.95/fin.mean(), 1.0)))
    return q, float(np.mean(2*half*q)), float(np.mean(need <= half)), 1-fin.mean()

print(f"problems={len(R)}  mean|K|={np.mean([r['nK'] for r in R]):.2f}  "
      f"mean wrong={np.mean([r['nwrong'] for r in R]):.3f}  "
      f"Z invalid in truth: {100*np.mean([r['zbad'] for r in R]):.1f}%")
print("\nmethod                                  raw95cov   scale     MEAN WIDTH   uncoverable")
for name, key in (("point estimate on the elicited set", "pt"),
                  ("selective hedge (phi_1 flagged claims)", "sel"),
                  ("hedge over EVERY claim, one revision", "alk"),
                  ("blanket hedge over the whole class", "bla")):
    q, w, raw, fi = calib(key)
    print(f"  {name:38s} {raw:6.3f}   "
          f"{'  n/a  ' if q is None else f'{q:7.2f}'}   "
          f"{'    n/a ' if w is None else f'{w:8.3f}'}   {fi:6.3f}")
mw = np.mean([r["sel"][1]-r["sel"][0] for r in R]); mb = np.mean([r["bla"][1]-r["bla"][0] for r in R])
ma = np.mean([r["alk"][1]-r["alk"][0] for r in R])
print(f"\nraw (uncalibrated) mean widths: selective {mw:.3f}  all-claims {ma:.3f}  blanket {mb:.3f}")
same = np.mean([r["sel"] == r["alk"] for r in R])
print(f"selective interval identical to the all-claims interval on {100*same:.1f}% of problems")
zero = np.mean([r["sel"][1]-r["sel"][0] == 0 for r in R])
print(f"selective interval has ZERO width on {100*zero:.1f}% of problems")
