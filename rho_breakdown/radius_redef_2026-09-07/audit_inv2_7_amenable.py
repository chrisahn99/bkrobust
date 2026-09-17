#!/usr/bin/env python3
"""inverse-2 audit, part 7: restrict to queries where y IS a descendant of x (O* is
only the HPM optimal set there; for y not in de(x) it returns the empty set, which the
back-door check then calls invalid for reasons that have nothing to do with K).
Re-measures r_val degeneracy, phi_1, retrieval, and the width table."""
import sys
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from collections import Counter
from harness import random_dag, ball, r_val
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, optimal_adjustment_set_dag,
                                    adjusted_estimand, is_valid_adjustment_set_mpdag,
                                    is_valid_adjustment_set_dag, random_sem)
P_ERR = 0.20

def flags(cpdag, K, x, y, z, g0):
    lb, den = [], 0
    for k in K:
        rest = [e for e in K if e != k]
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
    if key not in cache:
        cache[key] = [adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(d, x, y))
                      for d in enumerate_dag_extensions(g)]
    return cache[key]

def problem(rng, err):
    dag = random_dag(rng, 6, 0.35)
    if not dag.directed_edges: return None
    cpdag = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cpdag.undirected_edges)
    if not (2 <= len(und) <= 5): return None
    nodes = sorted(dag.nodes)
    cand = [(a, b) for a in nodes for b in nodes if a != b and b in dag.descendants(a)]
    rng.shuffle(cand)
    for x, y in cand:
        idx = rng.permutation(len(und))[:min(4, len(und))]
        K, wrong = [], []
        for t in idx:
            a, b = und[t]
            tr = (a, b) if (a, b) in dag.directed_edges else (b, a)
            if err and rng.random() < P_ERR:
                K.append((tr[1], tr[0])); wrong.append((tr[1], tr[0]))
            else: K.append(tr)
        g0 = apply_orientations(cpdag, K)
        if g0 is None: continue
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None: continue
        return dict(dag=dag, cpdag=cpdag, x=x, y=y, K=K, wrong=set(wrong),
                    g0=g0, z=frozenset(z), undirected=und)
    return None

# ---- A. truthful expert: r_val degeneracy and phi_1, on AMENABLE queries only
rng = np.random.default_rng(20260907)
rv, ph = [], []
n = tried = 0
while n < 250 and tried < 60000:
    tried += 1
    pr = problem(rng, err=False)
    if pr is None: continue
    rows = ball(pr)
    if not rows: continue
    n += 1
    rv.append(r_val(rows))
    ph.append(flags(pr["cpdag"], pr["K"], pr["x"], pr["y"], pr["z"], pr["g0"])[1])
c = Counter("UNREACHED" if r is None else r for r in rv)
print(f"A. TRUTHFUL EXPERT, y in de(x) only  (n={n})")
for k in sorted(c, key=lambda v: (v == "UNREACHED", v)):
    print(f"   r_val={str(k):<10} {c[k]:3d}  {100*c[k]/n:5.1f}%")
cp = Counter(ph)
print(f"   phi_1: {len(cp)} distinct, mass at 0 = {100*cp.get(0.0,0)/n:.1f}%, "
      f"mass at 1 = {100*cp.get(1.0,0)/n:.1f}%, mean = {np.mean(ph):.3f}")
viol = sum(1 for p, r in zip(ph, rv) if (p > 0) != (r == 1))
print(f"   violations of 'phi_1>0 <=> r_val==1': {viol}/{n}")

# ---- B. erring expert: retrieval + width
rng = np.random.default_rng(4242)
R = []; tried = 0
while len(R) < 250 and tried < 60000:
    tried += 1
    pr = problem(rng, err=True)
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
    bl = thetas(sem, cpdag, x, y, cache) or [point]
    R.append(dict(tau=tau, point=point, nL=len(L), nK=len(K), nwrong=len(pr["wrong"]),
                  hit=len(set(L) & pr["wrong"]),
                  zbad=not is_valid_adjustment_set_dag(dag, x, y, z),
                  pt=(point, point), sel=hedge(L), alk=hedge(K), bla=(min(bl), max(bl))))

tau = np.array([r["tau"] for r in R]); pts = np.array([r["point"] for r in R])
def calib(key):
    lo = np.minimum(np.array([r[key][0] for r in R]), pts)
    hi = np.maximum(np.array([r[key][1] for r in R]), pts)
    mid = (lo+hi)/2; half = (hi-lo)/2; need = np.abs(tau-mid)
    c = np.where(half > 1e-12, need/np.where(half > 1e-12, half, 1.0),
                 np.where(need < 1e-9, 0.0, np.inf))
    fin = np.isfinite(c)
    raw = float(np.mean(need <= half + 1e-12))
    if fin.mean() < 0.95: return None, None, raw, 1-fin.mean()
    q = float(np.quantile(np.sort(c[fin]), min(0.95/fin.mean(), 1.0)))
    return q, float(np.mean(2*half*q)), raw, 1-fin.mean()

nL = np.array([r["nL"] for r in R]); nw = np.array([r["nwrong"] for r in R])
hit = np.array([r["hit"] for r in R]); zbad = np.array([r["zbad"] for r in R])
print(f"\nB. ERRING EXPERT (p={P_ERR}), y in de(x)  (n={len(R)})  "
      f"mean wrong={nw.mean():.3f}  Z invalid in truth {100*zbad.mean():.1f}%")
print(f"   retrieval of wrong claims: precision {hit.sum()/max(nL.sum(),1):.3f}  "
      f"recall {hit.sum()/max(nw.sum(),1):.3f}  "
      f"base {nw.sum()/sum(r['nK'] for r in R):.3f}")
pos = nL > 0
print(f"   damage detection: TP {int((pos&zbad).sum())} FP {int((pos&~zbad).sum())} "
      f"FN {int((~pos&zbad).sum())} TN {int((~pos&~zbad).sum())}  "
      f"silent-failure rate among clean verdicts "
      f"{100*int((~pos&zbad).sum())/max(int((~pos).sum()),1):.1f}%")
print("\n   method                                  raw95cov   scale    MEAN WIDTH  uncoverable")
for name, key in (("point estimate on elicited set", "pt"),
                  ("selective hedge (phi_1 flags)", "sel"),
                  ("hedge over EVERY claim", "alk"),
                  ("blanket hedge over whole class", "bla")):
    q, w, raw, fi = calib(key)
    print(f"   {name:38s} {raw:6.3f}  {'  n/a ' if q is None else f'{q:6.2f}'}  "
          f"{'     n/a' if w is None else f'{w:8.3f}'}   {fi:6.3f}")
print(f"\n   selective == all-claims interval on "
      f"{100*np.mean([r['sel']==r['alk'] for r in R]):.1f}% of problems; "
      f"selective is a POINT on {100*np.mean([r['sel'][1]-r['sel'][0] < 1e-12 for r in R]):.1f}%")
print(f"   raw mean widths: sel {np.mean([r['sel'][1]-r['sel'][0] for r in R]):.3f}  "
      f"all {np.mean([r['alk'][1]-r['alk'][0] for r in R]):.3f}  "
      f"blanket {np.mean([r['bla'][1]-r['bla'][0] for r in R]):.3f}")
