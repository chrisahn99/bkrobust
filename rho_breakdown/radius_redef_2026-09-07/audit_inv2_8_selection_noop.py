#!/usr/bin/env python3
"""inverse-2 audit, part 8: does the SELECTION step do any work?
Hedge support = {point} u {revisions}, WITHOUT the extensions-of-G0 envelope, so the
only difference between the two hedges is which claims get revised."""
import sys
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, optimal_adjustment_set_dag,
                                    adjusted_estimand, is_valid_adjustment_set_mpdag,
                                    is_valid_adjustment_set_dag, random_sem)
from harness import random_dag
P_ERR = 0.20

def parts(cpdag, K, x, y, z, g0):
    lb, free, inert, dref, den = [], [], [], [], 0
    for k in K:
        rest = [e for e in K if e != k]
        gr = apply_orientations(cpdag, rest)
        gv = apply_orientations(cpdag, rest + [(k[1], k[0])])
        adm = [g for g in (gr, gv) if g is not None]
        if adm: den += 1
        if gv is None: dref.append(k)
        elif gr is not None and gr == g0: inert.append(k)
        elif any(not is_valid_adjustment_set_mpdag(g, x, y, z) for g in adm): lb.append(k)
        else: free.append(k)
    return lb, free, inert, dref, (len(lb)/den if den else 0.0)

def thetas(sem, g, x, y, cache):
    key = g.edge_string()
    if key not in cache:
        cache[key] = [adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(d, x, y))
                      for d in enumerate_dag_extensions(g)]
    return cache[key]

def problem(rng):
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
            if rng.random() < P_ERR: K.append((tr[1], tr[0])); wrong.append((tr[1], tr[0]))
            else: K.append(tr)
        g0 = apply_orientations(cpdag, K)
        if g0 is None: continue
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None: continue
        return dict(dag=dag, cpdag=cpdag, x=x, y=y, K=K, wrong=set(wrong), g0=g0,
                    z=frozenset(z))
    return None

rng = np.random.default_rng(4242)
R = []; tried = 0; pad_before = []; pad_after = []
while len(R) < 300 and tried < 60000:
    tried += 1
    pr = problem(rng)
    if pr is None: continue
    cpdag, K, x, y, z, g0, dag = (pr["cpdag"], pr["K"], pr["x"], pr["y"], pr["z"],
                                  pr["g0"], pr["dag"])
    sem = random_sem(dag, rng); cache = {}
    tau = sem.true_total_effect(x, y); point = adjusted_estimand(sem, x, y, z)
    L, F, I, D, phi = parts(cpdag, K, x, y, z, g0)
    def hedge(claims):
        vals = [point]
        for k in claims:
            rest = [e for e in K if e != k]
            for extra in ([], [(k[1], k[0])]):
                g = apply_orientations(cpdag, rest + extra)
                if g is not None: vals += thetas(sem, g, x, y, cache)
        return min(vals), max(vals)
    sel, alk = hedge(L), hedge(K)
    bl = thetas(sem, cpdag, x, y, cache) or [point]
    R.append(dict(tau=tau, point=point, sel=sel, alk=alk, bla=(min(bl), max(bl)),
                  nL=len(L), nF=len(F), nI=len(I), nD=len(D), phi=phi,
                  zbad=not is_valid_adjustment_set_dag(dag, x, y, z),
                  hit=len(set(L) & pr["wrong"]), nwrong=len(pr["wrong"]), nK=len(K)))
    # padding attack, amenable setting
    ent = [tuple(e) for e in sorted(g0.directed_edges) if tuple(e) not in set(K)]
    if ent:
        Kp = list(K) + ent
        if apply_orientations(cpdag, Kp) == g0:
            pad_before.append(phi); pad_after.append(parts(cpdag, Kp, x, y, z, g0)[4])

tau = np.array([r["tau"] for r in R]); pts = np.array([r["point"] for r in R])
def cov_width(key):
    lo = np.minimum(np.array([r[key][0] for r in R]), pts)
    hi = np.maximum(np.array([r[key][1] for r in R]), pts)
    half = (hi-lo)/2; mid = (lo+hi)/2; need = np.abs(tau-mid)
    raw = float(np.mean(need <= half + 1e-12))
    c = np.where(half > 1e-12, need/np.where(half > 1e-12, half, 1.0),
                 np.where(need < 1e-9, 0.0, np.inf))
    fin = np.isfinite(c)
    if fin.mean() < 0.95: return raw, None, float(np.mean(2*half)), 1-fin.mean()
    q = float(np.quantile(np.sort(c[fin]), min(0.95/fin.mean(), 1.0)))
    return raw, q, float(np.mean(2*half*q)), 1-fin.mean()

print(f"n={len(R)}   partition sizes: load-bearing {sum(r['nL'] for r in R)}, "
      f"free {sum(r['nF'] for r in R)}, inert {sum(r['nI'] for r in R)}, "
      f"data-refuted {sum(r['nD'] for r in R)}  (of {sum(r['nK'] for r in R)} claims)")
print("\nDOES THE SELECTION STEP CHANGE THE INTERVAL?")
diff = [i for i, r in enumerate(R) if r["sel"] != r["alk"]]
print(f"   selective interval differs from the all-claims interval on "
      f"{len(diff)}/{len(R)} problems ({100*len(diff)/len(R):.1f}%)")
ws = np.array([r["sel"][1]-r["sel"][0] for r in R])
wa = np.array([r["alk"][1]-r["alk"][0] for r in R])
print(f"   mean raw width  selective {ws.mean():.4f}   all-claims {wa.mean():.4f}   "
      f"saving {100*(1-ws.mean()/max(wa.mean(),1e-12)):.2f}%")
print("\nWIDTH AT HONEST 95% COVERAGE")
for name, key in (("point on elicited set", "pt" if False else None),):
    pass
for name, key in (("selective hedge (phi_1 flags)", "sel"),
                  ("hedge over EVERY claim (no phi_1 needed)", "alk"),
                  ("blanket hedge over whole class", "bla")):
    raw, q, w, fi = cov_width(key)
    print(f"   {name:42s} raw {raw:.3f}  scale "
          f"{'n/a' if q is None else f'{q:5.2f}'}  width {w:7.4f}  uncoverable {fi:.3f}")
print("\nPADDING ATTACK (amenable setting)")
if pad_before:
    b = np.array(pad_before); a = np.array(pad_after)
    print(f"   n={len(b)}  mean phi_1 {b.mean():.4f} -> {a.mean():.4f}   "
          f"changed on {100*np.mean(b!=a):.1f}% of problems")
    nz = b > 0
    if nz.any():
        print(f"   among phi_1>0: {b[nz].mean():.4f} -> {a[nz].mean():.4f}  "
              f"(mean ratio {np.mean(a[nz]/b[nz]):.3f}, min ratio {np.min(a[nz]/b[nz]):.3f})")
