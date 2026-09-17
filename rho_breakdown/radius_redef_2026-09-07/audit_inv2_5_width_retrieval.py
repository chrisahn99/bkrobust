#!/usr/bin/env python3
"""inverse-2 audit, part 5: with an ERRING expert, does the phi_1 flag list
(a) retrieve the wrong claims and (b) buy width at honest 95% coverage?"""
import sys
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import random_dag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.space import enumerate_space
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, adjusted_estimand,
                                    is_valid_adjustment_set_mpdag,
                                    is_valid_adjustment_set_dag, random_sem)

P_ERR = 0.20


def flags(cpdag, Kp, x, y, z, g0):
    lb, den = [], 0
    for k in Kp:
        rest = [e for e in Kp if e != k]
        g_ret = apply_orientations(cpdag, rest)
        g_rev = apply_orientations(cpdag, rest + [(k[1], k[0])])
        adm = [g for g in (g_ret, g_rev) if g is not None]
        if adm: den += 1
        if g_rev is None: continue
        if g_ret is not None and g_ret == g0: continue
        if any(not is_valid_adjustment_set_mpdag(g, x, y, z) for g in adm):
            lb.append(k)
    return lb, (len(lb)/den if den else 0.0), den


def theta(sem, g, x, y):
    zz = optimal_adjustment_set_mpdag(g, x, y)
    return None if zz is None else adjusted_estimand(sem, x, y, zz)


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
        Ktrue, Kp, wrong = [], [], []
        for t in idx:
            a, b = und[t]
            tr = (a, b) if (a, b) in dag.directed_edges else (b, a)
            Ktrue.append(tr)
            if rng.random() < P_ERR:
                Kp.append((tr[1], tr[0])); wrong.append((tr[1], tr[0]))
            else:
                Kp.append(tr)
        g0 = apply_orientations(cpdag, Kp)
        if g0 is None: continue
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None: continue
        return dict(dag=dag, cpdag=cpdag, x=x, y=y, K=Kp, wrong=set(wrong),
                    g0=g0, z=frozenset(z))
    return None


rng = np.random.default_rng(4242)
rows = []
n = tried = 0
while n < 250 and tried < 60000:
    tried += 1
    pr = problem(rng)
    if pr is None: continue
    cpdag, K, x, y, z, g0, dag = (pr["cpdag"], pr["K"], pr["x"], pr["y"], pr["z"],
                                  pr["g0"], pr["dag"])
    sem = random_sem(dag, rng)
    tau = sem.true_total_effect(x, y)
    point = adjusted_estimand(sem, x, y, z)
    L, phi, den = flags(cpdag, K, x, y, z, g0)

    def hedge(claims):
        vals = [point]
        for k in claims:
            rest = [e for e in K if e != k]
            for extra in ([], [(k[1], k[0])]):
                g = apply_orientations(cpdag, rest + extra)
                if g is None: continue
                t = theta(sem, g, x, y)
                if t is not None: vals.append(t)
        return min(vals), max(vals)

    sel_lo, sel_hi = hedge(L)                    # selective: flagged claims only
    all_lo, all_hi = hedge(K)                    # every claim, one revision
    space = enumerate_space(cpdag)
    bl = [t for t in (theta(sem, g, x, y) for g in space) if t is not None]
    bla_lo, bla_hi = (min(bl), max(bl)) if bl else (point, point)

    rows.append(dict(tau=tau, point=point, phi=phi,
                     nL=len(L), nK=len(K), nwrong=len(pr["wrong"]),
                     hit=len(set(L) & pr["wrong"]),
                     zbad=not is_valid_adjustment_set_dag(dag, x, y, z),
                     sel=(sel_lo, sel_hi), alk=(all_lo, all_hi), bla=(bla_lo, bla_hi)))
    n += 1

R = rows
tau = np.array([r["tau"] for r in R]); pt = np.array([r["point"] for r in R])

def calibrate(key):
    lo = np.array([r[key][0] for r in R]); hi = np.array([r[key][1] for r in R])
    mid = (lo + hi) / 2; half = (hi - lo) / 2
    need = np.abs(tau - mid)
    with np.errstate(divide="ignore", invalid="ignore"):
        c = np.where(half > 0, need / half, np.inf)
    c = np.where((half == 0) & (need == 0), 0.0, c)
    finite = np.sort(c[np.isfinite(c)])
    frac_inf = np.mean(~np.isfinite(c))
    if frac_inf > 0.05:
        return None, np.mean(2*half), np.mean(need <= half), frac_inf
    q = np.quantile(finite, min(0.95 / (1 - frac_inf), 1.0))
    return q, np.mean(2 * half * q), np.mean(need <= half), frac_inf

print(f"problems = {len(R)}   mean |K| = {np.mean([r['nK'] for r in R]):.2f}   "
      f"mean # genuinely wrong = {np.mean([r['nwrong'] for r in R]):.3f}")
print(f"Z actually invalid in the true DAG: {100*np.mean([r['zbad'] for r in R]):.1f}%")

nL = np.array([r["nL"] for r in R]); nw = np.array([r["nwrong"] for r in R])
hit = np.array([r["hit"] for r in R])
prec = hit.sum() / max(nL.sum(), 1); rec = hit.sum() / max(nw.sum(), 1)
print(f"\nRETRIEVAL of the genuinely wrong claims by the load-bearing flag")
print(f"   flagged total = {nL.sum()}   wrong total = {nw.sum()}   hits = {hit.sum()}")
print(f"   micro precision = {prec:.3f}    micro recall = {rec:.3f}")
base = nw.sum() / sum(r['nK'] for r in R)
print(f"   base rate (a random claim is wrong) = {base:.3f}   lift = {prec/base:.2f}x")

# does the flag detect DAMAGE (Z invalid), which is what the paper needs?
zbad = np.array([r["zbad"] for r in R])
print(f"\nDAMAGE DETECTION  (phi_1 > 0  as a test for 'Z is invalid in the truth')")
pos = nL > 0
tp = int((pos & zbad).sum()); fp = int((pos & ~zbad).sum())
fn = int((~pos & zbad).sum()); tn = int((~pos & ~zbad).sum())
print(f"   TP {tp}  FP {fp}  FN {fn}  TN {tn}")
print(f"   precision = {tp/max(tp+fp,1):.3f}   recall = {tp/max(tp+fn,1):.3f}")
print(f"   SILENT FAILURES (phi_1 == 0 yet Z invalid) = {fn}  "
      f"({100*fn/max((~pos).sum(),1):.1f}% of the 'clean' verdicts)")

print(f"\nWIDTH REQUIRED FOR HONEST 95% COVERAGE  (scale each interval to 95%)")
for name, key in (("selective hedge (phi_1 flags)", "sel"),
                  ("hedge over EVERY claim, 1 revision", "alk"),
                  ("blanket hedge over the whole class", "bla")):
    q, w, raw, fi = calibrate(key)
    qs = "UNCALIBRATABLE" if q is None else f"{q:8.2f}"
    print(f"   {name:36s} raw cov {raw:.3f}   scale {qs}   "
          f"width {w:8.3f}   zero-width-miss {fi:.3f}")
