#!/usr/bin/env python3
"""Audit 4: (a) is lambda==0 really 'no identified effect changes'?  (b) cost claim.
(c) gaming by over-claiming under a per-claim error RATE.  (d) retrieval quality."""
import sys, time
import numpy as np
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, optimal_adjustment_set_dag,
                                    is_valid_adjustment_set_dag, is_valid_adjustment_set_mpdag,
                                    adjusted_estimand, random_sem)
TOL = 1e-9
def random_dag(rng, n, p):
    o = list(rng.permutation(n)); e = set()
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < p: e.add((f"v{o[i]}", f"v{o[j]}"))
    return MPDAG(nodes=[f"v{i}" for i in range(n)], directed=e, undirected=[])
def theta_D(sem, d, x, y):
    o = optimal_adjustment_set_dag(d, x, y)
    if not is_valid_adjustment_set_dag(d, x, y, o):
        o = d.parents(x)
        if not is_valid_adjustment_set_dag(d, x, y, o): return None
    return adjusted_estimand(sem, x, y, o)
def theta_set(sem, g, x, y):
    return [t for t in (theta_D(sem, d, x, y) for d in enumerate_dag_extensions(g)) if t is not None]
def width(v): return (max(v)-min(v)) if v else 0.0
def fibre(cp, K, k):
    rest = [e for e in K if e != k]; out = []
    g = apply_orientations(cp, rest) if rest else cp
    if g is not None: out.append(g)
    g2 = apply_orientations(cp, rest + [(k[1], k[0])])
    if g2 is not None: out.append(g2)
    return out
def make(rng, n=9, p=0.30, k_claims=4, max_und=7):
    dag = random_dag(rng, n, p)
    if not dag.directed_edges: return None
    cp = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cp.undirected_edges)
    if not (2 <= len(und) <= max_und): return None
    sem = random_sem(dag, rng); nodes = sorted(dag.nodes)
    cand = [(x, y) for x in nodes for y in nodes if x != y and y in dag.descendants(x)]
    rng.shuffle(cand)
    for x, y in cand:
        thC = theta_set(sem, cp, x, y)
        if width(thC) <= TOL: continue
        idx = rng.permutation(len(und))[:min(k_claims, len(und))]
        K = [((a, b) if (a, b) in dag.directed_edges else (b, a)) for a, b in (und[t] for t in idx)]
        if not K: continue
        g0 = apply_orientations(cp, K)
        if g0 is None: continue
        if optimal_adjustment_set_mpdag(g0, x, y) is None: continue
        return dict(dag=dag, cpdag=cp, x=x, y=y, K=[tuple(e) for e in K], g0=g0,
                    z=frozenset(optimal_adjustment_set_mpdag(g0, x, y)), sem=sem,
                    thC=thC, und=und)
    return None

rng = np.random.default_rng(20260907)
P, tried = [], 0
N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
while len(P) < N and tried < 500000:
    tried += 1
    pr = make(rng)
    if pr is not None: P.append(pr)
print(f"# {len(P)} problems")

# ---- (a) what lambda==0 really certifies -----------------------------------
w0_pos = 0; lam0 = 0; lam0_moves = 0; lam0_zbreak = 0; lam0_and_w0pos = 0
for pr in P:
    sem, cp, K, x, y, g0 = pr["sem"], pr["cpdag"], pr["K"], pr["x"], pr["y"], pr["g0"]
    th0 = theta_set(sem, g0, x, y); w0 = width(th0)
    if w0 > TOL: w0_pos += 1
    lo0, hi0 = (min(th0), max(th0)) if th0 else (0, 0)
    s0 = set(round(v, 9) for v in th0)
    for k in K:
        fv, zb = [], False
        for g in fibre(cp, K, k):
            fv.extend(theta_set(sem, g, x, y))
            if not is_valid_adjustment_set_mpdag(g, x, y, pr["z"]): zb = True
        lam = width(th0 + fv) - w0
        if lam <= TOL:
            lam0 += 1
            if w0 > TOL: lam0_and_w0pos += 1
            if any(round(v, 9) not in s0 for v in fv): lam0_moves += 1
            if zb: lam0_zbreak += 1
print(f"\n== (a) lambda==0: what does it certify? ==")
print(f"problems with w(Theta(G0)) > 0 (analyst's own set NOT a point): {w0_pos}/{len(P)} "
      f"({100*w0_pos/len(P):.1f}%)")
print(f"lambda==0 claims: {lam0}")
print(f"   ... of which the fibre DOES contain an effect not in Theta(G0): {lam0_moves} "
      f"({100*lam0_moves/max(lam0,1):.1f}%)  <-- silent false negatives of the ranked list")
print(f"   ... of which the fibre invalidates Z: {lam0_zbreak} ({100*lam0_zbreak/max(lam0,1):.1f}%)")
print(f"   ... occurring on a problem where w(Theta(G0))>0: {lam0_and_w0pos}")

# ---- (b) cost -------------------------------------------------------------
cf, cb = [], []
for pr in P:
    cp, K = pr["cpdag"], pr["K"]
    c = sum(len(enumerate_dag_extensions(g)) for k in K for g in fibre(cp, K, k))
    cf.append(c); cb.append(len(enumerate_dag_extensions(cp)))
ex = sum(1 for a, b in zip(cf, cb) if a > b)
print(f"\n== (b) cost: 'strictly cheaper than the blanket-hedge baseline' ==")
print(f"mean sum|ext(fibres)| = {np.mean(cf):.1f}  vs  mean |ext(Chat)| = {np.mean(cb):.1f}")
print(f"fibre total EXCEEDS blanket on {ex}/{len(P)} problems ({100*ex/len(P):.1f}%); "
      f"mean ratio {np.mean([a/b for a,b in zip(cf,cb)]):.2f}x")

# ---- (c) gaming by over-claiming, per-claim error rate q -------------------
print(f"\n== (c) over-claiming: per-claim error rate q=0.20, delta=0 ==")
print(f"{'|K|':>4} {'cov':>6} {'width':>8} {'zero-w':>7} {'n':>5}")
for kc in (1, 2, 3, 4, 5, 6):
    rng2 = np.random.default_rng(101)
    cov = 0; ws = []; zw = 0; used = 0
    for pr in P:
        und, dag, cp, sem, x, y = pr["und"], pr["dag"], pr["cpdag"], pr["sem"], pr["x"], pr["y"]
        if len(und) < kc: continue
        idx = rng2.permutation(len(und))[:kc]
        Kt = [((a, b) if (a, b) in dag.directed_edges else (b, a)) for a, b in (und[t] for t in idx)]
        Kst = [((e[1], e[0]) if rng2.random() < 0.20 else e) for e in Kt]
        g0 = apply_orientations(cp, Kst)
        if g0 is None: continue
        if optimal_adjustment_set_mpdag(g0, x, y) is None: continue
        th0 = theta_set(sem, g0, x, y)
        if not th0: continue
        u = list(th0)
        for k in Kst:
            for g in fibre(cp, Kst, k): u.extend(theta_set(sem, g, x, y))
        lo, hi = min(u), max(u); truth = sem.true_total_effect(x, y)
        used += 1
        if lo - 1e-8 <= truth <= hi + 1e-8: cov += 1
        ws.append(hi - lo)
        if hi - lo <= TOL: zw += 1
    if used:
        print(f"{kc:>4} {cov/used:>6.3f} {np.mean(ws):>8.4f} {zw/used:>7.3f} {used:>5}")

# ---- (d) retrieval: does the lambda ranking find the WRONG claim? ----------
print(f"\n== (d) retrieval: rank of the single reversed claim under lambda, b=1 ==")
rng3 = np.random.default_rng(5)
ranks = []; ties_at_zero = 0; wrong_is_zero = 0; used = 0
for pr in P:
    K = pr["K"]; cp, sem, x, y = pr["cpdag"], pr["sem"], pr["x"], pr["y"]
    if len(K) < 2: continue
    j = int(rng3.integers(len(K)))
    Kst = [((e[1], e[0]) if i == j else e) for i, e in enumerate(K)]
    g0 = apply_orientations(cp, Kst)
    if g0 is None: continue
    if optimal_adjustment_set_mpdag(g0, x, y) is None: continue
    th0 = theta_set(sem, g0, x, y)
    if not th0: continue
    w0 = width(th0); lams = []
    for k in Kst:
        fv = []
        for g in fibre(cp, Kst, k): fv.extend(theta_set(sem, g, x, y))
        lams.append(width(th0 + fv) - w0)
    used += 1
    if lams[j] <= TOL: wrong_is_zero += 1
    if all(l <= TOL for l in lams): ties_at_zero += 1
    order = sorted(range(len(lams)), key=lambda i: -lams[i])
    ranks.append(order.index(j) + 1)
if used:
    print(f"cases {used}; the reversed claim scores lambda==0 (unrankable) on "
          f"{wrong_is_zero} ({100*wrong_is_zero/used:.1f}%)")
    print(f"whole ranking tied at zero on {ties_at_zero} ({100*ties_at_zero/used:.1f}%)")
    print(f"mean rank of the reversed claim: {np.mean(ranks):.2f} of {np.mean([len(pr['K']) for pr in P]):.2f}; "
          f"top-1 hit rate {sum(1 for r in ranks if r==1)/used:.3f}")
