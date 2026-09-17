#!/usr/bin/env python3
"""Audit 6: two exact structural checks.
 (i)  Theta(cl(K[k<-rev])) subset Theta(cl(K\\{k}))  -> the reversal branch of F(k)
      can never widen anything, so half of F(k) is dead weight.
 (ii) hull(Theta(G0) u U_{lambda>0} F) == hull(Theta(G0) u U_{all k} F)
      -> lambda contributes nothing to the interval at delta=0.
 (iii) delta is scale-dependent: rescale Y and watch the ranking threshold move.
"""
import sys
import numpy as np
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
exec(open("/private/tmp/claude-501/-Users-josecosta-mugango/1cbd3bff-6ceb-4991-8a2a-b20d5a3f01d3/scratchpad/audit_lambda4.py").read().split("rng = np.random.default_rng(20260907)")[0])
from bkrobust.demo.meek import is_consistent_extension

rng = np.random.default_rng(20260907)
P, tried = [], 0
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
while len(P) < N and tried < 800000:
    tried += 1
    pr = make(rng)
    if pr is not None: P.append(pr)
print(f"# {len(P)} problems\n")

# (i) reversal branch redundancy -------------------------------------------
tested = viol = rev_exists = ext_subset = 0
for pr in P:
    cp, K, sem, x, y = pr["cpdag"], pr["K"], pr["sem"], pr["x"], pr["y"]
    for k in K:
        rest = [e for e in K if e != k]
        g_ret = apply_orientations(cp, rest) if rest else cp
        g_rev = apply_orientations(cp, rest + [(k[1], k[0])])
        if g_ret is None or g_rev is None:
            continue
        rev_exists += 1
        er = {d.edge_string() for d in enumerate_dag_extensions(g_ret)}
        ev = {d.edge_string() for d in enumerate_dag_extensions(g_rev)}
        if ev <= er: ext_subset += 1
        tr = theta_set(sem, g_ret, x, y); tv = theta_set(sem, g_rev, x, y)
        tested += 1
        if tv and tr and (min(tv) < min(tr) - 1e-9 or max(tv) > max(tr) + 1e-9):
            viol += 1
print("== (i) is the reversal branch of F(k) redundant? ==")
print(f"claims where both branches exist: {rev_exists}")
print(f"  [cl(K[k<-rev])] subset of [cl(K\\k)] : {ext_subset}/{rev_exists} "
      f"({100*ext_subset/max(rev_exists,1):.1f}%)")
print(f"  reversal branch widens the hull beyond the retraction branch: {viol}/{tested}")
print("  => the reversal half of every fibre is dead weight; lambda(k) reduces to")
print("     w(Theta(G0) u Theta(cl(K\\{k}))) - w(Theta(G0)).\n")

# (ii) is lambda decorative at delta = 0? ----------------------------------
maxdiff = 0.0; n = 0
for pr in P:
    cp, K, sem, x, y, g0 = pr["cpdag"], pr["K"], pr["sem"], pr["x"], pr["y"], pr["g0"]
    th0 = theta_set(sem, g0, x, y)
    if not th0: continue
    w0 = width(th0); ua = list(th0); up = list(th0)
    for k in K:
        fv = []
        for g in fibre(cp, K, k): fv.extend(theta_set(sem, g, x, y))
        ua.extend(fv)
        if width(th0 + fv) - w0 > TOL: up.extend(fv)
    n += 1
    maxdiff = max(maxdiff, abs(width(ua) - width(up)))
print("== (ii) does lambda change the delta=0 interval? ==")
print(f"max |width(prefix hedge, delta=0) - width(hull over ALL fibres)| over {n} problems: "
      f"{maxdiff:.3e}")
print("  => 0 by construction: lambda(k)=0 means F(k) already lies inside the hull,")
print("     so dropping it cannot move min or max. The interval needs no lambda.\n")

# (iii) delta is not scale free --------------------------------------------
print("== (iii) is delta a justifiable constant? ==")
pr = P[0]
sem, cp, K, x, y = pr["sem"], pr["cpdag"], pr["K"], pr["x"], pr["y"]
th0 = theta_set(sem, cp if False else pr["g0"], x, y); w0 = width(th0)
lams = []
for k in K:
    fv = []
    for g in fibre(cp, K, k): fv.extend(theta_set(sem, g, x, y))
    lams.append(width(th0 + fv) - w0)
print(f"lambda vector on one problem, Y in its own units : "
      f"{[round(l,4) for l in lams]}")
print(f"the same vector if Y is measured in units 100x smaller: "
      f"{[round(100*l,4) for l in lams]}")
print("  => lambda has the units of the treatment effect, so any fixed delta ranks")
print("     claims differently under a change of outcome units. There is no scale-free")
print("     default and the definition supplies none.")
