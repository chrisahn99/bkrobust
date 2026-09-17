#!/usr/bin/env python3
"""witness-8 FINAL 2: (a) monotonicity of lambda, (b) the ACTUAL single-metric
number: width required for honest 95% coverage after scaling."""
import sys, itertools, collections
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import make_problem
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import optimal_adjustment_set_dag, adjusted_estimand, random_sem
TOL = 1e-9
def theta(cpdag, o, sem, x, y, c):
    k = frozenset(o)
    if k in c: return c[k]
    g = apply_orientations(cpdag, list(o))
    c[k] = frozenset() if g is None else frozenset(np.round(
        [adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(d, x, y))
         for d in enumerate_dag_extensions(g)], 12))
    return c[k]
def profiles(K, A):
    rest = [k for k in K if k not in A]; A = list(A)
    for b in itertools.product((0, 1), repeat=len(A)):
        yield rest + [(a[1], a[0]) for a, t in zip(A, b) if t == 1]
def hull(cpdag, K, A, sem, x, y, c):
    v = set(theta(cpdag, K, sem, x, y, c))
    for o in profiles(K, tuple(A)): v |= theta(cpdag, o, sem, x, y, c)
    return v
rng = np.random.default_rng(311)
mono_tests = mono_viol = 0
rows = []; got = tried = 0
while got < 300 and tried < 300000:
    tried += 1
    pr = make_problem(rng, n=7, p=0.35, k_claims=5)
    if pr is None: continue
    K0 = [tuple(e) for e in pr["K"]]
    if len(K0) < 4: continue
    cpdag, dag, x, y = pr["cpdag"], pr["dag"], pr["x"], pr["y"]
    sem = random_sem(dag, rng); truth = sem.true_total_effect(x, y)
    if abs(truth) < 1e-9: continue
    bad = [K0[i] for i in rng.permutation(len(K0))[:2]]
    K = [((k[1], k[0]) if k in bad else k) for k in K0]
    if apply_orientations(cpdag, K) is None: continue
    c = {}; base = theta(cpdag, K, sem, x, y, c)
    if not base: continue
    got += 1; w0 = max(base) - min(base)
    lam = {}
    for k in K:
        v = hull(cpdag, K, (k,), sem, x, y, c); lam[(k,)] = (max(v) - min(v)) - w0
    U1 = {k for k in K if lam[(k,)] > TOL}
    U2 = set(U1)
    for a, b in itertools.combinations(K, 2):
        v = hull(cpdag, K, (a, b), sem, x, y, c); L = (max(v) - min(v)) - w0
        lam[(a, b)] = L
        # monotonicity: lambda(pair) >= max of its singletons
        mono_tests += 2
        if L + 1e-12 < lam[(a,)]: mono_viol += 1
        if L + 1e-12 < lam[(b,)]: mono_viol += 1
        if L > TOL: U2 |= {a, b}
    def iv(S):
        v = hull(cpdag, K, tuple(sorted(S)), sem, x, y, c); return min(v), max(v)
    out = []
    for S in (U1, U2, set(K)):
        lo, hi = iv(S); out += [(lo + hi) / 2, (hi - lo) / 2]
    rows.append(out + [truth])
R = np.array(rows, float)
INF = R[:, 5] > TOL                      # blanket half-width > 0
print(f"monotonicity of lambda: {mono_viol} violations in {mono_tests} pair/singleton tests")
print(f"non-null {got}/{tried};  informative (blanket width>0) {int(INF.sum())}")
def calib(ci, hi, mask, name):
    ctr, hw, tr = R[mask, ci], R[mask, hi], R[mask, 6]
    grid = np.concatenate([np.linspace(0, 5, 2001), np.linspace(5, 200, 400)])
    for s in grid:
        cov = np.mean(np.abs(tr - ctr) <= s * hw + 1e-9)
        if cov >= 0.95:
            print(f"  {name:24s} s*={s:7.3f}  cov {cov:.3f}  "
                  f"calibrated mean width {np.mean(2*s*hw):.4f}")
            return
    cov = np.mean(np.abs(tr - ctr) <= 200 * hw + 1e-9)
    print(f"  {name:24s} INADMISSIBLE: coverage saturates at {cov:.3f} "
          f"(scaling cannot reach 0.95; zero-width intervals stay zero)")
print("\nWIDTH REQUIRED FOR HONEST 95% COVERAGE (informative stratum):")
for nm, ci, hi in (("b=1  (U_1)", 0, 1), ("b=2  (U_2)", 2, 3), ("blanket", 4, 5)):
    calib(ci, hi, INF, nm)
print("\nsame, over ALL non-null draws (blanket width may be 0):")
ALL = np.ones(len(R), bool)
for nm, ci, hi in (("b=1  (U_1)", 0, 1), ("b=2  (U_2)", 2, 3), ("blanket", 4, 5)):
    calib(ci, hi, ALL, nm)
print(f"\nzero-width rate: b=1 {np.mean(R[INF,1]<=TOL):.3f}   b=2 {np.mean(R[INF,3]<=TOL):.3f}"
      f"   blanket 0.000   (informative stratum)")

print("\ncoverage vs scale s (informative stratum):")
for nm, ci, hi in (("b=2 (U_2)", 2, 3), ("blanket", 4, 5)):
    line = "  ".join(f"s={s:.2f}:{np.mean(np.abs(R[INF,6]-R[INF,ci])<=s*R[INF,hi]+1e-9):.3f}"
                     for s in (0.25, 0.5, 0.75, 0.9, 0.99, 1.0))
    print(f"  {nm:12s} {line}")
