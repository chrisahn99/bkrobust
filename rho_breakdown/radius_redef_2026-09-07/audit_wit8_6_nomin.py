#!/usr/bin/env python3
"""witness-8 audit 6: the price of the ANTICHAIN.  Compare the proposal's
minimal family L_b against the un-minimised U_b = {k : k lies in SOME A with
|A|<=b and lambda(A)>0}.  Two false claims, non-null queries."""
import sys, itertools
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
rng = np.random.default_rng(4242)
rows = []; got = tried = 0
while got < 300 and tried < 250000:
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
    cache = {}; base = theta(cpdag, K, sem, x, y, cache)
    if not base: continue
    got += 1
    bw = max(base) - min(base)
    sing = {}
    for k in K:
        v = set(base)
        for o in profiles(K, (k,)): v |= theta(cpdag, o, sem, x, y, cache)
        sing[k] = (max(v) - min(v)) - bw
    L1 = {k for k in K if sing[k] > TOL}
    pairpos = []
    for a, b in itertools.combinations(K, 2):
        v = set(base)
        for o in profiles(K, (a, b)): v |= theta(cpdag, o, sem, x, y, cache)
        if (max(v) - min(v)) - bw > TOL: pairpos.append((a, b))
    syn = [p for p in pairpos if sing[p[0]] <= TOL and sing[p[1]] <= TOL]   # minimal only
    Umin = L1 | {c for p in syn for c in p}                                 # U L_2  (proposta)
    Uall = L1 | {c for p in pairpos for c in p}                             # sem antichain
    def hedge(S):
        v = set(base)
        for o in profiles(K, tuple(sorted(S))): v |= theta(cpdag, o, sem, x, y, cache)
        return (max(v) - min(v), min(v) - 1e-9 <= truth <= max(v) + 1e-9)
    w1, c1 = hedge(L1); wm, cm = hedge(Umin); wa, ca = hedge(Uall); wK, cK = hedge(set(K))
    rows.append((w1, c1, wm, cm, wa, ca, wK, cK, len(Umin), len(Uall), len(K)))
R = np.array(rows, float); I = R[:, 6] > TOL
print(f"{got} problemas nao-nulos de {tried};  estrato informativo {int(I.sum())}/{got}")
print(f"|U L_2| medio {R[I,8].mean():.2f}   |U_2 sem antichain| medio {R[I,9].mean():.2f}   |K| medio {R[I,10].mean():.2f}")
for name, wi, ci in (("b=1 (so singletons)", 0, 1), ("b=2 com antichain (proposta)", 2, 3),
                     ("b=2 SEM antichain", 4, 5), ("cobertor (classe toda)", 6, 7)):
    print(f"  {name:30s} cobertura {R[I,ci].mean():.3f}   largura {R[I,wi].mean():.4f}"
          f"   larg/cobertor {np.mean(R[I,wi]/R[I,6]):.3f}")
