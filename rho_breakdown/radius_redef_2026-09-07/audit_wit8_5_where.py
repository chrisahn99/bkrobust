#!/usr/bin/env python3
"""witness-8 audit 5: WHERE does b=1 fail with two false claims, and is that
exactly where the synergy pairs live?  Non-null queries only."""
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
rng = np.random.default_rng(777)
rows = []; got = tried = 0
while got < 300 and tried < 200000:
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
    badK = frozenset((k[1], k[0]) for k in bad)          # the false claims as the analyst wrote them
    cache = {}; base = theta(cpdag, K, sem, x, y, cache)
    if not base: continue
    got += 1
    bw = max(base) - min(base); sing = {}; parts = [set(base)]
    for k in K:
        v = set(base)
        for o in profiles(K, (k,)): v |= theta(cpdag, o, sem, x, y, cache)
        sing[k] = (max(v) - min(v)) - bw
        if sing[k] > TOL: parts.append(v)
    L1 = [k for k in K if sing[k] > TOL]; syn = []
    for a, b in itertools.combinations(K, 2):
        if sing[a] > TOL or sing[b] > TOL: continue
        v = set(base)
        for o in profiles(K, (a, b)): v |= theta(cpdag, o, sem, x, y, cache)
        if (max(v) - min(v)) - bw > TOL: syn.append((a, b)); parts.append(v)
    f1 = set().union(*parts[:1 + len(L1)]); f2 = set().union(*parts)
    vK = set(base)
    for o in profiles(K, tuple(K)): vK |= theta(cpdag, o, sem, x, y, cache)
    cov = lambda v: min(v) - 1e-9 <= truth <= max(v) + 1e-9
    rows.append(dict(nL1=len(L1), nsyn=len(syn), c1=cov(f1), c2=cov(f2), cK=cov(vK),
                     w1=max(f1)-min(f1), w2=max(f2)-min(f2), wK=max(vK)-min(vK),
                     rec1=len(badK & set(L1)), synhit=any(frozenset(p) == badK for p in syn),
                     info=(max(vK)-min(vK)) > TOL))
R = {k: np.array([r[k] for r in rows]) for k in rows[0]}
I = R["info"]
print(f"{got} problemas nao-nulos de {tried} sorteios;  estrato informativo {I.sum()}/{got}")
print(f"recall de L_1 sobre as 2 claims falsas (media de acertos): {R['rec1'][I].mean():.3f} de 2")
print(f"P(existe par sinergico) = {np.mean(R['nsyn'][I]>0):.3f}   "
      f"P(o par sinergico E exatamente o par falso) = {np.mean(R['synhit'][I]):.3f}")
for lab, m in (("L_1 vazio", I & (R['nL1']==0)), ("L_1 nao vazio", I & (R['nL1']>0)),
               ("sem par sinergico", I & (R['nsyn']==0)), ("com par sinergico", I & (R['nsyn']>0))):
    if m.sum()==0: continue
    print(f"  {lab:20s} n={m.sum():3d}  cob b=1 {R['c1'][m].mean():.3f}  cob b=2 {R['c2'][m].mean():.3f}"
          f"  cob cobertor {R['cK'][m].mean():.3f}  |  larg b=1 {R['w1'][m].mean():.3f}"
          f"  b=2 {R['w2'][m].mean():.3f}  cobertor {R['wK'][m].mean():.3f}")
f = I & ~R['c1']
print(f"\nquando b=1 FALHA (n={f.sum()}): P(par sinergico existe) = "
      f"{np.mean(R['nsyn'][f]>0) if f.sum() else float('nan'):.3f}   "
      f"P(b=2 conserta) = {np.mean(R['c2'][f]) if f.sum() else float('nan'):.3f}")
