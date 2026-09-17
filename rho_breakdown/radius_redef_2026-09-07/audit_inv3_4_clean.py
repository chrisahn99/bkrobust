#!/usr/bin/env python3
"""inverse-3 audit 4: FAIR re-run.  Keep only problems where Z = O*(G0) really is
certified valid in G0, so nu(0) = 0 holds as the definition asserts.  Allow more
claims (m up to 6) so the polynomial has room.  Then ask again:
  - what resolution does nu have inside the r_val = 1 stratum?
  - is it a rescaled |B1|?
  - does it track damage?
  - what does the 'exact probability' guarantee cost when q is misdeclared?
"""
import sys, itertools, math
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from scipy.stats import spearmanr, kendalltau
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, is_valid_adjustment_set_mpdag,
                                    optimal_adjustment_set_dag, random_sem, adjusted_estimand)
from harness import random_dag, ball, r_val
from audit_inv3_1_core import corruption_table, nu_profile, nu, nu_cond

QS = [0.01, 0.05, 0.10, 0.20, 0.40]


def clean_problem(rng, n=7, p=0.32, kmax=5):
    dag = random_dag(rng, n, p)
    if not dag.directed_edges:
        return None
    cp = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cp.undirected_edges)
    if not (2 <= len(und) <= 6):
        return None
    nodes = sorted(dag.nodes)
    cand = [(x, y) for x in nodes for y in nodes if x != y]
    rng.shuffle(cand)
    for x, y in cand:
        km = min(kmax, len(und))
        idx = rng.permutation(len(und))[:km]
        K = []
        for t in idx:
            a, b = und[t]
            K.append((a, b) if (a, b) in dag.directed_edges else (b, a))
        g0 = apply_orientations(cp, K)
        if g0 is None:
            continue
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None:
            continue
        if not is_valid_adjustment_set_mpdag(g0, x, y, z):   # <-- THE FILTER
            continue
        return dict(dag=dag, cpdag=cp, x=x, y=y, K=[tuple(e) for e in K],
                    g0=g0, z=frozenset(z), undirected=und)
    return None


if __name__ == "__main__":
    rng = np.random.default_rng(20260907)
    got, tried, R = 0, 0, []
    while got < 400 and tried < 200000:
        tried += 1
        pr = clean_problem(rng)
        if pr is None:
            continue
        rb = ball(pr)
        if not rb:
            continue
        m, rows = corruption_table(pr)
        if m == 0:
            continue
        nu_j, a_j = nu_profile(m, rows)
        assert nu_j[0] == 0.0
        got += 1
        # damage
        sem = random_sem(pr["dag"], rng)
        worst = 0.0
        for j, adm, brk, J in rows:
            if not (adm and brk):
                continue
            KJ = [(b, a) if i in J else (a, b) for i, (a, b) in enumerate(pr["K"])]
            g = apply_orientations(pr["cpdag"], KJ)
            for d in enumerate_dag_extensions(g, limit=30):
                s2 = random_sem(d, rng)
                zt = optimal_adjustment_set_dag(d, pr["x"], pr["y"])
                truth = adjusted_estimand(s2, pr["x"], pr["y"], zt)
                est = adjusted_estimand(s2, pr["x"], pr["y"], pr["z"])
                if abs(truth) > 1e-9:
                    worst = max(worst, abs(est - truth) / abs(truth))
        R.append(dict(m=m, nu_j=nu_j, a_j=a_j, rval=r_val(rb), worst=worst,
                      B1=sum(1 for j, a, b, _ in rows if j == 1 and a and b),
                      nu={q: nu(q, m, nu_j) for q in QS},
                      nuc={q: nu_cond(q, m, rows) for q in QS}))

    print(f"problemas LIMPOS: {got} de {tried} sorteios;  |K| = {dict(Counter(r['m'] for r in R))}")
    c = Counter("INF" if r["rval"] is None else r["rval"] for r in R)
    print(f"r_val: {dict(c)}   (r_val=1 em {100*c.get(1,0)/got:.1f}%)")

    print("\n[1] RESOLUCAO com nu(0) = 0 garantido")
    for q in QS:
        v = np.array([r["nu"][q] for r in R])
        print(f"   q={q:<5} distintos {len(set(np.round(v,12))):3d}  =0 em {100*np.mean(v==0):5.1f}%  "
              f"med {np.median(v):.4f}  max {v.max():.4f}  teto 1-(1-q)^5 = {1-(1-q)**5:.4f}")

    S = [r for r in R if r["rval"] == 1]
    print(f"\n[2] ESTRATO r_val = 1 LIMPO (n = {len(S)})")
    for q in QS:
        v = np.array([r["nu"][q] for r in S])
        print(f"   q={q:<5} distintos {len(set(np.round(v,12))):3d}  min {v.min():.4f} "
              f"med {np.median(v):.4f} max {v.max():.4f}")
    print(f"   celulas (m,|B1|): {len(set((r['m'],r['B1']) for r in S))}")
    v = np.array([r["nu"][0.10] for r in S])
    key = [(r["m"], r["B1"]) for r in S]
    grp = {}
    for k, val in zip(key, v):
        grp.setdefault(k, []).append(val)
    within = sum(np.var(g)*len(g) for g in grp.values())/len(v)
    print(f"   var(nu(0.10)) explicada por (m,|B1|): {100*(1-within/np.var(v)):.2f}%")

    print("\n[3] nu vs |B1| e invariancia de ordem em q (global)")
    B1 = np.array([r["B1"] for r in R], float)
    for q in QS:
        v = np.array([r["nu"][q] for r in R])
        print(f"   q={q:<5} spearman(nu,|B1|) = {spearmanr(v,B1).statistic:.4f}")
    print(f"   kendall tau ordem q=0.01 vs q=0.40: "
          f"{kendalltau([r['nu'][0.01] for r in R],[r['nu'][0.40] for r in R]).statistic:.4f}")

    print("\n[4] nu vs DANO (o que o target metric realmente cobra)")
    w = np.array([r["worst"] for r in R]); v = np.array([r["nu"][0.10] for r in R])
    ok = w > 0
    print(f"   spearman(nu(0.10), pior vies rel) | ha dano: {spearmanr(v[ok],w[ok]).statistic:+.4f} (n={ok.sum()})")
    print(f"   spearman(|B1|, pior vies rel)    | ha dano: {spearmanr(B1[ok],w[ok]).statistic:+.4f}")
    hi = w > np.median(w[ok]) if ok.any() else None
    print(f"   pior vies rel: mediana {np.median(w[ok]):.3f} p90 {np.percentile(w[ok],90):.3f} max {w[ok].max():.2f}")
    # is nu small where damage is absent?  (requirement (c))
    nod = ~ok
    print(f"   nu(0.10) onde NAO ha dano medido: media {v[nod].mean():.4f} max {v[nod].max():.4f} (n={nod.sum()})")
    print(f"   nu(0.10) onde HA dano:            media {v[ok].mean():.4f} max {v[ok].max():.4f}")

    print("\n[5] BERNSTEIN vs PROBABILIDADE DECLARADA (a_j decai)")
    maxm = max(r["m"] for r in R)
    for j in range(maxm+1):
        vals = [r["a_j"][j] for r in R if len(r["a_j"]) > j]
        print(f"   a_{j}: media {np.mean(vals):.4f}")
    for q in (0.05, 0.10, 0.20, 0.40):
        a = np.array([r["nu"][q] for r in R]); b = np.array([r["nuc"][q] for r in R])
        nz = b > 1e-12
        print(f"   q={q}: dif media {np.mean(a-b):+.5f}  max |dif| {np.max(np.abs(a-b)):.5f}  "
              f"erro rel mediano (onde >0) {np.median(np.abs(a-b)[nz]/b[nz]):.4f}  "
              f"p90 {np.percentile(np.abs(a-b)[nz]/b[nz],90):.4f}")
