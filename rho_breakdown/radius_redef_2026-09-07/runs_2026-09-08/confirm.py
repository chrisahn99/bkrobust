#!/usr/bin/env python3
"""Independent confirmation of the two numbers that decide the paper.

Q1  rho* : the width a SELECTIVE hedge (retraction hull, budget j) needs for
           95% coverage, relative to the BLANKET hedge over the whole class.
Q2  How much of that hedged ignorance is a dispute about whether X causes Y at
    all, rather than about which confounders to adjust for.

Everything is exact: population covariance, no data sampling. The interval is the
hull of possible total effects under the IDA convention.
"""
import sys, itertools
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
import numpy as np

from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import random_sem, adjusted_estimand
from harness import random_dag


def parents(dag, v):
    return {a for a, b in dag.directed_edges if b == v}


def theta(g, sem, x, y, cache):
    """IDA: the set of total effects over every DAG extension of g."""
    k = g.key() if hasattr(g, "key") else g
    if k in cache:
        return cache[k]
    vals = []
    for d in enumerate_dag_extensions(g):
        try:
            vals.append(adjusted_estimand(sem, x, y, parents(d, x)))
        except Exception:
            pass
    cache[k] = vals
    return vals


def hull(vals):
    return (min(vals), max(vals)) if vals else None


def problem(rng, n=7, p=0.3, kmax=4, n_false=0):
    """n_false claims are asserted REVERSED relative to the truth."""
    dag = random_dag(rng, n, p)
    if not dag.directed_edges:
        return None
    cpdag = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cpdag.undirected_edges)
    if not (2 <= len(und) <= 5):
        return None
    nodes = sorted(dag.nodes)
    pairs = [(a, b) for a in nodes for b in nodes if a != b]
    rng.shuffle(pairs)
    # a query whose endpoints are connected in the skeleton
    sk = {frozenset(e) for e in dag.directed_edges}
    for x, y in pairs:
        idx = rng.permutation(len(und))[:min(kmax, len(und))]
        K, truth = [], []
        for t in idx:
            a, b = und[t]
            e = (a, b) if (a, b) in dag.directed_edges else (b, a)
            truth.append(e); K.append(e)
        wrong = rng.permutation(len(K))[:min(n_false, len(K))]
        for w in wrong:
            K[w] = (K[w][1], K[w][0])            # afirma o inverso
        g0 = apply_orientations(cpdag, K)
        if g0 is None:
            continue                              # Meek rejeita: erro pego de graca
        return dict(dag=dag, cpdag=cpdag, x=x, y=y, K=K, g0=g0,
                    n_false=len(wrong))
    return None


def run(N=400, seed=20260907, n_false=0):
    rng = np.random.default_rng(seed)
    rows, tried = [], 0
    while len(rows) < N and tried < 40000:
        tried += 1
        pr = problem(rng, n_false=n_false)
        if pr is None:
            continue
        dag, cpdag, x, y, K = pr["dag"], pr["cpdag"], pr["x"], pr["y"], pr["K"]
        sem = random_sem(dag, rng)
        cache = {}
        try:
            tau = adjusted_estimand(sem, x, y, parents(dag, x))
        except Exception:
            continue

        blanket = hull(theta(cpdag, sem, x, y, cache))
        if blanket is None:
            continue
        wb = blanket[1] - blanket[0]

        # the retraction hull at budget j: union over subsets of K of size <= j
        hulls = {}
        for j in range(0, min(3, len(K)) + 1):
            vals = []
            for r in range(j + 1):
                for S in itertools.combinations(range(len(K)), r):
                    Kp = [K[i] for i in range(len(K)) if i not in S]
                    g = apply_orientations(cpdag, Kp) if Kp else cpdag
                    if g is None:
                        continue
                    vals += theta(g, sem, x, y, cache)
            hulls[j] = hull(vals)

        # Q2: does every extension of the CPDAG agree that Y descends from X?
        exts = enumerate_dag_extensions(cpdag)
        agree = all(y in d.descendants(x) for d in exts)
        # the width of the hull restricted to worlds where Y IS a descendant
        dvals = []
        for d in exts:
            if y in d.descendants(x):
                try: dvals.append(adjusted_estimand(sem, x, y, parents(d, x)))
                except Exception: pass
        dh = hull(dvals)
        wd = (dh[1] - dh[0]) if dh else 0.0

        rows.append(dict(tau=tau, wb=wb, blanket=blanket, hulls=hulls,
                         agree=agree, wdir=wd, nK=len(K)))
    return rows, tried


def report(rows, tried):
    N = len(rows)
    print(f"problemas: {N} (de {tried} sorteios)\n")
    inf = [r for r in rows if r["wb"] > 1e-9]        # a classe deixa alguma ignorancia
    print(f"informativos (largura blanket > 0): {len(inf)}  ({100*len(inf)/N:.1f}%)\n")

    print("=" * 66)
    print("Q1  COBERTURA E LARGURA POR ORCAMENTO DE RETRACAO")
    print("=" * 66)
    print(f"{'orcamento':>10} {'cobre':>8} {'largura/blanket':>17} {'= blanket':>11}")
    for j in range(4):
        hs = [r for r in rows if j in r["hulls"] and r["hulls"][j]]
        if not hs: continue
        cov = np.mean([h["hulls"][j][0] - 1e-12 <= h["tau"] <= h["hulls"][j][1] + 1e-12 for h in hs])
        ii = [h for h in hs if h["wb"] > 1e-9]
        ratio = np.mean([(h["hulls"][j][1]-h["hulls"][j][0]) / h["wb"] for h in ii]) if ii else float("nan")
        atom = np.mean([abs((h["hulls"][j][1]-h["hulls"][j][0]) - h["wb"]) < 1e-9 for h in ii]) if ii else float("nan")
        print(f"{j:>10} {cov:>8.3f} {ratio:>17.3f} {100*atom:>10.1f}%")
    covb = np.mean([r["blanket"][0]-1e-12 <= r["tau"] <= r["blanket"][1]+1e-12 for r in rows])
    print(f"{'blanket':>10} {covb:>8.3f} {1.0:>17.3f} {100.0:>10.1f}%")

    print("\n" + "=" * 66)
    print("Q2  QUANTO DA IGNORANCIA E DISPUTA DE DIRECAO, NAO DE AJUSTE")
    print("=" * 66)
    dis = [r for r in inf if not r["agree"]]
    print(f"informativos onde a classe DISCORDA se X causa Y: {len(dis)}/{len(inf)}"
          f"  ({100*len(dis)/len(inf):.1f}%)")
    share = [r["wdir"]/r["wb"] for r in inf]
    print(f"largura restrita aos mundos em que Y descende de X, "
          f"em fracao da blanket: media {np.mean(share):.3f}")
    print(f"   e EXATAMENTE ZERO em {100*np.mean([s < 1e-9 for s in share]):.1f}% dos informativos")


if __name__ == "__main__":
    import sys
    for nf in (0, 1, 2):
        print("\n" + "#" * 66)
        print(f"#  ARM: {nf} AFIRMACOES DO ESPECIALISTA ESTAO INVERTIDAS")
        print("#" * 66)
        rows, tried = run(N=300, n_false=nf)
        report(rows, tried)
