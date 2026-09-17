#!/usr/bin/env python3
"""O teste que nunca foi feito: o hedge LOCAL contra o hedge UNIFORME.

Toda medicao ate agora usou RH_j, que retrai QUALQUER subconjunto de j
afirmacoes. Mas a tese empirica da linha e que o dano e LOCAL: uma afirmacao
que nao toca a consulta quase nunca estraga o conjunto de ajuste. Se isso e
verdade, um hedge restrito as afirmacoes incidentes a consulta deveria comprar
a MESMA cobertura por MENOS largura. Isso nunca foi testado diretamente.

Politicas comparadas, todas reportando pelo mesmo caminho:
  PLAIN    nenhum hedge
  RH_1     retrai qualquer 1 afirmacao          (uniforme)
  RH_2     retrai qualquer 2                    (uniforme)
  LOCAL    retrai qualquer subconjunto das afirmacoes que TOCAM {X,Y}
  FAR      controle negativo: idem, mas as que NAO tocam
  BLANKET  a classe inteira
"""
import sys, itertools, json
from collections import defaultdict
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
import numpy as np
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.evaluate import random_sem, adjusted_estimand
from confirm import problem, theta, hull, parents


def hull_over(cpdag, K, subsets, sem, x, y, cache):
    vals = []
    for S in subsets:
        Kp = [K[i] for i in range(len(K)) if i not in S]
        g = apply_orientations(cpdag, Kp) if Kp else cpdag
        if g is None: continue
        vals += theta(g, sem, x, y, cache)
    return hull(vals)


def subsets_upto(n, j):
    return [S for r in range(j + 1) for S in itertools.combinations(range(n), r)]


def all_subsets_of(idx):
    idx = list(idx)
    return [tuple(sorted(S)) for r in range(len(idx) + 1)
            for S in itertools.combinations(idx, r)]


def run(N=600, seed=20260908, n_false=1):
    rng = np.random.default_rng(seed)
    rows, tried = [], 0
    while len(rows) < N and tried < 80000:
        tried += 1
        pr = problem(rng, n_false=n_false)
        if pr is None: continue
        dag, cpdag, x, y, K = pr["dag"], pr["cpdag"], pr["x"], pr["y"], pr["K"]
        sem = random_sem(dag, rng); cache = {}
        try:
            tau = adjusted_estimand(sem, x, y, parents(dag, x))
        except Exception:
            continue
        blanket = hull(theta(cpdag, sem, x, y, cache))
        if blanket is None or blanket[1] - blanket[0] <= 1e-9:
            continue                                  # so o estrato informativo
        wb = blanket[1] - blanket[0]

        near = [i for i, e in enumerate(K) if x in e or y in e]
        far = [i for i in range(len(K)) if i not in near]

        pol = {}
        pol["PLAIN"]   = hull_over(cpdag, K, [()], sem, x, y, cache)
        pol["RH_1"]    = hull_over(cpdag, K, subsets_upto(len(K), 1), sem, x, y, cache)
        pol["RH_2"]    = hull_over(cpdag, K, subsets_upto(len(K), 2), sem, x, y, cache)
        pol["LOCAL"]   = hull_over(cpdag, K, all_subsets_of(near), sem, x, y, cache)
        pol["FAR"]     = hull_over(cpdag, K, all_subsets_of(far), sem, x, y, cache)
        # === orcamento CASADO: retrair no maximo UMA, mas so as proximas / so as longes
        pol["RH1_near"] = hull_over(cpdag, K, [()] + [(i,) for i in near], sem, x, y, cache)
        pol["RH1_far"]  = hull_over(cpdag, K, [()] + [(i,) for i in far], sem, x, y, cache)
        pol["BLANKET"] = blanket
        if any(v is None for v in pol.values()):
            continue
        rows.append(dict(tau=tau, wb=wb, pol=pol, n_near=len(near), nK=len(K)))
    return rows, tried


if __name__ == "__main__":
    for nf in (1, 2):
        rows, tried = run(n_false=nf)
        n = len(rows)
        print(f"\n{'='*70}\n  {nf} afirmacao(oes) invertida(s)  |  {n} problemas informativos "
              f"(de {tried} sorteios)\n{'='*70}")
        print(f"{'politica':<9} {'cobre':>7} {'larg/blanket':>13} {'= blanket':>10} "
              f"{'retracoes':>10}")
        sizes = {"PLAIN": 1, "RH_1": None, "RH_2": None, "LOCAL": None, "FAR": None,
                 "BLANKET": None}
        for name in ("PLAIN", "RH1_far", "RH1_near", "RH_1", "LOCAL", "RH_2", "BLANKET"):
            cov = np.mean([r["pol"][name][0]-1e-12 <= r["tau"] <= r["pol"][name][1]+1e-12
                           for r in rows])
            w = np.mean([(r["pol"][name][1]-r["pol"][name][0]) / r["wb"] for r in rows])
            atom = np.mean([abs((r["pol"][name][1]-r["pol"][name][0]) - r["wb"]) < 1e-9
                            for r in rows])
            if name == "LOCAL":
                sz = np.mean([2**r["n_near"] for r in rows])
            elif name == "FAR":
                sz = np.mean([2**(r["nK"]-r["n_near"]) for r in rows])
            elif name == "RH_1":
                sz = np.mean([1+r["nK"] for r in rows])
            elif name == "RH1_near":
                sz = np.mean([1+r["n_near"] for r in rows])
            elif name == "RH1_far":
                sz = np.mean([1+r["nK"]-r["n_near"] for r in rows])
            elif name == "RH_2":
                sz = np.mean([1+r["nK"]+r["nK"]*(r["nK"]-1)/2 for r in rows])
            else:
                sz = 1.0
            mark = " <-- orcamento casado" if name in ("RH1_near","RH1_far") else ""
            print(f"{name:<9} {cov:>7.3f} {w:>13.3f} {100*atom:>9.1f}% {sz:>10.1f}{mark}")
        print(f"   afirmacoes que tocam a consulta: media "
              f"{np.mean([r['n_near'] for r in rows]):.2f} de {np.mean([r['nK'] for r in rows]):.2f}")
