#!/usr/bin/env python3
"""r_prod audit 1: what does escaping degeneracy actually COST in declarations,
and where does r_prod land once you have paid it?

r_prod = min{ m(d) : d in [Chat], C(d) subset F, Z invalid in d }.
Since the restriction only DELETES candidate breakers, r_prod >= r_K always.
To get r_prod > 1 the protocol must declare EVERY single-claim breaker unflippable.
So the price of escaping the r=1 artefact is exactly |B1|, the number of distinct
single-claim breakers, and the honest denominator is |K|.
"""
import sys
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import make_problem
from bkrobust.demo.meek import enumerate_dag_extensions
from bkrobust.demo.evaluate import is_valid_adjustment_set_dag


def profile(pr):
    """Per-DAG record: contradicted-claim set C(d), m(d), validity of Z."""
    K = [tuple(e) for e in pr["K"]]
    x, y, z = pr["x"], pr["y"], pr["z"]
    ext = enumerate_dag_extensions(pr["cpdag"])
    recs = []
    for d in ext:
        C = frozenset(c for c in K if not d.is_directed_edge(c[0], c[1]))
        recs.append((C, len(C), bool(is_valid_adjustment_set_dag(d, x, y, z)), d))
    return K, recs


def r_restricted(recs, F):
    """min m(d) over invalid d whose contradicted set lies inside F. None = infinity."""
    ms = [m for C, m, ok, _ in recs if not ok and C <= F]
    return min(ms) if ms else None


if __name__ == "__main__":
    rng = np.random.default_rng(20260907)
    got, tried = 0, 0
    rk_dist, b1_size, b1_frac, after_escape = Counter(), [], [], Counter()
    near_frac, k_sizes = [], []
    escape_all = 0
    while got < 400 and tried < 40000:
        tried += 1
        pr = make_problem(rng)
        if pr is None:
            continue
        K, recs = profile(pr)
        if not K:
            continue
        got += 1
        FULL = frozenset(K)
        rk = r_restricted(recs, FULL)
        rk_dist["INF" if rk is None else rk] += 1
        k_sizes.append(len(K))
        if rk != 1:
            continue
        # B1: the single claims whose reversal alone invalidates Z
        B1 = {next(iter(C)) for C, m, ok, _ in recs if not ok and m == 1}
        b1_size.append(len(B1))
        b1_frac.append(len(B1) / len(K))
        if len(B1) == len(K):
            escape_all += 1
        # touching the query?
        near = sum(1 for c in B1 if pr["x"] in c or pr["y"] in c)
        near_frac.append(near / len(B1))
        # pay the minimal price: F = K \ B1. Where does r_prod land?
        F = FULL - B1
        r2 = r_restricted(recs, F)
        after_escape["INF" if r2 is None else r2] += 1

    print(f"problemas: {got} de {tried} sorteios | |K| medio {np.mean(k_sizes):.2f}")
    print("\ndistribuicao de r_K (alfabeto completo, F = K):")
    for k in sorted(rk_dist, key=lambda v: (v == "INF", v)):
        print(f"   r_K = {str(k):<5} {rk_dist[k]:4d}  {100*rk_dist[k]/got:5.1f}%")
    n1 = len(b1_size)
    print(f"\nnos {n1} problemas com r_K = 1 ({100*n1/got:.1f}%):")
    print(f"   |B1| (breakers de 1 claim)  media {np.mean(b1_size):.2f}  mediana {np.median(b1_size):.0f}  max {max(b1_size)}")
    print(f"   |B1|/|K| = fracao de K que o protocolo TEM de declarar inflipavel:")
    print(f"        media {np.mean(b1_frac):.3f}   mediana {np.median(b1_frac):.3f}")
    print(f"        == 1.0 (todo K inflipavel) em {escape_all}/{n1} = {100*escape_all/n1:.1f}%")
    print(f"   fracao de B1 que toca X ou Y: media {np.mean(near_frac):.3f}")
    print("\n   r_prod DEPOIS de pagar o preco minimo (F = K \\ B1):")
    for k in sorted(after_escape, key=lambda v: (v == "INF", v)):
        print(f"        r_prod = {str(k):<5} {after_escape[k]:4d}  {100*after_escape[k]/n1:5.1f}%")
    inf = after_escape.get("INF", 0)
    print(f"\n=> escapar de r=1 leva a r_prod = INFINITO em {100*inf/n1:.1f}% dos casos.")
    print("   r_prod nao e continuo: e {1, ...} -> {1} ou {INF}. Ver audit 2 para a largura.")
