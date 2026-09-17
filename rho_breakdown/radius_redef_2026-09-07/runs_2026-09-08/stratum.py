#!/usr/bin/env python3
"""A degenerescencia e da ESTATISTICA ou da POPULACAO?

Ate agora todas as medidas foram sobre todos os problemas sorteados. Mas ~75%
deles nao tem ignorancia nenhuma: a classe de equivalencia concorda sobre o
efeito, o casco e um ponto, e nenhuma estatistica de robustez tem o que medir.
Este script repete r_val e o raio de nulidade CONDICIONANDO no estrato onde a
classe de fato deixa incerteza.
"""
import sys, itertools
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
import numpy as np
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.evaluate import (random_sem, adjusted_estimand,
                                    is_valid_adjustment_set_mpdag,
                                    optimal_adjustment_set_mpdag)
from bkrobust.demo.space import (enumerate_space, represented_dags, covering_pairs,
                                 neighbour_graph, bfs_distances)
from confirm import problem, theta, hull, parents


def dist(c, tot, keys):
    return {str(k): (c.get(k, 0), 100 * c.get(k, 0) / tot if tot else 0) for k in keys}


def run(N=1800, seed=20260908, n_false=1):
    rng = np.random.default_rng(seed)
    rows, tried = [], 0
    while len(rows) < N and tried < 200000:
        tried += 1
        pr = problem(rng, n_false=n_false)
        if pr is None: continue
        dag, cpdag, x, y, K, g0 = (pr["dag"], pr["cpdag"], pr["x"], pr["y"],
                                   pr["K"], pr["g0"])
        sem = random_sem(dag, rng); cache = {}
        try:
            tau = adjusted_estimand(sem, x, y, parents(dag, x))
        except Exception:
            continue
        blanket = hull(theta(cpdag, sem, x, y, cache))
        if blanket is None: continue
        wb = blanket[1] - blanket[0]

        # r_val, no espaco em torno de G0
        z = optimal_adjustment_set_mpdag(g0, x, y)
        rv = None
        if z is not None:
            space = enumerate_space(cpdag)
            nbrs = neighbour_graph(space, covering_pairs(space, represented_dags(space)))
            src = next((g for g in space if g == g0), None)
            if src is not None:
                dd = bfs_distances(nbrs, src)
                bad = [dd[g] for g in space if dd.get(g) is not None and dd[g] > 0
                       and not is_valid_adjustment_set_mpdag(g, x, y, z)]
                rv = min(bad) if bad else None

        # raio de nulidade
        rn = None
        for j in range(len(K) + 1):
            vals = []
            for r in range(j + 1):
                for S in itertools.combinations(range(len(K)), r):
                    Kp = [K[i] for i in range(len(K)) if i not in S]
                    g = apply_orientations(cpdag, Kp) if Kp else cpdag
                    if g is None: continue
                    vals += theta(g, sem, x, y, cache)
            h = hull(vals)
            if h and h[0] <= 0.0 <= h[1]:
                rn = j; break
        rows.append(dict(wb=wb, rval=rv, rnull=rn, tau=tau))
    return rows, tried


def show(rows, label):
    n = len(rows)
    if not n:
        print(f"{label}: vazio"); return
    print(f"\n--- {label}  (n = {n}) ---")
    for name, key in (("r_val", "rval"), ("raio de nulidade", "rnull")):
        c = Counter("nunca" if r[key] is None else r[key] for r in rows)
        ks = sorted([k for k in c if k != "nunca"]) + (["nunca"] if "nunca" in c else [])
        interior = sum(v for k, v in c.items() if k not in ("nunca", 0, 1))
        line = "  ".join(f"{k}:{100*c[k]/n:4.1f}%" for k in ks)
        print(f"  {name:<18} {line}")
        import math
        ps = [v/n for v in c.values()]
        H = -sum(p*math.log2(p) for p in ps if p > 0)
        eff = 2**H
        big = sum(1 for p in ps if p >= 0.05)
        print(f"  {'':<18} entropia {H:.2f} bits | valores efetivos {eff:.2f} | "
              f"niveis com >=5% da massa: {big}")


if __name__ == "__main__":
    rows, tried = run()
    inf = [r for r in rows if r["wb"] > 1e-9]
    non = [r for r in rows if r["wb"] <= 1e-9]
    print(f"{len(rows)} problemas de {tried} sorteios | "
          f"informativos {len(inf)} ({100*len(inf)/len(rows):.1f}%)")
    show(rows, "TODOS os problemas")
    show(non, "NAO informativos (a classe ja concorda: nada a medir)")
    show(inf, "INFORMATIVOS (a classe deixa incerteza)")
