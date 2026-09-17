#!/usr/bin/env python3
"""r_0 nas quatro redes gaussianas REAIS do corpus do Chris, a partir de Sigma-hat.

Tres versoes do mesmo raio, e a diferenca entre elas e' a metodologia:
  r_0 (oraculo)  -- Sigma populacional. Nao existe na pratica; e' a referencia.
  r_0 (pontual)  -- Sigma-hat, casco dos pontos. E' o que um analista ingenuo faz.
  r_0 (banda)    -- Sigma-hat, casco alargado por +-1,96 EP. E' o que ele deve fazer.

E uma decomposicao do mecanismo: o zero e' ESTRUTURAL (alguma extensao da efeito
exatamente nulo, isto e', nela X nem sequer e' causa de Y) ou de TRAVESSIA (o
casco cruza zero sem que nenhuma extensao o atinja)?
"""
import sys, json, time
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
import numpy as np
from rnull_real import (load, sem_of, dag_to_cpdag, apply_orientations, radii, first,
                        select_knowledge, component_of, optimal_adjustment_set_mpdag,
                        sample, pairs_for)

NETS = ("ecoli70", "magic-niab", "magic-irri")
NOBS = (200, 1000, 5000, 20000)
REPS = 5


def one(name, coverage, max_pairs=60, seed=20260909, jmax=4):
    rng = np.random.default_rng(seed)
    dag, W, var = load(name)
    sem = sem_of(dag, W, var)
    cpdag = dag_to_cpdag(dag)
    sigma_pop, idx_pop = sem.covariance(), {v: i for i, v in enumerate(dag.nodes)}
    K = sorted(select_knowledge(dag, cpdag, coverage))
    g0 = apply_orientations(cpdag, K) if K else cpdag
    if g0 is None:
        return []
    draws = {n: [np.cov(sample(dag, W, var, n, rng)[0], rowvar=False) for _ in range(REPS)]
             for n in NOBS}
    rows = []
    for (x, y) in pairs_for(dag, cpdag, max_pairs, rng=rng):
        Z = optimal_adjustment_set_mpdag(g0, x, y)
        if Z is None:
            continue
        comp = component_of(cpdag, x)
        Kloc = [e for e in K if component_of(cpdag, e[0]) == comp]
        if not Kloc:
            continue
        d = radii(cpdag, Kloc, x, y, Z, sigma_pop, idx_pop, jmax)
        if 0 not in d:
            continue
        r0 = first(d, "has_zero")
        rec = dict(network=name, coverage=coverage, x=x, y=y, comp=len(comp),
                   kloc=len(Kloc), r0=r0, neff0=d[0]["n_eff"],
                   structural=bool(d[r0]["structural_zero"]) if r0 is not None else None,
                   w1=(d[1]["hi"] - d[1]["lo"]) if 1 in d else None)
        for n in NOBS:
            pt, bd = [], []
            for sg in draws[n]:
                dn = radii(cpdag, Kloc, x, y, Z, sg, idx_pop, jmax, n=n)
                pt.append(first(dn, "has_zero"))
                bd.append(first(dn, "band_zero"))
            rec[f"pt_n{n}"] = pt
            rec[f"bd_n{n}"] = bd
        rows.append(rec)
    return rows


def agree(rows, key, r0key="r0"):
    """Fracao das REPLICAS em que o raio amostral concorda com o oraculo."""
    ok = tot = 0
    for r in rows:
        for v in r[key]:
            tot += 1
            ok += (v == r[r0key])
    return 100 * ok / tot if tot else float("nan")


def rate(rows, key):
    v = [a for r in rows for a in r[key]]
    return 100 * sum(1 for a in v if a is not None) / len(v) if v else float("nan")


if __name__ == "__main__":
    allrows = []
    for cov in (1.0, 0.5):
        print(f"\n{'='*104}\n  COBERTURA DO CONHECIMENTO = {cov}\n{'='*104}")
        print(f"{'rede':>12} {'inst':>5} {'r0<=orc':>8} {'estrut':>7} | "
              + " ".join(f"{'pt n='+str(n):>10}" for n in NOBS) + " | "
              + " ".join(f"{'bd n='+str(n):>10}" for n in NOBS))
        for name in NETS:
            t = time.time(); rows = one(name, cov)
            if not rows:
                print(f"{name:>12} — nenhuma instancia"); continue
            allrows += rows
            fin = [r for r in rows if r["r0"] is not None]
            st = 100 * sum(1 for r in fin if r["structural"]) / len(fin) if fin else float("nan")
            print(f"{name:>12} {len(rows):5d} {100*len(fin)/len(rows):7.1f}% {st:6.1f}% | "
                  + " ".join(f"{agree(rows,'pt_n'+str(n)):9.1f}%" for n in NOBS) + " | "
                  + " ".join(f"{agree(rows,'bd_n'+str(n)):9.1f}%" for n in NOBS)
                  + f"  [{time.time()-t:.0f}s]")
    json.dump(allrows, open("sweep_fast.json", "w"), indent=1)
    print("\nColunas pt/bd = % das replicas em que o raio calculado de Sigma-hat "
          "COINCIDE com o raio do oraculo.")
    print(f"{len(allrows)} instancias, {REPS} replicas por tamanho -> sweep_fast.json")
