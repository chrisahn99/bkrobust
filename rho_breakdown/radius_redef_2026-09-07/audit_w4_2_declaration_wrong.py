#!/usr/bin/env python3
"""r_prod audit 2: what an UNFLIPPABLE claim that is false does to coverage.

r_prod hedges over S = {d : C(d) subset F, m(d) <= t}. If the analyst's false
claim happens to be one the protocol declared unflippable, the true DAG has
C(dag) NOT subset F, so the true DAG is EXCLUDED FROM THE HEDGE SET AT EVERY t.
No widening of t recovers it. This is the failure the proposal defers to a
"separate sensitivity row"; it is measured here because the paper's target
metric is width at HONEST 95% coverage, and a structurally excluded truth
cannot be calibrated away by scaling a hull that does not contain it.

Protocol per problem: truthful K, then reverse ONE claim c* -> consistent-but-
false K'. Compare the stratum where c* is flippable against the stratum where
the protocol declared it unflippable.
"""
import sys
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import make_problem
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, optimal_adjustment_set_dag,
                                    is_valid_adjustment_set_dag, adjusted_estimand, random_sem)


def hedge(cpdag, Kp, F, t, sem, x, y):
    """Hull of the adjusted estimand over S = {d : C(d) subset F, m(d) <= t}."""
    vals, contains_truth = [], False
    for d in enumerate_dag_extensions(cpdag):
        C = frozenset(c for c in Kp if not d.is_directed_edge(c[0], c[1]))
        if len(C) > t or not (C <= F):
            continue
        try:
            vals.append(adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(d, x, y)))
        except Exception:
            continue
    return vals


if __name__ == "__main__":
    rng = np.random.default_rng(20260907)
    got = tried = 0
    rows = []
    while got < 600 and tried < 60000:
        tried += 1
        pr = make_problem(rng)
        if pr is None:
            continue
        K, cpdag, dag, x, y = [tuple(e) for e in pr["K"]], pr["cpdag"], pr["dag"], pr["x"], pr["y"]
        if not K:
            continue
        # precondition (P): the truthful analyst's Z must be valid in the truth
        if not is_valid_adjustment_set_dag(dag, x, y, pr["z"]):
            continue
        # reverse exactly one claim -> consistent-but-false knowledge
        cstar = K[int(rng.integers(len(K)))]
        Kp = [((c[1], c[0]) if c == cstar else c) for c in K]
        g0 = apply_orientations(cpdag, Kp)
        if g0 is None:
            continue
        zp = optimal_adjustment_set_mpdag(g0, x, y)
        if zp is None:
            continue
        got += 1
        sem = random_sem(dag, rng)
        truth = sem.true_total_effect(x, y)
        point = adjusted_estimand(sem, x, y, zp)        # the naive analyst estimate
        FULL = frozenset(Kp)
        cbad = (cstar[1], cstar[0])                      # the false claim as stated in K'
        # arm A: protocol calls the false claim flippable  (F = K')
        # arm B: protocol calls the false claim UNflippable (F = K' \ {cbad})
        for arm, F in (("A_flippable", FULL), ("B_unflippable", FULL - {cbad})):
            v = hedge(cpdag, Kp, F, t=1, sem=sem, x=x, y=y)
            if not v:
                continue
            lo, hi = min(v), max(v)
            tol = 1e-9 * max(1.0, abs(truth))
            miss = max(0.0, max(lo - truth, truth - hi))
            rows.append(dict(arm=arm, lo=lo, hi=hi, w=hi - lo, truth=truth, point=point,
                             cov=bool(miss <= tol), miss=(0.0 if miss <= tol else miss)))

    print(f"problemas usaveis: {got} de {tried} sorteios (precondicao P aplicada)\n")
    for arm in ("A_flippable", "B_unflippable"):
        r = [d for d in rows if d["arm"] == arm]
        cov = np.mean([d["cov"] for d in r])
        w = np.mean([d["w"] for d in r])
        misses = [d["miss"] for d in r if not d["cov"]]
        print(f"{arm:15s} n={len(r):4d}  cobertura do casco = {cov:.3f}   largura media = {w:.4f}")
        if misses:
            print(f"{'':15s}   erro quando NAO cobre: mediana {np.median(misses):.4f}"
                  f"  p90 {np.quantile(misses,0.9):.4f}  max {max(misses):.4f}")
    A = [d for d in rows if d["arm"] == "A_flippable"]
    B = [d for d in rows if d["arm"] == "B_unflippable"]
    if A and B:
        print(f"\n economia de largura de B sobre A: "
              f"{100*(1 - np.mean([d['w'] for d in B])/np.mean([d['w'] for d in A])):.1f}%")
        print(f" custo em cobertura:               "
              f"{np.mean([d['cov'] for d in A]) - np.mean([d['cov'] for d in B]):+.3f}")
        # what a naive scale-to-95 would need in arm B
        need = [d["miss"] for d in B]
        q = float(np.quantile(need, 0.95))
        qA = float(np.quantile([d["miss"] for d in A], 0.95))
        print(f" (braco A precisaria somar ~{qA:.4f})")
        wB = np.mean([d["w"] for d in B])
        print(f"\n para levar o braco B a 95% escalando o casco: precisa somar ~{q:.4f} de cada lado")
        print(f" largura calibrada de B ~= {wB + 2*q:.4f}  vs  largura de A {np.mean([d['w'] for d in A]):.4f}")
