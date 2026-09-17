#!/usr/bin/env python3
"""Confirm the no-op as a THEOREM: k not in L_flip => the reversed candidate's
adjusted estimand equals the point estimate, in EVERY DAG of the reversed class."""
import sys
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import random_dag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, is_valid_adjustment_set_mpdag,
                                    optimal_adjustment_set_dag, random_sem, adjusted_estimand)
from audit_wit1_2_coverage import elicit

rng = np.random.default_rng(2026)
checked, bad, n_ext = 0, 0, 0
tried = 0
while checked < 1200 and tried < 300000:
    tried += 1
    dag = random_dag(rng, 7, 0.32)
    if not dag.directed_edges: continue
    cpdag = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cpdag.undirected_edges)
    if not (2 <= len(und) <= 6): continue
    K = elicit(rng, dag, cpdag, 3)
    if not K: continue
    g0 = apply_orientations(cpdag, K)
    if g0 is None: continue
    nodes = sorted(dag.nodes); cand = [(a,b) for a in nodes for b in nodes if a != b]
    rng.shuffle(cand)
    sem = random_sem(dag, rng)
    for x, y in cand:
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None: continue
        z = frozenset(z)
        point = adjusted_estimand(sem, x, y, z)
        for k in K:
            rest = [e for e in K if e != k]
            g = apply_orientations(cpdag, rest + [(k[1], k[0])])
            if g is None: continue
            if not is_valid_adjustment_set_mpdag(g, x, y, z):
                continue          # k IS in L_flip; theorem says nothing
            checked += 1
            for d in enumerate_dag_extensions(g):
                n_ext += 1
                v = adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(d, x, y))
                if abs(v - point) > 1e-9:
                    bad += 1
        break

print("THEOREM: k not in L_flip and not in R  =>  every DAG-extension estimand of the")
print("         reversed class equals the point estimate (zero width contribution).")
print(f"  screened-out claims checked : {checked}")
print(f"  DAG extensions evaluated    : {n_ext}")
print(f"  counterexamples             : {bad}")
