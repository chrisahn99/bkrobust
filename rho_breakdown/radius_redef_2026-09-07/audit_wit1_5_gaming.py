#!/usr/bin/env python3
"""witness-1 audit 5: (a) can a verbose expert inflate R and empty L_flip?
   (b) does the budget-1 soundness lemma ever fail?"""
import sys
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import random_dag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, is_valid_adjustment_set_dag,
                                    optimal_adjustment_set_dag)
from audit_wit1_2_coverage import elicit, partition

# ---------- (a) VERBOSITY: same knowledge, more sentences ----------
print("=== (a) GAMING: does saying MORE (all of it true, no new information) shrink L_flip? ===")
rng = np.random.default_rng(4242)
res = {1: [], 2: [], 3: [], 4: [], 5: []}
tried, done = 0, 0
while done < 300 and tried < 200000:
    tried += 1
    dag = random_dag(rng, 7, 0.32)
    if not dag.directed_edges: continue
    cpdag = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cpdag.undirected_edges)
    if len(und) < 5: continue
    K_full = elicit(rng, dag, cpdag, len(und))   # every undirected edge, truthfully
    if K_full is None: continue
    g_full = apply_orientations(cpdag, K_full)
    if g_full is None: continue
    nodes = sorted(dag.nodes); cand = [(a,b) for a in nodes for b in nodes if a != b]
    rng.shuffle(cand)
    hit = None
    for x, y in cand:
        z = optimal_adjustment_set_mpdag(g_full, x, y)
        if z is None: continue
        hit = (x, y, frozenset(z)); break
    if hit is None: continue
    x, y, z = hit
    # G0 and Z are the SAME for every prefix that already pins the graph.
    ok = True
    row = {}
    for m in (1, 2, 3, 4, 5):
        if m > len(K_full): ok = False; break
        Km = K_full[:m]
        gm = apply_orientations(cpdag, Km)
        if gm is None: ok = False; break
        zm = optimal_adjustment_set_mpdag(gm, x, y)
        R, Lf, Lr, _, _ = partition(cpdag, Km, x, y, z)
        row[m] = (len(R), len(Lf), len(Lr), gm == g_full)
    if not ok: continue
    done += 1
    for m in row: res[m].append(row[m])

print(f"  {done} problems, K grown 1..5 claims over the SAME graph")
print(f"  {'|K|':>4s} {'mean|R|':>8s} {'mean|L_flip|':>13s} {'mean|L_ret|':>12s} {'%R of K':>9s} {'%G0 pinned':>11s}")
for m in (1,2,3,4,5):
    a = np.array(res[m])
    print(f"  {m:4d} {a[:,0].mean():8.2f} {a[:,1].mean():13.2f} {a[:,2].mean():12.2f} "
          f"{100*a[:,0].mean()/m:8.1f}% {100*a[:,3].mean():10.1f}%")
print("  -> R grows with sentence count on a FIXED graph: R measures verbosity, not refutability.")
print()

# ---------- (b) SOUNDNESS LEMMA STRESS TEST ----------
print("=== (b) SOUNDNESS: budget-1. If the false claim k is not in L_flip and not in R, is Z valid in the truth? ===")
rng = np.random.default_rng(99)
n, viol, in_R, in_Lf, clean_and_valid, clean_and_invalid = 0, 0, 0, 0, 0, 0
tried = 0
while n < 900 and tried < 400000:
    tried += 1
    dag = random_dag(rng, 7, 0.32)
    if not dag.directed_edges: continue
    cpdag = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cpdag.undirected_edges)
    if not (2 <= len(und) <= 6): continue
    K_true = elicit(rng, dag, cpdag, 3)
    if not K_true: continue
    i = int(rng.integers(len(K_true)))
    K_w = list(K_true); K_w[i] = (K_true[i][1], K_true[i][0])
    g0 = apply_orientations(cpdag, K_w)
    if g0 is None: continue
    nodes = sorted(dag.nodes); cand = [(a,b) for a in nodes for b in nodes if a != b]
    rng.shuffle(cand)
    for x, y in cand:
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None: continue
        z = frozenset(z); n += 1
        R, Lf, Lr, _, _ = partition(cpdag, K_w, x, y, z)
        k = K_w[i]
        zv = is_valid_adjustment_set_dag(dag, x, y, z)
        if k in R: in_R += 1
        elif k in Lf: in_Lf += 1
        else:
            if zv: clean_and_valid += 1
            else:  clean_and_invalid += 1; viol += 1
        break
print(f"  {n} budget-1 problems with exactly one false claim k")
print(f"    k landed in R      : {in_R:4d}  ({100*in_R/n:.1f}%)   [must be 0 for the R story to hold]")
print(f"    k landed in L_flip : {in_Lf:4d}  ({100*in_Lf/n:.1f}%)   [flagged: the audit fires]")
print(f"    k in neither, Z VALID   : {clean_and_valid:4d}  ({100*clean_and_valid/n:.1f}%)  [sound]")
print(f"    k in neither, Z INVALID : {clean_and_invalid:4d}  ({100*clean_and_invalid/n:.1f}%)  [LEMMA VIOLATION]")
print(f"  => soundness violations: {viol}")
