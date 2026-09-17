#!/usr/bin/env python3
"""witness-1 audit 4: is the L_flip screen a NO-OP for interval width, by theorem?

Claim under test: for k in L_ret \\ L_flip (dropped by the screen), the candidate
cl(K[k<-rev]) contributes the SAME estimand as the point estimate, hence zero width.
If so, L_flip vs L_ret is decorative for the target metric.
Also checks the demo hop counts asserted in the definition.
"""
import sys
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import random_dag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, is_valid_adjustment_set_mpdag,
                                    optimal_adjustment_set_dag, random_sem, adjusted_estimand)
from audit_wit1_2_coverage import elicit, partition, estimands

rng = np.random.default_rng(77)
same_z, same_val, diff_val, n_drop = 0, 0, 0, 0
tried, done = 0, 0
while done < 400 and tried < 200000:
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
    nodes = sorted(dag.nodes)
    cand = [(a,b) for a in nodes for b in nodes if a != b]
    rng.shuffle(cand)
    sem = random_sem(dag, rng)
    for x, y in cand:
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None: continue
        z = frozenset(z)
        R, Lf, Lr, grev, gret = partition(cpdag, K, x, y, z)
        dropped = [k for k in Lr if k not in Lf and k in grev]
        if not dropped: break
        done += 1
        point = adjusted_estimand(sem, x, y, z)
        for k in dropped:
            n_drop += 1
            g = grev[k]
            zz = optimal_adjustment_set_mpdag(g, x, y)
            if zz is not None and frozenset(zz) == z: same_z += 1
            vals = estimands(sem, x, y, [g])
            if all(abs(v - point) < 1e-9 for v in vals): same_val += 1
            else: diff_val += 1
        break

print("=== IS THE L_flip SCREEN A NO-OP FOR WIDTH? ===")
print(f"problems with a dropped candidate (k in L_ret \\ L_flip): {done}")
print(f"dropped candidates examined: {n_drop}")
print(f"  O*(cl(K[k<-rev])) identical to Z            : {same_z}/{n_drop}  ({100*same_z/max(n_drop,1):.1f}%)")
print(f"  candidate estimand identical to point est.  : {same_val}/{n_drop}  ({100*same_val/max(n_drop,1):.1f}%)")
print(f"  candidate estimand DIFFERENT (width gained)  : {diff_val}/{n_drop}  ({100*diff_val/max(n_drop,1):.1f}%)")
print()
print("=== DEMO HOP COUNTS ASSERTED IN THE DEFINITION ===")
from bkrobust.demo.scenario import scenarios, true_dag, TREATMENT, OUTCOME
from bkrobust.demo.space import (enumerate_space, represented_dags, covering_pairs,
                                 neighbour_graph, bfs_distances, atomic_moves)
dag = true_dag(); cp = dag_to_cpdag(dag)
space = enumerate_space(cp); reps = represented_dags(space)
covers = covering_pairs(space, reps); nbrs = neighbour_graph(space, covers)
for lbl in ("A", "B", "C"):
    K = [tuple(e) for e in scenarios()[lbl]["knowledge"]]
    g0 = apply_orientations(cp, K)
    g0m = next((g for g in space if g == g0), None)
    if g0m is None:
        print(f"  scenario {lbl}: G0 not in space"); continue
    d = bfs_distances(nbrs, g0m)
    print(f"  scenario {lbl}: K = {K}")
    for k in K:
        rest = [e for e in K if e != k]
        gr = apply_orientations(cp, rest)
        if gr is None:
            print(f"     retract {k[0]}->{k[1]}: FAIL"); continue
        grm = next((g for g in space if g == gr), None)
        hops = d.get(grm) if grm is not None else None
        rel = len(g0.directed_edges) - len(gr.directed_edges)
        print(f"     retract {k[0]}->{k[1]:6s}: hops={hops}  compelled edges released={rel}")
