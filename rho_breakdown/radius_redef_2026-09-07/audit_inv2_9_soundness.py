#!/usr/bin/env python3
"""inverse-2 audit, part 9: does the stated soundness theorem survive the PRIORITY ORDER?

Claimed: 'every claim not flagged load-bearing provably leaves Z valid in every DAG
extension of every admissible single revision of it.'
The partition is written as an ordered bullet list, so data-refuted and inert are
decided BEFORE the validity test ever runs.  Count the counterexamples.
"""
import sys
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import random_dag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag,
                                    is_valid_adjustment_set_mpdag)

def problem(rng):
    dag = random_dag(rng, 6, 0.35)
    if not dag.directed_edges: return None
    cpdag = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cpdag.undirected_edges)
    if not (2 <= len(und) <= 5): return None
    nodes = sorted(dag.nodes)
    cand = [(a,b) for a in nodes for b in nodes if a != b and b in dag.descendants(a)]
    rng.shuffle(cand)
    for x, y in cand:
        idx = rng.permutation(len(und))[:min(4, len(und))]
        K = []
        for t in idx:
            a, b = und[t]
            tr = (a,b) if (a,b) in dag.directed_edges else (b,a)
            K.append(tr if rng.random() > 0.2 else (tr[1], tr[0]))
        g0 = apply_orientations(cpdag, K)
        if g0 is None: continue
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None: continue
        return cpdag, K, x, y, frozenset(z), g0
    return None

rng = np.random.default_rng(999)
n = tried = 0
ce_dref = ce_inert = n_dref = n_inert = n_claims = 0
prob_with_ce = 0
while n < 300 and tried < 60000:
    tried += 1
    p = problem(rng)
    if p is None: continue
    cpdag, K, x, y, z, g0 = p
    n += 1; bad_here = False
    for k in K:
        n_claims += 1
        rest = [e for e in K if e != k]
        gr = apply_orientations(cpdag, rest)
        gv = apply_orientations(cpdag, rest + [(k[1], k[0])])
        adm = [g for g in (gr, gv) if g is not None]
        breaks = any(not is_valid_adjustment_set_mpdag(g, x, y, z) for g in adm)
        if gv is None:
            n_dref += 1
            if breaks: ce_dref += 1; bad_here = True     # NOT flagged, yet Z breaks
        elif gr is not None and gr == g0:
            n_inert += 1
            if breaks: ce_inert += 1; bad_here = True    # NOT flagged, yet Z breaks
    prob_with_ce += bad_here

print(f"problems {n}, claims {n_claims}")
print(f"data-refuted cell: {n_dref} claims; of those, an ADMISSIBLE single revision "
      f"(the retraction) still invalidates Z: {ce_dref}  "
      f"({100*ce_dref/max(n_dref,1):.1f}%)")
print(f"inert cell       : {n_inert} claims; of those, an ADMISSIBLE single revision "
      f"(the reversal) still invalidates Z: {ce_inert}  "
      f"({100*ce_inert/max(n_inert,1):.1f}%)")
print(f"problems where the stated soundness theorem is FALSE as written: "
      f"{prob_with_ce}/{n}  ({100*prob_with_ce/n:.1f}%)")
