#!/usr/bin/env python3
"""witness-1 audit 6: (a) retrieval quality of L_flip; (b) budget-2 silent failure;
   (c) does L_flip depend on the arbitrary choice of Z among valid sets of G0?"""
import sys
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import random_dag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, is_valid_adjustment_set_dag,
                                    all_valid_adjustment_sets_mpdag)
from audit_wit1_2_coverage import elicit, partition


def sweep(budget, target=700, seed=5):
    rng = np.random.default_rng(seed)
    out, tried = [], 0
    while len(out) < target and tried < 500000:
        tried += 1
        dag = random_dag(rng, 7, 0.32)
        if not dag.directed_edges: continue
        cpdag = dag_to_cpdag(dag)
        und = sorted(tuple(sorted(e)) for e in cpdag.undirected_edges)
        if not (2 <= len(und) <= 6): continue
        K_true = elicit(rng, dag, cpdag, 3)
        if not K_true or len(K_true) < budget: continue
        fi = set(rng.permutation(len(K_true))[:budget].tolist())
        K_w = [((k[1], k[0]) if i in fi else k) for i, k in enumerate(K_true)]
        false_claims = {K_w[i] for i in fi}
        g0 = apply_orientations(cpdag, K_w)
        if g0 is None: continue
        nodes = sorted(dag.nodes); cand = [(a,b) for a in nodes for b in nodes if a != b]
        rng.shuffle(cand)
        for x, y in cand:
            z = optimal_adjustment_set_mpdag(g0, x, y)
            if z is None: continue
            z = frozenset(z)
            R, Lf, Lr, _, _ = partition(cpdag, K_w, x, y, z)
            out.append(dict(zv=is_valid_adjustment_set_dag(dag, x, y, z),
                            Lf=set(Lf), R=set(R), Lr=set(Lr), F=false_claims,
                            m=len(K_w), g0=g0, cpdag=cpdag, x=x, y=y, K=K_w, dag=dag))
            break
    return out


for budget in (1, 2):
    rows = sweep(budget)
    n = len(rows)
    dmg = [r for r in rows if not r["zv"]]
    print(f"=== BUDGET {budget}: n={n}, real damage (Z invalid in truth) on {len(dmg)} ({100*len(dmg)/n:.1f}%) ===")
    # problem-level: does L_flip fire?
    fires = [r for r in rows if r["Lf"]]
    print(f"  L_flip fires on {100*len(fires)/n:.1f}% of problems; EMPTY on {100*(1-len(fires)/n):.1f}%")
    silent = [r for r in rows if not r["Lf"] and not r["zv"]]
    silent_full = [r for r in rows if not r["Lf"] and not r["R"] and not r["zv"]]
    print(f"  SILENT FAILURE (L_flip empty AND Z invalid) : {len(silent)}/{n} = {100*len(silent)/n:.2f}%")
    print(f"  ... of which R was also empty              : {len(silent_full)}  "
          f"(the rest were hidden inside R)")
    # retrieval scoring on the damaged stratum
    if dmg:
        rec = np.mean([len(r["Lf"] & r["F"]) / len(r["F"]) for r in dmg])
        rec_any = np.mean([bool(r["Lf"] & r["F"]) for r in dmg])
        print(f"  RECALL of false claims on damaged problems : {100*rec:.1f}% "
              f"(at least one caught: {100*rec_any:.1f}%)")
    prec = [len(r["Lf"] & r["F"]) / len(r["Lf"]) for r in rows if r["Lf"]]
    print(f"  PRECISION of L_flip (share of listed claims that are actually false): {100*np.mean(prec):.1f}%")
    fp = [r for r in rows if r["Lf"] and r["zv"]]
    print(f"  L_flip fires but Z is VALID in the truth  : {100*len(fp)/max(len(fires),1):.1f}% of firings")
    print()

# ---------- (c) Z-dependence ----------
print("=== (c) Is L_flip an artefact of the O* choice? Recompute over every valid Z of G0 ===")
rows = sweep(1, target=250, seed=11)
swings, sizes, agree = [], [], []
for r in rows:
    zs = all_valid_adjustment_sets_mpdag(r["g0"], r["x"], r["y"], max_size=3)
    zs = [frozenset(s) for s in zs]
    if len(zs) < 2: continue
    sets = []
    for z in zs[:8]:
        _, Lf, _, _, _ = partition(r["cpdag"], r["K"], r["x"], r["y"], z)
        sets.append(frozenset(Lf))
    swings.append(len(set(sets)))
    sizes.append(len(zs[:8]))
    agree.append(1.0 if len(set(sets)) == 1 else 0.0)
if swings:
    print(f"  {len(swings)} problems with >=2 valid adjustment sets")
    print(f"  mean distinct L_flip sets across the valid Z choices: {np.mean(swings):.2f} "
          f"(over {np.mean(sizes):.1f} sets tried)")
    print(f"  L_flip is the SAME for every valid Z on {100*np.mean(agree):.1f}% of problems")
    print(f"  => on {100*(1-np.mean(agree)):.1f}%, the 'load-bearing claims' the expert is shown "
          f"change if the analyst picks a different valid Z")
