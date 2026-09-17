#!/usr/bin/env python3
"""inverse-2 audit, part 1: reproduce phi_1 on the demo scenarios exactly as defined."""
import sys
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
from bkrobust.demo.scenario import true_dag, scenarios, TREATMENT, OUTCOME
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag,
                                    is_valid_adjustment_set_mpdag,
                                    is_valid_adjustment_set_dag)
from bkrobust.demo.space import (enumerate_space, represented_dags, covering_pairs,
                                 neighbour_graph, bfs_distances)


def partition(cpdag, K, x, y, z, g0):
    """The four-way partition exactly as written in the definition (bullet order = priority)."""
    out = {}
    for k in K:
        rest = [e for e in K if e != k]
        g_ret = apply_orientations(cpdag, rest)
        g_rev = apply_orientations(cpdag, rest + [(k[1], k[0])])
        adm = [g for g in (g_ret, g_rev) if g is not None]
        if g_rev is None:
            cell = "data-refuted"
        elif g_ret is not None and g_ret == g0:
            cell = "inert"
        elif any(not is_valid_adjustment_set_mpdag(g, x, y, z) for g in adm):
            cell = "load-bearing"
        else:
            cell = "free"
        out[k] = dict(cell=cell, g_ret=g_ret, g_rev=g_rev, n_admissible=len(adm))
    lb = sum(1 for v in out.values() if v["cell"] == "load-bearing")
    den = sum(1 for v in out.values() if v["n_admissible"] >= 1)
    return out, (lb / den if den else float("nan")), lb, den


dag = true_dag()
cpdag = dag_to_cpdag(dag)
X, Y = TREATMENT, OUTCOME
space = enumerate_space(cpdag)
reps = represented_dags(space)
nbrs = neighbour_graph(space, covering_pairs(space, reps))
print(f"space size = {len(space)}")

for lab, sc in scenarios().items():
    K = [tuple(e) for e in sc["knowledge"]]
    g0 = apply_orientations(cpdag, K)
    z = optimal_adjustment_set_mpdag(g0, X, Y)
    z = frozenset(z) if z is not None else None
    dist = bfs_distances(nbrs, next(g for g in space if g == g0))
    cells, phi, lb, den = partition(cpdag, K, X, Y, z, g0)
    z_true_ok = is_valid_adjustment_set_dag(dag, X, Y, z)
    print(f"\n=== scenario {lab} | K={K}")
    print(f"    Z = O*(G0) = {sorted(z)}   valid in the TRUE dag? {z_true_ok}")
    print(f"    phi_1 = {lb}/{den} = {phi:.4f}")
    for k, v in cells.items():
        hops = []
        for tag, g in (("retract", v["g_ret"]), ("reverse", v["g_rev"])):
            if g is None:
                hops.append(f"{tag}=INADMISSIBLE")
            else:
                gg = next((h for h in space if h == g), None)
                d = dist.get(gg) if gg is not None else None
                ok = is_valid_adjustment_set_mpdag(g, X, Y, z)
                hops.append(f"{tag}=hop{d}/Zvalid={ok}")
        print(f"      {str(k):24s} {v['cell']:14s} {'  '.join(hops)}")
