"""inverse-1: is r_alpha invariant to variables that have nothing to do with the query?
Add k isolated 2-node components to the demo DAG. They touch neither X, Y, Z nor any
elicited claim. r_val cannot move. Does r_alpha?"""
import sys, time
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
import numpy as np
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.scenario import TRUE_EDGES, K_TRUE, TREATMENT, OUTCOME, scenarios
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.space import (enumerate_space, represented_dags, covering_pairs,
                                 neighbour_graph, bfs_distances)
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag, is_valid_adjustment_set_mpdag

ALPHAS = [0.0, 0.05, 0.10, 0.25, 0.50]
K_A = list(scenarios()["A"]["knowledge"])

print(f"{'k isolated pairs':>17} {'|space|':>8} {'r_max':>6} {'r_val':>6}   " +
      "  ".join(f"r_a({a})" for a in ALPHAS) + "     sigma_global")
for k in range(0, 4):
    edges = list(TRUE_EDGES) + [(f"Iso{i}a", f"Iso{i}b") for i in range(k)]
    nodes = sorted({n for e in edges for n in e})
    dag = MPDAG(nodes, directed=edges)
    cp = dag_to_cpdag(dag)
    g0 = apply_orientations(cp, K_A)
    z = optimal_adjustment_set_mpdag(g0, TREATMENT, OUTCOME)
    t0 = time.time()
    space = enumerate_space(cp)
    reps = represented_dags(space); nbrs = neighbour_graph(space, covering_pairs(space, reps))
    g0k = next((g for g in space if g == g0), None)
    dist = bfs_distances(nbrs, g0k)
    rows = [(dist[g], bool(is_valid_adjustment_set_mpdag(g, TREATMENT, OUTCOME, z)))
            for g in space if g in dist]
    rmax = max(d for d, _ in rows)
    bad = [d for d, ok in rows if not ok and d > 0]
    rval = min(bad) if bad else None
    sig = {r: np.mean([ok for d, ok in rows if d <= r]) for r in range(rmax + 1)}
    ra = []
    for a in ALPHAS:
        S = [r for r in range(rmax + 1) if sig[r] >= 1 - a - 1e-12]
        ra.append(max(S) if S else -1)
    print(f"{k:17d} {len(rows):8d} {rmax:6d} {str(rval):>6}   " +
          "  ".join(f"{v:7d}" for v in ra) + f"     {sig[rmax]:.3f}   [{time.time()-t0:.0f}s]")
