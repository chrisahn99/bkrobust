"""inverse-1 (alpha-tolerant safe radius) : do the proposal's DEMO numbers reproduce?
Claimed: scenario A -> r_alpha = 2,2,3,4,9 for alpha = 0,.05,.10,.25,.50 against r_val = 3
         scenario C -> r_alpha = 1,1,2,3,9 against r_val = 2
         and r_0 = r_val - 1 exactly on all three scenarios.
"""
import sys
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
import numpy as np
from bkrobust.demo.scenario import true_dag, scenarios, TREATMENT, OUTCOME
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.space import (enumerate_space, represented_dags, covering_pairs,
                                 neighbour_graph, bfs_distances)
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag,
                                    is_valid_adjustment_set_mpdag)

dag = true_dag(); cp = dag_to_cpdag(dag)
space = enumerate_space(cp); reps = represented_dags(space)
covers = covering_pairs(space, reps); nbrs = neighbour_graph(space, covers)
print(f"space size = {len(space)}")

ALPHAS = [0.0, 0.05, 0.10, 0.25, 0.50]

for lab, sc in scenarios().items():
    K = list(sc["knowledge"])
    g0 = apply_orientations(cp, K)
    z = optimal_adjustment_set_mpdag(g0, TREATMENT, OUTCOME)
    g0k = next((g for g in space if g == g0), None)
    dist = bfs_distances(nbrs, g0k)
    rows = [(dist[g], bool(is_valid_adjustment_set_mpdag(g, TREATMENT, OUTCOME, z)))
            for g in space if g in dist]
    rmax = max(d for d, _ in rows)
    bad = [d for d, ok in rows if not ok and d > 0]
    rval = min(bad) if bad else None
    # cumulative safe mass over the ball
    sig = {}
    for r in range(rmax + 1):
        b = [ok for d, ok in rows if d <= r]
        sig[r] = sum(b) / len(b)
    ralpha = {}
    for a in ALPHAS:
        S = [r for r in range(rmax + 1) if sig[r] >= 1 - a - 1e-12]
        ralpha[a] = max(S) if S else -1
    # is the qualifying set a DOWN-SET?  (needed for "radius you can claim as safe")
    holes = {}
    for a in ALPHAS:
        S = set(r for r in range(rmax + 1) if sig[r] >= 1 - a - 1e-12)
        ra = ralpha[a]
        holes[a] = sorted(r for r in range(ra + 1) if r not in S) if ra >= 0 else []
    print(f"\n--- scenario {lab}: |Z|={len(z) if z is not None else None} Z={sorted(z) if z else z}")
    print(f"    reachable={len(rows)}  ecc(G0)=r_max={rmax}  r_val={rval}")
    print("    sigma(r) = " + "  ".join(f"{r}:{sig[r]:.3f}" for r in range(rmax + 1)))
    print("    r_alpha  = " + "  ".join(f"a={a}:{ralpha[a]}" for a in ALPHAS))
    print("    r_0 == r_val-1 ? " + (f"{ralpha[0.0]} vs {rval-1 if rval else 'NA'} -> "
          f"{'YES' if rval and ralpha[0.0]==rval-1 else 'NO'}"))
    for a in ALPHAS:
        if holes[a]:
            print(f"    NOT A DOWN-SET at alpha={a}: r_alpha={ralpha[a]} but sigma<1-a at r={holes[a]}")
    print(f"    shell sizes: " + "  ".join(
        f"{r}:{sum(1 for d,_ in rows if d==r)}" for r in range(rmax + 1)))
