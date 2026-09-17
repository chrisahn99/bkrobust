"""SEL_A audit 4 -- (a) is the screen's firing event the same event as r_val = 1,
(b) the redundancy ladder: how many copies of a claim defeat an order-r screen.
"""
import sys, itertools, time
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.evaluate import is_valid_adjustment_set_mpdag
from bkrobust.demo.space import (enumerate_space, represented_dags, covering_pairs,
                                 neighbour_graph, bfs_distances)
from audit_sel_2_coverage import build, hull_of, Z196

NOBS = 10000


def r_val_lattice(pr):
    """The original definition: smallest hop count at which Z is invalid."""
    space = enumerate_space(pr["cp"])
    reps = represented_dags(space)
    nbrs = neighbour_graph(space, covering_pairs(space, reps))
    g0 = next((g for g in space if g == pr["g0"]), None)
    if g0 is None:
        return None
    dist = bfs_distances(nbrs, g0)
    bad = [d for g, d in dist.items()
           if d is not None and d > 0
           and not is_valid_adjustment_set_mpdag(g, pr["x"], pr["y"], pr["z"])]
    return min(bad) if bad else None


def screen(pr, K, order, cache):
    """Active set of an order-`order` displacement screen on claim list K."""
    m = len(K)
    th0, delta = pr["theta0"], Z196*np.sqrt(pr["avar"]/NOBS)
    A = set()
    for r in range(1, order+1):
        for S in itertools.combinations(range(m), r):
            h = hull_of(pr["cp"], [K[i] for i in range(m) if i not in S],
                        pr["sem"], pr["x"], pr["y"], cache)
            if h is not None and max(abs(h[0]-th0), abs(h[1]-th0)) > delta:
                A |= set(S)
    return A


if __name__ == "__main__":
    rng = np.random.default_rng(31337)
    rows, tried, t0 = [], 0, time.time()
    while len(rows) < 250 and tried < 60000 and time.time()-t0 < 300:
        tried += 1
        pr = build(rng, 0)                      # truthful K: r_val is the paper's setting
        if pr is None:
            continue
        cache = {}
        try:
            A1 = screen(pr, pr["K"], 1, cache)
            rv = r_val_lattice(pr)
        except Exception:
            continue
        rows.append(dict(pr=pr, A1=A1, rv=rv, cache=cache))
    print(f"problems={len(rows)}  {time.time()-t0:.1f}s")

    fire = np.array([len(r["A1"]) > 0 for r in rows])
    rv1 = np.array([r["rv"] == 1 for r in rows])
    unreached = np.array([r["rv"] is None for r in rows])
    print(f"\n[a] r_val = 1 in {100*rv1.mean():.1f}% of problems  (the reported defect)")
    print(f"    the SEL_A screen fires (A1 non-empty) in {100*fire.mean():.1f}%")
    print(f"    agreement of the two indicators: {100*np.mean(fire == rv1):.1f}%  "
          f"(phi = {np.corrcoef(fire.astype(float), rv1.astype(float))[0,1]:.3f})")
    print(f"    fires & r_val=1: {100*np.mean(fire & rv1):.1f}% | "
          f"silent & r_val>1 or unreached: {100*np.mean(~fire & ~rv1):.1f}% | "
          f"fires & r_val>1: {100*np.mean(fire & ~rv1):.1f}% | "
          f"silent & r_val=1: {100*np.mean(~fire & rv1):.1f}%")
    print(f"    r_val UNREACHED (Z never breaks): {100*unreached.mean():.1f}%")

    # ---------------- (b) the redundancy ladder ----------------
    print("\n[b] REDUNDANCY LADDER: restate each claim as c logically equivalent copies")
    live = [r for r in rows if len(r["A1"]) > 0]
    print(f"    on the {len(live)} problems where the order-1 screen fires on the elicited K")
    for c in (1, 2, 3):
        res = {}
        for order in (1, 2, 3):
            if order >= 3 and c < 3:
                continue
            fires = 0
            for r in live:
                K = list(r["pr"]["K"]) * c
                if len(K) > 9:
                    continue
                A = screen(r["pr"], K, order, {})
                fires += len(A) > 0
            res[order] = 100*fires/max(len(live), 1)
        line = "  ".join(f"order-{o}: fires {v:5.1f}%" for o, v in res.items())
        print(f"      {c} cop{'y ' if c==1 else 'ies'}: {line}")
    print("    (a screen of order r is blind to a belief stated r+1 times; the cost of "
          "order r is C(m,r) closures on the inflated m)")
