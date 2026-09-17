"""Is beta(G) <= 1/2 at hop 1 a law for well-posed problems? Sweep configs."""
import sys
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import make_problem
from bkrobust.demo.space import (enumerate_space, represented_dags, covering_pairs,
                                 neighbour_graph, bfs_distances)
from bkrobust.demo.evaluate import is_valid_adjustment_set_dag

def run(n, p, k, target=120, seed=7):
    rng = np.random.default_rng(seed)
    got = tried = 0
    max1 = -1.0; over_half_1 = 0; hop1_pos = 0
    max2 = -1.0; over_half_2 = 0
    up_max = -1.0
    while got < target and tried < 30000:
        tried += 1
        pr = make_problem(rng, n=n, p=p, k_claims=k)
        if pr is None: continue
        space = enumerate_space(pr["cpdag"])
        if len(space) > 400: continue
        reps = represented_dags(space)
        g0 = next((g for g in space if g == pr["g0"]), None)
        if g0 is None: continue
        x,y,z = pr["x"], pr["y"], pr["z"]
        ext0 = reps[g0]
        if any(not is_valid_adjustment_set_dag(d,x,y,z) for d in ext0):
            continue                      # ill-posed: Z already broken at G0
        got += 1
        covers = covering_pairs(space, reps)
        nbrs = neighbour_graph(space, covers)
        dist = bfs_distances(nbrs, g0)
        for g, d in dist.items():
            if d == 0: continue
            ext = reps[g]
            b = sum(1 for dd in ext if not is_valid_adjustment_set_dag(dd,x,y,z))/len(ext)
            if d == 1:
                max1 = max(max1, b); hop1_pos += (b>0); over_half_1 += (b > 0.5+1e-12)
            elif d == 2:
                max2 = max(max2, b); over_half_2 += (b > 0.5+1e-12)
    return got, max1, over_half_1, hop1_pos, max2, over_half_2

print(f"{'n':>2} {'p':>5} {'k':>2} {'probs':>6} {'max beta hop1':>14} {'#hop1>0.5':>10} {'#hop1>0':>8} {'max beta hop2':>14} {'#hop2>0.5':>10}")
for n,p,k in [(5,0.4,2),(6,0.30,2),(6,0.35,3),(6,0.40,4),(7,0.30,3),(7,0.35,4),(8,0.30,4)]:
    got,m1,o1,h1,m2,o2 = run(n,p,k)
    print(f"{n:>2} {p:>5.2f} {k:>2} {got:>6} {m1:>14.4f} {o1:>10} {h1:>8} {m2:>14.4f} {o2:>10}")
