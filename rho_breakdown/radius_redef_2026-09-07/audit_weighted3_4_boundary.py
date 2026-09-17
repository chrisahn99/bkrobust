import sys
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from collections import Counter
from harness import make_problem
from bkrobust.demo.space import (enumerate_space, represented_dags, covering_pairs,
                                 neighbour_graph, bfs_distances)
from bkrobust.demo.evaluate import is_valid_adjustment_set_dag

rng = np.random.default_rng(4242)
got = tried = 0
strict_vs_not = Counter(); growth_agree = 0; growth_n = 0; both_inf = 0
while got < 200 and tried < 30000:
    tried += 1
    pr = make_problem(rng)
    if pr is None: continue
    space = enumerate_space(pr["cpdag"])
    if len(space) > 300: continue
    reps = represented_dags(space)
    g0 = next((g for g in space if g == pr["g0"]), None)
    if g0 is None: continue
    x,y,z = pr["x"], pr["y"], pr["z"]
    if any(not is_valid_adjustment_set_dag(d,x,y,z) for d in reps[g0]): continue
    got += 1
    covers = covering_pairs(space, reps); nbrs = neighbour_graph(space, covers)
    dist = bfs_distances(nbrs, g0); n0 = len(reps[g0])
    info = []
    for g, d in dist.items():
        if d == 0: continue
        ext = reps[g]
        b = sum(1 for dd in ext if not is_valid_adjustment_set_dag(dd,x,y,z))/len(ext)
        info.append((d, b, len(ext)))
    rs = min([d for d,b,_ in info if b > 0.5], default=None)
    rn = min([d for d,b,_ in info if b >= 0.5], default=None)
    strict_vs_not[(rs, rn)] += 1
    dg = min([d for d,_,n in info if n > 2*n0], default=None)
    if rs is None and dg is None: both_inf += 1
    if rs is not None:
        growth_n += 1; growth_agree += (dg == rs)
print(f"well-posed problems: {got}")
print("(r^0.5 with STRICT beta>alpha, r^0.5 with beta>=alpha):")
for k,v in sorted(strict_vs_not.items(), key=lambda kv: -kv[1]):
    print(f"   strict={k[0]!s:>4}  non-strict={k[1]!s:>4}   {v:3d}  ({100*v/got:.1f}%)")
flip = sum(v for (a,b),v in strict_vs_not.items() if a != b)
print(f"   the strict/non-strict choice changes r^0.5 on {flip}/{got} = {flip/got:.3f} of problems")
print(f"\nr^0.5 equals the pure class-growth distance min(d : |[G]|>2|[G0]|) on {growth_agree}/{growth_n} finite cases"
      + (f" = {growth_agree/growth_n:.3f}" if growth_n else ""))
print(f"both infinite together: {both_inf}/{got}")
