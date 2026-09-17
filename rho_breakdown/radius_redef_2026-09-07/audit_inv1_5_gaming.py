"""inverse-1: (a) how coarse is the alpha knob really, (b) can the analyst inflate r_alpha
by reporting a different -- still valid -- adjustment set at G0?"""
import sys, time
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.space import (enumerate_space, represented_dags, covering_pairs,
                                 neighbour_graph, bfs_distances)
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, is_valid_adjustment_set_mpdag,
                                    all_valid_adjustment_sets_mpdag, asymptotic_variance, random_sem)
from harness import make_problem

def ralpha_for(valid_by_dist, rmax, a):
    S = []
    for r in range(rmax + 1):
        v = [ok for d, ok in valid_by_dist if d <= r]
        if np.mean(v) >= 1 - a - 1e-12: S.append(r)
    return max(S) if S else -1

rng = np.random.default_rng(20260907)
GRID = np.round(np.arange(0, 0.96, 0.01), 2)
n_steps = []; act = []; games = []; var_price = []
t0 = time.time(); done = 0
while done < 250 and time.time() - t0 < 420:
    pr = make_problem(rng)
    if pr is None: continue
    space = enumerate_space(pr["cpdag"]); reps = represented_dags(space)
    nbrs = neighbour_graph(space, covering_pairs(space, reps))
    g0 = next((g for g in space if g == pr["g0"]), None)
    if g0 is None: continue
    dist = bfs_distances(nbrs, g0)
    if len(dist) < 4: continue
    x, y = pr["x"], pr["y"]; rmax = max(dist.values())
    def prof(z):
        return [(dist[g], bool(is_valid_adjustment_set_mpdag(g, x, y, z))) for g in space if g in dist]
    zstar = pr["z"]
    p0 = prof(zstar)
    if np.mean([ok for d, ok in p0 if d == 0]) < 1: continue     # keep only sigma(0)=1 problems
    done += 1
    vals = [ralpha_for(p0, rmax, a) for a in GRID]
    n_steps.append(len(set(vals)))
    first_change = next((a for a, v in zip(GRID, vals) if v != vals[0]), None)
    act.append(first_change)
    # gaming: any OTHER valid adjustment set at G0 with a strictly larger r_alpha at alpha=0.10
    alts = all_valid_adjustment_sets_mpdag(g0, x, y, max_size=4) or []
    base = ralpha_for(p0, rmax, 0.10)
    best = base; bz = zstar
    for z in alts:
        z = frozenset(z)
        if z == zstar: continue
        r = ralpha_for(prof(z), rmax, 0.10)
        if r > best: best, bz = r, z
    games.append((base, best))
    if best > base:
        sem = random_sem(pr["dag"], rng)
        try:
            var_price.append(asymptotic_variance(sem, x, y, bz) / asymptotic_variance(sem, x, y, zstar))
        except Exception: pass

print(f"problems (sigma(0)=1 only) = {done}   ({time.time()-t0:.0f}s)")
print(f"distinct r_alpha values over the WHOLE alpha grid 0..0.95: "
      f"mean {np.mean(n_steps):.2f}, median {np.median(n_steps):.0f}, "
      f"=1 (knob does nothing at all) in {100*np.mean([s==1 for s in n_steps]):.1f}% of problems")
have = [a for a in act if a is not None]
print(f"alpha at which r_alpha first moves off r_0: median {np.median(have):.2f}, "
      f"p25 {np.percentile(have,25):.2f}, p75 {np.percentile(have,75):.2f}; never moves in "
      f"{100*np.mean([a is None for a in act]):.1f}%")
up = [1 for b, g in games if g > b]
print(f"\nGAMING at alpha=0.10: a different VALID adjustment set at G0 raises r_alpha in "
      f"{100*len(up)/len(games):.1f}% of problems; mean gain "
      f"{np.mean([g-b for b,g in games if g>b]) if up else 0:.2f} hops "
      f"(max {max([g-b for b,g in games], default=0)})")
if var_price:
    print(f"asymptotic-variance price of the gamed set: median {np.median(var_price):.2f}x, "
          f"<=1.0x (free) in {100*np.mean([v<=1.0000001 for v in var_price]):.1f}% of the gamed cases")
