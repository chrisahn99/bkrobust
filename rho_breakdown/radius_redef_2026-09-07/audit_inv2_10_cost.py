#!/usr/bin/env python3
"""inverse-2 audit, part 10: the 'scales past exhaustive enumeration' claim.
phi_1 skips enumerate_space, but is_valid_adjustment_set_mpdag calls
enumerate_dag_extensions on graphs that are COARSER than G0 (retraction adds
undirected edges), so the exponential is still there."""
import sys, time
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import random_dag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo import meek as M
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.space import enumerate_space
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag, is_valid_adjustment_set_mpdag
import bkrobust.demo.evaluate as E

COUNT = {"ext": 0, "dags": 0}
_orig = M.enumerate_dag_extensions
def counting(g, limit=None):
    out = _orig(g, limit)
    COUNT["ext"] += 1; COUNT["dags"] += len(out)
    return out
M.enumerate_dag_extensions = counting

rng = np.random.default_rng(31337)
rows = []
tried = 0
while len(rows) < 40 and tried < 20000:
    tried += 1
    dag = random_dag(rng, 7, 0.32)
    if not dag.directed_edges: continue
    cpdag = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cpdag.undirected_edges)
    if not (4 <= len(und) <= 6): continue
    nodes = sorted(dag.nodes)
    cand = [(a,b) for a in nodes for b in nodes if a != b and b in dag.descendants(a)]
    rng.shuffle(cand)
    got = None
    for x, y in cand:
        K = []
        for t in rng.permutation(len(und))[:4]:
            a, b = und[t]
            K.append((a,b) if (a,b) in dag.directed_edges else (b,a))
        g0 = apply_orientations(cpdag, K)
        if g0 is None: continue
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None: continue
        got = (x, y, K, g0, frozenset(z)); break
    if got is None: continue
    x, y, K, g0, z = got
    COUNT["ext"] = COUNT["dags"] = 0
    t0 = time.time()
    for k in K:
        rest = [e for e in K if e != k]
        for extra in ([], [(k[1], k[0])]):
            g = apply_orientations(cpdag, rest + extra)
            if g is not None: is_valid_adjustment_set_mpdag(g, x, y, z)
    t_phi, d_phi = time.time()-t0, COUNT["dags"]
    COUNT["ext"] = COUNT["dags"] = 0
    t0 = time.time(); sp = enumerate_space(cpdag); t_sp, d_sp = time.time()-t0, COUNT["dags"]
    rows.append((len(und), len(sp), d_phi, d_sp, t_phi, t_sp))

a = np.array(rows, float)
print(f"n={len(a)}   mean |undirected(Chat)| = {a[:,0].mean():.1f}")
print(f"   |space| (MPDAGs enumerated by the space-based defs)  mean {a[:,1].mean():7.1f}")
print(f"   DAG extensions enumerated BY phi_1 alone             mean {a[:,2].mean():7.1f}"
      f"   max {a[:,2].max():.0f}")
print(f"   DAG extensions enumerated by enumerate_space          mean {a[:,3].mean():7.1f}")
print(f"   wall time  phi_1 {a[:,4].mean()*1000:6.1f} ms    enumerate_space "
      f"{a[:,5].mean()*1000:7.1f} ms   ratio {a[:,5].mean()/a[:,4].mean():.1f}x")
print("\n   phi_1 is cheaper by a constant/low-polynomial factor, but it still calls")
print("   enumerate_dag_extensions, so the exponential in the undirected component size")
print("   is NOT removed -- only the 3^k MPDAG loop around it is.")
