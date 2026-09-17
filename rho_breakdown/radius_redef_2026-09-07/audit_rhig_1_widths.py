import sys, itertools, time
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.space import (enumerate_space, represented_dags, covering_pairs,
                                 neighbour_graph, bfs_distances)
from bkrobust.demo.evaluate import (optimal_adjustment_set_dag, optimal_adjustment_set_mpdag,
                                    random_sem, adjusted_estimand)
from harness import random_dag, make_problem

def tau_range(g, sem, x, y, cache):
    key = g.edge_string()
    if key in cache: return cache[key]
    vals = []
    for D in enumerate_dag_extensions(g):
        try:
            z = optimal_adjustment_set_dag(D, x, y)
            vals.append(adjusted_estimand(sem, x, y, z))
        except Exception:
            pass
    out = (min(vals), max(vals)) if vals else None
    cache[key] = out
    return out

def rhig(cpdag, claims, j, sem, x, y, cache):
    m = len(claims); lo, hi = np.inf, -np.inf; ok = False
    for S in itertools.combinations(range(m), m - j):
        g = apply_orientations(cpdag, [claims[i] for i in S])
        if g is None: continue
        r = tau_range(g, sem, x, y, cache)
        if r is None: continue
        ok = True; lo = min(lo, r[0]); hi = max(hi, r[1])
    return (lo, hi) if ok else None

rng = np.random.default_rng(20260907)
rows = []
t0 = time.time(); got = 0; tried = 0
while got < 200 and tried < 20000 and time.time() - t0 < 240:
    tried += 1
    pr = make_problem(rng, n=6, p=0.35, k_claims=3)
    if pr is None: continue
    dag, cpdag, x, y, K, g0 = pr["dag"], pr["cpdag"], pr["x"], pr["y"], pr["K"], pr["g0"]
    if len(K) < 2: continue
    sem = random_sem(dag, rng)
    cache = {}
    m = len(K)
    try:
        widths = []
        for j in range(0, m + 1):
            r = rhig(cpdag, K, j, sem, x, y, cache)
            widths.append(None if r is None else r[1] - r[0])
        blanket = tau_range(cpdag, sem, x, y, cache)
        bw = blanket[1] - blanket[0]
        # granularity attack: restate K as the full closure orientation set
        und = set(tuple(sorted(e)) for e in cpdag.undirected_edges)
        Kclos = [e for e in g0.directed_edges if tuple(sorted(e)) in und]
        r1c = rhig(cpdag, Kclos, 1, sem, x, y, cache)
        w1c = None if r1c is None else r1c[1] - r1c[0]
        # true tau
        tau_true = adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(dag, x, y))
    except Exception as e:
        continue
    got += 1
    rows.append(dict(m=m, mclos=len(Kclos), k=len(pr["undirected"]), w=widths,
                     blanket=bw, w1_closure=w1c, tau=tau_true))

print(f"problems={got} tried={tried} time={time.time()-t0:.1f}s")
W1 = np.array([r["w"][1] for r in rows if r["w"][1] is not None])
B  = np.array([r["blanket"] for r in rows])
W1C= np.array([r["w1_closure"] for r in rows if r["w1_closure"] is not None])
print(f"mean width RHIG_1 = {W1.mean():.4f}   blanket = {B.mean():.4f}   sum-ratio = {W1.sum()/B.sum():.3f}")
print(f"zero-width at blanket: {100*np.mean(B<1e-9):.1f}%   zero-width RHIG_1: {100*np.mean(W1<1e-9):.1f}%")
nz = B > 1e-9
print(f"CONDITIONAL on blanket>0 (n={nz.sum()}): mean RHIG_1={W1[nz].mean():.4f} blanket={B[nz].mean():.4f} sum-ratio={W1[nz].sum()/B[nz].sum():.3f}")
print()
print("=== GRANULARITY / GAMING ATTACK: same G0, K restated as its own Meek closure ===")
print(f"mean m (as elicited) = {np.mean([r['m'] for r in rows]):.2f}   mean m (closure) = {np.mean([r['mclos'] for r in rows]):.2f}")
print(f"mean RHIG_1 width, K as elicited = {W1.mean():.4f}")
print(f"mean RHIG_1 width, K as closure  = {W1C.mean():.4f}")
shrunk = np.array([ (r['w'][1] or 0) > 1e-9 and (r['w1_closure'] or 0) < 1e-9 for r in rows])
print(f"problems where restating K as its closure sends RHIG_1 width to EXACTLY 0: {shrunk.sum()}/{len(rows)} = {100*shrunk.mean():.1f}%")
strict = np.array([ (r['w1_closure'] is not None and r['w'][1] is not None and r['w1_closure'] < r['w'][1]-1e-9) for r in rows])
print(f"problems where the restatement strictly NARROWS RHIG_1: {strict.sum()}/{len(rows)} = {100*strict.mean():.1f}%")
print()
print("=== width vs m (is j comparable across problems?) ===")
for mm in sorted(set(r['m'] for r in rows)):
    sub = [r for r in rows if r['m']==mm]
    w1 = np.array([r['w'][1] for r in sub]); b=np.array([r['blanket'] for r in sub])
    print(f"  m={mm}  n={len(sub)}  mean RHIG_1={w1.mean():.4f}  blanket={b.mean():.4f}  ratio={w1.sum()/max(b.sum(),1e-12):.3f}")
print()
print("=== saturation: at which j does RHIG_j equal the blanket? ===")
sat = Counter()
for r in rows:
    for j,w in enumerate(r['w']):
        if w is not None and abs(w - r['blanket']) < 1e-9:
            sat[j]+=1; break
    else: sat['never']+=1
print(dict(sat))
