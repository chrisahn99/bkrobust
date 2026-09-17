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
from bkrobust.demo.evaluate import (optimal_adjustment_set_dag, random_sem, adjusted_estimand)
from harness import random_dag

def tau_range(g, sem, x, y, cache):
    key = g.edge_string()
    if key in cache: return cache[key]
    vals = []
    for D in enumerate_dag_extensions(g):
        try: vals.append(adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(D, x, y)))
        except Exception: pass
    out = (min(vals), max(vals)) if vals else None
    cache[key] = out; return out

def rhig(cpdag, claims, j, sem, x, y, cache):
    m = len(claims); lo, hi = np.inf, -np.inf; ok=False
    for S in itertools.combinations(range(m), m-j):
        g = apply_orientations(cpdag, [claims[i] for i in S])
        if g is None: continue
        r = tau_range(g, sem, x, y, cache)
        if r is None: continue
        ok=True; lo=min(lo,r[0]); hi=max(hi,r[1])
    return (lo,hi) if ok else None

rng = np.random.default_rng(7)
rows=[]; t0=time.time(); tried=0
while len(rows) < 120 and tried < 60000 and time.time()-t0 < 300:
    tried += 1
    dag = random_dag(rng, 7, 0.32)
    if not dag.directed_edges: continue
    cpdag = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cpdag.undirected_edges)
    if not (2 <= len(und) <= 6): continue
    nodes = sorted(dag.nodes)
    cand = [(x,y) for x in nodes for y in nodes if x!=y]
    rng.shuffle(cand)
    sem = random_sem(dag, rng); cache={}
    for x,y in cand[:6]:
        blank = tau_range(cpdag, sem, x, y, cache)
        if blank is None or blank[1]-blank[0] < 1e-8:   # keep only NON-TRIVIAL problems
            continue
        km = int(rng.integers(1, min(4, len(und))+1))
        idx = rng.permutation(len(und))[:km]
        K=[]
        for t in idx:
            a,b = und[t]; K.append((a,b) if (a,b) in dag.directed_edges else (b,a))
        g0 = apply_orientations(cpdag, K)
        if g0 is None: continue
        try: tau_true = adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(dag, x, y))
        except Exception: continue
        m=len(K)
        w = [None]*(m+1)
        for j in range(m+1):
            r = rhig(cpdag, K, j, sem, x, y, cache)
            w[j] = None if r is None else (r[1]-r[0], r[0]<=tau_true+1e-9 and tau_true<=r[1]+1e-9)
        Kclos=[e for e in g0.directed_edges if tuple(sorted(e)) in set(und)]
        r1c = rhig(cpdag, Kclos, 1, sem, x, y, cache)
        # hop-1 ball hull (BIG_1) in the lattice
        space = enumerate_space(cpdag); reps = represented_dags(space)
        covers = covering_pairs(space, reps); nbrs = neighbour_graph(space, covers)
        g0s = next((g for g in space if g == g0), None)
        big1 = None
        if g0s is not None:
            dist = bfs_distances(nbrs, g0s)
            lo,hi=np.inf,-np.inf
            for g,d in dist.items():
                if d is not None and d<=1:
                    r = tau_range(g, sem, x, y, cache)
                    if r: lo=min(lo,r[0]); hi=max(hi,r[1])
            big1 = hi-lo if lo<np.inf else None
        rows.append(dict(m=m, mclos=len(Kclos), k=len(und), w=w, blanket=blank[1]-blank[0],
                         w1c=None if r1c is None else r1c[1]-r1c[0], big1=big1,
                         nspace=len(space)))
        break

print(f"NON-TRIVIAL problems only (blanket width > 0): n={len(rows)}  tried={tried}  {time.time()-t0:.1f}s")
B=np.array([r['blanket'] for r in rows]); W1=np.array([r['w'][1][0] for r in rows])
print(f"mean RHIG_1={W1.mean():.4f}  blanket={B.mean():.4f}  sum-ratio={W1.sum()/B.sum():.3f}")
print(f"RHIG_1 == blanket (saturates at j=1): {100*np.mean(np.abs(W1-B)<1e-9):.1f}%")
print(f"RHIG_1 width == 0 even though blanket>0: {100*np.mean(W1<1e-9):.1f}%")
cov1 = np.mean([r['w'][1][1] for r in rows]); cov0 = np.mean([r['w'][0][1] for r in rows])
print(f"coverage of tau_true (K is TRUTHFUL here): RHIG_0={cov0:.3f} RHIG_1={cov1:.3f}")
BG=np.array([r['big1'] if r['big1'] is not None else np.nan for r in rows])
ok=~np.isnan(BG)
print(f"\n=== RHIG_1 vs BIG_1 (hop-1 hull) on {ok.sum()} problems ===")
print(f"mean BIG_1={BG[ok].mean():.4f}  mean RHIG_1={W1[ok].mean():.4f}")
print(f"identical width: {100*np.mean(np.abs(BG[ok]-W1[ok])<1e-9):.1f}% | RHIG_1 wider: {100*np.mean(W1[ok]>BG[ok]+1e-9):.1f}% | BIG_1 wider: {100*np.mean(BG[ok]>W1[ok]+1e-9):.1f}%")
print(f"\n=== granularity attack (restate K as its Meek closure) ===")
casc=np.array([r['mclos']-r['m'] for r in rows])
print(f"cascade size (closure - elicited): mean {casc.mean():.2f} max {casc.max()}")
W1C=np.array([r['w1c'] if r['w1c'] is not None else np.nan for r in rows])
g=~np.isnan(W1C)
print(f"mean RHIG_1 as elicited={W1[g].mean():.4f}  as closure={W1C[g].mean():.4f}")
print(f"strictly narrower after restatement: {100*np.mean(W1C[g]<W1[g]-1e-9):.1f}%   driven to exactly 0: {100*np.mean((W1C[g]<1e-9)&(W1[g]>1e-9)):.1f}%")
sub = casc>0
print(f"  restricted to problems with a real cascade (n={sub.sum()}): narrower {100*np.mean(W1C[g&sub]<W1[g&sub]-1e-9):.1f}%  mean {W1[g&sub].mean():.4f} -> {W1C[g&sub].mean():.4f}")
print(f"\n=== saturation index ===")
sat=Counter()
for r in rows:
    for j,ww in enumerate(r['w']):
        if ww is not None and abs(ww[0]-r['blanket'])<1e-9: sat[j]+=1; break
    else: sat['never']+=1
print(dict(sat))
print(f"\nmean |space| (3^k enumeration) = {np.mean([r['nspace'] for r in rows]):.1f}; mean 2^k = {np.mean([2**r['k'] for r in rows]):.1f}; mean m = {np.mean([r['m'] for r in rows]):.2f}")
