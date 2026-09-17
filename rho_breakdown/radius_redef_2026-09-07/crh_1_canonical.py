#!/usr/bin/env python3
"""CRH_j: canonical retraction hull over inclusion-minimal generating sets.
Measures: rung count m*, saturation j_sat, width vs blanket, phrasing invariance,
containment vs the phrasing-dependent RHIG_j, and the (d* x j) coverage grid."""
import sys, itertools, time
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from bkrobust.demo.graph import canon
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import optimal_adjustment_set_dag, random_sem, adjusted_estimand
from harness import random_dag

# ---------- core ----------
def tau_range(g, sem, x, y, cache):
    key = (g.edge_string(), x, y)
    if key in cache: return cache[key]
    vals = []
    for D in enumerate_dag_extensions(g):
        try: vals.append(adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(D, x, y)))
        except Exception: pass
    out = (min(vals), max(vals)) if vals else None
    cache[key] = out; return out

def oriented_beyond(cpdag, g0):
    und = {canon(a,b) for a,b in cpdag.undirected_edges}
    return sorted((a,b) for (a,b) in g0.directed_edges if canon(a,b) in und)

def minimal_generating_sets(cpdag, g0, Or):
    gens = []
    for r in range(0, len(Or)+1):
        for S in itertools.combinations(Or, r):
            if any(set(G) <= set(S) for G in gens):   # not inclusion-minimal
                continue
            h = apply_orientations(cpdag, S)
            if h is not None and h == g0:
                gens.append(S)
    return gens

def crh(cpdag, gens, j, sem, x, y, cache):
    lo, hi, ok = np.inf, -np.inf, False
    for S in gens:
        t = max(len(S)-j, 0)
        for T in itertools.combinations(S, t):
            g = apply_orientations(cpdag, T)
            if g is None: continue
            r = tau_range(g, sem, x, y, cache)
            if r is None: continue
            ok = True; lo = min(lo, r[0]); hi = max(hi, r[1])
    return (lo, hi) if ok else None

def rhig(cpdag, K, j, sem, x, y, cache):        # phrasing-dependent original
    m = len(K); lo, hi, ok = np.inf, -np.inf, False
    for S in itertools.combinations(K, max(m-j,0)):
        g = apply_orientations(cpdag, S)
        if g is None: continue
        r = tau_range(g, sem, x, y, cache)
        if r is None: continue
        ok = True; lo=min(lo,r[0]); hi=max(hi,r[1])
    return (lo,hi) if ok else None

def dstar(gens, dag):
    """min over minimal generating sets of the number of facts reversed in the truth."""
    best = None
    for S in gens:
        f = sum(1 for (a,b) in S if (b,a) in dag.directed_edges)
        best = f if best is None else min(best, f)
    return best

# ---------- sweep ----------
rng = np.random.default_rng(11)  # cache key includes (x,y): a per-DAG cache keyed on the graph alone is WRONG
rows = []; t0 = time.time(); tried = 0
TARGET = 150
while len(rows) < TARGET and tried < 200000 and time.time()-t0 < 420:
    tried += 1
    dag = random_dag(rng, 7, 0.32)
    if not dag.directed_edges: continue
    cpdag = dag_to_cpdag(dag)
    und = sorted(canon(*e) for e in cpdag.undirected_edges)
    if not (2 <= len(und) <= 6): continue
    nodes = sorted(dag.nodes)
    cand = [(x,y) for x in nodes for y in nodes if x!=y]
    rng.shuffle(cand)
    sem = random_sem(dag, rng); cache = {}
    for x,y in cand[:6]:
        blank = tau_range(cpdag, sem, x, y, cache)
        if blank is None or blank[1]-blank[0] < 1e-8: continue     # NON-TRIVIAL only
        km = int(rng.integers(1, min(4, len(und))+1))
        idx = rng.permutation(len(und))[:km]
        nfalse = int(rng.integers(0, 3))                            # 0,1,2 reversed
        flip = set(rng.permutation(km)[:min(nfalse, km)].tolist())
        K = []
        for pos, t in enumerate(idx):
            a,b = und[t]
            true_dir = (a,b) if (a,b) in dag.directed_edges else (b,a)
            K.append(true_dir[::-1] if pos in flip else true_dir)
        g0 = apply_orientations(cpdag, K)
        if g0 is None: continue
        try: tau_true = adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(dag, x, y))
        except Exception: continue
        Or = oriented_beyond(cpdag, g0)
        gens = minimal_generating_sets(cpdag, g0, Or)
        if not gens: continue
        mstar = min(len(S) for S in gens); Mmax = max(len(S) for S in gens)
        ds = dstar(gens, dag)
        W, COV = [], []
        for j in range(0, mstar+1):
            r = crh(cpdag, gens, j, sem, x, y, cache)
            W.append(None if r is None else r[1]-r[0])
            COV.append(None if r is None else (r[0]<=tau_true+1e-9 <= 0 or (r[0]-1e-9<=tau_true<=r[1]+1e-9)))
        # phrasing variants of the SAME G0
        Kclos = Or                                     # K restated as its own closure
        r1_K   = rhig(cpdag, K, 1, sem, x, y, cache)
        r1_cl  = rhig(cpdag, list(Kclos), 1, sem, x, y, cache)
        crh1   = crh(cpdag, gens, 1, sem, x, y, cache)
        rows.append(dict(k=len(und), m=len(K), q=len(Or), mstar=mstar, Mmax=Mmax,
                         ngens=len(gens), dstar=ds, nfalse=nfalse,
                         blanket=blank[1]-blank[0], W=W, COV=COV,
                         w1_K=None if r1_K is None else r1_K[1]-r1_K[0],
                         w1_cl=None if r1_cl is None else r1_cl[1]-r1_cl[0],
                         w1_crh=None if crh1 is None else crh1[1]-crh1[0]))
        break

n=len(rows); print(f"n={n} non-trivial problems (blanket width>0), tried={tried}, {time.time()-t0:.1f}s\n")
B  = np.array([r['blanket'] for r in rows])
W1 = np.array([r['W'][1] if len(r['W'])>1 else r['W'][0] for r in rows])   # CRH_1 (=blanket if m*=1... W has m*+1 entries)
W1 = np.array([r['W'][min(1,len(r['W'])-1)] for r in rows])
print("--- ladder length m* (canonical rung count) ---")
print(" m* :", dict(sorted(Counter(r['mstar'] for r in rows).items())))
print(" |K|:", dict(sorted(Counter(r['m'] for r in rows).items())))
print(f" P(m*=1) = {np.mean([r['mstar']==1 for r in rows]):.3f}   median m* = {np.median([r['mstar'] for r in rows]):.1f}")
print(f" mean #minimal generating sets = {np.mean([r['ngens'] for r in rows]):.2f}  |Or|={np.mean([r['q'] for r in rows]):.2f}")

print("\n--- CRH_1 width vs blanket, on the non-trivial stratum ---")
print(f" mean CRH_1 = {W1.mean():.4f}   mean blanket = {B.mean():.4f}   sum-ratio = {W1.sum()/B.sum():.3f}")
print(f" P(CRH_1 == blanket) = {np.mean(np.abs(W1-B)<1e-9):.3f}")
print(f" P(CRH_1 width == 0) = {np.mean(W1<1e-9):.3f}")

sat = Counter()
for r in rows:
    for j,w in enumerate(r['W']):
        if w is not None and abs(w-r['blanket'])<1e-9: sat[j]+=1; break
    else: sat['never']+=1
print(" j_sat:", dict(sorted(sat.items(), key=lambda t: str(t[0]))))
interior = np.mean([sum(1 for a,b in zip(r['W'], r['W'][1:]) if b-a>1e-9) >= 2 for r in rows])
print(f" P(frontier strictly increases at >=2 rungs) = {interior:.3f}")

print("\n--- phrasing invariance (same G0, K vs its Meek closure) ---")
ok = np.array([r['w1_K'] is not None and r['w1_cl'] is not None for r in rows])
a = np.array([r['w1_K'] if r['w1_K'] is not None else np.nan for r in rows])
b = np.array([r['w1_cl'] if r['w1_cl'] is not None else np.nan for r in rows])
c = np.array([r['w1_crh'] if r['w1_crh'] is not None else np.nan for r in rows])
print(f" RHIG_1 changes under restatement: {np.mean(np.abs(a[ok]-b[ok])>1e-9):.3f}  (mean {np.nanmean(a[ok]):.4f} -> {np.nanmean(b[ok]):.4f})")
print(f" CRH_1 >= RHIG_1(as elicited): {np.mean(c[ok] >= a[ok]-1e-9):.3f}   strictly wider: {np.mean(c[ok] > a[ok]+1e-9):.3f}")
print(f" CRH_1 >= RHIG_1(as closure) : {np.mean(c[ok] >= b[ok]-1e-9):.3f}")

print("\n--- coverage grid: canonical defect d* x assumed budget j ---")
grid = {}
for r in rows:
    ds = r['dstar']
    for j,w in enumerate(r['W']):
        cov = r['COV'][j]
        grid.setdefault((ds,j), []).append(bool(cov))
for (ds,j) in sorted(grid):
    v = grid[(ds,j)]
    tag = "j>=d* (theorem)" if j>=ds else "j<d*  (informative)"
    print(f"  d*={ds} j={j}: coverage {np.mean(v):.3f}  (n={len(v):3d})  {tag}")
print("\n d* given #reversed claims:", {f: dict(sorted(Counter(r['dstar'] for r in rows if r['nfalse']==f).items())) for f in (0,1,2)})
