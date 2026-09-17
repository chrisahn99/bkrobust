#!/usr/bin/env python3
"""witness-8 FINAL: the b-bounded load-bearing cover U_b, un-minimised,
hop-restricted at b=2.  Measures: monotonicity of lambda, the b=1 zero-width
degeneracy, the b=2 non-degeneracy, cost, and the price of the hop restriction."""
import sys, itertools, collections
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import make_problem
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import optimal_adjustment_set_dag, adjusted_estimand, random_sem
TOL = 1e-9

def theta(cpdag, orient, sem, x, y, cache, ctr):
    """Theta(G[sigma]) = { value the analyst reports believing D }, D over extensions."""
    key = frozenset(orient)
    if key in cache: return cache[key]
    ctr[0] += 1
    g = apply_orientations(cpdag, list(orient))
    cache[key] = frozenset() if g is None else frozenset(np.round(
        [adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(d, x, y))
         for d in enumerate_dag_extensions(g)], 12))
    return cache[key]

def profiles(K, A):
    """All revision profiles sigma on A: each k in A either kept... no: retracted or reversed.
    Retraction = drop the claim; reversal = flip it.  2^|A| profiles."""
    rest = [k for k in K if k not in A]; A = list(A)
    for bits in itertools.product((0, 1), repeat=len(A)):
        yield rest + [(a[1], a[0]) for a, t in zip(A, bits) if t == 1]

def hull(cpdag, K, A, sem, x, y, cache, ctr):
    v = set(theta(cpdag, K, sem, x, y, cache, ctr))
    for o in profiles(K, tuple(A)): v |= theta(cpdag, o, sem, x, y, cache, ctr)
    return v

def skeleton_hops(cpdag, x, y):
    adj = collections.defaultdict(set)
    for a, b in cpdag.directed_edges: adj[a].add(b); adj[b].add(a)
    for e in cpdag.undirected_edges:
        a, b = tuple(e); adj[a].add(b); adj[b].add(a)
    d = {x: 0, y: 0}; q = collections.deque([x, y])
    while q:
        u = q.popleft()
        for w in adj[u]:
            if w not in d: d[w] = d[u] + 1; q.append(w)
    return d

rng = np.random.default_rng(90712)
rows = []; got = tried = 0
while got < 300 and tried < 300000:
    tried += 1
    pr = make_problem(rng, n=7, p=0.35, k_claims=5)
    if pr is None: continue
    K0 = [tuple(e) for e in pr["K"]]
    if len(K0) < 4: continue
    cpdag, dag, x, y = pr["cpdag"], pr["dag"], pr["x"], pr["y"]
    sem = random_sem(dag, rng); truth = sem.true_total_effect(x, y)
    if abs(truth) < 1e-9: continue                      # FIX 5: null queries out
    bad = [K0[i] for i in rng.permutation(len(K0))[:2]] # two false claims
    K = [((k[1], k[0]) if k in bad else k) for k in K0]
    if apply_orientations(cpdag, K) is None: continue
    cache = {}; ctr = [0]
    base = theta(cpdag, K, sem, x, y, cache, ctr)
    if not base: continue
    got += 1
    w0 = max(base) - min(base)
    hops = skeleton_hops(cpdag, x, y)
    lk = {k: min(hops.get(k[0], 99), hops.get(k[1], 99)) for k in K}

    lam1 = {}
    for k in K:
        v = hull(cpdag, K, (k,), sem, x, y, cache, ctr); lam1[k] = (max(v) - min(v)) - w0
    U1 = {k for k in K if lam1[k] > TOL}
    pairs_all = list(itertools.combinations(K, 2))
    pairs_hop = [(a, b) for a, b in pairs_all if min(lk[a], lk[b]) <= 1]
    def cover(pairs):
        U = set(U1)
        for a, b in pairs:
            if a in U and b in U: continue
            v = hull(cpdag, K, (a, b), sem, x, y, cache, ctr)
            if (max(v) - min(v)) - w0 > TOL: U |= {a, b}
        return U
    Ufull = cover(pairs_all); n_full = ctr[0]
    ctr2 = [0]; cache2 = dict(cache)
    Uhop = cover(pairs_hop)
    def rep(S):
        v = hull(cpdag, K, tuple(sorted(S)), sem, x, y, cache, ctr)
        return (max(v) - min(v), min(v) - 1e-9 <= truth <= max(v) + 1e-9)
    w1, c1 = rep(U1); wf, cf = rep(Ufull); wh, ch = rep(Uhop); wK, cK = rep(set(K))
    rows.append((w1, c1, wf, cf, wh, ch, wK, cK,
                 len(U1), len(Ufull), len(Uhop), len(K),
                 len(pairs_all), len(pairs_hop), w0))
R = np.array(rows, float)
I = R[:, 6] > TOL
print(f"non-null problems {got} / {tried} draws;  informative stratum {int(I.sum())}/{got}"
      f"  (blanket width > 0)")
print(f"|K| mean {R[I,11].mean():.2f}   |U_1| {R[I,8].mean():.2f}"
      f"   |U_2 full| {R[I,9].mean():.2f}   |U_2 hop<=1| {R[I,10].mean():.2f}")
print(f"pairs enumerated: all {R[I,12].mean():.2f}   hop-restricted {R[I,13].mean():.2f}"
      f"   ({100*R[I,13].sum()/R[I,12].sum():.1f}% of them)")
print(f"P(U_1 empty | informative) = {np.mean(R[I,8]==0):.3f}"
      f"    P(U_2 full empty | informative) = {np.mean(R[I,9]==0):.3f}")
print(f"P(Theta(G0) is a single point) = {np.mean(R[:,14] <= TOL):.3f}   (all non-null draws)")
for nm, wi, ci in (("b=1  (U_1)", 0, 1), ("b=2  full pairs", 2, 3),
                   ("b=2  hop<=1 pairs", 4, 5), ("blanket (whole class)", 6, 7)):
    print(f"  {nm:24s} coverage {R[I,ci].mean():.3f}   width {R[I,wi].mean():.4f}"
          f"   width/blanket {np.mean(R[I,wi]/R[I,6]):.3f}")
print(f"agreement U_2hop == U_2full : {np.mean([a==b for a,b in zip(R[I,10],R[I,9])]):.3f}")
