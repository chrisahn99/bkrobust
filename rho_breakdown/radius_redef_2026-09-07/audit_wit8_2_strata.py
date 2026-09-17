#!/usr/bin/env python3
"""witness-8 audit 2: stratified widths, the shape of lambda, and where the
synergy pairs sit relative to the query (the greedy's locality restriction)."""
import sys, itertools
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np, networkx as nx
from harness import make_problem
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import optimal_adjustment_set_dag, adjusted_estimand, random_sem

TOL = 1e-9

def theta(cpdag, orients, sem, x, y, cache):
    key = frozenset(orients)
    if key in cache: return cache[key]
    g = apply_orientations(cpdag, list(orients))
    vals = frozenset() if g is None else frozenset(
        np.round([adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(d, x, y))
                  for d in enumerate_dag_extensions(g)], 12))
    cache[key] = vals
    return vals

def profiles(K, A):
    rest = [k for k in K if k not in A]; A = list(A)
    for bits in itertools.product((0, 1), repeat=len(A)):
        yield rest + [(a[1], a[0]) for a, b in zip(A, bits) if b == 1]

def hull(cpdag, K, S, sem, x, y, cache, base):
    vals = set(base)
    for o in profiles(K, tuple(S)): vals |= theta(cpdag, o, sem, x, y, cache)
    return (min(vals), max(vals)) if vals else None

def skel_dist(cpdag, x, y):
    G = nx.Graph()
    G.add_nodes_from(cpdag.nodes)
    for a, b in cpdag.directed_edges: G.add_edge(a, b)
    for e in cpdag.undirected_edges: G.add_edge(*tuple(e))
    dx = nx.single_source_shortest_path_length(G, x)
    dy = nx.single_source_shortest_path_length(G, y)
    return lambda k: min(min(dx.get(k[0], 9), dx.get(k[1], 9)),
                         min(dy.get(k[0], 9), dy.get(k[1], 9)))

rng = np.random.default_rng(20260907)
got = tried = 0
rows = []; lam_pos = []; ell_L1 = []; ell_syn = []; base_pt = 0
while got < 400 and tried < 80000:
    tried += 1
    pr = make_problem(rng, n=7, p=0.35, k_claims=5)
    if pr is None: continue
    K = [tuple(e) for e in pr["K"]]
    if len(K) < 4: continue
    cpdag, dag, x, y = pr["cpdag"], pr["dag"], pr["x"], pr["y"]
    sem = random_sem(dag, rng); cache = {}
    base = theta(cpdag, K, sem, x, y, cache)
    if not base: continue
    got += 1
    base_w = max(base) - min(base); base_pt += (base_w <= TOL)
    ell = skel_dist(cpdag, x, y)
    sing = {}
    for k in K:
        h = hull(cpdag, K, (k,), sem, x, y, cache, base)
        sing[k] = (h[1] - h[0]) - base_w
    L1 = [k for k in K if sing[k] > TOL]
    lam_pos += [sing[k] for k in L1]; ell_L1 += [ell(k) for k in L1]
    syn = []
    for a, b in itertools.combinations(K, 2):
        if sing[a] > TOL or sing[b] > TOL: continue
        h = hull(cpdag, K, (a, b), sem, x, y, cache, base)
        if (h[1] - h[0]) - base_w > TOL:
            syn.append((a, b)); ell_syn += [ell(a), ell(b)]
    U1 = set(L1); U2 = U1 | {c for p in syn for c in p}
    hw = lambda S: (lambda h: h[1] - h[0])(hull(cpdag, K, sorted(S), sem, x, y, cache, base))
    rows.append(dict(nk=len(K), nL1=len(L1), nsyn=len(syn),
                     w0=base_w, w1=hw(U1), w2=hw(U2), wK=hw(K)))

R = {k: np.array([r[k] for r in rows]) for k in rows[0]}
print(f"problemas: {got} de {tried};  |K| medio {R['nk'].mean():.2f};  "
      f"Theta(G0) e um ponto em {100*base_pt/got:.1f}% dos casos")
print(f"\nlambda positivo: n={len(lam_pos)}  min={min(lam_pos):.3g}  "
      f"mediana={np.median(lam_pos):.3g}  max={max(lam_pos):.3g}")
print(f"  -> menor lambda positivo observado = {min(lam_pos):.3g}; a tolerancia dentro de "
      f"'lambda>0' so vira um botao acima disso")
print(f"\nhop ate a query:  claims de L_1 media {np.mean(ell_L1):.2f} (n={len(ell_L1)})  |  "
      f"claims dos pares sinergicos media {np.mean(ell_syn) if ell_syn else float('nan'):.2f} (n={len(ell_syn)})")

def show(name, m):
    if m.sum() == 0: print(f"{name:34s} n=0"); return
    print(f"{name:34s} n={m.sum():3d}  w1={R['w1'][m].mean():.4f}  w2={R['w2'][m].mean():.4f}  "
          f"wK={R['wK'][m].mean():.4f}  w1/wK={np.mean(R['w1'][m]/np.maximum(R['wK'][m],1e-12)):.3f}  "
          f"w2/wK={np.mean(R['w2'][m]/np.maximum(R['wK'][m],1e-12)):.3f}")
print()
show("todos", np.ones(len(rows), bool))
show("L_1 vazio", R['nL1'] == 0)
show("L_1 nao vazio", R['nL1'] > 0)
show("com par sinergico (L_2 > L_1)", R['nsyn'] > 0)
show("L_1 vazio E sinergia presente", (R['nL1'] == 0) & (R['nsyn'] > 0))
print(f"\nP(L_1 vazio) = {np.mean(R['nL1']==0):.3f}   P(sinergia) = {np.mean(R['nsyn']>0):.3f}   "
      f"P(L_1 vazio | sinergia) = {np.mean(R['nL1'][R['nsyn']>0]==0):.3f}")
print(f"P(wK == 0, ou seja nada e identificavel-sensivel) = {np.mean(R['wK']<=TOL):.3f}")
