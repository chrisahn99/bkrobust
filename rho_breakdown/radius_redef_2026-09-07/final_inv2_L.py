#!/usr/bin/env python3
"""FINAL inverse-2: the load-bearing SET L (validity-first partition), its degeneracy
profile, padding stability, the r_val cross-tab, the completeness check, and the
proximity-augmented hedge width table."""
import sys, itertools
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, optimal_adjustment_set_dag,
                                    adjusted_estimand, is_valid_adjustment_set_mpdag,
                                    is_valid_adjustment_set_dag, random_sem)
from bkrobust.demo.space import (enumerate_space, represented_dags, covering_pairs,
                                 neighbour_graph, bfs_distances)
from harness import random_dag
P_ERR = 0.20

def rev(k): return (k[1], k[0])

def admissible(cpdag, K, k):
    """single-claim revisions of k that the CPDAG + the rest of K still permit"""
    rest = [e for e in K if e != k]
    out = []
    for extra in ([], [rev(k)]):
        g = apply_orientations(cpdag, rest + extra)
        if g is not None:
            out.append(("retract" if not extra else "reverse", g))
    return out

def partition(cpdag, K, x, y, z, g0):
    """VALIDITY FIRST, residual labelled afterwards; data-refuted against Chat alone."""
    L, free, inert, dref = [], [], [], []
    for k in K:
        adm = admissible(cpdag, K, k)
        if any(not is_valid_adjustment_set_mpdag(g, x, y, z) for _, g in adm):
            L.append(k)                                    # load-bearing
        elif apply_orientations(cpdag, [rev(k)]) is None:
            dref.append(k)                                 # refuted by Chat ALONE
        elif any(tag == "retract" and g == g0 for tag, g in adm):
            inert.append(k)                                # Meek re-derives it
        else:
            free.append(k)
    return L, free, inert, dref

def pairs_revisions(cpdag, K, S):
    out = []
    for k1, k2 in itertools.combinations(sorted(set(S)), 2):
        rest = [e for e in K if e not in (k1, k2)]
        for e1 in ([], [rev(k1)]):
            for e2 in ([], [rev(k2)]):
                g = apply_orientations(cpdag, rest + e1 + e2)
                if g is not None:
                    out.append(g)
    return out

def problem(rng, n=6, p=0.35, kmax=4):
    dag = random_dag(rng, n, p)
    if not dag.directed_edges: return None
    cpdag = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cpdag.undirected_edges)
    if not (2 <= len(und) <= 5): return None
    nodes = sorted(dag.nodes)
    cand = [(a, b) for a in nodes for b in nodes if a != b and b in dag.descendants(a)]
    rng.shuffle(cand)
    for x, y in cand:
        idx = rng.permutation(len(und))[:min(kmax, len(und))]
        K, wrong = [], []
        for t in idx:
            a, b = und[t]
            tr = (a, b) if (a, b) in dag.directed_edges else (b, a)
            if rng.random() < P_ERR: K.append(rev(tr)); wrong.append(rev(tr))
            else: K.append(tr)
        g0 = apply_orientations(cpdag, K)
        if g0 is None: continue
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None: continue
        return dict(dag=dag, cpdag=cpdag, x=x, y=y, K=K, wrong=set(wrong), g0=g0,
                    z=frozenset(z), und=und)
    return None

def r_val_of(pr):
    space = enumerate_space(pr["cpdag"])
    reps = represented_dags(space); covers = covering_pairs(space, reps)
    nbrs = neighbour_graph(space, covers)
    g0 = next((g for g in space if g == pr["g0"]), None)
    if g0 is None: return "NA"
    dist = bfs_distances(nbrs, g0)
    bad = [dist[g] for g in space if dist.get(g) not in (None, 0)
           and not is_valid_adjustment_set_mpdag(g, pr["x"], pr["y"], pr["z"])]
    return min(bad) if bad else "UNREACHED"

rng = np.random.default_rng(20260907)
R, tried = [], 0
pad_L_before, pad_L_after, pad_phi_before, pad_phi_after = [], [], [], []
comp_viol = 0; comp_checked = 0
while len(R) < 300 and tried < 60000:
    tried += 1
    pr = problem(rng)
    if pr is None: continue
    cpdag, K, x, y, z, g0, dag = (pr["cpdag"], pr["K"], pr["x"], pr["y"], pr["z"],
                                  pr["g0"], pr["dag"])
    L, F, I, D = partition(cpdag, K, x, y, z, g0)
    P = [k for k in K if k[0] in (x, y) or k[1] in (x, y)]        # hop-0 claims
    S = sorted(set(L) | set(P))
    sem = random_sem(dag, rng); cache = {}
    def thetas(g):
        key = g.edge_string()
        if key not in cache:
            cache[key] = [adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(d, x, y))
                          for d in enumerate_dag_extensions(g)]
        return cache[key]
    point = adjusted_estimand(sem, x, y, z); tau = sem.true_total_effect(x, y)
    def hedge(graphs):
        vals = [point]
        for g in graphs: vals += thetas(g)
        return min(vals), max(vals)
    singles = lambda C: [g for k in C for _, g in admissible(cpdag, K, k)]
    h_L   = hedge(singles(L))
    h_K   = hedge(singles(K))
    h_S2  = hedge(singles(S) + pairs_revisions(cpdag, K, S))       # PROPOSED
    h_K2  = hedge(singles(K) + pairs_revisions(cpdag, K, K))       # all claims, order 2
    h_bla = hedge([cpdag])
    # completeness of the reordered partition: k not in L => every admissible revision
    # keeps Z valid in EVERY DAG extension
    for k in K:
        if k in L: continue
        for _, g in admissible(cpdag, K, k):
            comp_checked += 1
            if any(not is_valid_adjustment_set_dag(d, x, y, z)
                   for d in enumerate_dag_extensions(g)): comp_viol += 1
    # padding attack
    ent = [tuple(e) for e in sorted(g0.directed_edges) if tuple(e) not in set(K)]
    if ent and apply_orientations(cpdag, list(K) + ent) == g0:
        Kp = list(K) + ent
        L2, F2, I2, D2 = partition(cpdag, Kp, x, y, z, g0)
        den  = sum(1 for k in K  if admissible(cpdag, K,  k))
        den2 = sum(1 for k in Kp if admissible(cpdag, Kp, k))
        pad_L_before.append(len(L)); pad_L_after.append(len(L2))
        pad_phi_before.append(len(L)/den if den else 0.0)
        pad_phi_after.append(len(L2)/den2 if den2 else 0.0)
    R.append(dict(tau=tau, point=point, nL=len(L), nF=len(F), nI=len(I), nD=len(D),
                  nK=len(K), nP=len(P), nS=len(S), phi=len(L)/max(1,sum(1 for k in K if admissible(cpdag,K,k))),
                  hL=h_L, hK=h_K, hS2=h_S2, hK2=h_K2, hbl=h_bla,
                  rval=r_val_of(pr), hit=len(set(L) & pr["wrong"]), nwrong=len(pr["wrong"]),
                  zbad=not is_valid_adjustment_set_dag(dag, x, y, z)))

tau = np.array([r["tau"] for r in R]); pts = np.array([r["point"] for r in R])
def cov_width(key):
    lo = np.minimum(np.array([r[key][0] for r in R]), pts)
    hi = np.maximum(np.array([r[key][1] for r in R]), pts)
    half = (hi-lo)/2; mid = (lo+hi)/2; need = np.abs(tau-mid)
    raw = float(np.mean(need <= half + 1e-12))
    c = np.where(half > 1e-12, need/np.where(half > 1e-12, half, 1.0),
                 np.where(need < 1e-9, 0.0, np.inf))
    fin = np.isfinite(c)
    if fin.mean() < 0.95: return raw, None, float(np.mean(2*half)), 1-fin.mean()
    q = float(np.quantile(np.sort(c[fin]), min(0.95/fin.mean(), 1.0)))
    return raw, q, float(np.mean(2*half*q)), 1-fin.mean()

n = len(R)
print(f"n = {n} amenable problems (6 nodes, 2-5 undirected, |K|<=4, per-claim error 0.20)")
print("\n1. DEGENERACY PROFILE OF |L|")
c = Counter(r["nL"] for r in R)
for v in sorted(c): print(f"   |L| = {v}: {c[v]:3d}  {100*c[v]/n:5.1f}%")
print(f"   distinct values {len(c)}, mean {np.mean([r['nL'] for r in R]):.3f}, "
      f"max cell {100*max(c.values())/n:.1f}%")
print("\n2. CROSS-TAB WITH r_val")
ct = Counter((r["rval"], r["nL"] >= 1) for r in R)
for k in sorted(ct, key=lambda t: (str(t[0]), t[1])): print(f"   r_val={k[0]:<10} |L|>=1={k[1]}   {ct[k]}")
agree = sum(1 for r in R if (r["nL"] >= 1) == (r["rval"] == 1))
print(f"   (|L|>=1) == (r_val==1) on {agree}/{n} = {100*agree/n:.1f}%")
st = [r["nL"] for r in R if r["rval"] == 1]
if st:
    cs = Counter(st)
    print(f"   inside the r_val==1 stratum (n={len(st)}): |L| takes {len(cs)} values "
          f"{dict(sorted(cs.items()))}, max cell {100*max(cs.values())/len(st):.1f}%")
print("\n3. COMPLETENESS OF THE REORDERED PARTITION")
print(f"   revisions of non-flagged claims checked: {comp_checked}; violations: {comp_viol}")
print("\n4. PADDING ATTACK (append claims G0 already entails)")
if pad_L_before:
    b = np.array(pad_L_before); a = np.array(pad_L_after)
    pb = np.array(pad_phi_before); pa = np.array(pad_phi_after)
    print(f"   n={len(b)}   |L| changed on {100*np.mean(a!=b):.1f}%  (mean {b.mean():.3f} -> {a.mean():.3f})")
    print(f"              phi_1 changed on {100*np.mean(pa!=pb):.1f}%  (mean {pb.mean():.4f} -> {pa.mean():.4f})")
print("\n5. WIDTH AT HONEST 95% COVERAGE (lower is better)")
for name, key in (("hedge over L only (order 1)", "hL"),
                  ("hedge over every claim (order 1)", "hK"),
                  ("PROPOSED: order-1 on L u P, order-2 inside L u P", "hS2"),
                  ("every claim, order 2", "hK2"),
                  ("blanket hedge over the whole class", "hbl")):
    raw, q, w, fi = cov_width(key)
    print(f"   {name:52s} raw {raw:.3f}  scale {'n/a' if q is None else f'{q:5.2f}'}"
          f"  width {w:7.4f}  uncoverable {fi:.3f}")
pt_raw = float(np.mean(np.abs(tau-pts) < 1e-9))
print(f"   {'point estimate on the elicited set':52s} raw {pt_raw:.3f}  width 0.0000 (uncalibratable)")
print("\n6. RETRIEVAL SCORE OF L AGAINST THE ACTUALLY-WRONG CLAIMS")
tp = sum(r["hit"] for r in R); fp = sum(r["nL"] for r in R) - tp
fn = sum(r["nwrong"] for r in R) - tp
print(f"   micro precision {tp/max(1,tp+fp):.3f}  recall {tp/max(1,tp+fn):.3f}  "
      f"base rate {sum(r['nwrong'] for r in R)/sum(r['nK'] for r in R):.3f}")
print(f"\n   |P| mean {np.mean([r['nP'] for r in R]):.2f}, |S| mean {np.mean([r['nS'] for r in R]):.2f}, "
      f"|K| mean {np.mean([r['nK'] for r in R]):.2f}")
print(f"   partition totals: L {sum(r['nL'] for r in R)}, free {sum(r['nF'] for r in R)}, "
      f"inert {sum(r['nI'] for r in R)}, data-refuted(Chat-alone) {sum(r['nD'] for r in R)} "
      f"of {sum(r['nK'] for r in R)} claims")
