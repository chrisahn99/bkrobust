#!/usr/bin/env python3
"""FINAL DEFINITION -- the retraction hull ladder R_m, budget-calibrated.

Measures, on the regime where background knowledge does any work at all
(w(Theta(Chat)) > 0):
  1. budget-b soundness (exact, structural)   -- coverage at b >= b_true
  2. the ladder: coverage and width as a function of the admitted family
  3. the calibrated width vs the blanket hedge (the target-metric row)
  4. degeneracy: the distribution of nu = w(R)/w(Theta(Chat))
  5. cost, with the DAG cache: distinct theta evaluations vs the blanket's
"""
import sys, itertools
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np

from bkrobust.demo.graph import MPDAG
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_dag, optimal_adjustment_set_mpdag,
                                    is_valid_adjustment_set_dag, adjusted_estimand, random_sem)
from harness import random_dag

TOL = 1e-9


# ---------------------------------------------------------------- theta_D (audit fix 3)
def theta_D(sem, D, x, y, cache):
    key = D.key()
    if key in cache:
        return cache[key]
    # cn_D(x,y) empty  <=>  y is not a descendant of x in D  <=>  D asserts a zero
    # total effect.  This is the audit_rw_0_ostar_edge_case.py case: O*(D) is then
    # empty and empty is not a back-door set, so adjusting returns a confounded
    # coefficient instead of the zero D actually asserts.
    if y not in (D.descendants(x) & (D.ancestors(y) | {y})):
        cache[key] = 0.0
        return 0.0
    S = optimal_adjustment_set_dag(D, x, y)
    if not is_valid_adjustment_set_dag(D, x, y, S):
        S = D.parents(x)          # valid back-door set whenever y is a descendant of x
    v = adjusted_estimand(sem, x, y, S)
    cache[key] = v
    return v


def Theta(sem, G, x, y, cache):
    return [theta_D(sem, D, x, y, cache) for D in enumerate_dag_extensions(G)]


def hull(vals):
    return (min(vals), max(vals))


def w(iv):
    return iv[1] - iv[0]


def merge(a, b):
    return (min(a[0], b[0]), max(a[1], b[1]))


# ---------------------------------------------------------------- problem generator
def make(rng, n=9, p=0.30, kmax=4):
    dag = random_dag(rng, n, p)
    if not dag.directed_edges:
        return None
    cp = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cp.undirected_edges)
    if not (2 <= len(und) <= 7):
        return None
    nodes = sorted(dag.nodes)
    cand = [(x, y) for x in nodes for y in nodes if x != y]
    rng.shuffle(cand)
    for x, y in cand:
        km = min(kmax, len(und))
        idx = rng.permutation(len(und))[:km]
        K = []
        for t in idx:
            a, b = und[t]
            K.append((a, b) if (a, b) in dag.directed_edges else (b, a))
        g0 = apply_orientations(cp, K)
        if g0 is None:
            continue
        if optimal_adjustment_set_mpdag(g0, x, y) is None:
            continue
        sem = random_sem(dag, rng)
        cache = {}
        blanket = hull(Theta(sem, cp, x, y, cache))
        return dict(dag=dag, cpdag=cp, x=x, y=y, K=[tuple(e) for e in K], g0=g0,
                    sem=sem, blanket=blanket, cache=cache,
                    n_ext_blanket=len(enumerate_dag_extensions(cp)),
                    touches=any(a in (x, y) or b in (x, y) for a, b in K))
    return None


# ---------------------------------------------------------------- the ladder
def ladder(pr, Kstar, g0, bmax):
    """Ordered list of (S, interval, lam) for every retraction set S, |S| <= bmax.

    Order = (|S|, -lam(S)).  lam(S) = w(Theta(G0) u Theta(cl(K\\S))) - w(Theta(G0)).
    Returns (base_interval, entries, n_closures).
    """
    sem, x, y, cp, cache = pr["sem"], pr["x"], pr["y"], pr["cpdag"], pr["cache"]
    base = hull(Theta(sem, g0, x, y, cache))
    entries, ncl = [], 0
    for size in range(1, bmax + 1):
        for S in itertools.combinations(range(len(Kstar)), size):
            keep = [e for i, e in enumerate(Kstar) if i not in S]
            g = apply_orientations(cp, keep) if keep else cp
            ncl += 1
            if g is None:
                continue
            iv = hull(Theta(sem, g, x, y, cache))
            lam = w(merge(base, iv)) - w(base)
            entries.append((S, iv, lam, size))
    entries.sort(key=lambda t: (t[3], -t[2]))
    return base, entries, ncl


def run(seed=20260907, N=250, bmax=3, verbose=True):
    rng = np.random.default_rng(seed)
    P, tried = [], 0
    while len(P) < N and tried < 400000:
        tried += 1
        pr = make(rng)
        if pr is None:
            continue
        P.append(pr)
    nz = [pr for pr in P if w(pr["blanket"]) > TOL]
    if verbose:
        print(f"problems drawn: {len(P)}  (from {tried} draws)")
        print(f"REGIME FILTER w(Theta(Chat))>0 : {len(nz)}/{len(P)} = "
              f"{100*len(nz)/len(P):.1f}%  (the rest: background knowledge cannot matter)")
    return P, nz


if __name__ == "__main__":
    P, nz = run()
    # ---- 1. budget-b soundness, exact -------------------------------------------
    print("\n=== 1. BUDGET-b SOUNDNESS  (b_true false claims, R_b at b >= b_true) ===")
    rng = np.random.default_rng(1234)
    for b_true in (0, 1, 2):
        rows = {b: [0, 0] for b in range(0, 4)}   # [covered, used]
        widths = {b: [] for b in range(0, 4)}
        for pr in nz:
            K = pr["K"]
            if len(K) < b_true:
                continue
            pick = set(rng.permutation(len(K))[:b_true].tolist())
            Kst = [((e[1], e[0]) if i in pick else e) for i, e in enumerate(K)]
            g0 = apply_orientations(pr["cpdag"], Kst)
            if g0 is None:
                continue                       # a detectable (inconsistent) error
            if optimal_adjustment_set_mpdag(g0, pr["x"], pr["y"]) is None:
                continue
            t = pr["sem"].true_total_effect(pr["x"], pr["y"])
            base, entries, _ = ladder(pr, Kst, g0, min(3, len(Kst)))
            iv = base
            byb = {0: base}
            for S, e_iv, lam, size in sorted(entries, key=lambda z: z[3]):
                iv = merge(iv, e_iv)
                byb[size] = iv
            cur = base
            for b in range(0, 4):
                cur = byb.get(b, cur)
                rows[b][1] += 1
                if cur[0] - 1e-8 <= t <= cur[1] + 1e-8:
                    rows[b][0] += 1
                widths[b].append(w(cur))
        line = " | ".join(f"b={b}: cov={rows[b][0]/max(rows[b][1],1):.3f} "
                          f"wid={np.mean(widths[b]):.3f}" for b in range(0, 4))
        print(f"  b_true={b_true}  n={rows[0][1]:>3}  {line}")
