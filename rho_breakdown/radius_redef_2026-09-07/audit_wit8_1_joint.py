#!/usr/bin/env python3
"""witness-8 audit: the joint load-bearing cover L_b.

Measures, on the repo's own random-problem harness:
  (A) is lambda monotone in A?                      [the unstated load-bearing lemma]
  (B) |L_1| / |K|                                   [does the degeneracy migrate?]
  (C) is L_2 \\ L_1 empty?                           [the candidate's own expected failure]
  (D) what does greedy-from-positive-singletons reach?
  (E) width of the prefix hedge at b=1, b=2, blanket-over-K, blanket-over-class
  (F) sensitivity of |L_1| to the tolerance hidden inside "lambda(A) > 0"
"""
import sys, itertools
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from collections import Counter
from harness import make_problem
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import optimal_adjustment_set_dag, adjusted_estimand, random_sem

TOL = 1e-9

def theta(cpdag, orients, sem, x, y, cache):
    """Set of TRUE total effects over the DAG extensions of Meek(cpdag, orients)."""
    key = frozenset(orients)
    if key in cache:
        return cache[key]
    g = apply_orientations(cpdag, list(orients))
    if g is None:
        vals = frozenset()
    else:
        vs = []
        for d in enumerate_dag_extensions(g):
            vs.append(adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(d, x, y)))
        vals = frozenset(np.round(vs, 12))
    cache[key] = vals
    return vals

def profiles(K, A):
    """All sigma in {retract, reverse}^A -> revised orientation list."""
    rest = [k for k in K if k not in A]
    A = list(A)
    for bits in itertools.product((0, 1), repeat=len(A)):
        rev = [(a[1], a[0]) for a, b in zip(A, bits) if b == 1]
        yield rest + rev

def lam(cpdag, K, A, sem, x, y, cache, base):
    vals = set(base)
    for orients in profiles(K, A):
        vals |= theta(cpdag, orients, sem, x, y, cache)
    if not vals:
        return 0.0
    return (max(vals) - min(vals)) - (max(base) - min(base) if base else 0.0)

rng = np.random.default_rng(20260907)
got = tried = 0
mono_viol = 0; mono_tests = 0
frac_L1 = []; n_pairs_new = []; new_claims = []
w1 = []; w2 = []; wK = []; wCls = []; w0 = []
tolsweep = {t: [] for t in (0.0, 1e-12, 1e-9, 1e-6, 1e-3, 1e-2)}
K_sizes = []

while got < 250 and tried < 60000:
    tried += 1
    pr = make_problem(rng, n=7, p=0.35, k_claims=5)
    if pr is None:
        continue
    K = [tuple(e) for e in pr["K"]]
    if len(K) < 4:
        continue
    cpdag, dag, x, y = pr["cpdag"], pr["dag"], pr["x"], pr["y"]
    sem = random_sem(dag, rng)
    cache = {}
    base = theta(cpdag, K, sem, x, y, cache)          # Theta(G0)
    if not base:
        continue
    got += 1
    K_sizes.append(len(K))
    base_w = max(base) - min(base)
    w0.append(base_w)

    singles = {k: lam(cpdag, K, (k,), sem, x, y, cache, base) for k in K}
    pairs = {}
    for a, b in itertools.combinations(K, 2):
        pairs[(a, b)] = lam(cpdag, K, (a, b), sem, x, y, cache, base)
        mono_tests += 1
        if pairs[(a, b)] < max(singles[a], singles[b]) - 1e-9:
            mono_viol += 1

    L1 = [k for k in K if singles[k] > TOL]
    frac_L1.append(len(L1) / len(K))
    for t in tolsweep:
        tolsweep[t].append(sum(1 for k in K if singles[k] > t) / len(K))

    # L_2 \ L_1 : minimal positive pairs = both singletons zero, pair positive
    synergy = [(a, b) for (a, b), v in pairs.items()
               if v > TOL and singles[a] <= TOL and singles[b] <= TOL]
    n_pairs_new.append(len(synergy))
    U1 = set(L1)
    U2 = U1 | {c for p in synergy for c in p}
    new_claims.append(len(U2) - len(U1))

    # widths of the prefix hedge over a claim SET (all joint profiles on that set)
    def hedge_w(S):
        if not S:
            return base_w
        vals = set(base)
        for orients in profiles(K, tuple(S)):
            vals |= theta(cpdag, orients, sem, x, y, cache)
        return max(vals) - min(vals)
    w1.append(hedge_w(sorted(U1)))
    w2.append(hedge_w(sorted(U2)))
    wK.append(hedge_w(K))
    cls = theta(cpdag, [], sem, x, y, cache)           # blanket over the whole class
    wCls.append((max(cls) - min(cls)) if cls else float("nan"))

f = lambda v: f"{np.mean(v):.4f}"
print(f"problemas usaveis: {got} de {tried} sorteios;  |K| medio {np.mean(K_sizes):.2f}")
print()
print("(A) monotonicidade de lambda:", f"{mono_viol} violacoes em {mono_tests} pares testados")
print("(B) |L_1|/|K|: media", f(frac_L1), " P(L_1 = K) =",
      f"{np.mean([v == 1.0 for v in frac_L1]):.3f}",
      " P(L_1 vazio) =", f"{np.mean([v == 0.0 for v in frac_L1]):.3f}")
print("(C) pares sinergicos (L_2 \\ L_1): total", int(np.sum(n_pairs_new)),
      " em", got, "problemas;  P(>=1) =", f"{np.mean([v > 0 for v in n_pairs_new]):.3f}")
print("(D) claims que b=2 acrescenta a U L_1: total", int(np.sum(new_claims)),
      " media", f(new_claims))
print("(E) larguras medias:  Theta(G0)", f(w0), " hedge b=1", f(w1),
      " hedge b=2", f(w2), " hedge K inteiro", f(wK), " classe inteira", f(wCls))
print("    razao hedge_b1 / classe:", f"{np.mean(np.array(w1)/np.array(wCls)):.3f}",
      "   hedge_K / classe:", f"{np.mean(np.array(wK)/np.array(wCls)):.3f}")
print("(F) |L_1|/|K| por tolerancia:")
for t in sorted(tolsweep):
    print(f"      tol={t:<8g}  {np.mean(tolsweep[t]):.4f}")
