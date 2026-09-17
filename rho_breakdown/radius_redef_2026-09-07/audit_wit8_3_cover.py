#!/usr/bin/env python3
"""witness-8 audit 3: coverage and width when 1 or 2 elicited claims are FALSE.
Compares the hedge over the FLATTENED union (what the proposal literally says)
against the hedge over the FAMILY (only |A| <= b revised at a time)."""
import sys, itertools
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
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

def run(nfalse, ntarget=400, seed=20260907):
    rng = np.random.default_rng(seed)
    got = tried = 0
    acc = {k: [] for k in ("w1f", "w1x", "w2f", "w2x", "wK",
                           "c1f", "c1x", "c2f", "c2x", "cK", "nL1", "nsyn")}
    while got < ntarget and tried < 120000:
        tried += 1
        pr = make_problem(rng, n=7, p=0.35, k_claims=5)
        if pr is None: continue
        K0 = [tuple(e) for e in pr["K"]]
        if len(K0) < max(4, nfalse + 1): continue
        cpdag, dag, x, y = pr["cpdag"], pr["dag"], pr["x"], pr["y"]
        bad = [K0[i] for i in rng.permutation(len(K0))[:nfalse]]
        K = [((k[1], k[0]) if k in bad else k) for k in K0]     # analyst's FALSE knowledge
        if apply_orientations(cpdag, K) is None: continue
        sem = random_sem(dag, rng); truth = sem.true_total_effect(x, y)
        cache = {}
        base = theta(cpdag, K, sem, x, y, cache)
        if not base: continue
        got += 1
        bw = max(base) - min(base)
        sing, parts = {}, [set(base)]
        for k in K:
            v = set()
            for o in profiles(K, (k,)): v |= theta(cpdag, o, sem, x, y, cache)
            sing[k] = (max(v | set(base)) - min(v | set(base))) - bw
            if sing[k] > TOL: parts.append(v)
        L1 = [k for k in K if sing[k] > TOL]
        syn = []
        for a, b in itertools.combinations(K, 2):
            if sing[a] > TOL or sing[b] > TOL: continue
            v = set()
            for o in profiles(K, (a, b)): v |= theta(cpdag, o, sem, x, y, cache)
            if (max(v | set(base)) - min(v | set(base))) - bw > TOL:
                syn.append((a, b)); parts.append(v)
        U1, U2 = set(L1), set(L1) | {c for p in syn for c in p}
        fam1 = set().union(*parts[:1 + len(L1)])
        fam2 = set().union(*parts)
        def flat(S):
            v = set(base)
            for o in profiles(K, tuple(sorted(S))): v |= theta(cpdag, o, sem, x, y, cache)
            return v
        for tag, v in (("1f", fam1), ("1x", flat(U1)), ("2f", fam2),
                       ("2x", flat(U2)), ("K", flat(K))):
            acc["w" + tag].append(max(v) - min(v))
            acc["c" + tag].append(min(v) - 1e-9 <= truth <= max(v) + 1e-9)
        acc["nL1"].append(len(L1)); acc["nsyn"].append(len(syn))
    A = {k: np.array(v) for k, v in acc.items()}
    print(f"\n=== {nfalse} claim(s) FALSA(s):  {got} problemas de {tried} sorteios ===")
    print(f"P(L_1 vazio) = {np.mean(A['nL1']==0):.3f}   P(par sinergico) = {np.mean(A['nsyn']>0):.3f}")
    info = A["wK"] > TOL
    print(f"estrato informativo (largura do hedge-cobertor > 0): {info.sum()}/{got} = {info.mean():.3f}")
    hdr = (("hedge familia b=1", "1f"), ("hedge U L_1 achatado", "1x"),
           ("hedge familia b=2", "2f"), ("hedge U L_2 achatado", "2x"),
           ("cobertor (classe toda)", "K"))
    for name, t in hdr:
        w, c = A["w" + t], A["c" + t]
        print(f"  {name:24s} cobertura {c.mean():.3f}   largura media {w.mean():.4f}"
              f"   | no estrato informativo: cob {c[info].mean():.3f}  larg {w[info].mean():.4f}"
              f"  larg/cobertor {np.mean(w[info]/A['wK'][info]):.3f}")

run(1)
run(2)
