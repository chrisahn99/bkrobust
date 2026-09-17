#!/usr/bin/env python3
"""Audit 2: restrict to the regime the paper is about (background knowledge is
load-bearing: w(Theta(Chat)) > 0), diagnose fibre triviality, and test coverage."""
import sys, itertools, time
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
import numpy as np
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, optimal_adjustment_set_dag,
                                    is_valid_adjustment_set_mpdag, adjusted_estimand, random_sem)
TOL = 1e-9

def random_dag(rng, n, p):
    order = list(rng.permutation(n)); edges = set()
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < p: edges.add((f"v{order[i]}", f"v{order[j]}"))
    return MPDAG(nodes=[f"v{i}" for i in range(n)], directed=edges, undirected=[])

def theta_set(sem, g, x, y):
    exts = enumerate_dag_extensions(g)
    return [adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(d, x, y)) for d in exts], len(exts)

def width(v): return (max(v)-min(v)) if v else 0.0

def fibre(cpdag, K, k):
    rest = [e for e in K if e != k]; out = []
    g = apply_orientations(cpdag, rest) if rest else cpdag
    if g is not None: out.append(("retract", g))
    g2 = apply_orientations(cpdag, rest + [(k[1], k[0])])
    if g2 is not None: out.append(("reverse", g2))
    return out

def make(rng, n=9, p=0.30, k_claims=4, max_und=7, require_nonid=True):
    dag = random_dag(rng, n, p)
    if not dag.directed_edges: return None
    cpdag = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cpdag.undirected_edges)
    if not (2 <= len(und) <= max_und): return None
    nodes = sorted(dag.nodes)
    cand = [(x, y) for x in nodes for y in nodes if x != y]; rng.shuffle(cand)
    sem = random_sem(dag, rng)
    for x, y in cand:
        if require_nonid:
            thC, _ = theta_set(sem, cpdag, x, y)
            if width(thC) <= TOL:      # CPDAG already identifies the effect: BK irrelevant
                continue
        idx = rng.permutation(len(und))[:min(k_claims, len(und))]
        K = [((a, b) if (a, b) in dag.directed_edges else (b, a))
             for a, b in (und[t] for t in idx)]
        if not K: continue
        g0 = apply_orientations(cpdag, K)
        if g0 is None: continue
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None: continue
        return dict(dag=dag, cpdag=cpdag, x=x, y=y, K=[tuple(e) for e in K],
                    g0=g0, z=frozenset(z), und=und, sem=sem)
    return None

def analyse(P, label):
    print(f"\n########## {label}  (n_problems={len(P)}) ##########")
    lam0 = lampos = 0
    triv_retract = triv_reverse_none = fibre_trivial = n_claims = 0
    zero_hull = 0
    w_b1, w_bl, w_th0 = [], [], []
    lam0_zbreak = 0; lam0_tot = 0
    lams_all = []
    for pr in P:
        sem, x, y, cp, K, g0 = pr["sem"], pr["x"], pr["y"], pr["cpdag"], pr["K"], pr["g0"]
        th0, _ = theta_set(sem, g0, x, y); w0 = width(th0)
        thC, _ = theta_set(sem, cp, x, y)
        w_th0.append(w0); w_bl.append(width(thC))
        union = list(th0)
        for k in K:
            n_claims += 1
            fb = fibre(cp, K, k)
            tags = dict(fb)
            trivial = True
            if "retract" in tags:
                if tags["retract"] == g0: triv_retract += 1
                else: trivial = False
            if "reverse" not in tags: triv_reverse_none += 1
            else: trivial = False
            if trivial: fibre_trivial += 1
            fvals, zb = [], False
            for tag, g in fb:
                v, _ = theta_set(sem, g, x, y); fvals.extend(v)
                if not is_valid_adjustment_set_mpdag(g, x, y, pr["z"]): zb = True
            lam = width(th0 + fvals) - w0
            lams_all.append(lam)
            union.extend(fvals)
            if lam > TOL: lampos += 1
            else:
                lam0 += 1; lam0_tot += 1
                if zb: lam0_zbreak += 1
        wb1 = width(union); w_b1.append(wb1)
        if wb1 <= TOL: zero_hull += 1
    npb = len(P)
    print(f"lambda==0 : {lam0}/{n_claims} claims ({100*lam0/n_claims:.1f}%), "
          f"{lam0/npb:.2f} per problem of {n_claims/npb:.2f}")
    print(f"problems with an ALL-ZERO lambda vector (ranked list is fully tied): "
          f"{sum(1 for i in range(npb)):d} -> see zero_hull")
    print(f"budget-1 hull is a POINT on {zero_hull}/{npb} ({100*zero_hull/npb:.1f}%) problems")
    print(f"FIBRE TRIVIALITY: cl(K\\k)==G0 on {triv_retract}/{n_claims} "
          f"({100*triv_retract/n_claims:.1f}%);  cl(K[k<-rev])==None on "
          f"{triv_reverse_none}/{n_claims} ({100*triv_reverse_none/n_claims:.1f}%);  "
          f"BOTH (fibre adds nothing at all) on {fibre_trivial}/{n_claims} "
          f"({100*fibre_trivial/n_claims:.1f}%)")
    print(f"lambda==0 claims whose fibre nonetheless BREAKS Z: {lam0_zbreak}/{lam0_tot} "
          f"({100*lam0_zbreak/max(lam0_tot,1):.1f}%)")
    print(f"mean w(Theta(G0))={np.mean(w_th0):.5f}  mean w(budget-1 hull)={np.mean(w_b1):.5f}  "
          f"mean w(blanket)={np.mean(w_bl):.5f}")
    r = [a/b for a, b in zip(w_b1, w_bl) if b > TOL]
    print(f"budget-1 hull / blanket width ratio: mean {np.mean(r):.3f}, median {np.median(r):.3f} (n={len(r)})")
    nz = [l for l in lams_all if l > TOL]
    if nz: print(f"nonzero lambda values: n={len(nz)} min={min(nz):.4f} med={np.median(nz):.4f} max={max(nz):.4f}")
    return dict(w_b1=w_b1, w_bl=w_bl, zero_hull=zero_hull)

def coverage(P, budget):
    """Analyst states K with `budget` claims REVERSED (wrong). Does the hedge cover?"""
    rng = np.random.default_rng(7)
    n_ok = cov_b1 = cov_blanket = cov_point = usable = 0
    w_b1s, w_bls = [], []
    for pr in P:
        sem, x, y, cp, Ktrue = pr["sem"], pr["x"], pr["y"], pr["cpdag"], pr["K"]
        if len(Ktrue) < budget: continue
        pick = set(rng.permutation(len(Ktrue))[:budget].tolist())
        Kst = [((e[1], e[0]) if i in pick else e) for i, e in enumerate(Ktrue)]
        g0 = apply_orientations(cp, Kst)
        if g0 is None: continue                      # self-contradictory: analyst sees it
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None: continue
        usable += 1
        truth = sem.true_total_effect(x, y)
        th0, _ = theta_set(sem, g0, x, y)
        union = list(th0)
        for k in Kst:
            for tag, g in fibre(cp, Kst, k):
                v, _ = theta_set(sem, g, x, y); union.extend(v)
        thC, _ = theta_set(sem, cp, x, y)
        lo, hi = min(union), max(union)
        if lo - 1e-8 <= truth <= hi + 1e-8: cov_b1 += 1
        if min(thC) - 1e-8 <= truth <= max(thC) + 1e-8: cov_blanket += 1
        pt = adjusted_estimand(sem, x, y, z)
        if abs(pt - truth) < 1e-8: cov_point += 1
        w_b1s.append(hi - lo); w_bls.append(width(thC))
    print(f"budget={budget} wrong claim(s): usable {usable}")
    if usable:
        print(f"   coverage  point estimate    : {cov_point/usable:.3f}")
        print(f"   coverage  lambda budget-1 hull: {cov_b1/usable:.3f}   mean width {np.mean(w_b1s):.5f}")
        print(f"   coverage  blanket hedge     : {cov_blanket/usable:.3f}   mean width {np.mean(w_bls):.5f}")

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    for require in (False, True):
        rng = np.random.default_rng(20260907)
        P, tried = [], 0
        while len(P) < N and tried < 400000:
            tried += 1
            pr = make(rng, require_nonid=require)
            if pr is not None: P.append(pr)
        lab = ("REGIME B: only problems where the CPDAG does NOT already identify the effect"
               if require else "REGIME A: the generator as the harness writes it (no filter)")
        print(f"\n[{lab}] {len(P)} problems from {tried} draws")
        analyse(P, lab)
        for b in (0, 1, 2):
            coverage(P, b)
