#!/usr/bin/env python3
"""Audit of witness-2: claim leverage lambda(k).

Tests, in order:
  T1  sparsity of lambda (the claimed 2.97/3.42 zeros)
  T2  does lambda=0 really mean "revising k changes no identified effect"?
      and does lambda=0 coincide with "Z stays valid" (the r_val hop-1 event)?
  T3  width of the delta=0 prefix hedge vs the FULL budget-1 fibre hull
      (i.e. does lambda do any work at delta=0?)
  T4  zero-width hedges: how often, and can they be calibrated by scaling?
  T5  cost: sum over fibres of |extensions| vs |extensions(Chat)|
  T6  coverage under 0/1/2 wrong claims, lambda-hedge vs blanket hedge
  T7  is lambda just r_val renamed?
"""
import sys, itertools, time
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
import numpy as np

from bkrobust.demo.graph import MPDAG
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.space import (enumerate_space, represented_dags, covering_pairs,
                                 neighbour_graph, bfs_distances)
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag,
                                    optimal_adjustment_set_dag,
                                    is_valid_adjustment_set_mpdag,
                                    adjusted_estimand, random_sem)

TOL = 1e-9


def random_dag(rng, n, p):
    order = list(rng.permutation(n))
    edges = set()
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < p:
                edges.add((f"v{order[i]}", f"v{order[j]}"))
    return MPDAG(nodes=[f"v{i}" for i in range(n)], directed=edges, undirected=[])


def make_problem(rng, n=9, p=0.30, k_claims=4, max_und=6):
    dag = random_dag(rng, n, p)
    if not dag.directed_edges:
        return None
    cpdag = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cpdag.undirected_edges)
    if not (2 <= len(und) <= max_und):
        return None
    nodes = sorted(dag.nodes)
    cand = [(x, y) for x in nodes for y in nodes if x != y]
    rng.shuffle(cand)
    for x, y in cand:
        idx = rng.permutation(len(und))[:min(k_claims, len(und))]
        K = []
        for t in idx:
            a, b = und[t]
            K.append((a, b) if (a, b) in dag.directed_edges else (b, a))
        if not K:
            continue
        g0 = apply_orientations(cpdag, K)
        if g0 is None:
            continue
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None:
            continue
        sem = random_sem(dag, rng)
        return dict(dag=dag, cpdag=cpdag, x=x, y=y, K=[tuple(e) for e in K],
                    g0=g0, z=frozenset(z), und=und, sem=sem)
    return None


def theta_set(sem, g, x, y):
    """Theta(G) = { theta_D : D in extensions(G) }, plus the extension count."""
    exts = enumerate_dag_extensions(g)
    vals = []
    for d in exts:
        o = optimal_adjustment_set_dag(d, x, y)
        vals.append(adjusted_estimand(sem, x, y, o))
    return vals, len(exts)


def width(vals):
    return (max(vals) - min(vals)) if vals else 0.0


def fibre(cpdag, K, k):
    """F(k) = { cl(K \\ {k}), cl(K[k<-rev]) } minus None."""
    rest = [e for e in K if e != k]
    out = []
    g_ret = apply_orientations(cpdag, rest) if rest else cpdag
    if g_ret is not None:
        out.append(("retract", g_ret))
    rev = rest + [(k[1], k[0])]
    g_rev = apply_orientations(cpdag, rev)
    if g_rev is not None:
        out.append(("reverse", g_rev))
    return out


def run(n_problems=200, n=9, seed=20260907):
    rng = np.random.default_rng(seed)
    P = []
    tried = 0
    while len(P) < n_problems and tried < 200000:
        tried += 1
        pr = make_problem(rng, n=n)
        if pr is not None:
            P.append(pr)
    print(f"# {len(P)} usable problems from {tried} draws, n={n} nodes")

    lam_zero, lam_pos, per_problem_claims = 0, 0, []
    # T2 accumulators
    t2_lam0_theta_moves = 0        # lambda=0 but some fibre theta != theta_0
    t2_lam0_z_breaks = 0           # lambda=0 but Z invalid somewhere in fibre
    t2_lam0_total = 0
    t2_lampos_z_ok = 0             # lambda>0 but Z stays valid everywhere
    t2_lampos_total = 0
    # T3
    w_delta0, w_full_b1, w_blanket, w_theta0 = [], [], [], []
    # T5
    cost_fibres, cost_blanket = [], []
    # T4
    zero_width = 0
    # T7
    rval_is_1 = 0
    maxlam_pos = 0
    agree_rval_maxlam = 0

    t0 = time.time()
    for pr in P:
        sem, x, y, cpdag, K, g0 = pr["sem"], pr["x"], pr["y"], pr["cpdag"], pr["K"], pr["g0"]
        th0, n0 = theta_set(sem, g0, x, y)
        w0 = width(th0)
        w_theta0.append(w0)
        thC, nC = theta_set(sem, cpdag, x, y)
        w_blanket.append(width(thC))
        cost_blanket.append(nC)

        per_problem_claims.append(len(K))
        union_all = list(th0)
        union_pos = list(th0)
        cost = 0
        any_lam_pos = False
        for k in K:
            fvals = []
            zbreak = False
            for tag, g in fibre(cpdag, K, k):
                vals, ne = theta_set(sem, g, x, y)
                cost += ne
                fvals.extend(vals)
                if not is_valid_adjustment_set_mpdag(g, x, y, pr["z"]):
                    zbreak = True
            lam = width(th0 + fvals) - w0
            union_all.extend(fvals)
            if lam > TOL:
                lam_pos += 1
                any_lam_pos = True
                union_pos.extend(fvals)
                t2_lampos_total += 1
                if not zbreak:
                    t2_lampos_z_ok += 1
            else:
                lam_zero += 1
                t2_lam0_total += 1
                if zbreak:
                    t2_lam0_z_breaks += 1
                # does any fibre theta differ from theta_0 (which is a singleton here)?
                if fvals and (max(fvals) - min(th0) > TOL or min(th0) - min(fvals) > TOL
                              or max(th0) - min(fvals) > TOL or max(fvals) - max(th0) > TOL):
                    t2_lam0_theta_moves += 1
        cost_fibres.append(cost)
        w_delta0.append(width(union_pos))
        w_full_b1.append(width(union_all))
        if width(union_all) <= TOL:
            zero_width += 1
        if any_lam_pos:
            maxlam_pos += 1

        # T7: r_val
        space = enumerate_space(cpdag)
        reps = represented_dags(space)
        covers = covering_pairs(space, reps)
        nbrs = neighbour_graph(space, covers)
        g0m = next((g for g in space if g == g0), None)
        rv = None
        if g0m is not None:
            dist = bfs_distances(nbrs, g0m)
            bad = [dist[g] for g in space if g in dist and dist[g] > 0
                   and not is_valid_adjustment_set_mpdag(g, x, y, pr["z"])]
            rv = min(bad) if bad else None
        if rv == 1:
            rval_is_1 += 1
        if (rv == 1) == any_lam_pos:
            agree_rval_maxlam += 1
        pr["_rv"] = rv
        pr["_anylam"] = any_lam_pos

    el = time.time() - t0
    tot_claims = lam_zero + lam_pos
    npb = len(P)
    print(f"# sweep wall time {el:.1f} s")
    print()
    print("== T1 sparsity ==")
    print(f"claims per problem      : {tot_claims/npb:.2f}")
    print(f"lambda == 0 per problem : {lam_zero/npb:.2f}   ({100*lam_zero/tot_claims:.1f}% of claims)")
    print(f"lambda  > 0 per problem : {lam_pos/npb:.2f}")
    print()
    print("== T2 what lambda==0 actually means ==")
    print(f"lambda==0 claims                       : {t2_lam0_total}")
    print(f"  ... whose fibre theta differs from t0: {t2_lam0_theta_moves}  "
          f"({100*t2_lam0_theta_moves/max(t2_lam0_total,1):.1f}%)")
    print(f"  ... whose fibre BREAKS Z (r_val hop1): {t2_lam0_z_breaks}  "
          f"({100*t2_lam0_z_breaks/max(t2_lam0_total,1):.1f}%)")
    print(f"lambda>0 claims                        : {t2_lampos_total}")
    print(f"  ... whose fibre keeps Z valid        : {t2_lampos_z_ok}  "
          f"({100*t2_lampos_z_ok/max(t2_lampos_total,1):.1f}%)")
    print()
    print("== T3 does lambda do any work at delta=0? ==")
    d = [abs(a-b) for a, b in zip(w_delta0, w_full_b1)]
    print(f"mean width, delta=0 prefix hedge : {np.mean(w_delta0):.6f}")
    print(f"mean width, FULL budget-1 hull   : {np.mean(w_full_b1):.6f}")
    print(f"max |difference| over problems   : {max(d):.3e}   (0 => lambda is decorative at delta=0)")
    print(f"mean width, w(Theta(G0))         : {np.mean(w_theta0):.6f}  "
          f"(nonzero on {sum(1 for v in w_theta0 if v>TOL)}/{npb} problems)")
    print(f"mean width, blanket hedge Theta(Chat): {np.mean(w_blanket):.6f}")
    rat = [a/b for a, b in zip(w_full_b1, w_blanket) if b > TOL]
    print(f"budget-1 hull / blanket, mean ratio: {np.mean(rat):.3f}  (n={len(rat)})")
    print()
    print("== T4 zero-width hedges ==")
    print(f"problems where the budget-1 hull is a POINT: {zero_width}/{npb}  "
          f"({100*zero_width/npb:.1f}%)")
    print()
    print("== T5 cost claim ==")
    print(f"mean sum |ext(fibre)| : {np.mean(cost_fibres):.1f}")
    print(f"mean |ext(Chat)|      : {np.mean(cost_blanket):.1f}")
    frac = sum(1 for a, b in zip(cost_fibres, cost_blanket) if a > b)
    print(f"problems where fibre cost EXCEEDS blanket cost: {frac}/{npb} "
          f"({100*frac/npb:.1f}%)  [definition claims 'strictly cheaper']")
    print()
    print("== T7 lambda vs r_val ==")
    print(f"r_val == 1 on                : {rval_is_1}/{npb} ({100*rval_is_1/npb:.1f}%)")
    print(f"max_k lambda(k) > 0 on       : {maxlam_pos}/{npb} ({100*maxlam_pos/npb:.1f}%)")
    print(f"the two indicators agree on  : {agree_rval_maxlam}/{npb} "
          f"({100*agree_rval_maxlam/npb:.1f}%)")
    return P


if __name__ == "__main__":
    run(n_problems=int(sys.argv[1]) if len(sys.argv) > 1 else 200)
