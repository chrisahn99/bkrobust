#!/usr/bin/env python3
"""Audit of r^alpha (weighted-3): the alpha-substantial minimum.

For each random problem we compute, for every element G of the space:
  dist(G0,G), beta(G) = nu_G({d in [G] : Z invalid in d}), |[G]|, and whether
  G is a refinement of G0 (downward), a coarsening (upward), or neither.
Then r^alpha for a grid of alpha, plus the diagnostics the audit needs.
"""
import sys, time, json
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np

from harness import make_problem
from bkrobust.demo.space import (enumerate_space, represented_dags, covering_pairs,
                                 neighbour_graph, bfs_distances)
from bkrobust.demo.evaluate import is_valid_adjustment_set_dag

WELLPOSED = "--wellposed" in sys.argv
ALPHAS = [0.0, 0.01, 0.05, 0.10, 0.20, 0.25, 0.33, 0.40, 0.45, 0.49,
          0.50, 0.51, 0.55, 0.60, 0.75, 0.90]

def analyse(pr):
    space = enumerate_space(pr["cpdag"])
    reps = represented_dags(space)
    covers = covering_pairs(space, reps)
    nbrs = neighbour_graph(space, covers)
    g0 = next((g for g in space if g == pr["g0"]), None)
    if g0 is None:
        return None
    dist = bfs_distances(nbrs, g0)
    x, y, z = pr["x"], pr["y"], pr["z"]
    ext0 = reps[g0]
    beta0 = sum(1 for dd in ext0 if not is_valid_adjustment_set_dag(dd, x, y, z))/len(ext0)
    rows = []
    for g in space:
        d = dist.get(g)
        if d is None:
            continue
        ext = reps[g]
        bad = sum(1 for dd in ext if not is_valid_adjustment_set_dag(dd, x, y, z))
        beta = bad / len(ext)
        if g == g0:
            direction = "self"
        elif ext < ext0:
            direction = "down"      # refinement: strictly fewer worlds
        elif ext0 < ext:
            direction = "up"        # coarsening / retraction
        else:
            direction = "side"
        rows.append(dict(d=d, beta=beta, n=len(ext), dir=direction))
    return rows, beta0

def r_alpha(rows, a):
    cand = [r["d"] for r in rows if r["d"] > 0 and r["beta"] > a]
    return min(cand) if cand else None   # None = infinity

def main(n_problems=250, seed=20260907):
    rng = np.random.default_rng(seed)
    t0 = time.time()
    got = tried = illposed = 0
    per_alpha = {a: [] for a in ALPHAS}
    beta_all, beta_d1, down_beta_pos, side_count = [], [], 0, 0
    ball_sizes = []          # |{G : 0 < dist < r^alpha}| at a few alphas
    hedge = {a: [] for a in ALPHAS}
    knife = 0                # problems with a beta in [0.45,0.55] at dist 1
    space_sizes = []
    while got < n_problems and tried < 40000:
        tried += 1
        pr = make_problem(rng)
        if pr is None:
            continue
        out = analyse(pr)
        if out is None:
            continue
        rows, beta0 = out
        if WELLPOSED and beta0 > 0:
            continue
        if beta0 > 0:
            illposed += 1
        got += 1
        space_sizes.append(len(rows))
        for a in ALPHAS:
            r = r_alpha(rows, a)
            per_alpha[a].append(r)
            inside = [q for q in rows if q["d"] > 0 and (r is None or q["d"] < r)]
            hedge[a].append(len(inside))
        for q in rows:
            if q["dir"] != "self":
                beta_all.append(q["beta"])
            if q["d"] == 1:
                beta_d1.append(q["beta"])
            if q["dir"] == "down" and q["beta"] > 0:
                down_beta_pos += 1
            if q["dir"] == "side":
                side_count += 1
        if any(0.45 <= q["beta"] <= 0.55 for q in rows if q["d"] == 1):
            knife += 1
    el = time.time() - t0
    print(f"WELLPOSED_ONLY={WELLPOSED}  ill-posed (beta(G0)>0) kept in sample: {illposed}")
    print(f"problems={got} (from {tried} draws)  space size mean={np.mean(space_sizes):.1f} max={max(space_sizes)}  {el:.1f}s")
    print()
    print(f"{'alpha':>6} {'frac r=1':>9} {'frac r=2':>9} {'frac inf':>9} {'mean(finite)':>13} {'mean(inf->cap)':>15} {'mean hedge |ball|':>18}")
    for a in ALPHAS:
        rs = per_alpha[a]
        fin = [r for r in rs if r is not None]
        ninf = sum(1 for r in rs if r is None)
        cap = max(fin) + 1 if fin else 1
        capped = [r if r is not None else cap for r in rs]
        print(f"{a:>6.2f} {sum(1 for r in rs if r==1)/len(rs):>9.3f} {sum(1 for r in rs if r==2)/len(rs):>9.3f}"
              f" {ninf/len(rs):>9.3f} {np.mean(fin) if fin else float('nan'):>13.3f}"
              f" {np.mean(capped):>15.3f} {np.mean(hedge[a]):>18.2f}")
    print()
    b = np.array(beta_all)
    print(f"beta over all non-G0 elements (n={len(b)}): ==0: {np.mean(b==0):.3f}   in (0,0.25): {np.mean((b>0)&(b<0.25)):.3f}"
          f"   [0.25,0.5): {np.mean((b>=0.25)&(b<0.5)):.3f}   ==0.5: {np.mean(np.isclose(b,0.5)):.3f}"
          f"   (0.5,1): {np.mean((b>0.5)&(b<1)):.3f}   ==1: {np.mean(b==1):.3f}")
    b1 = np.array(beta_d1)
    print(f"beta at hop 1 (n={len(b1)}):            ==0: {np.mean(b1==0):.3f}   in (0,0.25): {np.mean((b1>0)&(b1<0.25)):.3f}"
          f"   [0.25,0.5): {np.mean((b1>=0.25)&(b1<0.5)):.3f}   ==0.5: {np.mean(np.isclose(b1,0.5)):.3f}"
          f"   (0.5,1): {np.mean((b1>0.5)&(b1<1)):.3f}   ==1: {np.mean(b1==1):.3f}")
    print(f"distinct beta values at hop 1: {sorted(set(round(v,4) for v in beta_d1))[:25]}")
    print(f"downward elements with beta>0: {down_beta_pos}   incomparable ('side') elements: {side_count}")
    print(f"problems with a hop-1 beta in [0.45,0.55] (knife-edge at alpha=0.5): {knife}/{got} = {knife/got:.3f}")
    # how often does a small alpha change the answer vs alpha=0
    r0 = per_alpha[0.0]
    for a in [0.05, 0.10, 0.25, 0.40, 0.49, 0.50]:
        ra = per_alpha[a]
        ch = sum(1 for u, v in zip(r0, ra) if u != v)
        print(f"  alpha={a:<5} changes r vs alpha=0 on {ch}/{got} = {ch/got:.3f} problems")

if __name__ == "__main__":
    args=[a for a in sys.argv[1:] if not a.startswith("--")]
    main(int(args[0]) if args else 250)
