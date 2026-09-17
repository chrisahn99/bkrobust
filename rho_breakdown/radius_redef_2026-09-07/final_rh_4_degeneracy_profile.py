#!/usr/bin/env python3
"""FINAL DEFINITION, part 4 -- the degeneracy profile, which is the thing to predict.

nu := w(R_{t*}) / w(Theta(Chat)) in [0,1].  Degenerate = all mass at 0 (the hedge is
a point) or at 1 (the hedge IS the blanket).  Non-degenerate = interior mass.
Hypothesis: interior mass is governed by f := |undirected edges of Chat| - |K|, the
ambiguity the claims do NOT touch.  Measured here as a function of f.
"""
import sys
import numpy as np
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag, random_sem
from harness import random_dag
from final_rh_1_ladder import Theta, hull, w, merge, ladder, TOL


def make_f(rng, n, p, nk, und_lo, und_hi):
    dag = random_dag(rng, n, p)
    if not dag.directed_edges:
        return None
    cp = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cp.undirected_edges)
    if not (und_lo <= len(und) <= und_hi) or len(und) < nk:
        return None
    nodes = sorted(dag.nodes)
    cand = [(x, y) for x in nodes for y in nodes if x != y]
    rng.shuffle(cand)
    for x, y in cand:
        idx = rng.permutation(len(und))[:nk]
        K = [((a, b) if (a, b) in dag.directed_edges else (b, a))
             for a, b in (und[t] for t in idx)]
        g0 = apply_orientations(cp, K)
        if g0 is None:
            continue
        if optimal_adjustment_set_mpdag(g0, x, y) is None:
            continue
        sem = random_sem(dag, rng)
        cache = {}
        blanket = hull(Theta(sem, cp, x, y, cache))
        if w(blanket) <= TOL:
            return None
        return dict(cpdag=cp, x=x, y=y, K=[tuple(e) for e in K], sem=sem,
                    blanket=blanket, cache=cache, f=len(und) - nk,
                    n_ext_blanket=len(enumerate_dag_extensions(cp)))
    return None


def run(seed, N, n, p, nk, und_lo, und_hi, q=0.20, tag=""):
    rng = np.random.default_rng(seed)
    rows, tried = [], 0
    while len(rows) < N and tried < 900000:
        tried += 1
        pr = make_f(rng, n, p, nk, und_lo, und_hi)
        if pr is None:
            continue
        K = pr["K"]
        flip = [bool(rng.random() < q) for _ in K]
        Kst = [((e[1], e[0]) if f else e) for e, f in zip(K, flip)]
        g0 = apply_orientations(pr["cpdag"], Kst)
        if g0 is None:
            continue
        if optimal_adjustment_set_mpdag(g0, pr["x"], pr["y"]) is None:
            continue
        base, entries, ncl = ladder(pr, Kst, g0, 2)
        if not entries:
            continue
        rows.append(dict(pr=pr, base=base, entries=entries, m=len(Kst),
                         b=sum(flip), t=pr["sem"].true_total_effect(pr["x"], pr["y"])))
    # calibrate on the ladder: t = 0, then top-j singletons, then all pairs
    grid = [("t=0", lambda r: r["base"])]
    for j in (1, 2, 3, 4):
        grid.append((f"top-{j}", lambda r, j=j: _take(r, j, 1)))
    grid.append(("all size<=1", lambda r: _take(r, 99, 1)))
    grid.append(("all size<=2", lambda r: _take(r, 99, 2)))
    star = None
    line = []
    for name, fn in grid:
        ivs = [fn(r) for r in rows]
        cov = np.mean([iv[0] - 1e-8 <= r["t"] <= iv[1] + 1e-8 for iv, r in zip(ivs, rows)])
        wid = np.mean([w(iv) for iv in ivs])
        line.append(f"{name}:{cov:.2f}/{wid:.3f}")
        if star is None and cov >= 0.95:
            star = (name, cov, wid, ivs)
    bl = np.mean([w(r["pr"]["blanket"]) for r in rows])
    print(f"\n{tag}  n={len(rows)} f={np.mean([r['pr']['f'] for r in rows]):.1f} "
          f"|K|={nk} blanket={bl:.4f} ext={np.mean([r['pr']['n_ext_blanket'] for r in rows]):.1f}")
    print("   ladder: " + "  ".join(line))
    if star is None:
        print("   NO ladder position reaches 0.95")
        return
    name, cov, wid, ivs = star
    nu = np.array([w(iv) / w(r["pr"]["blanket"]) for iv, r in zip(ivs, rows)])
    print(f"   t* = {name}: cov={cov:.3f} width={wid:.4f} = {100*wid/bl:.1f}% of blanket")
    print(f"   nu: mean={nu.mean():.3f}  at0={100*np.mean(nu<1e-6):.1f}%  "
          f"at1={100*np.mean(nu>1-1e-6):.1f}%  INTERIOR={100*np.mean((nu>=1e-6)&(nu<=1-1e-6)):.1f}%")


def _take(r, j, maxsize):
    iv, seen = r["base"], 0
    for S, e_iv, lam, size in r["entries"]:
        if size > maxsize:
            break
        if size < maxsize or seen < j:
            iv = merge(iv, e_iv)
            if size == maxsize:
                seen += 1
    return iv


if __name__ == "__main__":
    run(20260907, 200, 9, 0.30, 3, 3, 4, tag="[f small]  9 nodes, |und| 3-4, |K|=3")
    run(20260907, 200, 12, 0.22, 3, 6, 8, tag="[f large] 12 nodes, |und| 6-8, |K|=3")
    run(20260907, 200, 12, 0.22, 4, 7, 9, tag="[f large] 12 nodes, |und| 7-9, |K|=4")
