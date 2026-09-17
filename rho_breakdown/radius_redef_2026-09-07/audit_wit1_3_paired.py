#!/usr/bin/env python3
"""witness-1 audit 3: does the L_flip screen ever BUY anything?

Restricted to LIVE problems (the hedge actually moves the estimand). Paired
comparison L_flip vs L_ret vs no-screen, on width and on coverage.
Also: does dropping L_ret\\L_flip ever drop the EXTREME candidate?
"""
import sys
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import random_dag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag,
                                    is_valid_adjustment_set_mpdag,
                                    is_valid_adjustment_set_dag,
                                    optimal_adjustment_set_dag,
                                    random_sem, adjusted_estimand)
from audit_wit1_2_coverage import elicit, partition, estimands


def run(budget, n_live=250, seed=13):
    rng = np.random.default_rng(seed)
    rows, tried = [], 0
    while len(rows) < n_live and tried < 200000:
        tried += 1
        dag = random_dag(rng, 7, 0.32)
        if not dag.directed_edges:
            continue
        cpdag = dag_to_cpdag(dag)
        und = sorted(tuple(sorted(e)) for e in cpdag.undirected_edges)
        if not (2 <= len(und) <= 6):
            continue
        K_true = elicit(rng, dag, cpdag, 3)
        if K_true is None or len(K_true) < budget:
            continue
        fi = set(rng.permutation(len(K_true))[:budget].tolist())
        K_w = [((k[1], k[0]) if i in fi else k) for i, k in enumerate(K_true)]
        g0 = apply_orientations(cpdag, K_w)
        if g0 is None:
            continue
        nodes = sorted(dag.nodes)
        cand = [(a, b) for a in nodes for b in nodes if a != b]
        rng.shuffle(cand)
        sem = random_sem(dag, rng)
        for x, y in cand:
            z = optimal_adjustment_set_mpdag(g0, x, y)
            if z is None:
                continue
            z = frozenset(z)
            truth = adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(dag, x, y))
            if abs(truth) < 0.05:
                continue
            point = adjusted_estimand(sem, x, y, z)
            R, Lf, Lr, grev, gret = partition(cpdag, K_w, x, y, z)
            v_flip = [point] + estimands(sem, x, y, [grev[k] for k in Lf])
            v_ret = [point] + estimands(sem, x, y, [gret[k] for k in Lr])
            v_allf = [point] + estimands(sem, x, y, [grev[k] for k in K_w if k in grev])
            v_cls = [point] + estimands(sem, x, y, [cpdag])
            w = lambda v: max(v) - min(v)
            # LIVE = the blanket hedge is not a point
            if w(v_cls) < 1e-8:
                continue
            rows.append(dict(
                truth=truth, point=point, rel=abs(point - truth) / max(abs(truth), 1e-9),
                z_valid=is_valid_adjustment_set_dag(dag, x, y, z),
                nR=len(R), nLf=len(Lf), nLr=len(Lr), m=len(K_w),
                strict=set(Lf) < set(Lr),
                false_in_R=any(K_true[i] != K_w[i] and K_w[i] in R for i in fi),
                false_in_Lf=any(K_true[i] != K_w[i] and K_w[i] in Lf for i in fi),
                w_flip=w(v_flip), w_ret=w(v_ret), w_allf=w(v_allf), w_cls=w(v_cls),
                c_flip=min(v_flip) <= truth <= max(v_flip),
                c_ret=min(v_ret) <= truth <= max(v_ret),
                c_allf=min(v_allf) <= truth <= max(v_allf),
                c_cls=min(v_cls) <= truth <= max(v_cls),
                c_none=abs(point - truth) < 1e-8))
            break
    return rows


def rep(rows, budget):
    n = len(rows)
    print(f"===== LIVE problems, BUDGET {budget}: n={n} =====")
    print(f"Z valid in truth: {100*np.mean([r['z_valid'] for r in rows]):.1f}%   "
          f"mean |K|={np.mean([r['m'] for r in rows]):.2f} |R|={np.mean([r['nR'] for r in rows]):.2f} "
          f"|L_flip|={np.mean([r['nLf'] for r in rows]):.2f} |L_ret|={np.mean([r['nLr'] for r in rows]):.2f}")
    print(f"L_flip STRICTLY smaller than L_ret on {100*np.mean([r['strict'] for r in rows]):.1f}% of problems")
    print(f"false claim in R: {100*np.mean([r['false_in_R'] for r in rows]):.1f}%   "
          f"false claim in L_flip: {100*np.mean([r['false_in_Lf'] for r in rows]):.1f}%")
    print()
    print(f"{'method':22s} {'coverage':>9s} {'mean width':>11s} {'median':>9s}")
    for nm, cw, ww in [("point (no hedge)", "c_none", None),
                       ("L_flip hedge", "c_flip", "w_flip"),
                       ("L_ret hedge", "c_ret", "w_ret"),
                       ("all-flip, no screen", "c_allf", "w_allf"),
                       ("class (blanket)", "c_cls", "w_cls")]:
        cov = 100*np.mean([r[cw] for r in rows])
        if ww is None:
            print(f"{nm:22s} {cov:8.1f}% {0.0:11.3f} {0.0:9.3f}")
        else:
            w = np.array([r[ww] for r in rows])
            print(f"{nm:22s} {cov:8.1f}% {w.mean():11.3f} {np.median(w):9.3f}")
    print()
    st = [r for r in rows if r["strict"]]
    print(f"--- ON THE {len(st)} PROBLEMS WHERE THE SCREEN ACTUALLY FIRES (L_flip strictly smaller) ---")
    if st:
        wf = np.array([r["w_flip"] for r in st]); wr = np.array([r["w_ret"] for r in st])
        print(f"  width L_flip {wf.mean():.4f}  vs  L_ret {wr.mean():.4f}   "
              f"saving {100*(1-wf.mean()/max(wr.mean(),1e-12)):.1f}%")
        print(f"  identical width on {100*np.mean(np.abs(wf-wr)<1e-9):.1f}% of them "
              f"(= the dropped candidate was NOT extreme, so the screen bought nothing)")
        print(f"  coverage L_flip {100*np.mean([r['c_flip'] for r in st]):.1f}%  "
              f"vs L_ret {100*np.mean([r['c_ret'] for r in st]):.1f}%")
    print()


if __name__ == "__main__":
    for b in (1, 2):
        rep(run(b), b)
