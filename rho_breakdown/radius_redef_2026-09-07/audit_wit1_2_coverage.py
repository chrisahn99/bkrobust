#!/usr/bin/env python3
"""witness-1 audit 2: coverage and width of the L_flip hedge under budget-1 and budget-2.

The analyst asserts K_wrong = K_true with j claims reversed. Z = O*(cl(K_wrong)).
Hedge sets are built FROM K_wrong (the only thing the analyst has).
Interval = [min, max] of the adjusted estimand over {Z} u {O*(g) : g in hedge}.
"""
import sys
from collections import Counter, defaultdict
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


def elicit(rng, dag, cpdag, k_claims):
    und = sorted(tuple(sorted(e)) for e in cpdag.undirected_edges)
    if not und:
        return None
    idx = rng.permutation(len(und))[:min(k_claims, len(und))]
    K = []
    for t in idx:
        a, b = und[t]
        K.append((a, b) if (a, b) in dag.directed_edges else (b, a))
    return K


def partition(cpdag, K, x, y, z):
    R, Lf, Lr = [], [], []
    grev, gret = {}, {}
    for k in K:
        rest = [e for e in K if e != k]
        gr = apply_orientations(cpdag, rest + [(k[1], k[0])])
        if gr is None:
            R.append(k)
        else:
            grev[k] = gr
            if not is_valid_adjustment_set_mpdag(gr, x, y, z):
                Lf.append(k)
        gt = apply_orientations(cpdag, rest)
        if gt is not None:
            gret[k] = gt
            if not is_valid_adjustment_set_mpdag(gt, x, y, z):
                Lr.append(k)
    return R, Lf, Lr, grev, gret


def estimands(sem, x, y, graphs):
    """Adjusted estimand for O*(g) of each candidate graph; fall back to DAG-level."""
    vals = []
    for g in graphs:
        zz = optimal_adjustment_set_mpdag(g, x, y)
        if zz is not None:
            vals.append(adjusted_estimand(sem, x, y, zz))
        else:
            for d in enumerate_dag_extensions(g):
                vals.append(adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(d, x, y)))
    return vals


def run(budget, n_problems=300, seed=20260907):
    rng = np.random.default_rng(seed)
    rows = []
    tried = 0
    while len(rows) < n_problems and tried < 60000:
        tried += 1
        dag = random_dag(rng, 6, 0.35)
        if not dag.directed_edges:
            continue
        cpdag = dag_to_cpdag(dag)
        und = sorted(tuple(sorted(e)) for e in cpdag.undirected_edges)
        if not (2 <= len(und) <= 6):
            continue
        K_true = elicit(rng, dag, cpdag, 3)
        if K_true is None or len(K_true) < budget:
            continue
        flip_idx = rng.permutation(len(K_true))[:budget]
        K_w = [((k[1], k[0]) if i in flip_idx else k) for i, k in enumerate(K_true)]
        g0 = apply_orientations(cpdag, K_w)
        if g0 is None:
            continue
        nodes = sorted(dag.nodes)
        cand = [(a, b) for a in nodes for b in nodes if a != b]
        rng.shuffle(cand)
        picked = None
        for x, y in cand:
            z = optimal_adjustment_set_mpdag(g0, x, y)
            if z is None:
                continue
            picked = (x, y, frozenset(z))
            break
        if picked is None:
            continue
        x, y, z = picked
        sem = random_sem(dag, rng)
        truth = adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(dag, x, y))
        point = adjusted_estimand(sem, x, y, z)

        R, Lf, Lr, grev, gret = partition(cpdag, K_w, x, y, z)
        hedges = {
            "none": [],
            "L_flip": [grev[k] for k in Lf],
            "L_ret": [gret[k] for k in Lr],
            "all-flip (no screen)": [grev[k] for k in K_w if k in grev],
            "class (blanket)": [cpdag],
        }
        rec = dict(truth=truth, point=point, z_valid=is_valid_adjustment_set_dag(dag, x, y, z),
                   nR=len(R), nLf=len(Lf), nLr=len(Lr), m=len(K_w),
                   false_in_R=any(K_true[i] != K_w[i] and K_w[i] in R for i in range(len(K_w))),
                   false_in_Lf=any(K_true[i] != K_w[i] and K_w[i] in Lf for i in range(len(K_w))))
        for name, gs in hedges.items():
            vals = [point] + estimands(sem, x, y, gs)
            lo, hi = min(vals), max(vals)
            rec[name] = (lo, hi, hi - lo, lo <= truth <= hi)
        rows.append(rec)
    return rows


def report(rows, budget):
    print(f"===== BUDGET {budget}: {len(rows)} problems =====")
    zv = np.mean([r["z_valid"] for r in rows])
    print(f"Z actually valid in the true DAG: {100*zv:.1f}%  (so nominal point estimate "
          f"is already right on that share)")
    print(f"mean |K|={np.mean([r['m'] for r in rows]):.2f}  |R|={np.mean([r['nR'] for r in rows]):.2f}  "
          f"|L_flip|={np.mean([r['nLf'] for r in rows]):.2f}  |L_ret|={np.mean([r['nLr'] for r in rows]):.2f}")
    print(f"the FALSE claim landed in R      : {100*np.mean([r['false_in_R'] for r in rows]):.1f}%  "
          f"(R is sound only if this is 0 at budget 1)")
    print(f"the FALSE claim landed in L_flip : {100*np.mean([r['false_in_Lf'] for r in rows]):.1f}%")
    print()
    print(f"{'method':24s} {'raw cov':>8s} {'mean width':>11s} {'scale->95%':>11s} {'calib. width':>13s}")
    for name in ["none", "L_flip", "L_ret", "all-flip (no screen)", "class (blanket)"]:
        cov = np.mean([r[name][3] for r in rows])
        w = np.array([r[name][2] for r in rows])
        # calibrate: scale every interval about its centre until 95% coverage
        errs = []
        for r in rows:
            lo, hi, _, _ = r[name]
            c = 0.5 * (lo + hi)
            half = 0.5 * (hi - lo)
            d = abs(r["truth"] - c)
            errs.append(np.inf if half <= 1e-12 and d > 1e-12 else (d / half if half > 1e-12 else 0.0))
        errs = np.array(errs)
        s = np.quantile(errs[np.isfinite(errs)], 0.95) if np.isfinite(errs).mean() >= 0.95 else np.inf
        if np.isfinite(s):
            cw = np.mean(w * s)
            print(f"{name:24s} {100*cov:7.1f}% {w.mean():11.3f} {s:11.2f} {cw:13.3f}")
        else:
            frac = 100*np.mean(np.isinf(errs))
            print(f"{name:24s} {100*cov:7.1f}% {w.mean():11.3f} {'INF':>11s} "
                  f"{'UNREACHABLE':>13s}   ({frac:.1f}% zero-width misses)")


if __name__ == "__main__":
    for b in (1, 2):
        report(run(b), b)
        print()
