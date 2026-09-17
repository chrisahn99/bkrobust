#!/usr/bin/env python3
"""FINAL definition check v2: claim leverage lambda, its degeneracy profile, and
the width rule. Damage is measured as |theta_hat - tau_true|, NOT as the validity
predicate (which the v1 run showed fires on 66% of problems while the estimand
moves on 6.5%).
"""
import sys, itertools
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from final_leverage import problem, tau_range, calibrate, N_SAMPLE
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.evaluate import is_valid_adjustment_set_mpdag


def profile(p, budget_hedge=2):
    cp, K, x, y, z, sem, th = p["cpdag"], p["K"], p["x"], p["y"], p["z"], p["sem"], p["theta"]
    lam, inR, inLf, viol = {}, [], [], 0
    hull1 = [th]
    for k in K:
        rest = [e for e in K if e != k]
        gk = apply_orientations(cp, rest + [(k[1], k[0])])
        if gk is None:
            inR.append(k); lam[k] = 0.0; continue
        vals = tau_range(sem, x, y, gk)
        lam[k] = max(abs(v - th) for v in vals) if vals else 0.0
        hull1.extend(vals)
        if not is_valid_adjustment_set_mpdag(gk, x, y, z):
            inLf.append(k)
        elif lam[k] > 1e-9:
            viol += 1
    hull2 = list(hull1)
    for S in itertools.combinations(K, 2):
        rest = [e for e in K if e not in S]
        gS = apply_orientations(cp, rest + [(a[1], a[0]) for a in S])
        if gS is not None:
            hull2.extend(tau_range(sem, x, y, gS))
    return lam, inR, inLf, viol, (min(hull1), max(hull1)), (min(hull2), max(hull2))


def run(n_nodes, k_claims, budget, n, seed):
    rng = np.random.default_rng(seed)
    rows, tried = [], 0
    while len(rows) < n and tried < 500000:
        tried += 1
        p = problem(rng, n_nodes, k_claims, budget)
        if p is None:
            continue
        lam, inR, inLf, viol, (lo1, hi1), (lo2, hi2) = profile(p)
        cls = tau_range(p["sem"], p["x"], p["y"], p["cpdag"]) + [p["theta"]]
        v = sorted(lam.values(), reverse=True)
        rows.append(dict(theta=p["theta"], truth=p["truth"], se=p["se"], m=len(p["K"]),
                         nR=len(inR), nLf=len(inLf), viol=viol,
                         lmax=v[0] if v else 0.0, l2=v[1] if len(v) > 1 else 0.0,
                         npos=sum(1 for a in lam.values() if a > 1e-9),
                         lo1=lo1, hi1=hi1, lo2=lo2, hi2=hi2,
                         locls=min(cls), hicls=max(cls)))
    return rows


def report(rows, tag, budget):
    A = lambda k: np.array([r[k] for r in rows], dtype=float)
    n = len(rows)
    th, truth, se, lmax = A("theta"), A("truth"), A("se"), A("lmax")
    dmg = np.abs(th - truth)
    print(f"\n=========== {tag} | budget {budget} | n={n} ===========")
    print(f"mean |K|={A('m').mean():.2f}  |R|={A('nR').mean():.2f}  |L_flip|={A('nLf').mean():.2f}  "
          f"|{{k: lambda_k>0}}|={A('npos').mean():.2f}")
    print(f"(S) screen counterexamples (k not in L_flip, lambda_k>0): {int(A('viol').sum())}"
          f" / {int(A('m').sum())} claims")
    print(f"(B) lambda_max >= realised damage |theta-tau_true| violated on "
          f"{int(np.sum(lmax + 1e-12 < dmg))} / {n} problems")
    print(f"(Z) lambda_max == 0 on {100*np.mean(lmax<=1e-9):.1f}% ; of those, damage>0 on "
          f"{100*np.mean(dmg[lmax<=1e-9]>1e-9) if np.any(lmax<=1e-9) else 0:.1f}%")
    live = lmax > 1e-9
    if live.sum() >= 5:
        rel = lmax[live] / np.abs(th[live])
        q = np.round(rel, 2)
        vv, cc = np.unique(q, return_counts=True)
        print(f"(D) live stratum n={live.sum()} ({100*live.mean():.1f}%): largest atom of "
              f"lambda_max/|theta| at 2dp = {100*cc.max()/live.sum():.1f}% (value {vv[cc.argmax()]:.2f})")
        print(f"    p10={np.quantile(rel,.1):.3f} p50={np.quantile(rel,.5):.3f} "
              f"p90={np.quantile(rel,.9):.3f} max={rel.max():.2f}")
    rk = A("npos") >= 2
    print(f"(R) rankable (>=2 claims with lambda>0): {100*rk.mean():.1f}% of all, "
          f"{100*rk.sum()/max(live.sum(),1):.1f}% of live"
          + (f" ; top-1 beats top-2 by >10% on {100*np.mean((lmax>1.1*A('l2'))[rk]):.1f}%" if rk.sum() else ""))
    d1 = dmg > 1e-9
    if 0 < d1.sum() < n:
        r = np.empty(n)
        order = np.argsort(lmax); r[order] = np.arange(1, n+1)
        for v in np.unique(lmax):
            mm = lmax == v; r[mm] = r[mm].mean()
        n1, n0 = d1.sum(), n - d1.sum()
        print(f"(A) AUC of lambda_max for 'estimand actually damaged' "
              f"{(r[d1].sum()-n1*(n1+1)/2)/(n1*n0):.3f}  (prevalence {100*d1.mean():.1f}%)")
    for sub, nm2 in [(np.ones(n, bool), "ALL problems"),
                     (A("hicls") - A("locls") > 1e-9, "LIVE stratum (class hull > 0)")]:
        if sub.sum() < 5:
            continue
        print(f"(W) {nm2}: n={sub.sum()}   [lo-c*se, hi+c*se], c calibrated to 95%, N={N_SAMPLE}")
        print(f"    {'method':26s} {'c':>7s} {'cover':>7s} {'mean width':>11s} {'hull':>9s}")
        for nm, lo, hi in [("point (no hedge)", th, th),
                           ("flip hedge, budget 1", A("lo1"), A("hi1")),
                           ("flip hedge, budget 2", A("lo2"), A("hi2")),
                           ("blanket class hedge", A("locls"), A("hicls"))]:
            c, cov, w = calibrate(lo[sub], hi[sub], th[sub], se[sub], truth[sub])
            print(f"    {nm:26s} {c:7.2f} {100*cov:6.1f}% {w:11.4f} "
                  f"{np.mean((hi-lo)[sub]):9.4f}")


if __name__ == "__main__":
    for nn, kc, tag in [(7, 3, "7-node, K<=3"), (8, 5, "8-node, K<=5")]:
        for b in (1, 2):
            report(run(nn, kc, b, 500, 100 + b + nn), tag, b)
