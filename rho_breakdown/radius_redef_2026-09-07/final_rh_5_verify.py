#!/usr/bin/env python3
"""FINAL DEFINITION, part 5 -- the three structural checks and a seed sweep.

(i)  R at full budget IS the blanket hedge  =>  the family always calibrates.
(ii) cost at the calibrated position: distinct DAGs touched vs |ext(Chat)|.
(iii) headline saving across three seeds and three elicitation error rates.
"""
import sys, itertools
import numpy as np
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag
from final_rh_1_ladder import make, Theta, hull, w, merge, ladder, TOL


def cases(seed, N, q, kmax=4):
    rng = np.random.default_rng(seed)
    P, tried = [], 0
    while len(P) < N and tried < 600000:
        tried += 1
        pr = make(rng, kmax=kmax)
        if pr is None or w(pr["blanket"]) <= TOL:
            continue
        P.append(pr)
    out = []
    for pr in P:
        K = pr["K"]
        flip = [bool(rng.random() < q) for _ in K]
        Kst = [((e[1], e[0]) if f else e) for e, f in zip(K, flip)]
        g0 = apply_orientations(pr["cpdag"], Kst)
        if g0 is None:
            continue
        if optimal_adjustment_set_mpdag(g0, pr["x"], pr["y"]) is None:
            continue
        base, entries, ncl = ladder(pr, Kst, g0, len(Kst))
        out.append(dict(pr=pr, base=base, entries=entries, m=len(Kst), ncl=ncl,
                        b=sum(flip), t=pr["sem"].true_total_effect(pr["x"], pr["y"]),
                        Kst=Kst, g0=g0))
    return out


def hull_upto(c, maxsize, jtop=None):
    iv, seen = c["base"], 0
    for S, e_iv, lam, size in c["entries"]:
        if size < maxsize:
            iv = merge(iv, e_iv)
        elif size == maxsize:
            if jtop is None or seen < jtop:
                iv = merge(iv, e_iv)
                seen += 1
    return iv


if __name__ == "__main__":
    C = cases(20260907, 200, 0.20)
    # (i) full budget == blanket
    bad = 0
    for c in C:
        full = hull_upto(c, c["m"])
        bl = c["pr"]["blanket"]
        if abs(full[0] - bl[0]) > 1e-9 or abs(full[1] - bl[1]) > 1e-9:
            bad += 1
    print(f"(i)  R at full budget b=|K| equals Theta(Chat) exactly on {len(C)-bad}/{len(C)} "
          f"problems  (mismatches: {bad})")

    # (ii) cost at the calibrated position (base + every singleton, which the
    #      lambda-ordering needs anyway) vs the blanket's extension count
    ratios, over = [], 0
    for c in C:
        seen = set()
        for g in [c["g0"]] + [apply_orientations(c["pr"]["cpdag"],
                                                 [e for i, e in enumerate(c["Kst"]) if i != k])
                              for k in range(c["m"])]:
            if g is None:
                continue
            for d in enumerate_dag_extensions(g):
                seen.add(d.key())
        nb = c["pr"]["n_ext_blanket"]
        ratios.append(len(seen) / nb)
        over += len(seen) > nb
    print(f"(ii) distinct theta_D evaluations at t*<=1, as a fraction of the blanket's: "
          f"mean={np.mean(ratios):.3f} max={np.max(ratios):.3f}; exceeds the blanket on "
          f"{over}/{len(C)} problems")
    print(f"     Meek closures for the whole singleton tier: mean={np.mean([c['m'] for c in C]):.1f}")

    # (iii) headline across seeds and error rates
    print("\n(iii) seed x q sweep -- calibrated width vs the blanket hedge")
    for seed in (20260907, 4242, 99991):
        for q in (0.05, 0.10, 0.20, 0.35):
            C = cases(seed, 200, q)
            bl = np.mean([w(c["pr"]["blanket"]) for c in C])
            best = None
            for name, fn in [("t=0", lambda c: c["base"]),
                             ("top-1", lambda c: hull_upto(c, 1, 1)),
                             ("top-2", lambda c: hull_upto(c, 1, 2)),
                             ("all-1", lambda c: hull_upto(c, 1)),
                             ("all-2", lambda c: hull_upto(c, 2)),
                             ("all-3", lambda c: hull_upto(c, 3))]:
                ivs = [fn(c) for c in C]
                cov = np.mean([iv[0] - 1e-8 <= c["t"] <= iv[1] + 1e-8 for iv, c in zip(ivs, C)])
                wid = np.mean([w(iv) for iv in ivs])
                if best is None and cov >= 0.95:
                    best = (name, cov, wid, ivs)
            if best is None:
                print(f"  seed {seed} q={q}: no position reaches 0.95")
                continue
            name, cov, wid, ivs = best
            nu = np.array([w(iv) / w(c["pr"]["blanket"]) for iv, c in zip(ivs, C)])
            print(f"  seed {seed} q={q:<5} n={len(C):>3} t*={name:<6} cov={cov:.3f} "
                  f"width={wid:.4f} blanket={bl:.4f} saving={100*(1-wid/bl):5.1f}%  "
                  f"nu at0={100*np.mean(nu<1e-6):4.1f}% at1={100*np.mean(nu>1-1e-6):4.1f}% "
                  f"int={100*np.mean((nu>=1e-6)&(nu<=1-1e-6)):4.1f}%")
