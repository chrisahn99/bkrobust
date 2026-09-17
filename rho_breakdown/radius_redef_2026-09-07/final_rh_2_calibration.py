#!/usr/bin/env python3
"""FINAL DEFINITION, part 2 -- the fractional budget ladder, its calibration to
95% coverage, the width it needs, its degeneracy profile, and its cost.

The knob is t in [0, bmax]: admit every retraction set of size <= floor(t), plus
the top ceil(frac * C(m, floor(t)+1)) sets of size floor(t)+1 ordered by leverage
lambda.  R_t is nested and non-decreasing in t, so calibration is a line search.
"""
import sys, math, itertools
import numpy as np
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag
from final_rh_1_ladder import make, Theta, hull, w, merge, ladder, TOL

GRID = [0.0] + [b + f for b in range(0, 3) for f in (0.2, 0.4, 0.6, 0.8, 1.0)]


def R_at(base, entries, m_claims, t):
    """Hull at ladder position t.  entries are already sorted by (size, -lambda)."""
    b = int(math.floor(t + 1e-12))
    frac = t - b
    take_size = b + 1
    n_at = sum(1 for e in entries if e[3] == take_size)
    n_take = int(math.ceil(frac * n_at - 1e-12))
    iv, seen = base, 0
    for S, e_iv, lam, size in entries:
        if size <= b:
            iv = merge(iv, e_iv)
        elif size == take_size and seen < n_take:
            iv = merge(iv, e_iv)
            seen += 1
    return iv


def sweep(seed, N, q, bmax=3, kmax=4):
    rng = np.random.default_rng(seed)
    P, tried = [], 0
    while len(P) < N and tried < 500000:
        tried += 1
        pr = make(rng, kmax=kmax)
        if pr is None:
            continue
        if w(pr["blanket"]) <= TOL:
            continue                      # regime filter, stated in the paper
        P.append(pr)
    rng2 = np.random.default_rng(seed + 7)
    cases = []
    for pr in P:
        K = pr["K"]
        flip = [bool(rng2.random() < q) for _ in K]
        Kst = [((e[1], e[0]) if f else e) for e, f in zip(K, flip)]
        g0 = apply_orientations(pr["cpdag"], Kst)
        if g0 is None:
            continue                      # inconsistent -> the analyst sees it
        if optimal_adjustment_set_mpdag(g0, pr["x"], pr["y"]) is None:
            continue
        base, entries, ncl = ladder(pr, Kst, g0, min(bmax, len(Kst)))
        cases.append(dict(pr=pr, base=base, entries=entries, ncl=ncl,
                          m=len(Kst), b_true=sum(flip),
                          t=pr["sem"].true_total_effect(pr["x"], pr["y"])))
    return cases


def report(cases, q, label):
    n = len(cases)
    blanket_w = np.mean([w(c["pr"]["blanket"]) for c in cases])
    blanket_cov = np.mean([c["pr"]["blanket"][0] - 1e-8 <= c["t"] <= c["pr"]["blanket"][1] + 1e-8
                           for c in cases])
    print(f"\n--- {label}   q={q}  n={n}  mean|K|={np.mean([c['m'] for c in cases]):.2f} "
          f"  P(b_true>=1)={np.mean([c['b_true']>=1 for c in cases]):.2f}")
    print(f"    blanket hedge (Theta(Chat)) : cov={blanket_cov:.3f}  width={blanket_w:.4f}")
    tstar, curve = None, []
    for t in GRID:
        ivs = [R_at(c["base"], c["entries"], c["m"], t) for c in cases]
        cov = np.mean([iv[0] - 1e-8 <= c["t"] <= iv[1] + 1e-8 for iv, c in zip(ivs, cases)])
        wid = np.mean([w(iv) for iv in ivs])
        curve.append((t, cov, wid, ivs))
        if tstar is None and cov >= 0.95:
            tstar = (t, cov, wid, ivs)
    print("    ladder:  " + "  ".join(f"t={t:.1f}:{cov:.2f}/{wid:.3f}"
                                      for t, cov, wid, _ in curve[:11]))
    if tstar is None:
        print("    NO t reaches 0.95 -- fall back to the blanket")
        return
    t, cov, wid, ivs = tstar
    nu = [w(iv) / w(c["pr"]["blanket"]) for iv, c in zip(ivs, cases)]
    zero = np.mean([x < 1e-6 for x in nu])
    one = np.mean([x > 1 - 1e-6 for x in nu])
    print(f"    CALIBRATED t*={t:.1f}: cov={cov:.3f}  width={wid:.4f}  "
          f"= {100*wid/blanket_w:.1f}% of the blanket   (saving {100*(1-wid/blanket_w):.1f}%)")
    print(f"    DEGENERACY nu=w(R)/w(blanket): mean={np.mean(nu):.3f}  "
          f"at 0: {100*zero:.1f}%   at 1: {100*one:.1f}%   interior: {100*(1-zero-one):.1f}%")
    print(f"    zero-width-and-wrong: "
          f"{100*np.mean([w(iv)<1e-9 and not (iv[0]-1e-8<=c['t']<=iv[1]+1e-8) for iv,c in zip(ivs,cases)]):.1f}%")
    # locality split
    for name, sel in (("no claim touches X or Y", lambda c: not c["pr"]["touches"]),
                      ("some claim touches X/Y", lambda c: c["pr"]["touches"])):
        sub = [(iv, c) for iv, c in zip(ivs, cases) if sel(c)]
        if not sub:
            continue
        print(f"      [{name}] n={len(sub):>3} width={np.mean([w(i) for i,_ in sub]):.4f}  "
              f"nu={np.mean([w(i)/w(c['pr']['blanket']) for i,c in sub]):.3f}  "
              f"cov={np.mean([i[0]-1e-8<=c['t']<=i[1]+1e-8 for i,c in sub]):.3f}")
    # cost
    cost_ours = np.mean([len(c["pr"]["cache"]) for c in cases])
    cost_blanket = np.mean([c["pr"]["n_ext_blanket"] for c in cases])
    over = np.mean([len(c["pr"]["cache"]) > c["pr"]["n_ext_blanket"] for c in cases])
    print(f"    COST (distinct theta_D evaluations, DAG-cached): ours={cost_ours:.1f} "
          f"blanket={cost_blanket:.1f}  ours>blanket on {100*over:.1f}% of problems; "
          f"Meek closures={np.mean([c['ncl'] for c in cases]):.1f}")


if __name__ == "__main__":
    for q in (0.10, 0.20, 0.35):
        report(sweep(20260907, 300, q), q, "seed 20260907")
