"""weighted-5 audit 2: restrict to LIVE problems, test soundness + point-collapse."""
import sys, math
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from collections import Counter
from bkrobust.demo.evaluate import (optimal_adjustment_set_dag, is_valid_adjustment_set_dag,
                                    adjusted_estimand, asymptotic_variance, random_sem)
from audit_w0_rk import make_problem, m_of, r_K, gen
from audit_w5_rkdelta import profile, rk_delta, hull
INF = math.inf

def amenable(pr):
    d, x, y = pr['dag'], pr['x'], pr['y']
    O = optimal_adjustment_set_dag(d, x, y)
    return is_valid_adjustment_set_dag(d, x, y, O)

for nf in (0, 1, 2):
    rng = np.random.default_rng(4242 + nf)
    probs, tried = gen(rng, 600, n=7, p=0.40, k_claims=3, n_false=nf)
    rng2 = np.random.default_rng(31337 + nf)
    P = [profile(pr, rng2) for pr in probs]
    for p in P:
        p['spread'] = max(abs(b - p['beta_Z']) for _, b, _ in p['rec'])
        p['amen'] = amenable(p['pr'])
    live = [p for p in P if p['spread'] > 1e-9]
    print(f"\n{'='*86}\nf={nf}: {len(P)} problems, LIVE (fibre disagrees at all): "
          f"{len(live)} = {100*len(live)/len(P):.1f}%   amenable: "
          f"{100*np.mean([p['amen'] for p in P]):.1f}%")
    if not live: continue
    n = len(live)
    print(f"  on LIVE problems: median |beta_Z| {np.median([abs(p['beta_Z']) for p in live]):.4f}"
          f"  median spread {np.median([p['spread'] for p in live]):.4f}")
    print(f"\n  {'delta':>8s} {'P(r_K(d)=1)':>12s} {'P(hull is a POINT)':>19s} "
          f"{'mean hull w':>12s} {'2*delta':>8s} {'tau in hull':>12s} {'tau in hull|amen':>17s}"
          f" {'|bZ-tau|<=d':>12s}")
    for dl in (1e-12, .05, .10, .25, .50, 1.0, 2.0):
        deg = pt = 0; ws = []; cov = 0; cova = 0; na = 0; band = 0
        for p in live:
            r = rk_delta(p['rec'], p['beta_Z'], dl)
            deg += (r == 1)
            t = (p['nK'] if r is INF else r - 1)
            lo, hi = hull(p['rec'], p['beta_Z'], t)
            ws.append(hi - lo); pt += (hi - lo < 1e-12)
            ok = (lo - 1e-9 <= p['tau'] <= hi + 1e-9)
            cov += ok
            if p['amen']: na += 1; cova += ok
            band += (abs(p['beta_Z'] - p['tau']) <= dl)
        print(f"  {dl:8.2g} {deg/n:12.3f} {pt/n:19.3f} {np.mean(ws):12.4f} {2*dl:8.2f} "
              f"{cov/n:12.3f} {cova/max(na,1):17.3f} {band/n:12.3f}")

    # SOUNDNESS: whenever the budget holds (#false < r_K(delta)), is tau in the hull?
    print("\n  -- soundness test: problems where n_false < r_K(delta) --")
    for dl in (.05, .25, 1.0):
        tot = bad = tota = bada = 0
        for p in live:
            r = rk_delta(p['rec'], p['beta_Z'], dl)
            if not (nf < (99 if r is INF else r)):
                continue
            t = (p['nK'] if r is INF else r - 1)
            lo, hi = hull(p['rec'], p['beta_Z'], t)
            ok = (lo - 1e-9 <= p['tau'] <= hi + 1e-9)
            tot += 1; bad += (not ok)
            if p['amen']: tota += 1; bada += (not ok)
        print(f"    delta={dl:<5.2g} budget holds in {tot:4d} problems; "
              f"tau OUTSIDE the hull in {bad:3d} ({100*bad/max(tot,1):5.1f}%) "
              f"| amenable only: {bada}/{tota} ({100*bada/max(tota,1):5.1f}%)")

    # GAMING: what delta buys INF, and what it costs in honesty
    sp = np.array([p['spread'] for p in live])
    print(f"\n  -- gaming: quantiles of the observed spread (the ceiling delta must clear) --")
    print("     " + "  ".join(f"q{int(q*100)}={np.quantile(sp,q):.3f}" for q in (.5,.75,.9,.95,.99)))
    for dl in (.25, .5, 1.0, 2.0):
        inf = np.mean([rk_delta(p['rec'], p['beta_Z'], dl) is INF for p in live])
        miss = np.mean([abs(p['beta_Z']-p['tau']) > 0.10 for p in live])
        missinf = np.mean([(rk_delta(p['rec'], p['beta_Z'], dl) is INF)
                           and abs(p['beta_Z']-p['tau']) > 0.10 for p in live])
        print(f"     delta={dl:<4.2g} -> r_K(delta)=INF on {100*inf:5.1f}% of live problems; "
              f"of those, |beta_Z - tau| > 0.10 on {100*missinf/max(inf,1e-9):5.1f}% "
              f"(unconditional rate {100*miss:.1f}%)")
