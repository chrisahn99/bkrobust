"""weighted-5 audit 3: (a) is it r_val renamed? (b) redundancy gaming (c) can the
target metric even calibrate it? (d) does the soundness bug repair?"""
import sys, math
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from collections import Counter
from bkrobust.demo.evaluate import (optimal_adjustment_set_dag, is_valid_adjustment_set_dag,
                                    adjusted_estimand, asymptotic_variance, random_sem)
from audit_w0_rk import make_problem, m_of, r_K, r_val, gen
INF = math.inf

def parents(d, x):
    return {a for (a, b) in d.directed_edges if b == x}

def beta_ostar(sem, d, x, y):
    return adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(d, x, y))

def beta_pa(sem, d, x, y):
    """total effect implied by d: adjust for pa_d(X) -- ALWAYS valid in d."""
    return adjusted_estimand(sem, x, y, parents(d, x))

def rk_delta(rec, beta_Z, delta, key):
    best = INF
    for r in rec:
        if abs(r[key] - beta_Z) > delta:
            best = min(best, r['m'])
    return best

def hull(rec, beta_Z, t, key):
    v = [r[key] for r in rec if r['m'] <= t]
    return (min(v), max(v)) if v else (beta_Z, beta_Z)

NPROB = int(sys.argv[1]) if len(sys.argv) > 1 else 500
for nf in (1, 2):
    rng = np.random.default_rng(4242 + nf)
    probs, tried = gen(rng, NPROB, n=7, p=0.40, k_claims=3, n_false=nf)
    rng2 = np.random.default_rng(31337 + nf)
    P = []
    for pr in probs:
        dag, x, y, Z, K = pr['dag'], pr['x'], pr['y'], pr['z'], pr['K']
        sem = random_sem(dag, rng2)
        bZ = adjusted_estimand(sem, x, y, Z)
        tau = sem.true_total_effect(x, y)
        rec = [dict(m=m_of(d, K), bo=beta_ostar(sem, d, x, y), bp=beta_pa(sem, d, x, y),
                    ok=is_valid_adjustment_set_dag(d, x, y, Z)) for d in pr['fibre']]
        rv, _ = r_val(pr['cpdag'], pr['g0'], x, y, Z)
        rk = r_K(pr['fibre'], K, x, y, Z)
        Od = optimal_adjustment_set_dag(dag, x, y)
        P.append(dict(rec=rec, bZ=bZ, tau=tau, rv=rv, rk=rk, nK=len(K), pr=pr, sem=sem,
                      amen=is_valid_adjustment_set_dag(dag, x, y, Od),
                      av=asymptotic_variance(sem, x, y, Z),
                      spread=max(abs(r['bo'] - bZ) for r in rec)))
    live = [p for p in P if p['spread'] > 1e-9]
    n = len(live)
    print(f"\n{'#'*88}\nf={nf}  {len(P)} problems ({tried} draws), LIVE {n} ({100*n/len(P):.0f}%), "
          f"amenable {100*np.mean([p['amen'] for p in P]):.0f}%")

    # ---------- (a) IS IT r_val RENAMED? same problems, three statistics ----------
    print("\n(a) THREE STATISTICS ON THE SAME PROBLEMS (live only), P(stat = 1):")
    d05 = 0.05
    rows = [("r_val   (hop, validity)      ", [p['rv'] for p in live]),
            ("r_K     (contra, validity)   ", [p['rk'] for p in live]),
            ("r_K(0+) (contra, estimate)   ", [rk_delta(p['rec'], p['bZ'], 1e-12, 'bo') for p in live]),
            ("r_K(.05)(contra, estimate)   ", [rk_delta(p['rec'], p['bZ'], d05, 'bo') for p in live]),
            ("r_K(.25)                     ", [rk_delta(p['rec'], p['bZ'], .25, 'bo') for p in live])]
    for lab, v in rows:
        c = Counter('INF' if (x is INF or x is None) else x for x in v)
        print(f"  {lab} 1:{c.get(1,0)/n:.3f} 2:{c.get(2,0)/n:.3f} 3:{c.get(3,0)/n:.3f} "
              f"INF:{c.get('INF',0)/n:.3f}  |  P(in {{1,INF}})={(c.get(1,0)+c.get('INF',0))/n:.3f}")
    # agreement rate between r_K(delta) and r_val
    for dl in (1e-12, .05, .25):
        same = np.mean([ (rk_delta(p['rec'],p['bZ'],dl,'bo') == p['rv']) or
                         ((rk_delta(p['rec'],p['bZ'],dl,'bo') is INF) and p['rv'] is INF)
                         for p in live])
        samek = np.mean([ rk_delta(p['rec'],p['bZ'],dl,'bo') == p['rk'] for p in live])
        print(f"    delta={dl:<6.2g} r_K(delta)==r_val on {100*same:5.1f}% of live | "
              f"r_K(delta)==r_K(validity) on {100*samek:5.1f}%")

    # ---------- (b) REDUNDANCY GAMING: restate each claim c times ----------
    print("\n(b) GAMING BY RESTATEMENT: repeat every elicited claim c times (same G0, same Z)")
    for c in (1, 2, 3):
        vals = []
        for p in live:
            rec2 = [dict(m=c*r['m'], bo=r['bo']) for r in p['rec']]
            vals.append(rk_delta(rec2, p['bZ'], d05, 'bo'))
        cc = Counter('INF' if v is INF else v for v in vals)
        fin = [v for v in vals if v is not INF]
        print(f"    K restated x{c}: P(r_K(.05)=1)={cc.get(1,0)/n:.3f}  "
              f"median finite r_K(.05)={np.median(fin) if fin else float('nan'):.1f}  "
              f"mean finite={np.mean(fin) if fin else float('nan'):.2f}")

    # ---------- (c) CAN THE TARGET METRIC CALIBRATE IT? ----------
    print("\n(c) TARGET-METRIC CALIBRATION: scale each method's interval to 95% coverage of tau")
    NOBS = 2000
    for dl in (.05, .25, 1.0):
        # method: hull of S = {m < r_K(delta)} + sampling term
        half, err = [], []
        for p in live:
            r = rk_delta(p['rec'], p['bZ'], dl, 'bo')
            t = p['nK'] if r is INF else r - 1
            lo, hi = hull(p['rec'], p['bZ'], t, 'bo')
            se = math.sqrt(p['av'] / NOBS)
            centre = 0.5 * (lo + hi)
            half.append(0.5 * (hi - lo) + 1.96 * se)
            err.append(abs(p['tau'] - centre))
        half, err = np.array(half), np.array(err)
        ratio = err / np.maximum(half, 1e-300)
        cscale = np.quantile(ratio, 0.95)
        print(f"    delta={dl:<5.2g} raw coverage {np.mean(err<=half):.3f}  "
              f"mean raw width {2*np.mean(half):.4f}  -> scale c={cscale:8.2f}  "
              f"CALIBRATED mean width {2*cscale*np.mean(half):9.4f}")
    # blanket hedge over the whole class, same protocol
    half, err = [], []
    for p in live:
        lo, hi = hull(p['rec'], p['bZ'], p['nK'], 'bo')
        se = math.sqrt(p['av'] / NOBS)
        half.append(0.5*(hi-lo) + 1.96*se); err.append(abs(p['tau'] - 0.5*(lo+hi)))
    half, err = np.array(half), np.array(err)
    c2 = np.quantile(err/np.maximum(half,1e-300), .95)
    print(f"    BLANKET (whole class)  raw coverage {np.mean(err<=half):.3f}  "
          f"mean raw width {2*np.mean(half):.4f}  -> scale c={c2:8.2f}  "
          f"CALIBRATED mean width {2*c2*np.mean(half):9.4f}")
    # point-on-elicited
    half = np.array([1.96*math.sqrt(p['av']/NOBS) for p in live])
    err = np.array([abs(p['tau']-p['bZ']) for p in live])
    c3 = np.quantile(err/np.maximum(half,1e-300), .95)
    print(f"    POINT on elicited Z    raw coverage {np.mean(err<=half):.3f}  "
          f"mean raw width {2*np.mean(half):.4f}  -> scale c={c3:8.2f}  "
          f"CALIBRATED mean width {2*c3*np.mean(half):9.4f}")

    # ---------- (d) THE SOUNDNESS BUG AND ITS REPAIR ----------
    print("\n(d) SOUNDNESS: whenever m(d*)=n_false < r_K(delta), is |beta_Z - tau| <= delta?")
    for key, lab in (('bo', "beta(d)=O*(d)  [as defined]"), ('bp', "beta(d)=pa_d(X) [repaired]")):
        for dl in (.05, .25):
            tot = bad = 0; tota = bada = 0
            for p in live:
                r = rk_delta(p['rec'], p['bZ'], dl, key)
                if not (nf < (99 if r is INF else r)): continue
                ok = abs(p['bZ'] - p['tau']) <= dl + 1e-9
                tot += 1; bad += (not ok)
                if p['amen']: tota += 1; bada += (not ok)
            print(f"    {lab}  delta={dl:<5.2g} budget holds {tot:4d}x, VIOLATED {bad:3d} "
                  f"({100*bad/max(tot,1):5.1f}%) | amenable-only {bada}/{tota} "
                  f"({100*bada/max(tota,1):5.1f}%)")
