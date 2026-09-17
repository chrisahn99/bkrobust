"""Audit of weighted-5: r_K(delta), the damage-thresholded contradiction number."""
import sys, math, time
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from collections import Counter
from bkrobust.demo.evaluate import (optimal_adjustment_set_dag, is_valid_adjustment_set_dag,
                                    adjusted_estimand, asymptotic_variance, random_sem)
from audit_w0_rk import make_problem, m_of, r_K, r_val, gen

INF = math.inf

def profile(pr, rng):
    """population beta(d) for every d in the fibre, plus m(d)."""
    dag, x, y, Z, K = pr['dag'], pr['x'], pr['y'], pr['z'], pr['K']
    sem = random_sem(dag, rng)
    beta_Z = adjusted_estimand(sem, x, y, Z)
    tau    = sem.true_total_effect(x, y)
    rec = []
    for d in pr['fibre']:
        O = optimal_adjustment_set_dag(d, x, y)
        rec.append((m_of(d, K), adjusted_estimand(sem, x, y, O), frozenset(O)))
    av = asymptotic_variance(sem, x, y, Z)       # n-free part
    # r_K : validity flavour (delta -> 0 of the *validity* criterion)
    rk_val = r_K(pr['fibre'], K, x, y, Z)
    return dict(rec=rec, beta_Z=beta_Z, tau=tau, av=av, rk_val=rk_val,
                nK=len(K), n_false=pr['n_false'], sem=sem, pr=pr)

def rk_delta(rec, beta_Z, delta):
    best = INF
    for m, b, _ in rec:
        if abs(b - beta_Z) > delta:
            best = min(best, m)
    return best

def hull(rec, beta_Z, t):
    v = [b for m, b, _ in rec if m <= t]
    if not v: return (beta_Z, beta_Z)
    return (min(v), max(v))

if __name__ == "__main__":
    NPROB = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    for nf in (0, 1):
        rng = np.random.default_rng(4242 + nf)
        probs, tried = gen(rng, NPROB, n=7, p=0.40, k_claims=3, n_false=nf)
        rng2 = np.random.default_rng(31337 + nf)
        P = [profile(pr, rng2) for pr in probs]
        print(f"\n{'='*84}\nf = {nf} deliberately false claim(s) in K   "
              f"| {len(P)} problems from {tried} draws | |K| = 3")

        spread = np.array([max(abs(b - p['beta_Z']) for _, b, _ in p['rec']) for p in P])
        bz = np.array([abs(p['beta_Z']) for p in P])
        print(f"observed spread max|beta(d)-beta_Z| : median {np.median(spread):.4f} "
              f"q90 {np.quantile(spread,.9):.4f} max {spread.max():.4f}")
        print(f"|beta_Z|                            : median {np.median(bz):.4f}")

        print("\n-- distribution of r_K(delta), and P(=1), against r_K (validity) --")
        c = Counter('INF' if p['rk_val'] is INF else p['rk_val'] for p in P)
        n = len(P)
        print(f"  r_K (validity)      1:{c.get(1,0)/n:.3f}  2:{c.get(2,0)/n:.3f} "
              f" 3:{c.get(3,0)/n:.3f}  INF:{c.get('INF',0)/n:.3f}")
        for lab, dl in [("delta=0+   ", 1e-12), ("delta=0.01 ", .01), ("delta=0.05 ", .05),
                        ("delta=0.10 ", .10), ("delta=0.25 ", .25), ("delta=0.50 ", .50),
                        ("delta=1.00 ", 1.0)]:
            vals = [rk_delta(p['rec'], p['beta_Z'], dl) for p in P]
            cc = Counter('INF' if v is INF else v for v in vals)
            print(f"  r_K({lab})     1:{cc.get(1,0)/n:.3f}  2:{cc.get(2,0)/n:.3f} "
                  f" 3:{cc.get(3,0)/n:.3f}  INF:{cc.get('INF',0)/n:.3f}")
        # relative delta
        for frac in (0.05, 0.10, 0.25):
            vals = [rk_delta(p['rec'], p['beta_Z'], frac*abs(p['beta_Z'])) for p in P]
            cc = Counter('INF' if v is INF else v for v in vals)
            print(f"  r_K(delta={frac:.2f}|beta_Z|) 1:{cc.get(1,0)/n:.3f}  2:{cc.get(2,0)/n:.3f} "
                  f" 3:{cc.get(3,0)/n:.3f}  INF:{cc.get('INF',0)/n:.3f}")

        print("\n-- delta tied to the standard error: 2*SE(n) --")
        for nobs in (200, 2000, 20000, 200000, 2000000):
            vals, ds = [], []
            for p in P:
                se = math.sqrt(p['av']/nobs); ds.append(2*se)
                vals.append(rk_delta(p['rec'], p['beta_Z'], 2*se))
            cc = Counter('INF' if v is INF else v for v in vals)
            print(f"  n={nobs:>8d}  median delta={np.median(ds):.4f}  "
                  f"1:{cc.get(1,0)/n:.3f}  2:{cc.get(2,0)/n:.3f}  3:{cc.get(3,0)/n:.3f} "
                  f" INF:{cc.get('INF',0)/n:.3f}")

        print("\n-- the advertised regime: fragile identification, stable estimate --")
        for dl in (.05, .10, .25):
            gap = sum(1 for p in P
                      if (p['rk_val'] is not INF)
                      and rk_delta(p['rec'], p['beta_Z'], dl) > p['rk_val'])
            strict = sum(1 for p in P if p['rk_val'] is not INF)
            print(f"  delta={dl:.2f}: r_K(delta) > r_K in {gap}/{strict} of the problems "
                  f"where r_K is finite = {100*gap/max(strict,1):.1f}%")

        print("\n-- hull width of S = {m < r_K(delta)} (population, no sampling term) --")
        for dl in (.05, .10, .25, .50):
            w, wb, cov = [], [], 0
            for p in P:
                r = rk_delta(p['rec'], p['beta_Z'], dl)
                t = (p['nK'] if r is INF else r - 1)
                lo, hi = hull(p['rec'], p['beta_Z'], t)
                w.append(hi - lo)
                blo, bhi = hull(p['rec'], p['beta_Z'], p['nK'])
                wb.append(bhi - blo)
                cov += (lo - 1e-12 <= p['tau'] <= hi + 1e-12)
            print(f"  delta={dl:.2f}  mean width {np.mean(w):.4f} (cap 2*delta={2*dl:.2f}) "
                  f"| blanket {np.mean(wb):.4f} | structural coverage of tau {cov/n:.3f}")
