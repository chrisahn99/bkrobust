#!/usr/bin/env python3
"""Audit 5: THE TARGET METRIC. Calibrate each method to 95% empirical coverage by
scaling its interval about its centre; report the mean width it then needs."""
import sys
import numpy as np
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
exec(open("/private/tmp/claude-501/-Users-josecosta-mugango/1cbd3bff-6ceb-4991-8a2a-b20d5a3f01d3/scratchpad/audit_lambda4.py").read().split("rng = np.random.default_rng(20260907)")[0])

rng = np.random.default_rng(20260907)
P, tried = [], 0
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
while len(P) < N and tried < 800000:
    tried += 1
    pr = make(rng)
    if pr is not None: P.append(pr)
print(f"# {len(P)} problems (BK-relevant regime, theta_D corrected)\n")

def calibrated_width(rows, name):
    """rows = [(centre, halfwidth, truth)]. Scale s.t. 95% empirical coverage."""
    n = len(rows)
    dead = [1 for c, h, t in rows if h <= TOL and abs(t - c) > 1e-8]
    frac_dead = len(dead) / max(n, 1)
    if frac_dead > 0.05:
        print(f"  {name:<34} UNCALIBRATABLE: {100*frac_dead:.1f}% of cases have "
              f"zero width and miss the truth (>5%), so no scale reaches 95%")
        return None
    need = [0.0 if abs(t - c) <= 1e-9 else (abs(t - c) / h if h > TOL else np.inf)
            for c, h, t in rows]
    s = float(np.quantile(need, 0.95))
    w = float(np.mean([2 * s * h for c, h, t in rows]))
    cov = np.mean([1.0 if abs(t - c) <= s * h + 1e-9 else 0.0 for c, h, t in rows])
    print(f"  {name:<34} scale {s:>6.2f}   calibrated mean width {w:>8.4f}   "
          f"(achieved coverage {cov:.3f})")
    return w

for q in (0.0, 0.10, 0.20, 0.35):
    print(f"=== per-claim error rate q = {q:.2f} ===")
    rng2 = np.random.default_rng(2026)
    lam_rows, bl_rows, pt_rows, ret_rows = [], [], [], []
    for pr in P:
        cp, sem, x, y, K = pr["cpdag"], pr["sem"], pr["x"], pr["y"], pr["K"]
        Kst = [((e[1], e[0]) if rng2.random() < q else e) for e in K]
        g0 = apply_orientations(cp, Kst)
        if g0 is None: continue
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None: continue
        th0 = theta_set(sem, g0, x, y)
        if not th0: continue
        truth = sem.true_total_effect(x, y)
        # lambda / budget-1 fibre hull (delta = 0)
        u = list(th0)
        for k in Kst:
            for g in fibre(cp, Kst, k): u.extend(theta_set(sem, g, x, y))
        lam_rows.append((0.5*(min(u)+max(u)), 0.5*(max(u)-min(u)), truth))
        # blanket hedge over the whole equivalence class
        thC = pr["thC"]
        bl_rows.append((0.5*(min(thC)+max(thC)), 0.5*(max(thC)-min(thC)), truth))
        # point estimate on the elicited set (the "DoubleML on the elicited set" shape)
        pt = adjusted_estimand(sem, x, y, z)
        pt_rows.append((pt, 0.0, truth))
        # "retraction-only" hull: no reversal branch (a cheaper cousin)
        u2 = list(th0)
        for k in Kst:
            rest = [e for e in Kst if e != k]
            g = apply_orientations(cp, rest) if rest else cp
            if g is not None: u2.extend(theta_set(sem, g, x, y))
        ret_rows.append((0.5*(min(u2)+max(u2)), 0.5*(max(u2)-min(u2)), truth))
    print(f"  n = {len(lam_rows)}")
    calibrated_width(pt_rows,  "point estimate (elicited set)")
    calibrated_width(lam_rows, "lambda hedge, delta=0 (= budget-1 hull)")
    calibrated_width(ret_rows, "retraction-only hull (no reversals)")
    calibrated_width(bl_rows,  "blanket hedge (whole MEC)")
    print()
