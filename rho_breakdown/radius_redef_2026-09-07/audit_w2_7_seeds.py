#!/usr/bin/env python3
"""Audit 7: robustness of the two kill findings at larger n, and the
'zero-width and wrong' rate that makes the target metric unreachable."""
import sys
import numpy as np
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
exec(open("/private/tmp/claude-501/-Users-josecosta-mugango/1cbd3bff-6ceb-4991-8a2a-b20d5a3f01d3/scratchpad/audit_lambda4.py").read().split("rng = np.random.default_rng(20260907)")[0])

for seed in (20260907, 4242, 99991):
    rng = np.random.default_rng(seed)
    P, tried = [], 0
    while len(P) < 400 and tried < 900000:
        tried += 1
        pr = make(rng)
        if pr is not None: P.append(pr)
    out = []
    for b in (1, 2, 3):
        rng2 = np.random.default_rng(seed + 1)
        dead = used = cov = 0; ws = []
        for pr in P:
            cp, sem, x, y, K = pr["cpdag"], pr["sem"], pr["x"], pr["y"], pr["K"]
            if len(K) < b: continue
            pick = set(rng2.permutation(len(K))[:b].tolist())
            Kst = [((e[1], e[0]) if i in pick else e) for i, e in enumerate(K)]
            g0 = apply_orientations(cp, Kst)
            if g0 is None: continue
            if optimal_adjustment_set_mpdag(g0, x, y) is None: continue
            th0 = theta_set(sem, g0, x, y)
            if not th0: continue
            u = list(th0)
            for k in Kst:
                for g in fibre(cp, Kst, k): u.extend(theta_set(sem, g, x, y))
            lo, hi = min(u), max(u); t = sem.true_total_effect(x, y)
            used += 1; ws.append(hi - lo)
            ok = lo - 1e-8 <= t <= hi + 1e-8
            if ok: cov += 1
            if hi - lo <= TOL and not ok: dead += 1
        out.append(f"b={b}: n={used:>3} cov={cov/max(used,1):.3f} "
                   f"width={np.mean(ws):.4f} zero-width-and-wrong={100*dead/max(used,1):.1f}%")
    print(f"seed {seed}: {len(P)} problems | " + " | ".join(out))
