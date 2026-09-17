#!/usr/bin/env python3
"""inverse-2 audit, part 3: distribution of phi_1, the data-refuted laundering rate,
and the phi_1 = 0 / pair-fatal blind spot, on random problems."""
import sys, itertools
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import make_problem, ball, r_val
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.evaluate import is_valid_adjustment_set_mpdag


def cells_of(cpdag, K, x, y, z, g0):
    out = {}
    for k in K:
        rest = [e for e in K if e != k]
        g_ret = apply_orientations(cpdag, rest)
        g_rev = apply_orientations(cpdag, rest + [(k[1], k[0])])
        adm = [g for g in (g_ret, g_rev) if g is not None]
        if g_rev is None:
            cell = "data-refuted"
        elif g_ret is not None and g_ret == g0:
            cell = "inert"
        elif any(not is_valid_adjustment_set_mpdag(g, x, y, z) for g in adm):
            cell = "load-bearing"
        else:
            cell = "free"
        out[k] = (cell, len(adm))
    lb = sum(1 for c, _ in out.values() if c == "load-bearing")
    den = sum(1 for _, n in out.values() if n >= 1)
    return out, (lb / den if den else None), lb, den


def nu2_fatal(cpdag, K, x, y, z):
    """Is there an admissible revision of a PAIR of claims that invalidates Z?"""
    for a, b in itertools.combinations(K, 2):
        rest = [e for e in K if e not in (a, b)]
        for ra in (None, (a[1], a[0])):
            for rb in (None, (b[1], b[0])):
                extra = [r for r in (ra, rb) if r is not None]
                g = apply_orientations(cpdag, rest + extra)
                if g is None:
                    continue
                if not is_valid_adjustment_set_mpdag(g, x, y, z):
                    return True
    return False


rng = np.random.default_rng(20260907)
phis, rvals = [], []
laundered_tot = laundered_bad = drefuted_tot = 0
phi0_and_pairfatal = phi0 = 0
n = 0
tried = 0
while n < 300 and tried < 40000:
    tried += 1
    pr = make_problem(rng, k_claims=4)
    if pr is None:
        continue
    cpdag, K, x, y, z, g0 = pr["cpdag"], pr["K"], pr["x"], pr["y"], pr["z"], pr["g0"]
    if len(K) < 2:
        continue
    n += 1
    cells, phi, lb, den = cells_of(cpdag, K, x, y, z, g0)
    phis.append(phi)
    rows = ball(pr)
    rvals.append(r_val(rows) if rows else None)
    for k, (cell, _) in cells.items():
        if cell != "data-refuted":
            continue
        drefuted_tot += 1
        alone = apply_orientations(cpdag, [(k[1], k[0])])   # the data's own verdict
        if alone is not None:
            laundered_tot += 1
            if not is_valid_adjustment_set_mpdag(alone, x, y, z):
                laundered_bad += 1
    if phi == 0:
        phi0 += 1
        if nu2_fatal(cpdag, K, x, y, z):
            phi0_and_pairfatal += 1

print(f"problems = {n}  (mean |K| implied)")
c = Counter(phis)
print("\nphi_1 distribution (all problems):")
for v in sorted(c):
    print(f"   phi_1 = {v:.4f}   {c[v]:3d}  {100*c[v]/n:5.1f}%")
print(f"   distinct values = {len(c)}   min={min(c):.3f} max={max(c):.3f}")
print(f"   mass at 0 = {100*c.get(0.0,0)/n:.1f}%   mass at 1 = {100*c.get(1.0,0)/n:.1f}%")

sel = [p for p, r in zip(phis, rvals) if r == 1]
cs = Counter(sel)
print(f"\nphi_1 INSIDE the r_val == 1 stratum (n={len(sel)}):")
for v in sorted(cs):
    print(f"   phi_1 = {v:.4f}   {cs[v]:3d}  {100*cs[v]/len(sel):5.1f}%")
print(f"   distinct = {len(cs)}   mass at 0 = {100*cs.get(0.0,0)/len(sel):.1f}%  "
      f"mass at 1 = {100*cs.get(1.0,0)/len(sel):.1f}%")

print(f"\nDATA-REFUTED LAUNDERING")
print(f"   claims labelled data-refuted           : {drefuted_tot}")
print(f"   of those, CPDAG ALONE permits reversal : {laundered_tot}"
      f"  ({100*laundered_tot/max(drefuted_tot,1):.1f}%)")
print(f"   ... and that reversal INVALIDATES Z    : {laundered_bad}"
      f"  ({100*laundered_bad/max(drefuted_tot,1):.1f}% of the cell)")

print(f"\nBLIND SPOT")
print(f"   phi_1 == 0 problems                    : {phi0}  ({100*phi0/n:.1f}%)")
print(f"   ... of which SOME admissible PAIR revision kills Z: {phi0_and_pairfatal}"
      f"  ({100*phi0_and_pairfatal/max(phi0,1):.1f}%)")
