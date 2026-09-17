#!/usr/bin/env python3
"""inverse-2 audit, part 4: (a) exact cross-tab phi_1 vs r_val; (b) padding attack."""
import sys
from collections import Counter, defaultdict
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import make_problem, ball, r_val
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.evaluate import is_valid_adjustment_set_mpdag


def phi1(cpdag, K, x, y, z, g0):
    lb = den = 0
    flagged = []
    for k in K:
        rest = [e for e in K if e != k]
        g_ret = apply_orientations(cpdag, rest)
        g_rev = apply_orientations(cpdag, rest + [(k[1], k[0])])
        adm = [g for g in (g_ret, g_rev) if g is not None]
        if adm:
            den += 1
        if g_rev is None:
            continue                                    # data-refuted
        if g_ret is not None and g_ret == g0:
            continue                                    # inert
        if any(not is_valid_adjustment_set_mpdag(g, x, y, z) for g in adm):
            lb += 1; flagged.append(k)
    return (lb / den if den else None), lb, den, flagged


rng = np.random.default_rng(20260907)
tab = Counter(); n = 0; tried = 0
pad_rows = []
while n < 250 and tried < 40000:
    tried += 1
    pr = make_problem(rng, k_claims=4)
    if pr is None: continue
    cpdag, K, x, y, z, g0 = pr["cpdag"], pr["K"], pr["x"], pr["y"], pr["z"], pr["g0"]
    if len(K) < 2: continue
    rows = ball(pr)
    if not rows: continue
    n += 1
    phi, lb, den, flagged = phi1(cpdag, K, x, y, z, g0)
    r = r_val(rows)
    tab[(phi > 0, r)] += 1

    # --- padding attack: append claims already entailed by G0 (compelled edges of Chat,
    #     plus edges Meek derives from K).  G0, Z, the estimate, the truth: all unchanged.
    entailed = [e for e in sorted(g0.directed_edges) if tuple(e) not in set(K)]
    if entailed:
        Kpad = list(K) + [tuple(e) for e in entailed]
        g0p = apply_orientations(cpdag, Kpad)
        assert g0p == g0, "padding must not move G0"
        phip, lbp, denp, _ = phi1(cpdag, Kpad, x, y, z, g0)
        pad_rows.append((phi, phip, len(K), len(Kpad), lb, lbp))

print("CROSS-TAB  (phi_1 > 0)  x  r_val")
by_r = defaultdict(Counter)
for (pos, r), c in tab.items():
    by_r[r][pos] += c
for r in sorted(by_r, key=lambda v: (v is None, v)):
    cc = by_r[r]
    print(f"   r_val = {str(r):<10}  phi_1>0: {cc[True]:3d}   phi_1==0: {cc[False]:3d}")
eq = all((r == 1) == pos for (pos, r), c in tab.items() for _ in range(1))
viol = sum(c for (pos, r), c in tab.items() if (r == 1) != pos)
print(f"   violations of 'phi_1 > 0  <=>  r_val == 1' : {viol} / {n}"
      f"   ({100*viol/n:.1f}%)")

print("\nPADDING ATTACK  (append edges G0 already entails; G0/Z/estimate unchanged)")
if pad_rows:
    a = np.array([(p, q, k1, k2) for p, q, k1, k2, _, _ in pad_rows], float)
    moved = a[a[:,0] != a[:,1]]
    print(f"   problems with at least one paddable edge : {len(a)}")
    print(f"   phi_1 changed by padding                : {len(moved)}  "
          f"({100*len(moved)/len(a):.1f}%)")
    print(f"   mean phi_1 before = {a[:,0].mean():.4f}   after = {a[:,1].mean():.4f}")
    print(f"   mean |K| before = {a[:,2].mean():.1f}   after = {a[:,3].mean():.1f}")
    nz = a[a[:,0] > 0]
    if len(nz):
        print(f"   among phi_1 > 0 : mean {nz[:,0].mean():.4f} -> {nz[:,1].mean():.4f}"
              f"   (mean shrink factor {np.mean(nz[:,1]/nz[:,0]):.3f})")
        print(f"   worst single case: {nz[np.argmin(nz[:,1]/nz[:,0])][0]:.4f} -> "
              f"{nz[np.argmin(nz[:,1]/nz[:,0])][1]:.4f}")
    lbs = np.array([(l1, l2) for _, _, _, _, l1, l2 in pad_rows], float)
    print(f"   numerator (# load-bearing) changed?     : "
          f"{int((lbs[:,0] != lbs[:,1]).sum())} / {len(lbs)}")
