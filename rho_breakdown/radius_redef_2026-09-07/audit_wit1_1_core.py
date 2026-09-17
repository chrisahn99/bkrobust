#!/usr/bin/env python3
"""witness-1 audit 1: is L_flip non-degenerate, and is it r_val renamed?

Computes, per problem: R, L_flip, L_ret, r_val, and the hop-1 shell verdict.
K is elicited TRUTHFULLY here (the honest-expert case) so that the partition
is measured on its own terms before any lie is introduced.
"""
import sys
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import make_problem, ball, r_val
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.evaluate import is_valid_adjustment_set_mpdag


def partition(pr):
    """R, L_flip, L_ret for the problem's own K."""
    cpdag, K, x, y, z = pr["cpdag"], pr["K"], pr["x"], pr["y"], pr["z"]
    R, L_flip, L_ret = [], [], []
    for k in K:
        rest = [e for e in K if e != k]
        rev = rest + [(k[1], k[0])]
        g_rev = apply_orientations(cpdag, rev)
        if g_rev is None:
            R.append(k)
        elif not is_valid_adjustment_set_mpdag(g_rev, x, y, z):
            L_flip.append(k)
        g_ret = apply_orientations(cpdag, rest)
        if g_ret is None:
            L_ret.append(k)          # degenerate, shouldn't happen
        elif not is_valid_adjustment_set_mpdag(g_ret, x, y, z):
            L_ret.append(k)
    return R, L_flip, L_ret


if __name__ == "__main__":
    rng = np.random.default_rng(20260907)
    rows = []
    got, tried = 0, 0
    while got < 400 and tried < 40000:
        tried += 1
        pr = make_problem(rng)
        if pr is None:
            continue
        b = ball(pr)
        if not b:
            continue
        got += 1
        R, Lf, Lr = partition(pr)
        rows.append(dict(m=len(pr["K"]), R=len(R), Lf=len(Lf), Lr=len(Lr),
                         rval=r_val(b), Lf_set=set(Lf), Lr_set=set(Lr), R_set=set(R)))

    m = np.array([r["m"] for r in rows])
    nR = np.array([r["R"] for r in rows])
    nLf = np.array([r["Lf"] for r in rows])
    nLr = np.array([r["Lr"] for r in rows])
    rv = [r["rval"] for r in rows]

    print(f"problems: {len(rows)}")
    print(f"mean |K| = {m.mean():.2f}   mean |R| = {nR.mean():.2f}   "
          f"mean |L_flip| = {nLf.mean():.2f}   mean |L_ret| = {nLr.mean():.2f}")
    print()
    print("--- DEGENERACY RELOCATION: distribution of |L_flip| ---")
    c = Counter(nLf.tolist())
    for k in sorted(c):
        print(f"   |L_flip| = {k}: {c[k]:4d}  {100*c[k]/len(rows):5.1f}%")
    print(f"   => L_flip EMPTY on {100*(nLf==0).mean():.1f}% of problems")
    print(f"   => |L_flip| in {{0,1}} on {100*(nLf<=1).mean():.1f}% of problems")
    print()
    print("--- same question for |L_ret| (the screen L_flip claims to improve on) ---")
    c = Counter(nLr.tolist())
    for k in sorted(c):
        print(f"   |L_ret| = {k}: {c[k]:4d}  {100*c[k]/len(rows):5.1f}%")
    print(f"   => L_ret EMPTY on {100*(nLr==0).mean():.1f}%")
    print()
    print("--- LEMMA 2 CHECK: L_flip subset of L_ret ? ---")
    viol = [r for r in rows if not r["Lf_set"] <= r["Lr_set"]]
    print(f"   violations: {len(viol)} / {len(rows)}")
    print()
    print("--- IS IT r_val RENAMED? cross-tab (L_flip nonempty) x (r_val==1) ---")
    tab = Counter((r["Lf"] > 0, r["rval"] == 1) for r in rows)
    for (a, b_) in sorted(tab, key=lambda t: (t[0], t[1])):
        print(f"   L_flip!=0 {str(a):5s}  r_val==1 {str(b_):5s} : {tab[(a,b_)]:4d}")
    agree = sum(v for (a, b_), v in tab.items() if a == b_)
    print(f"   agreement: {100*agree/len(rows):.1f}%")
    print()
    print("--- and (L_ret nonempty) x (r_val==1) ---")
    tab = Counter((r["Lr"] > 0, r["rval"] == 1) for r in rows)
    for (a, b_) in sorted(tab, key=lambda t: (t[0], t[1])):
        print(f"   L_ret!=0 {str(a):5s}  r_val==1 {str(b_):5s} : {tab[(a,b_)]:4d}")
    agree = sum(v for (a, b_), v in tab.items() if a == b_)
    print(f"   agreement: {100*agree/len(rows):.1f}%")
    print()
    print(f"r_val distribution: {Counter('UNREACHED' if v is None else v for v in rv)}")
