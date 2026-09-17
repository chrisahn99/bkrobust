#!/usr/bin/env python3
"""FINAL DEFINITION, part 6 -- where the selective-widening bet stops paying.

The singleton tier (t <= 1) covers whenever at most one claim is false; its
coverage is therefore ~P(b_true <= 1) plus the cases where more errors happen not
to matter.  Sweep q and report where the tier falls below 0.95 and the calibrated
position jumps to the pairs tier, which costs the whole saving.
"""
import sys
import numpy as np
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
from final_rh_5_verify import cases, hull_upto
from final_rh_1_ladder import w

if __name__ == "__main__":
    print(" q     n   P(b>=2)  cov(t=0) cov(top1) cov(all1) cov(all2) | t*      width  blanket  saving")
    for q in (0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.45):
        C = cases(20260907, 400, q)
        bl = np.mean([w(c["pr"]["blanket"]) for c in C])
        cov, wid = {}, {}
        for name, fn in [("t=0", lambda c: c["base"]),
                         ("top-1", lambda c: hull_upto(c, 1, 1)),
                         ("all-1", lambda c: hull_upto(c, 1)),
                         ("all-2", lambda c: hull_upto(c, 2)),
                         ("all-3", lambda c: hull_upto(c, 3))]:
            ivs = [fn(c) for c in C]
            cov[name] = np.mean([iv[0] - 1e-8 <= c["t"] <= iv[1] + 1e-8 for iv, c in zip(ivs, C)])
            wid[name] = np.mean([w(iv) for iv in ivs])
        # INTEGER budget ladder only: b = 0, 1, 2, 3.  Keeps the exact guarantee.
        star = next((k for k in ("t=0", "all-1", "all-2", "all-3") if cov[k] >= 0.95), None)
        p2 = np.mean([c["b"] >= 2 for c in C])
        s = "  none" if star is None else f"{star:<7} {wid[star]:.4f} {bl:.4f}  {100*(1-wid[star]/bl):5.1f}%"
        print(f"{q:<5} {len(C):>4}  {p2:.3f}    {cov['t=0']:.3f}    {cov['top-1']:.3f}     "
              f"{cov['all-1']:.3f}     {cov['all-2']:.3f}   | {s}   "
              f"[w: b0={wid['t=0']:.4f} b1={wid['all-1']:.4f} b2={wid['all-2']:.4f}]")
