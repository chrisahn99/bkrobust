"""inverse-1: does the advertised identity r_0 = r_val - 1 hold on random problems,
and what does the UNREACHED convention do when sigma(0) < 1?"""
import sys, time
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import make_problem, ball, r_val

rng = np.random.default_rng(20260907)
P = []; t0 = time.time()
while len(P) < 500 and time.time() - t0 < 300:
    pr = make_problem(rng)
    if pr is None: continue
    rows = ball(pr)
    if not rows: continue
    rmax = max(d for d, _ in rows)
    if rmax < 1: continue
    P.append((rows, rmax))

n = len(P); s0bad = 0; agree = 0; disagree = []
for rows, rmax in P:
    sig = {r: np.mean([ok for d, ok in rows if d <= r]) for r in range(rmax + 1)}
    if sig[0] < 1 - 1e-12: s0bad += 1
    S = [r for r in range(rmax + 1) if sig[r] >= 1 - 1e-12]
    r0 = max(S) if S else -1
    rv = r_val(rows)
    tgt = (rv - 1) if rv is not None else rmax     # UNREACHED -> the whole ball is safe
    if r0 == tgt: agree += 1
    else: disagree.append((r0, rv, tgt, sig[0]))
print(f"problems = {n}")
print(f"Z = O*(G0) fails is_valid_adjustment_set_mpdag AT G0 itself in {100*s0bad/n:.1f}% of problems"
      f"  -> sigma(0) = 0 -> r_alpha = UNREACHED(-1) for every alpha < 1")
print(f"identity r_0 == r_val - 1 holds in {100*agree/n:.1f}% of problems")
c = Counter((a, b if b is not None else 'UNR') for a, b, _, _ in disagree)
for k, v in c.most_common(8):
    print(f"   counterexample r_0={k[0]}, r_val={k[1]}   x{v}")
