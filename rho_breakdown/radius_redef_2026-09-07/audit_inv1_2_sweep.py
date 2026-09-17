"""inverse-1: does the alpha knob buy a usable middle range, or does it swap a floor
collapse (r_val=1) for a ceiling collapse (r_alpha = eccentricity)?"""
import sys, time
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import make_problem, ball, r_val

ALPHAS = [0.0, 0.02, 0.05, 0.10, 0.15, 0.25, 0.40, 0.50]

def sigma(rows, rmax):
    return {r: np.mean([ok for d, ok in rows if d <= r]) for r in range(rmax + 1)}

def ralpha(sig, rmax, a):
    S = [r for r in range(rmax + 1) if sig[r] >= 1 - a - 1e-12]
    return max(S) if S else -1

rng = np.random.default_rng(20260907)
P = []
t0 = time.time()
while len(P) < 500 and time.time() - t0 < 420:
    pr = make_problem(rng)
    if pr is None: continue
    rows = ball(pr)
    if not rows: continue
    rmax = max(d for d, _ in rows)
    if rmax < 1: continue
    sig = sigma(rows, rmax)
    P.append(dict(rows=rows, rmax=rmax, sig=sig, rval=r_val(rows),
                  hop1_broken=any((not ok) for d, ok in rows if d == 1),
                  n=len(rows), sig_glob=sig[rmax]))
print(f"problems = {len(P)}   ({time.time()-t0:.0f}s)   mean space size {np.mean([p['n'] for p in P]):.1f}")
rv = Counter("UNREACHED" if p["rval"] is None else p["rval"] for p in P)
print("r_val: " + "  ".join(f"{k}:{100*rv[k]/len(P):.1f}%" for k in sorted(rv, key=lambda v: (v=='UNREACHED', v))))
print(f"mean eccentricity r_max = {np.mean([p['rmax'] for p in P]):.2f}")
print(f"global safe fraction sigma(r_max): mean {np.mean([p['sig_glob'] for p in P]):.3f}  "
      f"median {np.median([p['sig_glob'] for p in P]):.3f}")
sat_thr = [1 - p["sig_glob"] for p in P]
print(f"SATURATION THRESHOLD 1-sigma_global (alpha above it => r_alpha = eccentricity, no matter "
      f"what the local structure is): mean {np.mean(sat_thr):.3f}  median {np.median(sat_thr):.3f}  "
      f"p10 {np.percentile(sat_thr,10):.3f}  p90 {np.percentile(sat_thr,90):.3f}")

print(f"\n{'alpha':>6} {'r_a=0':>7} {'r_a=-1':>7} {'r_a=ecc':>8} {'mean r_a':>9} {'distinct':>9} "
      f"{'holes':>7} {'ra>=2 & hop1 broken':>21}")
for a in ALPHAS:
    ras = [ralpha(p["sig"], p["rmax"], a) for p in P]
    holes = 0; game = 0
    for p, ra in zip(P, ras):
        if ra >= 0 and any(p["sig"][r] < 1 - a - 1e-12 for r in range(ra + 1)):
            holes += 1
        if ra >= 2 and p["hop1_broken"]:
            game += 1
    print(f"{a:6.2f} {100*np.mean([r==0 for r in ras]):6.1f}% {100*np.mean([r==-1 for r in ras]):6.1f}% "
          f"{100*np.mean([r==p['rmax'] for r,p in zip(ras,P)]):7.1f}% {np.mean(ras):9.2f} "
          f"{len(set(ras)):9d} {100*holes/len(P):6.1f}% {100*game/len(P):20.1f}%")

# how much of the ball at r_alpha is actually the whole space?
print("\nfraction of the space swept into the hedge, |B_{r_alpha}|/|space|:")
for a in ALPHAS:
    fr = []
    for p in P:
        ra = ralpha(p["sig"], p["rmax"], a)
        fr.append(0.0 if ra < 0 else np.mean([d <= ra for d, _ in p["rows"]]))
    print(f"  alpha={a:4.2f}: mean {np.mean(fr):.3f}   median {np.median(fr):.3f}   "
          f"=1.0 in {100*np.mean([f>=0.999 for f in fr]):.1f}% of problems")

# per-problem: the alpha at which r_alpha first exceeds 0, and the alpha at which it saturates
first = []; sat = []
grid = np.linspace(0, 0.9, 91)
for p in P:
    f = next((a for a in grid if ralpha(p["sig"], p["rmax"], a) >= 1), None)
    s = next((a for a in grid if ralpha(p["sig"], p["rmax"], a) >= p["rmax"]), None)
    if f is not None: first.append(f)
    if s is not None: sat.append(s)
print(f"\nalpha needed for ANY widening (r_alpha>=1): median {np.median(first):.3f}  "
      f"p90 {np.percentile(first,90):.3f}   (undefined for {100*(1-len(first)/len(P)):.1f}% of problems)")
print(f"alpha at which r_alpha SATURATES to the eccentricity: median {np.median(sat):.3f}  "
      f"p10 {np.percentile(sat,10):.3f}")
gap = [s - f for p, f, s in zip(P, first, sat)] if len(first) == len(sat) == len(P) else None
if gap is not None:
    print(f"width of the usable alpha window (saturate - first widening): median {np.median(gap):.3f}, "
          f"zero or negative in {100*np.mean([g<=0 for g in gap]):.1f}% of problems")
