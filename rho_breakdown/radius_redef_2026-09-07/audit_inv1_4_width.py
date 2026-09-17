"""inverse-1 on the TARGET METRIC: width required for honest 95% coverage.
One elicited claim is REVERSED (false but data-consistent). Hedge over B_r(G0) with
r = r_alpha, versus the current definition (r_val - 1), versus the blanket hedge over
the whole space, versus no hedge. Calibrate every method by scaling its interval about
its midpoint until empirical coverage reaches 95%; report the mean calibrated width."""
import sys, time
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.space import (enumerate_space, represented_dags, covering_pairs,
                                 neighbour_graph, bfs_distances)
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, is_valid_adjustment_set_mpdag,
                                    optimal_adjustment_set_dag, is_valid_adjustment_set_dag,
                                    random_sem, adjusted_estimand)
from harness import random_dag

ALPHAS = [0.0, 0.05, 0.10, 0.25, 0.50]

def build(rng, n=6, p=0.35, kmax=3):
    dag = random_dag(rng, n, p)
    if not dag.directed_edges: return None
    cp = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cp.undirected_edges)
    if not (2 <= len(und) <= 5): return None
    nodes = sorted(dag.nodes); cand = [(x, y) for x in nodes for y in nodes if x != y]
    rng.shuffle(cand)
    for x, y in cand:
        zt = optimal_adjustment_set_dag(dag, x, y)
        if zt is None or not is_valid_adjustment_set_dag(dag, x, y, zt):
            continue          # the truth must itself admit adjustment, else nothing can cover
        km = min(kmax, len(und)); idx = rng.permutation(len(und))[:km]
        Kt = []
        for t in idx:
            a, b = und[t]; Kt.append((a, b) if (a, b) in dag.directed_edges else (b, a))
        jbad = int(rng.integers(len(Kt)))
        K = [(b, a) if i == jbad else (a, b) for i, (a, b) in enumerate(Kt)]
        g0 = apply_orientations(cp, K)
        if g0 is None: continue
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None: continue
        space = enumerate_space(cp)
        g0k = next((g for g in space if g == g0), None)
        if g0k is None: continue
        reps = represented_dags(space); covers = covering_pairs(space, reps)
        nbrs = neighbour_graph(space, covers); dist = bfs_distances(nbrs, g0k)
        if len(dist) < 4: continue
        sem = random_sem(dag, rng); tau = sem.true_total_effect(x, y)
        rows = []          # (distance, Z-valid-here, estimand of O*(g))
        for g in space:
            d = dist.get(g)
            if d is None: continue
            ok = bool(is_valid_adjustment_set_mpdag(g, x, y, z))
            zg = optimal_adjustment_set_mpdag(g, x, y)
            e = adjusted_estimand(sem, x, y, zg) if zg is not None else None
            rows.append((d, ok, e))
        if not any(e is not None for _, _, e in rows): continue
        return dict(rows=rows, tau=tau, est0=adjusted_estimand(sem, x, y, z),
                    rmax=max(d for d, _, _ in rows))
    return None

def sigma(rows, rmax):
    return {r: np.mean([ok for d, ok, _ in rows if d <= r]) for r in range(rmax + 1)}

def hedge(rows, r):
    v = [e for d, _, e in rows if d <= r and e is not None]
    return (min(v), max(v)) if v else None

def calibrated_width(iv, taus, target=0.95):
    """Smallest common scale k about the midpoint reaching >= target coverage; mean width."""
    lo = np.array([a for a, b in iv]); hi = np.array([b for a, b in iv]); t = np.array(taus)
    mid = (lo + hi) / 2; half = (hi - lo) / 2
    need = np.abs(t - mid)
    with np.errstate(divide='ignore', invalid='ignore'):
        k_need = np.where(half > 1e-12, need / np.maximum(half, 1e-300), np.inf)
    k_need = np.where((half <= 1e-12) & (need <= 1e-9), 0.0, k_need)
    reach = np.mean(np.isfinite(k_need))
    if reach < target - 1e-12:
        return None, reach, np.mean(2 * half)
    k = np.quantile(k_need[np.isfinite(k_need)], target / reach) if reach > 0 else np.inf
    k = max(k, 1.0)
    return float(np.mean(2 * half * k)), reach, float(np.mean(2 * half))

rng = np.random.default_rng(4242); P = []; t0 = time.time()
while len(P) < 400 and time.time() - t0 < 500:
    pr = build(rng)
    if pr: P.append(pr)
print(f"problems = {len(P)}   ({time.time()-t0:.0f}s)")

methods = {}
for p in P:
    p["sig"] = sigma(p["rows"], p["rmax"])
    bad = [d for d, ok, _ in p["rows"] if not ok and d > 0]
    p["rval"] = min(bad) if bad else None

def radius_of(p, kind, a=None):
    if kind == "blanket": return p["rmax"]
    if kind == "none": return 0
    if kind == "rval":   # the current definition's certified-safe ball
        return (p["rval"] - 1) if p["rval"] is not None else p["rmax"]
    S = [r for r in range(p["rmax"] + 1) if p["sig"][r] >= 1 - a - 1e-12]
    return max(S) if S else -1

specs = [("no hedge (point estimate)", "none", None),
         ("current defn: hedge over B_{r_val-1}", "rval", None)]
specs += [(f"r_alpha, alpha={a}", "alpha", a) for a in ALPHAS]
specs += [("blanket hedge over the whole space", "blanket", None)]

print(f"\n{'method':44s} {'raw cov':>8} {'raw wid':>9} {'max cov':>8} {'CAL WIDTH @95%':>15} {'mean r':>7} {'r=-1':>6}")
for label, kind, a in specs:
    iv = []; taus = []; rr = []; unr = 0
    for p in P:
        r = radius_of(p, kind, a); rr.append(r)
        if r < 0:
            unr += 1
            h = hedge(p["rows"], p["rmax"])   # UNREACHED: only honest fallback is the blanket
        else:
            h = hedge(p["rows"], r)
        if h is None: h = (p["est0"], p["est0"])
        iv.append(h); taus.append(p["tau"])
    cov = np.mean([lo - 1e-9 <= t <= hi + 1e-9 for (lo, hi), t in zip(iv, taus)])
    cw, reach, raw = calibrated_width(iv, taus)
    cws = "UNCALIBRATABLE" if cw is None else f"{cw:.4f}"
    print(f"{label:44s} {cov:8.3f} {raw:9.4f} {reach:8.3f} {cws:>15} {np.mean(rr):7.2f} {unr:6d}")
