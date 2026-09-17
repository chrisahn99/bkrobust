"""Criterion (c): what interval does the r^alpha hedge produce, and does it cover?

Problems are generated with ONE elicited claim deliberately REVERSED (the paper's
own failure mode: K consistent with the data but false). The true DAG is then
outside [G0]. We build the hedge exactly as the definition's
`how_it_yields_an_interval_width` states:
    S(alpha) = union over {G : 0 < dist(G0,G) < r^alpha} of {d in [G] : Z invalid in d}
    interval = hull of { point estimate } U { adjusted_estimand at O*(d) : d in S }
and compare against the blanket hedge over the whole class [Chat].
"""
import sys
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import random_dag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.space import (enumerate_space, represented_dags, covering_pairs,
                                 neighbour_graph, bfs_distances)
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, optimal_adjustment_set_dag,
                                    is_valid_adjustment_set_dag, random_sem, adjusted_estimand)

ALPHAS = [0.0, 0.25, 0.5]

def one(rng, n=6, p=0.35, k=3):
    dag = random_dag(rng, n, p)
    if not dag.directed_edges: return None
    cpdag = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cpdag.undirected_edges)
    if not (2 <= len(und) <= 6): return None
    nodes = sorted(dag.nodes)
    cand = [(x, y) for x in nodes for y in nodes if x != y]
    rng.shuffle(cand)
    for x, y in cand:
        idx = rng.permutation(len(und))[:min(k, len(und))]
        K = []
        for j, t in enumerate(idx):
            a, b = und[t]
            true_or = (a, b) if (a, b) in dag.directed_edges else (b, a)
            K.append((true_or[1], true_or[0]) if j == 0 else true_or)   # first claim REVERSED
        g0 = apply_orientations(cpdag, K)
        if g0 is None: continue
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None: continue
        z = frozenset(z)
        space = enumerate_space(cpdag)
        if not (2 <= len(space) <= 300): continue
        reps = represented_dags(space)
        g0s = next((g for g in space if g == g0), None)
        if g0s is None: continue
        if any(not is_valid_adjustment_set_dag(d, x, y, z) for d in reps[g0s]):
            continue                                   # ill-posed at its own origin
        covers = covering_pairs(space, reps)
        nbrs = neighbour_graph(space, covers)
        dist = bfs_distances(nbrs, g0s)
        beta = {}
        for g, d in dist.items():
            ext = reps[g]
            beta[g] = (d, [dd for dd in ext if not is_valid_adjustment_set_dag(dd, x, y, z)], len(ext))
        sem = random_sem(dag, rng)
        truth = adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(dag, x, y))
        point = adjusted_estimand(sem, x, y, z)
        damaged = not is_valid_adjustment_set_dag(dag, x, y, z)
        out = {}
        for a in ALPHAS:
            cands = [dd for g, (dd, bad, ntot) in beta.items() if dd > 0 and len(bad)/ntot > a]
            r = min(cands) if cands else None
            vals = [point]
            for g, (dd, bad, ntot) in beta.items():
                if dd > 0 and (r is None or dd < r):
                    for b in bad:
                        vals.append(adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(b, x, y)))
            out[a] = (min(vals), max(vals), r)
        allext = enumerate_dag_extensions(cpdag)
        bl = [adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(d, x, y)) for d in allext]
        out["blanket"] = (min(bl + [point]), max(bl + [point]), None)
        return dict(truth=truth, point=point, damaged=damaged, iv=out)
    return None

def calib_width(rows, key):
    """Smallest lambda scaling each interval about its centre that gives >=95% coverage."""
    lo = np.array([r["iv"][key][0] for r in rows]); hi = np.array([r["iv"][key][1] for r in rows])
    t = np.array([r["truth"] for r in rows])
    c = (lo + hi) / 2; h = (hi - lo) / 2
    need = np.abs(t - c)
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = np.where(h > 0, need / np.where(h > 0, h, 1), np.where(need <= 1e-12, 0.0, np.inf))
    lam = np.sort(lam)
    k = int(np.ceil(0.95 * len(lam))) - 1
    lam95 = lam[k]
    if not np.isfinite(lam95):
        return np.inf, np.inf, np.mean(2 * h)
    return lam95, np.mean(2 * h * lam95), np.mean(2 * h)

rng = np.random.default_rng(11)
rows = []
tried = 0
while len(rows) < 150 and tried < 20000:
    tried += 1
    r = one(rng)
    if r: rows.append(r)
dmg = [r for r in rows if r["damaged"]]
print(f"problems (one claim reversed): {len(rows)}   truth-damaged (Z invalid in the TRUE dag): {len(dmg)}")
print(f"{'method':>12} {'cov@raw':>8} {'mean width':>11} {'lambda for 95%':>15} {'calibrated width':>17}")
for key in ALPHAS + ["blanket"]:
    cov = np.mean([r["iv"][key][0] - 1e-9 <= r["truth"] <= r["iv"][key][1] + 1e-9 for r in rows])
    lam, cw, w = calib_width(rows, key)
    nm = f"r^{key}" if key != "blanket" else "blanket"
    print(f"{nm:>12} {cov:>8.3f} {w:>11.4f} {lam:>15.3f} {cw:>17.4f}")
print()
for key in ALPHAS + ["blanket"]:
    cov = np.mean([r["iv"][key][0] - 1e-9 <= r["truth"] <= r["iv"][key][1] + 1e-9 for r in dmg])
    w = np.mean([r["iv"][key][1] - r["iv"][key][0] for r in dmg])
    nm = f"r^{key}" if key != "blanket" else "blanket"
    print(f"  on the {len(dmg)} damaged problems: {nm:>10} coverage {cov:.3f}  mean width {w:.4f}")
rs = [r["iv"][0.0][2] for r in rows]
from collections import Counter
print("\nr^0 (= r_val) distribution:", Counter("inf" if v is None else v for v in rs))
print("r^0.5 distribution:", Counter("inf" if r["iv"][0.5][2] is None else r["iv"][0.5][2] for r in rows))
