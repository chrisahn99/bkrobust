#!/usr/bin/env python3
"""Audit 3: theta_D fixed (O* is invalid when cn(x,y) is empty), regime B only,
plus the delta sweep and the calibration-by-scaling test."""
import sys, time
import numpy as np
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag, optimal_adjustment_set_dag,
                                    is_valid_adjustment_set_dag, is_valid_adjustment_set_mpdag,
                                    adjusted_estimand, random_sem)
TOL = 1e-9

def random_dag(rng, n, p):
    order = list(rng.permutation(n)); e = set()
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < p: e.add((f"v{order[i]}", f"v{order[j]}"))
    return MPDAG(nodes=[f"v{i}" for i in range(n)], directed=e, undirected=[])

BAD_OSTAR = [0, 0]
def theta_D(sem, d, x, y):
    """The total effect DAG d implies, from the true observational law."""
    o = optimal_adjustment_set_dag(d, x, y)
    BAD_OSTAR[1] += 1
    if not is_valid_adjustment_set_dag(d, x, y, o):
        BAD_OSTAR[0] += 1
        o = d.parents(x)               # always a valid back-door set (no latents here)
        if not is_valid_adjustment_set_dag(d, x, y, o):
            return None
    return adjusted_estimand(sem, x, y, o)

def theta_set(sem, g, x, y):
    v = [theta_D(sem, d, x, y) for d in enumerate_dag_extensions(g)]
    return [t for t in v if t is not None]

def width(v): return (max(v)-min(v)) if v else 0.0

def fibre(cp, K, k):
    rest = [e for e in K if e != k]; out = []
    g = apply_orientations(cp, rest) if rest else cp
    if g is not None: out.append(g)
    g2 = apply_orientations(cp, rest + [(k[1], k[0])])
    if g2 is not None: out.append(g2)
    return out

def make(rng, n=9, p=0.30, k_claims=4, max_und=7):
    dag = random_dag(rng, n, p)
    if not dag.directed_edges: return None
    cp = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cp.undirected_edges)
    if not (2 <= len(und) <= max_und): return None
    sem = random_sem(dag, rng)
    nodes = sorted(dag.nodes)
    cand = [(x, y) for x in nodes for y in nodes if x != y and y in dag.descendants(x)]
    rng.shuffle(cand)
    for x, y in cand:
        thC = theta_set(sem, cp, x, y)
        if width(thC) <= TOL: continue          # BK irrelevant: CPDAG already identifies
        idx = rng.permutation(len(und))[:min(k_claims, len(und))]
        K = [((a, b) if (a, b) in dag.directed_edges else (b, a))
             for a, b in (und[t] for t in idx)]
        if not K: continue
        g0 = apply_orientations(cp, K)
        if g0 is None: continue
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None: continue
        return dict(dag=dag, cpdag=cp, x=x, y=y, K=[tuple(e) for e in K],
                    g0=g0, z=frozenset(z), sem=sem, thC=thC)
    return None

def hedge(sem, cp, K, x, y, delta):
    """Theta(G0) union the fibres of claims with lambda > delta. Returns (lo,hi,lams)."""
    g0 = apply_orientations(cp, K)
    th0 = theta_set(sem, g0, x, y); w0 = width(th0)
    lams, fib = [], {}
    for k in K:
        fv = []
        for g in fibre(cp, K, k): fv.extend(theta_set(sem, g, x, y))
        fib[k] = fv
        lams.append(width(th0 + fv) - w0)
    u = list(th0)
    for k, l in zip(K, lams):
        if l > delta: u.extend(fib[k])
    if not u: return None
    return min(u), max(u), lams, th0

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    rng = np.random.default_rng(20260907)
    P, tried = [], 0
    t0 = time.time()
    while len(P) < N and tried < 500000:
        tried += 1
        pr = make(rng)
        if pr is not None: P.append(pr)
    print(f"# {len(P)} problems, {tried} draws, {time.time()-t0:.1f}s")
    print(f"# theta_D fix fired: O*(D) was NOT a valid adjustment set on "
          f"{BAD_OSTAR[0]}/{BAD_OSTAR[1]} DAG evaluations ({100*BAD_OSTAR[0]/BAD_OSTAR[1]:.1f}%)"
          f"   <-- the definition's theta_D is wrong on these")

    # --- lambda profile on TRUE claim sets --------------------------------
    nz = z = 0; allzero = 0; lamvals = []
    for pr in P:
        h_ = hedge(pr["sem"], pr["cpdag"], pr["K"], pr["x"], pr["y"], 0.0)
        if h_ is None: continue
        lo, hi, lams, th0 = h_
        lamvals.extend(lams)
        for l in lams:
            if l > TOL: nz += 1
            else: z += 1
        if all(l <= TOL for l in lams): allzero += 1
    print(f"\n== T1' lambda sparsity (theta_D corrected, BK-relevant regime) ==")
    print(f"lambda==0 on {z}/{z+nz} claims ({100*z/(z+nz):.1f}%);  "
          f"{z/len(P):.2f} zeros per problem of {(z+nz)/len(P):.2f}")
    print(f"problems where EVERY claim has lambda==0 (ranked list fully tied, "
          f"hedge is a point): {allzero}/{len(P)} ({100*allzero/len(P):.1f}%)")

    # --- delta sweep, coverage + width under budget-b error ---------------
    print(f"\n== T3'/T6' delta sweep: coverage and width under b reversed claims ==")
    print(f"{'b':>2} {'delta':>7} {'cov':>6} {'width':>8} {'zero-w':>7} {'n':>5}")
    for b in (0, 1, 2):
        rng2 = np.random.default_rng(11)
        cases = []
        for pr in P:
            K = pr["K"]
            if len(K) < b: continue
            pick = set(rng2.permutation(len(K))[:b].tolist())
            Kst = [((e[1], e[0]) if i in pick else e) for i, e in enumerate(K)]
            g0 = apply_orientations(pr["cpdag"], Kst)
            if g0 is None: continue
            if optimal_adjustment_set_mpdag(g0, pr["x"], pr["y"]) is None: continue
            cases.append((pr, Kst))
        for delta in (0.0, 0.05, 0.2, 0.5):
            cov = zw = 0; ws = []; used = 0
            for pr, Kst in cases:
                h_ = hedge(pr["sem"], pr["cpdag"], Kst, pr["x"], pr["y"], delta)
                if h_ is None: continue
                lo, hi, lams, th0 = h_
                truth = pr["sem"].true_total_effect(pr["x"], pr["y"])
                if lo - 1e-8 <= truth <= hi + 1e-8: cov += 1
                ws.append(hi - lo); used += 1
                if hi - lo <= TOL: zw += 1
            print(f"{b:>2} {delta:>7.2f} {cov/max(used,1):>6.3f} {np.mean(ws):>8.4f} "
                  f"{zw/max(used,1):>7.3f} {used:>5}")
        # blanket baseline
        cov = 0; ws = []
        for pr, Kst in cases:
            truth = pr["sem"].true_total_effect(pr["x"], pr["y"])
            thC = pr["thC"]
            if min(thC) - 1e-8 <= truth <= max(thC) + 1e-8: cov += 1
            ws.append(width(thC))
        print(f"{b:>2} {'BLANKET':>7} {cov/max(len(cases),1):>6.3f} {np.mean(ws):>8.4f} "
              f"{0.0:>7.3f} {len(cases):>5}")

    # --- T4' calibration by scaling ---------------------------------------
    print(f"\n== T4' can the lambda hedge be calibrated to 95% by scaling? ==")
    for b in (1, 2):
        rng2 = np.random.default_rng(11)
        rows = []
        for pr in P:
            K = pr["K"]
            if len(K) < b: continue
            pick = set(rng2.permutation(len(K))[:b].tolist())
            Kst = [((e[1], e[0]) if i in pick else e) for i, e in enumerate(K)]
            g0 = apply_orientations(pr["cpdag"], Kst)
            if g0 is None: continue
            if optimal_adjustment_set_mpdag(g0, pr["x"], pr["y"]) is None: continue
            h_ = hedge(pr["sem"], pr["cpdag"], Kst, pr["x"], pr["y"], 0.0)
            if h_ is None: continue
            lo, hi, lams, th0 = h_
            c = 0.5*(lo+hi); h = 0.5*(hi-lo)
            rows.append((c, h, pr["sem"].true_total_effect(pr["x"], pr["y"])))
        # the scale s needed for 95% coverage of |truth-c| <= s*h
        need = []
        unreachable = 0
        for c, h, t in rows:
            if abs(t-c) <= 1e-9: need.append(0.0)
            elif h <= TOL: unreachable += 1     # zero half-width: NO scale ever covers
            else: need.append(abs(t-c)/h)
        print(f" b={b}: {len(rows)} cases; {unreachable} have zero width AND miss the truth "
              f"-> uncoverable by ANY scale ({100*unreachable/max(len(rows),1):.1f}%)")
        if unreachable/max(len(rows),1) > 0.05:
            print(f"        => 95% coverage by scaling is IMPOSSIBLE for this method at b={b}")
        else:
            s = np.quantile(need, 0.95) if need else float('nan')
            print(f"        scale for 95%: {s:.2f}; mean calibrated width "
                  f"{np.mean([2*s*h for c,h,t in rows]):.4f}")
