"""SEL_A audit 2 -- with FALSE claims: does the floor ever bind, does the ranking
retrieve the false claim better than a free baseline, and what does the honest-95%
width table look like against RHIG at the same screening order.
"""
import sys, itertools, time
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_dag, optimal_adjustment_set_mpdag,
                                    is_valid_adjustment_set_dag, random_sem,
                                    adjusted_estimand, asymptotic_variance)
from harness import random_dag

Z196 = 1.959963985


def tau_hull(g, sem, x, y, cache):
    key = g.edge_string()
    if key in cache:
        return cache[key]
    vals = []
    for D in enumerate_dag_extensions(g):
        try:
            vals.append(adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(D, x, y)))
        except Exception:
            pass
    out = (min(vals), max(vals)) if vals else None
    cache[key] = out
    return out


def hull_of(cp, claims, sem, x, y, cache):
    g = apply_orientations(cp, list(claims))
    return None if g is None else tau_hull(g, sem, x, y, cache)


def build(rng, n_false, n_nodes=6, p=0.35, kmax=3):
    """A problem whose elicited K contains exactly n_false REVERSED claims."""
    dag = random_dag(rng, n_nodes, p)
    if not dag.directed_edges:
        return None
    cp = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cp.undirected_edges)
    if not (2 <= len(und) <= 6):
        return None
    nodes = sorted(dag.nodes)
    cand = [(x, y) for x in nodes for y in nodes if x != y]
    rng.shuffle(cand)
    sem = random_sem(dag, rng)
    for x, y in cand:
        km = min(kmax, len(und))
        if km <= n_false:
            continue
        idx = rng.permutation(len(und))[:km]
        Ktrue = []
        for t in idx:
            a, b = und[t]
            Ktrue.append((a, b) if (a, b) in dag.directed_edges else (b, a))
        flip = list(rng.permutation(km)[:n_false])
        K = [(e[1], e[0]) if i in flip else e for i, e in enumerate(Ktrue)]
        g0 = apply_orientations(cp, K)
        if g0 is None:
            continue
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None:
            continue
        try:
            th0 = adjusted_estimand(sem, x, y, z)
            av = asymptotic_variance(sem, x, y, z)
        except Exception:
            continue
        return dict(dag=dag, cp=cp, x=x, y=y, K=K, false_idx=set(flip), g0=g0, z=set(z),
                    sem=sem, theta0=th0, avar=av,
                    tau_true=sem.true_total_effect(x, y))
    return None


def methods(pr, nobs, cache):
    cp, x, y, K, sem, th0 = pr["cp"], pr["x"], pr["y"], pr["K"], pr["sem"], pr["theta0"]
    m = len(K)
    se = np.sqrt(pr["avar"] / nobs)
    delta = Z196 * se
    D1 = []
    for i in range(m):
        h = hull_of(cp, [K[j] for j in range(m) if j != i], sem, x, y, cache)
        D1.append(0.0 if h is None else max(abs(h[0]-th0), abs(h[1]-th0)))
    D1 = np.array(D1)
    A1 = [i for i in range(m) if D1[i] > delta]
    D2 = {}
    A2 = set(A1)
    for i, j in itertools.combinations(range(m), 2):
        h = hull_of(cp, [K[t] for t in range(m) if t not in (i, j)], sem, x, y, cache)
        d = 0.0 if h is None else max(abs(h[0]-th0), abs(h[1]-th0))
        D2[(i, j)] = d
        if d > delta:
            A2 |= {i, j}
    A2 = sorted(A2)

    def hull_excl(S):
        h = hull_of(cp, [K[i] for i in range(m) if i not in S], sem, x, y, cache)
        return (th0, th0) if h is None else h

    def rhig(j):
        lo, hi = np.inf, -np.inf
        for S in itertools.combinations(range(m), m-j):
            h = hull_of(cp, [K[i] for i in S], sem, x, y, cache)
            if h is None:
                continue
            lo, hi = min(lo, h[0]), max(hi, h[1])
        return (th0, th0) if lo == np.inf else (lo, hi)

    out = {
        "point": (th0, th0),
        "SEL_A1": hull_excl(set(A1)),
        "SEL_A2": hull_excl(set(A2)),
        "RHIG_1": rhig(1),
        "RHIG_2": rhig(2),
        "blanket": tau_hull(cp, sem, x, y, cache) or (th0, th0),
    }
    return out, se, delta, D1, D2, A1, A2


def calibrated_width(intervals, ses, truths, th0s, target=0.95):
    """Scale ignorance half-widths and the noise term by lambda until coverage >= target."""
    lo = np.array([iv[0] for iv in intervals]); hi = np.array([iv[1] for iv in intervals])
    se = np.array(ses); tau = np.array(truths); th = np.array(th0s)
    best = None
    for lam in np.concatenate([np.arange(1.0, 20.01, 0.05), np.arange(21, 501, 1.0)]):
        L = th - lam*((th - lo) + Z196*se)
        H = th + lam*((hi - th) + Z196*se)
        cov = np.mean((tau >= L - 1e-12) & (tau <= H + 1e-12))
        if cov >= target:
            best = (lam, cov, np.mean(H - L))
            break
    if best is None:
        L = th - 500*((th - lo) + Z196*se); H = th + 500*((hi - th) + Z196*se)
        best = (np.inf, np.mean((tau >= L) & (tau <= H)), np.mean(H - L))
    return best


def run(n_false, nobs, nprob, seed):
    rng = np.random.default_rng(seed)
    rows, tried, t0 = [], 0, time.time()
    while len(rows) < nprob and tried < 100000 and time.time()-t0 < 300:
        tried += 1
        pr = build(rng, n_false)
        if pr is None:
            continue
        cache = {}
        try:
            iv, se, delta, D1, D2, A1, A2 = methods(pr, nobs, cache)
        except Exception:
            continue
        rows.append(dict(pr=pr, iv=iv, se=se, delta=delta, D1=D1, D2=D2, A1=A1, A2=A2))
    return rows


if __name__ == "__main__":
    NOBS = 10000
    for n_false in (1, 2):
        rows = run(n_false, NOBS, 300, 4242 + n_false)
        print(f"\n{'='*78}\n{n_false} FALSE claim(s) in K   n_obs={NOBS}   problems={len(rows)}")
        dam = np.array([abs(r["pr"]["tau_true"] - r["pr"]["theta0"]) > 1e-9 for r in rows])
        print(f"  the false claim(s) actually move the estimand in {100*dam.mean():.1f}% "
              f"(silent elsewhere)")

        # --- does the floor delta ever bind? ---
        allD = np.concatenate([r["D1"] for r in rows])
        alld = np.concatenate([np.full(len(r["D1"]), r["delta"]) for r in rows])
        inband = (allD > 1e-12) & (allD <= alld)
        print(f"  [floor] Delta_i values: exactly 0 in {100*np.mean(allD<=1e-12):.1f}%, "
              f"in (0, delta] in {100*inband.mean():.2f}%, > delta in "
              f"{100*np.mean(allD>alld):.1f}%   (n={len(allD)} claims)")
        pos = allD[allD > 1e-12]
        if len(pos):
            print(f"  [floor] smallest strictly positive Delta_i = {pos.min():.5f}; "
                  f"mean delta = {alld.mean():.5f}; "
                  f"ratio min(Delta+)/delta = {pos.min()/alld.mean():.2f}")

        # --- coverage and honest width ---
        print(f"  {'method':10s} {'raw cov':>8s} {'raw width':>10s} {'lambda':>8s} "
              f"{'cal cov':>8s} {'HONEST WIDTH':>13s}")
        for name in ("point", "SEL_A1", "SEL_A2", "RHIG_1", "RHIG_2", "blanket"):
            ivs = [r["iv"][name] for r in rows]
            ses = [r["se"] for r in rows]
            tau = [r["pr"]["tau_true"] for r in rows]
            th0 = [r["pr"]["theta0"] for r in rows]
            L = np.array([t - (t-iv[0]) - Z196*s for iv, s, t in zip(ivs, ses, th0)])
            H = np.array([t + (iv[1]-t) + Z196*s for iv, s, t in zip(ivs, ses, th0)])
            raw_cov = np.mean((np.array(tau) >= L) & (np.array(tau) <= H))
            raw_w = np.mean(H - L)
            lam, cov, wid = calibrated_width(ivs, ses, tau, th0)
            print(f"  {name:10s} {raw_cov:8.3f} {raw_w:10.4f} {lam:8.2f} {cov:8.3f} {wid:13.4f}")

        # --- SEL_A1 vs RHIG_1 width, head to head ---
        s1 = np.array([r["iv"]["SEL_A1"][1]-r["iv"]["SEL_A1"][0] for r in rows])
        r1 = np.array([r["iv"]["RHIG_1"][1]-r["iv"]["RHIG_1"][0] for r in rows])
        s2 = np.array([r["iv"]["SEL_A2"][1]-r["iv"]["SEL_A2"][0] for r in rows])
        r2 = np.array([r["iv"]["RHIG_2"][1]-r["iv"]["RHIG_2"][0] for r in rows])
        print(f"  [head-to-head] SEL_A1 narrower than RHIG_1: "
              f"{100*np.mean(s1 < r1-1e-9):.1f}%   wider: {100*np.mean(s1 > r1+1e-9):.1f}%   "
              f"identical: {100*np.mean(np.abs(s1-r1)<1e-9):.1f}%")
        print(f"                 SEL_A2 vs RHIG_2: narrower {100*np.mean(s2 < r2-1e-9):.1f}%  "
              f"wider {100*np.mean(s2 > r2+1e-9):.1f}%  identical {100*np.mean(np.abs(s2-r2)<1e-9):.1f}%")

        # --- retrieval: does the Delta ranking find the false claim? ---
        # scored only where the false claim actually damages
        hits_d = hits_b = hits_r = n_sc = 0
        for r in rows:
            pr = r["pr"]
            if abs(pr["tau_true"] - pr["theta0"]) <= 1e-9:
                continue
            n_sc += 1
            m = len(pr["K"])
            D = r["D1"]
            top = int(np.argmax(D))
            hits_d += top in pr["false_idx"]
            # free baseline: rank by whether the claim touches x or y
            touch = np.array([1.0 if (pr["x"] in e or pr["y"] in e) else 0.0 for e in pr["K"]])
            hits_b += int(np.argmax(touch)) in pr["false_idx"] if touch.max() > 0 else 0
            hits_r += len(pr["false_idx"]) / m          # random guess expectation
        if n_sc:
            print(f"  [retrieval, on the {n_sc} damaging problems] precision@1: "
                  f"Delta-ranking {hits_d/n_sc:.3f} | 'claim touches X or Y' "
                  f"{hits_b/n_sc:.3f} | random {hits_r/n_sc:.3f}")

    # --- floor vs sample size ---
    print(f"\n{'='*78}\n[floor vs n] does delta ever separate anything?")
    rows = run(1, 10000, 200, 99)
    allD = np.concatenate([r["D1"] for r in rows])
    base = np.concatenate([np.full(len(r["D1"]), r["delta"]*np.sqrt(10000)) for r in rows])
    for nobs in (100, 1000, 10**4, 10**5, 10**6, 10**8):
        d = base/np.sqrt(nobs)
        A_empty = []
        k = 0
        for r in rows:
            m = len(r["D1"])
            dd = r["delta"]*np.sqrt(10000/nobs)
            A_empty.append(not np.any(r["D1"] > dd))
            k += m
        inband = np.mean((allD > 1e-12) & (allD <= d))
        print(f"   n={nobs:>9}  A1 empty in {100*np.mean(A_empty):5.1f}%   "
              f"Delta_i in (0,delta] for {100*inband:5.2f}% of claims")
