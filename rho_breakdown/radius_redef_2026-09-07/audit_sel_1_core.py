"""SEL_A audit 1 -- structure of the screen, the 2-delta sandwich, bimodality,
one-sidedness, and the redundancy game.  Population functionals only.
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
    """(lo, hi, n_ext, n_nonamenable) over DAG extensions of g, tau via O*(D)."""
    key = g.edge_string()
    if key in cache:
        return cache[key]
    vals, bad = [], 0
    for D in enumerate_dag_extensions(g):
        try:
            z = optimal_adjustment_set_dag(D, x, y)
            if not is_valid_adjustment_set_dag(D, x, y, z):
                bad += 1
            vals.append(adjusted_estimand(sem, x, y, z))
        except Exception:
            pass
    out = (min(vals), max(vals), len(vals), bad) if vals else None
    cache[key] = out
    return out


def hull_of(cpdag, claims, sem, x, y, cache):
    g = apply_orientations(cpdag, list(claims))
    if g is None:
        return None
    return tau_hull(g, sem, x, y, cache)


def build(rng, n=6, p=0.35, kmax=3):
    dag = random_dag(rng, n, p)
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
        idx = rng.permutation(len(und))[:km]
        K = []
        for t in idx:
            a, b = und[t]
            K.append((a, b) if (a, b) in dag.directed_edges else (b, a))
        g0 = apply_orientations(cp, K)
        if g0 is None:
            continue
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None:
            continue
        try:
            theta0 = adjusted_estimand(sem, x, y, z)
            av = asymptotic_variance(sem, x, y, z)
        except Exception:
            continue
        return dict(dag=dag, cp=cp, x=x, y=y, K=K, g0=g0, z=set(z), sem=sem,
                    theta0=theta0, avar=av, und=und,
                    tau_true=sem.true_total_effect(x, y))
    return None


def analyse(pr, nobs, cache):
    cp, x, y, K, sem, th0 = pr["cp"], pr["x"], pr["y"], pr["K"], pr["sem"], pr["theta0"]
    m = len(K)
    delta = Z196 * np.sqrt(pr["avar"] / nobs)          # delta_noise
    # first-order screen
    D1, H1 = [], []
    none_closure = 0
    for i in range(m):
        h = hull_of(cp, [K[j] for j in range(m) if j != i], sem, x, y, cache)
        if h is None:
            none_closure += 1
            D1.append(0.0); H1.append((th0, th0, 0, 0)); continue
        D1.append(max(abs(h[0] - th0), abs(h[1] - th0)))
        H1.append(h)
    D1 = np.array(D1)
    A1 = [i for i in range(m) if D1[i] > delta]
    # second-order
    D2 = {}
    for i, j in itertools.combinations(range(m), 2):
        h = hull_of(cp, [K[t] for t in range(m) if t not in (i, j)], sem, x, y, cache)
        D2[(i, j)] = 0.0 if h is None else max(abs(h[0] - th0), abs(h[1] - th0))
    A2 = set(A1)
    for (i, j), d in D2.items():
        if d > delta:
            A2 |= {i, j}
    A2 = sorted(A2)
    sel1 = hull_of(cp, [K[i] for i in range(m) if i not in A1], sem, x, y, cache)
    sel2 = hull_of(cp, [K[i] for i in range(m) if i not in A2], sem, x, y, cache)
    # RHIG_1 / RHIG_2 = union of hulls over all (m-j)-subsets
    def rhig(j):
        lo, hi = np.inf, -np.inf
        for S in itertools.combinations(range(m), m - j):
            h = hull_of(cp, [K[i] for i in S], sem, x, y, cache)
            if h is None:
                continue
            lo, hi = min(lo, h[0]), max(hi, h[1])
        return None if lo == np.inf else (lo, hi)
    r1, r2 = rhig(1), rhig(2)
    blanket = tau_hull(cp, sem, x, y, cache)
    return dict(m=m, delta=delta, D1=D1, A1=A1, A2=A2, D2=D2,
                sel1=sel1, sel2=sel2, r1=r1, r2=r2, blanket=blanket,
                none_closure=none_closure, H1=H1, theta0=th0)


def w(iv):
    return None if iv is None else iv[1] - iv[0]


if __name__ == "__main__":
    NOBS = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
    NPROB = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    rng = np.random.default_rng(20260907)
    rows, tried, t0 = [], 0, time.time()
    while len(rows) < NPROB and tried < 60000 and time.time() - t0 < 420:
        tried += 1
        pr = build(rng)
        if pr is None:
            continue
        cache = {}
        try:
            a = analyse(pr, NOBS, cache)
        except Exception:
            continue
        if a["blanket"] is None:
            continue
        a["pr"] = pr
        rows.append(a)
    print(f"n_obs={NOBS}  problems={len(rows)}  tried={tried}  {time.time()-t0:.1f}s")
    print(f"mean |K|={np.mean([r['m'] for r in rows]):.2f}  "
          f"mean delta={np.mean([r['delta'] for r in rows]):.4f}")

    # ---------- 0. is the 'Delta_i := 0 if closure is None' branch ever taken? ----------
    nc = sum(r["none_closure"] for r in rows)
    print(f"\n[0] leave-one-out closures returning None: {nc} out of "
          f"{sum(r['m'] for r in rows)} -> the fail-OPEN default is "
          f"{'dead code' if nc == 0 else 'REACHABLE'}")

    # ---------- 1. how often does the screen fire, and how big is A ----------
    sizes = Counter(len(r["A1"]) for r in rows)
    print(f"\n[1] |A1| distribution: " +
          "  ".join(f"{k}:{100*v/len(rows):.1f}%" for k, v in sorted(sizes.items())))
    empty = 100 * np.mean([len(r["A1"]) == 0 for r in rows])
    print(f"    A1 empty in {empty:.1f}% of problems")
    print(f"    |A2| distribution: " +
          "  ".join(f"{k}:{100*v/len(rows):.1f}%"
                    for k, v in sorted(Counter(len(r['A2']) for r in rows).items())))

    # ---------- 2. bimodality: is SEL_A ever strictly between 0 and the blanket ----------
    S1 = np.array([w(r["sel1"]) for r in rows])
    B = np.array([w(r["blanket"]) for r in rows])
    R1 = np.array([w(r["r1"]) for r in rows])
    R2 = np.array([w(r["r2"]) if r["r2"] else np.nan for r in rows])
    S2 = np.array([w(r["sel2"]) for r in rows])
    eps = 1e-9
    zero = S1 < eps
    full = np.abs(S1 - B) < eps
    print(f"\n[2] SEL_A1 width: exactly 0 in {100*zero.mean():.1f}%, exactly the blanket in "
          f"{100*full.mean():.1f}%, strictly in between in "
          f"{100*np.mean(~zero & ~full):.1f}%")
    print(f"    (of the non-zero ones, {100*np.mean(full[~zero])if (~zero).sum() else float('nan'):.1f}% ARE the blanket)")
    print(f"    blanket itself is zero-width in {100*np.mean(B < eps):.1f}% of problems")

    # ---------- 3. the 2*delta sandwich ----------
    d = np.array([r["delta"] for r in rows])
    gain = R1 - S1                      # how much narrower SEL_A1 is than RHIG_1
    print(f"\n[3] width(RHIG_1) - width(SEL_A1): mean {gain.mean():+.4f}  "
          f"median {np.median(gain):+.4f}  max {gain.max():+.4f}")
    print(f"    violations of the bound gain <= 2*delta: "
          f"{int(np.sum(gain > 2*d + 1e-9))} of {len(rows)}   "
          f"(mean 2*delta = {2*d.mean():.4f})")
    print(f"    SEL_A1 STRICTLY WIDER than RHIG_1 in {100*np.mean(gain < -1e-9):.1f}% "
          f"(mean excess there {(-gain[gain < -1e-9]).mean() if np.any(gain < -1e-9) else 0:.4f})")
    g2 = R2 - S2
    ok = ~np.isnan(g2)
    print(f"    same for order 2: width(RHIG_2)-width(SEL_A2) mean {g2[ok].mean():+.4f}, "
          f"violations of 2*delta: {int(np.sum(g2[ok] > 2*d[ok] + 1e-9))}")
    nzB = B > eps
    print(f"    mean widths | blanket>0 (n={nzB.sum()}): SEL_A1={S1[nzB].mean():.4f} "
          f"RHIG_1={R1[nzB].mean():.4f} RHIG_2={R2[nzB].mean():.4f} blanket={B[nzB].mean():.4f}")
    print(f"    sum-ratio to blanket: SEL_A1={S1.sum()/B.sum():.3f}  RHIG_1={R1.sum()/B.sum():.3f}")

    # ---------- 4. one-sidedness: is 'Delta_max / width(RHIG_1) = 1.00' arithmetic? ----------
    onesided, ratios = [], []
    for r in rows:
        if r["r1"] is None or w(r["r1"]) < eps:
            continue
        lo, hi = r["r1"][0], r["r1"][1]
        th = r["theta0"]
        onesided.append(abs(lo - th) < eps or abs(hi - th) < eps)
        ratios.append(max(r["D1"]) / (hi - lo))
    onesided = np.array(onesided); ratios = np.array(ratios)
    print(f"\n[4] among problems with width(RHIG_1)>0 (n={len(ratios)}): "
          f"theta_0 is an ENDPOINT of the RHIG_1 hull in {100*onesided.mean():.1f}%")
    print(f"    median Delta_max/width(RHIG_1) = {np.median(ratios):.3f}; "
          f"= 1.000 exactly in {100*np.mean(np.abs(ratios-1) < 1e-9):.1f}% "
          f"-- and in {100*np.mean(np.abs(ratios[onesided]-1) < 1e-9) if onesided.sum() else float('nan'):.1f}% "
          f"of the one-sided ones")
    multi = np.array([len([i for i in r['A1']]) >= 2 for r in rows])
    print(f"    fraction with >=2 active claims: {100*multi.mean():.1f}% "
          f"(the ratio can only differ from 1 there or when two-sided)")

    # ---------- 5. non-amenable extensions inside the hull ----------
    tot = sum(r["blanket"][2] for r in rows)
    bad = sum(r["blanket"][3] for r in rows)
    print(f"\n[5] DAG extensions in the BLANKET hull whose O*(D) is NOT a valid "
          f"adjustment set in D: {bad}/{tot} = {100*bad/tot:.1f}%")
    frac = np.array([r["blanket"][3] / max(r["blanket"][2], 1) for r in rows])
    print(f"    problems with at least one such extension: {100*np.mean(frac>0):.1f}%; "
          f"mean share within a problem {100*frac.mean():.1f}%")

    # ---------- 6. the redundancy game: state every claim twice ----------
    print("\n[6] REDUNDANCY GAME -- the analyst restates each claim as two "
          "logically equivalent statements (same K, same G0, same Z)")
    fired1 = fired2 = 0; n_g = 0
    for r in rows:
        if len(r["A1"]) == 0:
            continue
        n_g += 1
        pr = r["pr"]; K2 = list(pr["K"]) + list(pr["K"])      # duplicate
        cache = {}
        m2 = len(K2)
        th0 = pr["theta0"]; delta = r["delta"]
        act1 = []
        for i in range(m2):
            h = hull_of(pr["cp"], [K2[j] for j in range(m2) if j != i],
                        pr["sem"], pr["x"], pr["y"], cache)
            if h is not None and max(abs(h[0]-th0), abs(h[1]-th0)) > delta:
                act1.append(i)
        act2 = set(act1)
        for i, j in itertools.combinations(range(m2), 2):
            h = hull_of(pr["cp"], [K2[t] for t in range(m2) if t not in (i, j)],
                        pr["sem"], pr["x"], pr["y"], cache)
            if h is not None and max(abs(h[0]-th0), abs(h[1]-th0)) > delta:
                act2 |= {i, j}
        fired1 += len(act1) > 0
        fired2 += len(act2) > 0
    print(f"    on the {n_g} problems where the screen DID fire on the elicited K:")
    print(f"      after duplication, order-1 screen still fires: {fired1}/{n_g} "
          f"= {100*fired1/max(n_g,1):.1f}%")
    print(f"      after duplication, order-2 screen still fires: {fired2}/{n_g} "
          f"= {100*fired2/max(n_g,1):.1f}%")
