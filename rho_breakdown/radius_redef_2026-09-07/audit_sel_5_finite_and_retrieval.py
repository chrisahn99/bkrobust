"""SEL_A audit 5 -- (a) the screen run on ESTIMATED moments, which is the only way
it can ever be run, and (b) the retrieval claim on a bigger sample.
"""
import sys, itertools, time
from dataclasses import dataclass
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_dag, optimal_adjustment_set_mpdag,
                                    is_valid_adjustment_set_dag, adjusted_estimand,
                                    asymptotic_variance)
from audit_sel_2_coverage import build, Z196

Z = Z196


class EmpiricalSEM:
    """Same interface adjusted_estimand/asymptotic_variance need, empirical Sigma."""
    def __init__(self, dag, sigma_hat):
        self.dag = dag
        self._s = sigma_hat

    def covariance(self):
        return self._s


def sample_sigma(sem, n, rng):
    nodes = list(sem.dag.nodes)
    idx = {v: i for i, v in enumerate(nodes)}
    data = np.zeros((n, len(nodes)))
    order = []
    seen = set()

    def visit(v):
        if v in seen:
            return
        seen.add(v)
        for p in sorted(sem.dag.parents(v)):
            visit(p)
        order.append(v)
    for v in sorted(nodes):
        visit(v)
    for v in order:
        col = rng.normal(0.0, np.sqrt(sem.noise_var[v]), n)
        for p in sorted(sem.dag.parents(v)):
            col = col + sem.weights[(p, v)] * data[:, idx[p]]
        data[:, idx[v]] = col
    return np.cov(data, rowvar=False, bias=False)


def hull(cp, claims, esem, x, y, cache):
    g = apply_orientations(cp, list(claims))
    if g is None:
        return None
    key = g.edge_string()
    if key in cache:
        return cache[key]
    vals = []
    for D in enumerate_dag_extensions(g):
        try:
            vals.append(adjusted_estimand(esem, x, y, optimal_adjustment_set_dag(D, x, y)))
        except Exception:
            pass
    out = (min(vals), max(vals)) if vals else None
    cache[key] = out
    return out


def one(pr, nobs, rng):
    """Everything computed from an n-sample covariance, as in a real analysis."""
    sig = sample_sigma(pr["sem"], nobs, rng)
    es = EmpiricalSEM(pr["dag"], sig)
    x, y, cp, K = pr["x"], pr["y"], pr["cp"], pr["K"]
    m = len(K)
    th0 = adjusted_estimand(es, x, y, pr["z"])
    se = np.sqrt(asymptotic_variance(es, x, y, pr["z"]) / nobs)
    delta = Z * se
    cache = {}
    D1 = []
    for i in range(m):
        h = hull(cp, [K[j] for j in range(m) if j != i], es, x, y, cache)
        D1.append(0.0 if h is None else max(abs(h[0]-th0), abs(h[1]-th0)))
    D1 = np.array(D1)
    A1 = [i for i in range(m) if D1[i] > delta]
    sel = hull(cp, [K[i] for i in range(m) if i not in A1], es, x, y, cache) or (th0, th0)
    lo, hi = np.inf, -np.inf
    for S in itertools.combinations(range(m), m-1):
        h = hull(cp, [K[i] for i in S], es, x, y, cache)
        if h:
            lo, hi = min(lo, h[0]), max(hi, h[1])
    r1 = (th0, th0) if lo == np.inf else (lo, hi)
    return dict(th0=th0, se=se, sel=sel, r1=r1, A1=A1, D1=D1, delta=delta)


def cover(iv, se, tau):
    return iv[0] - Z*se <= tau <= iv[1] + Z*se


if __name__ == "__main__":
    from audit_sel_3_fair import amenable
    print("[a] THE SCREEN ON ESTIMATED MOMENTS (one false claim, amenable problems)")
    for nobs in (200, 2000, 20000):
        rng = np.random.default_rng(5150)
        cs, cr, cp_, ws, wr, k, fire = 0, 0, 0, 0.0, 0.0, 0, 0
        t0 = time.time()
        while k < 200 and time.time()-t0 < 180:
            pr = build(rng, 1)
            if pr is None or not amenable(pr):
                continue
            try:
                o = one(pr, nobs, rng)
            except Exception:
                continue
            k += 1
            tau = pr["tau_true"]
            cs += cover(o["sel"], o["se"], tau)
            cr += cover(o["r1"], o["se"], tau)
            cp_ += cover((o["th0"], o["th0"]), o["se"], tau)
            ws += (o["sel"][1]-o["sel"][0]) + 2*Z*o["se"]
            wr += (o["r1"][1]-o["r1"][0]) + 2*Z*o["se"]
            fire += len(o["A1"]) > 0
        print(f"   n={nobs:>6}  problems={k}   screen fires {100*fire/k:.1f}%")
        print(f"      coverage: plain CI {cp_/k:.3f} | SEL_A1 {cs/k:.3f} | RHIG_1 {cr/k:.3f}")
        print(f"      mean width:          SEL_A1 {ws/k:.4f} | RHIG_1 {wr/k:.4f}  "
              f"(difference {100*(wr-ws)/max(wr,1e-12):+.2f}%)")

    print("\n[b] RETRIEVAL, bigger sample: does the Delta ranking find the FALSE claim?")
    from audit_sel_2_coverage import methods
    for n_false in (1, 2):
        rng = np.random.default_rng(909 + n_false)
        n_dam = 0; p1 = 0; base = 0; rand = 0.0; rec = 0; toptrue = 0; nfire = 0
        t0 = time.time()
        while n_dam < 200 and time.time()-t0 < 240:
            pr = build(rng, n_false)
            if pr is None or not amenable(pr):
                continue
            try:
                iv, se, delta, D1, D2, A1, A2 = methods(pr, 10000, {})
            except Exception:
                continue
            D1 = np.asarray(D1, float)
            if D1.max() > 0:
                nfire += 1
                toptrue += int(np.argmax(D1 + rng.random(len(D1))*1e-12)) not in pr["false_idx"]
            if abs(pr["tau_true"] - pr["theta0"]) <= 1e-9:
                continue
            n_dam += 1
            m = len(pr["K"])
            F = pr["false_idx"]
            p1 += int(np.argmax(D1 + rng.random(m)*1e-12)) in F
            touch = np.array([1.0 if (pr["x"] in e or pr["y"] in e) else 0.0 for e in pr["K"]])
            base += int(np.argmax(touch + rng.random(m)*1e-9)) in F
            rand += len(F)/m
            rec += any(D1[i] > 0 for i in F)
        print(f"   {n_false} false claim(s), {n_dam} damaging amenable problems "
              f"({time.time()-t0:.0f}s)")
        print(f"      precision@1: Delta {p1/n_dam:.3f} | 'touches X or Y' {base/n_dam:.3f} "
              f"| random {rand/n_dam:.3f}")
        print(f"      recall (some false claim has Delta>0): {rec/n_dam:.3f}")
        print(f"      among the {nfire} problems where the screen fires, the TOP-ranked "
              f"claim is a TRUE claim: {100*toptrue/max(nfire,1):.1f}%")
