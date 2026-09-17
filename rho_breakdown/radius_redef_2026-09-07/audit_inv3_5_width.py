#!/usr/bin/env python3
"""inverse-3 audit 5: THE TARGET METRIC.  Width required for honest 95% coverage.

Worlds are drawn from the SAME declared model nu(q) assumes: each claim independently
reversed w.p. q_true, restricted to admissible corruptions; then a DAG uniform in
that MPDAG's extensions; then a linear SEM.  The analyst reports the estimand
adjusted by Z = O*(G0).

Hedges compared (each an interval of candidate estimands, then SCALED to hit exactly
95% empirical coverage, then averaged -- lower is better):
  BLANKET   : range over the whole equivalence class of Chat
  NU(q_dec) : q-weighted central interval over the corruption-reachable graphs
  J<=1      : range over corruptions of at most ONE claim  (the truncation nu(q)
              already collapses to at small q; needs no q at all)
  POINT     : no hedge (the naive analyst)
"""
import sys, itertools, math
from collections import Counter
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_dag, random_sem,
                                    adjusted_estimand, is_valid_adjustment_set_mpdag)
from audit_inv3_4_clean import clean_problem
from audit_inv3_1_core import corruption_table, nu_profile, nu
from harness import ball, r_val

ALPHA = 0.05


def candidates(pr):
    """(J, G_J, [DAGs]) for every admissible corruption, plus the full class."""
    K, cp = pr["K"], pr["cpdag"]
    m = len(K)
    out = []
    for j in range(m + 1):
        for J in itertools.combinations(range(m), j):
            KJ = [(b, a) if i in J else (a, b) for i, (a, b) in enumerate(K)]
            g = apply_orientations(cp, KJ)
            if g is None:
                continue
            ext = enumerate_dag_extensions(g, limit=60)
            if ext:
                out.append((J, g, ext))
    full = enumerate_dag_extensions(cp, limit=400)
    return out, full


def weighted_interval(vals, ws, alpha):
    o = np.argsort(vals)
    v, w = np.array(vals)[o], np.array(ws)[o]
    w = w / w.sum()
    c = np.cumsum(w)
    lo = v[np.searchsorted(c, alpha / 2)]
    hi = v[min(np.searchsorted(c, 1 - alpha / 2), len(v) - 1)]
    return lo, hi


def run(q_true, q_dec, n_problems=220, draws=6, seed=20260907):
    rng = np.random.default_rng(seed)
    rows = []
    got = tried = 0
    while got < n_problems and tried < 200000:
        tried += 1
        pr = clean_problem(rng)
        if pr is None:
            continue
        cands, full = candidates(pr)
        if len(cands) < 2 or not full:
            continue
        got += 1
        m = len(pr["K"])
        x, y, Z = pr["x"], pr["y"], pr["z"]
        _, rws = corruption_table(pr)
        nu_j, _ = nu_profile(m, rws)
        rv = r_val(ball(pr))
        for _ in range(draws):
            # ---- draw a world from the declared model -----------------------
            for _try in range(200):
                J = tuple(i for i in range(m) if rng.random() < q_true)
                hit = [c for c in cands if c[0] == J]
                if hit:
                    break
            else:
                continue
            _, gJ, extJ = hit[0]
            d = extJ[rng.integers(len(extJ))]
            sem = random_sem(d, rng)
            truth = adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(d, x, y))
            point = adjusted_estimand(sem, x, y, Z)
            # ---- the hedges, all under the SAME sem -------------------------
            blanket = [adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(dd, x, y))
                       for dd in full]
            vals, ws, vals1 = [], [], []
            for Jc, gc, extc in cands:
                w = q_dec**len(Jc) * (1 - q_dec)**(m - len(Jc)) / len(extc)
                for dd in extc:
                    e = adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(dd, x, y))
                    vals.append(e); ws.append(w)
                    if len(Jc) <= 1:
                        vals1.append(e)
            nlo, nhi = weighted_interval(vals, ws, ALPHA)
            rows.append(dict(truth=truth, point=point, rval=rv, nu=nu(q_dec, m, nu_j),
                             blanket=(min(blanket), max(blanket)),
                             nuiv=(nlo, nhi), j1=(min(vals1), max(vals1)),
                             wsum=sum(q_dec**len(Jc)*(1-q_dec)**(m-len(Jc)) for Jc, _, _ in cands)))
    return rows


def calibrate(rows, key, target=0.95):
    """Smallest scale c on the half-width that reaches `target` coverage; mean width."""
    def cov(c):
        k = 0
        for r in rows:
            lo, hi = r[key]
            mid, half = (lo + hi) / 2, (hi - lo) / 2
            k += (mid - c*half - 1e-12) <= r["truth"] <= (mid + c*half + 1e-12)
        return k / len(rows)
    if cov(1e6) < target:
        return None, None, cov(1.0)
    lo, hi = 0.0, 1.0
    while cov(hi) < target:
        hi *= 2
        if hi > 1e7:
            return None, None, cov(1.0)
    for _ in range(60):
        mid = (lo + hi) / 2
        if cov(mid) >= target:
            hi = mid
        else:
            lo = mid
    c = hi
    w = np.mean([c * (r[key][1] - r[key][0]) for r in rows])
    return c, w, cov(1.0)


if __name__ == "__main__":
    for q_true, q_dec, tag in [(0.20, 0.20, "q declarado = q verdadeiro = 0.20  (o oraculo do autor)"),
                               (0.20, 0.05, "q SUBdeclarado: verdadeiro 0.20, declarado 0.05"),
                               (0.05, 0.20, "q SOBREdeclarado: verdadeiro 0.05, declarado 0.20")]:
        rows = run(q_true, q_dec)
        print(f"\n=== {tag}   (n = {len(rows)} mundos) ===")
        print(f"{'metodo':<10} {'cobertura crua':>15} {'escala p/ 95%':>14} {'LARGURA MEDIA':>15}")
        for key, name in [("blanket", "BLANKET"), ("nuiv", "NU(q)"), ("j1", "J<=1")]:
            c, w, raw = calibrate(rows, key)
            print(f"{name:<10} {raw:>15.4f} {('inf' if c is None else f'{c:.4f}'):>14} "
                  f"{('inf' if w is None else f'{w:.4f}'):>15}")
        S = [r for r in rows if r["rval"] == 1]
        if len(S) > 30:
            print(f"  -- estrato r_val = 1 (n = {len(S)}) --")
            for key, name in [("blanket", "BLANKET"), ("nuiv", "NU(q)"), ("j1", "J<=1")]:
                c, w, raw = calibrate(S, key)
                print(f"  {name:<8} {raw:>15.4f} {('inf' if c is None else f'{c:.4f}'):>14} "
                      f"{('inf' if w is None else f'{w:.4f}'):>15}")
        print(f"  massa total dos pesos do intervalo (deveria ser 1 se fosse a prob. declarada): "
              f"media {np.mean([r['wsum'] for r in rows]):.4f} min {min(r['wsum'] for r in rows):.4f}")
