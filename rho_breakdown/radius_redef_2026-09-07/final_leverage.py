#!/usr/bin/env python3
"""FINAL definition check: claim leverage lambda + the flip hedge.

lambda_k = max over DAG-extensions of cl(K[k<-rev]) of |tau_D(x,y) - theta_hat|,
with lambda_k := 0 when the reversal is inconsistent (k in R).

Checks, in order:
  (S)  screen soundness   : k not in L_flip  =>  lambda_k == 0        (counterexamples)
  (D)  degeneracy profile : atom mass of lambda_max/|theta_hat|, rankable fraction
  (A)  continuity payoff  : AUC of lambda_max for "Z invalid in the true DAG"
  (W)  width rule         : hull +/- c*se, c calibrated to 95% coverage, mean width
"""
import sys, itertools
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import random_dag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (optimal_adjustment_set_mpdag,
                                    optimal_adjustment_set_dag,
                                    is_valid_adjustment_set_mpdag,
                                    is_valid_adjustment_set_dag,
                                    random_sem, adjusted_estimand,
                                    asymptotic_variance)
from audit_wit1_2_coverage import elicit

N_SAMPLE = 2000          # nominal sample size behind the standard error


def tau_range(sem, x, y, g):
    """Every total effect realisable in the class g (over its DAG extensions)."""
    vals = []
    for d in enumerate_dag_extensions(g):
        zz = optimal_adjustment_set_dag(d, x, y)
        vals.append(adjusted_estimand(sem, x, y, zz))
    return vals


def problem(rng, n_nodes, k_claims, budget):
    dag = random_dag(rng, n_nodes, 0.32)
    if not dag.directed_edges:
        return None
    cpdag = dag_to_cpdag(dag)
    und = sorted(tuple(sorted(e)) for e in cpdag.undirected_edges)
    if not (2 <= len(und) <= 6):
        return None
    K_true = elicit(rng, dag, cpdag, k_claims)
    if K_true is None or len(K_true) < budget:
        return None
    fi = set(rng.permutation(len(K_true))[:budget].tolist())
    K = [((k[1], k[0]) if i in fi else k) for i, k in enumerate(K_true)]
    g0 = apply_orientations(cpdag, K)
    if g0 is None:
        return None
    sem = random_sem(dag, rng)
    nodes = sorted(dag.nodes)
    cand = [(a, b) for a in nodes for b in nodes if a != b]
    rng.shuffle(cand)
    for x, y in cand:
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None:
            continue
        z = frozenset(z)
        truth = adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(dag, x, y))
        if abs(truth) < 0.05:
            continue
        theta = adjusted_estimand(sem, x, y, z)
        se = float(np.sqrt(asymptotic_variance(sem, x, y, z) / N_SAMPLE))
        return dict(dag=dag, cpdag=cpdag, K=K, K_true=K_true, fi=fi, g0=g0,
                    x=x, y=y, z=z, sem=sem, truth=truth, theta=theta, se=se,
                    z_valid=is_valid_adjustment_set_dag(dag, x, y, z))
    return None


def leverage(p, budget_hedge=1):
    """lambda profile at budget 1 + the hedge hull at budget `budget_hedge`."""
    cp, K, x, y, z, sem, th = p["cpdag"], p["K"], p["x"], p["y"], p["z"], p["sem"], p["theta"]
    lam, inR, inLflip, viol = {}, [], [], 0
    hull1 = [th]
    for k in K:
        rest = [e for e in K if e != k]
        gk = apply_orientations(cp, rest + [(k[1], k[0])])
        if gk is None:
            inR.append(k); lam[k] = 0.0; continue
        valid = is_valid_adjustment_set_mpdag(gk, x, y, z)
        vals = tau_range(sem, x, y, gk)
        lam[k] = max(abs(v - th) for v in vals) if vals else 0.0
        hull1.extend(vals)
        if not valid:
            inLflip.append(k)
        elif lam[k] > 1e-9:
            viol += 1                      # screen soundness counterexample
    hull = list(hull1)
    if budget_hedge >= 2:
        for S in itertools.combinations(K, 2):
            rest = [e for e in K if e not in S]
            gS = apply_orientations(cp, rest + [(a[1], a[0]) for a in S])
            if gS is None:
                continue
            hull.extend(tau_range(sem, x, y, gS))
    return lam, inR, inLflip, viol, (min(hull), max(hull))


def calibrate(lo, hi, th, se, truth, target=0.95):
    """Smallest c >= 0 with empirical coverage >= target for [lo-c*se, hi+c*se]."""
    need = np.maximum((lo - truth) / se, (truth - hi) / se)   # c needed per problem
    need = np.maximum(need, 0.0)
    c = float(np.quantile(need, target))
    cov = float(np.mean(truth <= hi + c * se) * 1.0) if False else \
          float(np.mean((truth >= lo - c * se) & (truth <= hi + c * se)))
    width = float(np.mean((hi - lo) + 2 * c * se))
    return c, cov, width


def run(n_nodes=7, k_claims=3, budget=1, n=400, seed=7):
    rng = np.random.default_rng(seed)
    rows, tried = [], 0
    while len(rows) < n and tried < 400000:
        tried += 1
        p = problem(rng, n_nodes, k_claims, budget)
        if p is None:
            continue
        lam, inR, inLf, viol, (lo1, hi1) = leverage(p, 1)
        _, _, _, _, (lo2, hi2) = leverage(p, 2)
        cls = tau_range(p["sem"], p["x"], p["y"], p["cpdag"])
        vals = sorted(lam.values(), reverse=True)
        rows.append(dict(theta=p["theta"], truth=p["truth"], se=p["se"],
                         z_valid=p["z_valid"], m=len(p["K"]),
                         nR=len(inR), nLf=len(inLf), viol=viol,
                         lmax=vals[0] if vals else 0.0,
                         l2=vals[1] if len(vals) > 1 else 0.0,
                         npos=sum(1 for v in lam.values() if v > 1e-9),
                         lo1=lo1, hi1=hi1, lo2=lo2, hi2=hi2,
                         locls=min(cls + [p["theta"]]), hicls=max(cls + [p["theta"]])))
    return rows


def report(rows, tag):
    A = lambda k: np.array([r[k] for r in rows], dtype=float)
    n = len(rows)
    th, truth, se = A("theta"), A("truth"), A("se")
    lmax, l2, npos = A("lmax"), A("l2"), A("npos")
    print(f"\n================ {tag}  n={n} ================")
    print(f"mean |K|={A('m').mean():.2f}  mean |R|={A('nR').mean():.2f}  "
          f"mean |L_flip|={A('nLf').mean():.2f}  Z valid in truth: {100*A('z_valid').mean():.1f}%")
    print(f"(S) screen-soundness counterexamples (k not in L_flip but lambda_k>0): "
          f"{int(A('viol').sum())} over {int(A('m').sum())} claims")

    rel = lmax / np.abs(th)
    live = rel > 1e-9
    print(f"(D) lambda_max == 0 on {100*np.mean(~live):.1f}%   live stratum n={live.sum()}")
    if live.sum():
        q = np.round(rel[live], 2)
        vals, cnt = np.unique(q, return_counts=True)
        print(f"    largest atom of lambda_max/|theta| (2 dp) on the live stratum: "
              f"{100*cnt.max()/live.sum():.1f}% at value {vals[cnt.argmax()]:.2f}")
        print(f"    quantiles of lambda_max/|theta| on live: "
              f"p10={np.quantile(rel[live],.10):.3f} p50={np.quantile(rel[live],.50):.3f} "
              f"p90={np.quantile(rel[live],.90):.3f} max={rel[live].max():.2f}")
    rk = npos >= 2
    sep = (lmax > 1.1 * l2) & rk
    print(f"    rankable (>=2 claims with lambda>0): {100*np.mean(rk):.1f}%   "
          f"of those, top-1 separated by >10%: {100*(sep.sum()/max(rk.sum(),1)):.1f}%")

    lab = (~A("z_valid").astype(bool)).astype(float)
    if 0 < lab.sum() < n:
        order = np.argsort(lmax)
        r = np.empty(n); r[order] = np.arange(1, n + 1)
        # average ranks for ties
        for v in np.unique(lmax):
            m_ = lmax == v
            r[m_] = r[m_].mean()
        n1, n0 = lab.sum(), n - lab.sum()
        auc = (r[lab == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
        bin_ = (A("nLf") > 0).astype(float)
        rb = np.empty(n); rb[np.argsort(bin_)] = np.arange(1, n + 1)
        for v in np.unique(bin_):
            m_ = bin_ == v
            rb[m_] = rb[m_].mean()
        aucb = (rb[lab == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
        print(f"(A) AUC for 'Z invalid in the true DAG':  lambda_max {auc:.3f}   "
              f"vs binary 1[L_flip nonempty] {aucb:.3f}   (prevalence {100*lab.mean():.1f}%)")

    print(f"(W) width rule: [lo - c*se, hi + c*se], c calibrated to 95% coverage, N={N_SAMPLE}")
    print(f"    {'method':26s} {'c':>8s} {'coverage':>9s} {'mean width':>11s} {'hull width':>11s}")
    for nm, lo, hi in [("point (no hedge)", th, th),
                       ("flip hedge, budget 1", A("lo1"), A("hi1")),
                       ("flip hedge, budget 2", A("lo2"), A("hi2")),
                       ("blanket class hedge", A("locls"), A("hicls"))]:
        c, cov, w = calibrate(lo, hi, th, se, truth)
        print(f"    {nm:26s} {c:8.2f} {100*cov:8.1f}% {w:11.4f} {np.mean(hi-lo):11.4f}")


if __name__ == "__main__":
    for b in (1, 2):
        report(run(budget=b, n=400, seed=7 + b), f"7-node, truthful-then-{b}-reversed")
