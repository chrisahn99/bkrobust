"""Measured hedge width: certified-core hedge vs blanket hedge over the class."""
from __future__ import annotations
import itertools, json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import make_problem, Problem, minimal_unsafe, min_hitting_sets, powerset
from bkrobust.demo.meek import enumerate_dag_extensions, is_consistent_extension
from bkrobust.demo.evaluate import (
    random_sem, adjusted_estimand, optimal_adjustment_set_dag,
    is_valid_adjustment_set_mpdag,
)

def estimands(sem, g, x, y, cache):
    key = (g.directed_edges, g.undirected_edges)
    if key not in cache:
        vals = []
        for d in enumerate_dag_extensions(g):
            Zd = optimal_adjustment_set_dag(d, x, y)
            vals.append(adjusted_estimand(sem, x, y, Zd))
        cache[key] = vals
    return cache[key]

def run(seed, n_problems, kmin, kmax, tag):
    rng = np.random.default_rng(seed)
    rows = []
    t0 = time.time()
    tries = 0
    while len(rows) < n_problems and tries < 400000 and time.time() - t0 < 1500:
        tries += 1
        P = make_problem(rng, kmin, kmax)
        if P is None:
            continue
        cpdag, Kc, x, y, Z, dag = P["cpdag"], P["Kc"], P["x"], P["y"], P["Z"], P["dag"]
        m = len(Kc)
        prob = Problem(cpdag, Kc, x, y, Z)
        mins = minimal_unsafe(prob, m)
        c, hs = min_hitting_sets(mins, m)
        core = set(hs[0])
        sem = random_sem(dag, rng)
        theta_true = sem.true_total_effect(x, y)
        theta_point = adjusted_estimand(sem, x, y, Z)
        cache = {}
        TM, TW = [theta_point], [theta_point]
        nW = nM = 0
        for J in powerset(m):
            g = prob.world(J)
            if g is None:
                continue
            nW += 1
            vals = estimands(sem, g, x, y, cache)
            TW.extend(vals)
            if set(J) & core or not J:
                nM += 1
                TM.extend(vals)
        TC = estimands(sem, cpdag, x, y, cache)
        wC = max(TC) - min(TC)
        wM = max(TM) - min(TM)
        wW = max(TW) - min(TW)
        # structural coverage of the true effect
        tol = 1e-8 + 1e-6 * abs(theta_true)
        covM = (min(TM) - tol) <= theta_true <= (max(TM) + tol)
        covW = (min(TW) - tol) <= theta_true <= (max(TW) + tol)
        covPoint = abs(theta_point - theta_true) <= tol
        rows.append(dict(m=m, c=c, nW=nW, nM=nM, wM=wM, wW=wW, wC=wC,
                         nC=len(TC),
                         covM=bool(covM), covW=bool(covW), covPoint=bool(covPoint),
                         scale=abs(theta_true)))
    out = dict(tag=tag, seed=seed, secs=round(time.time()-t0,1), tries=tries, rows=rows)
    with open(f"width_{tag}.json","w") as f:
        json.dump(out, f)
    print(tag, "n=", len(rows), "secs=", out["secs"])

if __name__ == "__main__":
    run(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[1])
