"""Where does the estimand spread live: core-touching worlds or certified worlds?"""
from __future__ import annotations
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import make_problem, Problem, minimal_unsafe, min_hitting_sets, powerset
from bkrobust.demo.meek import enumerate_dag_extensions
from bkrobust.demo.evaluate import (random_sem, adjusted_estimand,
                                    optimal_adjustment_set_dag, optimal_adjustment_set_mpdag)

def run(seed, n_problems, kmin, kmax, tag):
    rng = np.random.default_rng(seed); rows = []; t0=time.time(); tries=0
    while len(rows) < n_problems and tries < 400000 and time.time()-t0 < 1200:
        tries += 1
        P = make_problem(rng, kmin, kmax)
        if P is None: continue
        cpdag, Kc, x, y, Z, dag = P["cpdag"], P["Kc"], P["x"], P["y"], P["Z"], P["dag"]
        m = len(Kc); prob = Problem(cpdag, Kc, x, y, Z)
        mins = minimal_unsafe(prob, m); c, hs = min_hitting_sets(mins, m); core = set(hs[0])
        sem = random_sem(dag, rng)
        th = sem.true_total_effect(x, y)
        pt = adjusted_estimand(sem, x, y, Z)
        F, M = [], [pt]
        nF = 0
        for J in powerset(m):
            g = prob.world(J)
            if g is None: continue
            vals = [adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(d, x, y))
                    for d in enumerate_dag_extensions(g)]
            if J and not (set(J) & core):
                F.extend(vals); nF += 1
            else:
                M.extend(vals)
        aM, bM = min(M), max(M)
        row = dict(m=m, c=c, nF=nF, th=th, pt=pt, aM=aM, bM=bM, wM=bM-aM)
        if F:
            aF, bF = min(F), max(F)
            row.update(aF=aF, bF=bF, wF=bF-aF,
                       inside=bool(aF >= aM - 1e-9 and bF <= bM + 1e-9),
                       excess_lo=max(0.0, aM-aF), excess_hi=max(0.0, bF-bM),
                       allTrue=bool(max(abs(v-th) for v in F) < 1e-9))
        rows.append(row)
    json.dump(dict(tag=tag, rows=rows, secs=round(time.time()-t0,1)),
              open(f"decomp_{tag}.json","w"))
    R = rows; n=len(R)
    F = [r for r in R if "wF" in r]
    print(f"{tag}: n={n} |K*|={np.mean([r['m'] for r in R]):.2f} c={np.mean([r['c'] for r in R]):.2f} secs={time.time()-t0:.1f}")
    print(f"  problems with >=1 certified non-trivial world: {len(F)}")
    if F:
        print(f"  certified-world estimands ALL equal theta_true on {np.mean([r['allTrue'] for r in F]):.4f}")
        print(f"  certified range inside must-hedge range on {np.mean([r['inside'] for r in F]):.4f}")
        print(f"  mean certified spread wF={np.mean([r['wF'] for r in F]):.4f} vs must-hedge wM={np.mean([r['wM'] for r in F]):.4f}")
        ex=[max(r['excess_lo'],r['excess_hi']) for r in F]
        print(f"  mean excess beyond must-hedge range={np.mean(ex):.6f} max={np.max(ex):.6f}")
if __name__=="__main__":
    run(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[1])
