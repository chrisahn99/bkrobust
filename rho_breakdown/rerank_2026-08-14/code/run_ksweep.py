"""
E1 STEP 3: the |K| sweep the pilot never ran.

MAX_K = 4 is hard-coded at run_linear.py:34, so rho* at |K| = 8 does NOT exist on
disk and cannot be recomputed from it. This script re-simulates it, with a
PAIRED control that is stronger than seed-matching:

  For each SCM we draw ONE uniform random ordering of the orientable statements
  K_all (the undirected edges of the CPDAG, oriented as in D*). Then
      K4 = ordering[:4]   K6 = ordering[:6]   K8 = ordering[:8]
  so K4 subset K6 subset K8 on an IDENTICAL graph, identical (x,y), identical
  Sigma, identical tau. K4 is still a uniform random 4-subset of K_all, i.e.
  distributionally identical to the original protocol, but now perfectly paired.
  This sidesteps the RNG-stream shift that makes a naive MAX_K change produce a
  different ensemble (MAX_K sits upstream of a conditional rng.choice at
  run_linear.py:65-67).

Everything downstream -- apply_background_knowledge, optimal_adjustment_set,
is_valid_adjustment_set, ols_coefficient, cov_linear, total_effect_linear,
make_linear_iscm -- is the pilot's validated machinery, imported unmodified.

Usage: python run_ksweep.py <n_scm> <ensemble> <out.json>
  ensemble = 'licensed'  : p in 5..10, deg in {1.5,2,2.5,4,6}  (claim's scope)
             'large'     : p in 15..25, deg in {1.5,2,2.5}     (outside scope)
"""
import json
import sys
from itertools import combinations
from multiprocessing import Pool

import numpy as np

from adjust import (cov_linear, is_amenable, is_valid_adjustment_set,
                    optimal_adjustment_set, ols_coefficient,
                    possibly_causal_paths, total_effect_linear)
from graphs import (MeekFail, apply_background_knowledge, dag_agrees_with,
                    dag_to_cpdag, directed_edges, random_dag, undirected_edges)
from scm import make_linear_iscm

K_SIZES = (4, 6, 8)
MAX_RHO = 3
NEED_U = max(K_SIZES)


def analyse_one(args):
    seed, p, deg = args
    rng = np.random.default_rng(seed)
    D = random_dag(p, deg, rng)
    C = dag_to_cpdag(D)
    U = undirected_edges(C)
    if len(U) < NEED_U:
        return None                      # rejection: cannot support |K| = 8

    pairs = [(x, y) for x in range(p) for y in range(p)
             if x != y and len(possibly_causal_paths(D, x, y)) > 0]
    if not pairs:
        return None
    x, y = pairs[int(rng.integers(len(pairs)))]

    scm = make_linear_iscm(D, rng)
    Sigma = cov_linear(scm["A"], scm["omega"])
    tau = total_effect_linear(scm["A"], x, y)
    if abs(tau) < 1e-3:
        return None

    # ONE ordering -> nested K4 subset K6 subset K8
    K_all = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U]
    order = rng.permutation(len(K_all))
    K_all = [K_all[i] for i in order]

    out = dict(seed=seed, p=p, deg=deg, x=x, y=y, tau=tau,
               n_undirected=len(U), n_edges=int(D.sum()),
               cpdag_amenable=bool(is_amenable(C, x, y)), arms={})

    n_cpdag_dir = len(directed_edges(C))
    for k in K_SIZES:
        K = K_all[:k]
        try:
            G0 = apply_background_knowledge(C, K)
        except MeekFail:
            continue                     # cannot happen for true K, but be safe
        O0 = optimal_adjustment_set(G0, x, y)
        if O0 is None:
            out["arms"][f"K{k}"] = dict(n_K=k, mpdag_amenable=False, members=[])
            continue
        est0 = ols_coefficient(Sigma, x, y, O0)

        members = []
        for rho in range(1, MAX_RHO + 1):
            for flip in combinations(range(k), rho):
                Kp = [(v, u) if idx in flip else (u, v)
                      for idx, (u, v) in enumerate(K)]
                rec = dict(rho=rho)
                try:
                    G = apply_background_knowledge(C, Kp)
                except MeekFail:
                    rec.update(consistent=False)
                    members.append(rec)
                    continue
                O = optimal_adjustment_set(G, x, y)
                rec.update(consistent=True,
                           expels_true_dag=not dag_agrees_with(D, G),
                           cascade=len(directed_edges(G)) - n_cpdag_dir - len(Kp),
                           amenable=O is not None)
                if O is not None:
                    rec["ostar_changed"] = (O != O0)
                    rec["ostar_valid"] = bool(is_valid_adjustment_set(D, x, y, O))
                    rec["est"] = ols_coefficient(Sigma, x, y, O)
                    rec["bias"] = rec["est"] - tau
                members.append(rec)

        out["arms"][f"K{k}"] = dict(K=[list(e) for e in K], n_K=k, O0=sorted(O0),
                                    est0=est0, mpdag_amenable=True, members=members)
    if not out["arms"]:
        return None
    return out


ENSEMBLES = {
    # exactly the pilot's primary-arm grid (run_linear.py:176-177), rejection-
    # sampled to n_undirected >= 8. The scope-faithful test.
    "original": (range(5, 9), [1.5, 2.0, 2.5]),
    "licensed": (range(5, 11), [1.5, 2.0, 2.5, 4.0, 6.0]),
    "large":    (range(15, 26), [1.5, 2.0, 2.5]),
}
OVERSAMPLE = {"original": 400, "licensed": 120, "large": 12}


def main():
    n_scm = int(sys.argv[1]) if len(sys.argv) > 1 else 600
    ens = sys.argv[2] if len(sys.argv) > 2 else "licensed"
    out_path = sys.argv[3] if len(sys.argv) > 3 else f"../results/ksweep_{ens}.json"
    ps, degs = ENSEMBLES[ens]
    ps = list(ps)

    # oversample: rejection rate is high (n_undirected >= 8 is rare on small graphs)
    over = OVERSAMPLE[ens]
    rng = np.random.default_rng(20260814)
    jobs = [(int(s), int(rng.choice(ps)), float(rng.choice(degs)))
            for s in range(n_scm * over)]

    # ORDERED imap so the kept set is the first n_scm passing seeds -> reproducible
    # (imap_unordered would make the kept set scheduler-dependent).
    with Pool(8) as pool:
        res = []
        it = pool.imap(analyse_one, jobs, chunksize=8)
        for r in it:
            if r:
                res.append(r)
                if len(res) >= n_scm:
                    break
        pool.terminate()
    with open(out_path, "w") as f:
        json.dump(res, f)
    print(f"[{ens}] kept {len(res)} SCMs from <= {len(jobs)} draws -> {out_path}")


if __name__ == "__main__":
    main()
