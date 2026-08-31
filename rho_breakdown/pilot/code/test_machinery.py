"""
Validation suite. A wrong Meek closure invalidates everything downstream, so the
implementation is checked against BRUTE-FORCE MEC enumeration, which is the
definition, not a re-implementation.

Run:  python test_machinery.py
"""
import sys
from itertools import combinations

import numpy as np

from graphs import (MeekFail, apply_background_knowledge, common_orientation_graph,
                    consistent_dag_extensions, dag_agrees_with, dag_to_cpdag,
                    directed_edges, is_consistent, meek_closure, random_dag,
                    skeleton, undirected_edges, v_structures)
from adjust import (causal_nodes, is_amenable, is_valid_adjustment_set,
                    optimal_adjustment_set, ols_coefficient, cov_linear,
                    total_effect_linear, possibly_causal_paths)
from scm import make_linear_iscm, sample_linear, sample_linear_do

FAILS = []


def check(name, cond, info=""):
    if cond:
        print(f"  PASS  {name}")
    else:
        print(f"  FAIL  {name}  {info}")
        FAILS.append(name)


# --------------------------------------------------------------- textbook Meek examples
def test_meek_textbook():
    print("\n[T0] Meek rules on textbook examples")
    # R1: 0 -> 1 -- 2, 0 not adj 2  =>  1 -> 2
    G = np.zeros((3, 3), dtype=np.int8)
    G[0, 1] = 1
    G[1, 2] = G[2, 1] = 1
    H = meek_closure(G)
    check("R1 orients 1->2", H[1, 2] == 1 and H[2, 1] == 0)

    # R2: 0 -> 2 -> 1, 0 -- 1  =>  0 -> 1
    G = np.zeros((3, 3), dtype=np.int8)
    G[0, 2] = 1
    G[2, 1] = 1
    G[0, 1] = G[1, 0] = 1
    H = meek_closure(G)
    check("R2 orients 0->1", H[0, 1] == 1 and H[1, 0] == 0)

    # R3: a=0,b=1,c=2,d=3 ; 0--1, 0--2, 0--3, 2->1, 3->1, 2 not adj 3 => 0->1
    G = np.zeros((4, 4), dtype=np.int8)
    for u, v in [(0, 1), (0, 2), (0, 3)]:
        G[u, v] = G[v, u] = 1
    G[2, 1] = 1
    G[3, 1] = 1
    H = meek_closure(G)
    check("R3 orients 0->1", H[0, 1] == 1 and H[1, 0] == 0)

    # R4: a=0,b=1,c=2,d=3 ; 0--1, 0--2, 0--3, 3->2, 2->1, 1 not adj 3 => 0->1
    G = np.zeros((4, 4), dtype=np.int8)
    for u, v in [(0, 1), (0, 2), (0, 3)]:
        G[u, v] = G[v, u] = 1
    G[3, 2] = 1
    G[2, 1] = 1
    H = meek_closure(G)
    check("R4 orients 0->1", H[0, 1] == 1 and H[1, 0] == 0)

    # v-structure: 0 -> 2 <- 1, 0,1 non-adjacent : CPDAG keeps both directed
    D = np.zeros((3, 3), dtype=np.int8)
    D[0, 2] = 1
    D[1, 2] = 1
    C = dag_to_cpdag(D)
    check("collider CPDAG keeps 0->2,1->2", np.array_equal(C, D))

    # chain 0->1->2 : CPDAG is fully undirected
    D = np.zeros((3, 3), dtype=np.int8)
    D[0, 1] = 1
    D[1, 2] = 1
    C = dag_to_cpdag(D)
    check("chain CPDAG fully undirected",
          C[0, 1] == 1 and C[1, 0] == 1 and C[1, 2] == 1 and C[2, 1] == 1)


# --------------------------------------------------------------- CPDAG vs brute force
def test_cpdag_bruteforce(n_graphs=200, seed=0):
    print(f"\n[T1] CPDAG == common orientations over the MEC (brute force), {n_graphs} graphs")
    rng = np.random.default_rng(seed)
    bad = 0
    for _ in range(n_graphs):
        p = int(rng.integers(3, 7))
        D = random_dag(p, rng.uniform(1.0, 3.0), rng)
        C = dag_to_cpdag(D)
        # brute force: all DAGs with skeleton(D) and v_structures(D)
        base = skeleton(D)
        mec = consistent_dag_extensions(base, ref_vstructs=v_structures(D))
        if not mec:
            bad += 1
            continue
        Cbf = common_orientation_graph(mec, skeleton(D))
        if not np.array_equal(C, Cbf):
            bad += 1
    check("CPDAG matches brute force on all graphs", bad == 0, f"{bad} mismatches")


# --------------------------------------------------------------- Meek + K vs brute force
def test_mpdag_bruteforce(n_graphs=200, seed=1):
    print(f"\n[T2] Meek+K consistency AND MPDAG == brute force, {n_graphs} graphs x all K")
    rng = np.random.default_rng(seed)
    bad_cons = bad_mpdag = n_K = 0
    for _ in range(n_graphs):
        p = int(rng.integers(3, 7))
        D = random_dag(p, rng.uniform(1.0, 3.0), rng)
        C = dag_to_cpdag(D)
        U = undirected_edges(C)
        if not U:
            continue
        mec = consistent_dag_extensions(C, ref_vstructs=v_structures(C))
        m = min(3, len(U))
        for edges in combinations(U, m):
            for bits in range(2 ** m):
                K = [(e[1], e[0]) if (bits >> k) & 1 else e for k, e in enumerate(edges)]
                n_K += 1
                # brute force: DAGs in the MEC that agree with every statement in K
                sub = [Dx for Dx in mec if all(Dx[i, j] == 1 and Dx[j, i] == 0 for i, j in K)]
                cons_bf = len(sub) > 0
                cons_me = is_consistent(C, K)
                if cons_bf != cons_me:
                    bad_cons += 1
                    continue
                if cons_bf:
                    G = apply_background_knowledge(C, K)
                    Gbf = common_orientation_graph(sub, skeleton(C))
                    if not np.array_equal(G, Gbf):
                        bad_mpdag += 1
    check(f"Meek consistency verdict matches brute force ({n_K} K's)", bad_cons == 0,
          f"{bad_cons} mismatches")
    check("Meek closure == common orientations over consistent extensions", bad_mpdag == 0,
          f"{bad_mpdag} mismatches")


# --------------------------------------------------------------- O* validity
def test_ostar_valid(n_graphs=300, seed=2):
    print(f"\n[T3] O*(X,Y,G_Ktrue) is a valid adjustment set in the true DAG, {n_graphs} graphs")
    rng = np.random.default_rng(seed)
    n_tested = bad = 0
    for _ in range(n_graphs):
        p = int(rng.integers(4, 9))
        D = random_dag(p, rng.uniform(1.5, 3.0), rng)
        C = dag_to_cpdag(D)
        pairs = [(x, y) for x in range(p) for y in range(p)
                 if x != y and len(possibly_causal_paths(D, x, y)) > 0]
        if not pairs:
            continue
        x, y = pairs[int(rng.integers(len(pairs)))]
        # K = true orientations of a random subset of undirected edges
        U = undirected_edges(C)
        K = [(i, j) if D[i, j] == 1 else (j, i)
             for (i, j) in U if rng.random() < 0.5]
        G = apply_background_knowledge(C, K)   # true K is always consistent
        O = optimal_adjustment_set(G, x, y)
        if O is None:
            continue
        n_tested += 1
        if not is_valid_adjustment_set(D, x, y, O):
            bad += 1
    check(f"O* valid in true DAG whenever amenable ({n_tested} cases)", bad == 0, f"{bad} invalid")


# --------------------------------------------------------------- analytic vs Monte-Carlo
def test_analytic_vs_mc(n_scm=8, seed=3):
    print(f"\n[T4] Monte-Carlo check of the analytic linear-iSCM quantities, {n_scm} SCMs")
    rng = np.random.default_rng(seed)
    err_sigma, err_beta, err_tau = [], [], []
    for _ in range(n_scm):
        p = int(rng.integers(4, 7))
        D = random_dag(p, 2.0, rng)
        scm = make_linear_iscm(D, rng)
        Sigma = cov_linear(scm["A"], scm["omega"])
        X = sample_linear(scm, 200_000, rng)
        err_sigma.append(np.abs(np.cov(X, rowvar=False) - Sigma).max())
        # a random regression coefficient, analytic vs empirical
        pairs = [(x, y) for x in range(p) for y in range(p)
                 if x != y and len(possibly_causal_paths(D, x, y)) > 0]
        if not pairs:
            continue
        x, y = pairs[int(rng.integers(len(pairs)))]
        Z = [k for k in range(p) if k not in (x, y) and rng.random() < 0.5]
        b_an = ols_coefficient(Sigma, x, y, Z)
        S = [x] + sorted(Z)
        b_mc = np.linalg.lstsq(np.c_[X[:, S], np.ones(len(X))], X[:, y], rcond=None)[0][0]
        err_beta.append(abs(b_an - b_mc))
        # true total effect: analytic path formula vs interventional MC
        tau_an = total_effect_linear(scm["A"], x, y)
        y1 = sample_linear_do(scm, 200_000, rng, {x: 1.0})[:, y].mean()
        y0 = sample_linear_do(scm, 200_000, rng, {x: -1.0})[:, y].mean()
        err_tau.append(abs(tau_an - (y1 - y0) / 2.0))
    check(f"Sigma analytic == empirical (max err {max(err_sigma):.4f})", max(err_sigma) < 0.03)
    check(f"OLS coef analytic == empirical (max err {max(err_beta):.4f})", max(err_beta) < 0.03)
    check(f"total effect analytic == do-MC (max err {max(err_tau):.4f})", max(err_tau) < 0.03)


def test_valid_set_recovers_effect(n_scm=30, seed=4):
    print(f"\n[T5] valid adjustment set => population OLS coef == true total effect, {n_scm} SCMs")
    rng = np.random.default_rng(seed)
    worst = 0.0
    n = 0
    for _ in range(n_scm):
        p = int(rng.integers(4, 9))
        D = random_dag(p, 2.0, rng)
        scm = make_linear_iscm(D, rng)
        Sigma = cov_linear(scm["A"], scm["omega"])
        pairs = [(x, y) for x in range(p) for y in range(p)
                 if x != y and len(possibly_causal_paths(D, x, y)) > 0]
        if not pairs:
            continue
        for (x, y) in pairs:
            tau = total_effect_linear(scm["A"], x, y)
            for _ in range(5):
                Z = frozenset(k for k in range(p) if k not in (x, y) and rng.random() < 0.5)
                if not is_valid_adjustment_set(D, x, y, Z):
                    continue
                b = ols_coefficient(Sigma, x, y, Z)
                worst = max(worst, abs(b - tau))
                n += 1
    check(f"every valid Z reproduces tau exactly ({n} sets, max err {worst:.2e})", worst < 1e-8)


if __name__ == "__main__":
    test_meek_textbook()
    test_cpdag_bruteforce()
    test_mpdag_bruteforce()
    test_ostar_valid()
    test_analytic_vs_mc()
    test_valid_set_recovers_effect()
    print("\n" + "=" * 60)
    if FAILS:
        print(f"FAILED: {FAILS}")
        sys.exit(1)
    print("ALL TESTS PASSED")
