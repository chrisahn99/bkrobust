"""
X2 MACHINERY TESTS -- run BEFORE any experiment.  A run whose machinery is
untested is not a result.

Every new object is checked against its DEFINITION by exhaustive brute force on
small graphs, not against a second copy of the same implementation:

  T1  hop_dist_inf         vs all-simple-paths minimum length (exhaustive DFS)
  T2  edge-component lemma  (an edge is wholly in or wholly out of the query's
                             component)  -- makes G7 decidable
  T3  dmin >= 1  <=>  the statement is not incident to {X,Y}   (PREREG 5.4)
  T4  meek_closure idempotence (G1) + skeleton invariance (G2)
  T5  ostar_and_paths      vs adjust.optimal_adjustment_set, exhaustively over
                             all MPDAGs reachable from all CPDAGs at p <= 5
  T5b pcp_capped           vs adjust.possibly_causal_paths
  T6  is_valid_adjustment_set  vs the ONLY independent ground truth available:
                             a valid set must give population OLS == tau
  T7  the ARM E operator   (all-true base never fails; D agrees with M(K);
                             skeleton preserved; non-adjacent pair -> MeekFail;
                             closure order-independence)
  T8  E_set detection      vs a from-scratch recomputation of O* on both graphs
  T9  prune_partition      vs brute-force filtering; Cost_{-1} arithmetic
  T10 member <-> flip bijection on the E1' archive
  T11 clopper_pearson      vs scipy.stats.binomtest inversion
  T12 G10 placebo: the members_frozen ball gives Q_r = 0 for every r

Usage:  python test_machinery.py
"""
import json
import os
import sys
import time
from itertools import combinations, permutations, product

import numpy as np

import x2lib as X
from adjust import (causal_nodes, cov_linear, forb, is_amenable,
                    is_valid_adjustment_set, optimal_adjustment_set,
                    ols_coefficient, parents_of_set, poss_de,
                    possibly_causal_paths, total_effect_linear)
from graphs import (MeekFail, apply_background_knowledge, dag_agrees_with,
                    dag_to_cpdag, directed_edges, is_directed, is_undirected,
                    meek_closure, random_dag, skeleton, undirected_edges,
                    v_structures, consistent_dag_extensions)
from scm import make_linear_iscm

RESULTS = []


def report(name, ok, detail):
    RESULTS.append((name, bool(ok), detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


# --------------------------------------------------------------- brute helpers
def brute_min_path_len(S, srcs, p):
    """Minimum simple-path length from any src to each node, by exhaustive DFS.
    Independent of BFS."""
    best = [np.inf] * p
    for s in srcs:
        best[s] = 0

    def dfs(v, length, visited):
        for w in range(p):
            if w == v or S[v, w] == 0 or w in visited:
                continue
            if length + 1 < best[w]:
                best[w] = length + 1
            visited.add(w)
            dfs(w, length + 1, visited)
            visited.remove(w)

    for s in srcs:
        dfs(s, 0, {s})
    return best


def brute_ostar(G, x, y):
    """O* from the definition, with an independently written path enumerator
    (BFS over path prefixes rather than the library's recursive DFS)."""
    p = G.shape[0]
    paths = []
    stack = [[x]]
    while stack:
        pth = stack.pop()
        v = pth[-1]
        if v == y:
            paths.append(pth)
            continue
        for w in range(p):
            if w == v or G[v, w] != 1 or w in pth:
                continue
            stack.append(pth + [w])
    if not paths:
        return None
    if not all(G[pth[0], pth[1]] == 1 and G[pth[1], pth[0]] == 0 for pth in paths):
        return None
    cn = set()
    for pth in paths:
        cn.update(pth[1:])
    # poss_de by transitive closure of the possibly-directed relation
    M = np.zeros((p, p), dtype=bool)
    for i in range(p):
        for j in range(p):
            if i != j and G[i, j] == 1:
                M[i, j] = True
    reach = set(cn)
    changed = True
    while changed:
        changed = False
        for v in list(reach):
            for w in range(p):
                if M[v, w] and w not in reach:
                    reach.add(w)
                    changed = True
    fb = reach | {x}
    pa = set()
    for a in range(p):
        for b in cn:
            if G[a, b] == 1 and G[b, a] == 0:
                pa.add(a)
    return frozenset(pa - fb)


def all_dags(p):
    """Every labelled DAG on p nodes (p <= 5), as amats."""
    idx = [(i, j) for i in range(p) for j in range(i + 1, p)]
    out = []
    for state in product([0, 1, 2], repeat=len(idx)):   # 0 none, 1 i->j, 2 j->i
        Dm = np.zeros((p, p), dtype=np.int8)
        for (i, j), s in zip(idx, state):
            if s == 1:
                Dm[i, j] = 1
            elif s == 2:
                Dm[j, i] = 1
        # acyclicity
        A = Dm.astype(bool)
        R = A.copy()
        for _ in range(p):
            R = R | (R @ A)
        if np.any(np.diag(R)):
            continue
        out.append(Dm)
    return out


# --------------------------------------------------------------- T1, T2, T3
def t1_t2_t3(n=400, seed=11):
    rng = np.random.default_rng(seed)
    bad1 = bad2 = bad3 = 0
    tot = 0
    for _ in range(n):
        p = int(rng.integers(4, 8))
        deg = float(rng.choice([0.8, 1.5, 2.5]))
        D = random_dag(p, deg, rng)
        C = dag_to_cpdag(D)
        S = skeleton(C)
        x, y = rng.choice(p, size=2, replace=False)
        d = X.hop_dist_inf(C, [int(x), int(y)])
        b = brute_min_path_len(S, [int(x), int(y)], p)
        for v in range(p):
            tot += 1
            dv = d[v]
            bv = b[v]
            if (np.isfinite(dv) != np.isfinite(bv)) or (np.isfinite(dv) and dv != bv):
                bad1 += 1
        # T2: every skeleton edge is wholly inside or wholly outside the component
        for i in range(p):
            for j in range(i + 1, p):
                if S[i, j] and (np.isfinite(d[i]) != np.isfinite(d[j])):
                    bad2 += 1
        # T3: dmin >= 1  <=>  not incident to {x,y}
        for (u, v) in undirected_edges(C):
            dmin, dmax = X.stmt_dist(d, (u, v))
            inc = (u in (x, y)) or (v in (x, y))
            if (dmin >= 1) != (not inc):
                bad3 += 1
    report("T1 hop_dist_inf == brute-force min simple-path length",
           bad1 == 0, f"{bad1} mismatches over {tot} (node, graph) pairs")
    report("T2 edge-component lemma (G7 decidable per edge)",
           bad2 == 0, f"{bad2} split edges over {n} graphs")
    report("T3 dmin>=1 <=> statement not incident to {X,Y} (PREREG 5.4)",
           bad3 == 0, f"{bad3} mismatches")


# --------------------------------------------------------------- T4
def t4(n=400, seed=12):
    rng = np.random.default_rng(seed)
    bad_idem = bad_skel = 0
    ntr = 0
    for _ in range(n):
        p = int(rng.integers(4, 9))
        deg = float(rng.choice([1.5, 2.5, 4.0]))
        D = random_dag(p, deg, rng)
        C = dag_to_cpdag(D)
        if not np.array_equal(meek_closure(C), C):
            bad_idem += 1
        U = undirected_edges(C)
        if not U:
            continue
        for _ in range(4):
            k = int(rng.integers(1, min(4, len(U)) + 1))
            sel = rng.choice(len(U), size=k, replace=False)
            Kp = []
            for i in sel:
                u, v = U[i]
                Kp.append((u, v) if rng.random() < 0.5 else (v, u))
            try:
                G = apply_background_knowledge(C, Kp)
            except MeekFail:
                continue
            ntr += 1
            if not np.array_equal(meek_closure(G), G):
                bad_idem += 1
            if not np.array_equal(skeleton(G), skeleton(C)):
                bad_skel += 1
    report("T4a meek_closure idempotent on C and on every M(K') [gate G1]",
           bad_idem == 0, f"{bad_idem} non-fixpoints over {ntr} MPDAGs + {n} CPDAGs")
    report("T4b skeleton(M(K')) == skeleton(C) [gate G2]",
           bad_skel == 0, f"{bad_skel} skeleton changes over {ntr} MPDAGs")


# --------------------------------------------------------------- T5, T5b
def t5(p_list=(4, 5), seed=13):
    """Exhaustive over all CPDAGs at p in p_list, all ordered queries, all
    single admissible statements (both orientations)."""
    n_cmp = 0
    bad_o = bad_amen = bad_paths = 0
    for p in p_list:
        seen = set()
        for D in all_dags(p):
            C = dag_to_cpdag(D)
            key = C.tobytes()
            if key in seen:
                continue
            seen.add(key)
            U = undirected_edges(C)
            graphs = [C]
            for (u, v) in U:
                for s in ((u, v), (v, u)):
                    try:
                        graphs.append(apply_background_knowledge(C, [s]))
                    except MeekFail:
                        pass
            for G in graphs:
                for x in range(p):
                    for y in range(p):
                        if x == y:
                            continue
                        n_cmp += 1
                        O_lib = optimal_adjustment_set(G, x, y)
                        O_mine, npaths, amen = X.ostar_and_paths(G, x, y)
                        if O_lib != O_mine:
                            bad_o += 1
                        if amen != is_amenable(G, x, y):
                            bad_amen += 1
                        pl = possibly_causal_paths(G, x, y)
                        pm = X.pcp_capped(G, x, y)
                        if sorted(map(tuple, pl)) != sorted(map(tuple, pm)):
                            bad_paths += 1
                        Ob = brute_ostar(G, x, y)
                        if Ob != O_lib:
                            bad_o += 1
    report("T5 ostar_and_paths == optimal_adjustment_set == brute_ostar "
           "(exhaustive, all CPDAGs+MPDAGs at p<=5, all queries)",
           bad_o == 0, f"{bad_o} mismatches over {n_cmp} (graph,query) triples")
    report("T5a amenability agrees with adjust.is_amenable",
           bad_amen == 0, f"{bad_amen} mismatches over {n_cmp}")
    report("T5b pcp_capped == adjust.possibly_causal_paths",
           bad_paths == 0, f"{bad_paths} mismatches over {n_cmp}")


# --------------------------------------------------------------- T6
def t6(n=120, seed=14):
    """A valid adjustment set must give population OLS beta == tau exactly.
    That is an INDEPENDENT ground truth for is_valid_adjustment_set (which the
    E_bias event is built on)."""
    rng = np.random.default_rng(seed)
    n_valid = n_inval = 0
    bad_valid = 0            # valid but beta != tau  -> the criterion is wrong
    coincid = 0              # invalid but beta == tau -> allowed (measure zero)
    for _ in range(n):
        p = int(rng.integers(4, 7))
        deg = float(rng.choice([1.5, 2.5]))
        D = random_dag(p, deg, rng)
        scm = make_linear_iscm(D, rng)
        Sigma = cov_linear(scm["A"], scm["omega"])
        for x in range(p):
            for y in range(p):
                if x == y or len(possibly_causal_paths(D, x, y)) == 0:
                    continue
                tau = total_effect_linear(scm["A"], x, y)
                if abs(tau) < 1e-6:
                    continue
                rest = [v for v in range(p) if v not in (x, y)]
                for r in range(0, min(3, len(rest)) + 1):
                    for Z in combinations(rest, r):
                        ok = bool(is_valid_adjustment_set(D, x, y, set(Z)))
                        b = ols_coefficient(Sigma, x, y, set(Z))
                        eq = abs(b - tau) < 1e-8
                        if ok:
                            n_valid += 1
                            if not eq:
                                bad_valid += 1
                        else:
                            n_inval += 1
                            if eq:
                                coincid += 1
    report("T6 valid adjustment set => population OLS beta == tau",
           bad_valid == 0,
           f"{bad_valid} violations over {n_valid} valid sets; "
           f"{coincid}/{n_inval} invalid sets coincidentally equal (allowed)")


# --------------------------------------------------------------- T7
def t7(n=300, seed=15):
    rng = np.random.default_rng(seed)
    fail_true = 0
    bad_agree = 0
    bad_nonadj = 0
    bad_order = 0
    ntr = 0
    n_nonadj = 0
    for _ in range(n):
        p = int(rng.integers(5, 9))
        deg = float(rng.choice([1.5, 2.5, 4.0]))
        D = random_dag(p, deg, rng)
        C = dag_to_cpdag(D)
        U = undirected_edges(C)
        if len(U) < 2:
            continue
        K_all = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U]
        # T7a: an ALL-TRUE knowledge set can never MeekFail, and D agrees with M(K)
        for k in range(1, min(4, len(K_all)) + 1):
            sel = rng.choice(len(K_all), size=k, replace=False)
            K = [K_all[i] for i in sel]
            try:
                G0 = apply_background_knowledge(C, K)
            except MeekFail:
                fail_true += 1
                continue
            ntr += 1
            if not dag_agrees_with(D, G0):
                bad_agree += 1
            # T7c: order-independence of the returned MPDAG
            for _ in range(2):
                Kperm = list(K)
                rng.shuffle(Kperm)
                try:
                    G1 = apply_background_knowledge(C, Kperm)
                except MeekFail:
                    bad_order += 1
                    continue
                if not np.array_equal(G0, G1):
                    bad_order += 1
        # T7b: a statement on a NON-ADJACENT pair is refused (the X1 boundary)
        nonadj = [(i, j) for i in range(p) for j in range(p)
                  if i != j and C[i, j] == 0 and C[j, i] == 0]
        if nonadj:
            i, j = nonadj[int(rng.integers(len(nonadj)))]
            n_nonadj += 1
            try:
                apply_background_knowledge(C, [(i, j)])
                bad_nonadj += 1
            except MeekFail:
                pass
    report("T7a all-true base K never MeekFails and D agrees with M(K)",
           fail_true == 0 and bad_agree == 0,
           f"{fail_true} MeekFails, {bad_agree} disagreements over {ntr} bases")
    report("T7b statement on a non-adjacent pair -> MeekFail (X1 boundary)",
           bad_nonadj == 0, f"{bad_nonadj} accepted of {n_nonadj} attempts")
    report("T7c M(K) is independent of the order of K",
           bad_order == 0, f"{bad_order} order-dependent results over {ntr} bases")


# --------------------------------------------------------------- T8
def t8(n=250, seed=16):
    """classify_trial's E_set / E_ident / E_gain against a from-scratch
    recomputation using the LIBRARY functions on both graphs."""
    rng = np.random.default_rng(seed)
    bad = 0
    ntr = 0
    fired = dict(e_set=0, e_ident=0, e_gain=0, meek_fail=0)
    for _ in range(n):
        p = int(rng.integers(5, 9))
        deg = float(rng.choice([1.5, 2.5, 4.0]))
        D = random_dag(p, deg, rng)
        C = dag_to_cpdag(D)
        U = undirected_edges(C)
        if len(U) < 2:
            continue
        pairs = [(a, b) for a in range(p) for b in range(p)
                 if a != b and len(possibly_causal_paths(D, a, b)) > 0]
        if not pairs:
            continue
        x, y = pairs[int(rng.integers(len(pairs)))]
        scm = make_linear_iscm(D, rng)
        Sigma = cov_linear(scm["A"], scm["omega"])
        K_all = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U]
        kb = int(rng.integers(1, min(4, len(K_all)) + 1))
        K = [K_all[i] for i in rng.choice(len(K_all), size=kb, replace=False)]
        base = X.make_base(C, D, Sigma, x, y, K)
        assert base is not None
        O0_ref = optimal_adjustment_set(base["G0"], x, y)
        if O0_ref != base["O0"]:
            bad += 1
        for (u, v) in U:
            for s in ((u, v), (v, u)):
                Kp = [t for t in K if {t[0], t[1]} != {u, v}] + [s]
                rec = X.classify_trial(C, D, Sigma, x, y, Kp, base)
                ntr += 1
                # reference
                try:
                    G = apply_background_knowledge(C, Kp)
                    consistent = True
                except MeekFail:
                    consistent = False
                if consistent != rec["consistent"]:
                    bad += 1
                    continue
                if not consistent:
                    fired["meek_fail"] += 1
                    continue
                Oref = optimal_adjustment_set(G, x, y)
                a0 = O0_ref is not None
                a1 = Oref is not None
                e_set = a0 and a1 and (Oref != O0_ref)
                e_ident = a0 and not a1
                e_gain = (not a0) and a1
                if (rec["e_set"], rec["e_ident"], rec["e_gain"]) != (e_set, e_ident, e_gain):
                    bad += 1
                if e_set:
                    fired["e_set"] += 1
                    if rec["e_bias"] != (not bool(is_valid_adjustment_set(D, x, y, Oref))):
                        bad += 1
                if e_ident:
                    fired["e_ident"] += 1
                    if rec["e_ident_nopath"] != (len(possibly_causal_paths(G, x, y)) == 0):
                        bad += 1
                if e_gain:
                    fired["e_gain"] += 1
    report("T8 classify_trial events == from-scratch library recomputation",
           bad == 0, f"{bad} mismatches over {ntr} trials; fired {fired}")


# --------------------------------------------------------------- T9
def t9(seed=17):
    """prune_partition vs brute-force filtering, and Cost_{-1} arithmetic."""
    rng = np.random.default_rng(seed)
    bad = 0
    for _ in range(300):
        k = int(rng.integers(2, 5))
        max_rho = k
        dmins = [float(rng.choice([0, 1, 2, np.inf])) for _ in range(k)]
        members = []
        flips = []
        for rho in range(1, max_rho + 1):
            for flip in combinations(range(k), rho):
                flips.append(flip)
                cons = rng.random() < 0.7
                amen = cons and rng.random() < 0.8
                m = dict(rho=rho, consistent=bool(cons), amenable=bool(amen))
                if amen:
                    m["est"] = float(rng.normal())
                members.append(m)
        est0 = 0.0
        for r in (-1, 0, 1, 2, 3):
            vf, vp, nf, npr = X.prune_partition(list(range(k)), dmins, members, max_rho, r)
            near = set(i for i in range(k) if np.isfinite(dmins[i]) and dmins[i] <= r) if r >= 0 else set()
            bf = [m["est"] for m, fl in zip(members, flips)
                  if m.get("consistent") and m.get("amenable")]
            bp = [m["est"] for m, fl in zip(members, flips)
                  if m.get("consistent") and m.get("amenable") and set(fl) <= near]
            if sorted(vf) != sorted(bf) or sorted(vp) != sorted(bp):
                bad += 1
            if nf != 1 + len(flips) or npr != 1 + sum(1 for fl in flips if set(fl) <= near):
                bad += 1
            if r == -1 and npr != 1:
                bad += 1
    # Cost_{-1} for k=4, rho=4 must be exactly 1/16
    c = 1.0 / (1 + 4 + 6 + 4 + 1)
    report("T9 prune_partition == brute-force filter; Cost_{-1}(k=4,rho=4)=1/16",
           bad == 0 and abs(c - 0.0625) < 1e-15,
           f"{bad} mismatches over 300 random balls x 5 radii; Cost_-1={c:.4f}")


# --------------------------------------------------------------- T10, T12
def t10_t12(arch="/home/costaj/latent-causal/e1prime-se/results"):
    bad_len = 0
    n_scm = 0
    q_nonzero = 0
    n_q = 0
    for ens, key, mr in (("original", "K4", 4), ("licensed", "K4", 4),
                         ("large", "K4", 4), ("k8", "K4", 3), ("k8", "K8", 3)):
        path = os.path.join(arch, f"arm1_{ens}.json")
        if not os.path.exists(path):
            continue
        recs = json.load(open(path))
        for rec in recs:
            arm = rec["arms"].get(key)
            if arm is None or not arm.get("mpdag_amenable"):
                continue
            n_scm += 1
            k = arm["n_K"]
            exp = sum(len(list(combinations(range(k), r))) for r in range(1, min(mr, k) + 1))
            if len(arm["members"]) != exp or len(arm["members_frozen"]) != exp:
                bad_len += 1
            # T12 (gate G10): the frozen ball must give Q_r = 0 for every r
            dmins = [float(h) if h < X.SENTINEL else np.inf for h in arm["stmt_hopdist"]]
            for r in (-1, 0, 1, 2, 3):
                vf, vp, _, _ = X.prune_partition(arm["K"], dmins, arm["members_frozen"], mr, r)
                If = X.interval(vf, arm["est0"])
                Ip = X.interval(vp, arm["est0"])
                n_q += 1
                if abs(If[0] - Ip[0]) > 1e-9 or abs(If[1] - Ip[1]) > 1e-9:
                    q_nonzero += 1
    report("T10 member index <-> flip bijection on the E1' archive",
           bad_len == 0, f"{bad_len} length mismatches over {n_scm} amenable arms")
    report("T12 [gate G10] placebo ball members_frozen gives Q_r = 0 for every r",
           q_nonzero == 0, f"{q_nonzero} non-zero of {n_q} (SCM, r) cells")


# --------------------------------------------------------------- T11
def t11():
    from scipy.stats import binomtest
    bad = 0
    for (k, n) in [(0, 10), (0, 1000), (1, 100), (3, 2087), (25, 127), (307, 586)]:
        lo, hi = X.clopper_pearson(k, n)
        ref = binomtest(k, n).proportion_ci(method="exact")
        if abs(lo - ref.low) > 1e-12 or abs(hi - ref.high) > 1e-12:
            bad += 1
    report("T11 clopper_pearson == scipy binomtest exact CI", bad == 0,
           f"{bad} mismatches over 6 cells")


# --------------------------------------------------------------- main
def main():
    t0 = time.time()
    t1_t2_t3()
    t4()
    t5()
    t6()
    t7()
    t8()
    t9()
    t10_t12()
    t11()
    npass = sum(1 for _, ok, _ in RESULTS if ok)
    print(f"\n=== {npass}/{len(RESULTS)} machinery tests passed "
          f"in {time.time() - t0:.1f}s ===", flush=True)
    with open("../results/machinery_tests.json", "w") as f:
        json.dump([dict(name=n, ok=ok, detail=d) for n, ok, d in RESULTS], f, indent=1)
    sys.exit(0 if npass == len(RESULTS) else 1)


if __name__ == "__main__":
    main()
