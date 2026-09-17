"""
E2 -- the empirical half of the R1-closed finiteness claim.

Three arms on the SAME 600 SCMs (the pilot's primary ensemble, regenerated
bitwise from the stored (seed, p, deg) triples):

  G-flip : generic K (true), reverse rho in {1,2,3} statements
  T-flip : tiered  K (true, FULL cross-tier set), reverse rho in {1,2,3}   <- destroys the tiering
  T-move : tiered  K (true, FULL cross-tier set), move rho in {1,2} nodes  <- KEEPS K' a tiering

Two pilot defects fixed (both load-bearing, see PREREGISTRATION.md sec.1):
  * run_linear.py:55  the `if D[u,v]==1` filter that DROPS statements disagreeing
    with the true DAG is REMOVED -- under a perturbed tiering those statements
    are the signal.
  * run_linear.py:34/56-58  MAX_K=4 subsampling of tiered K is REMOVED -- Bang &
    Didelez's R1-only property is about the FULL cross-tier constraint set.
    The generic arm is matched to |K_T| per SCM so the budget is equal.

Usage: python run_e2.py [out.json]
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
                    dag_to_cpdag, random_dag, topological_order,
                    undirected_edges)
from meek_rules import (ALL_RULES, apply_background_knowledge_rules,
                        cascade_size, n_orientation_changes,
                        sym_diff_statements)
from scm import make_linear_iscm

CAP_FLIP = 400     # max sampled flip-sets per rho (>= C(|K|,3) here: no subsampling occurs)
CAP_MOVE2 = 400    # max sampled 2-node tier moves (>= p*(T-1) choose 2 here: exhaustive)
KNOW_RNG_OFFSET = 99_000_000


# --------------------------------------------------------------- knowledge builders
def base_tiering(D, C, rng):
    """Bang & Didelez tiering: contiguous cuts of a topological order of D*.
    Returns (tier vector, K) with K the FULL cross-tier orientation set of C's
    undirected edges -- no MAX_K cap, no agreement filter."""
    order = topological_order(D)
    p = len(order)
    n_tiers = int(rng.integers(2, min(4, p) + 1))
    cuts = sorted(rng.choice(np.arange(1, p), size=n_tiers - 1, replace=False))
    tier = np.zeros(p, dtype=int)
    t = 0
    for pos, node in enumerate(order):
        while t < len(cuts) and pos >= cuts[t]:
            t += 1
        tier[node] = t
    return tier, n_tiers, K_from_tiering(tier, C)


def K_from_tiering(tier, C):
    """K = every undirected edge of C whose endpoints sit in different tiers,
    oriented lower-tier -> higher-tier."""
    return [((i, j) if tier[i] < tier[j] else (j, i))
            for (i, j) in undirected_edges(C) if tier[i] != tier[j]]


def generic_knowledge_matched(D, C, rng, size):
    """|K_G| = |K_T|, drawn from the undirected edges of C, oriented per D*."""
    U = undirected_edges(C)
    K_all = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U]
    if len(K_all) > size:
        idx = rng.choice(len(K_all), size=size, replace=False)
        K_all = [K_all[i] for i in sorted(idx)]
    return K_all


# --------------------------------------------------------------- ball iterators
def flip_ball(K, rng, cap=CAP_FLIP):
    """(rho, K') for every flip-set of size rho in {1,2,3} (sampled if > cap)."""
    out = []
    for rho in (1, 2, 3):
        if rho > len(K):
            break
        allf = list(combinations(range(len(K)), rho))
        if len(allf) > cap:
            idx = rng.choice(len(allf), size=cap, replace=False)
            allf = [allf[i] for i in sorted(idx)]
        for flip in allf:
            out.append((rho, [(v, u) if k in flip else (u, v)
                              for k, (u, v) in enumerate(K)]))
    return out


def move_ball(tier, n_tiers, C, rng, cap2=CAP_MOVE2):
    """(rho_move, K') for every reassignment of rho in {1,2} nodes to a DIFFERENT
    tier. K' is derived from the perturbed tiering, so it is STILL A TIERING --
    this is the whole point of the arm."""
    p = len(tier)
    out = []
    singles = [(v, t) for v in range(p) for t in range(n_tiers) if t != tier[v]]
    for (v, t) in singles:
        tp = tier.copy()
        tp[v] = t
        out.append((1, K_from_tiering(tp, C), [(int(v), int(t))]))
    doubles = [((v1, t1), (v2, t2))
               for i, (v1, t1) in enumerate(singles)
               for (v2, t2) in singles[i + 1:] if v1 != v2]
    if len(doubles) > cap2:
        idx = rng.choice(len(doubles), size=cap2, replace=False)
        doubles = [doubles[i] for i in sorted(idx)]
    for ((v1, t1), (v2, t2)) in doubles:
        tp = tier.copy()
        tp[v1] = t1
        tp[v2] = t2
        out.append((2, K_from_tiering(tp, C), [(int(v1), int(t1)), (int(v2), int(t2))]))
    return out


# --------------------------------------------------------------- member evaluation
def eval_member(C, U, D, Sigma, x, y, tau, K, Kp, G0, O0, rho, moves=None):
    rec = dict(rho=rho, n_K=len(K), n_Kp=len(Kp),
               n_changed=sym_diff_statements(K, Kp))
    if moves is not None:
        rec["moves"] = moves
    try:
        G = apply_background_knowledge(C, Kp)
    except MeekFail:
        rec.update(consistent=False)
        return rec
    # R1-ONLY closure of the SAME knowledge set -- attribution-free comparison
    try:
        G1 = apply_background_knowledge_rules(C, Kp, ("R1",))
        rec["r1_only_fail"] = False
        rec["r1_equals_full"] = bool(np.array_equal(G1, G))
        rec["n_oriented_r1"] = int(sum(1 for (i, j) in U
                                       if not (G1[i, j] == 1 and G1[j, i] == 1)))
    except MeekFail:
        rec["r1_only_fail"] = True
        rec["r1_equals_full"] = False
        rec["n_oriented_r1"] = None

    O = optimal_adjustment_set(G, x, y)
    rec.update(
        consistent=True,
        expels_true_dag=not dag_agrees_with(D, G),
        cascade=int(cascade_size(G, U, Kp)),
        d_orient=int(n_orientation_changes(G0, G, U)),
        amenable=O is not None,
    )
    if O is not None:
        rec["ostar_changed"] = bool(O != O0)
        rec["ostar_valid"] = bool(is_valid_adjustment_set(D, x, y, O))
        rec["ostar"] = sorted(int(v) for v in O)
        est = ols_coefficient(Sigma, x, y, O)
        rec["est"] = est
        rec["bias"] = est - tau
    return rec


# --------------------------------------------------------------- per-SCM driver
def analyse_one(args):
    seed, p, deg = args
    # ---- regenerate the pilot's SCM, consuming rng in EXACTLY analyse_one's order
    rng = np.random.default_rng(seed)
    D = random_dag(p, deg, rng)
    C = dag_to_cpdag(D)
    U = undirected_edges(C)
    if len(U) < 3:
        return None
    pairs = [(a, b) for a in range(p) for b in range(p)
             if a != b and len(possibly_causal_paths(D, a, b)) > 0]
    if not pairs:
        return None
    x, y = pairs[int(rng.integers(len(pairs)))]
    scm = make_linear_iscm(D, rng)
    Sigma = cov_linear(scm["A"], scm["omega"])
    tau = total_effect_linear(scm["A"], x, y)
    if abs(tau) < 1e-3:
        return None

    # ---- knowledge draws use a SEPARATE rng: this experiment's arms are not the
    #      pilot's draws, and must not be pooled with them (MAX_K sits upstream of
    #      a conditional rng.choice in the pilot, so the streams differ anyway).
    krng = np.random.default_rng(seed + KNOW_RNG_OFFSET)
    # A tiering with < 2 cross-tier orientable edges gives a degenerate ball, so
    # re-draw the ELICITATION (never an outcome) up to 20 times. Re-drawing and
    # draw-once-then-drop condition on the SAME event, so the conditional
    # distribution of tierings is identical; re-drawing merely retains more SCMs.
    for n_attempts in range(1, 21):
        tier, n_tiers, K_T = base_tiering(D, C, krng)
        if len(K_T) >= 2:
            break
    else:
        return None
    K_G = generic_knowledge_matched(D, C, krng, len(K_T))

    out = dict(seed=seed, p=p, deg=deg, x=x, y=y, tau=tau,
               n_undirected=len(U), n_edges=int(D.sum()), n_tiers=n_tiers,
               tier=[int(t) for t in tier],
               cpdag_amenable=bool(is_amenable(C, x, y)), n_attempts=n_attempts, arms={})

    for arm, K in (("G-flip", K_G), ("T-flip", K_T), ("T-move", K_T)):
        try:
            G0 = apply_background_knowledge(C, K)
        except MeekFail:
            continue                      # cannot happen for true K; be safe
        O0 = optimal_adjustment_set(G0, x, y)
        if O0 is None:
            continue                      # not amenable even under TRUE K -> excluded
        est0 = ols_coefficient(Sigma, x, y, O0)

        # sanity: TRUE tiered K must be R1-closed (Bang & Didelez, unperturbed)
        try:
            G0_r1 = apply_background_knowledge_rules(C, K, ("R1",))
            r1_base = bool(np.array_equal(G0_r1, G0))
        except MeekFail:
            r1_base = False

        members = []
        if arm == "T-move":
            for (rho, Kp, moves) in move_ball(tier, n_tiers, C,
                                              np.random.default_rng(seed + 7)):
                members.append(eval_member(C, U, D, Sigma, x, y, tau,
                                           K, Kp, G0, O0, rho, moves))
        else:
            for (rho, Kp) in flip_ball(K, np.random.default_rng(seed + 11)):
                members.append(eval_member(C, U, D, Sigma, x, y, tau,
                                           K, Kp, G0, O0, rho))

        out["arms"][arm] = dict(K=[list(e) for e in K], n_K=len(K),
                                O0=sorted(int(v) for v in O0), est0=est0,
                                r1_base_equals_full=r1_base, members=members)
    if not out["arms"]:
        return None
    return out


def main():
    out_path = sys.argv[1] if len(sys.argv) > 1 else "../results/e2_raw.json"
    ref_path = ("/Users/josecosta/research-pilots/"
                "latent-causal-rho-breakdown-knowledge/results/linear_raw.json")
    with open(ref_path) as f:
        ref = json.load(f)
    jobs = [(r["seed"], r["p"], r["deg"]) for r in ref]
    print(f"replaying {len(jobs)} pilot SCMs")
    with Pool(8) as pool:
        res = [r for r in pool.imap_unordered(analyse_one, jobs, chunksize=4) if r]

    # ---- MANDATORY: bitwise ensemble check against the pilot
    ref_tau = {r["seed"]: r["tau"] for r in ref}
    worst = 0.0
    for r in res:
        worst = max(worst, abs(r["tau"] - ref_tau[r["seed"]]))
    print(f"kept {len(res)} SCMs; worst |tau - tau_pilot| = {worst:.3e}")
    assert worst < 1e-12, "ENSEMBLE MISMATCH -- run is void"

    with open(out_path, "w") as f:
        json.dump(res, f)
    print(f"-> {out_path}")


if __name__ == "__main__":
    main()
