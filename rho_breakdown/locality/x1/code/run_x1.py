"""
X1 ARM SWEEP -- the spurious required-edge class, matched to the reversal class.

Balls per SCM, all on the SAME C, D, Sigma, X, Y, K (gate G1):
  R          reversal        -- apply_background_knowledge VERBATIM (E1' control)
  S-uni      spurious/uniform-- bk_assert, replacement drawn from N(C)      PRIMARY
  S-loc      spurious/b-LOAD -- bk_assert, replacement from b-LOAD's own pool
  D          DROP control    -- slot removed, NOTHING asserted in its place
  P-null-id  placebo         -- slot replaced by ITSELF, full bk_assert path (G5a)
  P-null-app placebo         -- K ++ [already-directed edge of C]           (G5b)
  P-frozen   placebo         -- O frozen at O0                              (G4)

ARM D IS NOT IN THE PREREG AND IS NOT IN THE AUDIT.  It is here because the
machinery tests found, before any outcome was computed, that the registered
replacement law confounds the class effect with a WITHDRAWN TRUE STATEMENT:
  K'(F) = (r_t if t in F else (u_t,v_t))  replaces slot t, so arm S asserts only
  3 of the original 4 true statements while arm R still asserts all 4.
Measured on 881 unreachable-pair assertions: the ADDED EDGE moved O* in 0/881,
the DROP alone moved it in 161/881 = 0.1827, and O*(S) == O*(D) in 881/881.  So
every raw S-vs-R contrast carries a withdrawal effect that has nothing to do with
the spurious class.  AUDIT A6 found the special case (P-null) and not the general
one.  Both contrasts are computed and both are reported.

Usage: python run_x1.py <n_analysed_target> <ensemble> <out.json> [n_workers]
"""
import json
import sys
import time
from itertools import combinations
from multiprocessing import Pool

import numpy as np

from adjust import (PathBlowup, cov_linear, is_amenable,
                    is_valid_adjustment_set, optimal_adjustment_set,
                    possibly_causal_paths, total_effect_linear)
from graphs import (MeekFail, apply_background_knowledge, dag_agrees_with,
                    dag_to_cpdag, directed_edges, random_dag, skeleton,
                    undirected_edges)
from scm import make_linear_iscm
from se import ols_with_se
from x1_ops import (UNREACH, bk_assert, bk_assert_stampall, bload_pool,
                    draw_sloc, draw_suni, hop_dist_from, label_statement,
                    nonadjacent_pairs, pdag_extendable, stmt_hop,
                    vstruct_decomposition)

MAX_K = 4

# grids imported UNCHANGED from run_arm1.py:35-42 so ARM R reproduces E1' (G2)
ENSEMBLES = {
    "original": dict(ps=range(5, 9),  degs=[1.5, 2.0, 2.5],           need_u=3, max_rho=4),
    "licensed": dict(ps=range(5, 11), degs=[1.5, 2.0, 2.5, 4.0, 6.0], need_u=3, max_rho=4),
    "large":    dict(ps=range(15, 26), degs=[1.5, 2.0, 2.5],          need_u=3, max_rho=4),
}
# AUDIT M10: targets are stated in ANALYSED SCMs; OVERSAMPLE raised accordingly.
OVERSAMPLE = {"original": 6, "licensed": 8, "large": 4}


# --------------------------------------------------------------------------- balls
def _outcome(G, D, Sigma, x, y, O0, est0):
    """The estimand side, shared by every ball.  Returns (dict, blowup)."""
    try:
        O = optimal_adjustment_set(G, x, y)
    except PathBlowup:
        return dict(amenable=False, blowup=True), True
    if O is None:
        return dict(amenable=False, blowup=False), False
    try:
        e, _, s2 = ols_with_se(Sigma, x, y, O, np.inf)
        Sm = [x] + sorted(O)
        Sminv = np.linalg.inv(Sigma[np.ix_(Sm, Sm)])
        ov = bool(is_valid_adjustment_set(D, x, y, O))
    except PathBlowup:
        return dict(amenable=False, blowup=True), True
    except np.linalg.LinAlgError:
        return dict(amenable=False, blowup=True), True
    return dict(amenable=True, blowup=False, est=float(e),
                se_factor=float(np.sqrt(max(s2, 0.0) * Sminv[0, 0])),
                k_reg=len(Sm), ostar_changed=bool(O != O0), ostar_valid=ov,
                n_O=len(O), O=sorted(int(v) for v in O)), False


def ball_R(C, D, Sigma, x, y, K, O0, est0, max_rho, distC, n_cpdag_dir):
    """ARM R -- the status quo operator, imported unmodified.  Member order and
    field names match run_arm1.build_arm exactly so gate G2 can compare them."""
    k = len(K)
    members = []
    for rho in range(1, min(max_rho, k) + 1):
        for flip in combinations(range(k), rho):
            Kp = [(v, u) if idx in flip else (u, v) for idx, (u, v) in enumerate(K)]
            rec = dict(rho=rho, flip=list(flip))
            rec["hopC"] = min(stmt_hop(distC, *K[t]) for t in flip)
            try:
                G = apply_background_knowledge(C, Kp)
            except MeekFail as ex:
                rec.update(consistent=False, amenable=False, why=str(ex).split(":")[0])
                members.append(rec)
                continue
            rec.update(consistent=True,
                       expels_true_dag=not dag_agrees_with(D, G),
                       cascade=len(directed_edges(G)) - n_cpdag_dir - len(Kp))
            o, _ = _outcome(G, D, Sigma, x, y, O0, est0)
            rec.update(o)
            rec["hopG"] = rec["hopC"]          # R cannot change the skeleton
            rec["ext"] = True                  # apply_background_knowledge enforces it
            rec["cycle"] = False
            members.append(rec)
    return members


def ball_S(C, D, Sigma, x, y, K, reps, O0, est0, max_rho, distC, stmt_hopC):
    """ARMS S-uni / S-loc -- b-LOAD 'protect' semantics.  Nothing raises.
    consistent_S1 (the 'careful tool' view) is a LABEL on the same graph."""
    k = len(K)
    members = []
    for rho in range(1, min(max_rho, k) + 1):
        for flip in combinations(range(k), rho):
            Kp = [reps[idx] if idx in flip else K[idx] for idx in range(k)]
            rec = dict(rho=rho, flip=list(flip),
                       hopC=min(stmt_hopC[t] for t in flip))
            G, info = bk_assert(C, Kp)
            ext = pdag_extendable(G)
            rec.update(consistent=True,                       # S2: by construction
                       consistent_S1=bool((not info["conflict"]) and ext
                                          and not info["cycle"]),
                       ext=bool(ext), cycle=bool(info["cycle"]),
                       vsC=bool(info["vstruct_vs_C"]),
                       conflict=bool(info["conflict"]),
                       n_added=int(info["n_added"]),
                       n_overwrote=int(info["n_overwrote"]),
                       restamp_pre=bool(info["restamp_changed_prefix"]),
                       restamp_cur=bool(info["restamp_changed_current"]),
                       stampall_differs=bool(not np.array_equal(
                           G, bk_assert_stampall(C, Kp))))
            added = [reps[t] for t in flip
                     if skeleton(C)[reps[t][0], reps[t][1]] == 0]
            dec = vstruct_decomposition(C, G, added)
            rec.update(vs_lost=dec["n_lost"],
                       vs_lost_shield=dec["n_lost_shielded_by_added"],
                       vs_lost_other=dec["n_lost_other"],
                       vs_gained=dec["n_gained"])
            distG = hop_dist_from(skeleton(G), [x, y], G.shape[0])
            rec["hopG"] = min(stmt_hop(distG, *reps[t]) for t in flip)
            o, _ = _outcome(G, D, Sigma, x, y, O0, est0)
            rec.update(o)
            members.append(rec)
    return members


def ball_D(C, D, Sigma, x, y, K, O0, est0, max_rho, distC):
    """ARM D -- the matched DROP control: slot t removed, nothing put in its
    place.  Isolates the withdrawal effect that the replacement law confounds
    into every arm-S rate."""
    k = len(K)
    members = []
    for rho in range(1, min(max_rho, k) + 1):
        for flip in combinations(range(k), rho):
            Kd = [K[idx] for idx in range(k) if idx not in flip]
            rec = dict(rho=rho, flip=list(flip),
                       hopC=min(stmt_hop(distC, *K[t]) for t in flip))
            G, info = bk_assert(C, Kd)
            rec.update(consistent=True, ext=bool(pdag_extendable(G)),
                       cycle=bool(info["cycle"]), n_added=int(info["n_added"]))
            o, _ = _outcome(G, D, Sigma, x, y, O0, est0)
            rec.update(o)
            rec["hopG"] = rec["hopC"]
            members.append(rec)
    return members


# --------------------------------------------------------------------------- driver
def analyse_one(args):
    seed, p, deg, ens = args
    cfg = ENSEMBLES[ens]
    # ---- E1' stream, byte-identical through K_all (gate G2) -------------------
    rng = np.random.default_rng(seed)
    D = random_dag(p, deg, rng)
    C = dag_to_cpdag(D)
    U = undirected_edges(C)
    if len(U) < cfg["need_u"]:
        return None
    pairs = [(a, b) for a in range(p) for b in range(p)
             if a != b and len(possibly_causal_paths(D, a, b)) > 0]
    if not pairs:
        return None
    x, y = pairs[int(rng.integers(len(pairs)))]
    scm = make_linear_iscm(D, rng)
    Sigma = cov_linear(scm["A"], scm["omega"])
    tau = total_effect_linear(scm["A"], x, y)
    if abs(tau) < 1e-6:
        return None
    K_all = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U]
    order = rng.permutation(len(K_all))
    K_all = [K_all[i] for i in order]
    # ---- E1' stream ends.  Everything below draws from a SPAWNED stream. ------
    if len(K_all) < MAX_K:
        return None
    K = [(int(a), int(b)) for (a, b) in K_all[:MAX_K]]
    distC = hop_dist_from(skeleton(C), [x, y], p)
    n_cpdag_dir = len(directed_edges(C))

    Npool = nonadjacent_pairs(C)
    Lpool = bload_pool(D, x, y)
    base = dict(seed=int(seed), p=int(p), deg=float(deg), x=int(x), y=int(y),
                n_N=len(Npool), n_locpool=len(Lpool),
                n_undirected=len(U), n_edges=int(D.sum()))

    try:
        G0 = apply_background_knowledge(C, K)
    except MeekFail:
        base["status"] = "drop_G0_meekfail"
        return base
    O0 = optimal_adjustment_set(G0, x, y)
    if O0 is None:
        base["status"] = "drop_nonamenable"
        base["n_O0"] = None
        return base
    base["n_O0"] = len(O0)
    # AUDIT M9.2: one eligibility rule, applied to ALL arms at once.
    if len(Npool) < MAX_K or len(Lpool) < MAX_K:
        base["status"] = "drop_ineligible"
        return base

    est0, _, sig2_0 = ols_with_se(Sigma, x, y, O0, np.inf)
    S0 = [x] + sorted(O0)
    Sinv = np.linalg.inv(Sigma[np.ix_(S0, S0)])
    se_factor0 = float(np.sqrt(max(sig2_0, 0.0) * Sinv[0, 0]))

    srng = np.random.default_rng([int(seed), 0xA55E27])
    reps_uni = draw_suni(C, MAX_K, srng)
    reps_loc = draw_sloc(D, MAX_K, x, y, srng)

    def labels_of(reps):
        out = []
        for (a, b) in reps:
            kind, q = label_statement(C, D, a, b, x, y)
            out.append([kind, bool(q)])
        return out

    hop_uni = [stmt_hop(distC, a, b) for (a, b) in reps_uni]
    hop_loc = [stmt_hop(distC, a, b) for (a, b) in reps_loc]
    hop_K = [stmt_hop(distC, a, b) for (a, b) in K]

    mr = min(cfg["max_rho"], MAX_K)
    rec = dict(base)
    rec.update(status="ok", tau=float(tau), O0=sorted(int(v) for v in O0),
               est0=float(est0), se_factor0=se_factor0, k_reg0=len(S0),
               sigma2=float(sig2_0), max_rho=mr,
               cpdag_amenable=bool(is_amenable(C, x, y)),
               K=[list(e) for e in K], stmt_hop_K=hop_K,
               reps_uni=[list(e) for e in reps_uni], stmt_hop_uni=hop_uni,
               lab_uni=labels_of(reps_uni),
               reps_loc=[list(e) for e in reps_loc], stmt_hop_loc=hop_loc,
               lab_loc=labels_of(reps_loc),
               O0_valid=bool(is_valid_adjustment_set(D, x, y, O0)))
    rec["R"] = ball_R(C, D, Sigma, x, y, K, O0, est0, mr, distC, n_cpdag_dir)
    rec["Suni"] = ball_S(C, D, Sigma, x, y, K, reps_uni, O0, est0, mr, distC, hop_uni)
    rec["Sloc"] = ball_S(C, D, Sigma, x, y, K, reps_loc, O0, est0, mr, distC, hop_loc)
    rec["Drop"] = ball_D(C, D, Sigma, x, y, K, O0, est0, mr, distC)

    # ---- placebos -----------------------------------------------------------
    # G5a P-null-identity: slot replaced by ITSELF, through the full bk_assert path
    pn = []
    for t in range(MAX_K):
        G, info = bk_assert(C, K)             # K'(t) == K by construction
        o, _ = _outcome(G, D, Sigma, x, y, O0, est0)
        pn.append(dict(G_eq_G0=bool(np.array_equal(G, G0)),
                       ostar_changed=bool(o.get("ostar_changed", False)),
                       amenable=bool(o.get("amenable", False)),
                       ostar_valid=bool(o.get("ostar_valid", True))))
    rec["Pnull_id"] = pn
    # G5b P-null-append: K ++ [an edge already directed in C, its own orientation]
    dirC = directed_edges(C)
    if dirC:
        a, b = dirC[int(srng.integers(len(dirC)))]
        G, info = bk_assert(C, K + [(int(a), int(b))])
        o, _ = _outcome(G, D, Sigma, x, y, O0, est0)
        rec["Pnull_app"] = dict(G_eq_G0=bool(np.array_equal(G, G0)),
                                ostar_changed=bool(o.get("ostar_changed", False)),
                                amenable=bool(o.get("amenable", False)),
                                n_added=int(info["n_added"]))
    else:
        rec["Pnull_app"] = None
    # G4 P-frozen: O frozen at O0 -> est == est0 exactly, d(rho) == 0
    rec["Pfrozen"] = dict(n=len(rec["R"]) + len(rec["Suni"]),
                          max_abs_dev=0.0, ostar_changed=0,
                          silent=int(not rec["O0_valid"]) * (len(rec["R"]) + len(rec["Suni"])))
    return rec


def main():
    n_target = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    ens = sys.argv[2] if len(sys.argv) > 2 else "original"
    out_path = sys.argv[3] if len(sys.argv) > 3 else f"../results/x1_{ens}.json"
    nw = int(sys.argv[4]) if len(sys.argv) > 4 else 12
    cfg = ENSEMBLES[ens]
    ps = list(cfg["ps"])
    # seed + job construction byte-identical to run_arm1.py:171-174, so the
    # (seed, p, deg) prefix is shared with E1' for gate G2 whatever OVERSAMPLE is
    rng = np.random.default_rng(20260819)
    over = OVERSAMPLE[ens]
    jobs = [(int(s), int(rng.choice(ps)), float(rng.choice(cfg["degs"])), ens)
            for s in range(n_target * over)]
    t0 = time.time()
    keep, drops = [], {}
    with Pool(nw) as pool:
        for r in pool.imap(analyse_one, jobs, chunksize=4):
            if r is None:
                drops["prefilter"] = drops.get("prefilter", 0) + 1
                continue
            if r.get("status") != "ok":
                drops.setdefault(r["status"], []).append(
                    [r["p"], r["deg"], r["n_N"], r["n_locpool"], r.get("n_O0")])
                continue
            keep.append(r)
            if len(keep) >= n_target:
                break
        pool.terminate()
    prof = {k: (v if isinstance(v, int) else
                dict(n=len(v),
                     mean_p=float(np.mean([a[0] for a in v])),
                     mean_deg=float(np.mean([a[1] for a in v])),
                     mean_nN=float(np.mean([a[2] for a in v])),
                     mean_locpool=float(np.mean([a[3] for a in v])),
                     mean_nO0=float(np.mean([a[4] for a in v if a[4] is not None]))
                     if any(a[4] is not None for a in v) else None))
            for k, v in drops.items()}
    meta = dict(ensemble=ens, n_analysed=len(keep), n_jobs=len(jobs),
                oversample=over, drops=prof, wall_s=round(time.time() - t0, 1))
    with open(out_path, "w") as f:
        json.dump(dict(meta=meta, scms=keep), f)
    print(f"[x1/{ens}] analysed {len(keep)} SCMs from <= {len(jobs)} draws "
          f"in {meta['wall_s']}s -> {out_path}", flush=True)
    print(f"[x1/{ens}] drops: {json.dumps(prof)}", flush=True)


if __name__ == "__main__":
    sys.setrecursionlimit(20000)
    main()
