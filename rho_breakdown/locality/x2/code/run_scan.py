"""
X2 ARM A (mechanism / E_gain, per M2) and ARM E (the primary S5 arm, per M3).

Both run on the SAME SCM draw, so the expensive part (random_dag ->
dag_to_cpdag -> make_linear_iscm with N_REF=100k) is paid once.

ARM A   base K = {} , so G0 = C.  Every undirected edge of C, both
        orientations.  This is PREREG 2.1 verbatim.  The AUDIT (A2) measured it
        to be null by construction; it is kept for (i) E_gain, which is the only
        event it can produce, and (ii) the mechanism statistics (amen(C), |cn|,
        undirected edges incident to cn) that explain why it is null.

ARM E   base K = k_base in 1..4 TRUE statements drawn as run_arm1.py:139, then
        EVERY undirected edge of C swept in BOTH orientations (replace in place
        if the edge is already in the base, else append).  This restores the
        dynamic range ARM A lacks while removing E1's |K| = 4 sampling of which
        edges may be misstated -- M3.

        A swept TRUE statement on an edge already in the base is a NULL TRIAL
        (K' is K_base reordered; T7c proves M(K') = M(K_base)).  Those trials
        must produce zero events and are used as an in-run placebo, gate G11.

Usage: python run_scan.py <n_scm> <ensemble> <out.json> [n_workers]
"""
import json
import sys
import time
from collections import defaultdict
from multiprocessing import Pool

import numpy as np

import x2lib as X
from adjust import causal_nodes, is_amenable, possibly_causal_paths
from graphs import is_undirected, undirected_edges

ENSEMBLES = {
    # id                 ps            degs                      need_u  query rule
    "licensed": dict(ps=range(5, 11),  degs=[1.5, 2.0, 2.5, 4.0, 6.0], need_u=3, q="uniform"),
    "large":    dict(ps=range(15, 26), degs=[1.5, 2.0, 2.5],           need_u=3, q="uniform"),
    "k8":       dict(ps=range(5, 11),  degs=[1.5, 2.0, 2.5, 4.0, 6.0], need_u=8, q="uniform"),
    "deep":     dict(ps=range(12, 21), degs=[1.5, 2.0, 2.5],           need_u=4, q="argmax_delta_C"),
    "deep8":    dict(ps=range(12, 21), degs=[1.5, 2.0, 2.5],           need_u=8, q="argmax_delta_C"),
}
JOB_SEED = {"licensed": 20260820, "large": 20260820, "k8": 20260820,
            "deep": 20260821, "deep8": 20260821}
OVERSAMPLE = {"licensed": 4, "large": 4, "k8": 70, "deep": 12, "deep8": 130}

COUNTERS = ("D_all", "D_con", "D_am", "D_s5", "N_meek", "N_set", "N_bias",
            "N_ident", "N_ident_nopath", "N_ident_unamen", "N_gain",
            "N_gain_paths", "N_gain_orient", "N_est_inf", "N_est_n",
            "cascade_pos", "N_null")
MAX_FAR_RECORDS = 8          # per SCM per arm, for E_ident / E_gain
MAX_SET_RECORDS = 50         # per SCM per arm, for E_set (the money case)


def sweep(sd, K_base, tag, acc, far_recs, ens):
    """One arm.  Mutates acc (counter dict) and far_recs (list)."""
    C, D, Sigma, x, y = sd["C"], sd["D"], sd["Sigma"], sd["x"], sd["y"]
    base = X.make_base(C, D, Sigma, x, y, K_base)
    if base is None:                       # impossible for an all-true base (T7a)
        return None
    dist = X.hop_dist_inf(C, [x, y])
    U = undirected_edges(C)
    base_pairs = {frozenset(t) for t in K_base}
    n_gain_far = n_id_far = n_set_far = 0
    near = dict(D_am=0, N_set=0)
    far = dict(D_am=0, N_set=0)
    for (u, v) in U:
        s_true = (u, v) if D[u, v] == 1 else (v, u)
        s_false = (s_true[1], s_true[0])
        dmin, dmax = X.stmt_dist(dist, (u, v))
        bmin, bmax = X.dbucket(dmin), X.dbucket(dmax)
        in_base = frozenset((u, v)) in base_pairs
        for s, tin in ((s_true, True), (s_false, False)):
            Kp = [t for t in K_base if frozenset(t) != frozenset((u, v))] + [s]
            null = in_base and tin
            rec = X.classify_trial(C, D, Sigma, x, y, Kp, base)
            op = "rep" if in_base else "add"
            kmin = (bmin, "T" if tin else "F", op)
            kmax = (bmax, "T" if tin else "F", op)
            a, ax = acc["dmin"][kmin], acc["dmax"][kmax]
            for h in (a, ax):
                h["D_all"] += 1
            if null:
                a["N_null"] += 1
            if rec["meek_fail"]:
                a["N_meek"] += 1
                ax["N_meek"] += 1
                continue
            for h in (a, ax):
                h["D_con"] += 1
            if rec["cascade_pos"]:
                a["cascade_pos"] += 1
            if rec["e_ident"]:
                a["N_ident"] += 1
                ax["N_ident"] += 1
                a["N_ident_nopath"] += int(rec["e_ident_nopath"])
                a["N_ident_unamen"] += int(rec["e_ident_unamen"])
                if X.far_finite(dmin) and n_id_far < MAX_FAR_RECORDS:
                    n_id_far += 1
                    far_recs.append(dict(arm=tag, ens=ens, seed=sd["seed"], p=sd["p"],
                                         deg=sd["deg"], x=x, y=y, s=list(s),
                                         dmin=dmin, dmax=dmax, true_in_D=tin,
                                         event="ident_nopath" if rec["e_ident_nopath"] else "ident_unamen",
                                         cascade=rec["cascade"], k_base=len(K_base)))
            elif rec["e_gain"]:
                a["N_gain"] += 1
                ax["N_gain"] += 1
                a["N_gain_paths"] += int(rec["e_gain_paths"])
                a["N_gain_orient"] += int(rec["e_gain_orient"])
                if X.far_finite(dmin) and n_gain_far < MAX_FAR_RECORDS:
                    n_gain_far += 1
                    far_recs.append(dict(arm=tag, ens=ens, seed=sd["seed"], p=sd["p"],
                                         deg=sd["deg"], x=x, y=y, s=list(s),
                                         dmin=dmin, dmax=dmax, true_in_D=tin,
                                         event="gain_paths" if rec["e_gain_paths"] else "gain_orient",
                                         cascade=rec["cascade"], k_base=len(K_base)))
            elif rec["amenable0"] and rec["amenable1"]:
                for h in (a, ax):
                    h["D_am"] += 1
                if (not tin) and np.isfinite(dmin):
                    a["D_s5"] += 1
                    ax["D_s5"] += 1
                if dmin == 0:
                    near["D_am"] += 1
                elif np.isfinite(dmin):
                    far["D_am"] += 1
                if rec["e_set"]:
                    a["N_set"] += 1
                    ax["N_set"] += 1
                    a["N_bias"] += int(rec["e_bias"])
                    ax["N_bias"] += int(rec["e_bias"])
                    if dmin == 0:
                        near["N_set"] += 1
                    elif np.isfinite(dmin):
                        far["N_set"] += 1
                    if X.far_finite(dmin) and n_set_far < MAX_SET_RECORDS:
                        n_set_far += 1
                        far_recs.append(dict(arm=tag, ens=ens, seed=sd["seed"], p=sd["p"],
                                             deg=sd["deg"], x=x, y=y, s=list(s),
                                             dmin=dmin, dmax=dmax, true_in_D=tin,
                                             event="set", e_bias=rec["e_bias"],
                                             O0=sorted(base["O0"]), O1=sorted(rec["ostar"]),
                                             est0=base["est0"], est=rec["est"],
                                             tau=sd["tau"], cascade=rec["cascade"],
                                             K_base=[list(t) for t in K_base],
                                             k_base=len(K_base)))
                a["N_est_inf"] += int(rec["e_est_inf"])
                a["N_est_n"] += int(rec["e_est_n"])
    return dict(near=near, far=far, amen0=base["amen0"],
                ncn=len(causal_nodes(base["G0"], x, y)) if base["amen0"] else 0)


def analyse_one(args):
    seed, p, deg, ens = args
    cfg = ENSEMBLES[ens]
    t0 = time.time()
    try:
        sd = X.draw_scm(seed, p, deg, cfg["need_u"], cfg["q"])
    except (X.PathCapExceeded, RecursionError):
        return dict(abort=True)
    if sd is None:
        return None
    C, D, x, y = sd["C"], sd["D"], sd["x"], sd["y"]
    acc = {"A": {"dmin": defaultdict(lambda: defaultdict(int)),
                 "dmax": defaultdict(lambda: defaultdict(int))},
           "E": {"dmin": defaultdict(lambda: defaultdict(int)),
                 "dmax": defaultdict(lambda: defaultdict(int))}}
    far_recs = []
    try:
        rA = sweep(sd, [], "A", acc["A"], far_recs, ens)
        rng = sd["rng"]
        kb = int(rng.integers(1, min(4, len(sd["K_all"])) + 1))
        K_base = sd["K_all"][:kb]
        rE = sweep(sd, K_base, "E", acc["E"], far_recs, ens)
    except (X.PathCapExceeded, RecursionError):
        return dict(abort=True)
    if time.time() - t0 > X.SCM_CAP:
        return dict(abort=True)
    # ARM A mechanism statistics (AUDIT A2)
    amenC = bool(is_amenable(C, x, y))
    ncn = uinc = 0
    if amenC:
        cn = causal_nodes(C, x, y)
        ncn = len(cn)
        uinc = sum(1 for (a, b) in undirected_edges(C) if a in cn or b in cn)
    out = dict(seed=int(seed), p=int(p), deg=float(deg), x=int(x), y=int(y),
               tau=float(sd["tau"]), nU=len(sd["U"]), n_edges=sd["n_edges"],
               amen_C=amenC, ncn_C=int(ncn), u_inc_cn_C=int(uinc),
               k_base=int(kb), overlap=sd["overlap"],
               amen0_A=bool(rA["amen0"]) if rA else None,
               amen0_E=bool(rE["amen0"]) if rE else None,
               near_A=rA["near"] if rA else None, far_A=rA["far"] if rA else None,
               near_E=rE["near"] if rE else None, far_E=rE["far"] if rE else None,
               far_recs=far_recs,
               acc={arm: {d: {"|".join(k): dict(v) for k, v in acc[arm][d].items()}
                          for d in ("dmin", "dmax")} for arm in ("A", "E")})
    return out


def main():
    n_scm = int(sys.argv[1])
    ens = sys.argv[2]
    out_path = sys.argv[3]
    nw = int(sys.argv[4]) if len(sys.argv) > 4 else 12
    cfg = ENSEMBLES[ens]
    ps = list(cfg["ps"])
    rng = np.random.default_rng(JOB_SEED[ens])
    jobs = [(int(s), int(rng.choice(ps)), float(rng.choice(cfg["degs"])), ens)
            for s in range(n_scm * OVERSAMPLE[ens])]
    t0 = time.time()
    res, n_abort = [], 0
    with Pool(nw) as pool:
        for r in pool.imap_unordered(analyse_one, jobs, chunksize=4):
            if r is None:
                continue
            if r.get("abort"):
                n_abort += 1
                continue
            res.append(r)
            if len(res) >= n_scm:
                break
        pool.terminate()
    with open(out_path, "w") as f:
        json.dump(dict(ens=ens, n_scm=len(res), n_abort=n_abort,
                       n_draws=len(jobs), recs=res), f, default=float)
    print(f"[scan/{ens}] kept {len(res)} SCMs, {n_abort} aborts, "
          f"<= {len(jobs)} draws, {time.time()-t0:.1f}s -> {out_path}", flush=True)


if __name__ == "__main__":
    main()
