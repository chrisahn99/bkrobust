"""
GATE G3, in its strongest form, plus the M15 split for ARM B.

E1' stores (seed, p, deg) per SCM, and `random_dag` consumes a deterministic
prefix of default_rng(seed), so D and C are regenerable WITHOUT re-running the
expensive iSCM standardisation.  The archive also stores x, y, K, O0 and every
member's flags.  So every structural claim in the archive can be re-derived by
an INDEPENDENT implementation and compared, member by member:

    n_edges, n_undirected, cpdag_amenable, stmt_hopdist, O0,
    consistent, amenable, ostar_changed, ostar_valid, cascade, expels_true_dag

Only `est` needs Sigma and is not re-derived here (Step 0 uses it directly).

At the same time this produces what the archive cannot: the M15 decomposition
of E_ident into `nopath` (the closure destroyed every proper possibly-causal
path -- the analyst reads "no causal effect", a WRONG ANSWER, not a lost one)
and `unamenable` (paths remain but one starts undirected).

Usage: python verify_armB.py <arch_dir> <out.json> [n_workers]
"""
import json
import sys
from collections import defaultdict
from itertools import combinations
from multiprocessing import Pool

import numpy as np

import x2lib as X
from adjust import is_valid_adjustment_set, optimal_adjustment_set, is_amenable
from graphs import (MeekFail, apply_background_knowledge, dag_agrees_with,
                    dag_to_cpdag, directed_edges, random_dag, undirected_edges)

ARCH = sys.argv[1] if len(sys.argv) > 1 else "/home/costaj/latent-causal/e1prime-se/results"
OUT = sys.argv[2] if len(sys.argv) > 2 else "../results/verify_armB.json"
NW = int(sys.argv[3]) if len(sys.argv) > 3 else 4
TAU_FLOOR = 1e-3
CELLS = [("original", "K4", 4), ("licensed", "K4", 4), ("large", "K4", 4),
         ("k8", "K4", 3), ("k8", "K6", 3), ("k8", "K8", 3)]


def one(args):
    rec, key, max_rho = args
    arm = rec["arms"].get(key)
    if arm is None or not arm.get("mpdag_amenable"):
        return None
    p, deg, seed = rec["p"], rec["deg"], rec["seed"]
    x, y = rec["x"], rec["y"]
    rng = np.random.default_rng(seed)
    D = random_dag(p, deg, rng)
    C = dag_to_cpdag(D)
    mism = defaultdict(int)
    if int(D.sum()) != rec["n_edges"]:
        mism["n_edges"] += 1
    if len(undirected_edges(C)) != rec["n_undirected"]:
        mism["n_undirected"] += 1
    if bool(is_amenable(C, x, y)) != bool(rec["cpdag_amenable"]):
        mism["cpdag_amenable"] += 1
    K = [tuple(e) for e in arm["K"]]
    dist = X.hop_dist_inf(C, [x, y])
    hop = []
    for (u, v) in K:
        dmin, _ = X.stmt_dist(dist, (u, v))
        hop.append(dmin)
    arch_hop = [float(h) if h < X.SENTINEL else np.inf for h in arm["stmt_hopdist"]]
    if [h if np.isfinite(h) else np.inf for h in hop] != arch_hop:
        mism["stmt_hopdist"] += 1
    G0 = apply_background_knowledge(C, K)
    O0, np0, amen0 = X.ostar_and_paths(G0, x, y)
    if O0 is None or sorted(O0) != arm["O0"]:
        mism["O0"] += 1
        return dict(mism=dict(mism), split=None, n_members=0)
    ndirC = len(directed_edges(C))
    split = defaultdict(lambda: defaultdict(int))
    nm = 0
    k = len(K)
    Oset_full, Oset_prune = {O0}, {r: {O0} for r in (-1, 0, 1, 2, 3)}
    for rho in range(1, min(max_rho, k) + 1):
        for flip in combinations(range(k), rho):
            m = arm["members"][nm]
            nm += 1
            Kp = [(v, u) if i in flip else (u, v) for i, (u, v) in enumerate(K)]
            try:
                G = apply_background_knowledge(C, Kp)
                cons = True
            except MeekFail:
                cons = False
            if cons != bool(m.get("consistent")):
                mism["consistent"] += 1
                continue
            if not cons:
                continue
            O, npth, amen = X.ostar_and_paths(G, x, y)
            if bool(amen) != bool(m.get("amenable")):
                mism["amenable"] += 1
            casc = len(directed_edges(G)) - ndirC - len(Kp)
            if casc != m.get("cascade"):
                mism["cascade"] += 1
            if (not dag_agrees_with(D, G)) != m.get("expels_true_dag"):
                mism["expels_true_dag"] += 1
            if amen:
                if (O != O0) != bool(m.get("ostar_changed")):
                    mism["ostar_changed"] += 1
                if bool(is_valid_adjustment_set(D, x, y, O)) != bool(m.get("ostar_valid")):
                    mism["ostar_valid"] += 1
                Oset_full.add(O)
                for r in (-1, 0, 1, 2, 3):
                    near = set(i for i in range(k)
                               if r >= 0 and np.isfinite(hop[i]) and hop[i] <= r)
                    if set(flip) <= near:
                        Oset_prune[r].add(O)
            else:
                # M15: the split the archive cannot carry
                hs = [hop[i] for i in flip]
                b = X.dbucket(min(hs))
                allfar = all(np.isfinite(h) and h >= 1 for h in hs)
                lab = f"rho{rho}|{b}|{'allfar' if allfar else 'mixed'}"
                split[lab]["N_ident"] += 1
                split[lab]["nopath" if npth == 0 else "unamen"] += 1
    qset = {str(r): bool(Oset_prune[r] != Oset_full) for r in (-1, 0, 1, 2, 3)}
    return dict(mism=dict(mism), split={a: dict(b) for a, b in split.items()},
                n_members=nm, qset=qset, n_Ostar=len(Oset_full))


def main():
    out = {}
    for ens, key, max_rho in CELLS:
        recs = json.load(open(f"{ARCH}/arm1_{ens}.json"))
        jobs = [(r, key, max_rho) for r in recs if abs(r["tau"]) >= TAU_FLOOR]
        tot_mism = defaultdict(int)
        tot_split = defaultdict(lambda: defaultdict(int))
        n_scm = n_mem = 0
        qs = defaultdict(int)
        nOs = []
        with Pool(NW) as pool:
            for r in pool.imap_unordered(one, jobs, chunksize=8):
                if r is None:
                    continue
                n_scm += 1
                n_mem += r["n_members"]
                nOs.append(r.get("n_Ostar", 1))
                for rr, ch in (r.get("qset") or {}).items():
                    qs[rr] += int(ch)
                for a, b in r["mism"].items():
                    tot_mism[a] += b
                for a, d in (r["split"] or {}).items():
                    for kk, vv in d.items():
                        tot_split[a][kk] += vv
        out[f"{ens}/{key}"] = dict(
            n_scm=n_scm, n_members=n_mem, mismatches=dict(tot_mism),
            ident_split={a: dict(b) for a, b in sorted(tot_split.items())},
            Qset={rr: dict(n=n_scm, n_changed=qs[rr], Q=qs[rr] / max(n_scm, 1),
                           Q_ci=list(X.clopper_pearson(qs[rr], n_scm)))
                  for rr in ("-1", "0", "1", "2", "3")},
            mean_n_Ostar=float(np.mean(nOs)) if nOs else 0.0)
        print(f"[G3] {ens}/{key}: {n_scm} SCMs, {n_mem} members, "
              f"mismatches={dict(tot_mism) or 0}", flush=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1, default=float)
    print(f"[G3] -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
