"""
G2  reversal-arm reproduction: on SCM seeds shared with E1', ARM R must match
    E1' member-for-member on consistent / amenable / ostar_changed / ostar_valid
    and est to 1e-12.  >= 500 shared SCMs, 100%.
G11 Dor-Tarsi in-run: on every member with |undirected_edges(G)| <= 12, check
    pdag_extendable(G) == (a brute-force consistent extension exists).  100%.
    Requires regenerating the graphs from their seeds -- the run stores outcomes,
    not graphs.
"""
import json
import sys
from itertools import combinations
from multiprocessing import Pool

import numpy as np

sys.setrecursionlimit(20000)
from adjust import optimal_adjustment_set, possibly_causal_paths
from graphs import (MeekFail, apply_background_knowledge, consistent_dag_extensions,
                    dag_to_cpdag, random_dag, skeleton, undirected_edges,
                    v_structures)
from scm import make_linear_iscm
from x1_ops import bk_assert, draw_suni, pdag_extendable


def g2(ens, e1path):
    x1 = json.load(open(f"../results/x1_{ens}.json"))["scms"]
    e1 = json.load(open(e1path))
    e1by = {r["seed"]: r for r in e1 if "K4" in r.get("arms", {})
            and r["arms"]["K4"].get("mpdag_amenable")}
    shared = ncmp = bad = 0
    badfields = {}
    for s in x1:
        r = e1by.get(s["seed"])
        if r is None:
            continue
        a = r["arms"]["K4"]
        if [list(e) for e in a["K"]] != s["K"]:
            badfields["K_mismatch"] = badfields.get("K_mismatch", 0) + 1
            continue
        shared += 1
        for me, mx in zip(a["members"], s["R"]):
            ncmp += 1
            for f in ("consistent", "amenable", "ostar_changed", "ostar_valid"):
                if bool(me.get(f, False)) != bool(mx.get(f, False)):
                    badfields[f] = badfields.get(f, 0) + 1
                    bad += 1
            if "est" in me and "est" in mx and abs(me["est"] - mx["est"]) > 1e-12:
                badfields["est"] = badfields.get("est", 0) + 1
                bad += 1
    return dict(ensemble=ens, shared_scms=shared, members_compared=ncmp,
                mismatches=bad, by_field=badfields,
                verdict="PASS" if (bad == 0 and shared >= 500) else "FAIL")


def _g11_one(args):
    seed, p, deg = args
    rng = np.random.default_rng(seed)
    D = random_dag(p, deg, rng)
    C = dag_to_cpdag(D)
    U = undirected_edges(C)
    if len(U) < 3:
        return (0, 0)
    pairs = [(a, b) for a in range(p) for b in range(p)
             if a != b and len(possibly_causal_paths(D, a, b)) > 0]
    if not pairs:
        return (0, 0)
    x, y = pairs[int(rng.integers(len(pairs)))]
    scm = make_linear_iscm(D, rng)
    K_all = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U]
    order = rng.permutation(len(K_all))
    K_all = [K_all[i] for i in order]
    if len(K_all) < 4:
        return (0, 0)
    K = [(int(a), int(b)) for (a, b) in K_all[:4]]
    srng = np.random.default_rng([int(seed), 0xA55E27])
    reps = draw_suni(C, 4, srng)
    if reps is None:
        return (0, 0)
    nchk = nbad = 0
    for rho in range(1, 5):
        for flip in combinations(range(4), rho):
            Kp = [reps[i] if i in flip else K[i] for i in range(4)]
            G, _ = bk_assert(C, Kp)
            if len(undirected_edges(G)) > 12:
                continue
            nchk += 1
            fast = pdag_extendable(G)
            slow = len(consistent_dag_extensions(G, ref_vstructs=v_structures(G),
                                                 limit=1)) > 0
            if fast != slow:
                nbad += 1
    return (nchk, nbad)


def g11(ens, n_scm, nw=12):
    from run_x1 import ENSEMBLES, OVERSAMPLE
    cfg = ENSEMBLES[ens]
    ps = list(cfg["ps"])
    rng = np.random.default_rng(20260819)
    jobs = [(int(s), int(rng.choice(ps)), float(rng.choice(cfg["degs"])))
            for s in range(n_scm)]
    with Pool(nw) as pool:
        res = pool.map(_g11_one, jobs, chunksize=4)
    nchk = sum(a for a, _ in res); nbad = sum(b for _, b in res)
    return dict(ensemble=ens, members_checked=nchk, disagreements=nbad,
                verdict="PASS" if nbad == 0 else "FAIL",
                note="|undirected_edges(G)| <= 12; oracle = "
                     "consistent_dag_extensions(G, ref_vstructs=v_structures(G), limit=1)")


if __name__ == "__main__":
    out = dict(
        G2=[g2("original", "/home/costaj/latent-causal/e1prime-se/results/arm1_original.json"),
            g2("licensed", "/home/costaj/latent-causal/e1prime-se/results/arm1_licensed.json")],
        G11=[g11("original", 3000), g11("large", 400)])
    print(json.dumps(out, indent=1))
    json.dump(out, open("../results/x1_gates_g2_g11.json", "w"), indent=1)
