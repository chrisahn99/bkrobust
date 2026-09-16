"""
X3 -- SKELETON PERTURBATION.  Does the locality law survive when the DISCOVERY is wrong?

Pre-registered in ../PREREG.md before the first number.  The research question is the one
posted to the Teams channel on 2026-08-31: the assumption that the data-driven method is
correct is very strong.  Every measurement in this line so far uses the ORACLE CPDAG.

Model of a skeleton error (an assumption, stated, not derived):
    the analyst's graph is  Chat = dag_to_cpdag(D')  with  D' = D minus one edge.

Arms:
    SK   G0 = Meek(Chat, K|)       -- PRIMARY: the same error inside our actual setting
    S0   G0 = Chat                 -- SECONDARY: the discovery error alone (see PREREG amendment 7)

Outcome, evaluated against the TRUE DAG D:
    meekfail / not_amenable / unchanged / changed_valid / silent_bias

Usage: python run_x3.py <n_analysed_target> <ensemble> <out.json> [n_workers]
"""
import json
import sys
import time
from multiprocessing import Pool

import numpy as np

from adjust import (PathBlowup, cov_linear, is_amenable,
                    is_valid_adjustment_set, optimal_adjustment_set,
                    ols_coefficient, possibly_causal_paths,
                    total_effect_linear)
from graphs import (MeekFail, apply_background_knowledge, dag_to_cpdag,
                    is_undirected, random_dag, skeleton, undirected_edges)
from scm import make_linear_iscm
from x1_ops import hop_dist_from, stmt_hop

MAX_K = 4

# grids imported UNCHANGED from run_x1.py:51-53 (gate G2)
ENSEMBLES = {
    "original": dict(ps=range(5, 9),   degs=[1.5, 2.0, 2.5],           need_u=3),
    "licensed": dict(ps=range(5, 11),  degs=[1.5, 2.0, 2.5, 4.0, 6.0], need_u=3),
    "large":    dict(ps=range(15, 26), degs=[1.5, 2.0, 2.5],           need_u=3),
}
OVERSAMPLE = {"original": 6, "licensed": 8, "large": 4}


def _outcome(G, D, Sigma, x, y, O0, tau):
    """The estimand side.  Returns a dict; never raises."""
    try:
        O = optimal_adjustment_set(G, x, y)
    except PathBlowup:
        return dict(outcome="blowup")
    if O is None:
        return dict(outcome="not_amenable")
    if O == O0:
        return dict(outcome="unchanged", n_O=len(O))
    try:
        valid = bool(is_valid_adjustment_set(D, x, y, O))
        est = float(ols_coefficient(Sigma, x, y, O))
    except (PathBlowup, np.linalg.LinAlgError):
        return dict(outcome="blowup")
    return dict(outcome="changed_valid" if valid else "silent_bias",
                n_O=len(O), est=est, abs_bias=float(abs(est - tau)),
                rel_bias=float(abs(est - tau) / max(abs(tau), 1e-12)))


def analyse_one(args):
    seed, p, deg, ens = args
    cfg = ENSEMBLES[ens]
    # ---- run_x1.py stream, byte-identical through K (gate G2) -----------------
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
    if len(K_all) < MAX_K:
        return None
    K = [(int(a), int(b)) for (a, b) in K_all[:MAX_K]]
    # ---- stream ends ----------------------------------------------------------

    base = dict(seed=int(seed), p=int(p), deg=float(deg), x=int(x), y=int(y),
                n_undirected=len(U), n_edges=int(D.sum()), tau=float(tau))

    # ---- eligibility (PREREG amendment 7): arm SK is PRIMARY -------------------
    # Identical to run_x1.py: G0 = Meek(C, K) must exist, O*(G0) defined, valid in D.
    try:
        G0_SK = apply_background_knowledge(C, K)
    except MeekFail:
        base["status"] = "drop_G0_meekfail"
        return base
    try:
        O0_SK = optimal_adjustment_set(G0_SK, x, y)
    except PathBlowup:
        base["status"] = "drop_blowup"
        return base
    if O0_SK is None:
        base["status"] = "drop_nonamenable"
        return base
    try:
        if not is_valid_adjustment_set(D, x, y, O0_SK):
            base["status"] = "drop_baseline_invalid"
            return base
    except PathBlowup:
        base["status"] = "drop_blowup"
        return base

    # secondary arm S0: only where the CPDAG alone is amenable
    cpdag_amenable = bool(is_amenable(C, x, y))
    O0_S0 = None
    if cpdag_amenable:
        try:
            O0_S0 = optimal_adjustment_set(C, x, y)
            if O0_S0 is not None and not is_valid_adjustment_set(D, x, y, O0_S0):
                O0_S0 = None
        except PathBlowup:
            O0_S0 = None

    distC = hop_dist_from(skeleton(C), [x, y], p)
    A = scm["A"]
    edges = [(int(i), int(j)) for i in range(p) for j in range(p) if D[i, j] == 1]

    perts = []
    for (i, j) in edges:
        Dp = D.copy()
        Dp[i, j] = 0
        try:
            Chat = dag_to_cpdag(Dp)
        except Exception:
            continue
        rec = dict(edge=[i, j], hop=int(stmt_hop(distC, i, j)),
                   w=float(abs(A[i, j])))
        # ---- arm SK (PRIMARY): the discovery error inside our setting --------
        Kr = [(a, b) for (a, b) in K if Chat[a, b] or Chat[b, a]]   # edge survives
        Kr = [(a, b) for (a, b) in Kr if is_undirected(Chat, a, b)] # still undirected
        rec["n_K_kept"] = len(Kr)
        try:
            Ghat = apply_background_knowledge(Chat, Kr) if Kr else Chat
            rec["SK"] = _outcome(Ghat, D, Sigma, x, y, O0_SK, tau)
        except MeekFail:
            rec["SK"] = dict(outcome="meekfail")
        except PathBlowup:
            rec["SK"] = dict(outcome="blowup")
        # ---- arm S0 (SECONDARY): no background knowledge ---------------------
        rec["S0"] = (_outcome(Chat, D, Sigma, x, y, O0_S0, tau)
                     if O0_S0 is not None else dict(outcome="baseline_unusable"))
        perts.append(rec)

    base.update(status="ok", n_O0_SK=len(O0_SK),
                cpdag_amenable=cpdag_amenable,
                n_O0_S0=(len(O0_S0) if O0_S0 is not None else None),
                K=[list(e) for e in K],
                hop_hist={str(h): sum(1 for r in perts if r["hop"] == h)
                          for h in sorted({r["hop"] for r in perts})},
                perts=perts)
    return base


def main():
    target = int(sys.argv[1])
    ens = sys.argv[2]
    out = sys.argv[3]
    nw = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    cfg = ENSEMBLES[ens]
    ps = list(cfg["ps"])
    rng = np.random.default_rng(20260831)
    n_draw = target * OVERSAMPLE[ens]
    jobs = [(int(s), int(rng.choice(ps)), float(rng.choice(cfg["degs"])), ens)
            for s in range(n_draw)]
    t0 = time.time()
    with Pool(nw) as pool:
        res = [r for r in pool.imap_unordered(analyse_one, jobs, chunksize=16)
               if r is not None]
    ok = [r for r in res if r.get("status") == "ok"]
    print(f"[{ens}] drawn={n_draw} returned={len(res)} analysed={len(ok)} "
          f"perturbations={sum(len(r['perts']) for r in ok)} "
          f"({time.time()-t0:.1f}s)", flush=True)
    json.dump(dict(ensemble=ens, target=target, n_drawn=n_draw,
                   prereg="../PREREG.md", scms=res), open(out, "w"))


if __name__ == "__main__":
    main()
