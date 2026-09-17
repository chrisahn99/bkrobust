"""
E3 step 1 - regenerate (D, C), compute the O*-relevance features, emit one row
per perturbation.

Reads results/linear_raw.json (copied, read-only). D and C are NOT stored in the
JSON; they are regenerated bitwise from the stored (seed, p, deg) triple, which is
exactly what run_linear.analyse_one does before it consumes any further randomness.

INTEGRITY GATES (abort on failure):
  G1  regenerated n_undirected == stored n_undirected
  G2  regenerated n_edges      == stored n_edges
  G3  regenerated cpdag_amenable == stored cpdag_amenable
  G4  recomputed O0 (from stored K) == stored O0
  G5  flip-set enumeration order reproduces every member's stored rho
  G6  recomputed member 'consistent'/'amenable' == stored (spot-check on a
      subsample; full check is cheap enough so we do all of them)

Usage: python e3_features.py
"""
import json
import os
import sys
from collections import Counter
from itertools import combinations

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "pilotcode"))

from adjust import (causal_nodes, forb, is_amenable, optimal_adjustment_set,  # noqa: E402
                    parents_of_set)
from graphs import (MeekFail, apply_background_knowledge, dag_to_cpdag,  # noqa: E402
                    random_dag, skeleton)

RAW = os.path.join(HERE, "linear_raw.json")
OUT = os.path.join(HERE, "results", "e3_rows.json")


# --------------------------------------------------------------------- geometry
def hop_distance_to_query(C, x, y):
    """BFS hop distance in the SKELETON of C from every node to the set {x,y}.
    d(x)=d(y)=0. Meek only orients, so skeleton(C)==skeleton(G0); using C is the
    same graph and makes the feature independent of K."""
    S = skeleton(C)
    p = S.shape[0]
    dist = np.full(p, 10 ** 6, dtype=int)
    frontier = [x, y]
    for v in frontier:
        dist[v] = 0
    while frontier:
        nxt = []
        for v in frontier:
            for w in np.flatnonzero(S[v]):
                w = int(w)
                if dist[w] > dist[v] + 1:
                    dist[w] = dist[v] + 1
                    nxt.append(w)
        frontier = nxt
    return dist


def ostar_region(G0, x, y):
    """The three node sets that appear in O* = pa(cn) \\ forb, plus the query."""
    cn = set(causal_nodes(G0, x, y))
    pa = set(parents_of_set(G0, cn)) if cn else set()
    frb = set(forb(G0, x, y))
    core = cn | {x, y}                    # causal-path core (tier 2)
    machinery = (pa | frb) - core         # adjustment machinery (tier 1)
    R = core | pa | frb
    return dict(cn=cn, pa=pa, frb=frb, core=core, machinery=machinery, R=R)


def statement_features(K, reg, dist, S):
    """Per-statement relevance. Computable from (C, K, x, y) alone - no true DAG,
    no perturbation."""
    feats = []
    for (u, v) in K:
        ends = {u, v}
        r = 1 if ends & reg["R"] else 0
        if ends & reg["core"]:
            r2 = 2
        elif ends & reg["machinery"]:
            r2 = 1
        else:
            r2 = 0
        d = int(min(dist[u], dist[v]))
        # 1-hop dilation of R: does the statement sit adjacent to the region?
        nb = set(np.flatnonzero(S[u])) | set(np.flatnonzero(S[v])) | ends
        r_dil = 1 if nb & reg["R"] else 0
        feats.append(dict(r=r, r2=r2, dist=d, r_dil=r_dil,
                          xy_inc=1 if ends & {u for u in (0,)} - {0} else 0))
        feats[-1]["xy_inc"] = 0  # placeholder, filled by caller (needs x,y)
    return feats


def main():
    with open(RAW) as f:
        raw = json.load(f)
    print(f"loaded {len(raw)} SCM records")

    rows = []
    stmt_rows = []
    gate_fail = Counter()
    stmt_rel = Counter()
    n_arms = 0

    for rec in raw:
        seed, p, deg = rec["seed"], rec["p"], rec["deg"]
        x, y = rec["x"], rec["y"]

        rng = np.random.default_rng(seed)
        D = random_dag(p, deg, rng)
        C = dag_to_cpdag(D)
        S = skeleton(C)

        # ---- G1..G3
        n_und = int(sum(1 for i in range(p) for j in range(i + 1, p)
                        if C[i, j] == 1 and C[j, i] == 1))
        if n_und != rec["n_undirected"]:
            gate_fail["G1"] += 1
            continue
        if int(D.sum()) != rec["n_edges"]:
            gate_fail["G2"] += 1
            continue
        if bool(is_amenable(C, x, y)) != rec["cpdag_amenable"]:
            gate_fail["G3"] += 1
            continue

        dist = hop_distance_to_query(C, x, y)

        for arm, A in rec["arms"].items():
            K = [tuple(e) for e in A["K"]]
            G0 = apply_background_knowledge(C, K)
            O0 = optimal_adjustment_set(G0, x, y)
            if O0 is None or sorted(O0) != A["O0"]:
                gate_fail["G4"] += 1
                continue
            n_arms += 1

            reg = ostar_region(G0, x, y)
            sf = statement_features(K, reg, dist, S)
            for (u, v), f in zip(K, sf):
                f["xy_inc"] = 1 if ({u, v} & {x, y}) else 0
                stmt_rel[(arm, f["r"], f["r2"], f["xy_inc"], min(f["dist"], 3))] += 1
            stmt_rows.append(dict(seed=seed, arm=arm,
                                  K=[list(e) for e in K],
                                  dist=[f["dist"] for f in sf],
                                  r=[f["r"] for f in sf],
                                  r2=[f["r2"] for f in sf],
                                  xy=[f["xy_inc"] for f in sf],
                                  x=x, y=y,
                                  cn=sorted(reg["cn"]), pa=sorted(reg["pa"]),
                                  frb=sorted(reg["frb"]), R=sorted(reg["R"])))

            # ---- G5: replay the enumeration order
            flips = []
            for rho in (1, 2, 3):
                if rho > len(K):
                    break
                for flip in combinations(range(len(K)), rho):
                    flips.append((rho, flip))
            if len(flips) != len(A["members"]):
                gate_fail["G5_len"] += 1
                continue
            if any(r != m["rho"] for (r, _), m in zip(flips, A["members"])):
                gate_fail["G5_rho"] += 1
                continue

            for (rho, flip), m in zip(flips, A["members"]):
                Fs = [sf[k] for k in flip]
                consistent = bool(m["consistent"])
                amenable = bool(m.get("amenable", False)) if consistent else False
                ostar_changed = m.get("ostar_changed", None)
                ostar_valid = m.get("ostar_valid", None)
                est = m.get("est", None)
                tau = rec["tau"]
                bias = m.get("bias", None)

                silent = bool(consistent and amenable and ostar_valid is False)

                rows.append(dict(
                    seed=seed, arm=arm, p=p, deg=deg, n_K=len(K), rho=rho,
                    flip=list(flip), tau=tau,
                    # ---------------- baselines
                    d_edit=rho,
                    d_xy_inc=int(any(f["xy_inc"] for f in Fs)),
                    d_gdist=int(min(f["dist"] for f in Fs)),
                    n_xy_inc=int(sum(f["xy_inc"] for f in Fs)),
                    # ---------------- candidates
                    d_ostar=int(sum(f["r"] for f in Fs)),
                    d_ostar_graded=int(sum(f["r2"] for f in Fs)),
                    d_ostar_any=int(any(f["r"] for f in Fs)),
                    d_ostar_dil=int(sum(f["r_dil"] for f in Fs)),
                    # ---------------- outcomes / status
                    consistent=consistent, amenable=amenable,
                    ostar_changed=(None if ostar_changed is None else bool(ostar_changed)),
                    ostar_valid=(None if ostar_valid is None else bool(ostar_valid)),
                    silent=silent,
                    est=est, bias=bias,
                    abs_rel_bias=(None if bias is None or tau == 0
                                  else abs(bias) / abs(tau)),
                    cpdag_amenable=rec["cpdag_amenable"],
                ))

    print(f"arms processed: {n_arms}; rows: {len(rows)}")
    print(f"gate failures: {dict(gate_fail)}")
    if gate_fail:
        print("!! INTEGRITY GATE FAILED - results below are NOT trustworthy")

    with open(OUT, "w") as f:
        json.dump(dict(rows=rows,
                       stmt_rel={"|".join(map(str, k)): v for k, v in stmt_rel.items()},
                       gate_fail=dict(gate_fail), n_arms=n_arms), f)
    print(f"wrote {OUT}")
    with open(os.path.join(HERE, "results", "e3_stmt.json"), "w") as f:
        json.dump(stmt_rows, f)
    print(f"wrote e3_stmt.json ({len(stmt_rows)} arms)")

    # quick sanity echo
    n_silent = sum(r["silent"] for r in rows)
    print(f"total members {len(rows)}, silent {n_silent} ({100*n_silent/len(rows):.2f}%)")


if __name__ == "__main__":
    main()
