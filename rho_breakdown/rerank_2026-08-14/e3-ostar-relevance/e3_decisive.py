"""
E3 step 3 - THE DECISIVE TEST for the honest risk.

The pre-registration's condition (3) compared d_ostar_graded (a 0..6 graded score)
against d_xy_inc (a BINARY feature). That comparison is granularity-handicapped:
any graded score beats a binary one on AUC whenever there is signal, because ties
inside the two bins are scored 0.5. This script builds MATCHED-GRANULARITY
comparators that carry ONLY hop-distance information - no cn/pa/forb, no O* -
and asks whether the estimand-induced region buys anything over them.

  n_xy_inc       = sum_k 1[dist(k)=0]                        (count, 0..3)
  d_dist_graded  = sum_k (2 if dist=0 else 1 if dist=1 else 0) (0..6, same range
                   and same functional form as d_ostar_graded)

Plus: nested logistic test, the rho>=2 counterexample hunt, and the both-arms
distance table that reconciles with the campaign's "0/791 statements".

Usage: python e3_decisive.py
"""
import json
import os
import sys
from collections import defaultdict
from itertools import combinations

import numpy as np
from scipy.stats import chi2, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "pilotcode"))
ROWS = os.path.join(HERE, "results", "e3_rows.json")
OUT = os.path.join(HERE, "results", "e3_decisive.json")

NBOOT = 1000
RNG = np.random.default_rng(20260814)


def auc(score, label):
    score = np.asarray(score, float)
    label = np.asarray(label, bool)
    n1, n0 = int(label.sum()), int((~label).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    order = np.argsort(score, kind="mergesort")
    s = score[order]
    ranks = np.empty(len(s), float)
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        ranks[i:j + 1] = (i + j) / 2.0 + 1.0
        i = j + 1
    r = np.empty(len(s), float)
    r[order] = ranks
    return float((r[label].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def logit_fit(X, y, iters=200):
    """Plain Newton-IRLS with tiny ridge for stability. Returns beta, loglik."""
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    n, k = X.shape
    b = np.zeros(k)
    for _ in range(iters):
        eta = X @ b
        pmu = 1.0 / (1.0 + np.exp(-np.clip(eta, -30, 30)))
        W = np.clip(pmu * (1 - pmu), 1e-9, None)
        g = X.T @ (y - pmu) - 1e-6 * b
        H = X.T @ (X * W[:, None]) + 1e-6 * np.eye(k)
        step = np.linalg.solve(H, g)
        b = b + step
        if np.max(np.abs(step)) < 1e-9:
            break
    eta = np.clip(X @ b, -30, 30)
    ll = float(np.sum(y * eta - np.log1p(np.exp(eta))))
    return b, ll


def main():
    with open(ROWS) as f:
        blob = json.load(f)
    rows = blob["rows"]

    gen = [r for r in rows if r["arm"] == "generic"]
    P3 = [r for r in gen if r["consistent"] and r["amenable"]]
    y = np.array([r["silent"] for r in P3], bool)
    seeds = np.array([r["seed"] for r in P3], int)

    # ---- matched-granularity distance comparators (need per-statement distance,
    # which we recover from the stored min over the flip-set only for d_gdist; so
    # rebuild the graded distance score from the per-statement table)
    # d_ostar_graded and n_xy_inc are already per-row. We add d_dist_graded by
    # re-deriving per-statement distances: they were folded into d_gdist (min) and
    # n_xy_inc (count of dist==0). To get the exact graded score we recompute.
    with open(os.path.join(HERE, "results", "e3_stmt.json")) as f:
        stmt = json.load(f)
    key = {(s["seed"], s["arm"]): s["dist"] for s in stmt}

    d_dist_graded, n_dist1 = [], []
    for r in P3:
        ds = key[(r["seed"], r["arm"])]
        F = r["flip"]
        d_dist_graded.append(sum(2 if ds[k] == 0 else 1 if ds[k] == 1 else 0 for k in F))
        n_dist1.append(sum(1 for k in F if ds[k] <= 1))
    d_dist_graded = np.array(d_dist_graded, float)
    n_dist1 = np.array(n_dist1, float)

    feats = {
        "d_edit": np.array([r["d_edit"] for r in P3], float),
        "d_xy_inc": np.array([r["d_xy_inc"] for r in P3], float),
        "n_xy_inc": np.array([r["n_xy_inc"] for r in P3], float),
        "d_gdist": -np.array([r["d_gdist"] for r in P3], float),
        "d_dist_graded": d_dist_graded,
        "n_dist_le1": n_dist1,
        "d_ostar": np.array([r["d_ostar"] for r in P3], float),
        "d_ostar_graded": np.array([r["d_ostar_graded"] for r in P3], float),
    }

    uniq = np.unique(seeds)
    s2i = {s: np.flatnonzero(seeds == s) for s in uniq}

    point = {f: auc(v, y) for f, v in feats.items()}
    boots = defaultdict(list)
    PAIRS = [("d_ostar_graded", "n_xy_inc"),
             ("d_ostar_graded", "d_dist_graded"),
             ("d_ostar_graded", "d_edit"),
             ("d_ostar_graded", "d_xy_inc"),
             ("d_ostar", "n_xy_inc")]
    for _ in range(NBOOT):
        pick = RNG.choice(len(uniq), size=len(uniq), replace=True)
        bi = np.concatenate([s2i[uniq[k]] for k in pick])
        yb = y[bi]
        if yb.sum() == 0 or (~yb).sum() == 0:
            continue
        cur = {f: auc(v[bi], yb) for f, v in feats.items()}
        for f, a in cur.items():
            boots[f].append(a)
        for a, b in PAIRS:
            boots[f"D_{a}_vs_{b}"].append(cur[a] - cur[b])
    ci = {k: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
          for k, v in boots.items()}

    out = {"n": len(P3), "n_silent": int(y.sum()), "n_scm": int(len(uniq)),
           "auc": point, "auc_ci": {f: ci[f] for f in feats},
           "deltas": {f"D_{a}_vs_{b}": dict(point=point[a] - point[b],
                                            ci=ci[f"D_{a}_vs_{b}"])
                      for a, b in PAIRS}}

    # ---- nested logistic: does O*-region info add over distance info?
    n = len(P3)
    one = np.ones((n, 1))
    Xd = np.column_stack([one, feats["d_dist_graded"], feats["d_edit"]])
    Xo = np.column_stack([one, feats["d_dist_graded"], feats["d_edit"],
                          feats["d_ostar_graded"]])
    _, ll_d = logit_fit(Xd, y.astype(float))
    bo, ll_o = logit_fit(Xo, y.astype(float))
    lr = 2 * (ll_o - ll_d)
    out["nested_logistic"] = dict(
        ll_distance_only=ll_d, ll_plus_ostar=ll_o, LR_stat=float(lr), df=1,
        p=float(chi2.sf(lr, 1)), beta_ostar=float(bo[-1]),
        note=("naive LR test: members are clustered in SCMs so the nominal p is "
              "anti-conservative; read the sign and magnitude, not the p"))

    # ---- Spearman on |bias|/|tau| with matched comparators
    m = np.array([r["abs_rel_bias"] is not None for r in P3])
    arb = np.array([r["abs_rel_bias"] if r["abs_rel_bias"] is not None else np.nan
                    for r in P3])
    sp = {f: float(spearmanr(v[m], arb[m])[0]) for f, v in feats.items()}
    out["spearman_abs_rel_bias"] = dict(n=int(m.sum()), point=sp)

    # ---- the rho>=2 counterexample hunt: silent bias with NO statement at dist 0
    ce = []
    for r in P3:
        if r["silent"] and r["d_gdist"] >= 1:
            ce.append(dict(seed=r["seed"], arm=r["arm"], p=r["p"], deg=r["deg"],
                           rho=r["rho"], flip=r["flip"], d_gdist=r["d_gdist"],
                           d_ostar=r["d_ostar"], d_ostar_graded=r["d_ostar_graded"],
                           tau=r["tau"], bias=r["bias"],
                           abs_rel_bias=r["abs_rel_bias"]))
    out["counterexamples_silent_at_gdist_ge1"] = ce

    # ---- both-arms per-statement distance table at rho=1 (reconcile with 0/791)
    tab = defaultdict(lambda: [0, 0])
    for r in rows:
        if r["rho"] != 1:
            continue
        d = min(r["d_gdist"], 4)
        tab[d][0] += 1
        tab[d][1] += int(r["silent"])
    out["rho1_both_arms_by_distance"] = {
        f"dist={d}": dict(n=v[0], n_silent=v[1], rate=v[1] / v[0])
        for d, v in sorted(tab.items())}

    # ---- silent rate by min-distance, all rho, both arms
    tab2 = defaultdict(lambda: [0, 0])
    for r in rows:
        if not (r["consistent"] and r["amenable"]):
            continue
        d = min(r["d_gdist"], 4)
        tab2[d][0] += 1
        tab2[d][1] += int(r["silent"])
    out["P3_both_arms_by_min_distance"] = {
        f"min_dist={d}": dict(n=v[0], n_silent=v[1], rate=v[1] / v[0])
        for d, v in sorted(tab2.items())}

    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)

    print("=== AUC on generic P3 (n=%d, silent=%d, SCMs=%d) ===" % (
        out["n"], out["n_silent"], out["n_scm"]))
    for f, a in sorted(point.items(), key=lambda kv: -kv[1]):
        print("  %-18s %.4f  [%.4f,%.4f]" % (f, a, *ci[f]))
    print("\n=== pre-registered + matched-granularity deltas ===")
    for k, v in out["deltas"].items():
        print("  %-42s %+.4f  [%+.4f,%+.4f]" % (k, v["point"], *v["ci"]))
    print("\n=== nested logistic (distance-only -> +O*) ===")
    print("  LR=%.2f  p=%.3g  beta_ostar=%+.4f" % (
        out["nested_logistic"]["LR_stat"], out["nested_logistic"]["p"],
        out["nested_logistic"]["beta_ostar"]))
    print("\n=== counterexamples: silent bias with min graph-distance >= 1 ===")
    print(json.dumps(ce, indent=1))
    print("\n=== rho=1, BOTH arms, by statement distance ===")
    print(json.dumps(out["rho1_both_arms_by_distance"], indent=1))
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
