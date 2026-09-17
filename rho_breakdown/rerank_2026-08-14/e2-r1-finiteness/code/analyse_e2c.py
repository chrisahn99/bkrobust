"""
E2 attack E (POST-HOC) -- the decisive one.

A tier move changes K in three ways at once: it can REVERSE a statement (an edge
cross-tier both before and after, but in the other direction), REMOVE one (the
edge becomes within-tier), or ADD one (a within-tier edge becomes cross-tier).

If the T-move ball is dominated by REMOVALS -- i.e. the "wrong expert" mostly
just says LESS rather than saying something FALSE -- then containment is an
artefact of a weaker perturbation and the SUPPORTED verdict is worthless. The
sharpest form of the question: how often does a tier move actually assert
something that CONTRADICTS the true DAG?

This re-derives D* per SCM (bitwise, from the stored (seed,p,deg)) so that each
statement of K' can be checked against ground truth.
"""
import json
import sys
from collections import Counter

import numpy as np
from scipy.stats import mannwhitneyu

from adjust import possibly_causal_paths
from graphs import dag_to_cpdag, random_dag, undirected_edges
from run_e2 import K_from_tiering

ARMS = ("G-flip", "T-flip", "T-move")


def desc(v):
    v = np.asarray(v, dtype=float)
    return (len(v), v.mean(), np.median(v), np.quantile(v, .90), v.max()) if len(v) else None


def main():
    data = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "../results/e2_raw.json"))
    matched = [r for r in data if all(a in r["arms"] for a in ARMS)]
    print("=" * 78)
    print("ATTACK E (POST-HOC) -- is a tier move a REAL error, or just LESS knowledge?")
    print("=" * 78)

    rows = []
    for r in matched:
        D = random_dag(r["p"], r["deg"], np.random.default_rng(r["seed"]))
        assert int(D.sum()) == r["n_edges"], "DAG replay mismatch"
        C = dag_to_cpdag(D)
        tier = np.array(r["tier"])
        K = [tuple(e) for e in r["arms"]["T-move"]["K"]]
        Kset = {frozenset(e): e for e in K}
        for m in r["arms"]["T-move"]["members"]:
            tp = tier.copy()
            for (v, t) in m["moves"]:
                tp[v] = t
            Kp = K_from_tiering(tp, C)
            Kpset = {frozenset(e): e for e in Kp}
            n_rev = sum(1 for k in Kset if k in Kpset and Kset[k] != Kpset[k])
            n_rem = sum(1 for k in Kset if k not in Kpset)
            n_add = sum(1 for k in Kpset if k not in Kset)
            n_false = sum(1 for e in Kp if D[e[0], e[1]] != 1)   # contradicts D*
            rows.append(dict(seed=r["seed"], rho=m["rho"], n_changed=m["n_changed"],
                             n_rev=n_rev, n_rem=n_rem, n_add=n_add, n_false=n_false,
                             consistent=m.get("consistent"),
                             r1eq=m.get("r1_equals_full"),
                             d_orient=m.get("d_orient"),
                             ostar_changed=m.get("ostar_changed"),
                             ostar_valid=m.get("ostar_valid")))

    R = [x for x in rows if x["n_changed"] > 0]
    print(f"\nT-move perturbations with n_changed>0: {len(R)}")
    print("\nCOMPOSITION OF A TIER MOVE (counts of statements per member):")
    for k in ("n_rev", "n_rem", "n_add", "n_false"):
        v = np.array([x[k] for x in R], dtype=float)
        print(f"  {k:<8} mean {v.mean():.3f}  median {np.median(v):.0f}  "
              f"max {v.max():.0f}  frac>0 {np.mean(v>0):.3f}")

    n_false_pos = sum(1 for x in R if x["n_false"] > 0)
    print(f"\n>>> {n_false_pos}/{len(R)} = {n_false_pos/len(R):.3f} of tier moves assert at "
          f"least one statement CONTRADICTING the true DAG.")
    print("    (If this were near 0 the arm would be measuring nothing.)")

    print("\nSPLIT: members that assert a FALSE statement vs those that only ADD/REMOVE true ones")
    print(f"{'subset':<26}{'n':>7}{'caught':>9}{'silent':>9}"
          f"{'d_orient med':>14}{'d_orient q90':>14}{'R1-suff':>9}")
    for lab, pred in (("asserts >=1 FALSE stmt", lambda x: x["n_false"] > 0),
                      ("only true add/remove", lambda x: x["n_false"] == 0)):
        S = [x for x in R if pred(x)]
        cons = [x for x in S if x["consistent"]]
        caught = sum(1 for x in S if not x["consistent"])
        silent = sum(1 for x in S if x["ostar_changed"] and not x["ostar_valid"])
        do = np.array([x["d_orient"] for x in cons], dtype=float)
        r1 = np.mean([x["r1eq"] for x in cons])
        print(f"{lab:<26}{len(S):>7}{caught/len(S):>9.3f}{silent/len(S):>9.3f}"
              f"{np.median(do):>14.2f}{np.quantile(do,.9):>14.2f}{r1:>9.4f}")

    # the containment claim restricted to REAL errors only
    print("\n" + "-" * 78)
    print("GATE 2 RE-RUN on the FALSE-ASSERTING SUBSET ONLY (the honest comparison:")
    print("both arms now definitely misstate the graph)")
    print("-" * 78)
    rt = np.array([x["d_orient"] / x["n_changed"] for x in R
                   if x["consistent"] and x["n_false"] > 0], dtype=float)
    raw = json.load(open("../results/e2_raw.json"))
    matchedmap = {r["seed"]: r for r in matched}
    rg, rf = [], []
    for r in matched:
        for a, acc in (("G-flip", rg), ("T-flip", rf)):
            for m in r["arms"][a]["members"]:
                if m.get("consistent") and m["n_changed"] > 0:
                    acc.append(m["d_orient"] / m["n_changed"])
    rg, rf = np.array(rg), np.array(rf)
    print(f"{'arm':<28}{'n':>7}{'mean':>9}{'median':>9}{'q90':>9}{'max':>9}")
    for lab, v in (("T-move (FALSE-asserting)", rt), ("T-flip (all, K=same)", rf),
                   ("G-flip (all)", rg)):
        d = desc(v)
        print(f"{lab:<28}{d[0]:>7}{d[1]:>9.3f}{d[2]:>9.3f}{d[3]:>9.3f}{d[4]:>9.3f}")
    print(f"  MWU T-move(false) < T-flip : p = "
          f"{mannwhitneyu(rt, rf, alternative='less').pvalue:.3e}")
    print(f"  MWU T-move(false) < G-flip : p = "
          f"{mannwhitneyu(rt, rg, alternative='less').pvalue:.3e}")
    print(f"  MWU T-move(false) > G-flip : p = "
          f"{mannwhitneyu(rt, rg, alternative='greater').pvalue:.3e}")

    eq = [x["r1eq"] for x in R if x["consistent"] and x["n_false"] > 0]
    print(f"\n  GATE 1 on the FALSE-asserting subset: {sum(eq)}/{len(eq)} = "
          f"{sum(eq)/len(eq):.4f}")


if __name__ == "__main__":
    main()
