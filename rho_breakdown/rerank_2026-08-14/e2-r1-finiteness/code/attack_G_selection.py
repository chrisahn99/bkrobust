"""
Second probe. The claim attributes containment to "the error model staying inside the
R1-closed sublattice". Its own null control (T-flip vs G-flip) says the FORM of K does
nothing. So the residual effect among matched operations must be either
  (a) R1-closure, or
  (b) WHICH statements a tier move happens to reverse (selection on topology).

Fully matched contrast: T-move restricted to PURE REVERSALS  vs  T-flip.
Same SCM, same K (identical tiered set), same operation (reverse rho statements),
same n_changed. The ONLY difference is which statements, and R1-closure.

Then: restrict T-flip to its R1-SUFFICIENT members. If R1-closure is the mechanism,
that restriction should close the gap. If it does not, R1-closure is not the mechanism.
"""
import gzip
import json
from collections import defaultdict

import numpy as np
from scipy.stats import mannwhitneyu, wilcoxon

from graphs import dag_to_cpdag, random_dag
from run_e2 import K_from_tiering

ARMS = ("G-flip", "T-flip", "T-move")


def desc(name, v):
    v = np.asarray(v, dtype=float)
    if not len(v):
        print(f"{name:<42}{'(empty)':>8}")
        return
    print(f"{name:<42}{len(v):>7}{v.mean():>9.3f}{np.median(v):>9.3f}"
          f"{np.quantile(v,.90):>8.3f}{v.max():>8.3f}")


def main():
    data = json.load(gzip.open("../results/e2_raw.json.gz"))
    matched = [r for r in data if all(a in r["arms"] for a in ARMS)]

    tmove, gflip, tflip = [], [], []
    for r in matched:
        D = random_dag(r["p"], r["deg"], np.random.default_rng(r["seed"]))
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
            rec = dict(m)
            rec.update(seed=r["seed"],
                       n_rev=sum(1 for k in Kset if k in Kpset and Kset[k] != Kpset[k]),
                       n_rem=sum(1 for k in Kset if k not in Kpset),
                       n_add=sum(1 for k in Kpset if k not in Kset),
                       n_false_new=sum(1 for e in Kp if D[e[0], e[1]] != 1
                                       and Kset.get(frozenset(e)) != e))
            tmove.append(rec)
        for arm, acc in (("G-flip", gflip), ("T-flip", tflip)):
            for m in r["arms"][arm]["members"]:
                rec = dict(m)
                rec.update(seed=r["seed"], n_rev=m["rho"], n_rem=0, n_add=0,
                           n_false_new=m["rho"])
                acc.append(rec)

    ok = lambda L: [x for x in L if x.get("consistent") and x["n_changed"] > 0]
    tm, gf, tf = ok(tmove), ok(gflip), ok(tflip)
    pure = [x for x in tm if x["n_add"] == 0 and x["n_rem"] == 0 and x["n_rev"] > 0]

    print("=" * 96)
    print("A. FULLY MATCHED: T-move PURE-REVERSAL vs T-flip  (same K, same operation)")
    print("=" * 96)
    print(f"{'stratum':<42}{'n':>7}{'mean':>9}{'median':>9}{'q90':>8}{'max':>8}")
    for k in (1, 2, 3):
        desc(f"  n_chg={k}  T-flip        d_orient",
             [x["d_orient"] for x in tf if x["n_changed"] == k])
        desc(f"  n_chg={k}  T-move PURE   d_orient",
             [x["d_orient"] for x in pure if x["n_changed"] == k])
    a = np.array([x["d_orient"] / x["n_changed"] for x in pure], float)
    b = np.array([x["d_orient"] / x["n_changed"] for x in tf], float)
    print(f"\n  pooled r:  T-move PURE {a.mean():.3f}   T-flip {b.mean():.3f}")
    print(f"  MWU H1: r(T-move PURE) < r(T-flip)  p = "
          f"{mannwhitneyu(a, b, alternative='less').pvalue:.3e}")

    print("\n" + "=" * 96)
    print("B. IS R1-CLOSURE THE MECHANISM? restrict T-flip to its R1-SUFFICIENT members")
    print("   (T-move is 100% R1-sufficient; if R1-closure causes containment, the")
    print("    R1-sufficient slice of T-flip should look like T-move)")
    print("=" * 96)
    tf_r1 = [x for x in tf if x.get("r1_equals_full")]
    tf_no = [x for x in tf if not x.get("r1_equals_full")]
    gf_r1 = [x for x in gf if x.get("r1_equals_full")]
    gf_no = [x for x in gf if not x.get("r1_equals_full")]
    print(f"{'slice':<42}{'n':>7}{'mean':>9}{'median':>9}{'q90':>8}{'max':>8}")
    for nm, S in (("T-flip  R1-sufficient   r", tf_r1),
                  ("T-flip  NOT R1-suff     r", tf_no),
                  ("G-flip  R1-sufficient   r", gf_r1),
                  ("G-flip  NOT R1-suff     r", gf_no),
                  ("T-move  (all R1-suff)   r", tm),
                  ("T-move PURE (R1-suff)   r", pure)):
        desc("  " + nm, [x["d_orient"] / x["n_changed"] for x in S])
    x1 = np.array([x["d_orient"] / x["n_changed"] for x in tf_r1], float)
    x0 = np.array([x["d_orient"] / x["n_changed"] for x in tf_no], float)
    print(f"\n  Within T-flip, does R1-sufficiency predict LESS change?")
    print(f"    MWU H1: r(R1-suff) < r(NOT R1-suff)  p = "
          f"{mannwhitneyu(x1, x0, alternative='less').pvalue:.3e}")
    print(f"    MWU H1: r(R1-suff) > r(NOT R1-suff)  p = "
          f"{mannwhitneyu(x1, x0, alternative='greater').pvalue:.3e}")
    y1 = np.array([x["d_orient"] / x["n_changed"] for x in gf_r1], float)
    y0 = np.array([x["d_orient"] / x["n_changed"] for x in gf_no], float)
    print(f"    G-flip same test: MWU H1 r(R1-suff) < r(NOT R1-suff) p = "
          f"{mannwhitneyu(y1, y0, alternative='less').pvalue:.3e}")

    print("\n" + "=" * 96)
    print("C. DAMAGE per FALSE statement asserted (the quantity a reader cares about)")
    print("=" * 96)
    print(f"{'arm':<42}{'n':>7}{'caught':>9}{'silent':>9}{'mean nfalse':>13}"
          f"{'silent/false':>13}")
    for nm, S_all in (("G-flip", gflip), ("T-flip", tflip), ("T-move", tmove)):
        S = [x for x in S_all if x["n_changed"] > 0 and x["n_false_new"] > 0]
        caught = np.mean([not x.get("consistent") for x in S])
        silent = np.mean([bool(x.get("ostar_changed")) and not x.get("ostar_valid")
                          for x in S])
        nf = np.mean([x["n_false_new"] for x in S])
        print(f"  {nm:<40}{len(S):>7}{caught:>9.3f}{silent:>9.3f}{nf:>13.3f}"
              f"{silent/nf:>13.4f}")

    print("\n" + "=" * 96)
    print("D. GATE-2 STATISTIC AT n_changed = 1 ONLY (denominator identical = 1,")
    print("   so r == d_orient and NO normalisation choice can be doing the work)")
    print("=" * 96)
    print(f"{'arm (n_changed==1)':<42}{'n':>7}{'mean':>9}{'median':>9}{'q90':>8}{'max':>8}")
    sub = {}
    for nm, S in (("G-flip", gf), ("T-flip", tf), ("T-move ALL", tm),
                  ("T-move: 1 reversal only", [x for x in pure if x["n_rev"] == 1]),
                  ("T-move: 1 removal only",
                   [x for x in tm if x["n_rem"] == 1 and x["n_add"] == 0 and x["n_rev"] == 0]),
                  ("T-move: 1 addition only",
                   [x for x in tm if x["n_add"] == 1 and x["n_rem"] == 0 and x["n_rev"] == 0])):
        S1 = [x for x in S if x["n_changed"] == 1]
        sub[nm] = S1
        desc("  " + nm, [x["d_orient"] for x in S1])
    print("\n  --> the n_changed==1 stratum decomposed by OPERATION shows whether the")
    print("      'containment' is about tiering or about removals/additions being cheap.")

    print("\n" + "=" * 96)
    print("E. PER-SCM PAIRED (355 clusters): T-move PURE-REVERSAL vs T-flip")
    print("=" * 96)
    gm, tmm = defaultdict(list), defaultdict(list)
    for x in tf:
        gm[x["seed"]].append(x["d_orient"] / x["n_changed"])
    for x in pure:
        tmm[x["seed"]].append(x["d_orient"] / x["n_changed"])
    seeds = sorted(set(gm) & set(tmm))
    A = np.array([np.mean(gm[s]) for s in seeds])
    B = np.array([np.mean(tmm[s]) for s in seeds])
    print(f"  SCMs with both: {len(seeds)} of 355   mean(T-movePURE - T-flip) = "
          f"{(B-A).mean():+.4f}  median {np.median(B-A):+.4f}")
    print(f"  Wilcoxon H1: T-movePURE < T-flip  p = "
          f"{wilcoxon(B, A, alternative='less').pvalue:.3e}")


if __name__ == "__main__":
    main()
