"""
REFUTATION PROBE: is Gate 2's containment an artefact of the DENOMINATOR?

n_changed = |K delta K'| counts a TRUE removal and a TRUE addition as one unit of
"perturbation", identically to a sign reversal. But by Meek monotonicity:
  - a removal can only SHRINK the oriented set  -> cannot inject a false orientation
  - adding a TRUE statement can only refine toward D* -> cannot inject a false one
  - only a REVERSAL / a FALSE addition injects error.

G-flip perturbations are 100% reversals => every unit of n_changed is error-injecting.
T-move perturbations are a MIX => harmless units pad the denominator.

So r = d_orient / n_changed is not commensurable across arms. The honest normaliser
is the number of FALSE statements asserted, which is exactly rho for G-flip.
"""
import gzip
import json

import numpy as np
from scipy.stats import mannwhitneyu

from graphs import dag_to_cpdag, random_dag
from run_e2 import K_from_tiering

ARMS = ("G-flip", "T-flip", "T-move")
RAW = "../results/e2_raw.json.gz"


def q(v, p):
    return float(np.quantile(np.asarray(v, dtype=float), p)) if len(v) else float("nan")


def desc(name, v):
    v = np.asarray(v, dtype=float)
    if not len(v):
        print(f"{name:<34}{'-':>8}")
        return
    print(f"{name:<34}{len(v):>8}{v.mean():>9.3f}{np.median(v):>9.3f}"
          f"{q(v,.75):>8.3f}{q(v,.90):>8.3f}{q(v,.99):>8.3f}{v.max():>8.3f}")


def main():
    data = json.load(gzip.open(RAW))
    matched = [r for r in data if all(a in r["arms"] for a in ARMS)]
    print(f"matched SCMs: {len(matched)}")

    tmove, gflip, tflip = [], [], []
    for r in matched:
        D = random_dag(r["p"], r["deg"], np.random.default_rng(r["seed"]))
        assert int(D.sum()) == r["n_edges"]
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
            n_false = sum(1 for e in Kp if D[e[0], e[1]] != 1)
            n_false_new = sum(1 for e in Kp if D[e[0], e[1]] != 1
                              and Kset.get(frozenset(e)) != e)
            rec = dict(m)
            rec.update(seed=r["seed"], n_rev=n_rev, n_rem=n_rem, n_add=n_add,
                       n_false=n_false, n_false_new=n_false_new)
            tmove.append(rec)
        for arm, acc in (("G-flip", gflip), ("T-flip", tflip)):
            for m in r["arms"][arm]["members"]:
                rec = dict(m)
                # a flip reverses rho TRUE statements -> rho false statements asserted
                rec.update(seed=r["seed"], n_rev=m["rho"], n_rem=0, n_add=0,
                           n_false=m["rho"], n_false_new=m["rho"])
                acc.append(rec)

    C_ = lambda L: [x for x in L if x.get("consistent") and x["n_changed"] > 0]
    tm, gf, tf = C_(tmove), C_(gflip), C_(tflip)

    print("\n" + "=" * 96)
    print("0. REPRODUCE THE REPORTED GATE 2 (sanity: my pipeline == theirs)")
    print("=" * 96)
    print(f"{'arm':<34}{'n':>8}{'mean':>9}{'median':>9}{'q75':>8}{'q90':>8}{'q99':>8}{'max':>8}")
    for nm, S in (("G-flip r=do/n_changed", gf), ("T-flip r=do/n_changed", tf),
                  ("T-move r=do/n_changed", tm)):
        desc(nm, [x["d_orient"] / x["n_changed"] for x in S])

    print("\n" + "=" * 96)
    print("1. COMPOSITION: what fraction of each arm's n_changed is ERROR-INJECTING?")
    print("=" * 96)
    for nm, S in (("G-flip", gf), ("T-flip", tf), ("T-move", tm)):
        nc = np.array([x["n_changed"] for x in S], float)
        nr = np.array([x["n_rev"] for x in S], float)
        na = np.array([x["n_add"] for x in S], float)
        nm_ = np.array([x["n_rem"] for x in S], float)
        nf = np.array([x["n_false_new"] for x in S], float)
        print(f"  {nm:<8} mean n_changed {nc.mean():.3f} | n_rev {nr.mean():.3f} "
              f"n_add {na.mean():.3f} n_rem {nm_.mean():.3f} | NEW-FALSE stmts "
              f"{nf.mean():.3f} | error-injecting share of n_changed "
              f"{(nf.sum()/nc.sum()):.3f}")

    print("\n" + "=" * 96)
    print("2. THE HONEST NORMALISER: d_orient per NEWLY-ASSERTED FALSE STATEMENT")
    print("   (for G-flip this IS n_changed, so G-flip's number is unchanged)")
    print("=" * 96)
    print(f"{'arm':<34}{'n':>8}{'mean':>9}{'median':>9}{'q75':>8}{'q90':>8}{'q99':>8}{'max':>8}")
    rr = {}
    for nm, S in (("G-flip  do/n_false_new", gf), ("T-flip  do/n_false_new", tf),
                  ("T-move  do/n_false_new", tm)):
        v = [x["d_orient"] / x["n_false_new"] for x in S if x["n_false_new"] > 0]
        rr[nm.split()[0]] = np.array(v, float)
        desc(nm, v)
    print(f"\n  MWU H1: r_false(T-move) > r_false(G-flip)  p = "
          f"{mannwhitneyu(rr['T-move'], rr['G-flip'], alternative='greater').pvalue:.3e}")
    print(f"  MWU H1: r_false(T-move) < r_false(G-flip)  p = "
          f"{mannwhitneyu(rr['T-move'], rr['G-flip'], alternative='less').pvalue:.3e}")
    print("  PRE-REGISTERED GATE 2 REFUTED-RULE: median or q90 of T-move GREATER than "
          "G-flip AND MWU(greater) p<0.05")

    print("\n" + "=" * 96)
    print("3. CASCADE PER ERROR-INJECTING STATEMENT (rule propagation, the theorem's object)")
    print("=" * 96)
    print(f"{'arm':<34}{'n':>8}{'mean':>9}{'median':>9}{'q75':>8}{'q90':>8}{'q99':>8}{'max':>8}")
    cc = {}
    for nm, S in (("G-flip  cascade/n_false_new", gf), ("T-flip  cascade/n_false_new", tf),
                  ("T-move  cascade/n_false_new", tm)):
        v = [x["cascade"] / x["n_false_new"] for x in S if x["n_false_new"] > 0]
        cc[nm.split()[0]] = np.array(v, float)
        desc(nm, v)
    print(f"\n  MWU H1: cascade_rate(T-move) > (G-flip) p = "
          f"{mannwhitneyu(cc['T-move'], cc['G-flip'], alternative='greater').pvalue:.3e}")
    for nm, S in (("G-flip", gf), ("T-flip", tf), ("T-move", tm)):
        Sf = [x for x in S if x["n_false_new"] > 0]
        fr = np.mean([x["cascade"] > 0 for x in Sf])
        print(f"  {nm:<8} frac cascading | >=1 false stmt asserted: {fr:.3f}  (n={len(Sf)})")

    print("\n" + "=" * 96)
    print("4. PURE-REVERSAL TIER MOVES vs G-FLIP at MATCHED n_changed")
    print("   (T-move restricted to n_add=n_rem=0 -> same operation as a G-flip)")
    print("=" * 96)
    pure = [x for x in tm if x["n_add"] == 0 and x["n_rem"] == 0 and x["n_rev"] > 0]
    print(f"  pure-reversal T-move members: {len(pure)} of {len(tm)} "
          f"({len(pure)/len(tm):.4f})")
    print(f"{'stratum':<34}{'n':>8}{'mean':>9}{'median':>9}{'q75':>8}{'q90':>8}{'q99':>8}{'max':>8}")
    for k in (1, 2, 3):
        desc(f"  n_chg={k} G-flip  d_orient",
             [x["d_orient"] for x in gf if x["n_changed"] == k])
        desc(f"  n_chg={k} T-move PURE d_orient",
             [x["d_orient"] for x in pure if x["n_changed"] == k])
    a = np.array([x["d_orient"] / x["n_changed"] for x in pure], float)
    b = np.array([x["d_orient"] / x["n_changed"] for x in gf], float)
    if len(a):
        print(f"\n  MWU H1: r(T-move PURE) > r(G-flip) p = "
              f"{mannwhitneyu(a, b, alternative='greater').pvalue:.3e}")
        print(f"  MWU H1: r(T-move PURE) < r(G-flip) p = "
              f"{mannwhitneyu(a, b, alternative='less').pvalue:.3e}")

    print("\n" + "=" * 96)
    print("5. IS THE 'PROVABLY HARMLESS' SUBSET A TAUTOLOGY? (n_false==0 members)")
    print("=" * 96)
    for nm, S in (("T-move", tmove),):
        Z = [x for x in S if x["n_changed"] > 0 and x["n_false"] == 0]
        caught = sum(1 for x in Z if not x.get("consistent"))
        silent = sum(1 for x in Z if x.get("ostar_changed") and not x.get("ostar_valid"))
        print(f"  {nm}: n={len(Z)}  caught={caught}  silent={silent}   "
              f"(K' all-true => D* agrees with G' => Alg.1 cannot FAIL and O* must be "
              f"valid; 0/0 is a THEOREM, not a measurement)")
        share = len(Z) / len([x for x in S if x["n_changed"] > 0])
        print(f"  share of the T-move ball that is provably harmless: {share:.4f}")

    print("\n" + "=" * 96)
    print("6. PER-SCM PAIRED TEST (members are nested in 355 SCMs; MWU over 29404 vs 1528")
    print("   pooled members is pseudo-replicated). Paired on the SCM mean.")
    print("=" * 96)
    from collections import defaultdict
    from scipy.stats import wilcoxon
    for lab, key in (("r = d_orient/n_changed", lambda x: x["d_orient"] / x["n_changed"]),
                     ("r_false = d_orient/n_false_new",
                      lambda x: x["d_orient"] / x["n_false_new"])):
        gm, tmm = defaultdict(list), defaultdict(list)
        for x in gf:
            if key is not None and (x["n_false_new"] > 0 or "n_changed" in x):
                try:
                    gm[x["seed"]].append(key(x))
                except ZeroDivisionError:
                    pass
        for x in tm:
            try:
                tmm[x["seed"]].append(key(x))
            except ZeroDivisionError:
                pass
        seeds = sorted(set(gm) & set(tmm))
        A = np.array([np.mean(gm[s]) for s in seeds])
        B = np.array([np.mean(tmm[s]) for s in seeds])
        d = B - A
        print(f"  {lab:<32} SCMs={len(seeds)}  mean(T-move - G-flip) = {d.mean():+.4f}  "
              f"median = {np.median(d):+.4f}  frac SCMs where T-move HIGHER = "
              f"{np.mean(d>0):.3f}")
        try:
            print(f"      Wilcoxon signed-rank H1: T-move > G-flip  p = "
                  f"{wilcoxon(B, A, alternative='greater').pvalue:.3e}")
            print(f"      Wilcoxon signed-rank H1: T-move < G-flip  p = "
                  f"{wilcoxon(B, A, alternative='less').pvalue:.3e}")
        except Exception as e:
            print("      wilcoxon failed:", e)


if __name__ == "__main__":
    main()
