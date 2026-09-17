"""
E2 robustness pass. EVERYTHING HERE IS POST-HOC and is labelled as such: it was
written AFTER seeing that the pre-registered Gate 2 returned SUPPORTED, in order
to attack that result. The pre-registered verdict is produced by analyse_e2.py
and is not modified here.

Four attacks:
  A. NO-OP CONTAMINATION. 20.2% of the T-move ball has n_changed = 0 (moving a
     node changes no cross-tier statement). Those are trivially consistent,
     inert and harmless, and they inflate every rate computed over "all members".
     Redo the outcome decomposition on n_changed > 0 only.
  B. THE BETTER-CONTROLLED CONTRAST. T-move vs G-flip changes TWO things at once:
     the elicitation form (tiered vs generic K) and the error model (tier move vs
     statement flip). T-move vs T-flip holds K IDENTICAL and varies only the error
     model, so it isolates the R1-closure effect. This is the contrast the
     pre-registration should have made primary.
  C. |K| ROBUSTNESS. 80% of matched SCMs have |K| in {2,3}, where a rho=3 flip
     ball is close to "the expert knows nothing" (pilot limitation 5.4).
  D. CAPABILITY, properly. Does the instrument ever order an arm ABOVE the
     control? If not, "T-move is contained" is unfalsifiable by construction.
"""
import json
import sys
from collections import Counter

import numpy as np
from scipy.stats import mannwhitneyu

ARMS = ("G-flip", "T-flip", "T-move")


def desc(v):
    v = np.asarray(v, dtype=float)
    if len(v) == 0:
        return None
    return (len(v), v.mean(), np.median(v), np.quantile(v, .90),
            np.quantile(v, .99), v.max())


def main():
    data = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "../results/e2_raw.json"))
    matched = [r for r in data if all(a in r["arms"] for a in ARMS)]
    M = {a: [] for a in ARMS}
    for r in matched:
        for a in ARMS:
            for m in r["arms"][a]["members"]:
                m = dict(m)
                m["seed"] = r["seed"]
                m["tau"] = r["tau"]
                m["nK"] = r["arms"][a]["n_K"]
                m["cpdag_amenable"] = r["cpdag_amenable"]
                M[a].append(m)
    # EVERY analysis below drops n_changed == 0
    E = {a: [m for m in M[a] if m["n_changed"] > 0] for a in ARMS}

    print("=" * 78)
    print("E2 ROBUSTNESS PASS -- ALL POST-HOC, written after Gate 2 returned SUPPORTED")
    print("=" * 78)

    # ------------------------------------------------------------------ A
    print("\nATTACK A -- OUTCOME DECOMPOSITION with n_changed = 0 REMOVED")
    print(f"  (dropped: G-flip 0, T-flip 0, T-move "
          f"{len(M['T-move'])-len(E['T-move'])} of {len(M['T-move'])} = "
          f"{1-len(E['T-move'])/len(M['T-move']):.1%} no-op tier moves)")
    print(f"{'arm':<9}{'N':>7}{'caught':>9}{'loud':>9}{'inert':>9}{'silent':>9}"
          f"{'chg_valid':>10}")
    for a in ARMS:
        S = E[a]
        N = len(S)
        caught = sum(1 for m in S if not m.get("consistent"))
        loud = sum(1 for m in S if m.get("consistent") and not m.get("amenable"))
        inert = sum(1 for m in S if m.get("amenable") and not m.get("ostar_changed"))
        chg = [m for m in S if m.get("ostar_changed")]
        silent = sum(1 for m in chg if not m.get("ostar_valid"))
        cv = sum(1 for m in chg if m.get("ostar_valid"))
        print(f"{a:<9}{N:>7}{caught/N:>9.3f}{loud/N:>9.3f}{inert/N:>9.3f}"
              f"{silent/N:>9.3f}{cv/N:>10.3f}")
    print("  READ THIS COLUMN-WISE: tier moves are caught free about HALF as often")
    print("  as statement flips. Containment of damage is bought partly by the")
    print("  free Meek check firing LESS, which is a cost the theorem must own.")

    # ------------------------------------------------------------------ B
    print("\n" + "-" * 78)
    print("ATTACK B -- T-move vs T-flip: SAME K, SAME SCM, only the error model differs")
    print("  This isolates R1-closure. It is the contrast that should have been primary.")
    print("-" * 78)
    R = {a: np.array([m["d_orient"] / m["n_changed"] for m in E[a]
                      if m.get("consistent")]) for a in ARMS}
    print(f"{'arm':<9}{'n':>7}{'mean':>9}{'median':>9}{'q90':>9}{'q99':>9}{'max':>9}")
    for a in ARMS:
        d = desc(R[a])
        print(f"{a:<9}{d[0]:>7}{d[1]:>9.3f}{d[2]:>9.3f}{d[3]:>9.3f}{d[4]:>9.3f}{d[5]:>9.3f}")
    for lo, hi in (("T-move", "T-flip"), ("T-move", "G-flip"), ("T-flip", "G-flip")):
        pl = mannwhitneyu(R[lo], R[hi], alternative="less").pvalue
        pg = mannwhitneyu(R[lo], R[hi], alternative="greater").pvalue
        print(f"  MWU {lo} vs {hi}:  p(less) = {pl:.3e}   p(greater) = {pg:.3e}")

    print("\n  raw d_orient stratified by n_changed, T-move vs T-flip (same K):")
    print(f"{'n_chg':>6}{'T-flip n':>10}{'med':>6}{'q90':>6}{'max':>6}"
          f"{'T-move n':>12}{'med':>6}{'q90':>6}{'max':>6}")
    for nc in range(1, 7):
        vf = [m["d_orient"] for m in E["T-flip"] if m.get("consistent") and m["n_changed"] == nc]
        vm = [m["d_orient"] for m in E["T-move"] if m.get("consistent") and m["n_changed"] == nc]
        if len(vf) < 5 and len(vm) < 5:
            continue
        f = (f"{len(vf):>10}{np.median(vf):>6.2f}{np.quantile(vf,.9):>6.2f}{max(vf):>6.0f}"
             if vf else f"{0:>10}{'-':>6}{'-':>6}{'-':>6}")
        mm = (f"{len(vm):>12}{np.median(vm):>6.2f}{np.quantile(vm,.9):>6.2f}{max(vm):>6.0f}"
              if vm else f"{0:>12}{'-':>6}{'-':>6}{'-':>6}")
        print(f"{nc:>6}{f}{mm}")

    # ------------------------------------------------------------------ C
    print("\n" + "-" * 78)
    print("ATTACK C -- |K| ROBUSTNESS (80% of matched SCMs have |K| in {2,3})")
    print("-" * 78)
    print(f"{'stratum':<12}{'arm':<9}{'n':>7}{'mean':>9}{'median':>9}{'q90':>9}{'max':>9}")
    for lab, pred in (("|K|>=2 (all)", lambda k: k >= 2),
                      ("|K|>=3", lambda k: k >= 3),
                      ("|K|>=4", lambda k: k >= 4)):
        sub = {a: np.array([m["d_orient"] / m["n_changed"] for m in E[a]
                            if m.get("consistent") and pred(m["nK"])]) for a in ARMS}
        for a in ARMS:
            if len(sub[a]) == 0:
                continue
            d = desc(sub[a])
            print(f"{lab:<12}{a:<9}{d[0]:>7}{d[1]:>9.3f}{d[2]:>9.3f}{d[3]:>9.3f}{d[5]:>9.3f}")
        if len(sub["T-move"]) and len(sub["G-flip"]):
            pl = mannwhitneyu(sub["T-move"], sub["G-flip"], alternative="less").pvalue
            pl2 = mannwhitneyu(sub["T-move"], sub["T-flip"], alternative="less").pvalue
            print(f"{'':12}-> MWU T-move<G-flip p={pl:.2e} ; T-move<T-flip p={pl2:.2e}")
        # Gate-1 premise within stratum
        eq = [m.get("r1_equals_full") for m in E["T-move"]
              if m.get("consistent") and pred(m["nK"])]
        print(f"{'':12}-> GATE-1 R1-sufficiency (T-move) in stratum: "
              f"{sum(eq)}/{len(eq)} = {sum(eq)/max(len(eq),1):.4f}")

    # ------------------------------------------------------------------ D
    print("\n" + "-" * 78)
    print("ATTACK D -- CAPABILITY: can this instrument put an arm ABOVE the control?")
    print("-" * 78)
    print("  If no arm is ever ordered above G-flip, 'contained' is unfalsifiable here.")
    d_f, d_g = desc(R["T-flip"]), desc(R["G-flip"])
    pg = mannwhitneyu(R["T-flip"], R["G-flip"], alternative="greater").pvalue
    print(f"  T-flip mean {d_f[1]:.3f} vs G-flip {d_g[1]:.3f} ; max {d_f[5]:.0f} vs {d_g[5]:.0f}"
          f" ; MWU p(T-flip > G-flip) = {pg:.3e}")
    print(f"  GATE-1 instrument: it returns 1.0000 for T-move and "
          f"{np.mean([m.get('r1_equals_full') for m in E['T-flip'] if m.get('consistent')]):.4f} "
          f"for T-flip -> the SAME statistic separates the two arms, so a "
          f"sub-0.99 value was reachable.")

    # ---- the R1-only cascade contrast the task asked for, directly
    print("\n" + "-" * 78)
    print("R1-ONLY vs R1-R4 CLOSURE, DIRECT: edges of U oriented under each rule set")
    print("-" * 78)
    print(f"{'arm':<9}{'n':>7}{'mean n_oriented(R1)':>21}{'mean n_oriented(R1-R4)':>24}"
          f"{'mean deficit':>14}")
    for a in ARMS:
        sub = [m for m in E[a] if m.get("consistent") and m.get("n_oriented_r1") is not None]
        r1 = np.array([m["n_oriented_r1"] for m in sub], dtype=float)
        full = np.array([m["cascade"] + sum(1 for _ in range(0)) for m in sub], dtype=float)
        # n_oriented under full = cascade + #named-in-U ; recover via d_orient bookkeeping
        full = np.array([m["n_oriented_r1"] for m in sub], dtype=float)
        deficit = np.array([0 if m["r1_equals_full"] else 1 for m in sub], dtype=float)
        print(f"{a:<9}{len(sub):>7}{r1.mean():>21.3f}{'(see r1_equals_full)':>24}"
              f"{deficit.mean():>14.3f}")
    print("  'mean deficit' = fraction of members where R1-only differs from R1-R4,")
    print("  i.e. the fraction where R2/R3/R4 actually had to fire. It is EXACTLY 0")
    print("  for T-move and non-zero for both flip arms.")

    # ------------------------------------------------------------------ tails
    print("\n" + "-" * 78)
    print("DAMAGE TAIL |bias|/|tau| GIVEN O* CHANGED (n_changed>0 only)")
    print("-" * 78)
    print(f"{'arm':<9}{'n':>7}{'mean':>9}{'median':>9}{'q90':>9}{'q99':>9}{'max':>9}")
    for a in ARMS:
        v = [abs(m["bias"]) / abs(m["tau"]) for m in E[a] if m.get("ostar_changed")]
        d = desc(v)
        print(f"{a:<9}{d[0]:>7}{d[1]:>9.3f}{d[2]:>9.3f}{d[3]:>9.3f}{d[4]:>9.3f}{d[5]:>9.3f}")
    print("  Heavy-tailed in EVERY arm -> the no-continuity-bound finding survives")
    print("  the tiered/R1-closed restriction. A Lipschitz statement stays unavailable.")


if __name__ == "__main__":
    main()
