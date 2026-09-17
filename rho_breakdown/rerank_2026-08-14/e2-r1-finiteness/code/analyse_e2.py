"""
E2 analysis -- executes the PRE-REGISTERED decision rule in PREREGISTRATION.md.
Nothing here chooses a threshold; the thresholds were fixed before the run.

Usage: python analyse_e2.py [raw.json]
"""
import json
import sys
from collections import Counter, defaultdict

import numpy as np
from scipy.stats import mannwhitneyu

RNG = np.random.default_rng(20260814)


def boot_ci(v, fn=np.median, B=4000, alpha=0.05):
    v = np.asarray(v, dtype=float)
    if len(v) == 0:
        return (float("nan"), float("nan"))
    idx = RNG.integers(0, len(v), size=(B, len(v)))
    stats = np.array([fn(v[i]) for i in idx])
    return (float(np.quantile(stats, alpha / 2)), float(np.quantile(stats, 1 - alpha / 2)))


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def desc(v):
    v = np.asarray(v, dtype=float)
    if len(v) == 0:
        return dict(n=0)
    return dict(n=int(len(v)), mean=float(v.mean()), median=float(np.median(v)),
                q75=float(np.quantile(v, .75)), q90=float(np.quantile(v, .90)),
                q99=float(np.quantile(v, .99)), max=float(v.max()))


def main():
    raw_path = sys.argv[1] if len(sys.argv) > 1 else "../results/e2_raw.json"
    data = json.load(open(raw_path))
    ARMS = ("G-flip", "T-flip", "T-move")

    # ---------------- matched SCM set: every arm present
    matched = [r for r in data if all(a in r["arms"] for a in ARMS)]
    print("=" * 78)
    print("E2 -- R1-CLOSED FINITENESS, EMPIRICAL HALF")
    print("=" * 78)
    print(f"SCMs replayed from pilot ensemble : 600")
    print(f"SCMs retained (|K_T|>=2, O* exists): {len(data)}")
    print(f"SCMs with ALL THREE arms (matched) : {len(matched)}")
    nk = Counter()
    for r in matched:
        nk[r["arms"]["T-move"]["n_K"]] += 1
    print(f"|K| distribution (matched, both arms share it): {dict(sorted(nk.items()))}")
    print(f"tiering re-draws needed: {dict(sorted(Counter(r['n_attempts'] for r in matched).items()))}")

    # ---------------- collect members
    M = {a: [] for a in ARMS}
    for r in matched:
        for a in ARMS:
            for m in r["arms"][a]["members"]:
                m = dict(m)
                m["seed"] = r["seed"]
                m["tau"] = r["tau"]
                m["cpdag_amenable"] = r["cpdag_amenable"]
                M[a].append(m)

    print("\n" + "-" * 78)
    print("MEMBER COUNTS AND CONSISTENCY (all members of the matched set)")
    print("-" * 78)
    print(f"{'arm':<9}{'members':>9}{'consistent':>12}{'cons.rate':>11}{'n_chg=0':>9}")
    for a in ARMS:
        cons = [m for m in M[a] if m.get("consistent")]
        noop = [m for m in M[a] if m["n_changed"] == 0]
        print(f"{a:<9}{len(M[a]):>9}{len(cons):>12}{len(cons)/max(len(M[a]),1):>11.3f}"
              f"{len(noop):>9}")

    # ================================================================= GATE 1
    print("\n" + "=" * 78)
    print("GATE 1 -- PREMISE: is a PERTURBED tiering still R1-closed?")
    print("  pre-registered: s >= 0.99 -> premise holds ; s < 0.99 -> REFUTED")
    print("=" * 78)
    base_ok = sum(1 for r in matched if r["arms"]["T-move"]["r1_base_equals_full"])
    print(f"UNPERTURBED tiered K, R1-only closure == R1-R4 closure: "
          f"{base_ok}/{len(matched)} = {base_ok/len(matched):.4f}   "
          f"(Bang & Didelez sanity check)")
    base_ok_g = sum(1 for r in matched if r["arms"]["G-flip"]["r1_base_equals_full"])
    print(f"UNPERTURBED generic K, same check                    : "
          f"{base_ok_g}/{len(matched)} = {base_ok_g/len(matched):.4f}")

    gate1 = {}
    for a in ARMS:
        cons = [m for m in M[a] if m.get("consistent")]
        eq = sum(1 for m in cons if m.get("r1_equals_full"))
        lo, hi = wilson(eq, len(cons))
        gate1[a] = eq / max(len(cons), 1)
        print(f"PERTURBED {a:<7} R1-only == R1-R4 : {eq}/{len(cons)} = "
              f"{eq/max(len(cons),1):.4f}  [{lo:.4f},{hi:.4f}]")
    s = gate1["T-move"]
    GATE1 = "REFUTED" if s < 0.99 else "PASS"
    print(f"\n>>> GATE 1 = {GATE1}   (s = {s:.4f}, threshold 0.99)")

    # ================================================================= GATE 2
    print("\n" + "=" * 78)
    print("GATE 2 -- CONTAINMENT: orientation changes per unit perturbation")
    print("  r = d_orient / n_changed, consistent members with n_changed > 0")
    print("=" * 78)
    R = {}
    for a in ARMS:
        R[a] = [m["d_orient"] / m["n_changed"] for m in M[a]
                if m.get("consistent") and m["n_changed"] > 0]
    print(f"{'arm':<9}{'n':>7}{'mean':>9}{'median':>9}{'q75':>7}{'q90':>7}{'q99':>7}{'max':>7}")
    for a in ARMS:
        d = desc(R[a])
        print(f"{a:<9}{d['n']:>7}{d['mean']:>9.3f}{d['median']:>9.3f}"
              f"{d['q75']:>7.3f}{d['q90']:>7.3f}{d['q99']:>7.3f}{d['max']:>7.3f}")
    for a in ARMS:
        lo, hi = boot_ci(R[a])
        print(f"  median[{a}] 95% CI = [{lo:.3f}, {hi:.3f}]")

    rt, rg = np.array(R["T-move"]), np.array(R["G-flip"])
    p_greater = mannwhitneyu(rt, rg, alternative="greater").pvalue
    p_less = mannwhitneyu(rt, rg, alternative="less").pvalue
    print(f"\nMann-Whitney U  H1: r(T-move) >  r(G-flip)  p = {p_greater:.3e}")
    print(f"Mann-Whitney U  H1: r(T-move) <  r(G-flip)  p = {p_less:.3e}")

    med_t, med_g = np.median(rt), np.median(rg)
    q90_t, q90_g = np.quantile(rt, .90), np.quantile(rg, .90)
    max_t, max_g = rt.max(), rg.max()
    refuted = ((med_t > med_g) or (q90_t > q90_g)) and p_greater < 0.05
    supported = (med_t <= med_g and q90_t <= q90_g and max_t <= max_g and p_less < 0.05)
    GATE2 = "REFUTED" if refuted else ("SUPPORTED" if supported else "INCONCLUSIVE")
    print(f"\n  median  T-move {med_t:.3f}  vs  G-flip {med_g:.3f}   -> T-move "
          f"{'GREATER' if med_t > med_g else 'not greater'}")
    print(f"  q90     T-move {q90_t:.3f}  vs  G-flip {q90_g:.3f}   -> T-move "
          f"{'GREATER' if q90_t > q90_g else 'not greater'}")
    print(f"  max     T-move {max_t:.3f}  vs  G-flip {max_g:.3f}")
    print(f"\n>>> GATE 2 = {GATE2}")

    # ---- raw d_orient stratified by n_changed (rate normalisation is a design
    #      choice; this is the same comparison without it)
    print("\nSTRATIFIED BY n_changed (raw d_orient, no normalisation):")
    print(f"{'n_chg':>6}{'':2}{'G-flip n':>9}{'med':>6}{'q90':>6}{'max':>6}"
          f"{'':3}{'T-move n':>9}{'med':>6}{'q90':>6}{'max':>6}")
    for nc in range(1, 7):
        row = [nc]
        cells = []
        for a in ("G-flip", "T-move"):
            v = [m["d_orient"] for m in M[a]
                 if m.get("consistent") and m["n_changed"] == nc]
            cells.append(v)
        if len(cells[0]) < 10 and len(cells[1]) < 10:
            continue
        s_out = f"{nc:>6}{'':2}"
        for v in cells:
            if v:
                s_out += (f"{len(v):>9}{np.median(v):>6.2f}"
                          f"{np.quantile(v,.9):>6.2f}{max(v):>6.0f}{'':3}")
            else:
                s_out += f"{0:>9}{'-':>6}{'-':>6}{'-':>6}{'':3}"
        print(s_out)

    # ---- cascade, and the pilot's own metric for continuity
    print("\nCASCADE (rule-propagated orientations beyond those named in K'):")
    print(f"{'arm':<9}{'n':>7}{'mean':>9}{'median':>9}{'q90':>7}{'max':>7}"
          f"{'frac>0':>9}")
    for a in ARMS:
        v = [m["cascade"] for m in M[a] if m.get("consistent")]
        d = desc(v)
        f0 = np.mean([c > 0 for c in v]) if v else float("nan")
        print(f"{a:<9}{d['n']:>7}{d['mean']:>9.3f}{d['median']:>9.3f}"
              f"{d['q90']:>7.3f}{d['max']:>7.0f}{f0:>9.3f}")

    # ================================ CAPABILITY CHECK (C13, pre-registered)
    print("\n" + "=" * 78)
    print("C13 CAPABILITY CHECK -- could the instrument have returned the other verdict?")
    print("=" * 78)
    for a in ARMS:
        v = np.array(R[a])
        vals = Counter(np.round(v, 6))
        big = [(k, c / len(v)) for k, c in vals.items() if c / len(v) > 0.05]
        print(f"{a:<9} distinct r values with >5% mass: {len(big)}  "
              f"-> {sorted([(float(k), round(m,3)) for k, m in big])[:6]}")
    print("\nPOSITIVE CONTROL (pilot reported tiered cascade 0.760 vs generic 0.278,")
    print("50.1% vs 19.9% of flips cascading -- pipeline must reproduce the SIGN):")
    ct = [m["cascade"] for m in M["T-flip"] if m.get("consistent")]
    cg = [m["cascade"] for m in M["G-flip"] if m.get("consistent")]
    print(f"  cascade mean  T-flip {np.mean(ct):.3f}  vs  G-flip {np.mean(cg):.3f}  "
          f"-> {'REPRODUCED' if np.mean(ct) > np.mean(cg) else 'NOT REPRODUCED'}")
    print(f"  frac cascading T-flip {np.mean([c>0 for c in ct]):.3f}  vs  "
          f"G-flip {np.mean([c>0 for c in cg]):.3f}")

    # ================================================================= GATE 3
    print("\n" + "=" * 78)
    print("GATE 3 -- FINITENESS FLAVOUR (secondary, no threshold, distributions only)")
    print("=" * 78)
    print(f"{'arm':<9}{'SCMs':>6}{'distinct O* per ball: mean':>28}{'median':>8}{'q90':>6}{'max':>6}")
    for a in ARMS:
        counts = []
        for r in matched:
            os_ = {tuple(m["ostar"]) for m in r["arms"][a]["members"] if "ostar" in m}
            counts.append(len(os_))
        d = desc(counts)
        print(f"{a:<9}{d['n']:>6}{d['mean']:>28.3f}{d['median']:>8.1f}"
              f"{d['q90']:>6.1f}{d['max']:>6.0f}")

    print("\nDAMAGE TAIL |bias|/|tau| CONDITIONAL ON O* CHANGING "
          "(the campaign was burned by a mean):")
    print(f"{'arm':<9}{'n':>7}{'mean':>9}{'median':>9}{'q90':>9}{'q99':>9}{'max':>9}")
    for a in ARMS:
        v = [abs(m["bias"]) / abs(m["tau"]) for m in M[a]
             if m.get("ostar_changed")]
        d = desc(v)
        if d["n"]:
            print(f"{a:<9}{d['n']:>7}{d['mean']:>9.3f}{d['median']:>9.3f}"
                  f"{d['q90']:>9.3f}{d['q99']:>9.3f}{d['max']:>9.3f}")

    print("\nOUTCOME DECOMPOSITION (consistent members; 'silent' = O* changed & invalid):")
    print(f"{'arm':<9}{'N':>7}{'caught':>9}{'loud':>9}{'inert':>9}{'silent':>9}"
          f"{'chg_valid':>10}")
    for a in ARMS:
        allm = M[a]
        N = len(allm)
        caught = sum(1 for m in allm if not m.get("consistent"))
        loud = sum(1 for m in allm if m.get("consistent") and not m.get("amenable"))
        inert = sum(1 for m in allm if m.get("amenable") and not m.get("ostar_changed"))
        chg = [m for m in allm if m.get("ostar_changed")]
        silent = sum(1 for m in chg if not m.get("ostar_valid"))
        cv = sum(1 for m in chg if m.get("ostar_valid"))
        print(f"{a:<9}{N:>7}{caught/N:>9.3f}{loud/N:>9.3f}{inert/N:>9.3f}"
              f"{silent/N:>9.3f}{cv/N:>10.3f}")

    # ---- utility vs safety stratification (pilot's Theorem-1 candidate)
    print("\nUTILITY vs SAFETY (does damage live only where K is load-bearing?):")
    for a in ARMS:
        for lab, want in (("CPDAG already amenable (K = efficiency only)", True),
                          ("CPDAG NOT amenable (K load-bearing)", False)):
            sub = [m for m in M[a] if m["cpdag_amenable"] == want]
            if not sub:
                continue
            dmg = sum(1 for m in sub if m.get("ostar_changed") and not m.get("ostar_valid"))
            print(f"  {a:<8} {lab:<45} {dmg}/{len(sub)}")

    # ================================================================= VERDICT
    verdict = "REFUTED" if (GATE1 == "REFUTED" or GATE2 == "REFUTED") else (
        "SUPPORTED" if GATE2 == "SUPPORTED" else "INCONCLUSIVE")
    print("\n" + "=" * 78)
    print(f"GATE 1 = {GATE1}   GATE 2 = {GATE2}   =>   OVERALL VERDICT = {verdict}")
    print("=" * 78)

    json.dump(dict(gate1=gate1, gate1_verdict=GATE1, gate2_verdict=GATE2,
                   verdict=verdict,
                   rate={a: desc(R[a]) for a in ARMS},
                   p_greater=float(p_greater), p_less=float(p_less),
                   n_matched=len(matched)),
              open("../results/e2_summary.json", "w"), indent=1)


if __name__ == "__main__":
    main()
