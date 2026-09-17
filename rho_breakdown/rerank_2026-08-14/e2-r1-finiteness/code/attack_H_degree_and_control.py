"""
Third probe.
(i) Is the residual pure-reversal containment a DEGREE-SELECTION artefact?
    A pure-reversal tier move requires every undirected neighbour of the moved node
    to lie STRICTLY BETWEEN the old and the new tier -> selects low-degree nodes.
    R1 chains along chordless undirected paths, so a reversal on a low-degree node
    cascades less FOR REASONS OF DEGREE, not R1-closure.
(ii) Does the PRE-REGISTERED POSITIVE CONTROL (T-flip must cascade MORE than G-flip)
     actually pass? The pre-registration makes the whole verdict unreadable if not.
"""
import gzip
import json
from collections import defaultdict

import numpy as np
from scipy.stats import mannwhitneyu

from graphs import dag_to_cpdag, random_dag, undirected_edges
from run_e2 import K_from_tiering

ARMS = ("G-flip", "T-flip", "T-move")


def main():
    data = json.load(gzip.open("../results/e2_raw.json.gz"))
    matched = [r for r in data if all(a in r["arms"] for a in ARMS)]

    tmove, gflip, tflip = [], [], []
    for r in matched:
        D = random_dag(r["p"], r["deg"], np.random.default_rng(r["seed"]))
        C = dag_to_cpdag(D)
        U = undirected_edges(C)
        udeg = defaultdict(int)
        for (i, j) in U:
            udeg[i] += 1
            udeg[j] += 1
        tier = np.array(r["tier"])
        K = [tuple(e) for e in r["arms"]["T-move"]["K"]]
        Kset = {frozenset(e): e for e in K}
        for m in r["arms"]["T-move"]["members"]:
            tp = tier.copy()
            for (v, t) in m["moves"]:
                tp[v] = t
            Kp = K_from_tiering(tp, C)
            Kpset = {frozenset(e): e for e in Kp}
            rev = [k for k in Kset if k in Kpset and Kset[k] != Kpset[k]]
            rec = dict(m)
            rec.update(seed=r["seed"], n_rev=len(rev),
                       n_rem=sum(1 for k in Kset if k not in Kpset),
                       n_add=sum(1 for k in Kpset if k not in Kset))
            # degree signature of the reversed edge(s)
            rec["revdeg"] = sorted(
                sum(udeg[n] for n in k) for k in rev) if rev else []
            tmove.append(rec)
        for arm, acc in (("G-flip", gflip), ("T-flip", tflip)):
            Ka = [tuple(e) for e in r["arms"][arm]["K"]]
            for m in r["arms"][arm]["members"]:
                rec = dict(m)
                rec.update(seed=r["seed"], n_rev=m["rho"], n_rem=0, n_add=0)
                acc.append(rec)
        # recover which statements each flip reversed: rebuild the ball deterministically
        from run_e2 import flip_ball
        for arm, acc in (("G-flip", gflip), ("T-flip", tflip)):
            Ka = [tuple(e) for e in r["arms"][arm]["K"]]
            ball = flip_ball(Ka, np.random.default_rng(r["seed"] + 11))
            base = len(acc) - len(ball)
            for idx, (rho, Kp) in enumerate(ball):
                Kps = {frozenset(e): e for e in Kp}
                Kas = {frozenset(e): e for e in Ka}
                rev = [k for k in Kas if Kas[k] != Kps[k]]
                acc[base + idx]["revdeg"] = sorted(
                    sum(udeg[n] for n in k) for k in rev)

    ok = lambda L: [x for x in L if x.get("consistent") and x["n_changed"] > 0]
    tm, gf, tf = ok(tmove), ok(gflip), ok(tflip)
    pure1 = [x for x in tm if x["n_changed"] == 1 and x["n_rev"] == 1]
    gf1 = [x for x in gf if x["n_changed"] == 1]
    tf1 = [x for x in tf if x["n_changed"] == 1]

    print("=" * 96)
    print("(i) SINGLE REVERSAL, STRATIFIED BY UNDIRECTED DEGREE OF THE REVERSED EDGE")
    print("    deg = udeg(u)+udeg(v) in the undirected subgraph of the CPDAG")
    print("=" * 96)
    print(f"{'deg':>5} | {'T-flip n':>9}{'mean':>8}{'q90':>7} | "
          f"{'G-flip n':>9}{'mean':>8}{'q90':>7} | {'T-movePURE n':>13}{'mean':>8}{'q90':>7}")
    allsets = {}
    for lab, S in (("tf", tf1), ("gf", gf1), ("tm", pure1)):
        d = defaultdict(list)
        for x in S:
            if len(x["revdeg"]) == 1:
                d[x["revdeg"][0]].append(x["d_orient"])
        allsets[lab] = d
    degs = sorted(set().union(*[set(v) for v in allsets.values()]))
    tot = {k: [0, 0.0] for k in allsets}
    for dg in degs:
        row = f"{dg:>5} |"
        for lab in ("tf", "gf", "tm"):
            v = np.array(allsets[lab].get(dg, []), float)
            if len(v):
                row += f"{len(v):>9}{v.mean():>8.3f}{np.quantile(v,.9):>7.2f} |"
            else:
                row += f"{'-':>9}{'-':>8}{'-':>7} |"
        print(row)
    # degree-reweighted: reweight T-move PURE to the G-flip degree distribution
    print("\n  DEGREE-REWEIGHTED single-reversal mean d_orient")
    gfw = {k: len(v) for k, v in allsets["gf"].items()}
    n = sum(gfw.values())
    for lab in ("tf", "gf", "tm"):
        num, den = 0.0, 0.0
        for dg, w in gfw.items():
            v = allsets[lab].get(dg, [])
            if v:
                num += w * np.mean(v)
                den += w
        print(f"    {lab:<4} reweighted to G-flip degree profile: {num/den:.3f} "
              f"(coverage {den/n:.3f})")
    print("\n  RAW (unweighted) single-reversal mean d_orient")
    for lab in ("tf", "gf", "tm"):
        v = np.concatenate([np.array(x, float) for x in allsets[lab].values()])
        print(f"    {lab:<4} {v.mean():.3f}  n={len(v)}  q90={np.quantile(v,.9):.2f}")

    print("\n" + "=" * 96)
    print("(ii) PRE-REGISTERED POSITIVE CONTROL: 'T-flip must cascade MORE than G-flip'")
    print("     PREREGISTRATION.md sec.4: if the pipeline cannot reproduce it, the")
    print("     instrument is broken and NO VERDICT IS READABLE.")
    print("=" * 96)
    cg = np.array([x["cascade"] for x in gf], float)
    ct = np.array([x["cascade"] for x in tf], float)
    print(f"  cascade  G-flip mean {cg.mean():.3f} (n={len(cg)})   "
          f"T-flip mean {ct.mean():.3f} (n={len(ct)})")
    print(f"  MWU H1: cascade(T-flip) > cascade(G-flip)  p = "
          f"{mannwhitneyu(ct, cg, alternative='greater').pvalue:.4f}")
    rg = np.array([x["d_orient"] / x["n_changed"] for x in gf], float)
    rt = np.array([x["d_orient"] / x["n_changed"] for x in tf], float)
    print(f"  MWU H1: r(T-flip) > r(G-flip)              p = "
          f"{mannwhitneyu(rt, rg, alternative='greater').pvalue:.4f}")
    print(f"  pilot ratio to reproduce: 0.760/0.278 = {0.760/0.278:.2f}x ; "
          f"measured here: {ct.mean()/cg.mean():.2f}x")

    print("\n" + "=" * 96)
    print("(iii) THE HEADLINE NUMBER: q90 at n_changed==1, by operation")
    print("=" * 96)
    for lab, S in (("G-flip (1 reversal)", gf1), ("T-flip (1 reversal)", tf1),
                   ("T-move ALL", [x for x in tm if x["n_changed"] == 1]),
                   ("T-move 1 reversal", pure1),
                   ("T-move 1 removal",
                    [x for x in tm if x["n_changed"] == 1 and x["n_rem"] == 1
                     and x["n_rev"] == 0 and x["n_add"] == 0]),
                   ("T-move 1 addition",
                    [x for x in tm if x["n_changed"] == 1 and x["n_add"] == 1
                     and x["n_rev"] == 0 and x["n_rem"] == 0])):
        v = np.array([x["d_orient"] for x in S], float)
        print(f"  {lab:<24} n={len(v):>6}  mean {v.mean():.3f}  med "
              f"{np.median(v):.2f}  q90 {np.quantile(v,.9):.2f}  max {v.max():.0f}")


if __name__ == "__main__":
    main()
