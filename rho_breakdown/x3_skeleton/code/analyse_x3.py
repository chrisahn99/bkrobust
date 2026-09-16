"""X3 analysis -- exactly the statistics pre-registered in ../PREREG.md."""
import json
import sys
from collections import Counter, defaultdict

import numpy as np

Z = 1.959963984540054


def wilson(k, n, z=Z):
    if n == 0:
        return (float("nan"), float("nan"))
    ph = k / n
    d = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / d
    h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def cluster_boot(per_scm, B=4000, seed=20260831):
    """per_scm: list of (k, n) per SCM.  Paired SCM-level cluster bootstrap."""
    if not per_scm:
        return (float("nan"), float("nan"), float("nan"))
    arr = np.array(per_scm, dtype=float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(arr), size=(B, len(arr)))
    ks = arr[idx, 0].sum(axis=1)
    ns = arr[idx, 1].sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        rates = np.where(ns > 0, ks / ns, np.nan)
    lo, hi = np.nanpercentile(rates, [2.5, 97.5])
    return (float(lo), float(hi), float(np.nanstd(rates)))


def design_effect(k, n, per_scm):
    if n == 0 or k == 0:
        return float("nan")
    p = k / n
    se_srs = np.sqrt(p * (1 - p) / n)
    _, _, se_cl = cluster_boot(per_scm)
    return float((se_cl / se_srs) ** 2) if se_srs > 0 else float("nan")


def collect(scms, arm, hop_key):
    """Returns outcome Counter and per-SCM (silent, denom) for a hop bucket."""
    cnt = Counter()
    per_scm = []
    for s in scms:
        if s.get("status") != "ok":
            continue
        k = n = 0
        for r in s["perts"]:
            if not hop_key(r["hop"]):
                continue
            o = r[arm]["outcome"]
            cnt[o] += 1
            if o in ("unchanged", "changed_valid", "silent_bias"):
                n += 1
                if o == "silent_bias":
                    k += 1
        if n:
            per_scm.append((k, n))
    return cnt, per_scm


def rate_block(scms, arm, hop_key, label):
    cnt, per_scm = collect(scms, arm, hop_key)
    n = sum(cnt[o] for o in ("unchanged", "changed_valid", "silent_bias"))
    k = cnt["silent_bias"]
    wl, wh = wilson(k, n)
    bl, bh, _ = cluster_boot(per_scm)
    de = design_effect(k, n, per_scm)
    return dict(label=label, k=k, n=n, rate=(k / n if n else float("nan")),
                wilson=[wl, wh], cluster=[bl, bh], design_effect=de,
                n_scm=len(per_scm), outcomes=dict(cnt))


def main():
    out = {}
    for ens in ("original", "licensed", "large"):
        d = json.load(open(f"../results/x3_{ens}.json"))
        scms = d["scms"]
        ok = [s for s in scms if s.get("status") == "ok"]
        e = {"n_drawn": d["n_drawn"], "n_returned": len(scms), "n_analysed": len(ok),
             "status": dict(Counter(s.get("status", "?") for s in scms)),
             "cpdag_amenable": dict(k=sum(1 for s in ok if s["cpdag_amenable"]),
                                    n=len(ok)),
             "n_perturbations": sum(len(s["perts"]) for s in ok)}
        for arm in ("SK", "S0"):
            e[arm] = {
                "hop0":  rate_block(ok, arm, lambda h: h == 0, "hop 0"),
                "hop1":  rate_block(ok, arm, lambda h: h == 1, "hop 1"),
                "hop2p": rate_block(ok, arm, lambda h: h >= 2, "hop >= 2"),
                "hop1p": rate_block(ok, arm, lambda h: h >= 1, "hop >= 1"),
                "all":   rate_block(ok, arm, lambda h: True,   "all"),
            }
            r0 = e[arm]["hop0"]["rate"]
            r1 = e[arm]["hop1p"]["rate"]
            e[arm]["ratio_hop0_over_hop1p"] = (r0 / r1) if r1 else float("inf")
        # bias magnitude among silent-bias events, arm SK
        rel = defaultdict(list)
        for s in ok:
            for r in s["perts"]:
                if r["SK"]["outcome"] == "silent_bias":
                    rel["h0" if r["hop"] == 0 else "h1p"].append(r["SK"]["rel_bias"])
        e["rel_bias_SK"] = {k: dict(n=len(v), median=float(np.median(v)),
                                    q90=float(np.quantile(v, 0.9)),
                                    max=float(np.max(v))) for k, v in rel.items() if v}
        # weight of the deleted edge, silent vs not (arm SK)
        ws, wn = [], []
        for s in ok:
            for r in s["perts"]:
                o = r["SK"]["outcome"]
                if o == "silent_bias":
                    ws.append(r["w"])
                elif o in ("unchanged", "changed_valid"):
                    wn.append(r["w"])
        if ws and wn:
            e["weight_deleted_edge"] = dict(
                silent=dict(n=len(ws), median=float(np.median(ws))),
                not_silent=dict(n=len(wn), median=float(np.median(wn))))
        out[ens] = e

    # ---- pre-registered decision rule: arm SK, ensemble licensed --------------
    L = out["licensed"]["SK"]
    r0, r1 = L["hop0"]["rate"], L["hop1p"]["rate"]
    ratio = L["ratio_hop0_over_hop1p"]
    if ratio >= 5 and r1 < 0.02:
        verdict = "SUPPORTED"
    elif ratio < 2 or r1 >= 0.05:
        verdict = "REFUTED"
    else:
        verdict = "INCONCLUSIVE"
    out["_verdict"] = dict(rule="arm SK, ensemble licensed",
                           r_hop0=r0, r_hop1p=r1, ratio=ratio, verdict=verdict)
    # ---- gates ---------------------------------------------------------------
    out["_gates"] = dict(
        G1_placebo="0.0000 by construction (eligibility requires the unperturbed O* to be valid in D)",
        G3_dynamic_range=dict(overall_rate=out["licensed"]["SK"]["all"]["rate"],
                              ok=0.0 < out["licensed"]["SK"]["all"]["rate"] < 1.0),
        G4_hop_mass=dict(n_hop1p=L["hop1p"]["n"], ok=L["hop1p"]["n"] >= 200),
        G5_not_amenable_reported=True,
    )
    json.dump(out, open("../results/x3_analysis.json", "w"), indent=1)

    # ---- print ---------------------------------------------------------------
    for ens in ("original", "licensed", "large"):
        e = out[ens]
        print(f"\n════════ {ens}  analysed={e['n_analysed']}  perturbations={e['n_perturbations']}")
        ca = e["cpdag_amenable"]
        print(f"  CPDAG amenable without K: {ca['k']}/{ca['n']} = {ca['k']/ca['n']:.4f}")
        for arm in ("SK", "S0"):
            print(f"  --- arm {arm} ---")
            for key in ("hop0", "hop1", "hop2p", "hop1p", "all"):
                b = e[arm][key]
                if b["n"] == 0:
                    print(f"    {b['label']:>9}  (no mass)"); continue
                print(f"    {b['label']:>9}  silent {b['k']:6}/{b['n']:6} = {b['rate']:.4f} "
                      f" wilson[{b['wilson'][0]:.4f},{b['wilson'][1]:.4f}] "
                      f" cluster[{b['cluster'][0]:.4f},{b['cluster'][1]:.4f}] "
                      f" deff={b['design_effect']:.2f}")
            print(f"    ratio hop0 / hop>=1 = {e[arm]['ratio_hop0_over_hop1p']:.2f}")
        if "weight_deleted_edge" in e:
            w = e["weight_deleted_edge"]
            print(f"  |w| of deleted edge: silent median {w['silent']['median']:.4f} "
                  f"(n={w['silent']['n']})  vs not-silent {w['not_silent']['median']:.4f} "
                  f"(n={w['not_silent']['n']})")
        if e.get("rel_bias_SK"):
            for k, v in e["rel_bias_SK"].items():
                print(f"  rel|bias| {k}: n={v['n']} median={v['median']:.4f} "
                      f"q90={v['q90']:.4f} max={v['max']:.2f}")
    print("\n════════ OUTCOME TAXONOMY, arm SK, licensed, all hops")
    for k, v in sorted(out["licensed"]["SK"]["all"]["outcomes"].items(),
                       key=lambda kv: -kv[1]):
        tot = sum(out["licensed"]["SK"]["all"]["outcomes"].values())
        print(f"    {k:18} {v:7}  {v/tot:.4f}")
    print("\n════════ PRE-REGISTERED VERDICT")
    v = out["_verdict"]
    print(f"    r(hop 0) = {v['r_hop0']:.4f}   r(hop>=1) = {v['r_hop1p']:.4f}   "
          f"ratio = {v['ratio']:.2f}")
    print(f"    ==> {v['verdict']}")
    print("\n════════ GATES"); print("   ", json.dumps(out["_gates"]))


if __name__ == "__main__":
    main()
