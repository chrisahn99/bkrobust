"""
E1 analysis: rho* under sign / rel50 / ident, at |K| = 4, 6, 8, paired on the
same SCMs, plus ignorance-interval widths.

Usage: python analyse_ksweep.py ../results/ksweep_licensed.json licensed
"""
import json
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, ".")
from analyse import boot_ci, breakdown_radius, fmt  # noqa: E402

MAX_RHO = 3


def interval_at(arm, rho, est0):
    ests = [est0] + [m["est"] for m in arm["members"]
                     if m["rho"] <= rho and m.get("consistent") and m.get("amenable")]
    n_na = sum(1 for m in arm["members"]
               if m["rho"] <= rho and m.get("consistent") and not m.get("amenable"))
    return min(ests), max(ests), n_na


def arm_table(rs, key):
    """rs = SCMs that have an amenable MPDAG under EVERY arm (paired subset)."""
    out = {"n_scm": len(rs)}
    for crit in ("sign", "rel50", "ident"):
        rr = np.array([breakdown_radius(r["arms"][key], r["tau"],
                                        r["arms"][key]["est0"], MAX_RHO, crit)
                       for r in rs])
        n = len(rr)
        bins = {str(b): round(float((rr == b).sum()) / n, 4) for b in (1, 2, 3, 4)}
        out[crit] = {
            "bins": bins,
            "counts": {str(b): int((rr == b).sum()) for b in (1, 2, 3, 4)},
            "censored": round(float((rr > 3).mean()), 4),
            "censored_ci": fmt(boot_ci((rr > 3).astype(float))),
            "max_bin": round(max(bins.values()), 4),
            "n_bins_ge_5pct": sum(1 for v in bins.values() if v >= 0.05),
            "median": float(np.median(rr)),
        }
    out["widths"] = {}
    for rho in (1, 2, 3):
        wn = np.array([(lambda t: (t[1] - t[0]) / abs(r["tau"]))(
            interval_at(r["arms"][key], rho, r["arms"][key]["est0"])) for r in rs])
        out["widths"][f"rho{rho}"] = {
            "mean": round(float(wn.mean()), 4),
            "mean_ci": fmt(boot_ci(wn)),
            "median": round(float(np.median(wn)), 4),
            "q90": round(float(np.percentile(wn, 90)), 4),
            "frac_zero_width": round(float((wn < 1e-9).mean()), 4),
        }
    # member-level bookkeeping
    ms = [m for r in rs for m in r["arms"][key]["members"]]
    out["members_total"] = len(ms)
    out["frac_meek_inconsistent"] = round(
        float(np.mean([not m["consistent"] for m in ms])), 4)
    r1 = [m for m in ms if m["rho"] == 1]
    cons1 = [m for m in r1 if m["consistent"]]
    amen1 = [m for m in cons1 if m["amenable"]]
    out["rho1"] = {
        "n": len(r1),
        "frac_meek_inconsistent": round(float(np.mean([not m["consistent"] for m in r1])), 4),
        "frac_not_amenable_of_consistent": round(
            float(np.mean([not m["amenable"] for m in cons1])), 4) if cons1 else None,
        "frac_silent_bias_of_all": round(
            float(sum(1 for m in amen1 if m["ostar_changed"] and not m["ostar_valid"]) / len(r1)), 4),
        "changed_still_valid": int(sum(1 for m in amen1
                                       if m["ostar_changed"] and m["ostar_valid"])),
    }
    return out


def main():
    path = sys.argv[1]
    label = sys.argv[2]
    res = json.load(open(path))
    keys = ["K4", "K6", "K8"]

    S = {"source": path, "label": label, "n_scm_raw": len(res)}
    S["p_dist"] = dict(sorted(Counter(r["p"] for r in res).items()))
    S["deg_dist"] = dict(sorted(Counter(r["deg"] for r in res).items()))
    S["n_undirected"] = {
        "mean": round(float(np.mean([r["n_undirected"] for r in res])), 3),
        "median": float(np.median([r["n_undirected"] for r in res])),
        "min": int(min(r["n_undirected"] for r in res)),
        "max": int(max(r["n_undirected"] for r in res)),
    }
    S["mpdag_amenable_under_true_K"] = {
        k: round(float(np.mean([r["arms"].get(k, {}).get("mpdag_amenable", False)
                                for r in res])), 4) for k in keys}

    paired = [r for r in res
              if all(r["arms"].get(k, {}).get("mpdag_amenable") for k in keys)]
    S["n_paired"] = len(paired)
    S["arms"] = {k: arm_table(paired, k) for k in keys}
    print(json.dumps(S, indent=1))
    with open(f"../results/ksweep_analysis_{label}.json", "w") as f:
        json.dump(S, f, indent=1)


if __name__ == "__main__":
    main()
