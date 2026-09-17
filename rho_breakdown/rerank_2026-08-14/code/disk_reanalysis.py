"""
Step 1/2 of E1: recompute rho* under all three criteria from data ALREADY ON DISK
(|K| <= 4). No simulation. Reproduces the published table and adds the full
distribution + interval widths that the published summary only bucketed.

Usage: python disk_reanalysis.py
"""
import json
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, ".")
from analyse import boot_ci, breakdown_radius, fmt, interval_at  # noqa: E402

RAW = "/Users/josecosta/research-pilots/latent-causal-rho-breakdown-knowledge/results/linear_raw.json"


def table(res, arm_name):
    rs = [r for r in res if arm_name in r["arms"]]
    out = {"n_scm": len(rs), "n_K": dict(sorted(Counter(
        r["arms"][arm_name]["n_K"] for r in rs).items()))}
    for crit in ("sign", "rel50", "ident"):
        rr = np.array([breakdown_radius(r["arms"][arm_name], r["tau"],
                                        r["arms"][arm_name]["est0"],
                                        min(3, r["arms"][arm_name]["n_K"]), crit)
                       for r in rs])
        n = len(rr)
        bins = {str(k): round(float((rr == k).sum()) / n, 4) for k in (1, 2, 3, 4)}
        out[crit] = {
            "bins_rho*_1_2_3_censored": bins,
            "counts": {str(k): int((rr == k).sum()) for k in (1, 2, 3, 4)},
            "censored_frac_ci": fmt(boot_ci((rr > 3).astype(float))),
            "max_bin": round(max(bins.values()), 4),
            "n_bins_ge_5pct": sum(1 for v in bins.values() if v >= 0.05),
            "median": float(np.median(rr)),
        }
    # interval widths
    out["widths"] = {}
    for rho in (1, 2, 3):
        wn = []
        for r in rs:
            arm = r["arms"][arm_name]
            if arm["n_K"] < rho:
                continue
            lo, hi, _, _ = interval_at(arm, rho, arm["est0"])
            wn.append((hi - lo) / abs(r["tau"]))
        wn = np.array(wn)
        out["widths"][f"rho{rho}"] = {
            "n": len(wn), "mean": round(float(wn.mean()), 4),
            "median": round(float(np.median(wn)), 4),
            "q90": round(float(np.percentile(wn, 90)), 4),
            "frac_zero_width": round(float((wn < 1e-9).mean()), 4),
        }
    return out


def main():
    res = json.load(open(RAW))
    S = {"source": RAW, "n_scm_total": len(res),
         "n_undirected_dist": dict(sorted(Counter(r["n_undirected"] for r in res).items()))}
    for arm in ("generic", "tiered"):
        S[arm] = table(res, arm)
    print(json.dumps(S, indent=1))
    with open("../results/disk_reanalysis.json", "w") as f:
        json.dump(S, f, indent=1)


if __name__ == "__main__":
    main()
