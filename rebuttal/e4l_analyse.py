"""E4L analysis: the registered statistics of PREREGISTRATION.md (section E4L)."""

from __future__ import annotations

import json
import math
from collections import defaultdict

import numpy as np

from common import OUT

D = OUT / "e4l"
BINS = [(0, 0, "1"), (1, 10, "2-10"), (11, 100, "11-100"), (101, 10**3, "101-10^3"),
        (10**3 + 1, 10**4, "10^3-10^4"), (10**4 + 1, 10**6, "10^4-10^6"),
        (10**6 + 1, 10**9, "10^6-10^9"), (10**9 + 1, None, ">10^9")]
ORDER = [b[2] for b in BINS] + ["unknown"]


def sbin(r):
    s = r.get("class_size")
    if s is None:
        return "unknown"
    s = int(s)
    if s == 1:
        return "1"
    for lo, hi, lab in BINS[1:]:
        if s >= lo and (hi is None or s <= hi):
            return lab
    return ">10^9"


def wilson(k, n):
    if n == 0:
        return (None, None)
    z = 1.96
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (round(c - h, 3), round(c + h, 3))


def faster(r, leg):
    L = r.get(leg)
    return bool(L and L["complete"] and L["seconds"] < r["hyb_s"])


def main():
    rows = [json.loads(line) for line in open(D / "rows.jsonl")]
    ok = [r for r in rows if r.get("status") == "ok"]
    out = {"n_rows": len(rows),
           "status": {s: sum(1 for r in rows if r.get("status", "").split(":")[0] == s)
                      for s in {r.get("status", "").split(":")[0] for r in rows}},
           "consistency_failures": [r["uid"] for r in ok if r.get("consistency_fail")],
           "hybrid_not_exact": sum(1 for r in ok if not r.get("hyb_exact"))}
    tab = defaultdict(lambda: defaultdict(list))
    for r in ok:
        tab[r["pop"]][sbin(r)].append(r)
        tab["ALL"][sbin(r)].append(r)
    res = {}
    for pop, bins in tab.items():
        res[pop] = {}
        for b in ORDER:
            rs = bins.get(b, [])
            if not rs:
                continue
            n = len(rs)
            kl = sum(faster(r, "list") for r in rs)
            kb = sum(faster(r, "bnb") for r in rs)
            res[pop][b] = {
                "n": n,
                "list_faster": round(kl / n, 3), "list_faster_ci": wilson(kl, n),
                "bnb_faster": round(kb / n, 3), "bnb_faster_ci": wilson(kb, n),
                "list_complete": round(sum(r["list"]["complete"] for r in rs) / n, 3),
                "bnb_complete": round(sum(r["bnb"]["complete"] for r in rs) / n, 3),
                "hyb_median_s": float(np.median([r["hyb_s"] for r in rs])),
                "hyb_max_s": float(max(r["hyb_s"] for r in rs)),
                "list_projected_median_s": float(np.median([r["list_projected_s"] for r in rs
                                                             if r.get("list_projected_s")]))
                if any(r.get("list_projected_s") for r in rs) else None,
            }
    out["by_pop_bin"] = res
    proj = {}
    for pop in sorted({r["pop"] for r in ok}) + ["ALL"]:
        rs = [r for r in ok if pop == "ALL" or r["pop"] == pop]
        pj = [r["list_projected_s"] for r in rs if r.get("list_projected_s")]
        proj[pop] = {
            "n": len(rs), "n_projected": len(pj),
            "share_list_over_1h": round(np.mean([p > 3600 for p in pj]), 3) if pj else None,
            "share_list_over_1y": round(np.mean([p > 3.15e7 for p in pj]), 3) if pj else None,
            "share_hybrid_under_1s": round(np.mean([r["hyb_s"] < 1 for r in rs]), 3),
            "hyb_max_s": float(max(r["hyb_s"] for r in rs)),
            "list_faster_overall": round(np.mean([faster(r, "list") for r in rs]), 3),
            "bnb_faster_overall": round(np.mean([faster(r, "bnb") for r in rs]), 3),
            "list_complete_overall": round(np.mean([r["list"]["complete"] for r in rs]), 3),
            "bnb_complete_overall": round(np.mean([r["bnb"]["complete"] for r in rs]), 3),
            "class_size_log10_quantiles": [round(float(q), 1) for q in np.percentile(
                [r["class_size_log10"] for r in rs if r.get("class_size_log10") is not None],
                [10, 50, 90, 99])] if rs else None,
        }
    out["projection"] = proj
    # where does each method win on large classes (>10^4)
    big = [r for r in ok if r.get("class_size") and int(r["class_size"]) > 10**4]
    out["large_classes"] = {
        "n": len(big),
        "bnb_complete": round(np.mean([r["bnb"]["complete"] for r in big]), 3) if big else None,
        "bnb_faster": round(np.mean([faster(r, "bnb") for r in big]), 3) if big else None,
        "bnb_median_s_when_complete": float(np.median([r["bnb"]["seconds"] for r in big
                                                       if r["bnb"]["complete"]]))
        if any(r["bnb"]["complete"] for r in big) else None,
        "hyb_median_s": float(np.median([r["hyb_s"] for r in big])) if big else None,
    }
    # long-cap validation
    lp = D / "rows_long.jsonl"
    if lp.exists():
        main_by = {r["uid"]: r for r in ok}
        val = []
        for line in open(lp):
            r = json.loads(line)
            if r.get("status") != "ok" or r["uid"] not in main_by:
                continue
            m = main_by[r["uid"]]
            probe = m["list"]["dags_per_s"]
            longr = r["list"]["dags_per_s"]
            val.append({"uid": r["uid"], "size_log10": r["class_size_log10"],
                        "probe_rate": probe, "long_rate": longr, "ratio": longr / probe,
                        "long_complete": r["list"]["complete"], "long_seconds": r["list"]["seconds"],
                        "projected_s": m.get("list_projected_s")})
        rat = [v["ratio"] for v in val]
        out["long_cap_validation"] = {
            "n": len(val), "rate_ratio_median": float(np.median(rat)) if rat else None,
            "rate_ratio_min": float(min(rat)) if rat else None,
            "rate_ratio_max": float(max(rat)) if rat else None,
            "share_within_factor_2": round(np.mean([0.5 <= x <= 2 for x in rat]), 3) if rat else None,
            "completed": [v for v in val if v["long_complete"]],
        }
    (D / "analysis.json").write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps({k: v for k, v in out.items() if k != "by_pop_bin"}, indent=1, default=str))
    for pop in ("P-E2", "P-REL", "P-BIG", "P-SYN", "ALL"):
        if pop in res:
            print("==", pop)
            for b, v in res[pop].items():
                print(f"  {b:10s} n={v['n']:5d} list_faster={v['list_faster']} {v['list_faster_ci']} "
                      f"bnb_faster={v['bnb_faster']} list_complete={v['list_complete']} "
                      f"bnb_complete={v['bnb_complete']} hyb_med={v['hyb_median_s']:.4f} "
                      f"proj_med={v['list_projected_median_s']}")


if __name__ == "__main__":
    main()
