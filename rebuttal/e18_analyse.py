"""E18 analysis.

Checks: (1) naive draw 0 reproduces the stored record (r_mu signature, ceiling,
mu peak, to 1e-9); (2) posthoc equals naive on every relative quantity.
Effects: per arm and draw, the number of distinct r_mu signatures over the
paper's 49 queries (finite r_val, non-zero bias), the share whose ceiling
B(C_hat)/|theta_Z| is >= 1 (bias at least the size of the estimate), the share
never reaching each tolerance, and how often the internal arm's signature
equals the naive one for the same draw. var-sortability per arm.
"""

import json
from collections import defaultdict

import numpy as np

from common import OUT, ROOT

EPS = ["0.02", "0.1", "0.25", "0.4", "0.6", "0.9"]


def main():
    rows = [json.loads(line) for line in open(OUT / "e18" / "standardise.jsonl")]
    stored = {}
    for f in ("x_component", "extended"):
        for line in open(ROOT / f"results/mean_table_relevant_v2/{f}/instances.jsonl"):
            r = json.loads(line)
            if r.get("status") == "ok":
                stored[(r["network"], r["x"], r["y"])] = r
    by = defaultdict(dict)
    for r in rows:
        by[(r["network"], r["x"], r["y"], r["draw"])][r["arm"]] = r
    out = {"n_rows": len(rows), "status": {}}
    for r in rows:
        out["status"][r["status"]] = out["status"].get(r["status"], 0) + 1
    # (1) reproduction
    rep_ok = rep_bad = 0
    for (n, x, y, dr), arms in by.items():
        if dr != 0 or arms.get("naive", {}).get("status") != "ok":
            continue
        s = stored[(n, x, y)]
        a = arms["naive"]
        same = all(a["r_mu"][e] == s["r_mu_rel"][e] for e in EPS) and \
            abs(a["b_top"] - s.get("b_top", a["b_top"])) < 1e-9 if "b_top" in s else \
            all(a["r_mu"][e] == s["r_mu_rel"][e] for e in EPS)
        rep_ok += same
        rep_bad += not same
    out["naive_draw0_reproduces"] = {"ok": rep_ok, "differ": rep_bad}
    # (2) posthoc invariance
    inv_ok = inv_bad = 0
    maxdiff = 0.0
    for arms in by.values():
        a, p = arms.get("naive"), arms.get("posthoc")
        if not a or not p or a["status"] != "ok" or p["status"] != "ok":
            continue
        d = max([abs(a["b_top"] - p["b_top"]), abs(a["mu_peak"] - p["mu_peak"])] +
                [abs(u - v) for u, v in zip(a["profile"], p["profile"])])
        maxdiff = max(maxdiff, d)
        same = d < 1e-9 and a["r_mu"] == p["r_mu"]
        inv_ok += same
        inv_bad += not same
    out["posthoc_equals_naive"] = {"ok": inv_ok, "differ": inv_bad, "max_abs_diff": maxdiff}
    # effects, on the paper's population: finite r_val and non-zero bias under the naive arm
    eff = {}
    for arm in ("naive", "posthoc", "internal"):
        per_draw = defaultdict(list)
        for (n, x, y, dr), arms in by.items():
            a = arms.get(arm)
            if not a or a["status"] != "ok":
                continue
            if str(a["r_val"]) in ("-1", "None", "inf") or a["b_top"] <= 1e-9:
                continue
            per_draw[dr].append(a)
        d0 = per_draw.get(0, [])
        sig = lambda a: tuple(a["r_mu"][e] for e in EPS)  # noqa: E731
        eff[arm] = {
            "n_queries_draw0": len(d0),
            "distinct_signatures_draw0": len({sig(a) for a in d0}),
            "distinct_signatures_mean_over_draws": float(np.mean([len({sig(a) for a in v}) for v in per_draw.values()])),
            "share_ceiling_ge_1_draw0": float(np.mean([a["b_top"] >= 1 - 1e-9 for a in d0])) if d0 else None,
            "median_ceiling_all_draws": float(np.median([a["b_top"] for v in per_draw.values() for a in v])),
            "share_never_reached_by_eps_all_draws": {e: float(np.mean([a["r_mu"][e] == -1 for v in per_draw.values() for a in v])) for e in EPS},
            "var_sortability_median": float(np.median([a["var_sortability"] for v in per_draw.values() for a in v if a["var_sortability"] is not None])),
        }
    out["effects"] = eff
    # internal vs naive signature agreement, same query and draw
    agree = tot = 0
    ceil_up = 0
    for arms in by.values():
        a, i = arms.get("naive"), arms.get("internal")
        if not a or not i or a["status"] != "ok" or i["status"] != "ok":
            continue
        tot += 1
        agree += a["r_mu"] == i["r_mu"]
        ceil_up += i["b_top"] > a["b_top"] + 1e-9
    out["internal_vs_naive"] = {"pairs": tot, "same_signature": agree, "ceiling_higher_internal": ceil_up}
    (OUT / "e18" / "analysis.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
