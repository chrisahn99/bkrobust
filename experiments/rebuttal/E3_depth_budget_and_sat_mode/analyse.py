"""E3 analysis: the three registered statistics of PREREG.md, from the run.py output.

    .venv/bin/python experiments/rebuttal/E3_depth_budget_and_sat_mode/analyse.py            # full run
    .venv/bin/python experiments/rebuttal/E3_depth_budget_and_sat_mode/analyse.py --smoke    # smoke run

Writes results.json next to the rows.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
BUDGETS = (0, 1, 2, 3, 4, 6)
MODES = ("ladder", "minimise")
REF_B = 3
N_BOOT = 10_000
EPS = 1e-6  # a SAT time is never 0; this only guards a log of a clock tick


def load(dir_: Path):
    rows = {}
    for f in sorted(dir_.glob("rows_shard*.jsonl")):
        for line in f.read_text().splitlines():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            rows[r["unit"]] = r  # the last line of a unit wins if a shard was rerun
    return [rows[k] for k in sorted(rows)]


def arm(row, b, mode):
    for a in row["arms"]:
        if a["b"] == b and a["mode"] == mode:
            return a
    raise KeyError((b, mode))


def boot_ci(stat, n, rng):
    if n == 0:
        return [None, None]
    vals = [stat(rng.integers(0, n, n)) for _ in range(N_BOOT)]
    vals = [v for v in vals if v is not None and math.isfinite(v)]
    return [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))] if vals else [None, None]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "results/rebuttal/E3"), help="run.py --out directory")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    dir_ = Path(args.out) / ("smoke" if args.smoke else "full")
    rows = load(dir_)
    ok = [r for r in rows if r.get("arms")]
    errors = [{"unit": r["unit"], "error": r.get("error")} for r in rows if not r.get("arms")]
    res = {"n_rows": len(rows), "n_with_arms": len(ok), "errors": errors,
           "criterion": sorted({r["criterion"] for r in rows}), "commit": sorted({r["commit"] for r in rows})}

    # 1. consistency
    bad = [{"unit": r["unit"], "failures": r["failures"]} for r in ok if not r["consistent"]]
    res["consistency"] = {"n_inconsistent": len(bad), "inconsistent": bad,
                          "n_capped_arms": sum(r["n_capped"] for r in ok),
                          "n_inexact_arms": sum(r["n_inexact"] for r in ok)}
    match_keys = ("m_undirected", "k_g0", "largest_component", "largest_component_density")
    res["reconstruction"] = {
        "graph_matches_file": sum(all(r["rebuilt_matches_file"][k] for k in match_keys) for r in ok),
        "z_size_matches_file": sum(r["query"]["z_size_matches_file"] for r in ok),
        "query_fallback": sum(r["query"]["source"] == "fallback" for r in ok),
        "of": len(ok),
    }
    radii = [r["radius"] for r in ok]
    res["radius_counts"] = {str(k): radii.count(k) for k in sorted(set(radii), key=lambda v: (v is None, v))}

    # 2. budget: T(b) per mode, ratio to b = 3, paired bootstrap over instances
    rng = np.random.default_rng(0)
    n = len(ok)
    tot = {(b, m): np.array([arm(r, b, m)["total_s"] for r in ok]) for b in BUDGETS for m in MODES}
    budget = {}
    for m in MODES:
        per_b = {}
        for b in BUDGETS:
            t, t3 = tot[(b, m)], tot[(REF_B, m)]
            ratio = float(t.sum() / t3.sum()) if n and t3.sum() > 0 else None
            methods = {}
            for r in ok:
                k = arm(r, b, m)["method"]
                methods[k] = methods.get(k, 0) + 1
            per_b[str(b)] = {
                "T_s": float(t.sum()), "median_s": float(np.median(t)) if n else None,
                "ratio_to_b3": ratio,
                "ratio_ci95": boot_ci(lambda i, t=t, t3=t3: float(t[i].sum() / t3[i].sum())
                                      if t3[i].sum() > 0 else None, n, rng),
                "n_capped": sum(arm(r, b, m)["capped"] for r in ok), "methods": methods,
            }
        best = min(BUDGETS, key=lambda b: per_b[str(b)]["T_s"]) if n else None
        budget[m] = {"per_b": per_b, "b_star": best}
    res["budget"] = budget

    # 3. SAT mode at b = 0, where every finite-radius instance reaches the SAT stage
    pairs = []
    for r in ok:
        la, mn = arm(r, 0, "ladder"), arm(r, 0, "minimise")
        if la["method"] == "e1_ladder" and mn["method"] == "e1_minimise" and la["exact"] and mn["exact"]:
            pairs.append((r["unit"], la["sat_s"], mn["sat_s"], r["radius"]))
    lr = np.array([math.log((p[2] + EPS) / (p[1] + EPS)) for p in pairs])
    gm = float(math.exp(lr.mean())) if len(lr) else None
    ci = boot_ci(lambda i: float(math.exp(lr[i].mean())), len(lr), rng)
    by_r = {}
    for p, v in zip(pairs, lr):
        by_r.setdefault(str(p[3]), []).append(v)
    res["sat_mode_b0"] = {
        "n_pairs": len(pairs), "geomean_minimise_over_ladder": gm, "ci95": ci,
        "sum_ladder_sat_s": float(sum(p[1] for p in pairs)), "sum_minimise_sat_s": float(sum(p[2] for p in pairs)),
        "geomean_by_radius": {k: float(math.exp(np.mean(v))) for k, v in sorted(by_r.items())},
        "n_by_radius": {k: len(v) for k, v in sorted(by_r.items())},
    }
    t3l, t3m = tot[(REF_B, "ladder")], tot[(REF_B, "minimise")]
    res["total_ratio_minimise_over_ladder_b3"] = float(t3m.sum() / t3l.sum()) if n and t3l.sum() > 0 else None

    out = dir_ / "results.json"
    out.write_text(json.dumps(res, indent=1))

    def f(v, w=8, d=3):
        return f"{v:{w}.{d}f}" if isinstance(v, float) else f"{str(v):>{w}}"

    print(f"E3  rows={len(rows)}  with arms={n}  errors={len(errors)}  criterion={res['criterion']}")
    print(f"consistency: inconsistent instances = {len(bad)}  capped arms = {res['consistency']['n_capped_arms']}"
          f"  inexact arms = {res['consistency']['n_inexact_arms']}")
    rc = res["reconstruction"]
    print(f"reconstruction: graph matches file {rc['graph_matches_file']}/{rc['of']}, |Z| matches "
          f"{rc['z_size_matches_file']}/{rc['of']}, fallback queries {rc['query_fallback']}")
    print(f"radii: {res['radius_counts']}")
    print("\n   b   mode        T(s)   median   T/T(3)  95% CI               capped")
    for m in MODES:
        for b in BUDGETS:
            p = budget[m]["per_b"][str(b)]
            lo, hi = p["ratio_ci95"]
            print(f"{b:4d}   {m:9s} {f(p['T_s'], 8, 2)} {f(p['median_s'], 8, 4)} {f(p['ratio_to_b3'], 7, 3)}"
                  f"  [{f(lo, 6, 3)}, {f(hi, 6, 3)}]  {p['n_capped']:6d}")
        print(f"     b* ({m}) = {budget[m]['b_star']}")
    s = res["sat_mode_b0"]
    lo, hi = s["ci95"]
    print(f"\nSAT stage at b=0, minimise / ladder (geometric mean over {s['n_pairs']} instances): "
          f"{f(gm, 0, 3)}  95% CI [{f(lo, 0, 3)}, {f(hi, 0, 3)}]")
    print(f"  by radius: {s['geomean_by_radius']}  n: {s['n_by_radius']}")
    print(f"total time at b=3, minimise / ladder: {f(res['total_ratio_minimise_over_ladder_b3'], 0, 3)}")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
