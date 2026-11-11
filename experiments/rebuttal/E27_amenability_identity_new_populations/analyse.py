"""E27 analysis: the registered share r_val == r_amen, timing, radius distribution (see PREREG.md).

    .venv/bin/python experiments/rebuttal/E27_amenability_identity_new_populations/analyse.py --out results/rebuttal/E27
Writes <out>/results.json and prints a short report.
"""

from __future__ import annotations

import argparse
import collections
import glob
import json
import statistics
from pathlib import Path

from scipy.stats import beta

UNREACHED = -1


def cp(k: int, n: int) -> list[float | None]:
    if n == 0:
        return [None, None]
    lo = 0.0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(0.975, k + 1, n - k))
    return [round(lo, 4), round(hi, 4)]


def share(rows) -> dict:
    n = len(rows)
    k = sum(1 for r in rows if r["equal"])
    return {"n": n, "equal": k, "unequal": n - k, "share": round(k / n, 4) if n else None, "cp95": cp(k, n)}


def q(xs, p):
    if not xs:
        return None
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(p * len(xs)))]


def lab(r):
    return "inf" if r == UNREACHED else str(r)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    rows = {}
    for f in sorted(glob.glob(str(out / "units*.jsonl"))):
        for line in open(f):
            r = json.loads(line)
            rows[r["key"]] = r
    rows = list(rows.values())
    fails = [json.loads(l) for f in glob.glob(str(out / "consistency_failures*.jsonl")) for l in open(f)]
    errs = sum(len(open(f).read().splitlines()) for f in glob.glob(str(out / "errors*.jsonl")))
    census = json.loads((out / "census_population_A.json").read_text()) if (out / "census_population_A.json").exists() else None

    status = collections.Counter(r["status"] for r in rows)
    ok = [r for r in rows if r["status"] == "ok"]
    decided = [r for r in ok if r.get("both_exact")]
    censored = [r for r in ok if not r.get("both_exact")]
    fin = [r for r in decided if r["r_val"] != UNREACHED or r["r_amen"] != UNREACHED]
    new = [r for r in decided if not r.get("in_table1_frame")]
    by_cond = {c: share([r for r in decided if r["condition"] == c])
               for c in sorted({r["condition"] for r in decided})}
    by_net = {n: share([r for r in decided if r["network"] == n])
              for n in sorted({r["network"] for r in decided})}
    worst = share(decided + [{**r, "equal": False} for r in censored])

    g0rows = [r for r in rows if r.get("r_amen_call") is not None]
    t_val = [r["t_val_s"] for r in ok if r.get("r_val_call") == "ok"]
    t_amen = [r["t_amen_s"] for r in g0rows if r.get("r_amen_call") == "ok"]
    paired = [r["t_val_s"] / r["t_amen_s"] for r in ok
              if r.get("r_val_call") == "ok" and r.get("r_amen_call") == "ok" and r["t_amen_s"] > 0]
    enum_rows = [r for r in rows if r.get("enum_call") == "ok"]
    enum_amen_agree = sum(1 for r in enum_rows if r.get("r_amen_exact") and r["enum_r_amen"] == r["r_amen"])
    enum_amen_n = sum(1 for r in enum_rows if r.get("r_amen_exact"))
    enum_val = [r for r in enum_rows if r["status"] == "ok" and r.get("r_val_exact")]
    enum_val_agree = sum(1 for r in enum_val if r["enum_r_val"] == r["r_val"])

    res = {
        "registered": share(decided),
        "prediction": "share 1.000, zero disagreements",
        "prediction_met": bool(decided) and all(r["equal"] for r in decided) and not fails,
        "finite_units": share(fin),
        "outside_table1_frame": share(new),
        "censored_counted_as_disagreement": worst,
        "n_units": len(rows), "status_counts": dict(status), "n_ok": len(ok),
        "n_censored_ok": len(censored),
        "censor_calls": dict(collections.Counter((r.get("r_val_call"), r.get("r_amen_call")) for r in censored)),
        "by_condition": by_cond, "by_network": by_net,
        "unequal_units": [{k: r.get(k) for k in ("key", "r_val", "r_amen", "z_star", "k_g0", "r_val_method",
                                                   "amenable_at_g0", "enum_r_val", "enum_r_amen")}
                          for r in decided if not r["equal"]],
        "r_val_distribution": dict(sorted(collections.Counter(lab(r["r_val"]) for r in decided).items())),
        "r_amen_distribution_all_g0_units": dict(sorted(collections.Counter(
            lab(r["r_amen"]) for r in g0rows if r.get("r_amen_exact")).items())),
        "timing_s": {
            "r_val": {"n": len(t_val), "median": q(t_val, .5), "p90": q(t_val, .9), "max": max(t_val, default=None),
                      "total": round(sum(t_val), 2)},
            "r_amen": {"n": len(t_amen), "median": q(t_amen, .5), "p90": q(t_amen, .9),
                       "max": max(t_amen, default=None), "total": round(sum(t_amen), 2)},
            "paired_ratio_t_val_over_t_amen": {"n": len(paired),
                                               "median": round(statistics.median(paired), 3) if paired else None,
                                               "p10": q(paired, .1), "p90": q(paired, .9)},
            "speedup_bar_2_met": bool(paired) and statistics.median(paired) >= 2,
        },
        "enumeration_check": {"r_amen": [enum_amen_agree, enum_amen_n], "r_val": [enum_val_agree, len(enum_val)]},
        "consistency_failures": len(fails), "failure_kinds": dict(collections.Counter(f["kind"] for f in fails)),
        "errors": errs,
        "population_A": None if census is None else {k: census[k] for k in
                                                       ("n_files_manifest", "n_networks_outside_33",
                                                        "n_sampled_pairs", "n_scored_units")}
        | {"networks": [(n["network"], n["status"], n.get("undirected_edges"), n["sampled_pairs"])
                        for n in census["networks"]]},
    }
    (out / "results.json").write_text(json.dumps(res, indent=1, default=str))
    R = res["registered"]
    print(f"E27 registered: r_val == r_amen on {R['equal']} of {R['n']} decided units "
          f"(share {R['share']}, CP95 {R['cp95']}); unequal {R['unequal']}; censored {len(censored)}")
    print(f"  finite: {res['finite_units']}  outside Table 1 frame: {res['outside_table1_frame']}")
    print(f"  statuses: {dict(status)}")
    print(f"  r_val distribution: {res['r_val_distribution']}")
    print(f"  timing: {json.dumps(res['timing_s'])}")
    print(f"  enumeration agree: {res['enumeration_check']}; failures {len(fails)}; errors {errs}")
    print(f"  population A: {res['population_A']}")
    print(f"  prediction met: {res['prediction_met']}")


if __name__ == "__main__":
    main()
