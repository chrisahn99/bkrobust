"""E16 analysis: registered statistic, table, results.json. Works on full and smoke output.

    python experiments/rebuttal/E16_radius_where_separation_undefined/analyse.py --out results/rebuttal/E16
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
REGISTERED = "submitted/pooled_tau_radius_undefined"
SECONDARY = ("corrected/pooled_tau_radius_undefined", "submitted/within_tau_radius_undefined",
             "corrected/within_tau_radius_undefined",
             "corrected/within_paired_radius_minus_separation_matched")


def fmt(res: dict) -> str:
    p, lo, hi = res.get("point"), res.get("ci_lo_2p5"), res.get("ci_hi_97p5")
    if p is None:
        return "n/a"
    ci = f"[{lo:+.3f}, {hi:+.3f}]" if lo is not None and hi is not None else "[n/a]"
    return f"{p:+.3f} {ci}"


def n_of(res: dict):
    return res.get("n_units", res.get("n", res.get("n_separation_undefined")))


def verdict(res: dict) -> str:
    lo, hi = res.get("ci_lo_2p5"), res.get("ci_hi_97p5")
    if lo is None or hi is None:
        return "no interval"
    if lo > 0:
        return "interval above zero"
    if hi < 0:
        return "interval below zero"
    return "interval contains zero"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(REPO_ROOT / "results/rebuttal/E16"))
    args = ap.parse_args(argv)
    out = Path(args.out)
    rows = {}
    for f in sorted(glob.glob(str(out / "rows*.jsonl"))):
        for line in open(f, encoding="utf-8"):
            line = line.strip()
            if line:
                r = json.loads(line)
                rows[r["key"]] = r  # the last line written for a cell wins (resumed runs)
    if not rows:
        print(f"no rows in {out}", file=sys.stderr)
        return 1

    failures = {k: r["consistency_failures"] for k, r in rows.items() if r["consistency_failures"]}
    ci_mismatch = {k: r["reproduction"]["ci_published"] for k, r in rows.items()
                   if r["reproduction"].get("ci_ok") is False}
    print(f"{'cell':58s} {'n':>6s}  tau_b or lead [95% CI]")
    for k in sorted(rows):
        r = rows[k]
        if r["cell"] == "describe_undefined":
            continue
        print(f"{k:58s} {str(n_of(r['result'])):>6s}  {fmt(r['result'])}")
    for c in ("submitted", "corrected"):
        d = rows.get(f"{c}/describe_undefined")
        if d:
            x = d["result"]
            print(f"{c}: separation undefined on {x['n_separation_undefined']} of {x['n_scorable']} scorable "
                  f"units, {len(x['networks_undefined'])} networks, {x['n_states_undefined']} states; "
                  f"radius counts {x['radius_counts_undefined']}")

    reg = rows.get(REGISTERED)
    summary = {
        "experiment": "E16",
        "registered_statistic": "Kendall tau_b(GAC radius, AUC_frac) on pooled separation-undefined units, "
                                "submitted cohort, network-cluster bootstrap 95% CI",
        "registered": None if reg is None else {
            "point": reg["result"].get("point"), "ci": [reg["result"].get("ci_lo_2p5"), reg["result"].get("ci_hi_97p5")],
            "n": n_of(reg["result"]), "n_networks": reg["result"].get("n_networks"),
            "loo_tau_min": reg["result"].get("loo_tau_min"), "loo_tau_max": reg["result"].get("loo_tau_max"),
            "verdict": verdict(reg["result"]), "n_boot": reg["n_boot"],
        },
        "secondary": {k: {"point": rows[k]["result"].get("point"),
                          "ci": [rows[k]["result"].get("ci_lo_2p5"), rows[k]["result"].get("ci_hi_97p5")],
                          "n": n_of(rows[k]["result"]), "verdict": verdict(rows[k]["result"])}
                      for k in SECONDARY if k in rows},
        "cells": {k: r["result"] for k, r in rows.items()},
        "consistency_failures": failures,
        "published_interval_mismatches": ci_mismatch,
        "smoke": any(r.get("smoke") for r in rows.values()),
        "n_boot": sorted({r["n_boot"] for r in rows.values()}),
        "criterion": sorted({r["criterion"] for r in rows.values()}),
        "commit": sorted({r["commit"] for r in rows.values()}),
        "n_cells": len(rows),
    }
    if reg:
        print(f"\nREGISTERED {REGISTERED}: {fmt(reg['result'])}, n={n_of(reg['result'])} -> {verdict(reg['result'])}")
        if verdict(reg["result"]) != "interval above zero":
            print("  per PREREG: no ranking claim on these units; say only that the radius is defined there.")
    if summary["smoke"]:
        print("SMOKE output: intervals use few resamples and are not results.")
    if ci_mismatch:
        print(f"published intervals not reproduced within 0.02: {ci_mismatch}")
    if failures:
        print(f"CONSISTENCY FAILURES in {len(failures)} cells: {failures}")
    (out / "results.json").write_text(json.dumps(summary, indent=2, default=float))
    print(f"wrote {out / 'results.json'}")
    return 2 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
