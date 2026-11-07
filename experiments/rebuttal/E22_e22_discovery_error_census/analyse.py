#!/usr/bin/env python3
"""E22 analysis: the registered statistic, the census tables and results.json.

    .venv/bin/python experiments/rebuttal/E22_e22_discovery_error_census/analyse.py --out results/rebuttal/E22

Exit 2 when a row carries a consistency failure or the population is incomplete (full run only), 0 otherwise.
"""

from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
B, BOOT_SEED = 2000, 20261106
LOCAL_LIMIT = 0.20
ERROR_TYPES = ("query_direction_wrong", "missing_adjacency", "reversed_compelled", "false_v_structure",
               "missed_v_structure", "extra_adjacency", "under_oriented", "over_oriented")
CONDITIONS = ("y_not_true_descendant", "forbidden_in_z", "open_noncausal_path")


def load_rows(out: Path) -> list[dict]:
    rows: dict[str, dict] = {}
    for p in sorted(Path(out).glob("rows.*.jsonl")):
        for line in p.open():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            rows.setdefault(r["unit"], r)
    return list(rows.values())


def cluster_boot(rows: list[dict], pred, n_boot: int = B, seed: int = BOOT_SEED) -> dict:
    """Share of rows with pred, and a percentile interval resampling learned graphs."""
    by = collections.defaultdict(list)
    for r in rows:
        by[r["gid"]].append(1 if pred(r) else 0)
    gids = sorted(by)
    num = np.array([sum(by[g]) for g in gids], float)
    den = np.array([len(by[g]) for g in gids], float)
    k, n = int(num.sum()), int(den.sum())
    if n == 0:
        return {"count": 0, "n": 0, "share": None, "lo": None, "hi": None}
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(gids), size=(n_boot, len(gids)))
    s = num[idx].sum(1) / den[idx].sum(1)
    return {"count": k, "n": n, "share": k / n, "lo": float(np.percentile(s, 2.5)),
            "hi": float(np.percentile(s, 97.5)), "n_graphs": len(gids)}


def fmt(c: dict) -> str:
    if c["share"] is None:
        return "  n/a"
    return f"{c['count']:4d}/{c['n']:<4d} {100 * c['share']:5.1f} % [{100 * c['lo']:5.1f}, {100 * c['hi']:5.1f}]"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="results/rebuttal/E22")
    ap.add_argument("--n-boot", type=int, default=B)
    args = ap.parse_args(argv)
    out = Path(args.out)
    out = out if out.is_absolute() else REPO / out
    rows = load_rows(out)
    cfgs = [json.loads(p.read_text()) for p in sorted(out.glob("config.*.json"))]
    smoke = any(c.get("smoke") for c in cfgs)
    pop = cfgs[0]["n_rows_population"] if cfgs else None
    bad = [(r["unit"], r["consistency_fail"]) for r in rows if r.get("consistency_fail")]
    scored = [r for r in rows if "errors" in r]
    inv = [r for r in scored if r["z_valid_at_truth"] == 0]
    val = [r for r in scored if r["z_valid_at_truth"] == 1]
    nb = args.n_boot

    reg = cluster_boot(inv, lambda r: r["s3_hypothesis"] == 1, nb)
    lit = cluster_boot(inv, lambda r: r["census_only_extra_in_A"] == 1, nb)
    nonlocal_ = cluster_boot(inv, lambda r: r["any_error_within2"] == 0, nb)
    by_cond = {c: cluster_boot(inv, lambda r, c=c: r["first_failing"] == c, nb) for c in CONDITIONS}
    by_near = {t: cluster_boot(inv, lambda r, t=t: r["nearest_type"] == t, nb) for t in (*ERROR_TYPES, "none")}
    near2 = {t: {"invalid": cluster_boot(inv, lambda r, t=t: r["n_errors_within2"][t] > 0, nb),
                 "valid": cluster_boot(val, lambda r, t=t: r["n_errors_within2"][t] > 0, nb)} for t in ERROR_TYPES}
    s3_valid = cluster_boot(val, lambda r: r["s3_hypothesis"] == 1, nb)
    per_net = {}
    for net in sorted({r["network"] for r in scored}):
        iv = [r for r in inv if r["network"] == net]
        per_net[net] = {"n_valid": sum(r["network"] == net for r in val), "n_invalid": len(iv),
                        "registered_count": sum(r["s3_hypothesis"] == 1 for r in iv),
                        "no_error_within2": sum(r["any_error_within2"] == 0 for r in iv),
                        "first_failing": dict(collections.Counter(r["first_failing"] for r in iv)),
                        "nearest_type": dict(collections.Counter(r["nearest_type"] for r in iv))}
    # why the hypothesis fails on invalid rows
    why = collections.Counter(
        "missing_in_A" if r["s3_missing_in_A"] else "reversed_in_A" if r["s3_reversed_in_A"]
        else "closure_fails" if r["s3_closure_ok"] == 0 else "holds" for r in inv)

    against = []
    if reg["count"] > 0:
        against.append(f"{reg['count']} invalid rows satisfy the local S3 hypothesis (predicted 0)")
    if nonlocal_["share"] is not None and nonlocal_["share"] > LOCAL_LIMIT:
        against.append(f"{100 * nonlocal_['share']:.1f} % of invalid rows have no error within two hops (limit 20 %)")
    incomplete = (not smoke) and pop is not None and len(rows) != pop

    lines = [f"E22 {'SMOKE ' if smoke else ''}rows {len(rows)} (population {pop}), scored {len(scored)}: "
             f"{len(inv)} invalid, {len(val)} valid, on {len({r['gid'] for r in scored})} learned graphs",
             f"consistency failures: {len(bad)}" + ("; INCOMPLETE population" if incomplete else ""),
             "",
             "REGISTERED: invalid rows whose differences from the truth inside An({x,y} u Z) are only",
             f"  orientation-compatible extra adjacencies (S3 predicts 0): {fmt(reg)}",
             f"  same hypothesis among valid rows (not registered): {fmt(s3_valid)}",
             f"  literal census variant (census errors in A all extra adjacencies): {fmt(lit)}",
             f"  why the hypothesis fails on invalid rows: {dict(why)}",
             f"against the authors: {'; '.join(against) if against else 'no'}",
             "",
             f"invalid rows with no error of any kind within two hops of x or y: {fmt(nonlocal_)}  (against us above 20 %)",
             "",
             "first failing condition at the truth (invalid rows)"]
    lines += [f"  {c:24s} {fmt(by_cond[c])}" for c in CONDITIONS]
    lines += ["", "nearest error type (invalid rows; ties by the order below)"]
    lines += [f"  {t:24s} {fmt(by_near[t])}" for t in (*ERROR_TYPES, "none")]
    lines += ["", "error type within two hops of x or y: invalid rows | valid rows"]
    lines += [f"  {t:24s} {fmt(near2[t]['invalid'])} | {fmt(near2[t]['valid'])}" for t in ERROR_TYPES]
    lines += ["", "per network: n_invalid, n_valid, registered count, no error within two hops"]
    lines += [f"  {n:12s} {v['n_invalid']:4d} {v['n_valid']:4d} {v['registered_count']:4d} {v['no_error_within2']:4d}"
              for n, v in per_net.items()]
    print("\n".join(lines))

    res = {"experiment": "E22_e22_discovery_error_census", "smoke": smoke, "n_rows": len(rows),
           "n_population": pop, "n_scored": len(scored), "n_invalid": len(inv), "n_valid": len(val),
           "consistency_ok": not bad and not incomplete, "n_consistency_failures": len(bad),
           "registered_statistic": {"name": "invalid rows satisfying the local S3 hypothesis", "prediction": 0, **reg},
           "registered_statistic_valid_rows": s3_valid, "literal_census_variant": lit,
           "hypothesis_failure_reasons_invalid": dict(why),
           "no_error_within_two_hops": {**nonlocal_, "limit": LOCAL_LIMIT},
           "first_failing_condition": by_cond, "nearest_error_type": by_near, "error_within_two_hops": near2,
           "per_network": per_net, "against_authors": against, "bootstrap": {"draws": nb, "seed": BOOT_SEED,
                                                                           "cluster": "learned graph"},
           "repo_commit": sorted({r.get("repo_commit") for r in rows if r.get("repo_commit")})}
    (out / "results.json").write_text(json.dumps(res, indent=1))
    if bad:
        for u, f in bad[:10]:
            print(f"FAIL {u}: {f}")
    return 2 if (bad or incomplete) else 0


if __name__ == "__main__":
    sys.exit(main())
