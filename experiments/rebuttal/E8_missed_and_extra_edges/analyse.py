"""E8 analysis: the two registered statistics of PREREG.md, plus the secondary tables.

    .venv/bin/python experiments/rebuttal/E8_missed_and_extra_edges/analyse.py            # full run
    .venv/bin/python experiments/rebuttal/E8_missed_and_extra_edges/analyse.py --smoke    # smoke output

Reads E8_shard*.jsonl (or E8_smoke.jsonl) in --out, writes results.json (or results_smoke.json) there.
"""
import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
Z = 1.959963984540054
BAR = 0.002
X3_REFERENCE = dict(k=15, n=35501, rate=0.000423, source="pilot x3 skeleton-perturbation study, 2026-08-31")
USABLE = ("unchanged", "changed_valid", "silent_bias")
ENS = ("original", "licensed", "large")


def wilson(k, n):
    if n == 0:
        return [float("nan"), float("nan")]
    ph = k / n
    d = 1 + Z * Z / n
    c = (ph + Z * Z / (2 * n)) / d
    h = Z * np.sqrt(ph * (1 - ph) / n + Z * Z / (4 * n * n)) / d
    return [max(0.0, c - h), min(1.0, c + h)]


def cluster_boot(per_scm, B=4000, seed=20261106):
    """SCM-level cluster bootstrap of a pooled ratio: perturbations of one SCM share D and the query."""
    if not per_scm:
        return [float("nan"), float("nan")]
    arr = np.array(per_scm, dtype=float)
    idx = np.random.default_rng(seed).integers(0, len(arr), size=(B, len(arr)))
    ks, ns = arr[idx, 0].sum(axis=1), arr[idx, 1].sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        rates = np.where(ns > 0, ks / ns, np.nan)
    if np.all(np.isnan(rates)):
        return [float("nan"), float("nan")]
    lo, hi = np.nanpercentile(rates, [2.5, 97.5])
    return [float(lo), float(hi)]


def block(rows, field, arm, keep):
    """Silent-bias rate over the perturbations of `field` ('dels' or 'adds') that pass `keep`."""
    cnt, per_scm = Counter(), []
    for r in rows:
        k = n = 0
        for q in r[field]:
            if not keep(q):
                continue
            o = q[arm]["o"]
            cnt[o] += 1
            if o in USABLE:
                n += 1
                k += o == "silent_bias"
        if n:
            per_scm.append((k, n))
    n = sum(cnt[o] for o in USABLE)
    k = cnt["silent_bias"]
    return dict(k=k, n=n, rate=(k / n if n else float("nan")), wilson=wilson(k, n),
                cluster=cluster_boot(per_scm), n_scm=len(per_scm), outcomes=dict(cnt))


def load(out, smoke):
    files = [out / "E8_smoke.jsonl"] if smoke else sorted(out.glob("E8_shard*.jsonl"))
    rows, bad = {}, 0
    for f in files:
        if not f.exists():
            continue
        for line in open(f):
            try:
                r = json.loads(line)
            except Exception:
                bad += 1
                continue
            rows[r["key"]] = r  # a unit written twice (resharding) counts once
    return files, list(rows.values()), bad


def fmt(b):
    if b["n"] == 0:
        return "no mass"
    return (f"{b['k']:6d}/{b['n']:7d} = {b['rate']:.6f}  wilson[{b['wilson'][0]:.6f},{b['wilson'][1]:.6f}]"
            f"  cluster[{b['cluster'][0]:.6f},{b['cluster'][1]:.6f}]")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(REPO / "results/rebuttal/E8"))
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    out = Path(args.out)
    files, rows, bad = load(out, args.smoke)
    if not rows:
        sys.exit(f"no rows in {out}")
    ok = [r for r in rows if r.get("status") == "ok"]
    res = dict(experiment="E8", smoke=args.smoke, files=[str(f) for f in files], n_units=len(rows),
               unparsable_lines=bad, status=dict(Counter(r.get("status") for r in rows)),
               criteria=sorted({r.get("criterion") for r in rows}),
               commits=sorted({r.get("commit") for r in rows}))
    anomalies = Counter(a["kind"] for r in rows for a in r.get("anomalies", []))
    res["anomalies"] = dict(anomalies)
    res["oracle_disagreements"] = res["status"].get("oracle_disagree", 0)

    hop = {"hop0": lambda q: q["hop"] == 0, "hop1": lambda q: q["hop"] == 1,
           "hop2p": lambda q: q["hop"] >= 2, "all": lambda q: True}
    # ---- deletion arm ----------------------------------------------------------------------
    dels = {}
    for scope, rs in [("pooled", ok)] + [(e, [r for r in ok if r["ens"] == e]) for e in ENS]:
        dels[scope] = {arm: {h: block(rs, "dels", arm, f) for h, f in hop.items()} for arm in ("SK", "S0")}
    res["deletion"] = dels
    res["deletion_hop2p_unreachable"] = sum(1 for r in ok for q in r["dels"] if q["hop"] >= 1 << 20)

    # ---- registered statistic 1 ------------------------------------------------------------
    s1 = dels["pooled"]["SK"]["hop2p"]
    if s1["n"] == 0:
        v1 = "NO MASS"
    elif s1["rate"] >= BAR:
        v1 = "AGAINST US"
    elif s1["cluster"][1] < BAR:
        v1 = "SUPPORTED"
    else:
        v1 = "INCONCLUSIVE"
    res["registered_1"] = dict(statistic="silent-bias rate, deletions at hop >= 2, arm SK, pooled",
                               bar=BAR, **{k: s1[k] for k in ("k", "n", "rate", "wilson", "cluster", "n_scm")},
                               underpowered=s1["n"] < 5000, verdict=v1, x3_reference=X3_REFERENCE)

    # ---- addition arm ----------------------------------------------------------------------
    adds = {}
    for scope, rs in [("pooled", ok)] + [(e, [r for r in ok if r["ens"] == e]) for e in ENS]:
        adds[scope] = {arm: block(rs, "adds", arm, lambda q: True) for arm in ("SK", "S0")}
    for kind in ("forced", "free"):
        adds[f"pooled_{kind}"] = {arm: block(ok, "adds", arm, lambda q, k=kind: q["kind"] == k)
                                  for arm in ("SK", "S0")}
    for h, f in hop.items():
        adds[f"pooled_{h}"] = {"SK": block(ok, "adds", "SK", f)}
    z0 = Counter(q.get("Z0_valid_in_D_plus") for r in ok for q in r["adds"])
    res["addition"] = adds
    res["addition_Z0_valid_in_D_plus"] = {str(k): v for k, v in z0.items()}
    res["addition_tested"] = sum(r["n_add_tested"] for r in ok)
    res["addition_candidates"] = sum(r["n_add_candidates"] for r in ok)
    failures = [dict(key=r["key"], edge=q["e"], kind=q["kind"], hop=q["hop"])
                for r in ok for q in r["adds"] if q["SK"]["o"] == "silent_bias"]
    F = len(failures)
    res["registered_2"] = dict(statistic="additions, arm SK, O* defined and invalid in the true DAG",
                               F=F, n_usable=adds["pooled"]["SK"]["n"],
                               wilson_upper_rate=adds["pooled"]["SK"]["wilson"][1],
                               verdict=("AGAINST US" if F else "SUPPORTED (F = 0)"),
                               failures=failures[:50])

    # ---- print -----------------------------------------------------------------------------
    print(f"E8 {'SMOKE ' if args.smoke else ''}units={len(rows)} ok={len(ok)} "
          f"criteria={res['criteria']} commits={[c[:10] for c in res['commits'] if c]}")
    print("status:", res["status"])
    print(f"oracle disagreements: {res['oracle_disagreements']}   anomalies: {dict(anomalies) or 0}"
          f"   unparsable lines: {bad}")
    for arm in ("SK", "S0"):
        print(f"\nDELETION, arm {arm}      silent bias k/n = rate")
        for scope in ("pooled",) + ENS:
            for h in hop:
                print(f"  {scope:9} {h:5}  {fmt(dels[scope][arm][h])}")
    print(f"\nADDITION ({res['addition_tested']} tested of {res['addition_candidates']} candidates)")
    for scope in ("pooled",) + ENS + ("pooled_forced", "pooled_free"):
        for arm in adds[scope]:
            b = adds[scope][arm]
            print(f"  {scope:14} {arm}  {fmt(b)}  outcomes={b['outcomes']}")
    print(f"  committed set Z0 still valid in D plus the edge: {res['addition_Z0_valid_in_D_plus']}")
    r1, r2 = res["registered_1"], res["registered_2"]
    print("\nREGISTERED 1  silent-bias rate at hop >= 2, arm SK, pooled")
    print(f"  {r1['k']}/{r1['n']} = {r1['rate']:.6f}  cluster 95% [{r1['cluster'][0]:.6f}, "
          f"{r1['cluster'][1]:.6f}]  bar {BAR}  -> {v1}{'  (UNDERPOWERED: n < 5000)' if r1['underpowered'] else ''}")
    print(f"  X3 reference (post hoc, other seeds): {X3_REFERENCE['k']}/{X3_REFERENCE['n']} = {X3_REFERENCE['rate']}")
    print("REGISTERED 2  additions where O* is invalid in the true DAG, arm SK")
    print(f"  F = {F} of {r2['n_usable']} usable  -> {r2['verdict']}")
    if res["oracle_disagreements"] or anomalies:
        print("\nWARNING: oracle disagreements or theory anomalies present; PREREG says fix and redo before reading.")
    if args.smoke:
        print("\n(smoke output: the numbers above are a code-path check, not results)")

    dest = out / ("results_smoke.json" if args.smoke else "results.json")
    dest.write_text(json.dumps(res, indent=1) + "\n")
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
