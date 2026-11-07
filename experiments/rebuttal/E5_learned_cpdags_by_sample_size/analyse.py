#!/usr/bin/env python3
"""E5 analysis: invalidity at truth and the radius on learned CPDAGs, by n.

    .venv/bin/python experiments/rebuttal/E5_learned_cpdags_by_sample_size/analyse.py --out results/rebuttal/E5

Reads every rows.*.jsonl in --out, runs the consistency checks of PREREG.md, prints the tables and
writes results.json. Exit 2 when a consistency check fails, 0 otherwise.
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
PRIMARY_ARM, PRIMARY_METHOD = "A_ORACLE_R", "pc"
N_LO, N_HI = 200, 5000
B, BOOT_SEED = 2000, 20261106
INF = float("inf")


def load_rows(out: Path) -> list[dict]:
    rows: dict[str, dict] = {}
    for p in sorted(Path(out).glob("rows.*.jsonl")):
        for line in p.open():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            rows.setdefault(r["unit"], r)  # resumed shards never rewrite a unit; first copy wins
    return list(rows.values())


def interval(r, exact, lower):
    if r is None:
        return (1.0, INF)
    r = int(r)
    if int(exact):
        return (INF, INF) if r == -1 else (float(r), float(r))
    if r == -1:
        return (float(lower) if lower is not None else 1.0, INF)
    return (1.0, INF)


def cross_phase_failures(rows: list[dict]) -> list[str]:
    """Check 4: Phase B on the reference graph equals Phase A wherever the knowledge is the same."""
    a = {(r["network"], r["arm"], r["rep"], r["x"], r["y"]): r for r in rows if r["phase"] == "A"}
    ref = {r["network"]: (r["method"], r["n"], r["seed"]) for r in rows if r["phase"] == "A"}
    bad = []
    for r in rows:
        if r["phase"] != "B" or ref.get(r["network"]) != (r["method"], r["n"], r["seed"]):
            continue
        s = a.get((r["network"], r["arm"], r["rep"], r["x"], r["y"]))
        if s is None or s["k_sha256"] != r["k_sha256"]:
            continue  # A_ORACLE_R, or knowledge the transfer changed: not the same row
        for f in ("chat_source_sha256", "status", "z_size", "z_valid_at_truth", "r_claim", "r_claim_status"):
            if s.get(f) != r.get(f) and not (r.get("claim_wall_timeout") or s.get("claim_wall_timeout")):
                bad.append(f"{r['unit']}: {f} {r.get(f)!r} vs phase A {s.get(f)!r}")
        if r.get("status") == "ok" and not (r.get("hop_wall_timeout") or s.get("hop_wall_timeout")):
            i, j = (interval(x.get("r_hop"), x.get("hop_exact"), x.get("hop_lower_bound")) for x in (r, s))
            if max(i[0], j[0]) > min(i[1], j[1]):
                bad.append(f"{r['unit']}: gac radius {i} vs phase A {j}")
    return bad


def nkey(n) -> str:
    return "true" if n is None else str(n)


def per_graph(rows: list[dict], arm: str, method: str) -> dict:
    """(network, n) -> list of (invalid, certifiable, total) per graph."""
    g: dict = collections.defaultdict(lambda: [0, 0, 0])
    for r in rows:
        if r["phase"] != "B" or r["arm"] != arm or r["method"] != method:
            continue
        c = g[(r["network"], nkey(r["n"]), r["gid"])]
        c[2] += 1
        if r["status"] == "ok":
            c[1] += 1
            c[0] += int(r["z_valid_at_truth"] == 0)
    cells: dict = collections.defaultdict(list)
    for (net, n, _), c in sorted(g.items()):
        cells[(net, n)].append(c)
    return cells


def pooled(cells: dict, n: str, idx: dict | None = None) -> float:
    inv = cert = 0
    for (net, nn), gs in cells.items():
        if nn != n:
            continue
        pick = gs if idx is None else [gs[i] for i in idx[(net, nn)]]
        inv += sum(c[0] for c in pick)
        cert += sum(c[1] for c in pick)
    return inv / cert if cert else float("nan")


def boot_delta(cells: dict, lo: str, hi: str) -> dict:
    point_lo, point_hi = pooled(cells, lo), pooled(cells, hi)
    rng = np.random.default_rng(BOOT_SEED)
    keys = sorted(k for k in cells if k[1] in (lo, hi))
    draws, d_lo, d_hi = [], [], []
    for _ in range(B):
        idx = {k: rng.integers(len(cells[k]), size=len(cells[k])) for k in keys}
        a, b = pooled(cells, lo, idx), pooled(cells, hi, idx)
        d_lo.append(a)
        d_hi.append(b)
        draws.append(b - a)
    draws = np.array(draws)

    def ci(v):
        v = np.array(v)
        v = v[~np.isnan(v)]
        return [round(float(np.percentile(v, 2.5)), 4), round(float(np.percentile(v, 97.5)), 4)] if len(v) else [None, None]

    return {"share_lo": round(point_lo, 4), "share_hi": round(point_hi, 4),
            "delta": round(point_hi - point_lo, 4), "delta_ci95": ci(draws),
            "share_lo_ci95": ci(d_lo), "share_hi_ci95": ci(d_hi),
            "n_graphs_lo": sum(len(cells[k]) for k in keys if k[1] == lo),
            "n_graphs_hi": sum(len(cells[k]) for k in keys if k[1] == hi),
            "boot_nan": int(np.isnan(draws).sum())}


def table(rows: list[dict]) -> list[dict]:
    g: dict = collections.defaultdict(list)
    for r in rows:
        if r["phase"] == "B":
            g[(r["arm"], r["method"], nkey(r["n"]))].append(r)
    out = []
    order = {"200": 0, "1000": 1, "5000": 2, "true": 3}
    for (arm, method, n), rs in sorted(g.items(), key=lambda kv: (kv[0][0], kv[0][1], order.get(kv[0][2], 9))):
        ok = [r for r in rs if r["status"] == "ok"]
        inv = sum(1 for r in ok if r["z_valid_at_truth"] == 0)
        timeouts = sum(1 for r in ok if r.get("hop_wall_timeout") or r.get("claim_wall_timeout"))
        exact = [r for r in ok if r.get("hop_exact") == 1 and not r.get("hop_wall_timeout")]
        finite = [r["r_hop"] for r in exact if r["r_hop"] != -1]
        out.append({
            "arm": arm, "method": method, "n": n, "graphs": len({r["gid"] for r in rs}), "rows": len(rs),
            "certifiable": len(ok), "certifiable_share": round(len(ok) / len(rs), 4) if rs else None,
            "invalid_at_truth": inv, "invalid_share": round(inv / len(ok), 4) if ok else None,
            "invalid_share_all_rows": round(inv / len(rs), 4) if rs else None,
            "knowledge_inconsistent": sum(1 for r in rs if r["status"] == "knowledge_inconsistent"),
            "class_correct_share": round(np.mean([r["g_class_correct"] for r in rs]), 4),
            "radius_finite": len(finite), "radius_unbounded": len(exact) - len(finite),
            "radius_lower_only": len(ok) - len(exact) - timeouts,
            "radius_median_finite": float(np.median(finite)) if finite else None,
            "radius_eq_1": sum(1 for v in finite if v == 1),
            "B6_dangerous": sum(1 for r in ok if r.get("bucket_B6_hop") == "dangerous"),
            "B7_dangerous": sum(1 for r in ok if r.get("bucket_B7_claim") == "dangerous"),
            "wall_timeouts": timeouts,
            "timeout_flag": bool(ok) and timeouts / len(ok) > 0.05,
        })
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/rebuttal/E5")
    args = ap.parse_args(argv)
    out = Path(args.out)
    out = out if out.is_absolute() else REPO / out
    rows = load_rows(out)
    if not rows:
        print(f"no rows under {out}")
        return 1

    phase_a = [r for r in rows if r["phase"] == "A"]
    fails = {r["unit"]: r["consistency_fail"] for r in rows if r.get("consistency_fail")}
    cross = cross_phase_failures(rows)
    criteria = sorted({r.get("criterion") for r in rows})
    commits = sorted({r.get("repo_commit") for r in rows})
    a_ok = [r for r in phase_a if r["status"] == "ok"]
    checks = {
        "phase_a_rows": len(phase_a), "phase_a_certifiable": len(a_ok),
        "phase_a_disagreements": sum(1 for r in phase_a if r.get("consistency_fail")),
        "phase_a_gac_ne_backdoor_exact": sum(
            1 for r in a_ok if r.get("hop_exact") == 1 and r.get("bd_hop_exact") == 1 and r.get("r_hop") != r.get("bd_r_hop")),
        "phase_b_row_failures": sum(1 for r in rows if r["phase"] == "B" and r.get("consistency_fail")),
        "crosscheck_rows": sum(1 for r in rows if "xc_r_hop" in r),
        "cross_phase_failures": len(cross), "criteria": criteria, "repo_commits": commits,
    }
    print("consistency:", json.dumps(checks))
    for u, f in list(fails.items())[:15]:
        print(f"  FAIL {u}: {f}")
    for line in cross[:15]:
        print(f"  CROSS {line}")

    tab = table(rows)
    hdr = ("arm", "method", "n", "graphs", "rows", "certifiable", "invalid_at_truth", "invalid_share",
           "invalid_share_all_rows", "knowledge_inconsistent", "radius_finite", "radius_unbounded",
           "radius_lower_only", "B7_dangerous", "wall_timeouts")
    print("\n" + " | ".join(hdr))
    for t in tab:
        print(" | ".join(str(t[h]) + ("!" if h == "wall_timeouts" and t["timeout_flag"] else "") for h in hdr))

    deltas = {}
    for arm in sorted({r["arm"] for r in rows if r["phase"] == "B"}):
        for method in sorted({r["method"] for r in rows if r["phase"] == "B"} - {"true_cpdag"}):
            cells = per_graph(rows, arm, method)
            if any(k[1] == str(N_LO) for k in cells) and any(k[1] == str(N_HI) for k in cells):
                deltas[f"{arm}|{method}"] = boot_delta(cells, str(N_LO), str(N_HI))
    reg = deltas.get(f"{PRIMARY_ARM}|{PRIMARY_METHOD}")
    verdict = None
    if reg and reg["delta_ci95"][1] is not None:
        verdict = "falls_with_n" if reg["delta_ci95"][1] < 0 else "does_not_fall"
    anchor = [t for t in tab if t["n"] == "true" and t["arm"] == PRIMARY_ARM]
    print(f"\nregistered: {PRIMARY_ARM} {PRIMARY_METHOD}, invalid share n={N_HI} minus n={N_LO}: "
          f"{json.dumps(reg)} -> {verdict}")
    print(f"true-CPDAG anchor ({PRIMARY_ARM}): {[(t['certifiable'], t['invalid_share']) for t in anchor]}")
    for k, v in deltas.items():
        print(f"  delta {k}: {v['share_lo']} -> {v['share_hi']}, {v['delta']} {v['delta_ci95']}")

    per_net = collections.defaultdict(dict)
    for (net, n), gs in per_graph(rows, PRIMARY_ARM, PRIMARY_METHOD).items():
        inv, cert = sum(c[0] for c in gs), sum(c[1] for c in gs)
        per_net[net][n] = {"graphs": len(gs), "certifiable": cert, "invalid_share": round(inv / cert, 4) if cert else None}

    any_fail = bool(fails) or bool(cross)
    res = {"experiment": "E5_learned_cpdags_by_sample_size", "registered_statistic": {
        "arm": PRIMARY_ARM, "method": PRIMARY_METHOD, "n_lo": N_LO, "n_hi": N_HI,
        "bootstrap": {"draws": B, "seed": BOOT_SEED, "resample": "graphs within (network, n)"},
        **(reg or {}), "verdict": verdict}, "consistency": checks, "consistency_ok": not any_fail,
        "deltas": deltas, "table": tab, "primary_per_network": per_net,
        "failures": {"rows": dict(list(fails.items())[:200]), "cross_phase": cross[:200]}}
    (out / "results.json").write_text(json.dumps(res, indent=1, default=str))
    print(f"\nwrote {out / 'results.json'}")
    return 2 if any_fail else 0


if __name__ == "__main__":
    sys.exit(main())
