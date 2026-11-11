"""E12 analysis: N_net and the censored fraction C (PREREG.md), then the audit and the timings.

    python experiments/rebuttal/E12_networks_above_500_nodes/analyse.py [--out results/rebuttal/E12]
"""

from __future__ import annotations

import argparse
import collections
import json
import statistics
import sys
from pathlib import Path

from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[3]
REGISTERED = ("link", "munin", "munin2", "munin3", "munin4")


def cp(k: int, n: int, a: float = 0.05) -> list:
    if n == 0:
        return [None, None]
    lo = 0.0 if k == 0 else float(beta.ppf(a / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - a / 2, k + 1, n - k))
    return [round(lo, 4), round(hi, 4)]


def group(arm: str) -> str:
    return arm if arm in ("D_RAND", "D_DEGEN", "A_TRUE", "A_ORACLE") else "D_LLM*"


def q(xs: list) -> dict:
    xs = [x for x in xs if x is not None]
    if not xs:
        return {"n": 0, "median": None, "max": None}
    return {"n": len(xs), "median": round(statistics.median(xs), 4), "max": round(max(xs), 4)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "results/rebuttal/E12"))
    args = ap.parse_args(argv)
    out = Path(args.out)
    rows: dict = {}
    for f in sorted(out.glob("rows*.jsonl")):
        for line in f.open():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            rows[r["key"]] = r  # the last write of a unit wins (an error row redone)
    rows = list(rows.values())
    if not rows:
        print(f"no rows in {out}")
        return 1
    smoke = any(r.get("smoke") for r in rows)
    nets = sorted({r["network"] for r in rows})
    cert = [r for r in rows if r.get("status") == "ok"]
    status = collections.Counter(r.get("status") for r in rows)
    n_cens = sum(r["censored"] for r in cert)
    exact_nets = {r["network"] for r in cert if not r["censored"]}
    n_net = len(exact_nets & set(REGISTERED))
    C = n_cens / len(cert) if cert else None
    against = (n_net < 2) or (C is not None and C > 0.5)
    cons = [r for r in rows if r.get("consistency")]

    per_net = {}
    for n in nets:
        rs = [r for r in rows if r["network"] == n]
        cs = [r for r in rs if r.get("status") == "ok"]
        hyb = collections.Counter((r.get("hybrid") or {}).get("outcome", "not_called") for r in cs)
        per_net[n] = {
            "nodes": rs[0].get("nodes"), "cpdag_undirected": rs[0].get("cpdag_undirected"),
            "units": len(rs), "certifiable": len(cs),
            "top_state": sum(r.get("hop_method") == "top_state" for r in cs),
            "finite_exact": sum(not r["censored"] and r["r_hop"] not in (-1, None) for r in cs),
            "censored": sum(r["censored"] for r in cs),
            "no_retractable": sum(r.get("hop_method") == "no_retractable_edges" for r in cs),
            "hybrid_outcomes": dict(hyb),
            "top_seconds": q([r.get("top_seconds") for r in cs]),
            "search_seconds": q([(r.get("search") or {}).get("seconds") for r in cs]),
            "hybrid_wall_seconds": q([(r.get("hybrid") or {}).get("wall_seconds") for r in cs]),
            "claim_censored": sum(r.get("r_claim_status") == "censored_time" for r in cs),
            "unit_seconds": q([r.get("unit_seconds") for r in rs]),
        }

    audit = {}
    for g in ("D_LLM*", "D_RAND", "D_DEGEN", "A_TRUE", "A_ORACLE"):
        cs = [r for r in cert if group(r["arm"]) == g and "z_valid_at_truth" in r]
        k = sum(r["z_valid_at_truth"] for r in cs)
        audit[g] = {
            "certifiable": len(cs), "z_valid_at_truth": k,
            "fraction": round(k / len(cs), 4) if cs else None, "ci95": cp(k, len(cs)),
            "bucket_B6_hop": dict(collections.Counter(r["bucket_B6_hop"] for r in cs)),
            "bucket_B7_claim": dict(collections.Counter(r["bucket_B7_claim"] for r in cs)),
        }
    radii = collections.Counter(
        ("inf" if r["r_hop"] == -1 else r["r_hop"]) if not r["censored"]
        else f">={r['hop_lower_bound']}" for r in cert)

    res = {
        "smoke": smoke, "criterion": sorted({r.get("criterion") for r in rows}),
        "commit": sorted({r.get("commit") for r in rows}),
        "units": len(rows), "status": dict(status), "certifiable": len(cert),
        "registered": {"N_net": n_net, "networks_with_exact_radius": sorted(exact_nets),
                       "C_censored_fraction": round(C, 4) if C is not None else None,
                       "C_ci95": cp(n_cens, len(cert)), "censored": n_cens,
                       "against_authors": against},
        "consistency_failures": len(cons),
        "consistency_examples": [{"key": r["key"], "what": r["consistency"]} for r in cons[:10]],
        "radius_distribution": {str(k): v for k, v in sorted(radii.items(), key=str)},
        "per_network": per_net, "audit_by_condition": audit,
    }
    (out / "results.json").write_text(json.dumps(res, indent=1))

    print(f"E12{' (SMOKE: not the experiment)' if smoke else ''}: {len(rows)} units, "
          f"{len(cert)} certifiable, criterion {res['criterion']}")
    print(f"REGISTERED  N_net = {n_net} of 5   C = {res['registered']['C_censored_fraction']} "
          f"[{res['registered']['C_ci95'][0]}, {res['registered']['C_ci95'][1]}] "
          f"({n_cens}/{len(cert)})   against the authors: {against}")
    print(f"consistency failures: {len(cons)}")
    print(f"{'network':10}{'nodes':>6}{'units':>6}{'cert':>6}{'top':>6}{'finite':>7}"
          f"{'cens':>6}{'top s med':>11}{'hyb s max':>11}")
    for n, v in per_net.items():
        print(f"{n:10}{v['nodes'] or 0:>6}{v['units']:>6}{v['certifiable']:>6}{v['top_state']:>6}"
              f"{v['finite_exact']:>7}{v['censored']:>6}"
              f"{v['top_seconds']['median'] or 0:>11.4f}{v['hybrid_wall_seconds']['max'] or 0:>11.2f}")
    print(f"{'condition':10}{'cert':>6}{'valid@truth':>13}  B6 buckets / B7 buckets")
    for g, v in audit.items():
        print(f"{g:10}{v['certifiable']:>6}{str(v['fraction']):>13}  {v['bucket_B6_hop']} / "
              f"{v['bucket_B7_claim']}")
    print(f"status: {dict(status)}")
    print(f"wrote {out / 'results.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
