"""E18 analysis: signatures, ceilings and the Table 2 contrasts per arm and draw.

    PYTHONPATH=src .venv/bin/python experiments/rebuttal/E18_section7_standardised_models/analyse.py [--smoke | --out DIR]

Reads every units.*.jsonl in the output folder, prints a table and writes results.json there.
The registered statistic (PREREG.md): distinct r_mu signatures over the 49 queries per
(arm, draw); primary = mean over draws of internal minus naive, with a 95 % t-interval.
"""

from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
EPS_KEYS = ["0.02", "0.1", "0.25", "0.4", "0.6", "0.9"]
ARMS = ["naive", "posthoc", "internal"]
INF = float("inf")
Q_DG25 = ("barley", "dg25", "s2225")
Q_AU = ("Schipf_2010", "A", "U")
Q_CDR = ("Kampen_2014", "AFF", "CDR")
STORED = {"signatures": 22, "signatures_with_ceiling": 32, "radius_values": 3}  # App. E.5


def load_units(out: Path) -> dict[tuple, dict[str, Any]]:
    units: dict[tuple, dict[str, Any]] = {}
    for path in sorted(out.glob("units.*.jsonl")):
        for line in path.open():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            units[(r["network"], r["x"], r["y"], r["draw"], r["arm"])] = r
    return units


def sig(r: dict[str, Any]) -> tuple:
    return tuple(r["r_mu_rel"][k] for k in EPS_KEYS)


def ceil3(r: dict[str, Any]) -> float:
    return round(r["b_top"], 3)


def as_dist(v: int) -> float:
    return INF if v == -1 else float(v)  # -1 is UNREACHED (UNR): never crosses


def contrasts(by_q: dict[tuple, dict[str, Any]]) -> dict[str, bool | None]:
    dg = by_q.get(Q_DG25)
    c1 = None
    if dg is not None:
        r = [as_dist(v) for v in sig(dg)]
        c1 = len(set(r[:5])) == 1 and r[0] < INF and r[5] > r[0]
    c2 = None
    if dg is not None and Q_AU in by_q and Q_CDR in by_q:
        base = as_dist(sig(dg)[5])
        c2 = as_dist(sig(by_q[Q_AU])[5]) > base and as_dist(sig(by_q[Q_CDR])[5]) > base
    return {"C1_dg25_crosses_2_to_60_at_once": c1, "C2_AU_and_CDR_reach_90_later": c2}


def t_interval(xs: list[float]) -> tuple[float | None, float | None, float | None]:
    if not xs:
        return None, None, None
    m = statistics.fmean(xs)
    if len(xs) < 2:
        return m, None, None
    from scipy.stats import t
    half = t.ppf(0.975, len(xs) - 1) * statistics.stdev(xs) / len(xs) ** 0.5
    return m, m - half, m + half


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None, help="default results/rebuttal/E18 (smoke: .../E18/smoke)")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    out = Path(args.out or ("results/rebuttal/E18/smoke" if args.smoke else "results/rebuttal/E18"))
    out = out if out.is_absolute() else ROOT / out
    units = load_units(out)
    if not units:
        raise SystemExit(f"no units in {out}")
    failures = sum(1 for p in out.glob("failures.*.jsonl") for _ in p.open())
    draws = sorted({k[3] for k in units})
    queries = sorted({k[:3] for k in units})
    status: dict[str, int] = {}
    for r in units.values():
        status[r["status"]] = status.get(r["status"], 0) + 1
    excluded = sorted({f"{r['network']} {r['x']}->{r['y']} draw {r['draw']} {r['arm']}: {r['status']}"
                       for r in units.values() if r["status"] != "ok"})

    per: dict[str, dict[int, dict[str, Any]]] = {a: {} for a in ARMS}
    ok = {k: r for k, r in units.items() if r["status"] == "ok"}
    for arm in ARMS:
        for d in draws:
            by_q = {k[:3]: r for k, r in ok.items() if k[3] == d and k[4] == arm}
            if not by_q:
                continue
            rows = list(by_q.values())
            naive = {k[:3]: r for k, r in ok.items() if k[3] == d and k[4] == "naive"}
            naive0 = {k[:3]: r for k, r in ok.items() if k[3] == 0 and k[4] == arm}
            same_naive = [sig(r) == sig(naive[q]) for q, r in by_q.items() if q in naive]
            same_draw0 = [sig(r) == sig(naive0[q]) for q, r in by_q.items() if q in naive0]
            vs = [r["sortability"]["var_sortability"] for r in rows if r["sortability"]["var_sortability"] is not None]
            r2 = [r["sortability"]["r2_sortability"] for r in rows if r["sortability"]["r2_sortability"] is not None]
            per[arm][d] = {
                "n_queries": len(rows),
                "signatures": len({sig(r) for r in rows}),
                "signatures_with_ceiling": len({(sig(r), ceil3(r)) for r in rows}),
                "ceiling_values": len({ceil3(r) for r in rows}),
                "radius_and_ceiling": len({(r["radius_gac"], ceil3(r)) for r in rows}),
                "radius_values": len({r["radius_gac"] for r in rows}),
                "frac_signature_as_naive_same_draw": statistics.fmean(same_naive) if same_naive else None,
                "frac_signature_as_draw0_same_arm": statistics.fmean(same_draw0) if same_draw0 else None,
                "mean_var_sortability": statistics.fmean(vs) if vs else None,
                "mean_r2_sortability": statistics.fmean(r2) if r2 else None,
                **contrasts(by_q),
            }

    def paired(key: str) -> dict[str, Any]:
        ds = [d for d in draws if d in per["internal"] and d in per["naive"]]
        diffs = [per["internal"][d][key] - per["naive"][d][key] for d in ds]
        m, lo, hi = t_interval(diffs)
        return {"per_draw": diffs, "mean": m, "ci95": [lo, hi], "n_draws": len(ds)}

    def med(arm: str, key: str) -> float | None:
        vals = [v[key] for v in per[arm].values()]
        return statistics.median(vals) if vals else None

    def frac(arm: str, key: str) -> float | None:
        vals = [v[key] for v in per[arm].values() if v[key] is not None]
        return statistics.fmean(vals) if vals else None

    collapse = None
    if med("internal", "signatures") is not None and med("naive", "signatures") is not None:
        collapse = med("internal", "signatures") <= 0.5 * med("naive", "signatures")
    posthoc_units = [r for k, r in ok.items() if k[4] == "posthoc"]
    posthoc_same = sum(1 for r in posthoc_units if r.get("checks", {}).get("posthoc_match"))
    repro = [r["checks"].get("repro_stored") for k, r in ok.items() if k[3] == 0 and k[4] == "naive"]
    naive0 = per["naive"].get(0, {})
    results = {
        "experiment": "E18", "out": str(out), "n_units": len(units), "status": status,
        "excluded": excluded, "consistency_failures": failures, "draws": draws,
        "n_queries": len(queries),
        "criterion": sorted({r.get("criterion") for r in units.values()}),
        "commit": sorted({r.get("commit") for r in units.values()}),
        "registered": {
            "statistic": "distinct r_mu signatures over the queries, internal minus naive, per draw",
            "signatures": paired("signatures"),
            "signatures_with_ceiling": paired("signatures_with_ceiling"),
            "median_signatures": {a: med(a, "signatures") for a in ARMS},
            "median_signatures_with_ceiling": {a: med(a, "signatures_with_ceiling") for a in ARMS},
            "collapse_under_internal": collapse,
            "frac_draws_C1": {a: frac(a, "C1_dg25_crosses_2_to_60_at_once") for a in ARMS},
            "frac_draws_C2": {a: frac(a, "C2_AU_and_CDR_reach_90_later") for a in ARMS},
        },
        "posthoc_units_matching_naive": [posthoc_same, len(posthoc_units)],
        "naive_draw0_reproduces_stored_rows": [sum(1 for v in repro if v), len(repro)],
        "naive_draw0_vs_appendix_E5": {k: [naive0.get(k), v] for k, v in STORED.items()},
        "per_arm_draw": per,
    }
    (out / "results.json").write_text(json.dumps(results, indent=2) + "\n")

    cols = ["n_queries", "signatures", "signatures_with_ceiling", "ceiling_values", "radius_values",
            "frac_signature_as_naive_same_draw", "C1_dg25_crosses_2_to_60_at_once", "C2_AU_and_CDR_reach_90_later"]
    short = ["n", "sig", "sig+ceil", "ceil", "r_val", "=naive", "C1", "C2"]
    print(f"E18 {out}  units={len(units)} status={status} consistency_failures={failures}")
    print(f"{'arm':9}{'draw':>5} " + " ".join(f"{s:>8}" for s in short))
    for arm in ARMS:
        for d in draws:
            if d not in per[arm]:
                continue
            v = per[arm][d]
            cells = []
            for c in cols:
                x = v[c]
                cells.append(f"{x:8.2f}" if isinstance(x, float) else f"{str(x):>8}")
            print(f"{arm:9}{d:>5} " + " ".join(cells))
    reg = results["registered"]
    s = reg["signatures"]
    print(f"\nregistered: internal - naive distinct signatures, mean {s['mean']} 95% CI {s['ci95']} "
          f"over {s['n_draws']} draws; medians {reg['median_signatures']}; collapse={collapse}")
    print(f"posthoc units matching naive: {posthoc_same}/{len(posthoc_units)}; "
          f"naive draw 0 reproduces stored rows: {results['naive_draw0_reproduces_stored_rows']}; "
          f"vs App. E.5 [this run, stored]: {results['naive_draw0_vs_appendix_E5']}")
    if excluded:
        print(f"excluded units: {len(excluded)} (listed in results.json)")
    print(f"wrote {out / 'results.json'}")


if __name__ == "__main__":
    main()
