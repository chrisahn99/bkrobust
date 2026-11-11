"""E2 analysis: the registered statistic (PREREG.md) from the rows run.py wrote.

    .venv/bin/python experiments/rebuttal/E2_running_time_by_graph_size/analyse.py --out results/rebuttal/E2

Censored values (over the cap, out of memory, killed, or skipped by the per-series rule) are
+infinity: a median is finite only when more than half of the runs finished, and its interval
comes from order statistics, so no censored value is ever imputed.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
INF = math.inf
CENSORED = ("timeout", "mem", "over_cap", "skipped")
PRACTICAL_S = 60.0
K_BINS = ((1, 5), (6, 10), (11, 20), (21, 50), (51, 100), (101, 200), (201, 10**9))


def load(out: Path) -> list[dict]:
    last: dict[str, dict] = {}
    for f in sorted(out.glob("rows_*.jsonl")):
        for line in f.read_text().splitlines():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            last[r["uid"]] = r  # a retried unit replaces its earlier row
    return list(last.values())


def stage_value(row: dict, stage: str):
    """Seconds for one stage; INF if censored; None if the stage did not apply to this unit."""
    p, legs = row["prep"], row["legs"]
    if stage == "cpdag_build":
        return p["cpdag_s"]
    if stage == "top_check":
        return p["top_median_s"]
    if stage == "validity_call":
        return p["g0_call_median_s"]
    name = {"search_b3": "search", "sat_build": "sat", "sat_ladder": "sat",
            "hybrid": "hybrid", "search_full": "search_full"}[stage]
    leg = legs.get(name, {})
    st = leg.get("status")
    if st in (None, "not_needed"):
        return None
    if stage == "sat_ladder" and p["top_valid"]:
        return None  # the ladder alone runs only where the top check fails, as in the hybrid
    if stage == "sat_build":
        b = leg.get("build_seconds", leg.get("partial", {}).get("build_seconds"))
        return b if b is not None else (INF if st in CENSORED else None)
    if st in CENSORED:
        return INF
    if st != "ok":
        return None  # infrastructure error: excluded and counted
    if stage == "sat_ladder":
        return leg["ladder_total_seconds"]
    return leg["seconds"]


def median_ci(vals: list[float], level: float = 0.95):
    """Sample median and a distribution-free interval from order statistics (binomial(m, 1/2))."""
    v = sorted(vals)
    m = len(v)
    if m == 0:
        return None, None, None
    med = v[(m - 1) // 2] if m % 2 else (v[m // 2 - 1] + v[m // 2]) / 2
    pmf = [math.comb(m, i) / 2**m for i in range(m + 1)]
    best = (0, m - 1)
    for lo in range(m // 2 + 1):
        hi = m - 1 - lo
        if hi < lo:
            break
        if sum(pmf[lo + 1:hi + 1]) >= level:  # P(v[lo] <= median <= v[hi])
            best = (lo, hi)
    return med, v[best[0]], v[best[1]]


def summarise(vals: list[float]) -> dict:
    med, lo, hi = median_ci(vals)
    fin = [x for x in vals if x != INF]
    return {"m": len(vals), "censored": len(vals) - len(fin), "median_s": med,
            "ci95_lo_s": lo, "ci95_hi_s": hi, "max_s": max(vals) if vals else None,
            "max_finished_s": max(fin) if fin else None}


def fmt(x) -> str:
    if x is None:
        return "-"
    if x == INF:
        return ">cap"
    return f"{x * 1000:.3g}ms" if x < 1 else f"{x:.3g}s"


def clean(o):
    """Censored infinities as the string "censored", so results.json is standard JSON."""
    if isinstance(o, float) and o == INF:
        return "censored"
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    return o


STAGES = ("cpdag_build", "top_check", "validity_call", "search_b3", "sat_build", "sat_ladder",
          "hybrid", "search_full")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="results/rebuttal/E2", help="directory holding rows_*.jsonl")
    args = ap.parse_args()
    out = Path(args.out)
    if not out.is_absolute():
        out = ROOT / out
    rows = load(out)
    if not rows:
        print(f"no rows in {out}")
        return 1

    status = Counter(r.get("status") for r in rows)
    used = [r for r in rows if r.get("status") == "ok"]
    criteria = sorted({r.get("criterion") for r in rows})
    commits = sorted({r.get("commit") for r in rows})
    cons = [(r["uid"], r["consistency_failures"]) for r in used if r.get("consistency_failures")]
    errors = Counter((leg, v.get("status")) for r in used for leg, v in r["legs"].items()
                     if v.get("status") == "error")

    by_fn: dict[tuple, dict[str, list]] = defaultdict(lambda: defaultdict(list))
    by_k: dict[tuple, dict[str, list]] = defaultdict(lambda: defaultdict(list))
    methods: dict[tuple, Counter] = defaultdict(Counter)
    k_range: dict[tuple, list] = defaultdict(list)
    for r in used:
        key = (r["family"], r["n"])
        k = r["prep"]["k_undirected"]
        kb = next(b for b in K_BINS if b[0] <= k <= b[1])
        k_range[key].append(k)
        h = r["legs"].get("hybrid", {})
        methods[key][h.get("method") if h.get("status") == "ok" else h.get("status")] += 1
        for s in STAGES:
            v = stage_value(r, s)
            if v is not None:
                by_fn[key][s].append(v)
                by_k[(r["family"], kb)][s].append(v)

    table = []
    for key in sorted(by_fn):
        entry = {"family": key[0], "n": key[1], "k_min": min(k_range[key]),
                 "k_max": max(k_range[key]), "hybrid_methods": dict(methods[key])}
        for s in STAGES:
            entry[s] = summarise(by_fn[key][s])
        table.append(entry)
    ktable = []
    for key in sorted(by_k):
        entry = {"family": key[0], "k_bin": f"{key[1][0]}-{key[1][1] if key[1][1] < 10**9 else 'inf'}"}
        for s in STAGES:
            entry[s] = summarise(by_k[key][s])
        ktable.append(entry)

    # Headline: largest sparse n up to which every hybrid run (at that n and below) took <= 60 s.
    practical = None
    for e in sorted((e for e in table if e["family"] == "sparse"), key=lambda e: e["n"]):
        mx = e["hybrid"]["max_s"]
        if mx is None or mx > PRACTICAL_S:
            break
        practical = e["n"]
    small_slow = [e["n"] for e in table if e["n"] <= 48 and (e["hybrid"]["max_s"] or 0) > PRACTICAL_S]

    print(f"E2 running time  criterion={criteria}  commit={commits}  units={dict(status)}")
    cols = ("n", "k", "top_check", "validity_call", "search_b3", "sat_build", "hybrid med [95% CI]",
            "hybrid max", "cens", "methods")
    for fam in ("appD", "sparse"):
        print(f"\nfamily {fam}")
        print(" | ".join(cols))
        for e in (e for e in table if e["family"] == fam):
            h = e["hybrid"]
            print(" | ".join([
                str(e["n"]), f"{e['k_min']}-{e['k_max']}", fmt(e["top_check"]["median_s"]),
                fmt(e["validity_call"]["median_s"]),
                f"{fmt(e['search_b3']['median_s'])}/{fmt(e['search_b3']['max_s'])}",
                f"{fmt(e['sat_build']['median_s'])}/{fmt(e['sat_build']['max_s'])}",
                f"{fmt(h['median_s'])} [{fmt(h['ci95_lo_s'])}, {fmt(h['ci95_hi_s'])}]",
                fmt(h["max_s"]), f"{h['censored']}/{h['m']}",
                ",".join(f"{k}:{v}" for k, v in sorted(e["hybrid_methods"].items(), key=str)),
            ]))
    print("\nby undirected edges k (hybrid median/max, search_b3 median/max, sat_build median/max)")
    for e in ktable:
        print(f"{e['family']:6s} k={e['k_bin']:8s} hybrid {fmt(e['hybrid']['median_s'])}/"
              f"{fmt(e['hybrid']['max_s'])}  search {fmt(e['search_b3']['median_s'])}/"
              f"{fmt(e['search_b3']['max_s'])}  sat_build {fmt(e['sat_build']['median_s'])}/"
              f"{fmt(e['sat_build']['max_s'])}  (m={e['hybrid']['m']})")
    print(f"\npractical size (sparse, every hybrid run <= {PRACTICAL_S:g}s): n <= {practical}")
    print(f"n <= 48 with a hybrid run > {PRACTICAL_S:g}s: {small_slow or 'none'}")
    print(f"consistency failures: {len(cons)}  leg errors: {dict(errors) or 'none'}")

    res = {"criterion": criteria, "commits": commits, "unit_status": dict(status),
           "n_used": len(used), "practical_size_sparse_n": practical,
           "practical_threshold_s": PRACTICAL_S, "n_le_48_with_hybrid_over_60s": small_slow,
           "against_authors": bool(small_slow) or (practical is None or practical < 64),
           "consistency_failures": cons,
           "leg_errors": {f"{a}:{b}": c for (a, b), c in errors.items()},
           "by_family_n": table, "by_family_k": ktable}
    (out / "results.json").write_text(json.dumps(clean(res), indent=1))
    print(f"wrote {out / 'results.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
