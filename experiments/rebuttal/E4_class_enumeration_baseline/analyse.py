"""E4 analysis: when is listing [C^] faster than the hybrid? Registered statistic in PREREG.md.

    .venv/bin/python experiments/rebuttal/E4_class_enumeration_baseline/analyse.py [DIR]
    .venv/bin/python experiments/rebuttal/E4_class_enumeration_baseline/analyse.py DIR --check NAME   # FIDELITY.yaml

DIR defaults to results/rebuttal/E4. Writes DIR/results.json.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import random
import statistics
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
BINS = [(1, 1, "1"), (2, 10, "2-10"), (11, 100, "11-100"), (101, 1_000, "101-1e3"), (1_001, 10_000, "1e3-1e4"),
        (10_001, 100_000, "1e4-1e5"), (100_001, math.inf, ">1e5")]
FAMILIES = ("four", "highk", "timing")


def load(d: Path) -> list[dict]:
    rows: dict[str, dict] = {}
    for f in sorted(glob.glob(str(d / "rows_*.jsonl"))):
        for line in Path(f).read_text().splitlines():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            rows[r["unit_id"]] = r  # a unit written twice (resumed shard) counts once
    return list(rows.values())


def size_bin(r: dict) -> str:
    s = r.get("class_size")
    if s is None:
        return "unknown"
    return next(lab for lo, hi, lab in BINS if lo <= s <= hi)


def wilson(k: int, n: int, z: float = 1.959964) -> tuple[float | None, float | None]:
    if n == 0:
        return None, None
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, c - h), min(1.0, c + h)


def faster(r: dict, leg: str) -> bool:
    """A listing counts as faster only if it finished and beat the hybrid's wall time."""
    e = r.get(leg)
    return bool(e) and e["complete"] and e["seconds"] < r["hybrid"]["wall_seconds"]


def boot_median(xs: list[float], b: int = 2000, seed: int = 0) -> tuple[float | None, float | None, float | None]:
    if not xs:
        return None, None, None
    rng = random.Random(seed)
    meds = sorted(statistics.median(rng.choices(xs, k=len(xs))) for _ in range(b))
    return statistics.median(xs), meds[int(0.025 * b)], meds[int(0.975 * b) - 1]


def cell(rows: list[dict], leg: str) -> dict:
    n = len(rows)
    k = sum(faster(r, leg) for r in rows)
    lo, hi = wilson(k, n)
    ratios = [r[leg]["seconds"] / max(r["hybrid"]["wall_seconds"], 1e-9) for r in rows if r.get(leg)]
    med, mlo, mhi = boot_median(ratios)
    return {"n": n, "n_faster": k, "frac_faster": k / n if n else None, "wilson95": [lo, hi],
            "n_incomplete": sum(1 for r in rows if r.get(leg) and not r[leg]["complete"]),
            "median_ratio_to_hybrid": med, "median_ratio_boot95": [mlo, mhi],
            "max_class_size_faster": max((r["class_size"] for r in rows if faster(r, leg) and r.get("class_size")), default=None)}


def analyse(d: Path) -> dict:
    rows = load(d)
    ok = [r for r in rows if r.get("status") == "ok"]
    res: dict = {"dir": str(d), "n_rows": len(rows), "n_ok": len(ok),
                 "n_excluded": sum(r.get("status") == "excluded" for r in rows),
                 "n_error": sum(r.get("status") == "error" for r in rows),
                 "n_inconsistent": sum(not r.get("consistent", False) for r in rows),
                 "inconsistent_units": [r["unit_id"] for r in rows if not r.get("consistent", False)][:50],
                 "n_hybrid_inexact": sum(bool(r.get("flag_hybrid_inexact")) for r in ok),
                 "criteria": sorted({r.get("criterion") for r in rows}), "commits": sorted({r.get("commit") for r in rows}),
                 "families": {}, "pooled": {}}
    labels = [b[2] for b in BINS] + ["unknown"]
    for leg in ("enum", "bnb"):
        res["pooled"][leg] = {lab: cell([r for r in ok if size_bin(r) == lab], leg) for lab in labels}
        cross = next((lab for lab in labels[:-1] if res["pooled"][leg][lab]["n"] and res["pooled"][leg][lab]["frac_faster"] < 0.5), None)
        res["pooled"][leg]["crossover_bin"] = cross
    for fam in FAMILIES:
        fr = [r for r in ok if r["family"] == fam]
        res["families"][fam] = {
            leg: {"overall": cell(fr, leg), "by_bin": {lab: cell([r for r in fr if size_bin(r) == lab], leg) for lab in labels}}
            for leg in ("enum", "bnb")
        }
        res["families"][fam]["hybrid_methods"] = {m: sum(r["hybrid"]["method"] == m for r in fr)
                                                  for m in sorted({r["hybrid"]["method"] for r in fr})}
    return res


def fmt(x, nd=2):
    return "-" if x is None else f"{x:.{nd}f}"


def print_table(res: dict) -> None:
    print(f"E4 {res['dir']}: {res['n_ok']} ok, {res['n_excluded']} excluded, {res['n_error']} errors, "
          f"{res['n_inconsistent']} inconsistent, {res['n_hybrid_inexact']} hybrid inexact; criterion {res['criteria']}")
    print("\nRegistered: fraction of instances where listing [C^] (enum) finishes faster than the hybrid")
    print(f"{'family':8} {'bin':9} {'n':>5} {'enum':>6} {'95% CI':>13} {'ratio':>9} {'bnb':>6}")
    for fam in FAMILIES:
        f = res["families"][fam]
        if not f["enum"]["overall"]["n"]:
            continue
        for lab, c in list(f["enum"]["by_bin"].items()) + [("ALL", f["enum"]["overall"])]:
            if not c["n"]:
                continue
            b = f["bnb"]["overall"] if lab == "ALL" else f["bnb"]["by_bin"][lab]
            lo, hi = c["wilson95"]
            print(f"{fam:8} {lab:9} {c['n']:5d} {fmt(c['frac_faster']):>6} [{fmt(lo)}, {fmt(hi)}] "
                  f"{fmt(c['median_ratio_to_hybrid'], 1):>9} {fmt(b['frac_faster']):>6}")
    print(f"\nPooled crossover bin (lowest bin with fraction < 0.5): enum {res['pooled']['enum']['crossover_bin']}, "
          f"bnb {res['pooled']['bnb']['crossover_bin']}")
    print("ratio = median enum seconds / hybrid seconds (incomplete listings enter with their cut-off time, a lower bound)")


# --------------------------------------------------------------------------- fidelity checks


def check(name: str, d: Path) -> int:
    rows = load(d)
    ok = [r for r in rows if r.get("status") == "ok"]
    mans = [json.loads(Path(f).read_text()) for f in glob.glob(str(d / "manifest_*.json"))]
    if name == "criterion_gac":
        good = bool(rows) and all(r.get("criterion") == "gac" for r in rows)
    elif name == "commit_recorded":
        good = bool(rows) and all(r.get("commit") not in (None, "", "unknown") for r in rows)
    elif name == "no_inconsistency":
        good = bool(rows) and all(r.get("consistent") for r in rows)
    elif name == "all_units":
        want = sum(m["n_units_shard"] for m in mans)
        good = bool(mans) and len(rows) == want and mans[0]["n_shards"] == len({m["shard"] for m in mans})
    elif name == "families_complete":
        n = {f: sum(r["family"] == f for r in rows) for f in FAMILIES}
        four_n = next((m["n_four_node"] for m in mans if m.get("n_four_node")), None)
        good = n["highk"] == 144 and n["timing"] == 256 and four_n is not None and n["four"] == four_n
        print(n, "four-node selection:", four_n)
    elif name == "cap_1e5":
        good = bool(rows) and all(r.get("cap") == 100_000 for r in rows)
    elif name == "listing_is_class":
        # complete listings must have exactly the counted number of DAGs
        good = bool(ok) and all(r["enum"]["n_dags"] == r["class_size"] for r in ok if r["enum"]["complete"] and r.get("class_size"))
    elif name == "lib_cross_checked":
        good = any(r.get("lib", {}).get("same_dag_set") for r in ok) and all(r["lib"]["same_dag_set"] for r in ok if "lib" in r)
    elif name == "not_smoke":
        good = bool(rows) and not any(r.get("smoke") for r in rows)
    else:
        print(f"unknown check {name}")
        return 3
    print(f"{name}: {'pass' if good else 'FAIL'}")
    return 0 if good else 1


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("dir", nargs="?", default="results/rebuttal/E4")
    p.add_argument("--check", default=None)
    a = p.parse_args()
    d = Path(a.dir)
    if not d.is_absolute():
        d = REPO / d
    if a.check:
        return check(a.check, d)
    res = analyse(d)
    (d / "results.json").write_text(json.dumps(res, indent=1, default=str))
    print_table(res)
    print(f"\nwrote {d / 'results.json'}")
    return 1 if res["n_inconsistent"] else 0


if __name__ == "__main__":
    sys.exit(main())
