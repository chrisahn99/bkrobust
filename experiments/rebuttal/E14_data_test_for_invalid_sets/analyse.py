"""E14 analysis: detection minus false-alarm rate of the coefficient-agreement test (PREREG.md).

    python experiments/rebuttal/E14_data_test_for_invalid_sets/analyse.py --out results/rebuttal/E14
"""

from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
REGISTERED_N = 10000
THRESHOLD = 0.1
BOOT = 2000
SEED = 20261106
STRATA = ("1", "2", "3+", "no_failure_found")


def load_units(out: Path) -> list[dict]:
    by_key: dict[str, dict] = {}
    for f in sorted(out.glob("units*.jsonl")):
        for line in f.open():
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            if rec.get("status") != "error" or rec["key"] not in by_key:
                by_key[rec["key"]] = rec
    return list(by_key.values())


def rates(units: list[dict], n: str, field: str) -> np.ndarray:
    vals = [u["summary"][n][field] for u in units if u["summary"].get(n, {}).get(field) is not None]
    return np.array(vals, dtype=float)


def mean_ci(a: np.ndarray, rng) -> dict:
    if a.size == 0:
        return {"n_units": 0, "rate": None, "ci": None}
    boot = a[rng.integers(a.size, size=(BOOT, a.size))].mean(axis=1)
    return {"n_units": int(a.size), "rate": float(a.mean()),
            "ci": [float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))]}


def delta_ci(bad: np.ndarray, good: np.ndarray, rng) -> dict:
    """D - F with a cluster bootstrap over units, each group resampled on its own."""
    if bad.size == 0 or good.size == 0:
        return {"delta": None, "ci": None}
    bd = bad[rng.integers(bad.size, size=(BOOT, bad.size))].mean(axis=1)
    gd = good[rng.integers(good.size, size=(BOOT, good.size))].mean(axis=1)
    diff = bd - gd
    return {"delta": float(bad.mean() - good.mean()),
            "ci": [float(np.quantile(diff, 0.025)), float(np.quantile(diff, 0.975))]}


def block(units: list[dict], n: str, field: str, rng) -> dict:
    bad = [u for u in units if u["z_valid_at_truth"] is False]
    good = [u for u in units if u["z_valid_at_truth"] is True]
    rb, rg = rates(bad, n, field), rates(good, n, field)
    return {"detection": mean_ci(rb, rng), "false_alarm": mean_ci(rg, rng), **delta_ci(rb, rg, rng)}


def fmt(x) -> str:
    return "  -  " if x is None else f"{x:.3f}"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, default=ROOT / "results/rebuttal/E14")
    args = ap.parse_args(argv)
    units = load_units(args.out)
    if not units:
        print(f"no units in {args.out}", file=sys.stderr)
        return 1
    configs = [json.loads(f.read_text()) for f in sorted(args.out.glob("config*.json"))]
    status = collections.Counter(u.get("status") for u in units)
    ok = [u for u in units if u.get("status") == "ok"]
    n_values = sorted({n for u in ok for n in u["summary"]}, key=int)
    reg_n = str(REGISTERED_N) if str(REGISTERED_N) in n_values else (n_values[-1] if n_values else None)
    rng = np.random.default_rng(SEED)

    res: dict = {"experiment": "E14", "status_counts": dict(status), "n_ok_units": len(ok),
                 "registered_n": reg_n, "registered_is_10000": reg_n == str(REGISTERED_N),
                 "hard_failure_units": sum(bool(u.get("hard_failures")) for u in units),
                 "released_truth_mismatches": sum(u.get("matches_released_truth") is False for u in units),
                 "commits": sorted({u.get("commit") for u in units}),
                 "criteria": sorted({u.get("criterion") for u in units}),
                 "smoke": sorted({c.get("smoke") for c in configs}),
                 "n_invalid_z": sum(u["z_valid_at_truth"] is False for u in ok),
                 "n_valid_z": sum(u["z_valid_at_truth"] is True for u in ok),
                 "by_n": {}}
    for n in n_values:
        entry = {"all": block(ok, n, "reject_rate", rng),
                 "pairwise": block(ok, n, "pairwise_reject_rate", rng),
                 "by_radius": {s: block([u for u in ok if u["radius_stratum"] == s], n, "reject_rate", rng)
                               for s in STRATA}}
        clean = [u for u in ok if u["all_valid_at_truth"]]
        mixed = [u for u in ok if u["z_valid_at_truth"] and not u["all_valid_at_truth"]]
        entry["size_all_sets_valid"] = mean_ci(rates(clean, n, "reject_rate"), rng)
        entry["alarm_z_valid_other_invalid"] = mean_ci(rates(mixed, n, "reject_rate"), rng)
        res["by_n"][n] = entry
    bad = [u for u in ok if u["z_valid_at_truth"] is False]
    res["population_ceiling_invalid_z"] = mean_ci(np.array([u["pop_detectable_rate"] for u in bad]), rng)
    res["population_ceiling_by_radius"] = {
        s: mean_ci(np.array([u["pop_detectable_rate"] for u in bad if u["radius_stratum"] == s]), rng)
        for s in STRATA}
    reg = res["by_n"].get(reg_n, {}).get("all", {}) if reg_n else {}
    res["registered"] = {"statistic": "Delta = detection - false alarm, joint Wald test, alpha 0.05",
                         "n": reg_n, "delta": reg.get("delta"), "ci": reg.get("ci"),
                         "against_authors": None if reg.get("delta") is None else reg["delta"] < THRESHOLD}
    (args.out / "results.json").write_text(json.dumps(res, indent=1))

    print(f"E14  units: {dict(status)}  hard-failure units: {res['hard_failure_units']}  "
          f"released-truth mismatches: {res['released_truth_mismatches']}")
    print(f"     committed set invalid at truth: {res['n_invalid_z']}, valid: {res['n_valid_z']}")
    print(f"{'n':>6} {'stratum':>17} {'units D/F':>10} {'detect':>7} {'f.alarm':>7} {'delta':>7}  95% CI")
    for n in n_values:
        rows = [("all", res["by_n"][n]["all"])] + list(res["by_n"][n]["by_radius"].items())
        for name, b in rows:
            ci = "" if b["ci"] is None else f"[{b['ci'][0]:.3f}, {b['ci'][1]:.3f}]"
            print(f"{n:>6} {name:>17} {b['detection']['n_units']:>4}/{b['false_alarm']['n_units']:<5} "
                  f"{fmt(b['detection']['rate']):>7} {fmt(b['false_alarm']['rate']):>7} {fmt(b['delta']):>7}  {ci}")
        s, m = res["by_n"][n]["size_all_sets_valid"], res["by_n"][n]["alarm_z_valid_other_invalid"]
        pw = res["by_n"][n]["pairwise"]
        print(f"{n:>6}   size (all sets valid at truth) {fmt(s['rate'])} on {s['n_units']}; "
              f"alarm (Z valid, another set invalid) {fmt(m['rate'])} on {m['n_units']}; "
              f"pairwise delta {fmt(pw['delta'])}")
    pc = res["population_ceiling_invalid_z"]
    print(f"population ceiling (invalid Z, sets disagree in population): {fmt(pc['rate'])} on {pc['n_units']}")
    r = res["registered"]
    ci = "" if r["ci"] is None else f"[{r['ci'][0]:.3f}, {r['ci'][1]:.3f}]"
    print(f"REGISTERED n={r['n']}: delta {fmt(r['delta'])} {ci}  "
          f"against the authors (delta < {THRESHOLD}): {r['against_authors']}")
    print(f"wrote {args.out / 'results.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
