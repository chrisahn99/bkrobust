"""E1 analysis: how many of the released lower-bound rows the CPDAG check settles.

Reads every ``rows*.jsonl`` in the output folder (shards and resumed runs alike),
keeps one record per row, and writes ``results.json`` beside them.

    python experiments/rebuttal/E1_cpdag_check_on_lower_bound_rows/analyse.py --out results/rebuttal/E1
"""

from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[3]
N_RELEASED = 372


def clopper_pearson(k: int, n: int, level: float = 0.95) -> tuple[float | None, float | None]:
    if n == 0:
        return None, None
    a = (1 - level) / 2
    lo = 0.0 if k == 0 else float(beta.ppf(a, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - a, k + 1, n - k))
    return lo, hi


def load(out: Path) -> list[dict]:
    """One record per row; a completed record beats an error record for the same row."""
    best: dict[str, dict] = {}
    for f in sorted(out.glob("rows*.jsonl")):
        for line in f.open():
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            old = best.get(rec["key"])
            if old is None or old.get("status") == "error":
                best[rec["key"]] = rec
    return [best[k] for k in sorted(best)]


def arm_class(arm: str) -> str:
    return arm if arm in ("D_RAND", "A_TRUE", "A_ORACLE") else "deployment"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=ROOT / "results/rebuttal/E1")
    ap.add_argument("--results", type=Path, default=None, help="default <out>/results.json")
    args = ap.parse_args(argv)
    recs = load(args.out)
    if not recs:
        print(f"no rows in {args.out}")
        return 1
    ok = [r for r in recs if r.get("status") == "ok"]
    mism = [r for r in recs if r.get("status") == "rebuild_mismatch"]
    errs = [r for r in recs if r.get("status") == "error"]
    hard = [r for r in recs if r.get("hard_failures")]
    smoke = any(r.get("smoke") for r in recs)

    res: dict = {
        "experiment": "E1",
        "smoke": smoke,
        "criterion": "gac",
        "commits": sorted({r.get("commit", "unknown") for r in recs}),
        "n_released_population": N_RELEASED,
        "n_records": len(recs),
        "n_rebuilt": len(ok),
        "n_rebuild_mismatch": len(mism),
        "rebuild_mismatch_fields": dict(collections.Counter(f for r in mism for f in r.get("mismatch", []))),
        "n_error": len(errs),
        "n_hard_failure_rows": len(hard),
        "hard_failures": {r["key"]: r["hard_failures"] for r in hard},
    }
    n = len(ok)
    for crit in ("gac", "backdoor"):
        s = sum(r[crit]["valid_at_cpdag"] for r in ok)
        lo, hi = clopper_pearson(s, n)
        unsettled = [r for r in ok if not r[crit]["valid_at_cpdag"]]
        res[crit] = {
            "settled_infinite": s,
            "n": n,
            "proportion": s / n if n else None,
            "ci95_clopper_pearson": [lo, hi],
            "outcome_unsettled": dict(collections.Counter(r[crit]["outcome"] for r in unsettled)),
            "radius_found_values": dict(collections.Counter(
                r[crit]["large_search"]["radius"] for r in unsettled if r[crit]["outcome"] == "radius_found")),
            "best_lower_bound_values": dict(collections.Counter(
                r[crit]["best_lb"] for r in unsettled if r[crit]["best_lb"] is not None)),
            "reproduces_released_lb": sum(r[crit]["reproduces_released_lb"] for r in ok),
            "released_lb_overstated_by_formula": sum(r[crit]["released_lb_overstated"] for r in ok),
            "top_seconds_median": sorted(r[crit]["top_seconds"] for r in ok)[n // 2] if n else None,
        }
    res["criteria_disagree_at_cpdag"] = sum(r["criteria_disagree_at_cpdag"] for r in ok)

    by = collections.defaultdict(lambda: [0, 0, 0, 0])
    for r in ok:
        for key in (("network", r["network"]), ("arm", arm_class(r["arm"]))):
            c = by[key]
            c[0] += 1
            c[1] += r["gac"]["valid_at_cpdag"]
            c[2] += r["backdoor"]["valid_at_cpdag"]
            c[3] += r["criteria_disagree_at_cpdag"]
    res["by_stratum"] = {f"{k[0]}={k[1]}": dict(zip(("n", "gac_settled", "backdoor_settled", "disagree"), v))
                         for k, v in sorted(by.items())}

    target = args.results or args.out / "results.json"
    target.write_text(json.dumps(res, indent=1, default=str))

    g = res["gac"]
    print(f"E1{' (SMOKE)' if smoke else ''}: rows read {len(recs)} of {N_RELEASED}; rebuilt {n}, "
          f"mismatch {len(mism)}, error {len(errs)}, hard-failure rows {len(hard)}")
    if n:
        print(f"REGISTERED  settled as infinite under GAC: {g['settled_infinite']}/{n} = "
              f"{g['proportion']:.3f}  95% CI [{g['ci95_clopper_pearson'][0]:.3f}, "
              f"{g['ci95_clopper_pearson'][1]:.3f}]")
    print(f"{'stratum':<24}{'n':>5}{'gac':>6}{'bd':>6}{'disagree':>10}")
    for name, v in res["by_stratum"].items():
        print(f"{name:<24}{v['n']:>5}{v['gac_settled']:>6}{v['backdoor_settled']:>6}{v['disagree']:>10}")
    for crit in ("gac", "backdoor"):
        c = res[crit]
        print(f"{crit:<9} settled {c['settled_infinite']}/{n}; unsettled outcomes {c['outcome_unsettled']}; "
              f"release lb reproduced {c['reproduces_released_lb']}/{n}, "
              f"overstated by one {c['released_lb_overstated_by_formula']}")
    print(f"criteria disagree at the CPDAG: {res['criteria_disagree_at_cpdag']}")
    for k, v in res["hard_failures"].items():
        print(f"HARD {k}: {v}")
    print(f"wrote {target}")
    return 1 if hard else 0


if __name__ == "__main__":
    sys.exit(main())
