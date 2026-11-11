"""Read-only checks named in FIDELITY.yaml. Exit 0 = pass. Run from the repository root.

    python experiments/rebuttal/E19_bias_radius_internal_standardisation/fidelity_checks.py <check> [--out DIR]
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]


def load(out: Path) -> list[dict]:
    rows: dict[str, dict] = {}
    for f in sorted(out.glob("units.shard*of*.jsonl")):
        for line in f.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                rows[r["instance_id"]] = r
    return list(rows.values())


def fails_with(rows: list[dict], prefix: str) -> list[str]:
    return [r["instance_id"] for r in rows for f in r.get("consistency_failures", []) if f.startswith(prefix)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("check")
    ap.add_argument("--out", default=str(ROOT / "results" / "rebuttal" / "E19"))
    ap.add_argument("--baseline-csv", default=str(ROOT / "results" / "epsilon" / "study" / "results.csv"))
    a = ap.parse_args()
    rows = load(Path(a.out))
    if not rows:
        print("no rows")
        return 1
    c = a.check
    if c == "covers_900":
        with open(a.baseline_csv, newline="") as fh:
            want = {r["instance_id"] for r in csv.DictReader(fh) if r["status"] == "ok"}
        got = {r["instance_id"] for r in rows}
        print(f"baseline ok {len(want)}, written {len(got)}, missing {len(want - got)}, extra {len(got - want)}")
        return 0 if want == got and len(want) == 900 else 1
    if c == "identity":
        bad = fails_with(rows, "identity_") + fails_with(rows, "rebuild")
        print(f"rows whose rebuilt instance differs from the stored one: {len(bad)}")
        return 0 if not bad else 1
    if c == "naive_reproduces":
        bad = fails_with(rows, "naive vs baseline")
        ok = sum(1 for r in rows if r["arms"].get("naive", {}).get("status") == "ok")
        print(f"naive ok {ok}/{len(rows)}, mismatches with the stored study {len(set(bad))}")
        return 0 if not bad and ok >= 0.95 * len(rows) else 1
    if c == "posthoc_equals_naive":
        bad = fails_with(rows, "posthoc vs naive")
        print(f"posthoc differs from naive on {len(set(bad))} rows")
        return 0 if not bad else 1
    if c == "internal_is_standardised":
        ia = [r["arms"]["internal"] for r in rows if r["arms"].get("internal", {}).get("status") == "ok"]
        na = [r["arms"]["naive"] for r in rows if r["arms"].get("naive", {}).get("status") == "ok"]
        worst = max((x["max_unit_var_error"] for x in ia), default=float("inf"))
        naive_dev = sum(x["max_unit_var_error"] for x in na) / max(1, len(na))
        # the naive arm must NOT be standardised, or the two arms are the same model
        print(f"internal max |var-1| = {worst:.3g} over {len(ia)}; naive mean max |var-1| = {naive_dev:.3g}")
        return 0 if ia and worst <= 1e-9 and naive_dev > 0.05 else 1
    if c == "same_r_val":
        bad = fails_with(rows, "r_val differs")
        print(f"rows whose r_val differs across arms: {len(bad)}")
        return 0 if not bad else 1
    if c == "strategies_agree":
        bad = [x for x in (fails_with(rows, "naive: incremental") + fails_with(rows, "posthoc: incremental")
                           + fails_with(rows, "internal: incremental"))]
        print(f"incremental/bisection disagreements: {len(bad)}")
        return 0 if not bad else 1
    if c == "criterion_and_commit":
        crit = {r.get("criterion") for r in rows}
        commits = {r.get("repo_commit") for r in rows}
        dirty = any(r.get("repo_code_dirty") for r in rows)
        print(f"criteria {crit}, commits {commits}, code dirty {dirty}")
        return 0 if crit == {"gac"} and len(commits) == 1 and "nogit" not in commits and not dirty else 1
    if c == "exclusions_small":
        bad = sum(1 for r in rows if not r.get("status_all_ok"))
        print(f"excluded (timeout or error in an arm) {bad}/{len(rows)}")
        return 0 if bad <= 0.05 * len(rows) else 1
    print(f"unknown check {c}")
    return 2


if __name__ == "__main__":
    sys.exit(main())
