"""E13 analysis: the three registered counts (PREREG.md) with Clopper-Pearson intervals.

    python experiments/rebuttal/E13_committed_set_consistency_checks/analyse.py [--out results/rebuttal/E13]
"""

from __future__ import annotations

import argparse
import glob
import json
from collections import Counter
from pathlib import Path

from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[3]


def cp(k: int, n: int, level: float = 0.95) -> list:
    """Clopper-Pearson interval for k successes in n; exact, and honest at k = 0."""
    if n == 0:
        return [None, None]
    a = (1 - level) / 2
    lo = 0.0 if k == 0 else float(beta.ppf(a, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - a, k + 1, n - k))
    return [round(lo, 4), round(hi, 4)]


def share(k: int, n: int) -> dict:
    return {"k": k, "n": n, "share": None if n == 0 else round(k / n, 4), "ci95": cp(k, n)}


def load(out: Path) -> list[dict]:
    rows: dict[str, dict] = {}
    for f in sorted(glob.glob(str(out / "rows*.jsonl"))):
        for line in open(f):
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            if r.get("status") != "error" or r["key"] not in rows:
                rows[r["key"]] = r  # later complete lines replace earlier error lines
    return list(rows.values())


def part_a(rows: list[dict]) -> dict:
    res = {}
    for label, pop in (("real", [r for r in rows if r.get("real_condition")]),
                       ("scrambled", [r for r in rows if not r.get("real_condition")])):
        ok = [r for r in pop if r["status"] == "ok"]
        amen = [r for r in ok if r["amenable"]]
        res[label] = {
            "n_units": len(pop),
            "status": dict(Counter(r["status"] for r in pop)),
            "n_reproduced_none": len(ok),
            "A_amenable": share(len(amen), len(ok)),
            "amenable_with_valid_closed_form": sum(1 for r in amen if r.get("z_valid_at_g0")),
            "y_possible_descendant": sum(1 for r in ok if r["y_possible_descendant"]),
            "closed_form_verdicts": dict(Counter(r["closed_form_verdict"] for r in ok)),
            "closed_form_in_some_O_D": sum(1 for r in amen if r.get("z_in_some_O_D")),
        }
    return res


def part_b(rows: list[dict]) -> dict:
    ok = [r for r in rows if r["status"] == "ok"]
    out = {"n_rows": len(rows), "status": dict(Counter(r["status"] for r in rows))}
    for verdict in ("closed_form", "not_amenable", "y_not_possible_descendant"):
        v = [r for r in ok if r["verdict"] == verdict]
        comp = [r for r in v if r["b"]["enumeration"] == "complete" and r["b"]["n_extensions"] > 0]
        d = {"n": len(v), "n_enumerated": len(comp),
             "excluded": dict(Counter(r["b"]["enumeration"] for r in v
                                      if not (r["b"]["enumeration"] == "complete" and r["b"]["n_extensions"] > 0))),
             "extensions_agree": sum(1 for r in comp if r["b"]["n_distinct_O_D"] == 1)}
        if verdict == "closed_form":
            dis = [r for r in comp if not r["b"].get("agree")]
            d["B_disagree"] = share(len(dis), len(comp))
            d["B_split"] = {"extensions_disagree": sum(1 for r in dis if r["b"]["n_distinct_O_D"] > 1),
                            "agree_on_other_set": sum(1 for r in dis if r["b"]["n_distinct_O_D"] == 1)}
            d["validity_in_every_extension_checked"] = sum(
                1 for r in comp if r["b"].get("z_valid_in_every_extension") is not None)
            d["examples"] = [r["key"] for r in dis[:10]]
        out[verdict] = d
    return out


def part_c(rows: list[dict]) -> dict:
    ok = [r for r in rows if r["status"] == "ok" and "c" in r]
    inv = [r for r in ok if not r["c"]["valid_at_truth"]]
    checked = [r for r in inv if r["c"].get("violation") is not None]
    viol = [r for r in checked if r["c"]["violation"]]
    slack = Counter(r["c"].get("slack") for r in checked if r["c"].get("slack") is not None)
    cor = [r for r in ok if r.get("corollary", {}).get("match") is not None]
    return {
        "n_ok_rows": len(ok),
        "n_invalid_at_truth": len(inv),
        "bound_kind": dict(Counter(r["c"]["bound_kind"] for r in inv)),
        "C_violations": share(len(viol), len(checked)),
        "violation_examples": [r["key"] for r in viol[:10]],
        "slack_distribution": {str(k): v for k, v in sorted(slack.items())},
        "secondary_rval_by_enumeration": {
            "n_checked": len(cor),
            "mismatches": sum(1 for r in cor if not r["corollary"]["match"]),
            "by_kind": dict(Counter(r["c"]["bound_kind"] for r in cor)),
            "not_enumerated": dict(Counter(r.get("corollary", {}).get("enumeration") for r in ok
                                           if r.get("corollary", {}).get("match") is None)),
        },
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=ROOT / "results/rebuttal/E13")
    args = ap.parse_args()
    rows = load(args.out)
    a_rows = [r for r in rows if r["part"] == "a"]
    bc_rows = [r for r in rows if r["part"] == "bc"]
    configs = [json.loads(Path(f).read_text()) for f in sorted(glob.glob(str(args.out / "config*.json")))]
    res = {
        "experiment": "E13",
        "criterion": sorted({r["criterion"] for r in rows}),
        "commits": sorted({r["commit"] for r in rows}),
        "smoke": any(r.get("smoke") for r in rows),
        "n_errors": sum(1 for r in rows if r["status"] == "error"),
        "n_hard": sum(1 for r in rows if r.get("hard")),
        "hard_examples": [[r["key"], r["hard"]] for r in rows if r.get("hard")][:20],
        "configs": configs,
        "a": part_a(a_rows),
        "b": part_b(bc_rows),
        "c": part_c(bc_rows),
    }
    (args.out / "results.json").write_text(json.dumps(res, indent=1, default=str))

    def fmt(s):
        return f"{s['k']}/{s['n']}  [{s['ci95'][0]}, {s['ci95'][1]}]" if s["n"] else "0/0"

    a, b, c = res["a"], res["b"], res["c"]
    print(f"E13  criterion {res['criterion']}  commits {[x[:12] for x in res['commits']]}  smoke {res['smoke']}")
    print(f"errors {res['n_errors']}   hard consistency failures {res['n_hard']}")
    print("-" * 72)
    print(f"(a) real-name units           {a['real']['n_units']:>6}   reproduced 'none' {a['real']['n_reproduced_none']}")
    print(f"    A amenable at G0          {fmt(a['real']['A_amenable'])}   with valid closed form "
          f"{a['real']['amenable_with_valid_closed_form']}")
    print(f"    scrambled-name units      {a['scrambled']['n_units']:>6}   amenable {fmt(a['scrambled']['A_amenable'])}")
    cf = b["closed_form"]
    print(f"(b) closed_form rows          {cf['n']:>6}   enumerated {cf['n_enumerated']}   excluded {cf['excluded']}")
    print(f"    B disagree                {fmt(cf['B_disagree'])}   split {cf['B_split']}")
    for v in ("not_amenable", "y_not_possible_descendant"):
        print(f"    {v:<26}{b[v]['n']:>6}   enumerated {b[v]['n_enumerated']}   extensions agree {b[v]['extensions_agree']}")
    print(f"(c) invalid at truth          {c['n_invalid_at_truth']:>6}   bounds {c['bound_kind']}")
    print(f"    C violations              {fmt(c['C_violations'])}")
    print(f"    slack (reversed - bound)  {c['slack_distribution']}")
    s = c["secondary_rval_by_enumeration"]
    print(f"    r_val by Cor. 1(c) enum   checked {s['n_checked']}  mismatches {s['mismatches']}")
    print(f"wrote {args.out / 'results.json'}")


if __name__ == "__main__":
    main()
