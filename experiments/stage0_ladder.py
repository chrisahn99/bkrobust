"""Stage 0: the free hop-count rule against the committed radius.

Reads ``results/axisa3/instances.jsonl`` and nothing else. For every
admissible row the rule predicts the separation where it is measured and 1
where it is undefined, and the direction of every miss is recorded: a
prediction above the committed radius over-certifies, one below it is
conservative. The script also checks the direction of the session-5 law on
the measurable rows, the pairing of instances across coverage levels, and the
intraclass correlation of the ``r = 1`` indicator by network.

Writes ``results/stage0/ladder_rows.csv`` and ``results/stage0/ladder_summary.json``.

    python experiments/stage0_ladder.py
"""

from __future__ import annotations

import collections
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INSTANCES = ROOT / "results" / "axisa3" / "instances.jsonl"
OUT = ROOT / "results" / "stage0"

ILLUSTRATION = {"confounding", "mediator", "paths", "M-bias"}


def free_rule(row: dict) -> int:
    """Separation where measured, 1 where undefined."""
    if row["separation_status"] == "measured":
        return int(row["separation"])
    return 1


def icc_one_way(groups: dict[str, list[float]]) -> tuple[float, float]:
    """One-way ANOVA intraclass correlation and the m0 group-size constant."""
    k = len(groups)
    n = sum(len(g) for g in groups.values())
    grand = sum(sum(g) for g in groups.values()) / n
    msb = sum(len(g) * (sum(g) / len(g) - grand) ** 2 for g in groups.values()) / (k - 1)
    msw = sum(sum((x - sum(g) / len(g)) ** 2 for x in g) for g in groups.values()) / (n - k)
    m0 = (n - sum(len(g) ** 2 for g in groups.values()) / n) / (k - 1)
    return (msb - msw) / (msb + (m0 - 1) * msw), m0


def main() -> None:
    """Score the free rule against every admissible row and write the summary."""
    rows = [json.loads(line) for line in INSTANCES.open()]
    adm = [r for r in rows if r.get("admissible")]
    OUT.mkdir(parents=True, exist_ok=True)

    with (OUT / "ladder_rows.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "network",
                "x",
                "y",
                "coverage",
                "radius",
                "separation",
                "separation_status",
                "k_g0",
                "method",
                "free_rule",
                "diff",
            ]
        )
        for r in adm:
            b = free_rule(r)
            w.writerow(
                [
                    r["network"],
                    r["x"],
                    r["y"],
                    r["coverage"],
                    r["radius"],
                    r["separation"],
                    r["separation_status"],
                    r["k_g0"],
                    r["method"],
                    b,
                    b - r["radius"],
                ]
            )

    over = [r for r in adm if free_rule(r) > r["radius"]]
    under = [r for r in adm if free_rule(r) < r["radius"]]
    undefined = [r for r in adm if r["separation_status"] != "measured"]
    measurable = [r for r in adm if r["separation_status"] == "measured"]
    outside = [r for r in adm if r["network"] not in ILLUSTRATION]

    def law(r: dict) -> int:
        return min(int(r["separation"]), int(r["k_g0"]))

    pairs: dict[tuple, dict] = collections.defaultdict(dict)
    for r in adm:
        pairs[(r["network"], r["x"], r["y"])][r["coverage"]] = r["radius"]
    three = [v for v in pairs.values() if len(v) == 3]

    icc, m0 = icc_one_way(
        {
            n: [1.0 if r["radius"] == 1 else 0.0 for r in adm if r["network"] == n]
            for n in {r["network"] for r in adm}
        }
    )

    summary = {
        "rows": len(adm),
        "distinct_pairs": len(pairs),
        "pairs_admissible_at_all_three_coverages": len(three),
        "paired_radius_across_coverage": {
            "unchanged": sum(1 for v in three if len(set(v.values())) == 1),
            "lower_at_lower_coverage": sum(1 for v in three if v[0.25] < v[1.0]),
            "higher_at_lower_coverage": sum(1 for v in three if v[0.25] > v[1.0]),
        },
        "free_rule": {
            "exact": sum(1 for r in adm if free_rule(r) == r["radius"]),
            "within_one": sum(1 for r in adm if abs(free_rule(r) - r["radius"]) <= 1),
            "over_certifies": len(over),
            "over_certifies_by_network": dict(collections.Counter(r["network"] for r in over)),
            "under_certifies": len(under),
            "under_certifies_by_network": dict(
                collections.Counter(r["network"] for r in under).most_common()
            ),
            "under_certifies_in_undefined_stratum": sum(
                1 for r in under if r["separation_status"] != "measured"
            ),
            "abs_diff_ge_2": dict(
                collections.Counter(
                    r["network"] for r in adm if abs(free_rule(r) - r["radius"]) >= 2
                )
            ),
        },
        "undefined_stratum": {
            "rows": len(undefined),
            "radius_distribution": dict(
                sorted(collections.Counter(r["radius"] for r in undefined).items())
            ),
        },
        "law_on_measurable_rows": {
            "rows": len(measurable),
            "radius_equal_min": sum(1 for r in measurable if r["radius"] == law(r)),
            "radius_greater_than_min": sum(1 for r in measurable if r["radius"] > law(r)),
            "radius_less_than_min": sum(1 for r in measurable if r["radius"] < law(r)),
            "radius_greater_than_k_g0_anywhere": sum(1 for r in adm if r["radius"] > r["k_g0"]),
            "radius_equal_k_g0": sum(1 for r in adm if r["radius"] == r["k_g0"]),
        },
        "outside_illustration": {
            "rows": len(outside),
            "max_radius": max(r["radius"] for r in outside),
            "free_rule_exact": sum(1 for r in outside if free_rule(r) == r["radius"]),
            "free_rule_over_certifies": sum(1 for r in outside if free_rule(r) > r["radius"]),
        },
        "method_by_radius": {
            f"{m}:{rad}": c
            for (m, rad), c in sorted(
                collections.Counter((r["method"], r["radius"]) for r in adm).items()
            )
        },
        "by_coverage": {
            str(c): {
                "rows": len(s),
                "radius_1": sum(1 for r in s if r["radius"] == 1),
                "free_rule_exact": sum(1 for r in s if free_rule(r) == r["radius"]),
                "free_rule_over_certifies": sum(1 for r in s if free_rule(r) > r["radius"]),
                "undefined": sum(1 for r in s if r["separation_status"] != "measured"),
                "median_k_g0": sorted(r["k_g0"] for r in s)[len(s) // 2],
            }
            for c, s in ((c, [r for r in adm if r["coverage"] == c]) for c in (1.0, 0.5, 0.25))
        },
        "icc_radius_is_one": {
            "icc": round(icc, 4),
            "networks": len({r["network"] for r in adm}),
            "m0": round(m0, 2),
        },
    }
    (OUT / "ladder_summary.json").write_text(json.dumps(summary, indent=1))
    for key in (
        "rows",
        "distinct_pairs",
        "free_rule",
        "law_on_measurable_rows",
        "undefined_stratum",
        "outside_illustration",
        "icc_radius_is_one",
    ):
        print(key, json.dumps(summary[key]))


if __name__ == "__main__":
    main()
