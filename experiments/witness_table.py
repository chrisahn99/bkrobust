"""Phase B: the first-failing retraction set, named.

A radius is a number. What an analyst can act on is the sentence the number
comes from: the estimate holds unless these specific claims are taken back. This
finds that set for every certifiable row, in the analyst's own variable names,
and prints beside it what each cheap statistic said about the same row.

A count baseline can return a number. It cannot return which claim. That is the
column this table exists for, and it is the last thing a referee reads.

Note what the last three columns are and are not. They print what each cheap
statistic said **in its own unit**, beside the retraction depth **in claims**.
The hop radius and the closure size read high against a claim depth by
construction, since one asserted claim can carry many covering steps, so the gap
between them is the unit mismatch and not an error either instrument made. It is
reported here because that mismatch is the reason the certificate has to be
re-denominated before any of these numbers is put in front of an analyst.

The witness is the lexicographically first subset at the smallest depth that
breaks validity, which makes it reproducible rather than arbitrary; where
several subsets at that depth break it, the count is reported too, since an
analyst facing six single claims that each break the estimate is in a different
position from one facing one.

Writes ``results/frame/witness_rows.csv`` and ``results/frame/witness_summary.json``.

    python experiments/witness_table.py --coverage 1.0
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import itertools
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from bkrobust.benchmarks.describe import parse_file, to_mpdag  # noqa: E402
from bkrobust.benchmarks.measure import separation  # noqa: E402
from bkrobust.demo.evaluate import (  # noqa: E402
    is_valid_adjustment_set_dag,
    optimal_adjustment_set_mpdag,
)
from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag  # noqa: E402
from bkrobust.hybrid import breakdown_radius  # noqa: E402
from bkrobust.mpdag_criterion.criterion import is_amenable  # noqa: E402

MODELS = ROOT / "results" / "axisa3" / "networks" / "example_models"
INSTANCES = ROOT / "results" / "axisa3" / "instances.jsonl"
OUT = ROOT / "results" / "frame"

MAX_DEPTH = 2
SUBSET_BUDGET = 400


def load(name: str) -> tuple[MPDAG, MPDAG, dict[str, str]]:
    """Parse one network and return its DAG, its CPDAG and the raw-name map."""
    path = next(p for p in sorted(MODELS.iterdir()) if p.name.startswith(name))
    parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
    dag = to_mpdag(parsed)
    return dag, dag_to_cpdag(dag), dict(parsed.name_map)


def witness(
    cpdag: MPDAG, k: list, x: str, y: str, z: frozenset
) -> tuple[int | None, tuple, int, str]:
    """The smallest retraction that invalidates ``z``, and how many there are at that depth.

    Returns ``(depth, the first such subset, how many subsets break it at that
    depth, what the retraction costs)``. The cost is ``validity`` when the set
    merely stops being valid and ``identifiability`` when the query also stops
    being amenable, which is the worse of the two for an analyst: there is then
    no adjustment set to fall back on.
    """
    used = 0
    for depth in range(1, MAX_DEPTH + 1):
        if depth > len(k):
            break
        hits: list[tuple] = []
        cost = "validity"
        for subset in itertools.combinations(range(len(k)), depth):
            if used >= SUBSET_BUDGET:
                break
            drop = set(subset)
            kept = [e for i, e in enumerate(k) if i not in drop]
            g = apply_orientations(cpdag, kept)
            used += 1
            if g is None:
                continue
            if not is_gac_valid_mpdag(g, x, y, z):
                hits.append(tuple(k[i] for i in subset))
                if not hits[:-1] and not is_amenable(g, x, y):
                    cost = "identifiability"
        if hits:
            return depth, hits[0], len(hits), cost
    return None, (), 0, ""


def where(edges: tuple, x: str, y: str) -> str:
    """Where the retraction sits relative to the query.

    The locality the project measured for orientation statements, restated in
    the unit an analyst asserts: a claim pointing into the treatment is the one
    that opens a back-door path when it is taken back.
    """
    if len(edges) != 1:
        return "multiple claims"
    a, b = edges[0]
    if b == x:
        return "into the treatment"
    if a == x:
        return "out of the treatment"
    if b == y:
        return "into the outcome"
    if a == y:
        return "out of the outcome"
    return "elsewhere in the component"


def render(edges: tuple, names: dict[str, str]) -> str:
    """The retraction set as an analyst would read it."""
    return " and ".join(f"{names.get(a, a)} -> {names.get(b, b)}" for a, b in edges)


def main(argv: list[str] | None = None) -> None:
    """Build the decision table over the committed admissible pairs."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--coverage", type=float, default=1.0)
    ap.add_argument("--skip", default="pathfinder")
    args = ap.parse_args(argv)
    skip = {s for s in args.skip.split(",") if s}

    rows = [json.loads(line) for line in INSTANCES.open()]
    adm = [
        r
        for r in rows
        if r.get("admissible") and r["coverage"] == args.coverage and r["network"] not in skip
    ]

    graphs: dict[str, tuple[MPDAG, MPDAG, dict[str, str]]] = {}
    recovering: dict[str, list] = {}
    for name in sorted({r["network"] for r in adm}):
        graphs[name] = load(name)
        recovering[name] = sorted(knowledge_to_recover(graphs[name][0], graphs[name][1]))

    OUT.mkdir(parents=True, exist_ok=True)
    fields = [
        "network",
        "treatment",
        "outcome",
        "n_claims",
        "retraction_depth",
        "retraction_set",
        "retraction_sits",
        "n_sets_at_depth",
        "what_it_costs",
        "free_rule_said",
        "hop_radius_said",
        "closure_size_said",
        "claim_count_said",
        "free_rule_unit_gap",
        "hop_unit_gap",
        "closure_unit_gap",
        "z_valid_at_truth",
    ]
    out_rows: list[dict] = []
    t0 = time.perf_counter()
    for r in adm:
        dag, cpdag, names = graphs[r["network"]]
        k = recovering[r["network"]]
        x, y = r["x"], r["y"]
        g0 = apply_orientations(cpdag, k)
        o = optimal_adjustment_set_mpdag(g0, x, y)
        if not o:
            continue
        z = frozenset(o)
        depth, subset, n_sets, cost = witness(cpdag, k, x, y, z)
        if depth is None:
            continue
        sep, sep_status = separation(cpdag, x, z)
        free = sep if sep_status == "measured" else 1
        hop = breakdown_radius(cpdag, None, x, y, z, g0=g0, time_limit_s=60.0)
        k_g0 = sum(
            1 for (a, b) in g0.directed_edges if tuple(sorted((a, b))) in cpdag.undirected_edges
        )
        out_rows.append(
            {
                "network": r["network"],
                "treatment": names.get(x, x),
                "outcome": names.get(y, y),
                "n_claims": len(k),
                "retraction_depth": depth,
                "retraction_set": render(subset, names),
                "retraction_sits": where(subset, x, y),
                "n_sets_at_depth": n_sets,
                "what_it_costs": cost,
                "free_rule_said": free,
                "hop_radius_said": hop.radius,
                "closure_size_said": k_g0,
                "claim_count_said": len(k),
                "free_rule_unit_gap": free - depth,
                "hop_unit_gap": hop.radius - depth,
                "closure_unit_gap": k_g0 - depth,
                "z_valid_at_truth": int(is_valid_adjustment_set_dag(dag, x, y, z)),
            }
        )

    with (OUT / "witness_rows.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(out_rows)

    def reads_high(field: str) -> dict:
        """How far an instrument's own unit sits above the claim depth."""
        gaps = [r[field] for r in out_rows]
        return {
            "reads_high": sum(1 for g in gaps if g > 0),
            "exact": sum(1 for g in gaps if g == 0),
            "reads_low": sum(1 for g in gaps if g < 0),
            "max_gap": max(gaps),
        }

    over = [r for r in out_rows if r["free_rule_unit_gap"] > 0]
    summary = {
        "coverage": args.coverage,
        "rows": len(out_rows),
        "networks": len({r["network"] for r in out_rows}),
        "seconds": round(time.perf_counter() - t0, 1),
        "retraction_depth": dict(
            sorted(collections.Counter(r["retraction_depth"] for r in out_rows).items())
        ),
        "what_it_costs": dict(collections.Counter(r["what_it_costs"] for r in out_rows)),
        "sets_at_depth": {
            "one_set_only": sum(1 for r in out_rows if r["n_sets_at_depth"] == 1),
            "median": sorted(r["n_sets_at_depth"] for r in out_rows)[len(out_rows) // 2],
            "max": max(r["n_sets_at_depth"] for r in out_rows),
        },
        "unit_gap_against_claim_depth": {
            "free_rule": reads_high("free_rule_unit_gap"),
            "hop_radius": reads_high("hop_unit_gap"),
            "closure_size": reads_high("closure_unit_gap"),
            "note": (
                "A positive gap is the unit mismatch, not an error: a hop count and a closure "
                "size are not counts of retractable claims. Every instrument denominated in "
                "anything but asserted claims reads high against the depth an analyst can act on."
            ),
        },
        "free_rule_reads_high": {
            "rows": len(over),
            "by_network": dict(collections.Counter(r["network"] for r in over)),
        },
        "where_the_retraction_sits": {
            "all_rows": dict(
                collections.Counter(r["retraction_sits"] for r in out_rows).most_common()
            ),
            "where_free_rule_reads_high": dict(
                collections.Counter(
                    r["retraction_sits"] for r in out_rows if r["free_rule_unit_gap"] > 0
                ).most_common()
            ),
            "where_free_rule_is_exact": dict(
                collections.Counter(
                    r["retraction_sits"] for r in out_rows if r["free_rule_unit_gap"] == 0
                ).most_common()
            ),
        },
        "single_claim_networks": {
            r["network"]: r["n_claims"] for r in out_rows if r["n_claims"] == 1
        },
        "widest_gap_between_hop_and_claim": sorted(
            (
                {
                    "network": r["network"],
                    "treatment": r["treatment"],
                    "outcome": r["outcome"],
                    "hop": r["hop_radius_said"],
                    "claim_depth": r["retraction_depth"],
                    "retraction_set": r["retraction_set"],
                }
                for r in out_rows
            ),
            key=lambda d: -(d["hop"] - d["claim_depth"]),
        )[:5],
    }
    (OUT / "witness_summary.json").write_text(json.dumps(summary, indent=1))
    head = {k: v for k, v in summary.items() if k != "free_rule_reads_high"}
    print(json.dumps(head, indent=1))
    print("\nthe sentence the table exists for, one row per network:")
    seen: set[str] = set()
    for r in out_rows:
        if r["network"] in seen:
            continue
        seen.add(r["network"])
        if len(seen) > 8:
            break
        print(
            f"  {r['network']}: the estimate for {r['treatment']} on {r['outcome']} holds "
            f"unless {r['retraction_set']} is retracted"
        )


if __name__ == "__main__":
    main()
