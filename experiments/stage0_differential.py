"""Stage 0: do the patched measurement path and the committed rows still agree?

Five changes landed in ``benchmarks/measure.py``: a nested selection mode, a
three-valued verdict for the optimal set, ``g0_undirected_edges`` assigned on
every branch, the dispatch leg and the search budget as columns, and the
rejection label ``not_amenable`` split out of ``o_g0_not_identified``. Four are
additive. Two change something a reader could quote, and this is the test that
says exactly what changed.

Under the default selection mode every committed field must be reproduced,
except the rejection label, whose new values decompose the old one. Under the
nested mode the coverage-1.0 rows must still be reproduced, because the two
modes coincide there, and the partial-coverage delta is the population table.

Writes ``results/stage0/differential.json``.

    python experiments/stage0_differential.py --per-network 5
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from bkrobust.benchmarks.describe import parse_file  # noqa: E402
from bkrobust.benchmarks.measure import evaluate  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402

INSTANCES = ROOT / "results" / "axisa3" / "instances.jsonl"
MODELS = ROOT / "results" / "axisa3" / "networks" / "example_models"
OUT = ROOT / "results" / "stage0"

#: The old label and the labels that now decompose it.
LABEL_SPLIT = {"o_g0_not_identified": {"o_g0_not_identified", "not_amenable", "o_g0_empty"}}
#: Fields that must be identical under the default mode.
COMPARED = (
    "admissible",
    "radius",
    "method",
    "exact",
    "separation",
    "separation_status",
    "k_g0",
    "z_size",
    "component_size",
)


def load(name: str) -> tuple[MPDAG, MPDAG]:
    """Parse one network file and return its DAG and CPDAG."""
    path = next(p for p in sorted(MODELS.iterdir()) if p.name.startswith(name))
    parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
    dag = MPDAG(parsed.nodes, directed=parsed.edges)
    return dag, dag_to_cpdag(dag)


def compare(committed: dict, got: dict) -> list[str]:
    """Field names on which the recomputed row differs from the committed one."""
    bad = []
    for f in COMPARED:
        if committed.get(f) != got.get(f):
            bad.append(f)
    old, new = committed.get("reject_reason", ""), got.get("reject_reason", "")
    if old != new and new not in LABEL_SPLIT.get(old, {old}):
        bad.append("reject_reason")
    return bad


def main(argv: list[str] | None = None) -> None:
    """Recompute a sample of committed rows under both selection modes."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-network", type=int, default=5)
    ap.add_argument("--seed", type=int, default=20260912)
    ap.add_argument("--skip", default="pathfinder")
    args = ap.parse_args(argv)
    skip = {s for s in args.skip.split(",") if s}

    rows = [json.loads(line) for line in INSTANCES.open()]
    committed = {
        (r["network"], r["x"], r["y"], r["coverage"]): r
        for r in rows
        if r.get("x") and r["network"] not in skip
    }
    adm = [r for r in rows if r.get("admissible") and r["network"] not in skip]

    by_net = collections.defaultdict(list)
    for r in adm:
        by_net[r["network"]].append((r["x"], r["y"]))
    rng = np.random.default_rng(args.seed)
    picked: list[tuple[str, str, str]] = []
    for net in sorted(by_net):
        pairs = sorted(set(by_net[net]))
        take = min(len(pairs), args.per_network)
        idx = rng.choice(len(pairs), size=take, replace=False)
        picked.extend((net, *pairs[int(i)]) for i in sorted(idx))

    graphs: dict[str, tuple[MPDAG, MPDAG]] = {}
    report: dict = {
        "seed": args.seed,
        "skipped_networks": sorted(skip),
        "pairs": len(picked),
        "stride": {"checked": 0, "mismatched": 0, "mismatch_rows": [], "label_decomposition": {}},
        "nested": {
            "coverage_1_checked": 0,
            "coverage_1_mismatched": 0,
            "partial_coverage_changed": 0,
            "partial_coverage_checked": 0,
            "radius_delta": {},
            "status_delta": {},
        },
    }
    t0 = time.perf_counter()
    for net, x, y in picked:
        if net not in graphs:
            graphs[net] = load(net)
        dag, cpdag = graphs[net]
        for coverage in (1.0, 0.5, 0.25):
            key = (net, x, y, coverage)
            if key not in committed:
                continue
            old = committed[key]

            got = evaluate(net, dag, cpdag, x, y, coverage).as_row()
            report["stride"]["checked"] += 1
            bad = compare(old, got)
            if bad:
                report["stride"]["mismatched"] += 1
                report["stride"]["mismatch_rows"].append(
                    {
                        "key": list(key),
                        "fields": bad,
                        "committed": {f: old.get(f) for f in bad},
                        "recomputed": {f: got.get(f) for f in bad},
                    }
                )
            if old.get("reject_reason") in LABEL_SPLIT:
                d = report["stride"]["label_decomposition"]
                d[got["reject_reason"]] = d.get(got["reject_reason"], 0) + 1

            nested = evaluate(net, dag, cpdag, x, y, coverage, selection_mode="nested").as_row()
            if coverage >= 1.0:
                report["nested"]["coverage_1_checked"] += 1
                if compare(old, nested):
                    report["nested"]["coverage_1_mismatched"] += 1
            else:
                report["nested"]["partial_coverage_checked"] += 1
                if (nested.get("admissible"), nested.get("radius")) != (
                    got.get("admissible"),
                    got.get("radius"),
                ):
                    report["nested"]["partial_coverage_changed"] += 1
                    delta = f"{got.get('radius')}->{nested.get('radius')}"
                    report["nested"]["radius_delta"][delta] = (
                        report["nested"]["radius_delta"].get(delta, 0) + 1
                    )
                    sd = f"{got.get('reject_reason')}->{nested.get('reject_reason')}"
                    report["nested"]["status_delta"][sd] = (
                        report["nested"]["status_delta"].get(sd, 0) + 1
                    )

    report["seconds"] = round(time.perf_counter() - t0, 1)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "differential.json").write_text(json.dumps(report, indent=1))
    print(json.dumps({k: v for k, v in report.items() if k != "stride"}, indent=1))
    print(
        "stride:", json.dumps({k: v for k, v in report["stride"].items() if k != "mismatch_rows"})
    )
    for row in report["stride"]["mismatch_rows"][:10]:
        print("  MISMATCH", json.dumps(row))


if __name__ == "__main__":
    main()
