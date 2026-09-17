"""Stage 0: do the committed radii reproduce from the committed inputs?

Nothing downstream is worth running until they do. For a seeded sample of the
admissible rows in ``results/axisa3/instances.jsonl`` this rebuilds the inputs
the way ``benchmarks.measure.evaluate`` does, calls ``hybrid.breakdown_radius``,
and compares the radius, the dispatch leg and the exactness flag against the
committed row.

Writes ``results/stage0/reproduce.json``.

    python experiments/stage0_reproduce.py --sample 30 --seed 20260912
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
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.hybrid import breakdown_radius  # noqa: E402

INSTANCES = ROOT / "results" / "axisa3" / "instances.jsonl"
MODELS = ROOT / "results" / "axisa3" / "networks" / "example_models"
OUT = ROOT / "results" / "stage0"


def select_knowledge(k_true: list, coverage: float) -> list:
    """The stride of ``benchmarks.measure.select_knowledge`` on a precomputed recovering set."""
    if coverage >= 1.0:
        return list(k_true)
    keep = round(len(k_true) * coverage)
    if keep <= 0:
        return []
    step = len(k_true) / keep
    return [k_true[min(len(k_true) - 1, int(i * step))] for i in range(keep)]


def load(name: str) -> tuple[MPDAG, MPDAG]:
    """Parse one network file and return its DAG and CPDAG."""
    path = next(p for p in sorted(MODELS.iterdir()) if p.name.startswith(name))
    parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
    dag = MPDAG(parsed.nodes, directed=parsed.edges)
    return dag, dag_to_cpdag(dag)


def main(argv: list[str] | None = None) -> None:
    """Recompute a sample of committed radii and report every disagreement."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", type=int, default=30)
    ap.add_argument("--seed", type=int, default=20260912)
    ap.add_argument("--skip", default="pathfinder")
    args = ap.parse_args(argv)
    skip = {s for s in args.skip.split(",") if s}

    rows = [json.loads(line) for line in INSTANCES.open()]
    adm = [r for r in rows if r.get("admissible") and r["network"] not in skip]

    # Stratify by dispatch leg so the ladder branch is exercised, not only the fast one.
    by_leg = collections.defaultdict(list)
    for r in adm:
        by_leg[r["method"]].append(r)
    rng = np.random.default_rng(args.seed)
    picked = []
    for _, group in sorted(by_leg.items()):
        take = min(len(group), max(1, args.sample // len(by_leg)))
        idx = rng.choice(len(group), size=take, replace=False)
        picked.extend(group[int(i)] for i in sorted(idx))

    graphs: dict[str, tuple[MPDAG, MPDAG]] = {}
    recovering: dict[str, list] = {}
    checked, disagreements = [], []
    t0 = time.perf_counter()
    for r in picked:
        net = r["network"]
        if net not in graphs:
            graphs[net] = load(net)
            recovering[net] = sorted(knowledge_to_recover(*graphs[net]))
        _, cpdag = graphs[net]
        k = select_knowledge(recovering[net], r["coverage"])
        g0 = apply_orientations(cpdag, k)
        z = frozenset(optimal_adjustment_set_mpdag(g0, r["x"], r["y"]))
        out = breakdown_radius(cpdag, None, r["x"], r["y"], z, g0=g0, time_limit_s=300.0)
        row = {
            "network": net,
            "x": r["x"],
            "y": r["y"],
            "coverage": r["coverage"],
            "committed": {"radius": r["radius"], "method": r["method"], "exact": r["exact"]},
            "recomputed": {"radius": out.radius, "method": out.method, "exact": out.exact},
            "z_size_committed": r["z_size"],
            "z_size_recomputed": len(z),
        }
        checked.append(row)
        if (out.radius, out.method, out.exact, len(z)) != (
            r["radius"],
            r["method"],
            r["exact"],
            r["z_size"],
        ):
            disagreements.append(row)

    OUT.mkdir(parents=True, exist_ok=True)
    report = {
        "seed": args.seed,
        "skipped_networks": sorted(skip),
        "rows_checked": len(checked),
        "by_dispatch_leg": dict(collections.Counter(r["committed"]["method"] for r in checked)),
        "disagreements": len(disagreements),
        "disagreement_rows": disagreements,
        "seconds": round(time.perf_counter() - t0, 1),
        "rows": checked,
    }
    (OUT / "reproduce.json").write_text(json.dumps(report, indent=1))
    print(
        f"checked {len(checked)} rows across {len(report['by_dispatch_leg'])} dispatch legs "
        f"in {report['seconds']} s: {len(disagreements)} disagreements"
    )
    for row in disagreements:
        print("  ", json.dumps(row))


if __name__ == "__main__":
    main()
