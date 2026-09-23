"""Stage 0: the radius in asserted claims, on the committed rows.

Three checks, all read-only against ``results/axisa3/``.

1. For every admissible row at coverage 1.0, retract each asserted claim in
   turn, re-close under Meek, and test whether the committed optimal set is
   still GAC-valid. A row where some single retraction breaks validity has a
   radius of 1 in claim units, whatever its hop-count radius is. The fraction
   of single retractions that break validity (``phi_1``) is recorded as well.
2. Whether the knowledge sets produced by ``select_knowledge`` at coverage
   0.25, 0.5 and 1.0 are nested.
3. Whether the ``o_g0_not_identified`` rejections are amenable at ``G0``. A
   non-amenable query has no valid adjustment set, so nothing can be certified
   on it whatever the optimal-set definition.

``pathfinder`` is skipped by default: its 85-node chain component makes the
greedy recovering set slow to build, and it yields three rows, all censored.

Writes ``results/stage0/claim_radius_rows.csv`` and
``results/stage0/claim_radius_summary.json``.

    python experiments/stage0_claim_radius.py [--skip pathfinder]
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from bkrobust.benchmarks.describe import parse_file  # noqa: E402
from bkrobust.demo.evaluate import optimal_adjustment_set_dag  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag  # noqa: E402
from bkrobust.mpdag_criterion.criterion import is_amenable  # noqa: E402

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
    """Run the three checks and write the rows and the summary."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip", default="pathfinder", help="comma-separated networks to skip")
    args = ap.parse_args(argv)
    skip = {s for s in args.skip.split(",") if s}

    rows = [json.loads(line) for line in INSTANCES.open()]
    adm = [r for r in rows if r.get("admissible") and r["network"] not in skip]
    rejected = [
        r
        for r in rows
        if r.get("reject_reason") == "o_g0_not_identified" and r["network"] not in skip
    ]
    networks = sorted({r["network"] for r in adm} | {r["network"] for r in rejected})

    graphs: dict[str, tuple[MPDAG, MPDAG]] = {}
    recovering: dict[str, list] = {}
    for n in networks:
        t = time.perf_counter()
        graphs[n] = load(n)
        recovering[n] = sorted(knowledge_to_recover(*graphs[n]))
        print(
            f"loaded {n}: {len(recovering[n])} asserted claims, {time.perf_counter() - t:.1f} s",
            flush=True,
        )

    # 2. nesting of the coverage sweep
    nesting = {"networks": len(networks), "k05_in_k1": 0, "k025_in_k05": 0, "violations": []}
    for n in networks:
        k = {c: set(select_knowledge(recovering[n], c)) for c in (1.0, 0.5, 0.25)}
        nesting["k05_in_k1"] += k[0.5] <= k[1.0]
        nesting["k025_in_k05"] += k[0.25] <= k[0.5]
        if not k[0.25] <= k[0.5]:
            nesting["violations"].append(
                {"network": n, "sizes": [len(k[1.0]), len(k[0.5]), len(k[0.25])]}
            )

    # 3. amenability of the o_g0_not_identified rejections
    amen = collections.Counter()
    for r in rejected:
        dag, cpdag = graphs[r["network"]]
        g0 = apply_orientations(cpdag, select_knowledge(recovering[r["network"]], r["coverage"]))
        amen[
            "amenable" if g0 is not None and is_amenable(g0, r["x"], r["y"]) else "not_amenable"
        ] += 1

    # 1. single-claim retraction at coverage 1.0
    OUT.mkdir(parents=True, exist_ok=True)
    cov1 = [r for r in adm if r["coverage"] == 1.0]
    results = []
    t0 = time.perf_counter()
    with (OUT / "claim_radius_rows.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "network",
                "x",
                "y",
                "hop_radius",
                "n_asserted",
                "k_g0",
                "single_retractions_breaking_validity",
                "single_retractions_inconsistent",
                "claim_radius_is_one",
                "phi_1",
            ]
        )
        for i, r in enumerate(cov1):
            dag, cpdag = graphs[r["network"]]
            k = recovering[r["network"]]
            z = frozenset(optimal_adjustment_set_dag(dag, r["x"], r["y"]))
            broken = inconsistent = 0
            for claim in k:
                g = apply_orientations(cpdag, [e for e in k if e != claim])
                if g is None:
                    inconsistent += 1
                elif not is_gac_valid_mpdag(g, r["x"], r["y"], z):
                    broken += 1
            phi = broken / len(k) if k else None
            results.append((r["radius"], len(k), r["k_g0"], broken, phi))
            w.writerow(
                [
                    r["network"],
                    r["x"],
                    r["y"],
                    r["radius"],
                    len(k),
                    r["k_g0"],
                    broken,
                    inconsistent,
                    int(broken > 0),
                    phi,
                ]
            )
            if i % 100 == 0:
                print(f"  {i}/{len(cov1)} rows, {time.perf_counter() - t0:.0f} s", flush=True)

    phis = sorted(p for *_, p in results if p is not None)
    summary = {
        "skipped_networks": sorted(skip),
        "asserted_claims_per_network": {n: len(recovering[n]) for n in networks},
        "coverage_1_rows": len(results),
        "claim_radius_is_one": sum(1 for *_, b, _ in results if b > 0),
        "no_single_retraction_breaks_validity": sum(1 for *_, b, _ in results if b == 0),
        "hop_radius_where_claim_radius_is_one": dict(
            sorted(collections.Counter(h for h, *_, b, _ in results if b > 0).items())
        ),
        "phi_1": {"min": phis[0], "median": phis[len(phis) // 2], "max": phis[-1]},
        "median_asserted": sorted(n for _, n, *_ in results)[len(results) // 2],
        "median_k_g0": sorted(kg for _, _, kg, *_ in results)[len(results) // 2],
        "coverage_sweep_nesting": nesting,
        "o_g0_not_identified": {"rows": len(rejected), **amen},
    }
    (OUT / "claim_radius_summary.json").write_text(json.dumps(summary, indent=1))
    for key in (
        "coverage_1_rows",
        "claim_radius_is_one",
        "no_single_retraction_breaks_validity",
        "hop_radius_where_claim_radius_is_one",
        "phi_1",
        "median_asserted",
        "median_k_g0",
        "coverage_sweep_nesting",
        "o_g0_not_identified",
    ):
        print(key, json.dumps(summary[key]))


if __name__ == "__main__":
    main()
