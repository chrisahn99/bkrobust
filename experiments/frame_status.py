"""Phase B: what the frozen frame yields under each knowledge source.

Three sources answer on the same frozen rows, so the arms differ in status
composition rather than in population.

``A-TRUE`` is the audit control: a truthful subset of the recovering set, taken
as a prefix of one seeded permutation so that a smaller coverage is a subset of
a larger one. ``D-RAND`` is the vacuity control and the reason any knowledge
result is a measurement of knowledge rather than of cardinality: the same edges
as ``A-TRUE``, so size, granularity and distance from the treatment match
exactly, with each orientation drawn uniformly among the choices Meek closure
accepts. ``D-DEGEN`` asserts nothing.

Every quantity on the computation path is a function of the CPDAG, the
knowledge and the query. The truth is read once at the end of each row, to score
the commission count and whether the committed set is valid at the truth, and
never before.

Writes ``results/frame/status_rows.csv`` and ``results/frame/status_summary.json``.

    python experiments/frame_status.py --max-nodes 500 --reps 3
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import itertools
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

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
FRAME = ROOT / "results" / "frame" / "frame.jsonl"
OUT = ROOT / "results" / "frame"

MAX_DEPTH = 3
SUBSET_BUDGET = 200
#: Undirected edges in the analyst graph above which the extension enumeration
#: that defines the committed optimal set is not entered.
MAX_UNDIRECTED = 8
ARMS = ("A_TRUE", "D_RAND", "D_DEGEN")


def load(name: str) -> tuple[MPDAG, MPDAG]:
    """Parse one network and return its DAG and CPDAG."""
    path = next(p for p in sorted(MODELS.iterdir()) if p.name.startswith(name))
    parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
    dag = to_mpdag(parsed)
    return dag, dag_to_cpdag(dag)


def nested_prefix(k_true: list, coverage: float, seed: int) -> list:
    """A prefix of one seeded permutation, so coverage is a retraction sequence."""
    if coverage >= 1.0:
        return list(k_true)
    keep = round(len(k_true) * coverage)
    if keep <= 0:
        return []
    order = np.random.default_rng(seed).permutation(len(k_true))
    return [k_true[int(i)] for i in sorted(order[:keep])]


def chance_level(
    cpdag: MPDAG, base: list, rng: np.random.Generator, retries: int = 12
) -> list | None:
    """The same edges as ``base``, each oriented uniformly among Meek-consistent choices."""
    for _ in range(retries):
        drawn = [(a, b) if rng.integers(2) else (b, a) for (a, b) in base]
        if apply_orientations(cpdag, drawn) is not None:
            return drawn
    return None


def breaking_depths(
    cpdag: MPDAG, k: list, x: str, y: str, z: frozenset
) -> tuple[int | None, str, int | None]:
    """Smallest retraction breaking validity, and the one breaking identifiability."""
    r_claim = r_id = None
    status = f"gt_{min(len(k), MAX_DEPTH)}"
    used = 0
    for depth in range(1, MAX_DEPTH + 1):
        if depth > len(k):
            break
        width = math.comb(len(k), depth)
        if used + width > SUBSET_BUDGET:
            if r_claim is None:
                status = f"censored_at_depth_{depth}"
            break
        for subset in itertools.combinations(range(len(k)), depth):
            drop = set(subset)
            g = apply_orientations(cpdag, [e for i, e in enumerate(k) if i not in drop])
            used += 1
            if g is None:
                continue
            if r_claim is None and not is_gac_valid_mpdag(g, x, y, z):
                r_claim, status = depth, "exact"
            if r_id is None and not is_amenable(g, x, y):
                r_id = depth
            if r_claim is not None and r_id is not None:
                break
        if r_claim is not None and r_id is not None:
            break
    return r_claim, status, r_id


def score(
    dag: MPDAG, cpdag: MPDAG, x: str, y: str, k: list, radius_limit: float
) -> tuple[str, dict]:
    """Status of one (row, knowledge) pair and, where it is certifiable, its radii."""
    g0 = apply_orientations(cpdag, k)
    if g0 is None:
        return "knowledge_inconsistent", {}
    n_und = len(g0.undirected_edges)
    if not is_amenable(g0, x, y):
        return "not_amenable", {"g0_undirected": n_und}
    if n_und > MAX_UNDIRECTED:
        return "extensions_intractable", {"g0_undirected": n_und}
    o = optimal_adjustment_set_mpdag(g0, x, y)
    if o is None:
        return "optimal_not_identified", {"g0_undirected": n_und}
    if not o:
        return "optimal_empty", {"g0_undirected": n_und}
    z = frozenset(o)
    if not is_gac_valid_mpdag(g0, x, y, z):
        return "z_invalid_at_g0", {"g0_undirected": n_und}

    sep, sep_status = separation(cpdag, x, z)
    k_g0 = sum(1 for (a, b) in g0.directed_edges if tuple(sorted((a, b))) in cpdag.undirected_edges)
    # The bounded search answers in microseconds when a failure is near. Proving
    # that none exists is the expensive direction, and on this arm a third of the
    # rows are in it, so the ladder is given a short wall and its exactness flag
    # is carried rather than a censored answer being read as a measurement.
    hop = breakdown_radius(cpdag, None, x, y, z, g0=g0, time_limit_s=radius_limit)
    r_claim, claim_status, r_id = breaking_depths(cpdag, k, x, y, z)
    n_break = sum(
        1
        for i in range(len(k))
        if (gi := apply_orientations(cpdag, [e for j, e in enumerate(k) if j != i])) is not None
        and not is_gac_valid_mpdag(gi, x, y, z)
    )
    return "ok", {
        "g0_undirected": n_und,
        "z_size": len(z),
        "separation": sep,
        "separation_status": sep_status,
        "k_g0": k_g0,
        "r_hop": hop.radius,
        "hop_method": hop.method,
        "hop_exact": int(hop.exact),
        "r_claim": r_claim,
        "r_claim_status": claim_status,
        "r_id": r_id,
        "phi_1": round(n_break / len(k), 4) if k else None,
        # read once, after every observable quantity is fixed
        "n_wrong_claims": sum(1 for (a, b) in k if not dag.is_directed_edge(a, b)),
        "z_valid_at_truth": int(is_valid_adjustment_set_dag(dag, x, y, z)),
    }


def _stable(key: str) -> int:
    """A seed component that does not depend on the interpreter's hash salt."""
    return int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)


def main(argv: list[str] | None = None) -> None:
    """Score every frame row under each arm and publish the status composition."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-nodes", type=int, default=500)
    ap.add_argument("--coverages", default="1.0,0.5")
    ap.add_argument("--reps", type=int, default=3, help="chance-level draws per row")
    ap.add_argument("--seed", type=int, default=20260913)
    ap.add_argument(
        "--radius-limit",
        type=float,
        default=5.0,
        help="seconds the hop-radius ladder may spend before its answer is inexact",
    )
    args = ap.parse_args(argv)
    coverages = [float(v) for v in args.coverages.split(",")]

    frame = [json.loads(line) for line in FRAME.open()]
    by_net = collections.defaultdict(list)
    for row in frame:
        by_net[row["network"]].append(row)

    graphs: dict[str, tuple[MPDAG, MPDAG]] = {}
    recovering: dict[str, list] = {}
    excluded: list[dict] = []
    for name in sorted(by_net):
        dag, cpdag = load(name)
        if len(dag.nodes) > args.max_nodes:
            excluded.append({"network": name, "nodes": len(dag.nodes)})
            continue
        graphs[name] = (dag, cpdag)
        recovering[name] = sorted(knowledge_to_recover(dag, cpdag))
        print(f"loaded {name}: {len(dag.nodes)} nodes, {len(recovering[name])} claims", flush=True)
    print(f"excluded for size: {excluded}", flush=True)

    OUT.mkdir(parents=True, exist_ok=True)
    fields = [
        "network",
        "x",
        "y",
        "arm",
        "coverage",
        "rep",
        "n_k",
        "status",
        "g0_undirected",
        "z_size",
        "separation",
        "separation_status",
        "k_g0",
        "r_hop",
        "hop_method",
        "hop_exact",
        "r_claim",
        "r_claim_status",
        "r_id",
        "phi_1",
        "n_wrong_claims",
        "z_valid_at_truth",
    ]
    fh = (OUT / "status_rows.csv").open("w", newline="")
    writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
    writer.writeheader()

    counts: dict[tuple[str, float], collections.Counter] = collections.defaultdict(
        collections.Counter
    )
    radii: dict[tuple[str, float], list[dict]] = collections.defaultdict(list)
    t0 = time.perf_counter()
    done = 0
    for name in sorted(graphs):
        dag, cpdag = graphs[name]
        k_true = recovering[name]
        for row in by_net[name]:
            x, y = row["x"], row["y"]
            for coverage in coverages:
                base = nested_prefix(k_true, coverage, args.seed + _stable(name) % 9973)
                plans: list[tuple[str, int, list | None]] = [("A_TRUE", 0, base)]
                for rep in range(args.reps):
                    rng = np.random.default_rng(
                        [args.seed, _stable(f"{name}|{x}|{y}"), int(coverage * 100), rep]
                    )
                    plans.append(("D_RAND", rep, chance_level(cpdag, base, rng)))
                if coverage == coverages[0]:
                    plans.append(("D_DEGEN", 0, []))
                for arm, rep, k in plans:
                    done += 1
                    if k is None:
                        counts[(arm, coverage)]["chance_draw_failed"] += 1
                        continue
                    status, extra = score(dag, cpdag, x, y, k, args.radius_limit)
                    counts[(arm, coverage)][status] += 1
                    writer.writerow(
                        {
                            "network": name,
                            "x": x,
                            "y": y,
                            "arm": arm,
                            "coverage": coverage,
                            "rep": rep,
                            "n_k": len(k),
                            "status": status,
                            **extra,
                        }
                    )
                    if status == "ok":
                        radii[(arm, coverage)].append({"network": name, **extra})
            if done % 200 < 6:
                print(f"  {done} evaluations, {time.perf_counter() - t0:.0f} s", flush=True)
    fh.close()

    def digest(rows: list[dict]) -> dict:
        if not rows:
            return {}
        hops = collections.Counter(r["r_hop"] for r in rows)
        claims = collections.Counter(
            r["r_claim"] if r["r_claim"] is not None else r["r_claim_status"] for r in rows
        )
        wrong = [r["n_wrong_claims"] for r in rows]
        return {
            "rows": len(rows),
            "networks": len({r["network"] for r in rows}),
            "r_hop": dict(sorted(hops.items(), key=lambda kv: str(kv[0]))),
            "r_claim": dict(sorted(claims.items(), key=lambda kv: str(kv[0]))),
            "median_k_g0": sorted(r["k_g0"] for r in rows)[len(rows) // 2],
            "mean_wrong_claims": round(sum(wrong) / len(wrong), 3),
            "z_invalid_at_truth": sum(1 for r in rows if not r["z_valid_at_truth"]),
            "separation_undefined": sum(1 for r in rows if r["separation_status"] != "measured"),
            "hop_unreached": sum(1 for r in rows if r["r_hop"] == -1),
            "hop_inexact": sum(1 for r in rows if not r["hop_exact"]),
        }

    summary = {
        "frame_rows": len(frame),
        "frame_rows_evaluated": sum(len(by_net[n]) for n in graphs),
        "networks_evaluated": len(graphs),
        "excluded_for_size": excluded,
        "evaluations": done,
        "seconds": round(time.perf_counter() - t0, 1),
        "settings": {
            "coverages": coverages,
            "reps": args.reps,
            "seed": args.seed,
            "max_nodes": args.max_nodes,
            "max_undirected": MAX_UNDIRECTED,
        },
        "status_by_arm": {f"{a}@{c}": dict(v) for (a, c), v in sorted(counts.items())},
        "radii_by_arm": {f"{a}@{c}": digest(v) for (a, c), v in sorted(radii.items())},
    }
    (OUT / "status_summary.json").write_text(json.dumps(summary, indent=1))
    print("\nstatus composition on the frozen frame:")
    for key, v in summary["status_by_arm"].items():
        total = sum(v.values())
        ok = v.get("ok", 0)
        print(f"  {key:16s} n={total:5d}  ok={ok:5d} ({ok / total:6.1%})  {dict(v)}")
    print("\nradii where certifiable:")
    for key, v in summary["radii_by_arm"].items():
        if v:
            print(f"  {key:16s} {json.dumps(v)}")


if __name__ == "__main__":
    main()
