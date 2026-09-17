"""Stage 1: does the closed-form set agree with the enumerating one where both exist?

The enumerating definition computes the optimal set in every DAG extension of
the analyst's graph and reports non-identification when they disagree. It is
exponential in the graph's undirected edges, so it caps what the pipeline can
evaluate at a small number of them. The closed form is polynomial and always
returns something on an amenable graph. This measures, on rows where the
enumerating definition terminates, how often the two agree, and it checks the
closed form's validity on every row including the ones enumeration cannot reach.

Three populations, because the answer may differ across them: the committed
admissible rows with the true recovering set (fully oriented analyst graph at
full coverage, partly oriented below it); the same rows under a truthful
subset of half the claims, where the analyst graph keeps free edges; and the
same rows under uniformly corrupted claims, where the graph is both partly open
and wrong.

Writes ``results/stage0/lever0_differential.json``.

    python experiments/lever0_differential.py --max-undirected 8
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

from bkrobust.benchmarks.describe import parse_file, to_mpdag  # noqa: E402
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag  # noqa: E402
from bkrobust.mpdag_criterion.optimal import committed_adjustment_set  # noqa: E402
from bkrobust.synth.knowledge import flip  # noqa: E402

INSTANCES = ROOT / "results" / "axisa3" / "instances.jsonl"
MODELS = ROOT / "results" / "axisa3" / "networks" / "example_models"
OUT = ROOT / "results" / "stage0"


def load(name: str) -> tuple[MPDAG, MPDAG]:
    """Parse one network and return its DAG and CPDAG."""
    path = next(p for p in sorted(MODELS.iterdir()) if p.name.startswith(name))
    parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
    dag = to_mpdag(parsed)
    return dag, dag_to_cpdag(dag)


def compare(g0: MPDAG, x: str, y: str, max_und: int) -> dict:
    """One row: the closed form against enumeration, and the closed form's validity."""
    cs = committed_adjustment_set(g0, x, y)
    rec = {
        "closed_verdict": cs.verdict,
        "closed_valid": bool(cs.z is not None and is_gac_valid_mpdag(g0, x, y, cs.z)),
        "g0_undirected": len(g0.undirected_edges),
        "enum": "not_run",
    }
    if len(g0.undirected_edges) <= max_und:
        o = optimal_adjustment_set_mpdag(g0, x, y)
        if o is None:
            rec["enum"] = "not_identified"
        elif not o:
            rec["enum"] = "empty"
        else:
            rec["enum"] = "set"
            rec["agree"] = cs.z is not None and set(cs.z) == set(o)
            rec["closed_subset_of_enum"] = cs.z is not None and set(cs.z) <= set(o)
            rec["enum_subset_of_closed"] = cs.z is not None and set(o) <= set(cs.z)
    return rec


def _stable(key: str) -> int:
    """A seed component that does not depend on the interpreter's hash salt."""
    return int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)


def main(argv: list[str] | None = None) -> None:
    """Run the three populations and write the differential."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-undirected", type=int, default=8)
    ap.add_argument("--seed", type=int, default=20260913)
    ap.add_argument("--skip", default="pathfinder")
    args = ap.parse_args(argv)
    skip = {s for s in args.skip.split(",") if s}

    rows = [json.loads(line) for line in INSTANCES.open()]
    pairs = sorted(
        {
            (r["network"], r["x"], r["y"])
            for r in rows
            if r.get("admissible") and r["network"] not in skip
        }
    )
    graphs: dict[str, tuple[MPDAG, MPDAG]] = {}
    recovering: dict[str, list] = {}
    for n in sorted({p[0] for p in pairs}):
        graphs[n] = load(n)
        recovering[n] = sorted(knowledge_to_recover(*graphs[n]))

    pops = {"true_full": [], "true_half": [], "corrupted_half": []}
    t0 = time.perf_counter()
    for i, (net, x, y) in enumerate(pairs):
        _, cpdag = graphs[net]
        k = recovering[net]
        rng = np.random.default_rng([args.seed, _stable(f"{net}|{x}|{y}")])
        half = [k[int(j)] for j in sorted(rng.permutation(len(k))[: max(1, len(k) // 2)])]
        for name, kk in (
            ("true_full", k),
            ("true_half", half),
            ("corrupted_half", flip(half, rng, 0.5)),
        ):
            g0 = apply_orientations(cpdag, kk)
            if g0 is None:
                pops[name].append({"closed_verdict": "knowledge_inconsistent", "enum": "not_run"})
                continue
            pops[name].append(compare(g0, x, y, args.max_undirected))
        if i % 100 == 0:
            print(f"  {i}/{len(pairs)} pairs, {time.perf_counter() - t0:.0f} s", flush=True)

    def digest(recs: list[dict]) -> dict:
        both = [
            r
            for r in recs
            if r.get("enum") == "set" and r["closed_verdict"] in ("closed_form", "canonical")
        ]
        return {
            "rows": len(recs),
            "closed_verdict": dict(collections.Counter(r["closed_verdict"] for r in recs)),
            "closed_valid_where_set": sum(1 for r in recs if r.get("closed_valid")),
            "closed_invalid_where_set": sum(
                1
                for r in recs
                if r["closed_verdict"] in ("closed_form", "canonical") and not r["closed_valid"]
            ),
            "enumeration": dict(collections.Counter(r["enum"] for r in recs)),
            "both_defined": len(both),
            "agree": sum(1 for r in both if r.get("agree")),
            "closed_strict_subset_of_enum": sum(
                1 for r in both if r.get("closed_subset_of_enum") and not r.get("agree")
            ),
            "enum_strict_subset_of_closed": sum(
                1 for r in both if r.get("enum_subset_of_closed") and not r.get("agree")
            ),
            "closed_returns_set_where_enum_not_identified": sum(
                1
                for r in recs
                if r.get("enum") == "not_identified"
                and r["closed_verdict"] in ("closed_form", "canonical")
            ),
            "closed_returns_set_where_enum_not_run": sum(
                1
                for r in recs
                if r.get("enum") == "not_run"
                and r["closed_verdict"] in ("closed_form", "canonical")
            ),
            "undirected_edges_max": max((r.get("g0_undirected", 0) for r in recs), default=0),
        }

    summary = {
        "pairs": len(pairs),
        "max_undirected_for_enumeration": args.max_undirected,
        "seconds": round(time.perf_counter() - t0, 1),
        "populations": {k: digest(v) for k, v in pops.items()},
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "lever0_differential.json").write_text(json.dumps(summary, indent=1))
    for k, v in summary["populations"].items():
        print(f"\n{k}:")
        for kk, vv in v.items():
            print(f"  {kk}: {vv}")


if __name__ == "__main__":
    main()
