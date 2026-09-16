"""Post-hoc diagnostic: is the contradiction rate a fact about structure or about
how the analyst's knowledge was written down?

Registered in ``results/axis_robustness_real/PREREGISTRATION.md`` Appendix D,
**with its directional prediction, before it was run**. It is a diagnostic, not a
stratum: it produces no τ, it enters no anchor table, and it is labelled post-hoc
wherever it is reported.

The question. On the committed real corpus the flip-arm contradiction rate at
depth 1 is ~0.013; on synthetic structure it was 0.627. Two things differ at
once: the **graphs**, and the **knowledge model** --
``demo/example.py::knowledge_to_recover`` returns a *greedy-minimal generator
set*, while ``synth/knowledge.py::draw_k_true`` returns a random subset of the
CPDAG's undirected edges. This module holds the graphs fixed and moves only the
knowledge model.

Method. **Exhaustive, never sampled**: every claim in turn is reversed and
``apply_orientations`` is asked whether the result admits a consistent MPDAG. For
a set of size ``|K|`` that is exactly ``|K|`` closures, so the answer is a count
rather than an estimate and carries no Monte-Carlo error at all.
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from pathlib import Path
from typing import Any

import numpy as np

from bkrobust.benchmarks.measure import select_knowledge
from bkrobust.demo.meek import apply_orientations
from bkrobust.robustness import real_survival as rs
from bkrobust.synth.knowledge import draw_k_true

FIELDS = [
    "network", "coverage", "knowledge_model", "n_k", "n_single_reversals",
    "n_contradictory", "contradiction_rate", "n_undirected_cpdag",
    "g0_is_true_dag", "g0_undirected_edges",
]


def exhaustive_single_reversal(cpdag: Any, k: list[tuple[str, str]]) -> int:
    """Count the single-claim reversals of ``k`` that Meek closure rejects.

    Exhaustive: every claim is reversed in turn. No sampling, so the result is a
    count and not an estimate.

    Args:
        cpdag: The oracle CPDAG.
        k: The claim set, as ``(tail, head)`` pairs.

    Returns:
        The number of reversals for which ``apply_orientations`` returns
        ``None``.
    """
    k = sorted(k)
    n = 0
    for i in range(len(k)):
        perturbed = [(b, a) if j == i else (a, b) for j, (a, b) in enumerate(k)]
        if apply_orientations(cpdag, perturbed) is None:
            n += 1
    return n


def rows_for_network(network: str, coverages: list[float]) -> list[dict[str, Any]]:
    """Both knowledge models for one network.

    The minimal-generator model is evaluated at every coverage the frame carries
    for this network, because that is what the main sweep uses. The
    all-undirected-edges model has no coverage parameter -- it asserts every
    undirected edge -- so it is evaluated once and labelled ``coverage = 1.0``.

    Args:
        network: Network name.
        coverages: The coverage levels the frame carries for this network.

    Returns:
        One row per (knowledge model, coverage).
    """
    dag, cpdag = rs.load_network(network)
    n_und = len(cpdag.undirected_edges)
    out: list[dict[str, Any]] = []
    models: dict[tuple[str, float], list[tuple[str, str]]] = {
        ("minimal_generator", cov): sorted(select_knowledge(dag, cpdag, cov))
        for cov in coverages
    }
    models[("all_undirected_edges", 1.0)] = sorted(
        draw_k_true(dag, cpdag, np.random.default_rng(rs.derived_seed(network, "kfull")), 1.0)
    )
    for (name, coverage), k in sorted(models.items()):
        g0 = apply_orientations(cpdag, list(k))
        # Computed once and reused: the exhaustive scan is |K| Meek closures, and
        # on pathfinder one closure costs about a second.
        n_contra = exhaustive_single_reversal(cpdag, list(k))
        out.append({
            "network": network,
            "coverage": coverage,
            "knowledge_model": name,
            "n_k": len(k),
            "n_single_reversals": len(k),
            "n_contradictory": n_contra,
            "contradiction_rate": (n_contra / len(k) if k else None),
            "n_undirected_cpdag": n_und,
            "g0_is_true_dag": (
                g0 is not None and g0.directed_edges == dag.directed_edges
            ),
            "g0_undirected_edges": len(g0.undirected_edges) if g0 is not None else None,
        })
    return out


def main(argv: list[str] | None = None) -> int:
    """Entry point.

    Args:
        argv: Argument vector, or ``None`` for ``sys.argv``.

    Returns:
        Process exit status.
    """
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", default="results/axis_robustness_real")
    p.add_argument("--frame", default="results/axis_robustness_real/frame.jsonl")
    args = p.parse_args(argv)

    frame = [json.loads(l) for l in Path(args.frame).read_text().splitlines() if l.strip()]
    networks = sorted({r["network"] for r in frame})
    cov_by_net: dict[str, list[float]] = {}
    for r in frame:
        cov_by_net.setdefault(r["network"], [])
        if r["coverage"] not in cov_by_net[r["network"]]:
            cov_by_net[r["network"]].append(r["coverage"])

    rows: list[dict[str, Any]] = []
    for net in networks:
        rows.extend(rows_for_network(net, sorted(cov_by_net[net])))

    out_dir = Path(args.out_dir)
    path = out_dir / "knowledge_model_diagnostic.csv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    summary: dict[str, Any] = {"n_networks": len(networks)}
    for model, label in (("minimal_generator", "minimal_generator_all_coverages"),
                         ("minimal_generator_cov1", "minimal_generator_coverage_1.0"),
                         ("all_undirected_edges", "all_undirected_edges")):
        if model == "minimal_generator_cov1":
            sel = [r for r in rows if r["knowledge_model"] == "minimal_generator"
                   and r["coverage"] == 1.0 and r["n_k"] > 0]
        else:
            sel = [r for r in rows if r["knowledge_model"] == model and r["n_k"] > 0]
        tot_k = sum(r["n_k"] for r in sel)
        tot_c = sum(r["n_contradictory"] for r in sel)
        summary[label] = {
            "n_cells_with_claims": len(sel),
            "total_single_reversals": tot_k,
            "total_contradictory": tot_c,
            "pooled_contradiction_rate": (tot_c / tot_k) if tot_k else None,
            "network_weighted_rate": (
                statistics.mean(r["contradiction_rate"] for r in sel) if sel else None
            ),
            "median_over_networks": (
                statistics.median(r["contradiction_rate"] for r in sel) if sel else None
            ),
            "cells_at_exactly_zero": sum(1 for r in sel if r["contradiction_rate"] == 0.0),
            "median_n_k": statistics.median(r["n_k"] for r in sel) if sel else None,
        }
    (out_dir / "knowledge_model_diagnostic.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
