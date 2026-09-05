"""The exhaustive differential sweep of the graphical criterion against the oracle.

Why this is a module and not a notebook cell: the claim the criterion rests on is
"it agrees with :func:`bkrobust.core.oracle.is_valid` everywhere in the scope we
can enumerate", and a claim like that is only worth what its *stated scope* is
worth. This script fixes the scope in code, writes it into the results file
alongside the counts, and can be re-run to reproduce them.

Scope, exhaustive within it:

* every CPDAG on 3 and on 4 labelled nodes that has at least one undirected edge
  (:func:`bkrobust.search.conjecture_study.all_cpdags`), together with
* every element of that CPDAG's corrected perturbation space
  (:func:`bkrobust.search.space_fixed.build_space_correct`), deduplicated across
  CPDAGs since spaces overlap;
* every ordered pair ``(x, y)`` with ``x != y``;
* every subset ``z`` of the remaining nodes of size at most 3 (which, at
  ``n <= 4``, is every subset).

Two strata are reported separately, because "0 disagreements" is only
informative if the hard cases were actually present: the non-amenable graphs,
where the criterion refuses every ``z`` outright, and the graphs whose ``[G]`` is
empty, where the oracle's convention rather than the graph theory decides the
answer.

Run with::

    PYTHONPATH=src python3 -m bkrobust.mpdag_criterion.sweep

Written to run on Python 3.9 as well as the repository's 3.11 target; ``zip`` is
never called with ``strict=`` (3.10+).
"""

from __future__ import annotations

import itertools
import json
import time
from collections import Counter
from pathlib import Path as FilePath
from typing import Any

from bkrobust.core.oracle import extensions, is_valid
from bkrobust.demo.graph import MPDAG
from bkrobust.mpdag_criterion.criterion import is_amenable, why_invalid
from bkrobust.search.conjecture_study import all_cpdags
from bkrobust.search.space_fixed import build_space_correct

Node = str

#: Where the sweep's numbers are written, relative to the repository root.
DEFAULT_OUT = FilePath("results/axisb4/criterion_agreement.json")

SCOPE = (
    "Exhaustive. Graphs: every CPDAG on 3 and 4 labelled nodes with at least one "
    "undirected edge, plus every element of each such CPDAG's corrected space "
    "(build_space_correct), deduplicated across CPDAGs. Queries: every ordered "
    "(x, y) with x != y, and every subset z of the remaining nodes of size <= 3 "
    "(at n <= 4 that is every subset). Compared against "
    "bkrobust.core.oracle.is_valid, which enumerates [G]."
)


def sweep_graphs(sizes: tuple[int, ...] = (3, 4)) -> list[MPDAG]:
    """Every graph in the sweep's scope, deduplicated and deterministically ordered.

    Args:
        sizes: Node counts to enumerate CPDAGs for.

    Returns:
        The CPDAGs and their corrected-space elements, sorted by
        ``(number of nodes, edge string)``.
    """
    seen: dict[tuple[tuple[Node, ...], str], MPDAG] = {}
    for n in sizes:
        for cpdag in all_cpdags(n):
            if not cpdag.undirected_edges:
                continue
            seen[(cpdag.nodes, cpdag.edge_string())] = cpdag
            for element in build_space_correct(cpdag).elements:
                seen[(element.nodes, element.edge_string())] = element
    return [seen[k] for k in sorted(seen)]


def candidate_sets(nodes: tuple[Node, ...], x: Node, y: Node, max_size: int) -> list[frozenset]:
    """Every subset of ``nodes - {x, y}`` of size at most ``max_size``.

    Args:
        nodes: The graph's vertices.
        x: The treatment node.
        y: The outcome node.
        max_size: Largest subset size to emit.

    Returns:
        The subsets, ordered by size then sorted labels.
    """
    rest = sorted(set(nodes) - {x, y})
    out: list[frozenset] = []
    for size in range(min(max_size, len(rest)) + 1):
        for combo in itertools.combinations(rest, size):
            out.append(frozenset(combo))
    return out


def run_sweep(sizes: tuple[int, ...] = (3, 4), max_z: int = 3) -> dict[str, Any]:
    """Compare the graphical criterion to the enumerating oracle over the whole scope.

    Args:
        sizes: Node counts to enumerate CPDAGs for.
        max_z: Largest candidate adjustment set to test.

    Returns:
        A JSON-serialisable record of the scope, the counts, the two reported
        strata, the distribution of failure reasons, and up to five minimised
        disagreements (smallest graph, then smallest ``z``).
    """
    started = time.time()
    graphs = sweep_graphs(sizes)

    n_cases = n_agree = n_disagree = 0
    non_amenable_cases = non_amenable_agree = 0
    empty_extension_graphs = 0
    empty_extension_cases = empty_extension_agree = 0
    reasons: Counter = Counter()
    disagreements: list[dict[str, Any]] = []

    for g in graphs:
        has_ext = bool(extensions(g))
        if not has_ext:
            empty_extension_graphs += 1
        for x, y in itertools.permutations(sorted(g.nodes), 2):
            amenable = is_amenable(g, x, y)
            for z in candidate_sets(g.nodes, x, y, max_z):
                reason = why_invalid(g, x, y, z)
                mine = reason == ""
                truth = is_valid(z, g, x, y)
                n_cases += 1
                reasons[reason or "valid"] += 1
                if not amenable:
                    non_amenable_cases += 1
                if not has_ext:
                    empty_extension_cases += 1
                if mine == truth:
                    n_agree += 1
                    if not amenable:
                        non_amenable_agree += 1
                    if not has_ext:
                        empty_extension_agree += 1
                else:
                    n_disagree += 1
                    disagreements.append(
                        {
                            "graph": g.edge_string(),
                            "nodes": list(g.nodes),
                            "x": x,
                            "y": y,
                            "z": sorted(z),
                            "oracle": truth,
                            "criterion": mine,
                            "reason": reason,
                            "amenable": amenable,
                            "n_extensions": len(extensions(g)),
                        }
                    )

    disagreements.sort(key=lambda d: (len(d["nodes"]), len(d["graph"]), len(d["z"]), d["graph"]))
    return {
        "scope": SCOPE,
        "sizes": list(sizes),
        "max_z": max_z,
        "n_graphs": len(graphs),
        "n_cases": n_cases,
        "n_agree": n_agree,
        "n_disagree": n_disagree,
        "non_amenable_cases": non_amenable_cases,
        "non_amenable_agree": non_amenable_agree,
        "empty_extension_graphs": empty_extension_graphs,
        "empty_extension_cases": empty_extension_cases,
        "empty_extension_agree": empty_extension_agree,
        "reason_counts": dict(sorted(reasons.items())),
        "disagreements": disagreements[:5],
        "seconds": round(time.time() - started, 2),
    }


N5_SCOPE = (
    "Out-of-scope stress pass, not part of the headline numbers. Every 5th CPDAG "
    "on 5 labelled nodes that has between 1 and 6 undirected edges, plus every "
    "element of each one's corrected space; every ordered (x, y) with x != y; "
    "every z of size <= 3. Sampled deterministically by list slicing, not by any "
    "RNG."
)


def run_n5_sample(step: int = 5, max_undirected: int = 6, max_z: int = 3) -> dict[str, Any]:
    """A larger, deterministically sampled stress pass at n=5.

    The n<=4 sweep is exhaustive but small, and the two definitions that had to be
    repaired (``possde`` and amenability) both only misbehave on graphs dense
    enough to have shielded paths. This pass buys roughly thirty times the case
    count at the price of exhaustiveness, and its sampling is a fixed list slice
    so it reproduces exactly.

    Args:
        step: Take every ``step``-th CPDAG.
        max_undirected: Skip CPDAGs with more undirected edges than this, whose
            spaces and extension sets explode.
        max_z: Largest candidate adjustment set to test.

    Returns:
        A JSON-serialisable record with the same counts as :func:`run_sweep`.
    """
    started = time.time()
    cpdags = [
        c for c in all_cpdags(5) if c.undirected_edges and len(c.undirected_edges) <= max_undirected
    ]
    seen: dict[tuple[tuple[Node, ...], str], MPDAG] = {}
    for cpdag in cpdags[::step]:
        seen[(cpdag.nodes, cpdag.edge_string())] = cpdag
        for element in build_space_correct(cpdag).elements:
            seen[(element.nodes, element.edge_string())] = element
    graphs = [seen[k] for k in sorted(seen)]

    n_cases = n_agree = n_disagree = 0
    non_amenable_cases = non_amenable_agree = 0
    for g in graphs:
        for x, y in itertools.permutations(sorted(g.nodes), 2):
            amenable = is_amenable(g, x, y)
            for z in candidate_sets(g.nodes, x, y, max_z):
                mine = why_invalid(g, x, y, z) == ""
                truth = is_valid(z, g, x, y)
                n_cases += 1
                if not amenable:
                    non_amenable_cases += 1
                if mine == truth:
                    n_agree += 1
                    if not amenable:
                        non_amenable_agree += 1
                else:
                    n_disagree += 1
    return {
        "scope": N5_SCOPE,
        "step": step,
        "n_graphs": len(graphs),
        "n_cases": n_cases,
        "n_agree": n_agree,
        "n_disagree": n_disagree,
        "non_amenable_cases": non_amenable_cases,
        "non_amenable_agree": non_amenable_agree,
        "seconds": round(time.time() - started, 2),
    }


def main(out_path: str | FilePath = DEFAULT_OUT, *, with_n5: bool = True) -> dict[str, Any]:
    """Run the sweep and write its record to ``out_path``.

    Args:
        out_path: Destination JSON file. Parent directories are created.
        with_n5: Also run :func:`run_n5_sample` and record it under
            ``"n5_stress_sample"``. Adds roughly a minute and a half.

    Returns:
        The record that was written.
    """
    record = run_sweep()
    if with_n5:
        record["n5_stress_sample"] = run_n5_sample()
    path = FilePath(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in record.items() if k != "scope"}, indent=2, sort_keys=True))
    return record


if __name__ == "__main__":
    main()
