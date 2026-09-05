"""The differential sweep behind every number claimed for :mod:`bkrobust.gac`.

Why this is a module and not a notebook cell: the claims this package rests on
are all of the form "it agrees with *that* implementation everywhere in the
scope we can enumerate", and a claim like that is worth exactly what its stated
scope is worth. This script fixes each scope in code, writes it into the results
file next to the counts, and can be re-run to reproduce them.

Four checks, matching the four sections of ``results/axisa2/gac_agreement.json``:

1. **DAG level, GAC vs back-door, classified.** Over exhaustive small DAGs,
   :func:`~bkrobust.gac.dag_level.is_gac_valid_dag` is compared with
   :func:`bkrobust.demo.evaluate.is_valid_adjustment_set_dag`. Disagreement is
   *expected* -- the GAC is strictly weaker -- so the informative output is not
   the count but the **classification**: every disagreement must be a case where
   back-door rejects and the GAC accepts because some member of ``z`` is a
   descendant of ``x`` that is not on or below a proper causal path to ``y``.
   Anything else would be a bug, and is counted separately as unexplainable.
   The one-sided implication (back-door valid => GAC valid) is checked on every
   case.
2. **MPDAG level vs enumeration -- the semantic anchor.** Exhaustive over the
   same graph scope :mod:`bkrobust.mpdag_criterion.sweep` uses,
   :func:`~bkrobust.gac.mpdag_level.is_gac_valid_mpdag` is compared with
   ``all(is_gac_valid_dag(d, ...) for d in extensions(g))``. Disagreements must
   be **zero**; the non-amenable stratum is reported separately, because a
   "0 disagreements" headline is only informative if the hard cases were
   present.
3. **The O-set claim, tested rather than assumed.** Most prior results in this
   repository fix ``Z = O(G0)``, the optimal adjustment set. If the GAC and
   back-door never disagree *on that particular ``Z``*, then switching criteria
   changes nothing for those results -- which is a finding either way, so it is
   measured.
4. **Performance.** Cold-cache per-call timings of ``is_gac_valid_mpdag`` on
   dense Erdos-Renyi instances built exactly as :mod:`bkrobust.hybrid` builds
   them. The bar this package was written against is under 5 ms per call at
   ``n = 20``.

Run with::

    PYTHONPATH=src python3 -m bkrobust.gac.sweep

Written to run on Python 3.9 as well as the repository's 3.11 target: neither
``zip(strict=)`` nor ``itertools.pairwise`` (3.10+) is used, and randomness is
drawn only from explicitly passed generators.
"""

from __future__ import annotations

import itertools
import json
import platform
import statistics
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from bkrobust.core.oracle import extensions
from bkrobust.demo.evaluate import (
    is_valid_adjustment_set_dag,
    optimal_adjustment_set_mpdag,
)
from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.gac.dag_level import forbidden_set_dag, is_gac_valid_dag
from bkrobust.gac.mpdag_level import clear_cache, is_gac_valid_mpdag, why_invalid_gac
from bkrobust.mpdag_criterion.criterion import is_amenable, is_valid_mpdag
from bkrobust.search.conjecture_study import all_cpdags, all_dags
from bkrobust.search.space_fixed import build_space_correct
from bkrobust.synth.generators import erdos_renyi_dag

Node = str
Edge = tuple[str, str]

#: Where the sweep's numbers are written, relative to the repository root.
DEFAULT_OUT = Path("results/axisa2/gac_agreement.json")

CHECK1_SCOPE = (
    "Graphs: every labelled DAG on 4 nodes (all_dags(4)) and {five} of the "
    "labelled DAGs on 5 nodes (all_dags(5)[::{step5}]); when the step is 1 that "
    "is exhaustive, and any larger step is a deterministic list slice, never an "
    "RNG sample. Queries: every ordered (x, y) with x != y, and every subset z "
    "of the remaining nodes of size <= 3 (at n <= 5 with x and y removed that is "
    "every subset). Compares bkrobust.gac.dag_level.is_gac_valid_dag with "
    "bkrobust.demo.evaluate.is_valid_adjustment_set_dag (Pearl back-door)."
)

CHECK2_SCOPE = (
    "Exhaustive. Graphs: every CPDAG on 3 and 4 labelled nodes with at least one "
    "undirected edge, plus every element of each such CPDAG's corrected space "
    "(bkrobust.search.space_fixed.build_space_correct), deduplicated across "
    "CPDAGs -- the same graph scope as bkrobust.mpdag_criterion.sweep. Queries: "
    "every ordered (x, y) with x != y, and every subset z of the remaining nodes "
    "of size <= 3 (at n <= 4 that is every subset). Compares "
    "bkrobust.gac.mpdag_level.is_gac_valid_mpdag with the enumerated semantics "
    "all(is_gac_valid_dag(d, x, y, z) for d in extensions(g)), False when [g] is "
    "empty."
)

CHECK3_SCOPE = (
    "Same graph scope as check 2, but Z is not swept: for each CPDAG G0 and each "
    "ordered (x, y), Z = optimal_adjustment_set_mpdag(G0, x, y) is fixed once "
    "(skipped when it is None, i.e. when the extensions of G0 disagree on the "
    "optimal set) and then evaluated at G0 itself and at every element of G0's "
    "corrected space -- mirroring how the repository actually uses O(G0). At each "
    "graph, is_gac_valid_mpdag is compared with "
    "bkrobust.mpdag_criterion.is_valid_mpdag (the back-door criterion)."
)

CHECK4_SCOPE = (
    "G = apply_orientations(dag_to_cpdag(erdos_renyi_dag(n, rng, p)), "
    "sorted(knowledge_to_recover(dag, cpdag))), exactly as bkrobust.hybrid builds "
    "its instances; seeds 0-2 per (n, p); n in {12, 20, 30}, p in {0.70, 0.85}. "
    "The 'half_knowledge' rows keep only the first half of that sorted knowledge "
    "list, so the graph stays partially oriented and amenability has real work to "
    "do. Every call is preceded by clear_cache(), so each one pays for its own "
    "Meek closure, adjacency tables and possible-descendant sets -- the "
    "pessimistic figure, and the relevant one when a search walks distinct graphs."
)


# --------------------------------------------------------------------------------
# Shared scope helpers
# --------------------------------------------------------------------------------


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


def mpdag_scope(sizes: tuple[int, ...] = (3, 4)) -> list[MPDAG]:
    """Every graph in the MPDAG scope, deduplicated and deterministically ordered.

    Args:
        sizes: Node counts to enumerate CPDAGs for.

    Returns:
        The CPDAGs with at least one undirected edge and their corrected-space
        elements, sorted by ``(nodes, edge string)``.
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


# --------------------------------------------------------------------------------
# Check 1: DAG level, GAC versus back-door
# --------------------------------------------------------------------------------


def classify_disagreement(dag: MPDAG, x: Node, y: Node, z: frozenset) -> str:
    """Name the reason the GAC and back-door differ on one DAG-level case.

    The only legitimate difference is the one the GAC was introduced to make:
    back-door bars every descendant of ``x``, the GAC bars only what is on or
    below a proper causal path from ``x`` to ``y``. So a legitimate disagreement
    has the GAC accepting, back-door rejecting, and some member of ``z`` a
    descendant of ``x`` outside ``forb(x, y, dag)``.

    Args:
        dag: A fully oriented DAG.
        x: The treatment node.
        y: The outcome node.
        z: The candidate adjustment set the two criteria disagreed on.

    Returns:
        ``"z_has_descendant_of_x_off_the_causal_route"`` for the expected class,
        ``"backdoor_valid_but_gac_invalid"`` if the one-sided implication broke,
        or ``"unexplained"``.
    """
    gac = is_gac_valid_dag(dag, x, y, z)
    if not gac:
        return "backdoor_valid_but_gac_invalid"
    off_route = (z & dag.descendants(x)) - forbidden_set_dag(dag, x, y)
    if off_route:
        return "z_has_descendant_of_x_off_the_causal_route"
    return "unexplained"


def run_check1(step5: int = 1, max_z: int = 3) -> dict[str, Any]:
    """Compare the DAG-level GAC with the back-door criterion and classify every difference.

    Args:
        step5: Take every ``step5``-th DAG on 5 nodes. ``1`` is exhaustive and is
            the default; larger steps exist so the test suite can run a bounded
            slice of the same code.
        max_z: Largest candidate adjustment set to test.

    Returns:
        A JSON-serialisable record of the scope, the counts, the disagreement
        classification and up to five minimised examples.
    """
    started = time.time()
    graphs = list(all_dags(4)) + list(all_dags(5))[::step5]

    n_cases = n_agree = n_disagree = 0
    classes: Counter = Counter()
    implication_failures = 0
    examples: list[dict[str, Any]] = []

    for dag in graphs:
        for x, y in itertools.permutations(sorted(dag.nodes), 2):
            for z in candidate_sets(dag.nodes, x, y, max_z):
                gac = is_gac_valid_dag(dag, x, y, z)
                backdoor = is_valid_adjustment_set_dag(dag, x, y, z)
                n_cases += 1
                if backdoor and not gac:
                    implication_failures += 1
                if gac == backdoor:
                    n_agree += 1
                    continue
                n_disagree += 1
                kind = classify_disagreement(dag, x, y, z)
                classes[kind] += 1
                if kind != "z_has_descendant_of_x_off_the_causal_route" or len(examples) < 5:
                    examples.append(
                        {
                            "graph": dag.edge_string(),
                            "n_nodes": len(dag.nodes),
                            "x": x,
                            "y": y,
                            "z": sorted(z),
                            "gac": gac,
                            "backdoor": backdoor,
                            "class": kind,
                        }
                    )

    return {
        "scope": CHECK1_SCOPE.format(
            step5=step5, five="every one" if step5 == 1 else f"every {step5}th"
        ),
        "n_graphs": len(graphs),
        "n_cases": n_cases,
        "n_agree": n_agree,
        "n_disagree": n_disagree,
        "disagreement_classes": dict(sorted(classes.items())),
        "n_unexplainable_disagreements": classes["unexplained"]
        + classes["backdoor_valid_but_gac_invalid"],
        "backdoor_valid_implies_gac_valid": implication_failures == 0,
        "implication_failures": implication_failures,
        "examples": examples[:5],
        "seconds": round(time.time() - started, 2),
    }


# --------------------------------------------------------------------------------
# Check 2: the semantic anchor
# --------------------------------------------------------------------------------


def enumerated_gac_semantics(g: MPDAG, x: Node, y: Node, z: frozenset) -> bool:
    """The predicate :func:`~bkrobust.gac.mpdag_level.is_gac_valid_mpdag` must equal.

    ``z`` is GAC-valid in every DAG extension of ``g``, with the repository's
    convention that an empty ``[g]`` certifies nothing.

    Args:
        g: The MPDAG.
        x: The treatment node.
        y: The outcome node.
        z: The candidate adjustment set.

    Returns:
        True iff ``[g]`` is non-empty and ``z`` is GAC-valid in each of its DAGs.
    """
    exts = extensions(g)
    if not exts:
        return False
    return all(is_gac_valid_dag(d, x, y, z) for d in exts)


def run_check2(sizes: tuple[int, ...] = (3, 4), max_z: int = 3) -> dict[str, Any]:
    """Compare the MPDAG-level GAC with the enumerated semantics over the whole scope.

    Args:
        sizes: Node counts to enumerate CPDAGs for.
        max_z: Largest candidate adjustment set to test.

    Returns:
        A JSON-serialisable record of the scope, the counts, the non-amenable
        and empty-``[G]`` strata, the reason distribution and up to five
        minimised disagreements.
    """
    started = time.time()
    graphs = mpdag_scope(sizes)

    n_cases = n_agree = n_disagree = 0
    non_amenable_cases = non_amenable_agree = 0
    empty_extension_graphs = empty_extension_cases = 0
    reasons: Counter = Counter()
    disagreements: list[dict[str, Any]] = []

    for g in graphs:
        has_ext = bool(extensions(g))
        if not has_ext:
            empty_extension_graphs += 1
        for x, y in itertools.permutations(sorted(g.nodes), 2):
            amenable = is_amenable(g, x, y)
            for z in candidate_sets(g.nodes, x, y, max_z):
                reason = why_invalid_gac(g, x, y, z)
                mine = reason == ""
                truth = enumerated_gac_semantics(g, x, y, z)
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
                else:
                    n_disagree += 1
                    disagreements.append(
                        {
                            "graph": g.edge_string(),
                            "nodes": list(g.nodes),
                            "x": x,
                            "y": y,
                            "z": sorted(z),
                            "semantics": truth,
                            "criterion": mine,
                            "reason": reason,
                            "amenable": amenable,
                            "n_extensions": len(extensions(g)),
                        }
                    )

    disagreements.sort(key=lambda d: (len(d["nodes"]), len(d["graph"]), len(d["z"]), d["graph"]))
    return {
        "scope": CHECK2_SCOPE,
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
        "reason_counts": dict(sorted(reasons.items())),
        "disagreements": disagreements[:5],
        "seconds": round(time.time() - started, 2),
    }


CHECK2_N5_SCOPE = (
    "Out-of-scope stress pass on the same comparison as check 2, not part of the "
    "headline numbers. Every {step}th CPDAG on 5 labelled nodes with between 1 and "
    "{max_undirected} undirected edges, plus every element of each one's corrected "
    "space; every ordered (x, y) with x != y; every z of size <= 3. Sampled "
    "deterministically by list slicing, never by an RNG. It exists because the one "
    "part of the MPDAG criterion that is measured rather than proved -- that "
    "gac_forbidden_set equals the union of the extensions' forbidden sets -- can "
    "only fail on graphs dense enough for two possible-descendant witnesses to "
    "need incompatible extensions, and n <= 4 is thin on those."
)


def run_check2_n5_sample(step: int = 5, max_undirected: int = 6, max_z: int = 3) -> dict[str, Any]:
    """A larger, deterministically sampled stress pass of check 2 at n=5.

    Args:
        step: Take every ``step``-th qualifying CPDAG.
        max_undirected: Skip CPDAGs with more undirected edges than this, whose
            spaces and extension sets explode.
        max_z: Largest candidate adjustment set to test.

    Returns:
        A JSON-serialisable record with the same headline counts as
        :func:`run_check2`.
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
    non_amenable_cases = 0
    for g in graphs:
        for x, y in itertools.permutations(sorted(g.nodes), 2):
            amenable = is_amenable(g, x, y)
            for z in candidate_sets(g.nodes, x, y, max_z):
                mine = why_invalid_gac(g, x, y, z) == ""
                truth = enumerated_gac_semantics(g, x, y, z)
                n_cases += 1
                if not amenable:
                    non_amenable_cases += 1
                if mine == truth:
                    n_agree += 1
                else:
                    n_disagree += 1
    return {
        "scope": CHECK2_N5_SCOPE.format(step=step, max_undirected=max_undirected),
        "step": step,
        "n_graphs": len(graphs),
        "n_cases": n_cases,
        "n_agree": n_agree,
        "n_disagree": n_disagree,
        "non_amenable_cases": non_amenable_cases,
        "seconds": round(time.time() - started, 2),
    }


# --------------------------------------------------------------------------------
# Check 3: the O-set
# --------------------------------------------------------------------------------


#: Why check 3 comes out the way it does, recorded next to the number so the
#: reader does not have to re-derive it. ``O = pa(cn) \\ (cn u {x})`` can never
#: contain a descendant of ``x``: if ``w in de(x)`` were a parent of some causal
#: node ``c``, then ``w`` would be an ancestor of ``y`` as well as a descendant
#: of ``x``, hence itself in ``cn`` and excluded. So ``O ∩ de(x) = {}`` always,
#: and on a set disjoint from ``de(x)`` the two criteria provably coincide --
#: the GAC's forbidden set is contained in ``de(x) u {x, y}``, and blocking every
#: back-door path is equivalent to blocking every proper non-causal path once
#: ``Z`` misses ``de(x)``. The measured agreement is a check on that argument,
#: not a surprise.
CHECK3_NOTE = (
    "O = pa(cn) \\ (cn u {x}) is always disjoint from de(x) -- a descendant of x "
    "that parents a causal node would itself be a causal node -- and on sets "
    "disjoint from de(x) the GAC and the back-door criterion provably coincide. "
    "The count below tests that argument rather than discovering something new; "
    "n_both_valid / n_both_invalid are reported so the reader can see the "
    "comparison is not vacuous."
)


def run_check3(sizes: tuple[int, ...] = (3, 4)) -> dict[str, Any]:
    """Ask whether the GAC and back-door ever differ on ``Z = O(G0)`` specifically.

    Args:
        sizes: Node counts to enumerate CPDAGs for.

    Returns:
        A JSON-serialisable record of the scope, how often ``O(G0)`` was
        identified at all, the case count, the disagreement count, the
        valid/invalid split and up to five examples.
    """
    started = time.time()
    n_cases = n_agree = n_disagree = 0
    n_queries = n_optimal_none = 0
    n_both_valid = n_both_invalid = 0
    examples: list[dict[str, Any]] = []

    for n in sizes:
        for cpdag in all_cpdags(n):
            if not cpdag.undirected_edges:
                continue
            elements = [cpdag, *build_space_correct(cpdag).elements]
            seen: dict[str, MPDAG] = {}
            for element in elements:
                seen.setdefault(element.edge_string(), element)
            visited = [seen[k] for k in sorted(seen)]
            for x, y in itertools.permutations(sorted(cpdag.nodes), 2):
                n_queries += 1
                optimal = optimal_adjustment_set_mpdag(cpdag, x, y)
                if optimal is None:
                    n_optimal_none += 1
                    continue
                z = frozenset(optimal)
                for g in visited:
                    gac = is_gac_valid_mpdag(g, x, y, z)
                    backdoor = is_valid_mpdag(g, x, y, z)
                    n_cases += 1
                    if gac == backdoor:
                        n_agree += 1
                        if gac:
                            n_both_valid += 1
                        else:
                            n_both_invalid += 1
                    else:
                        n_disagree += 1
                        if len(examples) < 5:
                            examples.append(
                                {
                                    "g0": cpdag.edge_string(),
                                    "graph": g.edge_string(),
                                    "x": x,
                                    "y": y,
                                    "z": sorted(z),
                                    "gac": gac,
                                    "backdoor": backdoor,
                                }
                            )

    return {
        "scope": CHECK3_SCOPE,
        "note": CHECK3_NOTE,
        "sizes": list(sizes),
        "n_queries": n_queries,
        "n_optimal_not_identified": n_optimal_none,
        "n_cases": n_cases,
        "n_agree": n_agree,
        "n_disagree": n_disagree,
        "n_both_valid": n_both_valid,
        "n_both_invalid": n_both_invalid,
        "examples": examples,
        "seconds": round(time.time() - started, 2),
    }


# --------------------------------------------------------------------------------
# Check 4: performance
# --------------------------------------------------------------------------------


def dense_instance(n: int, seed: int, edge_prob: float, *, half_knowledge: bool = False) -> MPDAG:
    """One dense Erdos-Renyi instance, built the way :mod:`bkrobust.hybrid` builds them.

    Args:
        n: Number of nodes.
        seed: Seed for the sole generator; no other randomness is used.
        edge_prob: Erdos-Renyi edge probability.
        half_knowledge: Assert only the first half of the recovering knowledge,
            leaving the graph partially oriented.

    Returns:
        The resulting MPDAG.

    Raises:
        RuntimeError: If the orientations FAIL to close, which would mean the
            knowledge extracted from the DAG is inconsistent with its own CPDAG.
    """
    rng = np.random.default_rng(seed)
    dag = erdos_renyi_dag(n, rng, edge_prob)
    cpdag = dag_to_cpdag(dag)
    knowledge = sorted(knowledge_to_recover(dag, cpdag))
    if half_knowledge:
        knowledge = knowledge[: len(knowledge) // 2]
    g = apply_orientations(cpdag, knowledge)
    if g is None:
        raise RuntimeError(f"orientations FAILed on n={n}, seed={seed}, p={edge_prob}")
    return g


def _queries(g: MPDAG, n_treatments: int = 6) -> list[tuple[Node, Node, frozenset]]:
    """A deterministic spread of ``(x, y, z)`` queries on one graph.

    Args:
        g: The graph to query.
        n_treatments: How many treatments to take from the front of the sorted
            node list.

    Returns:
        Queries covering several treatments and ``|z|`` in ``{0, 2, 4}``.
    """
    nodes = sorted(g.nodes)
    out: list[tuple[Node, Node, frozenset]] = []
    for x in nodes[:n_treatments]:
        y = nodes[-1] if x != nodes[-1] else nodes[0]
        rest = [v for v in nodes if v not in (x, y)]
        for size in (0, 2, 4):
            out.append((x, y, frozenset(rest[:size])))
    return out


def time_calls(g: MPDAG, queries: list[tuple[Node, Node, frozenset]]) -> list[float]:
    """Time :func:`~bkrobust.gac.mpdag_level.is_gac_valid_mpdag` once per query, cold.

    Every call is preceded by :func:`~bkrobust.gac.mpdag_level.clear_cache`, so
    no call benefits from a memo another call filled.

    Args:
        g: The graph to query.
        queries: The ``(x, y, z)`` triples.

    Returns:
        Milliseconds per call, in query order.
    """
    out: list[float] = []
    for x, y, z in queries:
        clear_cache()
        started = time.perf_counter()
        is_gac_valid_mpdag(g, x, y, z)
        out.append((time.perf_counter() - started) * 1000.0)
    return out


def _summarise(samples: list[float]) -> dict[str, Any]:
    """Reduce a list of per-call milliseconds to the reported statistics.

    Args:
        samples: Milliseconds per call.

    Returns:
        Call count, mean, median, p95 and worst, each rounded to 4 decimals.
    """
    ordered = sorted(samples)
    p95 = ordered[min(len(ordered) - 1, int(0.95 * len(ordered)))]
    return {
        "calls": len(samples),
        "mean_ms": round(statistics.fmean(samples), 4),
        "median_ms": round(statistics.median(samples), 4),
        "p95_ms": round(p95, 4),
        "worst_ms": round(max(samples), 4),
    }


def run_check4(
    ns: tuple[int, ...] = (12, 20, 30),
    probs: tuple[float, ...] = (0.7, 0.85),
    seeds: tuple[int, ...] = (0, 1, 2),
) -> dict[str, Any]:
    """Time ``is_gac_valid_mpdag`` cold-cache on dense instances.

    Args:
        ns: Node counts.
        probs: Erdos-Renyi edge probabilities.
        seeds: Seeds per ``(n, p)``.

    Returns:
        A JSON-serialisable record with one row per ``(n, p, knowledge)``, the
        platform, and whether the 5 ms bar at ``n = 20`` was met.
    """
    started = time.time()
    rows: list[dict[str, Any]] = []
    for n in ns:
        for p in probs:
            for half in (False, True):
                samples: list[float] = []
                edges: list[int] = []
                undirected: list[int] = []
                for seed in seeds:
                    g = dense_instance(n, seed, p, half_knowledge=half)
                    edges.append(len(g.directed_edges) + len(g.undirected_edges))
                    undirected.append(len(g.undirected_edges))
                    samples.extend(time_calls(g, _queries(g)))
                rows.append(
                    {
                        "n": n,
                        "edge_prob": p,
                        "knowledge": "half" if half else "full",
                        "mean_edges": round(statistics.fmean(edges), 1),
                        "mean_undirected_edges": round(statistics.fmean(undirected), 1),
                        "cold_cache_is_gac_valid_mpdag": _summarise(samples),
                    }
                )
    worst_at_20 = max(r["cold_cache_is_gac_valid_mpdag"]["worst_ms"] for r in rows if r["n"] == 20)
    return {
        "scope": CHECK4_SCOPE,
        "platform": {
            "system": platform.system(),
            "machine": platform.machine(),
            "python": platform.python_version(),
        },
        "rows": rows,
        "bar_ms_per_call_at_n20": 5.0,
        "worst_ms_at_n20": worst_at_20,
        "bar_met": worst_at_20 < 5.0,
        "seconds": round(time.time() - started, 2),
    }


# --------------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------------


def main(out_path: str | Path = DEFAULT_OUT, *, with_n5: bool = True) -> dict[str, Any]:
    """Run all four checks and write the record to ``out_path``.

    Args:
        out_path: Destination JSON file. Parent directories are created.
        with_n5: Also run :func:`run_check2_n5_sample`, recorded under
            ``"check2_n5_stress_sample"``. Adds roughly a minute.

    Returns:
        The record that was written.
    """
    record: dict[str, Any] = {
        "generated_by": "PYTHONPATH=src python3 -m bkrobust.gac.sweep",
        "python": sys.version.split()[0],
        "check1_dag_gac_vs_backdoor": run_check1(),
        "check2_mpdag_vs_enumeration": run_check2(),
        "check3_optimal_set": run_check3(),
        "check4_performance": run_check4(),
    }
    if with_n5:
        record["check2_n5_stress_sample"] = run_check2_n5_sample()
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps(record, indent=2, sort_keys=True))
    return record


if __name__ == "__main__":
    main()
