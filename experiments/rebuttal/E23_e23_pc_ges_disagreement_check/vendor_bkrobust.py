# Copied verbatim from experiments/rebuttal/E5_learned_cpdags_by_sample_size/vendor_bkrobust.py
# (sha256 20133f9aaf103e27006616698b84df0c4390180e5ad6ce67da9987d7e81d06e3, copied 2026-11-06). Its own header follows.
#
"""Functions E5 needs that the submitted tree lacks, copied from the second author's tree.

Source repository: bkrobust (main branch)
    commit ba7ddc12c446f1087e58dc428720c84d33b1b622
Code is copied verbatim; long docstrings are shortened. The source files below are STAGED in that tree but not in HEAD, so the commit alone does
not identify them; their sha256 at copy time (2026-11-06) is given instead.

    src/bkrobust/benchmarks/audit.py          fccfd536d28cc48ad8afe5cc73a393c5904b1690c458188e298d0d3cb2782146
        -> is_valid_adjustment_set_gac_dag (verbatim, plus a networkx < 3.3 fallback)
    src/bkrobust/mpdag_criterion/optimal.py   538cbb6f63b4eea978f433b8b7c5b4fb63c2bf683f69ca2d6b781224cd9ed734
        -> CommittedSet, possible_ancestors, causal_nodes_poly, committed_adjustment_set (verbatim)
    experiments/ledger_sweep.py               f84eccd7fd5930bf9ab6e05158c110ca3675cb3eb62fb870a6581f8c2324ea3c
        -> sha, breaking_depths, claim_lower_bound, bucket, chance_level, audit_row and the
           constants they use; TIER is cut to the four T3 networks. instrument_row is NOT copied: run.py has its own
           version, because the stored one passes use_ladder / search_closure_budget, which the
           submitted breakdown_radius does not accept.
    experiments/pest_panel.py                 2a18c4e51b09942f87df30422ad3f33ee931912e5b5c22e7da57b77315e2a624
        -> graph_from_json (verbatim)
    src/bkrobust/core/conventions.py          9e923d5a86244445564fee79caeea35f52f0f3c80cfc4eaaf455f88eb44d23ce
        -> NO_RETRACTABLE_EDGES (the constant only)

Every other dependency is imported from the submitted tree (src/bkrobust on the rebuttal branch).
gac/mpdag_level.py, demo/meek.py, demo/graph.py, demo/evaluate.py, demo/example.py,
mpdag_criterion/criterion.py and paths.py are byte-identical between the two trees; hybrid.py,
benchmarks/measure.py, search/exact.py and core/conventions.py differ.
"""

from __future__ import annotations

import hashlib
import itertools
import math
from collections.abc import Iterable
from dataclasses import dataclass

import networkx as nx
import numpy as np

from bkrobust.benchmarks.measure import component_of
from bkrobust.demo.evaluate import is_valid_adjustment_set_dag, optimal_adjustment_set_dag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.gac.mpdag_level import gac_forbidden_set, is_gac_valid_mpdag
from bkrobust.mpdag_criterion.criterion import is_amenable
from bkrobust.mpdag_criterion.paths import possible_descendants

Node = str

# --- core/conventions.py ---------------------------------------------------------

NO_RETRACTABLE_EDGES: str = "no_retractable_edges"

# --- benchmarks/audit.py ---------------------------------------------------------


def _d_separated(g: nx.DiGraph, x: set, y: set, z: set) -> bool:
    # is_d_separator appeared in networkx 3.3; d_separated is the older name of the same test
    if hasattr(nx, "is_d_separator"):
        return bool(nx.is_d_separator(g, x, y, z))
    return bool(nx.d_separated(g, x, y, z))


def is_valid_adjustment_set_gac_dag(dag: MPDAG, x: str, y: str, z: Iterable[str]) -> bool:
    """Whether ``z`` satisfies the generalised adjustment criterion in ``dag`` (networkx only)."""
    zs = frozenset(z)
    if x in zs or y in zs:
        return False
    d = nx.DiGraph()
    d.add_nodes_from(dag.nodes)
    d.add_edges_from(dag.directed_edges)
    causal = nx.descendants(d, x) & (nx.ancestors(d, y) | {y})
    forbidden = {x}
    for c in causal:
        forbidden |= nx.descendants(d, c) | {c}
    if zs & forbidden:
        return False
    proper = d.copy()
    proper.remove_edges_from([(x, c) for c in list(d.successors(x)) if c in causal])
    return _d_separated(proper, {x}, {y}, set(zs))


# --- mpdag_criterion/optimal.py --------------------------------------------------


@dataclass(frozen=True)
class CommittedSet:
    """The set an analyst commits to on their own graph, and how it was reached."""

    z: frozenset[Node] | None
    verdict: str
    causal_nodes: frozenset[Node]
    forbidden: frozenset[Node]


def possible_ancestors(g: MPDAG, node: Node) -> set[Node]:
    """Nodes that can reach ``node`` along a possibly causal path, plus ``node``."""
    return {u for u in g.nodes if node in possible_descendants(g, u)} | {node}


def causal_nodes_poly(g: MPDAG, x: Node, y: Node) -> set[Node]:
    """Nodes on a possibly causal path from ``x`` to ``y``, by reachability."""
    down = possible_descendants(g, x)
    return {v for v in down if v != x and (v == y or y in possible_descendants(g, v))}


def committed_adjustment_set(g: MPDAG, x: Node, y: Node) -> CommittedSet:
    """The set the analyst commits to on ``g``, polynomial, with its provenance."""
    if y not in possible_descendants(g, x):
        return CommittedSet(None, "y_not_possible_descendant", frozenset(), frozenset())
    if not is_amenable(g, x, y):
        return CommittedSet(None, "not_amenable", frozenset(), frozenset())
    cn = causal_nodes_poly(g, x, y)
    forb = gac_forbidden_set(g, x, y)
    parents: set[Node] = set()
    for v in cn:
        parents |= g.parents(v)
    z = frozenset(parents - forb - cn - {x})
    if is_gac_valid_mpdag(g, x, y, z):
        return CommittedSet(z, "closed_form", frozenset(cn), frozenset(forb))
    canonical = frozenset((possible_ancestors(g, x) | possible_ancestors(g, y)) - forb - {x, y})
    if is_gac_valid_mpdag(g, x, y, canonical):
        return CommittedSet(canonical, "canonical", frozenset(cn), frozenset(forb))
    return CommittedSet(None, "no_valid_set", frozenset(cn), frozenset(forb))


# --- experiments/ledger_sweep.py -------------------------------------------------

MAX_DEPTH = 3
SUBSET_BUDGET = 300
EXHAUSTIVE_MAX_RETRACTABLE = 12
LADDER_MAX_NODES = 60
TIER = {n: "T3" for n in ("arth150", "ecoli70", "magic-irri", "magic-niab")}
RUNGS = ("B2_free_rule", "B4_kg0", "B5_k", "B6_hop", "B7_claim")


def sha(g: MPDAG) -> str:
    """A hash of a graph's edge string, for the provenance stamp."""
    return hashlib.sha256(g.edge_string().encode()).hexdigest()[:16]


def breaking_depths(
    cpdag: MPDAG, k: list, x: str, y: str, z: frozenset
) -> tuple[int | None, str, int | None, tuple]:
    """Smallest retraction breaking validity (with its witness) and identifiability."""
    r_claim = r_id = None
    witness: tuple = ()
    status = f"gt_{min(len(k), MAX_DEPTH)}" if k else NO_RETRACTABLE_EDGES
    used = 0
    for depth in range(1, MAX_DEPTH + 1):
        if depth > len(k):
            break
        if used + math.comb(len(k), depth) > SUBSET_BUDGET:
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
                r_claim, status, witness = depth, "exact", tuple(k[i] for i in subset)
            if r_id is None and not is_amenable(g, x, y):
                r_id = depth
            if r_claim is not None and r_id is not None:
                break
        if r_claim is not None and r_id is not None:
            break
    if r_claim is None and status.startswith("gt_") and len(k) <= MAX_DEPTH:
        status = "unreached"  # every non-empty subset was tried
    return r_claim, status, r_id, witness


def claim_lower_bound(r_claim: int | None, status: str) -> int | None:
    """What the claim instrument certifies: the radius, or its censored floor."""
    if r_claim is not None:
        return r_claim
    if status in ("unreached", NO_RETRACTABLE_EDGES):
        return None
    if status.startswith("censored_at_depth_"):
        return int(status.rsplit("_", 1)[1])
    if status.startswith("gt_"):
        return int(status.rsplit("_", 1)[1]) + 1
    raise ValueError(status)


def bucket(r: int | None, w_c: int, z_valid: bool) -> str:
    """The four-way audit verdict for one instrument's certified radius."""
    covered = True if r is None else w_c <= r - 1
    if covered:
        return "dangerous" if not z_valid else "held"
    return "safe" if not z_valid else "slack"


def chance_level(
    cpdag: MPDAG, base: list, rng: np.random.Generator, retries: int = 20
) -> list | None:
    """Same edges as ``base``, orientations uniform among Meek-consistent choices."""
    for _ in range(retries):
        drawn = [(a, b) if rng.integers(2) else (b, a) for (a, b) in base]
        if apply_orientations(cpdag, drawn) is not None:
            return drawn
    return None


def audit_row(dag: MPDAG, cpdag: MPDAG, k: list, x: str, y: str, z: frozenset, cert: dict) -> dict:
    """The truth, revealed after the freeze: commission, omission, validity, buckets."""
    rec: dict = {}
    comp = component_of(cpdag, x)
    w_c = sum(1 for (a, b) in k if not dag.is_directed_edge(a, b))
    comp_edges = {tuple(sorted(e)) for e in cpdag.undirected_edges if comp and e[0] in comp}
    asserted = {tuple(sorted(e)) for e in k}
    rec["n_omitted"] = len(comp_edges - asserted)
    rec["n_wrong_claims"] = w_c
    z_ok = is_valid_adjustment_set_gac_dag(dag, x, y, z)
    rec["z_valid_at_truth"] = int(z_ok)
    rec["z_backdoor_at_truth"] = int(is_valid_adjustment_set_dag(dag, x, y, z))
    if not z_ok:
        rec["severity"] = "validity"
    elif set(z) == set(optimal_adjustment_set_dag(dag, x, y)):
        rec["severity"] = "benign"
    else:
        rec["severity"] = "efficiency"
    for rung in RUNGS:
        rec[f"bucket_{rung}"] = bucket(cert[rung], w_c, z_ok)
    return rec


# --- experiments/pest_panel.py ---------------------------------------------------


def graph_from_json(e: dict) -> MPDAG:
    """An MPDAG from the ``nodes``/``directed``/``undirected`` fields of a record."""
    return MPDAG(
        e["nodes"],
        directed=[tuple(x) for x in e["directed"]],
        undirected=[tuple(x) for x in e["undirected"]],
    )
