"""Functions E24 needs from E5, copied verbatim so that E24 never imports another experiment folder.

Sources (copied 2026-11-06):
    experiments/rebuttal/E5_learned_cpdags_by_sample_size/vendor_bkrobust.py
        sha256 20133f9aaf103e27006616698b84df0c4390180e5ad6ce67da9987d7e81d06e3
        -> is_valid_adjustment_set_gac_dag (+ _d_separated), CommittedSet, possible_ancestors,
           causal_nodes_poly, committed_adjustment_set, sha, graph_from_json
           (those were themselves copied from bkrobust (main branch), see that file's header)
    experiments/rebuttal/E5_learned_cpdags_by_sample_size/run.py
        sha256 0e8d66763a8434f13ae2fe7f34c9e93ed821fc307145bd9088e4f08125d99df6
        -> canon, greedy_consistent, truthful, file_sha, repo_commit
Everything else is imported from the submitted tree (src/bkrobust).
"""

from __future__ import annotations

import hashlib
import subprocess
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import networkx as nx

from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.gac.mpdag_level import gac_forbidden_set, is_gac_valid_mpdag
from bkrobust.mpdag_criterion.criterion import is_amenable
from bkrobust.mpdag_criterion.paths import possible_descendants

Node = str


# --- E5 vendor_bkrobust.py (benchmarks/audit.py) ----------------------------------


def _d_separated(g: nx.DiGraph, x: set, y: set, z: set) -> bool:
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


# --- E5 vendor_bkrobust.py (mpdag_criterion/optimal.py) ---------------------------


@dataclass(frozen=True)
class CommittedSet:
    z: frozenset[Node] | None
    verdict: str
    causal_nodes: frozenset[Node]
    forbidden: frozenset[Node]


def possible_ancestors(g: MPDAG, node: Node) -> set[Node]:
    return {u for u in g.nodes if node in possible_descendants(g, u)} | {node}


def causal_nodes_poly(g: MPDAG, x: Node, y: Node) -> set[Node]:
    down = possible_descendants(g, x)
    return {v for v in down if v != x and (v == y or y in possible_descendants(g, v))}


def committed_adjustment_set(g: MPDAG, x: Node, y: Node) -> CommittedSet:
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


# --- E5 vendor_bkrobust.py (ledger_sweep.py, pest_panel.py) -----------------------


def sha(g: MPDAG) -> str:
    return hashlib.sha256(g.edge_string().encode()).hexdigest()[:16]


def graph_from_json(e: dict) -> MPDAG:
    return MPDAG(
        e["nodes"],
        directed=[tuple(x) for x in e["directed"]],
        undirected=[tuple(x) for x in e["undirected"]],
    )


# --- E5 run.py --------------------------------------------------------------------


def canon(e) -> tuple[str, str]:
    return tuple(sorted(e))


def greedy_consistent(g: MPDAG, claims: list) -> tuple[list, int]:
    """Keep claims in order, skipping any whose addition makes the Meek closure fail."""
    kept, dropped = [], 0
    for e in claims:
        if apply_orientations(g, [*kept, e]) is None:
            dropped += 1
        else:
            kept.append(e)
    return kept, dropped


def truthful(dag: MPDAG, edges) -> list:
    return [(a, b) if dag.is_directed_edge(a, b) else (b, a) for (a, b) in edges if dag.has_edge(a, b)]


def file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repo_commit(repo: Path) -> str:
    try:
        out = subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "HEAD"], capture_output=True, text=True, timeout=20
        )
        return out.stdout.strip() or "unknown"
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"
