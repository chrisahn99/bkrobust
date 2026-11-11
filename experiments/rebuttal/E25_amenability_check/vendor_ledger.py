"""Functions vendored from the second author's working tree, which the submitted tree lacks.

Repository: bkrobust (main branch)
commit ba7ddc12c446f1087e58dc428720c84d33b1b622
Both source files below were staged but not committed at copy time ('A' in git status), so the
copy is from the working tree as it stood on 2026-11-06.

* src/bkrobust/mpdag_criterion/optimal.py: CommittedSet, possible_ancestors, causal_nodes_poly,
  committed_adjustment_set. Code copied unchanged (docstrings shortened): the released rows
  read their committed set off it.
* experiments/ledger_sweep.py: sha, chance_level, and the per-arm knowledge plan of main()
  (rewritten as rebuild_knowledge, same arithmetic, same seed derivation).
* src/bkrobust/benchmarks/audit.py: is_valid_adjustment_set_gac_dag, the release's truth oracle
  over networkx. Kept beside the submitted gac.dag_level.is_gac_valid_dag as an independent check.
  One change: nx.is_d_separator falls back to nx.d_separated on networkx < 3.3.

Every bkrobust name imported here exists in the submitted tree.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable
from dataclasses import dataclass

import networkx as nx
import numpy as np

from bkrobust.demo.example import knowledge_to_recover
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.gac.mpdag_level import gac_forbidden_set, is_gac_valid_mpdag
from bkrobust.mpdag_criterion.criterion import is_amenable
from bkrobust.mpdag_criterion.paths import possible_descendants

Node = str


# ---------------------------------------------------------------- optimal.py (unchanged)


@dataclass(frozen=True)
class CommittedSet:
    z: frozenset[Node] | None
    verdict: str
    causal_nodes: frozenset[Node]
    forbidden: frozenset[Node]


def possible_ancestors(g: MPDAG, node: Node) -> set[Node]:
    """Nodes that can reach ``node`` along a possibly causal path, plus ``node``."""
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


# ---------------------------------------------------------------- ledger_sweep.py

RELEASE_SEED = 20260913


def sha(g: MPDAG) -> str:
    """A hash of a graph's edge string, for the provenance stamp."""
    return hashlib.sha256(g.edge_string().encode()).hexdigest()[:16]


def k_sha(k: list) -> str:
    return hashlib.sha256(json.dumps(sorted(k)).encode()).hexdigest()[:16]


def chance_level(cpdag: MPDAG, base: list, rng: np.random.Generator, retries: int = 20) -> list | None:
    """Same edges as ``base``, orientations uniform among Meek-consistent choices."""
    for _ in range(retries):
        drawn = [(a, b) if rng.integers(2) else (b, a) for (a, b) in base]
        if apply_orientations(cpdag, drawn) is not None:
            return drawn
    return None


def rebuild_knowledge(net: str, arm: str, rep: int, dag: MPDAG, cpdag: MPDAG, know: dict,
                      seed: int = RELEASE_SEED) -> list | None:
    """The arm's knowledge on one network, as ledger_sweep.main builds its plans."""
    base_llm = [tuple(e) for e in (know["D_LLM"]["networks"].get(net) or {}).get("k", [])]
    if arm == "D_RAND":
        net_seed = int(hashlib.sha256(net.encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng([seed, net_seed, rep])
        return chance_level(cpdag, base_llm, rng)
    if arm == "A_TRUE":
        return [(a, b) if dag.is_directed_edge(a, b) else (b, a) for (a, b) in base_llm]
    if arm == "A_ORACLE":
        return sorted(knowledge_to_recover(dag, cpdag))
    if arm == "D_DEGEN":
        return []
    entry = know[arm]["networks"].get(net) or {}
    return [tuple(e) for e in entry.get("k", [])]


# ---------------------------------------------------------------- audit.py


def _d_separated(g: nx.DiGraph, x: set, y: set, z: set) -> bool:
    if hasattr(nx, "is_d_separator"):
        return bool(nx.is_d_separator(g, x, y, z))
    return bool(nx.d_separated(g, x, y, z))


def is_valid_adjustment_set_gac_dag(dag: MPDAG, x: str, y: str, z: Iterable[str]) -> bool:
    """Perković et al. (2018) GAC on a fully oriented graph, over networkx."""
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
