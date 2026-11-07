"""Functions copied from the second author's working tree, unchanged in behaviour.

Source repository: bkrobust (main branch)
commit ba7ddc12c446f1087e58dc428720c84d33b1b622
The source files below are staged in that repository's index but NOT in HEAD, so
their content is pinned by sha256 instead:

  src/bkrobust/benchmarks/discovery.py      2702a484e83ca9e9846a87a199d14119f0dda45aa26a65fda6dcb9f8569dabf5
      gaussian_parameters, topological_order, sample_ancestral, learn_graph,
      from_endpoint_matrix
  experiments/pest_panel.py                 2a18c4e51b09942f87df30422ad3f33ee931912e5b5c22e7da57b77315e2a624
      PANEL_SEED, ALPHA, WALL_S, load_gaussian, sample_for, data_sha, graph_from_json
  src/bkrobust/mpdag_criterion/optimal.py   538cbb6f63b4eea978f433b8b7c5b4fb63c2bf683f69ca2d6b781224cd9ed734
      CommittedSet, possible_ancestors, causal_nodes_poly, committed_adjustment_set
  experiments/ledger_sweep.py               f84eccd7fd5930bf9ab6e05158c110ca3675cb3eb62fb870a6581f8c2324ea3c
      sha, chance_level
  experiments/questionnaire_build.py        67e572f5654a542b3e43ed873994ca410ff94e804c8a151e2b1826fabd1873f0
      stable_seed

Everything they import (MPDAG, Meek closure, GAC, amenability, possible
descendants, the network parser) is taken from the submitted repository's
``src/``; those modules were checked byte-identical between the two trees on
2026-11-06. ``causal-learn`` is imported lazily, only by ``learn_graph``.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from bkrobust.benchmarks.describe import parse_file, to_mpdag
from bkrobust.demo.graph import MPDAG, canon
from bkrobust.demo.meek import apply_orientations
from bkrobust.gac.mpdag_level import gac_forbidden_set, is_gac_valid_mpdag
from bkrobust.mpdag_criterion.criterion import is_amenable
from bkrobust.mpdag_criterion.paths import possible_descendants

Node = str
Edge = tuple[str, str]

PANEL_SEED = 20260912
ALPHA = 0.01
WALL_S = 900


# --- questionnaire_build.py / ledger_sweep.py ----------------------------------


def stable_seed(name: str) -> int:
    """A per-network seed that does not depend on the interpreter's hash salt."""
    return int(hashlib.sha256(name.encode()).hexdigest()[:8], 16)


def sha(g: MPDAG) -> str:
    """A hash of a graph's edge string, for the provenance stamp."""
    return hashlib.sha256(g.edge_string().encode()).hexdigest()[:16]


def chance_level(
    cpdag: MPDAG, base: list, rng: np.random.Generator, retries: int = 20
) -> list | None:
    """Same edges as ``base``, orientations uniform among Meek-consistent choices."""
    for _ in range(retries):
        drawn = [(a, b) if rng.integers(2) else (b, a) for (a, b) in base]
        if apply_orientations(cpdag, drawn) is not None:
            return drawn
    return None


# --- discovery.py ----------------------------------------------------------------


def gaussian_parameters(path: Path) -> tuple[dict[Edge, float], dict[Node, float]]:
    """Fitted coefficients and residual variances of a bnlearn Gaussian JSON network."""
    data = json.loads(path.read_text())
    weights: dict[Edge, float] = {}
    noise_var: dict[Node, float] = {}
    for child, cpd in data["cpds"].items():
        noise_var[child] = float(cpd["variance"][0])
        for parent, coef in cpd["coefficients"].items():
            if parent != "(Intercept)":
                weights[(parent, child)] = float(coef[0])
    return weights, noise_var


def topological_order(dag: MPDAG) -> list[Node]:
    """A parents-first order of ``dag``'s nodes, deterministic for a given graph."""
    remaining = {v: len(dag.parents(v)) for v in dag.nodes}
    ready = sorted(v for v, d in remaining.items() if d == 0)
    order: list[Node] = []
    while ready:
        v = ready.pop(0)
        order.append(v)
        for c in sorted(dag.children(v)):
            remaining[c] -= 1
            if remaining[c] == 0:
                ready.append(c)
        ready.sort()
    if len(order) != len(dag.nodes):
        raise ValueError("graph is not acyclic")
    return order


def sample_ancestral(
    dag: MPDAG,
    weights: dict[Edge, float],
    noise_var: dict[Node, float],
    n: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """``n`` draws from the linear-Gaussian SEM, columns in ``dag.nodes`` order."""
    idx = {v: i for i, v in enumerate(dag.nodes)}
    data = np.zeros((n, len(idx)))
    for v in topological_order(dag):
        col = rng.normal(0.0, float(np.sqrt(noise_var[v])), size=n)
        for p in sorted(dag.parents(v)):
            col = col + weights[(p, v)] * data[:, idx[p]]
        data[:, idx[v]] = col
    return data


def learn_graph(data: np.ndarray, method: str, alpha: float = 0.01) -> np.ndarray:
    """Run one discovery algorithm and return causal-learn's endpoint matrix."""
    if method == "pc":
        from causallearn.search.ConstraintBased.PC import pc

        return np.asarray(
            pc(data, alpha=alpha, indep_test="fisherz", stable=True, show_progress=False).G.graph
        )
    if method == "ges":
        from causallearn.search.ScoreBased.GES import ges

        return np.asarray(ges(data, score_func="local_score_BIC")["G"].graph)
    raise ValueError(f"unknown discovery method {method!r}")


def from_endpoint_matrix(matrix: np.ndarray, nodes: Sequence[Node]) -> tuple[MPDAG, dict[str, int]]:
    """Convert causal-learn's endpoint matrix into an :class:`MPDAG`."""
    directed: list[Edge] = []
    undirected: list[Edge] = []
    counts = {"directed": 0, "undirected": 0, "conflict": 0, "other": 0}
    p = len(nodes)
    for i in range(p):
        for j in range(i + 1, p):
            a, b = int(matrix[i, j]), int(matrix[j, i])
            if a == 0 and b == 0:
                continue
            if a == -1 and b == 1:
                directed.append((nodes[i], nodes[j]))
                counts["directed"] += 1
            elif a == 1 and b == -1:
                directed.append((nodes[j], nodes[i]))
                counts["directed"] += 1
            elif a == -1 and b == -1:
                undirected.append(canon(nodes[i], nodes[j]))
                counts["undirected"] += 1
            elif a == 1 and b == 1:
                undirected.append(canon(nodes[i], nodes[j]))
                counts["conflict"] += 1
            else:
                undirected.append(canon(nodes[i], nodes[j]))
                counts["other"] += 1
    return MPDAG(nodes, directed=directed, undirected=undirected), counts


# --- pest_panel.py ---------------------------------------------------------------


def load_gaussian(models: Path, name: str) -> tuple[MPDAG, dict, dict]:
    """The true DAG, its fitted coefficients and residual variances.

    Same as the panel's ``load_gaussian`` with the models folder as an argument
    instead of a module constant.
    """
    path = models / f"{name}.json"
    parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
    dag = to_mpdag(parsed)
    weights, noise = gaussian_parameters(path)
    if set(weights) != set(dag.directed_edges):
        raise ValueError(f"{name}: coefficients and arcs disagree")
    return dag, weights, noise


def sample_for(dag: MPDAG, weights: dict, noise: dict, name: str, n: int, seed: int) -> np.ndarray:
    """The panel's sample for one cell, reproducible from its four seeds."""
    rng = np.random.default_rng([PANEL_SEED, stable_seed(name), n, seed])
    return sample_ancestral(dag, weights, noise, n, rng)


def data_sha(data: np.ndarray) -> str:
    """A short hash of a sample, so a later stage can prove it regenerated the same one."""
    return hashlib.sha256(np.ascontiguousarray(data).tobytes()).hexdigest()[:16]


def graph_from_json(e: dict) -> MPDAG:
    """An MPDAG from the ``nodes``/``directed``/``undirected`` fields of a record."""
    return MPDAG(
        e["nodes"],
        directed=[tuple(x) for x in e["directed"]],
        undirected=[tuple(x) for x in e["undirected"]],
    )


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
