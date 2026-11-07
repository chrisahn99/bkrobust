"""Functions E7 needs that the submitted tree lacks, copied from the second author's tree.

Source repository: bkrobust (main branch)
    commit ba7ddc12c446f1087e58dc428720c84d33b1b622
The source files below are STAGED in that repository's index but not in HEAD, so the commit
alone does not identify them; their sha256 at copy time (2026-11-06) is given instead.

    experiments/questionnaire_build.py        67e572f5654a542b3e43ed873994ca410ff94e804c8a151e2b1826fabd1873f0
        stable_seed
    experiments/ledger_sweep.py               f84eccd7fd5930bf9ab6e05158c110ca3675cb3eb62fb870a6581f8c2324ea3c
        sha
    src/bkrobust/benchmarks/discovery.py      2702a484e83ca9e9846a87a199d14119f0dda45aa26a65fda6dcb9f8569dabf5
        gaussian_parameters, topological_order, sample_ancestral, learn_graph, from_endpoint_matrix
    experiments/pest_panel.py                 2a18c4e51b09942f87df30422ad3f33ee931912e5b5c22e7da57b77315e2a624
        ALPHA, PANEL_SEED, load_gaussian, sample_for, data_sha, graph_from_json
    src/bkrobust/mpdag_criterion/optimal.py   538cbb6f63b4eea978f433b8b7c5b4fb63c2bf683f69ca2d6b781224cd9ed734
        CommittedSet, possible_ancestors, causal_nodes_poly, committed_adjustment_set
    src/bkrobust/benchmarks/audit.py          fccfd536d28cc48ad8afe5cc73a393c5904b1690c458188e298d0d3cb2782146
        is_valid_adjustment_set_gac_dag

The function bodies were extracted with the ast module, so they are verbatim except for two
changes marked "E7 change": load_gaussian takes the models folder as an argument, and the
networkx d-separation call falls back to the pre-3.3 name. Everything they import (MPDAG,
Meek closure, GAC, amenability, possible descendants, the network parser) comes from the
submitted repository's src/; those modules are byte-identical between the two trees
(checked 2026-11-06: demo/graph.py, demo/meek.py, gac/*.py, mpdag_criterion/criterion.py,
mpdag_criterion/paths.py, benchmarks/describe.py). causal-learn is imported lazily, only by
learn_graph.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import networkx as nx
import numpy as np

from bkrobust.benchmarks.describe import parse_file, to_mpdag
from bkrobust.demo.graph import MPDAG, canon
from bkrobust.gac.mpdag_level import gac_forbidden_set, is_gac_valid_mpdag
from bkrobust.mpdag_criterion.criterion import is_amenable
from bkrobust.mpdag_criterion.paths import possible_descendants

Node = str
Edge = tuple[str, str]


# --- experiments/questionnaire_build.py ------------------------------------
def stable_seed(name: str) -> int:
    """A per-network seed that does not depend on the interpreter's hash salt."""
    return int(hashlib.sha256(name.encode()).hexdigest()[:8], 16)



# --- experiments/ledger_sweep.py -------------------------------------------
def sha(g: MPDAG) -> str:
    """A hash of a graph's edge string, for the provenance stamp."""
    return hashlib.sha256(g.edge_string().encode()).hexdigest()[:16]



# --- src/bkrobust/benchmarks/discovery.py ----------------------------------
def gaussian_parameters(path: Path) -> tuple[dict[Edge, float], dict[Node, float]]:
    """Fitted coefficients and residual variances of a bnlearn Gaussian JSON network.

    The four Gaussian networks of the corpus carry, beside ``nodes`` and
    ``arcs``, a ``cpds`` block with one linear equation per node: an intercept,
    a coefficient per parent and a residual variance, all fitted to real data.
    The intercept is dropped, since a CPDAG and an adjustment estimand are
    invariant to it.

    Args:
        path: The ``.json`` file.

    Returns:
        ``(weights, noise_var)`` with ``weights[(parent, child)]`` the
        coefficient of ``parent`` in ``child``'s equation.
    """
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
    """``n`` draws from the linear-Gaussian SEM, columns in ``dag.nodes`` order.

    Each node is its parents' weighted sum plus Gaussian noise, generated in a
    topological order. Parents are visited sorted, so the draw sequence does not
    depend on set iteration order.

    Args:
        dag: The generating DAG.
        weights: Structural coefficients keyed ``(parent, child)``.
        noise_var: Residual variance per node.
        n: Number of rows.
        rng: The generator to draw from.

    Returns:
        An ``(n, p)`` array with column ``i`` the node ``dag.nodes[i]``.
    """
    idx = {v: i for i, v in enumerate(dag.nodes)}
    data = np.zeros((n, len(idx)))
    for v in topological_order(dag):
        col = rng.normal(0.0, float(np.sqrt(noise_var[v])), size=n)
        for p in sorted(dag.parents(v)):
            col = col + weights[(p, v)] * data[:, idx[p]]
        data[:, idx[v]] = col
    return data

def learn_graph(data: np.ndarray, method: str, alpha: float = 0.01) -> np.ndarray:
    """Run one discovery algorithm and return causal-learn's endpoint matrix.

    Args:
        data: An ``(n, p)`` sample.
        method: ``"pc"`` for PC with the Fisher-z test at level ``alpha``, the
            stable skeleton and the package's default collider handling;
            ``"ges"`` for GES with the BIC score.
        alpha: The test level for PC; ignored by GES.

    Returns:
        The ``(p, p)`` integer endpoint matrix.

    Raises:
        ValueError: On an unknown method.
    """
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
    """Convert causal-learn's endpoint matrix into an :class:`MPDAG`.

    Args:
        matrix: The ``(p, p)`` endpoint matrix, in the order of ``nodes``.
        nodes: Node names, one per row of ``matrix``.

    Returns:
        The graph and a count of what the matrix held: ``directed``,
        ``undirected``, ``conflict`` (two arrowheads, carried as undirected) and
        ``other`` (any endpoint pair outside the CPDAG vocabulary, carried as
        undirected as well).
    """
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



# --- experiments/pest_panel.py ---------------------------------------------
ALPHA = 0.01

PANEL_SEED = 20260912

def load_gaussian(models: Path, name: str) -> tuple[MPDAG, dict, dict, dict[str, str]]:
    """The true DAG, its fitted coefficients and residual variances, and the raw names."""
    # E7 change: the models folder is an argument instead of the module constant MODELS
    path = models / f"{name}.json"
    parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
    dag = to_mpdag(parsed)
    weights, noise = gaussian_parameters(path)
    if set(weights) != set(dag.directed_edges):
        raise ValueError(f"{name}: coefficients and arcs disagree")
    return dag, weights, noise, dict(parsed.name_map)

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



# --- src/bkrobust/mpdag_criterion/optimal.py -------------------------------
@dataclass(frozen=True)
class CommittedSet:
    """The set an analyst commits to on their own graph, and how it was reached.

    Attributes:
        z: The adjustment set, or ``None`` when no valid set exists.
        verdict: One of ``closed_form`` (the optimal-set formula, valid on
            ``G0``), ``canonical`` (the formula's set was not valid and the
            canonical set was committed to instead), ``y_not_possible_descendant``
            (the analyst's own graph says the treatment cannot cause the outcome,
            so there is no effect to adjust for; the frame drew the pair from a
            CPDAG in which it was possible and the knowledge oriented it away),
            ``not_amenable`` (no adjustment set can identify the effect on
            ``G0``), or ``no_valid_set`` (amenable, the outcome possibly caused,
            and yet neither candidate is valid, which on an amenable maximal PDAG
            should not happen and is reported rather than swallowed).
        causal_nodes: The polynomial causal-node approximation used.
        forbidden: The GAC forbidden set used.
    """

    z: frozenset[Node] | None
    verdict: str
    causal_nodes: frozenset[Node]
    forbidden: frozenset[Node]

def possible_ancestors(g: MPDAG, node: Node) -> set[Node]:
    """Nodes that can reach ``node`` along a possibly causal path, plus ``node``."""
    return {u for u in g.nodes if node in possible_descendants(g, u)} | {node}

def causal_nodes_poly(g: MPDAG, x: Node, y: Node) -> set[Node]:
    """Nodes on a possibly causal path from ``x`` to ``y``, by reachability.

    ``possde(x)`` intersected with ``possan(y)`` plus ``y``, minus ``x``. On a DAG this is
    exact. On a maximal PDAG it is the reachability shortcut that
    :func:`bkrobust.mpdag_criterion.criterion.causal_nodes` documents as
    inexact, because two possibly causal sub-paths need not concatenate into a
    simple one. It is used here because the alternative enumerates paths and is
    exponential, and because the set it feeds is checked for validity on ``G0``
    before anything is committed to.
    """
    down = possible_descendants(g, x)
    return {v for v in down if v != x and (v == y or y in possible_descendants(g, v))}

def committed_adjustment_set(g: MPDAG, x: Node, y: Node) -> CommittedSet:
    """The set the analyst commits to on ``g``, polynomial, with its provenance.

    Args:
        g: The analyst's graph, Meek-closed.
        x: Treatment.
        y: Outcome.

    Returns:
        A :class:`CommittedSet`. ``z`` is ``None`` only when the verdict is
        ``not_amenable`` or ``no_valid_set``.
    """
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



# --- src/bkrobust/benchmarks/audit.py --------------------------------------
def is_valid_adjustment_set_gac_dag(dag: MPDAG, x: str, y: str, z: Iterable[str]) -> bool:
    """Whether ``z`` satisfies the generalised adjustment criterion in ``dag``.

    Perković, Textor, Kalisch and Maathuis (2018) on a fully oriented graph:
    ``z`` contains neither endpoint, avoids the forbidden set, which is the
    treatment together with every descendant of a node on a causal path from
    treatment to outcome, and d-separates treatment from outcome in the proper
    back-door graph, where the first edge of every such causal path is removed.

    Args:
        dag: A fully oriented graph.
        x: The treatment.
        y: The outcome.
        z: The candidate adjustment set.

    Returns:
        True iff the adjustment formula over ``z`` identifies the effect.
    """
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
    # E7 change: is_d_separator appeared in networkx 3.3; d_separated is the older name
    if hasattr(nx, "is_d_separator"):
        return bool(nx.is_d_separator(proper, {x}, {y}, set(zs)))
    return bool(nx.d_separated(proper, {x}, {y}, set(zs)))
