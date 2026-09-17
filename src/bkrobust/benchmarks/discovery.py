"""Estimated CPDAGs for the deployment panel: sampling, discovery, conversion, scoring.

The ledger's CPDAG is the true DAG's. The estimated panel replaces it with the
output of a discovery algorithm run on samples the true DAG generated, which is
the practitioner's situation. This module holds what that panel needs and the
ledger does not: the Gaussian parameters the four bnlearn networks ship with,
ancestral sampling from them, the two discovery algorithms through the
``causal-learn`` package, the conversion of its graph encoding into
:class:`~bkrobust.demo.graph.MPDAG`, and the structure statistics a learned
graph is scored on against the truth.

``causal-learn`` is imported lazily and only by :func:`learn_graph`; nothing
else here needs it, and it is not a declared dependency of the package. Its
endpoint encoding, on the adjacency matrix ``M`` of a ``GeneralGraph``:
``M[i, j] = -1`` with ``M[j, i] = 1`` is ``i -> j``; both ``-1`` is ``i -- j``;
both ``1`` is a conflict, an arrowhead at each end, which PC produces when two
v-structures disagree about one edge. A conflict carries no direction an
analyst could act on, so it is converted to an undirected edge and counted.
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path

import numpy as np
from sklearn.metrics import adjusted_rand_score

from bkrobust.demo.graph import MPDAG, canon, undirected_components
from bkrobust.demo.meek import is_meek_closed, meek_closure

Node = str
Edge = tuple[str, str]

#: The two discovery methods the panel runs.
METHODS: tuple[str, ...] = ("pc", "ges")


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


def to_endpoint_matrix(g: MPDAG) -> np.ndarray:
    """The inverse of :func:`from_endpoint_matrix` on a graph without conflicts."""
    idx = {v: i for i, v in enumerate(g.nodes)}
    m = np.zeros((len(idx), len(idx)), dtype=int)
    for a, b in g.directed_edges:
        m[idx[a], idx[b]] = -1
        m[idx[b], idx[a]] = 1
    for a, b in g.undirected_edges:
        m[idx[a], idx[b]] = -1
        m[idx[b], idx[a]] = -1
    return m


def partition_labels(g: MPDAG) -> list[int]:
    """A chain-component label per node of ``g``, in ``g.nodes`` order.

    Nodes in the same chain component share a label; every node outside any
    component of size two or more gets a label of its own, so two graphs whose
    components coincide get an adjusted Rand index of one.
    """
    label = {v: -1 for v in g.nodes}
    for k, comp in enumerate(undirected_components(g)):
        for v in comp:
            label[v] = k
    nxt = len(undirected_components(g))
    out = []
    for v in g.nodes:
        if label[v] < 0:
            label[v] = nxt
            nxt += 1
        out.append(label[v])
    return out


def structure_scores(g: MPDAG, dag: MPDAG, oracle: MPDAG) -> tuple[dict, MPDAG]:
    """Score a learned graph against the truth and the oracle CPDAG.

    The graph every later stage uses is the Meek closure of ``g`` where the
    closure succeeds and ``g`` itself where it does not; both facts are
    recorded. Skeleton precision and recall are against the true DAG's
    skeleton; orientation accuracy is over the learned directed edges that
    exist in the truth; the adjusted Rand index compares chain-component
    partitions with the oracle's.

    Args:
        g: The converted learned graph.
        dag: The generating DAG.
        oracle: ``dag``'s CPDAG.

    Returns:
        The scores and the graph to use downstream.
    """
    closed_on_arrival = is_meek_closed(g)
    closure = meek_closure(g)
    used = closure if closure is not None else g
    true_skel = dag.skeleton()
    skel = used.skeleton()
    tp = len(true_skel & skel)
    n_dir, n_und = len(used.directed_edges), len(used.undirected_edges)
    present = [e for e in used.directed_edges if canon(*e) in true_skel]
    comps = undirected_components(used)
    scores = {
        "meek_closed_on_arrival": int(closed_on_arrival),
        "closure_ok": int(closure is not None),
        "n_forced_by_closure": len(g.undirected_edges) - n_und,
        "acyclic": int(used.is_acyclic()),
        "n_learned_edges": n_dir + n_und,
        "n_directed": n_dir,
        "n_undirected": n_und,
        "undirected_fraction": round(n_und / (n_dir + n_und), 6) if n_dir + n_und else None,
        "skeleton_tp": tp,
        "skeleton_precision": round(tp / len(skel), 6) if skel else None,
        "skeleton_recall": round(tp / len(true_skel), 6) if true_skel else None,
        "n_directed_in_truth": len(present),
        "oriented_correct": (
            round(sum(1 for a, b in present if dag.is_directed_edge(a, b)) / len(present), 6)
            if present
            else None
        ),
        "n_components_ge2": len(comps),
        "max_component": max((len(c) for c in comps), default=0),
        "components": [sorted(c) for c in comps],
        "partition": partition_labels(used),
        "ari_vs_oracle": round(
            float(adjusted_rand_score(partition_labels(oracle), partition_labels(used))), 6
        ),
    }
    return scores, used
