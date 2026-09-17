"""The audit's truth oracle, written apart from the instrument it scores.

The certificate decides validity on the analyst's MPDAG with the generalised
adjustment criterion. When the audit later asks whether the committed set was
valid at the generating DAG, it must ask the same question, so this module
implements that criterion on a DAG over networkx and nothing from the
criterion package. The repository's ``is_valid_adjustment_set_dag`` is Pearl's
back-door criterion, which forbids every descendant of the treatment; it is
sufficient for validity and not necessary, and it rejects sets the adjustment
formula accepts. On the ledger it rejected four of 1,169 committed sets the
criterion accepts, each containing a descendant of the treatment that lies on
no causal path to the outcome.
"""

from __future__ import annotations

from collections.abc import Iterable

import networkx as nx

from bkrobust.demo.graph import MPDAG


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
    return bool(nx.is_d_separator(proper, {x}, {y}, set(zs)))
