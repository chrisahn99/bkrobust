"""Vendored from src/bkrobust/benchmarks/audit.py.

commit ba7ddc12c446f1087e58dc428720c84d33b1b622
git status of that file at copy time: 'A  src/bkrobust/benchmarks/audit.py' (staged, not committed; copied from the
working tree as it stood on 2026-11-06). Code below this docstring is copied unchanged from
the file's 'from __future__' line on. The submitted tree lacks this module; the released column z_valid_at_truth was computed with it.
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
