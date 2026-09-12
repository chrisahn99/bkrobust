"""The Henckel--Perkovic--Maathuis optimal adjustment set, at MPDAG level.

Why this module exists
----------------------

:func:`bkrobust.demo.evaluate.optimal_adjustment_set_mpdag` returns the optimal
set only when *every* DAG extension of ``G`` agrees on it, and ``None``
otherwise. That is a defensible reading of "identified by ``G`` alone", but it
is **not** the object the adjustment literature calls the optimal set of an
MPDAG, and the difference is not cosmetic. Under extension-agreement, any ``G``
that admits the true DAG and returns a set at all must return the true DAG's own
optimal set, so achievable asymptotic variance cannot vary with what the analyst
asserts: the conclusion restates the definition.

The graphical object is

.. math::

    \\mathbf{O}(X, Y, G) \;=\; \\mathrm{pa}\\bigl(\\mathrm{cn}(X, Y, G)\\bigr)
                              \;\\setminus\; \\mathrm{forb}(X, Y, G),

defined on ``G`` directly and returning a set whenever ``G`` is amenable. It can
move as knowledge grows: orienting edges adds parents to the causal nodes and
can drop nodes from the forbidden set, so more knowledge can buy a different --
and possibly better -- set. That is the behaviour the extension-agreement
reading hides.

Both readings are useful and this module replaces neither. It exists so a claim
about efficiency can be checked under the standard definition rather than under
a project-local one.

Cost and scope
--------------

``cn`` and ``forb`` come from :mod:`bkrobust.mpdag_criterion.criterion`, which
states them as the literature does and is **exponential** in the number of
vertices; that module's own caveats about shielded possibly-causal paths apply
here unchanged. Instances in the efficiency sweep have components of 6--12
nodes, where the cost is irrelevant. Do not call this at scale.

Optimality of ``O`` among GAC-valid sets is a result of the cited literature and
is not re-derived here. :func:`check_agrees_with_dag_level` is the acceptance
test actually run: on a fully oriented graph this must reproduce
:func:`bkrobust.demo.evaluate.optimal_adjustment_set_dag` exactly.
"""

from __future__ import annotations

from bkrobust.demo.graph import MPDAG, Node
from bkrobust.demo.evaluate import optimal_adjustment_set_dag
from bkrobust.mpdag_criterion.criterion import causal_nodes, forbidden_set, is_amenable


def optimal_adjustment_set_mpdag_hpm(g: MPDAG, x: Node, y: Node) -> set[Node] | None:
    """``O(x, y, g) = pa(cn(x, y, g)) \\ forb(x, y, g)``, or ``None`` if not amenable.

    A function of ``g`` alone: no extension enumeration, and a set is returned
    whenever the effect is identified by adjustment at all.

    Args:
        g: The MPDAG, assumed Meek-closed.
        x: The treatment node.
        y: The outcome node.

    Returns:
        The optimal adjustment set, or ``None`` when ``g`` is not amenable
        relative to ``(x, y)``, or when no possibly causal path from ``x`` to
        ``y`` exists.
    """
    if not is_amenable(g, x, y):
        return None
    cn = causal_nodes(g, x, y)
    if not cn:
        return None
    pa_cn: set[Node] = set()
    for n in sorted(cn):
        pa_cn |= g.parents(n)
    return pa_cn - forbidden_set(g, x, y)


def check_agrees_with_dag_level(dag: MPDAG, x: Node, y: Node) -> bool:
    """Acceptance test: on a fully oriented graph the MPDAG formula must agree.

    Args:
        dag: A fully oriented DAG.
        x: The treatment node.
        y: The outcome node.

    Returns:
        ``True`` if :func:`optimal_adjustment_set_mpdag_hpm` reproduces
        :func:`~bkrobust.demo.evaluate.optimal_adjustment_set_dag` on ``dag``.

    Raises:
        ValueError: If ``dag`` still carries undirected edges.
    """
    if dag.undirected_edges:
        raise ValueError("check_agrees_with_dag_level needs a fully oriented graph")
    got = optimal_adjustment_set_mpdag_hpm(dag, x, y)
    want = set(optimal_adjustment_set_dag(dag, x, y))
    if got is None:
        return not want and y not in dag.descendants(x)
    return got == want
