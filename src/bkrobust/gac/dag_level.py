r"""The generalised adjustment criterion on a fully oriented DAG.

Why this layer exists
---------------------
:func:`bkrobust.demo.evaluate.is_valid_adjustment_set_dag` implements Pearl's
**back-door** criterion, which is sufficient but *not necessary* for a set to be
a valid adjustment set. The gap is not academic: back-door bans **every**
descendant of the treatment, while what actually has to be banned is only what
sits on, or below, the causal route from the treatment to the outcome. The
smallest witness is ``X -> Y``, ``X -> W`` with ``W`` unattached to ``Y``. Here
``Z = {W}`` gives ``sum_w P(y | x, w) P(w) = P(y | x) = P(y | do(x))`` -- a
perfectly valid adjustment set -- yet back-door rejects it because ``W`` is a
descendant of ``X``.

The **generalised adjustment criterion** (GAC) of Shpitser, VanderWeele and
Robins (2010) and Perkovic, Textor, Kalisch and Maathuis (2018) closes that gap:
it is sound *and complete* for adjustment. This module states it at the DAG
level, where every definition is unambiguous, and it is the reference against
which the MPDAG lift in :mod:`bkrobust.gac.mpdag_level` is validated.

The criterion
-------------
``Z`` is a valid adjustment set for ``(x, y)`` in the DAG ``D`` iff

1. ``Z & forb(x, y, D) = {}``, where ``forb`` is everything on or downstream of
   a proper causal path from ``x`` to ``y``, together with ``x`` and ``y``; and
2. ``Z`` blocks every proper **non-causal** path from ``x`` to ``y``.

A *proper* path from ``x`` to ``y`` meets ``{x}`` only at its first node, which
is automatic for a single-node treatment; the word is kept because the
literature's definitions are stated for node **sets** ``X`` and the distinction
is real there.

The reduction used for condition 2, derived
-------------------------------------------
Condition 2 is a quantifier over paths. It is decided here by a single
d-separation test on a modified graph, exactly as the back-door criterion is:

    Let ``D'`` be ``D`` with every edge ``x -> v`` deleted for which ``v`` lies
    on a proper causal path from ``x`` to ``y`` -- i.e. the *first edges of the
    proper causal paths*, and only those. **Given condition 1**, ``Z`` blocks
    every proper non-causal path from ``x`` to ``y`` in ``D`` iff ``x`` and
    ``y`` are d-separated by ``Z`` in ``D'``.

Note the hypothesis: the equivalence is not unconditional, it holds once
condition 1 has been checked. That is why this module tests condition 1 first
and why the order is not cosmetic. The three steps:

* *Every proper causal path is destroyed.* Such a path starts ``x -> v`` with
  ``v`` a child of ``x`` on it, hence ``v in cn(x, y, D)``, hence that edge is
  deleted. So no path of ``D'`` from ``x`` to ``y`` is causal, and the paths
  that remain are exactly the proper non-causal paths of ``D`` that do not begin
  with a deleted edge.
* *The non-causal paths that are destroyed were blocked anyway.* Let ``p`` be a
  proper non-causal path starting ``x -> v`` with ``v in cn(x, y, D)``. Walk
  along ``p`` from ``x``; since ``p`` is non-causal some edge eventually points
  backwards, so there is a first node ``c`` entered forwards and left by an edge
  pointing back into it. The segment ``x -> v -> ... -> c`` is directed, so
  ``c`` is in ``de(v) u {v}``, hence in ``forb``, and so is every descendant of
  ``c``. ``c`` is a collider on ``p`` and condition 1 puts none of its
  descendants (nor ``c`` itself) in ``Z``, so ``c`` blocks ``p``. Deleting the
  edge therefore cannot hide an *open* path.
* *Blocking does not change between ``D`` and ``D'``.* The only way it could is
  a collider ``c`` opened in ``D`` by a directed path ``c ~> z_i`` that runs
  through a deleted edge ``x -> v``. Then ``z_i in de(v) subset forb``, which
  condition 1 forbids. So every collider open in ``D`` is open in ``D'``, and
  conversely ``D'`` is a subgraph of ``D``.

The path-enumerating :func:`is_gac_valid_dag_by_paths` is kept as an independent
second implementation of the same predicate: it checks condition 2 directly,
path by path, with no reduction at all. The two are asserted equal over every
DAG on 4 nodes in ``tests/gac/test_gac.py``, which is what turns the derivation
above from an argument into a measured claim.

Conventions fixed here
----------------------
* ``forb`` always contains ``x`` **and** ``y``, whether or not a causal path
  from ``x`` to ``y`` exists. The literature's ``forb`` picks ``y`` up via
  ``cn`` whenever a causal path exists; adding it unconditionally only affects
  the causally disconnected case, where adjusting for the outcome itself is
  meaningless anyway, and it keeps ``Z ∩ {x, y} = {}`` a stated part of the
  criterion rather than an accident of d-separation's behaviour on a
  conditioned endpoint.
* ``x == y`` is a degenerate query and is answered False, matching
  :func:`bkrobust.core.oracle.is_valid`.

Written to run on Python 3.9 as well as the repository's 3.11 target: builtin
generics and ``X | Y`` unions appear only in annotations, which
``from __future__ import annotations`` defers to strings, and neither
``zip(strict=)`` nor ``itertools.pairwise`` (3.10+) is used.
"""

from __future__ import annotations

from collections.abc import Iterable

from bkrobust.demo.evaluate import is_dseparated
from bkrobust.demo.graph import MPDAG
from bkrobust.mpdag_criterion.paths import is_blocked, is_non_causal, simple_paths

Node = str


def causal_nodes_dag(dag: MPDAG, x: Node, y: Node) -> set[Node]:
    """``cn(x, y, D)``: the nodes on proper causal paths from ``x`` to ``y``, minus ``x``.

    Computed as ``de(x) & (an(y) u {y})``: a node ``v != x`` lies on some
    directed path ``x -> ... -> v -> ... -> y`` iff ``x`` reaches it and it
    reaches ``y``. In a DAG the two segments concatenate into a *simple* path,
    because a node cannot be both an ancestor and a descendant of another
    without closing a cycle -- this is the step that fails on an MPDAG and is
    why :mod:`bkrobust.gac.mpdag_level` cannot reuse this formula.

    ``y`` itself is in the result whenever any causal path exists; ``x`` never
    is, since :meth:`~bkrobust.demo.graph.MPDAG.descendants` is exclusive.

    Args:
        dag: A fully oriented DAG.
        x: The treatment node.
        y: The outcome node.

    Returns:
        The causal nodes, empty when no directed path from ``x`` to ``y`` exists.
    """
    return dag.descendants(x) & (dag.ancestors(y) | {y})


def forbidden_set_dag(dag: MPDAG, x: Node, y: Node) -> set[Node]:
    """``forb(x, y, D)``: nodes on or downstream of a proper causal path, plus ``x`` and ``y``.

    ``cn(x, y, D) u de(cn(x, y, D)) u {x, y}``. This is the set the GAC bars
    from ``Z``, and it is in general a **strict subset** of the back-door bar
    ``de(x) u {x, y}``: a descendant of ``x`` that neither lies on nor hangs off
    the causal route to ``y`` is forbidden by back-door and permitted here. That
    difference is the entire behavioural gap between the two criteria; see the
    module docstring for the minimal witness.

    On the inclusion of ``y``: see the module docstring's conventions. It only
    matters when no causal path from ``x`` to ``y`` exists, in which case the
    literature's ``forb`` would be ``{x}`` alone.

    Args:
        dag: A fully oriented DAG.
        x: The treatment node.
        y: The outcome node.

    Returns:
        The forbidden set, always containing ``x`` and ``y``.
    """
    cn = causal_nodes_dag(dag, x, y)
    out: set[Node] = set(cn) | {x, y}
    for node in sorted(cn):
        out |= dag.descendants(node)
    return out


def _causal_first_edges(dag: MPDAG, x: Node, cn: set[Node]) -> set[tuple[Node, Node]]:
    """The edges out of ``x`` that begin a proper causal path to the outcome.

    An edge ``x -> v`` begins one iff ``v`` is itself a causal node, i.e. lies on
    a directed path from ``x`` to the outcome.

    Args:
        dag: A fully oriented DAG.
        x: The treatment node.
        cn: ``causal_nodes_dag(dag, x, y)``, passed in so it is computed once.

    Returns:
        The qualifying ``(tail, head)`` pairs.
    """
    return {(x, v) for v in dag.children(x) & cn}


def is_gac_valid_dag(dag: MPDAG, x: Node, y: Node, z: Iterable[Node]) -> bool:
    """Whether ``z`` satisfies the generalised adjustment criterion for ``(x, y)`` in ``dag``.

    The two conditions of the module docstring: ``z`` misses
    :func:`forbidden_set_dag`, and ``z`` blocks every proper non-causal path
    from ``x`` to ``y`` -- the latter decided by one d-separation test on the
    graph with the first edges of the proper causal paths deleted, which is
    equivalent *given* the first condition (derivation in the module docstring).

    This predicate is strictly weaker than
    :func:`bkrobust.demo.evaluate.is_valid_adjustment_set_dag`: every back-door
    valid set is GAC-valid, and some GAC-valid sets are not back-door valid.

    Args:
        dag: A fully oriented DAG.
        x: The treatment node.
        y: The outcome node.
        z: The candidate adjustment set.

    Returns:
        True iff ``z`` is a valid adjustment set for ``(x, y)`` in ``dag``.
    """
    if x == y:
        return False
    z_set = frozenset(z)
    if z_set & forbidden_set_dag(dag, x, y):
        return False
    cn = causal_nodes_dag(dag, x, y)
    drop = _causal_first_edges(dag, x, cn)
    trimmed = MPDAG(
        dag.nodes,
        [e for e in dag.directed_edges if e not in drop],
        dag.undirected_edges,
    )
    return is_dseparated(trimmed, x, y, z_set)


def is_gac_valid_dag_by_paths(dag: MPDAG, x: Node, y: Node, z: Iterable[Node]) -> bool:
    """The same predicate, deciding condition 2 by enumerating paths instead.

    .. warning::
       **Exponential in the number of vertices.** A dense graph on ``n`` nodes
       has ``Theta((n-2)!)`` simple paths between two of them. This is a
       reference implementation for tests at small ``n``, never a decision
       procedure -- call :func:`is_gac_valid_dag` instead.

    It exists because the reduction behind :func:`is_gac_valid_dag` is an
    argument (module docstring) rather than a proof-checked theorem, and an
    argument is worth what its differential test is worth. This function makes
    no use of the reduction: it walks every proper non-causal path and asks
    :func:`bkrobust.mpdag_criterion.paths.is_blocked` directly. Every node of a
    path in a DAG is of definite status, so the definite-status qualification
    that the MPDAG-level criterion needs is vacuous here.

    Args:
        dag: A fully oriented DAG.
        x: The treatment node.
        y: The outcome node.
        z: The candidate adjustment set.

    Returns:
        True iff ``z`` is a valid adjustment set for ``(x, y)`` in ``dag``.
    """
    if x == y:
        return False
    z_set = frozenset(z)
    if z_set & forbidden_set_dag(dag, x, y):
        return False
    for path in simple_paths(dag, x, y):
        if is_non_causal(dag, path) and not is_blocked(dag, path, z_set):
            return False
    return True


__all__ = [
    "causal_nodes_dag",
    "forbidden_set_dag",
    "is_gac_valid_dag",
    "is_gac_valid_dag_by_paths",
]
