"""The adjustment set an analyst commits to, in polynomial time.

``demo.evaluate.optimal_adjustment_set_mpdag`` computes the DAG-level optimal set
in every extension of the analyst's graph and returns ``None`` the moment two of
them disagree. That has two consequences neither of which is about the query.
It is exponential in the number of undirected edges the analyst's graph retains,
so a knowledge set that leaves a component partly open cannot be evaluated at
all above a small cap. And it deletes exactly the instances where the analyst's
knowledge is incomplete, which is every real deployment.

This module gives the analyst a set to commit to by the closed form of Henckel,
Perkovic and Maathuis lifted to a maximal PDAG: the parents of the causal nodes,
minus the forbidden set. Every piece is polynomial and already in the tree.
Where that set is not valid on the analyst's graph the canonical adjustment set
of Perkovic et al. is tried, which is valid whenever any valid set exists on an
amenable graph. The caller learns which of the two was committed to, because
the second is a valid set and not an optimal one, and any claim about
efficiency must drop on rows where the fallback fired.

Whether the closed form is *optimal* on maximal PDAGs, or only valid, is
decided by ``experiments/lever0_differential.py`` against the enumerating
definition on the rows where the latter terminates, not asserted here.
"""

from __future__ import annotations

from dataclasses import dataclass

from bkrobust.demo.graph import MPDAG
from bkrobust.gac.mpdag_level import gac_forbidden_set, is_gac_valid_mpdag
from bkrobust.mpdag_criterion.criterion import is_amenable
from bkrobust.mpdag_criterion.paths import possible_descendants

Node = str


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
