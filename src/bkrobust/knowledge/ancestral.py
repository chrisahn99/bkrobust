"""Reduce an ancestral claim to the orientation claims it entails, and to nothing more.

An ancestral claim, "``a`` causes ``b``, possibly through other variables", says
that ``a`` is an ancestor of ``b`` in the true DAG. It is the one granularity of
background knowledge that a causal order over adjacent pairs does not already
express: it names a directed path, not an edge, and on a CPDAG it may pin down
one edge, several, or none at all.

The instruments downstream consume orientation claims on undirected edges of
the CPDAG. So an ancestral claim has to be converted, and the conversion must
be *sound*: an orientation is emitted only when every DAG of the equivalence
class in which ``a`` is an ancestor of ``b`` contains it. The rule implemented
here is sound and deliberately incomplete, so that no orientation is ever a
guess.

**The rule.** Let ``P`` be the set of possibly directed paths from ``a`` to
``b`` in the CPDAG: simple paths along which every edge is either undirected
or directed in the direction of travel. A directed path of any DAG in the
class is, read on the CPDAG, a member of ``P``, so if ``a`` is an ancestor of
``b`` the true path is one of them.

* If every path in ``P`` starts with the same edge ``a - c`` and that edge is
  undirected in the CPDAG, the claim entails ``a -> c``.
* If every path in ``P`` ends with the same edge ``c - b`` and that edge is
  undirected in the CPDAG, the claim entails ``c -> b``.
* Otherwise the claim is ``not_reducible``: the paths diverge at both ends,
  or the unique end edges are already compelled, or ``P`` is empty, in which
  case the claim contradicts the CPDAG. Nothing is emitted.

Both end edges can fire at once, and Meek's rules then close whatever lies
between them. The rule does not chase entailments that only appear once the
divergent paths are jointly considered; the unit test records one such graph,
where two divergent last edges are in fact both entailed and the rule stays
silent. Silence costs coverage, never validity.

The first-edge set is computed without enumerating paths: ``a - c`` starts a
path in ``P`` if and only if ``c`` is a possibly directed neighbour of ``a``
and ``b`` is reachable from ``c`` along possibly directed edges in the graph
with ``a`` removed. The last-edge set is the mirror image. Linear in the size
of the graph per claim.
"""

from __future__ import annotations

from dataclasses import dataclass

from bkrobust.demo.graph import MPDAG, Edge, Node

REDUCED = "reduced"
NOT_REDUCIBLE = "not_reducible"


@dataclass(frozen=True)
class AncestralReduction:
    """What one ancestral claim ``a`` causes ``b`` entails on the CPDAG.

    Attributes:
        a: The claimed ancestor.
        b: The claimed descendant.
        orientations: The entailed orientation claims as ``(tail, head)``; empty
            when the claim is not reducible.
        status: ``"reduced"`` or ``"not_reducible"``.
        reason: Why: ``unique_first_edge``, ``unique_last_edge``,
            ``unique_first_and_last_edge``, ``no_possibly_directed_path``,
            ``end_edges_diverge`` or ``end_edges_already_directed``.
        n_first_edges: How many distinct first edges the paths in ``P`` use.
        n_last_edges: How many distinct last edges they use.
    """

    a: Node
    b: Node
    orientations: tuple[Edge, ...]
    status: str
    reason: str
    n_first_edges: int
    n_last_edges: int


def _forward(g: MPDAG, node: Node) -> set[Node]:
    """Nodes one possibly directed step from ``node``: undirected neighbours and children."""
    return g.neighbors(node) | g.children(node)


def _backward(g: MPDAG, node: Node) -> set[Node]:
    """Nodes one possibly directed step *into* ``node``: undirected neighbours and parents."""
    return g.neighbors(node) | g.parents(node)


def _reaches(g: MPDAG, start: Node, target: Node, avoid: Node) -> bool:
    """Whether ``target`` is reachable from ``start`` possibly-directedly, never via ``avoid``."""
    if start == avoid:
        return False
    seen = {start}
    stack = [start]
    while stack:
        u = stack.pop()
        if u == target:
            return True
        for v in _forward(g, u):
            if v != avoid and v not in seen:
                seen.add(v)
                stack.append(v)
    return False


def first_edges(g: MPDAG, a: Node, b: Node) -> set[Node]:
    """Nodes ``c`` such that some possibly directed path from ``a`` to ``b`` starts ``a - c``."""
    return {c for c in _forward(g, a) if _reaches(g, c, b, avoid=a)}


def last_edges(g: MPDAG, a: Node, b: Node) -> set[Node]:
    """Nodes ``c`` such that some possibly directed path from ``a`` to ``b`` ends ``c - b``."""
    return {c for c in _backward(g, b) if _reaches(g, a, c, avoid=b)}


def reduce_ancestral_claim(g: MPDAG, a: Node, b: Node) -> AncestralReduction:
    """Apply the sound reduction rule to the claim that ``a`` is an ancestor of ``b``.

    Args:
        g: The CPDAG the claim is read against. Not mutated.
        a: The claimed ancestor.
        b: The claimed descendant.

    Returns:
        The entailed orientations, or none, with the reason.

    Raises:
        ValueError: If ``a`` and ``b`` are the same node or not in ``g``.
    """
    if a == b or a not in g.nodes or b not in g.nodes:
        raise ValueError(f"ancestral claim needs two distinct nodes of the graph: {(a, b)}")
    firsts = first_edges(g, a, b)
    lasts = last_edges(g, a, b)
    if not firsts:
        return AncestralReduction(a, b, (), NOT_REDUCIBLE, "no_possibly_directed_path", 0, 0)
    out: list[Edge] = []
    reasons: list[str] = []
    if len(firsts) == 1:
        (c,) = firsts
        if g.is_undirected_edge(a, c):
            out.append((a, c))
            reasons.append("first")
    if len(lasts) == 1:
        (c,) = lasts
        if g.is_undirected_edge(c, b) and (c, b) not in out:
            out.append((c, b))
            reasons.append("last")
    if out:
        reason = "unique_" + "_and_".join(reasons) + "_edge"
        return AncestralReduction(a, b, tuple(out), REDUCED, reason, len(firsts), len(lasts))
    reason = (
        "end_edges_already_directed" if len(firsts) == 1 or len(lasts) == 1 else "end_edges_diverge"
    )
    return AncestralReduction(a, b, (), NOT_REDUCIBLE, reason, len(firsts), len(lasts))
