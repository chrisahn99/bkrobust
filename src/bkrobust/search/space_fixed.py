r"""The CORRECTED perturbation space, and the resolution of Task 0.

Summary of the defect
---------------------

The inherited :func:`bkrobust.demo.space.enumerate_space` filters candidates
through ``is_valid_mpdag``, which requires each connected component of the
**undirected subgraph** to be chordal. That condition characterises *CPDAGs*
(essential graphs are chain graphs with chordal chain components). It is **not**
a property of an MPDAG obtained by adding background knowledge to a CPDAG.

Adding knowledge can place a directed edge *inside* what would otherwise be one
undirected component. The component then contains a partially directed cycle, so
the graph is not a chain graph, and its undirected subgraph can look chordless
because the chord that would fix it is present but **directed**. Such a graph is
nonetheless a perfectly ordinary knowledge state.

Evidence that the excluded graphs are legitimate (all verified computationally,
see ``results/axisb2/task0_*.json``):

* each is **reachable**: it equals ``Meek(cpdag, K)`` for an explicit consistent
  ``K``;
* each is **Meek-closed** under R1-R4;
* each represents a **non-empty** set of DAGs;
* each is **maximally oriented** -- for every undirected edge, both orientations
  occur among its DAG extensions, which is the defining property of an MPDAG.

So the chordality test was wrong for this purpose, not the notion of validity.

The corrected membership test
-----------------------------

``G`` is a knowledge state of ``Ĉ`` exactly when closing ``Ĉ`` under ``G``'s own
non-compelled orientations returns ``G``:

    ``apply_orientations(cpdag, G.directed \\ cpdag.directed) == G``

This fixpoint predicate was verified to hold for **every** reachable state and
to characterise reachability exactly (1,588 states over all CPDAGs on 3 and 4
nodes).

Scale of the effect
-------------------

The previous session reported 0.09% *of Meek closures*. Measured per element of
the space, the loss is far larger and **grows with density**: 0% for k <= 4
undirected edges, 3.85% at k=5 (affecting 100% of CPDAGs), 8.25% at k=7 and
10.19% at k=8. Every number in the previous session computed against
``enumerate_space`` therefore moves, and is re-run here side by side.

Why the old function is left alone
----------------------------------

:func:`bkrobust.demo.space.enumerate_space` is **not** modified, so the previous
session's committed results remain exactly reproducible from the code that
produced them. New work uses the functions here.
"""

from __future__ import annotations

import itertools

from bkrobust.core.oracle import extensions
from bkrobust.core.spacelib import Space
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.space import covering_pairs, neighbour_graph

Edge = tuple[str, str]


def knowledge_of(cpdag: MPDAG, g: MPDAG) -> list[Edge]:
    """The orientations ``g`` adds to ``cpdag`` -- its background knowledge."""
    return sorted(set(g.directed_edges) - set(cpdag.directed_edges))


def is_knowledge_state(cpdag: MPDAG, g: MPDAG) -> bool:
    """Whether ``g`` is a reachable knowledge state of ``cpdag``.

    The fixpoint test described in the module docstring: closing ``cpdag`` under
    ``g``'s own non-compelled orientations must return ``g`` itself.
    """
    h = apply_orientations(cpdag, knowledge_of(cpdag, g))
    return h is not None and h == g


def enumerate_space_correct(cpdag: MPDAG) -> list[MPDAG]:
    """Every knowledge state of ``cpdag``, deterministically ordered.

    Enumerates the ``3^k`` assignments of orientation state to the ``k``
    undirected edges, Meek-closes each, and deduplicates. This is reachability by
    construction, so it needs no validity predicate at all -- which is what makes
    it immune to the chordality defect.

    Args:
        cpdag: The CPDAG to refine.

    Returns:
        The knowledge states, sorted by (number of directed edges, edge string).
    """
    und = sorted(cpdag.undirected_edges)
    seen: dict[str, MPDAG] = {}
    for assign in itertools.product((0, 1, 2), repeat=len(und)):
        ors: list[Edge] = []
        for (a, b), v in zip(und, assign):  # noqa: B905 - equal by construction
            if v == 1:
                ors.append((a, b))
            elif v == 2:
                ors.append((b, a))
        g = apply_orientations(cpdag, ors)
        if g is not None:
            seen.setdefault(g.edge_string(), g)
    return sorted(seen.values(), key=lambda g: (len(g.directed_edges), g.edge_string()))


def build_space_correct(cpdag: MPDAG) -> Space:
    """Build the corrected space with its order, covers and neighbour graph.

    Mirrors :func:`bkrobust.core.spacelib.build_space` but over the corrected
    element set. The covering relation is still computed from model inclusion by
    enumerating represented DAGs -- no conjectured characterisation.
    """
    elements = tuple(enumerate_space_correct(cpdag))
    reps = {g: frozenset(extensions(g)) for g in elements}
    covers = frozenset(covering_pairs(list(elements), reps))
    nbrs = neighbour_graph(list(elements), covers)
    return Space(cpdag=cpdag, elements=elements, reps=reps, covers=covers, neighbours=nbrs)


def is_maximally_oriented(g: MPDAG) -> bool | None:
    """Whether every undirected edge of ``g`` is genuinely undecided in ``[g]``.

    The defining property of a maximally oriented PDAG, and the check that
    settled Task 0: the graphs the old filter excluded all satisfy it.

    Returns:
        ``True`` / ``False``, or ``None`` when ``[g]`` is empty and the question
        is vacuous.
    """
    exts = enumerate_dag_extensions(g)
    if not exts:
        return None
    for a, b in sorted(g.undirected_edges):
        fwd = any(d.is_directed_edge(a, b) for d in exts)
        bwd = any(d.is_directed_edge(b, a) for d in exts)
        if not (fwd and bwd):
            return False
    return True
