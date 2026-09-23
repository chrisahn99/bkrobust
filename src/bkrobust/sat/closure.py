r"""Meek closure as a set of clauses, and the correctness argument for it.

The encoding represents an MPDAG ``G`` refining the CPDAG ``Ĉ`` by one boolean
per **ordered** adjacent pair:

    ``d[u, v] = 1``  iff  ``u -> v`` is a directed edge of ``G``.

Edges already compelled in ``Ĉ`` are fixed; the remaining ("free") edges are the
ones background knowledge may orient. ``d[u,v]`` and ``d[v,u]`` are mutually
exclusive; both false means the edge is left undirected.

Meek closure as forbidden configurations
---------------------------------------

A graph is Meek-closed exactly when **no rule fires**, so each rule instance is
encoded as a clause forbidding its firing configuration. The applicable vertex
tuples are fixed by the skeleton and enumerated once at build time. Writing
``und[u,v]`` for "the edge is present and left undirected":

* **R1** forbid ``d[a,b] & und[b,c]``           (``a``, ``c`` non-adjacent)
* **R2** forbid ``d[a,c] & d[c,b] & und[a,b]``
* **R3** forbid ``und[a,b] & und[a,c] & und[a,e] & d[c,b] & d[e,b]``
  (``c``, ``e`` non-adjacent)
* **R4** forbid ``und[a,b] & und[a,c] & d[c,w] & d[w,b]``
  (``c``, ``b`` non-adjacent)

An earlier draft encoded R1-R4 as bare implications over the ``d`` variables,
dropping the premises that certain edges are **undirected**. That is sound for R1
and R2 -- there the only other option creates a new v-structure or a directed
cycle, both independently forbidden -- but **wrong for R3 and R4**, whose
premises genuinely require ``a-c`` and ``a-d`` to be undirected. The over-strong
clause set forbade legitimate knowledge states; the stage-1 validation below
caught it immediately, with 66 of 133 CPDAGs mismatching. The rule forms above
match ``bkrobust.demo.meek`` line for line.

Two further families are needed, since Meek-closedness alone does not pin down
membership of the space:

* **no new v-structure**: forbid ``d[a,b] & d[c,b]`` for non-adjacent ``a``,
  ``c``, except where that v-structure is already present in ``Ĉ``;
* **acyclicity**, via a position variable per vertex.

**None of this is taken on trust.** :func:`solutions` enumerates every model and
the accompanying test checks it equals the corrected space exactly, exhaustively
at small ``n``.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from typing import Any

from ortools.sat.python import cp_model

from bkrobust.demo.graph import MPDAG, canon

Node = str
Edge = tuple[str, str]


@dataclass
class OrientationVars:
    """The ``d[u,v]`` booleans and the structural facts derived from ``Ĉ``.

    Attributes:
        cpdag: The CPDAG being refined.
        nodes: Vertices in sorted order (determinism).
        d: ``(u, v) -> BoolVar`` for every ordered adjacent pair.
        free: The undirected edges of ``Ĉ``, canonically ordered.
        adjacent: Unordered adjacency lookup.
    """

    cpdag: MPDAG
    nodes: tuple[str, ...]
    d: dict[Edge, Any]
    free: tuple[Edge, ...]
    adjacent: set[Edge] = field(default_factory=set)

    def is_adjacent(self, u: str, v: str) -> bool:
        """Whether ``u`` and ``v`` are adjacent in the CPDAG."""
        return canon(u, v) in self.adjacent


def add_orientation_vars(model: cp_model.CpModel, cpdag: MPDAG) -> OrientationVars:
    """Create the ``d[u,v]`` booleans and pin the compelled edges.

    Args:
        model: The CP-SAT model.
        cpdag: The CPDAG.

    Returns:
        An :class:`OrientationVars`.
    """
    nodes = tuple(cpdag.nodes)
    adjacent = {canon(a, b) for a, b in cpdag.skeleton()}
    d: dict[Edge, Any] = {}
    for a, b in sorted(adjacent):
        d[(a, b)] = model.NewBoolVar(f"d_{a}_{b}")
        d[(b, a)] = model.NewBoolVar(f"d_{b}_{a}")
        # an edge is oriented at most one way
        model.AddBoolOr([d[(a, b)].Not(), d[(b, a)].Not()])

    for a, b in sorted(cpdag.directed_edges):
        model.Add(d[(a, b)] == 1)
        model.Add(d[(b, a)] == 0)

    free = tuple(sorted(cpdag.undirected_edges))
    return OrientationVars(cpdag=cpdag, nodes=nodes, d=d, free=free, adjacent=adjacent)


def add_undirected_vars(model: cp_model.CpModel, ov: OrientationVars) -> dict[Edge, Any]:
    """One boolean per edge: the edge is present and left undirected."""
    und: dict[Edge, Any] = {}
    for a, b in sorted(ov.adjacent):
        u = model.NewBoolVar(f"und_{a}_{b}")
        # u <-> (not d[a,b]) and (not d[b,a])
        model.AddBoolOr([ov.d[(a, b)], ov.d[(b, a)], u])
        model.AddImplication(u, ov.d[(a, b)].Not())
        model.AddImplication(u, ov.d[(b, a)].Not())
        und[canon(a, b)] = u
    return und


def add_meek_closure(model: cp_model.CpModel, ov: OrientationVars) -> int:
    """Forbid every Meek-rule firing configuration, plus new v-structures and cycles.

    Returns:
        The number of clauses added.
    """
    d, nodes = ov.d, ov.nodes
    und = add_undirected_vars(model, ov)
    n_clauses = 0

    def forbid(lits: list[Any]) -> None:
        nonlocal n_clauses
        model.AddBoolOr([lit.Not() for lit in lits])
        n_clauses += 1

    # R1 and R2
    for a, b, c in itertools.permutations(nodes, 3):
        ab, bc, ac = ov.is_adjacent(a, b), ov.is_adjacent(b, c), ov.is_adjacent(a, c)
        if ab and bc and not ac:
            forbid([d[(a, b)], und[canon(b, c)]])
        if ac and bc and ab:
            forbid([d[(a, c)], d[(c, b)], und[canon(a, b)]])

    for a, b in itertools.permutations(nodes, 2):
        if not ov.is_adjacent(a, b):
            continue
        others = [w for w in nodes if w not in (a, b)]
        # R3
        for c, e in itertools.combinations(others, 2):
            if ov.is_adjacent(c, e):
                continue
            if (
                ov.is_adjacent(a, c)
                and ov.is_adjacent(a, e)
                and ov.is_adjacent(c, b)
                and ov.is_adjacent(e, b)
            ):
                forbid([und[canon(a, b)], und[canon(a, c)], und[canon(a, e)], d[(c, b)], d[(e, b)]])
        # R4
        for c in others:
            if ov.is_adjacent(c, b) or not ov.is_adjacent(a, c):
                continue
            for w in others:
                if w == c or not (ov.is_adjacent(c, w) and ov.is_adjacent(w, b)):
                    continue
                forbid([und[canon(a, b)], und[canon(a, c)], d[(c, w)], d[(w, b)]])

    # no NEW v-structure: a->b<-c with a,c non-adjacent, unless already in the CPDAG
    cpdag_v = set()
    for b in nodes:
        ps = sorted(ov.cpdag.parents(b))
        for i, a in enumerate(ps):
            for c in ps[i + 1 :]:
                if not ov.is_adjacent(a, c):
                    cpdag_v.add((canon(a, c), b))
    for b in nodes:
        nb = [w for w in nodes if ov.is_adjacent(w, b)]
        for a, c in itertools.combinations(nb, 2):
            if ov.is_adjacent(a, c):
                continue
            if (canon(a, c), b) in cpdag_v:
                continue
            forbid([d[(a, b)], d[(c, b)]])

    # acyclicity via positions
    pos = {v: model.NewIntVar(0, len(nodes) - 1, f"pos_{v}") for v in nodes}
    for (u, v), var in d.items():
        model.Add(pos[u] < pos[v]).OnlyEnforceIf(var)
    return n_clauses


def solutions(cpdag: MPDAG, limit: int | None = None) -> list[frozenset[Edge]]:
    """Enumerate every model of the closure clauses, as sets of directed edges.

    Used only for validation -- it is the brute-force check that the clause set
    describes exactly the corrected space.

    Args:
        cpdag: The CPDAG.
        limit: Stop after this many solutions.

    Returns:
        One frozenset of directed ``(tail, head)`` pairs per model, sorted.
    """
    model = cp_model.CpModel()
    ov = add_orientation_vars(model, cpdag)
    add_meek_closure(model, ov)

    collected: list[frozenset[Edge]] = []
    keys = sorted(ov.d)

    class _Collect(cp_model.CpSolverSolutionCallback):
        def __init__(self) -> None:
            super().__init__()

        def on_solution_callback(self) -> None:
            got = frozenset(k for k in keys if self.Value(ov.d[k]))
            collected.append(got)
            if limit is not None and len(collected) >= limit:
                self.StopSearch()

    solver = cp_model.CpSolver()
    solver.parameters.num_workers = 1
    solver.parameters.random_seed = 0
    solver.parameters.enumerate_all_solutions = True
    solver.Solve(model, _Collect())
    return sorted(collected, key=lambda s: sorted(s))
