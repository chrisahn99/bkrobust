r"""The failure predicate: a witness DAG, encoded existentially.

Why a witness DAG at all
------------------------

``Z`` is *valid* in an MPDAG ``G`` iff it satisfies the adjustment criterion in
**every** DAG of ``[G]`` -- a universal quantifier, which alone would put the
problem in QBF. But the radius is defined by *failure*, and failure is
``∃ D ∈ [G]`` with ``Z`` invalid in ``D`` -- a plain existential. So the model
carries the witness DAG as first-class variables beside the MPDAG orientation
variables and asks the solver for both at once. That is what keeps this in SAT,
and the witness the solver returns is the falsifying scenario the project wants
to report anyway.

Matching the reference oracle exactly
-------------------------------------

The encoding must agree with ``bkrobust.demo.evaluate.is_valid_adjustment_set_dag``,
which is the definition everything else in the repository is calibrated against.
That function declares ``Z`` invalid in ``D`` when either

* **(A)** ``X`` or ``Y`` lies in ``Z``, or some ``z ∈ Z`` is a descendant of
  ``X`` in ``D``; or
* **(B)** ``X`` and ``Y`` are d-connected given ``Z`` in ``D_X̄``, the graph with
  the edges *out of* ``X`` deleted.

Both are encoded below, over ``D`` for (A) and over ``D_X̄`` for (B). The
distinction matters: in ``D_X̄`` the treatment has no descendants at all, so
using the wrong graph for (A) would silently disable that mode.

Soundness of the path encoding, and why not a fixpoint
------------------------------------------------------

d-connection is naturally expressed as reachability over (vertex, direction)
states -- the Bayes-ball recursion. Encoded as a boolean fixpoint that recursion
is **unsafe in both directions**: the state graph has cycles (descend to a
collider, ascend again), so a set of mutually-supporting states can be true with
no grounding in ``X`` (over-reporting failure, radii too small), while adding
level variables to forbid that lets the solver instead leave a genuinely
reachable state false (under-reporting failure, radii too large). Under-reporting
is the more dangerous of the two here, because the ladder's UNSAT answers are
what certify robustness.

So the d-connecting path is encoded **explicitly**, as a simple path of bounded
length with one position-indexed variable block. Every d-connecting path can be
taken simple, so bounding the length by the vertex count loses nothing, and the
encoding is exact for SAT *and* UNSAT.

Two fixpoints do remain, and both are safe because they range over a graph the
model already constrains to be **acyclic**, where the recursion has a unique
solution evaluable in topological order:

* ``rx[v]``   -- ``v`` is a descendant of ``X`` in ``D``;
* ``ancz[v]`` -- ``v`` is in ``Z`` or has a descendant in ``Z``, in ``D_X̄``.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import Any

from ortools.sat.python import cp_model

from bkrobust.demo.graph import canon
from bkrobust.sat.closure import OrientationVars

Edge = tuple[str, str]


@dataclass
class WitnessVars:
    """The witness DAG and the literal asserting ``Z`` fails in it."""

    e: dict[Edge, Any]
    fails: Any


def add_witness_dag(model: cp_model.CpModel, ov: OrientationVars) -> dict[Edge, Any]:
    """A DAG ``D`` extending the MPDAG ``G``, with no v-structure absent from ``Ĉ``.

    Args:
        model: The CP-SAT model.
        ov: The MPDAG orientation variables.

    Returns:
        ``(u, v) -> BoolVar`` meaning ``u -> v`` in ``D``.
    """
    nodes = ov.nodes
    e: dict[Edge, Any] = {}
    for a, b in sorted(ov.adjacent):
        e[(a, b)] = model.NewBoolVar(f"e_{a}_{b}")
        e[(b, a)] = model.NewBoolVar(f"e_{b}_{a}")
        model.AddExactlyOne([e[(a, b)], e[(b, a)]])  # D orients every edge
        model.AddImplication(ov.d[(a, b)], e[(a, b)])  # D extends G
        model.AddImplication(ov.d[(b, a)], e[(b, a)])

    # acyclicity
    epos = {v: model.NewIntVar(0, len(nodes) - 1, f"epos_{v}") for v in nodes}
    for (u, v), var in e.items():
        model.Add(epos[u] < epos[v]).OnlyEnforceIf(var)

    # no v-structure absent from the CPDAG
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
            if ov.is_adjacent(a, c) or (canon(a, c), b) in cpdag_v:
                continue
            model.AddBoolOr([e[(a, b)].Not(), e[(c, b)].Not()])
    return e


def add_failure(
    model: cp_model.CpModel,
    ov: OrientationVars,
    e: dict[Edge, Any],
    x: str,
    y: str,
    z: frozenset[str],
) -> Any:
    """Constrain the witness so that ``Z`` is an INVALID adjustment set in it.

    Returns:
        A literal that is true exactly when ``Z`` fails in the witness DAG.
        (The model asserts it, so it is true in any solution.)
    """
    nodes = ov.nodes
    n = len(nodes)

    # ---- (A) some z in Z is a descendant of X in D -------------------------
    rx = {v: model.NewBoolVar(f"rx_{v}") for v in nodes}
    model.Add(rx[x] == 0)
    for v in nodes:
        if v == x:
            continue
        supports = []
        for u in nodes:
            if u == v or not ov.is_adjacent(u, v):
                continue
            t = model.NewBoolVar(f"rxsup_{u}_{v}")
            if u == x:
                model.AddBoolAnd([e[(u, v)]]).OnlyEnforceIf(t)
                model.AddBoolOr([e[(u, v)].Not()]).OnlyEnforceIf(t.Not())
            else:
                model.AddBoolAnd([e[(u, v)], rx[u]]).OnlyEnforceIf(t)
                model.AddBoolOr([e[(u, v)].Not(), rx[u].Not()]).OnlyEnforceIf(t.Not())
            supports.append(t)
        if supports:
            model.AddMaxEquality(rx[v], supports)
        else:
            model.Add(rx[v] == 0)
    mode_a = model.NewBoolVar("mode_a")
    zs = [rx[v] for v in sorted(z) if v in rx]
    if zs:
        model.AddMaxEquality(mode_a, zs)
    else:
        model.Add(mode_a == 0)

    # ---- D_Xbar: edges out of X deleted -----------------------------------
    def ep(u: str, v: str) -> Any | None:
        """``u -> v`` present in ``D_Xbar``; None when the edge cannot be there."""
        if u == x or not ov.is_adjacent(u, v):
            return None
        return e[(u, v)]

    # ancz[v]: v in Z, or v has a descendant in Z, in D_Xbar
    ancz = {v: model.NewBoolVar(f"ancz_{v}") for v in nodes}
    for v in nodes:
        if v in z:
            model.Add(ancz[v] == 1)
            continue
        supports = []
        for w in nodes:
            lit = ep(v, w)
            if lit is None:
                continue
            t = model.NewBoolVar(f"anczsup_{v}_{w}")
            model.AddBoolAnd([lit, ancz[w]]).OnlyEnforceIf(t)
            model.AddBoolOr([lit.Not(), ancz[w].Not()]).OnlyEnforceIf(t.Not())
            supports.append(t)
        if supports:
            model.AddMaxEquality(ancz[v], supports)
        else:
            model.Add(ancz[v] == 0)

    # ---- (B) an explicit simple d-connecting path from X to Y in D_Xbar ----
    q = {(i, v): model.NewBoolVar(f"q_{i}_{v}") for i in range(n) for v in nodes}
    act = [model.NewBoolVar(f"act_{i}") for i in range(n)]
    mode_b = model.NewBoolVar("mode_b")

    model.AddImplication(mode_b, act[0])
    model.AddImplication(mode_b, act[1])  # at least X and Y
    for i in range(n - 1):
        model.AddImplication(act[i + 1], act[i])  # positions form a prefix
    for i in range(n):
        model.Add(sum(q[(i, v)] for v in nodes) == 1).OnlyEnforceIf(act[i])
        model.Add(sum(q[(i, v)] for v in nodes) == 0).OnlyEnforceIf(act[i].Not())
    for v in nodes:  # simple path
        model.Add(sum(q[(i, v)] for i in range(n)) <= 1)
    model.AddImplication(mode_b, q[(0, x)])
    # the last active position is Y
    for i in range(1, n):
        end = model.NewBoolVar(f"end_{i}")
        model.AddBoolAnd([act[i]] + ([act[i + 1].Not()] if i + 1 < n else [])).OnlyEnforceIf(end)
        model.AddImplication(end, q[(i, y)])
        model.AddBoolOr([act[i].Not()] + ([act[i + 1]] if i + 1 < n else []) + [end])
    # consecutive positions joined by an edge PRESENT in D_Xbar
    for i in range(n - 1):
        for u, v in itertools.permutations(nodes, 2):
            both = [q[(i, u)], q[(i + 1, v)], act[i + 1]]
            fwd, bwd = ep(u, v), ep(v, u)
            opts = [lit for lit in (fwd, bwd) if lit is not None]
            if not opts:
                model.AddBoolOr([lit.Not() for lit in both])
            else:
                model.AddBoolOr([lit.Not() for lit in both] + opts)
    # collider / non-collider conditions at internal positions
    for i in range(1, n - 1):
        for u, v, w in itertools.permutations(nodes, 3):
            if not (ov.is_adjacent(u, v) and ov.is_adjacent(v, w)):
                continue
            ctx = [q[(i - 1, u)], q[(i, v)], q[(i + 1, w)], act[i + 1], mode_b]
            into_l, into_r = ep(u, v), ep(w, v)
            neg = [lit.Not() for lit in ctx]
            if into_l is not None and into_r is not None:
                # collider at v  =>  v must have a descendant in Z
                model.AddBoolOr([*neg, into_l.Not(), into_r.Not(), ancz[v]])
                if v in z:
                    # v is in Z, so the path is open at v ONLY as a collider,
                    # which needs BOTH arrows pointing in. Requiring only their
                    # disjunction lets a chain through Z pass as a collider --
                    # that was a real bug, caught by the stage-2 predicate test
                    # (144 of 26,304 cases, all over-reporting failure).
                    model.AddBoolOr([*neg, into_l])
                    model.AddBoolOr([*neg, into_r])
            elif v in z:
                # v cannot be a collider here, so v in Z blocks the path
                model.AddBoolOr(neg)

    model.AddBoolOr([mode_a, mode_b])
    fails = model.NewBoolVar("fails")
    model.Add(fails == 1)
    return fails
