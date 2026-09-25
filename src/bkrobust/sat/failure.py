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

GAC mode
--------

``add_failure`` also encodes GAC failure, selected with ``criterion="gac"``
(the default) versus ``criterion="backdoor"`` for the mode described above.
GAC failure in the witness DAG ``D`` mirrors
:func:`bkrobust.gac.dag_level.is_gac_valid_dag` exactly:

* ``cn[v]`` (``v != X``) -- ``v`` is a descendant of ``X`` in ``D`` (``rx[v]``)
  **and** ``v`` is ``Y`` or an ancestor of ``Y`` in ``D``. The second conjunct
  is a new topological fixpoint, ``ancY[v]``, built exactly like ``ancz`` but
  over the *full* ``D`` (not ``D_X̄``) and targeting ``{Y}``.
* ``forb[v]`` -- ``v`` is ``X``, ``Y``, in ``cn``, or a descendant (in ``D``)
  of a node already in ``forb``. Also a topological fixpoint, propagated
  forward along ``D``'s edges exactly as ``rx`` is.
* ``mode_a' = OR_{z in Z} forb[z]``.
* ``D'`` is ``D`` with each edge ``X -> v`` replaced by ``e[(X,v)] AND NOT
  cn[v]`` -- only the *first edges of proper causal paths* are cut, not every
  edge out of ``X`` as in the back-door mode. ``ancz'`` and the d-connecting
  path of mode ``B'`` are built exactly as their back-door counterparts, just
  over ``D'`` instead of ``D_X̄``.
* ``fails = mode_a' OR mode_b'``, sound by the same argument as
  :mod:`bkrobust.gac.dag_level`'s module docstring: the equivalence "``Z``
  blocks every proper non-causal path in ``D`` iff ``X ⟂ Y | Z`` in ``D'``"
  needs ``mode_a'`` false, and asserting the disjunction is correct either way
  (if ``mode_a'`` holds the set already fails; if not, the reduction applies).

Both modes reuse the same simple-path machinery (:func:`_add_dconn_mode`), so
model size does not double -- only the selected criterion's fixpoints and path
variables are built.
"""

from __future__ import annotations

import itertools
from collections.abc import Iterable
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


def _membership_vars(
    model: cp_model.CpModel, nodes: tuple[str, ...], members: Iterable[str], tag: str
) -> dict[str, Any]:
    """One fixed ``BoolVar`` per node, true iff the node is in ``members``.

    A plain ``model.Add(v == 0/1)`` rather than ``NewConstant``, so these
    behave identically to every other literal in the model (consistent with
    the rest of this module's style).
    """
    members_set = set(members)
    out = {}
    for v in nodes:
        bv = model.NewBoolVar(f"{tag}_{v}")
        model.Add(bv == (1 if v in members_set else 0))
        out[v] = bv
    return out


def _add_fixpoint(
    model: cp_model.CpModel,
    ov: OrientationVars,
    base: dict[str, Any],
    edge_of: Any,
    tag: str,
) -> dict[str, Any]:
    """A generic forward topological fixpoint.

    ``out[v] = base[v] OR OR_u out[u]`` over edges ``u -> v`` given by
    ``edge_of(u, v)`` (``None`` if absent).

    Safe because the model already constrains the underlying digraph acyclic
    (see the module docstring): there is a unique solution, evaluable in
    topological order, and the CP-SAT constraints below pin exactly that
    solution regardless of the order the Python loop visits vertices in.
    """
    nodes = ov.nodes
    out = {v: model.NewBoolVar(f"{tag}_{v}") for v in nodes}
    for v in nodes:
        supports = [base[v]]
        for u in nodes:
            if u == v:
                continue
            lit = edge_of(u, v)
            if lit is None:
                continue
            t = model.NewBoolVar(f"{tag}sup_{u}_{v}")
            model.AddBoolAnd([lit, out[u]]).OnlyEnforceIf(t)
            model.AddBoolOr([lit.Not(), out[u].Not()]).OnlyEnforceIf(t.Not())
            supports.append(t)
        model.AddMaxEquality(out[v], supports)
    return out


def _add_reverse_fixpoint(
    model: cp_model.CpModel,
    ov: OrientationVars,
    base: dict[str, Any],
    edge_of: Any,
    tag: str,
) -> dict[str, Any]:
    """Like :func:`_add_fixpoint` but propagating *backward* along edges.

    ``out[v] = base[v] OR OR_w out[w]`` over edges ``v -> w`` given by
    ``edge_of(v, w)``. Used for "has a descendant with property P" fixpoints
    (``ancz``, ``ancY``), which is why the recursion looks the other way
    around from ``rx``/``forb`` but is equally safe on an acyclic digraph.
    """
    nodes = ov.nodes
    out = {v: model.NewBoolVar(f"{tag}_{v}") for v in nodes}
    for v in nodes:
        supports = [base[v]]
        for w in nodes:
            if w == v:
                continue
            lit = edge_of(v, w)
            if lit is None:
                continue
            t = model.NewBoolVar(f"{tag}sup_{v}_{w}")
            model.AddBoolAnd([lit, out[w]]).OnlyEnforceIf(t)
            model.AddBoolOr([lit.Not(), out[w].Not()]).OnlyEnforceIf(t.Not())
            supports.append(t)
        model.AddMaxEquality(out[v], supports)
    return out


def _add_dconn_mode(
    model: cp_model.CpModel,
    ov: OrientationVars,
    x: str,
    y: str,
    z: frozenset[str],
    edge_of: Any,
    ancz: dict[str, Any],
    tag: str,
) -> Any:
    """An explicit simple ``X``-``Y`` path, d-connecting given ``Z``.

    Built over the graph whose ``u -> v`` presence is given by
    ``edge_of(u, v)`` (``None`` if absent). ``ancz[v]`` must already hold
    "``v`` in ``Z`` or has a descendant in ``Z``" over that same graph. Shared
    by both the back-door mode (over ``D_X̄``) and the GAC mode (over ``D'``)
    so the path machinery is written once.

    Returns:
        A literal true iff such a path exists.
    """
    nodes = ov.nodes
    n = len(nodes)
    q = {(i, v): model.NewBoolVar(f"q{tag}_{i}_{v}") for i in range(n) for v in nodes}
    act = [model.NewBoolVar(f"act{tag}_{i}") for i in range(n)]
    mode_b = model.NewBoolVar(f"mode_b{tag}")

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
        end = model.NewBoolVar(f"end{tag}_{i}")
        model.AddBoolAnd([act[i]] + ([act[i + 1].Not()] if i + 1 < n else [])).OnlyEnforceIf(end)
        model.AddImplication(end, q[(i, y)])
        model.AddBoolOr([act[i].Not()] + ([act[i + 1]] if i + 1 < n else []) + [end])
    # consecutive positions joined by an edge PRESENT in the graph
    for i in range(n - 1):
        for u, v in itertools.permutations(nodes, 2):
            both = [q[(i, u)], q[(i + 1, v)], act[i + 1]]
            fwd, bwd = edge_of(u, v), edge_of(v, u)
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
            into_l, into_r = edge_of(u, v), edge_of(w, v)
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
    return mode_b


def add_failure(
    model: cp_model.CpModel,
    ov: OrientationVars,
    e: dict[Edge, Any],
    x: str,
    y: str,
    z: frozenset[str],
    *,
    criterion: str = "gac",
) -> Any:
    """Constrain the witness so that ``Z`` is an INVALID adjustment set in it.

    Args:
        model: The CP-SAT model.
        ov: The MPDAG orientation variables.
        e: The witness DAG's edge variables, from :func:`add_witness_dag`.
        x: Treatment.
        y: Outcome.
        z: The candidate adjustment set.
        criterion: ``"gac"`` (default) encodes generalised-adjustment-criterion
            failure, matching :func:`bkrobust.gac.dag_level.is_gac_valid_dag`.
            ``"backdoor"`` encodes Pearl back-door failure, matching
            :func:`bkrobust.demo.evaluate.is_valid_adjustment_set_dag` (the
            original mode; see the module docstring). Model size is not
            affected by the choice not taken -- only the selected criterion's
            variables and constraints are built.

    Returns:
        A literal that is true exactly when ``Z`` fails in the witness DAG.
        (The model asserts it, so it is true in any solution.)
    """
    if criterion not in ("gac", "backdoor"):
        raise ValueError(f"unknown criterion: {criterion!r}")
    nodes = ov.nodes

    # ---- rx[v]: v is a descendant of X in D --------------------------------
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

    if criterion == "backdoor":
        # ---- (A) some z in Z is a descendant of X in D ---------------------
        mode_a = model.NewBoolVar("mode_a")
        zs = [rx[v] for v in sorted(z) if v in rx]
        if zs:
            model.AddMaxEquality(mode_a, zs)
        else:
            model.Add(mode_a == 0)

        # ---- D_Xbar: edges out of X deleted ---------------------------------
        def ep(u: str, v: str) -> Any | None:
            """``u -> v`` present in ``D_Xbar``; None when it cannot be there."""
            if u == x or not ov.is_adjacent(u, v):
                return None
            return e[(u, v)]

        # ancz[v]: v in Z, or v has a descendant in Z, in D_Xbar
        base_z = _membership_vars(model, nodes, z, "basez")
        ancz = _add_reverse_fixpoint(model, ov, base_z, ep, "ancz")

        # ---- (B) an explicit simple d-connecting path in D_Xbar -------------
        mode_b = _add_dconn_mode(model, ov, x, y, z, ep, ancz, "")

        model.AddBoolOr([mode_a, mode_b])
        fails = model.NewBoolVar("fails")
        model.Add(fails == 1)
        return fails

    # ==== GAC mode ===========================================================
    # ancY[v]: v is Y, or v is an ancestor of Y in D (full D, not D_Xbar).
    base_y = _membership_vars(model, nodes, {y}, "basey")

    def efull(u: str, v: str) -> Any | None:
        if not ov.is_adjacent(u, v):
            return None
        return e[(u, v)]

    ancy = _add_reverse_fixpoint(model, ov, base_y, efull, "ancy")

    # cn[v] (v != x): v descendant of X, and (v == y or v ancestor of Y).
    cn: dict[str, Any] = {}
    for v in nodes:
        t = model.NewBoolVar(f"cn_{v}")
        if v == x:
            model.Add(t == 0)
        elif v == y:
            model.Add(t == rx[v])
        else:
            model.AddBoolAnd([rx[v], ancy[v]]).OnlyEnforceIf(t)
            model.AddBoolOr([rx[v].Not(), ancy[v].Not()]).OnlyEnforceIf(t.Not())
        cn[v] = t

    # decn[v]: v is in cn, or a descendant (in D) of a node in cn -- i.e.
    # de(cn) union cn, seeded ONLY at cn (not at x or y, which are forbidden
    # as isolated points below, not as propagation sources -- forb is
    # cn u de(cn) u {x, y}, not de(x) u de(y) u {x, y}, which would silently
    # collapse this mode back to the back-door bar).
    decn = _add_fixpoint(model, ov, cn, efull, "decn")
    forb: dict[str, Any] = {}
    for v in nodes:
        t = model.NewBoolVar(f"forb_{v}")
        if v == x or v == y:
            model.Add(t == 1)
        else:
            model.Add(t == decn[v])
        forb[v] = t

    mode_a_p = model.NewBoolVar("mode_a_gac")
    zs = [forb[v] for v in sorted(z) if v in forb]
    if zs:
        model.AddMaxEquality(mode_a_p, zs)
    else:
        model.Add(mode_a_p == 0)

    # D': D with each edge X -> v replaced by e[(X,v)] AND NOT cn[v] -- only
    # the first edges of proper causal paths are cut.
    eprime: dict[Edge, Any] = {}
    for (u, v), var in e.items():
        if u == x:
            t = model.NewBoolVar(f"eprime_{u}_{v}")
            model.AddBoolAnd([var, cn[v].Not()]).OnlyEnforceIf(t)
            model.AddBoolOr([var.Not(), cn[v]]).OnlyEnforceIf(t.Not())
            eprime[(u, v)] = t
        else:
            eprime[(u, v)] = var

    def epp(u: str, v: str) -> Any | None:
        if not ov.is_adjacent(u, v):
            return None
        return eprime[(u, v)]

    base_zp = _membership_vars(model, nodes, z, "basezp")
    anczp = _add_reverse_fixpoint(model, ov, base_zp, epp, "anczp")

    mode_b_p = _add_dconn_mode(model, ov, x, y, z, epp, anczp, "_gac")

    model.AddBoolOr([mode_a_p, mode_b_p])
    fails = model.NewBoolVar("fails")
    model.Add(fails == 1)
    return fails
