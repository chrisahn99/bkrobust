r"""The missing back-door lower bound, and its exhaustive validation.

Background
----------

An adjustment set ``Z``, valid at ``G0``, stops being valid at some graph ``G``
above it for exactly one of two reasons (Pearl's back-door criterion, split
into its two clauses):

(A) some ``z in Z`` becomes a possible **descendant** of ``X`` -- handled by
    :func:`bkrobust.search.exact.descendant_lower_bound`;
(B) some **back-door path** from ``X`` to ``Y`` becomes **unblocked** given
    ``Z`` -- handled here.

:func:`bkrobust.search.exact.descendant_lower_bound` is admissible but
incomplete on its own: at n=4 (see ``results/search/bounds_n4.json``, commit
``b315704``) it was UNDEFINED on 60 of 792 finite-radius instances, precisely
those whose binding failure mode is (B). This module closes that gap.

The bound implemented here
---------------------------

A back-door path from ``X`` to ``Y`` is a simple path in the skeleton whose
first edge points *into* ``X``. Given ``Z``, it is blocked at an internal
(non-endpoint) vertex ``m`` iff:

* ``m`` is a **non-collider** on the path (chain or fork) and ``m in Z``, or
* ``m`` is a **collider** on the path and ``Z`` contains neither ``m`` nor any
  descendant of ``m``.

To *unblock* the path we therefore need, at every internal vertex ``m``:

* if ``m in Z``: ``m`` must be a **collider** (both path edges point into
  ``m``) -- this is necessary and sufficient, since ``m`` being in ``Z``
  itself discharges the "``m`` or a descendant of ``m`` is in ``Z``" clause;
* if ``m not in Z``: ``m`` must be a **non-collider**.

That second bullet is a deliberate simplification, stated here rather than
hidden. The exact criterion also allows ``m not in Z`` to work as a collider,
*provided some descendant of m is in Z* -- realising that needs a second
directed path (from ``m`` down to some ``z``), which is exactly the kind of
recursive sub-problem :func:`descendant_lower_bound` solves, but chaining it
here would require jointly verifying two paths share a consistent extension,
well beyond what either bound currently attempts. Dropping that option costs
completeness (some paths that could unblock via a collider-plus-descendant at
a non-``Z`` vertex are missed, i.e. some sub-cases of this failure mode
contribute no candidate), but costs nothing on **admissibility**: a path
excluded because we didn't model one of its unblocking strategies simply
contributes no candidate to the ``min`` below, rather than a wrong number.
Every candidate this function *does* return corresponds to a strategy it can
actually exhibit, so nothing is ever invented.

Under that restriction, orientations along a fixed simple path decompose into
a chain of local constraints: the desired state of edge ``i`` interacts only
with the desired state of edge ``i-1`` and ``i+1``, through the collider
condition at the vertex between them. That is solved by a small dynamic
program over the two possible states of each edge (see ``_path_cost``), the
same shortest-path idiom :func:`descendant_lower_bound` uses, generalised from
a chain of nodes to a chain of *edge states*.

Edge cost (identical convention to :func:`descendant_lower_bound`): orienting
an edge the way a path step demands costs 0 if it is already that way or is
still undirected in ``g0`` (free -- some extension can supply it), 1 if it
currently points the *other* way and is retractable, and is infeasible if it
points the other way and is data-compelled (frozen by the CPDAG).

``UNREACHED`` (see :mod:`bkrobust.core.conventions`) means this mode yields no
bound -- not that the cost is infinite. It is returned when no back-door path
in the skeleton admits a feasible orientation/blocking assignment at all.
"""

from __future__ import annotations

import itertools
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from bkrobust.core.conventions import UNREACHED
from bkrobust.core.oracle import is_valid
from bkrobust.core.resultsio import ResultWriter, write_manifest
from bkrobust.core.spacelib import Space, distances_from, radius
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.search.conjecture_study import all_cpdags
from bkrobust.search.exact import (
    descendant_lower_bound,
    guided_upper_bound,
    retractable_edges,
)
from bkrobust.search.space_fixed import build_space_correct

Node = str
Edge = tuple[str, str]

# --------------------------------------------------------------------------------
# Skeleton path enumeration
# --------------------------------------------------------------------------------


def _skeleton_adjacency(g: MPDAG) -> dict[Node, set[Node]]:
    """Undirected adjacency of ``g``'s skeleton, as a dict of sets."""
    adj: dict[Node, set[Node]] = {n: set() for n in g.nodes}
    for a, b in g.skeleton():
        adj[a].add(b)
        adj[b].add(a)
    return adj


def enumerate_simple_paths(g: MPDAG, x: Node, y: Node) -> list[list[Node]]:
    """Every simple path from ``x`` to ``y`` in ``g``'s skeleton, deterministically ordered.

    A path is a node list ``[x, ..., y]`` with no repeats. Includes the direct
    edge (a two-node path) when ``x`` and ``y`` are adjacent. Deterministic
    (sorted neighbour expansion) so results never depend on set iteration
    order or ``PYTHONHASHSEED``.

    At the node counts this module targets (``n <= 6``) brute-force DFS
    enumeration is cheap; nothing here is optimised for larger graphs.
    """
    if x == y:
        return []
    adj = _skeleton_adjacency(g)
    paths: list[list[Node]] = []
    visited = {x}
    path = [x]

    def dfs(node: Node) -> None:
        if node == y:
            paths.append(list(path))
            return
        for nxt in sorted(adj.get(node, ())):
            if nxt in visited:
                continue
            visited.add(nxt)
            path.append(nxt)
            dfs(nxt)
            path.pop()
            visited.discard(nxt)

    dfs(x)
    return paths


# --------------------------------------------------------------------------------
# Per-path cost: a small DP over edge orientation states
# --------------------------------------------------------------------------------


def _cost_to_orient(g0: MPDAG, retractable: set[Edge], tail: Node, head: Node) -> int | None:
    """Cost to have the ``tail``-``head`` edge of ``g0`` oriented ``tail -> head``.

    Same convention as the edge weights in
    :func:`bkrobust.search.exact.descendant_lower_bound`: 0 if already that way
    or still undirected, 1 if reversed and retractable, ``None`` (infeasible)
    if reversed and data-compelled.
    """
    if g0.is_directed_edge(tail, head) or g0.is_undirected_edge(tail, head):
        return 0
    if g0.is_directed_edge(head, tail):
        return 1 if (head, tail) in retractable else None
    return None  # not adjacent -- should not be reached on a skeleton path


def _path_cost(
    g0: MPDAG, retractable: set[Edge], z: frozenset[Node], path: list[Node]
) -> int | None:
    """Minimum retractions to make ``path`` an unblocked back-door path, or ``None``.

    ``path[0]`` is ``X``, ``path[-1]`` is ``Y``. Edge ``i`` joins ``path[i]``
    and ``path[i + 1]``; its state is ``'F'`` (``path[i] -> path[i+1]``) or
    ``'B'`` (``path[i+1] -> path[i]``). Edge 0 is forced to ``'B'`` (arrow into
    ``X``, the back-door requirement). At each internal vertex ``path[i]``
    (``1 <= i <= len(path) - 2``) the collider/non-collider requirement from
    the module docstring prunes the transition from edge ``i-1``'s state to
    edge ``i``'s state. Returns ``None`` if no state sequence survives.
    """
    m = len(path) - 1  # number of edges
    cost0 = _cost_to_orient(g0, retractable, path[1], path[0])  # edge 0 forced 'B'
    if cost0 is None:
        return None
    dp: dict[str, int] = {"B": cost0}

    for i in range(1, m):
        vi = path[i]
        must_be_z_collider = vi in z
        new_dp: dict[str, int] = {}
        for s_cur in ("F", "B"):
            if s_cur == "F":
                cost_edge = _cost_to_orient(g0, retractable, path[i], path[i + 1])
            else:
                cost_edge = _cost_to_orient(g0, retractable, path[i + 1], path[i])
            if cost_edge is None:
                continue
            best_prev: int | None = None
            for s_prev, prev_cost in dp.items():
                collider = s_prev == "F" and s_cur == "B"
                if collider != must_be_z_collider:
                    continue
                total = prev_cost + cost_edge
                if best_prev is None or total < best_prev:
                    best_prev = total
            if best_prev is not None:
                new_dp[s_cur] = best_prev
        dp = new_dp
        if not dp:
            return None

    return min(dp.values())


def backdoor_lower_bound(cpdag: MPDAG, g0: MPDAG, z: frozenset[Node], x: Node, y: Node) -> int:
    """Minimum retractions before some back-door path from ``X`` to ``Y`` can unblock given ``Z``.

    See the module docstring for the exact definition and the deliberate
    simplification (colliders at ``Z``-members only; non-``Z`` vertices must be
    non-colliders). Enumerates every simple skeleton path from ``x`` to ``y``,
    computes each path's minimum retraction cost by :func:`_path_cost`, and
    returns the minimum over paths that admit a feasible assignment at all.

    Args:
        cpdag: The CPDAG (fixes which edges are frozen).
        g0: The analyst's graph.
        z: The adjustment set.
        x: Treatment.
        y: Outcome.

    Returns:
        The minimum retraction count for this failure mode, or
        :data:`~bkrobust.core.conventions.UNREACHED` if no back-door path can
        be unblocked under the strategies this function models.

    Note:
        Like :func:`~bkrobust.search.exact.descendant_lower_bound`, this is a
        bound for **one** failure mode: admissible on its own, but only
        complete for the radius in combination with the other mode. Combine
        with :func:`combined_lower_bound`, never use alone as a certificate.
    """
    retractable = set(retractable_edges(cpdag, g0))
    best = UNREACHED
    for path in enumerate_simple_paths(cpdag, x, y):
        cost = _path_cost(g0, retractable, z, path)
        if cost is not None and (best == UNREACHED or cost < best):
            best = cost
    return best


def combined_lower_bound(cpdag: MPDAG, g0: MPDAG, z: frozenset[Node], x: Node, y: Node) -> int:
    """``min`` of the descendant- and back-door-mode lower bounds, over the modes that are defined.

    Both :func:`~bkrobust.search.exact.descendant_lower_bound` and
    :func:`backdoor_lower_bound` can independently return
    :data:`~bkrobust.core.conventions.UNREACHED`, meaning *that mode* gives no
    information -- never that its cost is 0 or infinite. An undefined mode is
    therefore dropped from the ``min``, not substituted with a value: treating
    it as 0 would make the combined bound trivial and wrong (0 is always
    "admissible" but useless), and treating it as infinity would make it
    dominate and silently discard a real bound from the other mode. If both
    modes are undefined the combined bound is undefined too.

    Returns:
        ``min(descendant_lower_bound, backdoor_lower_bound)`` over whichever of
        the two are not :data:`~bkrobust.core.conventions.UNREACHED`, or
        :data:`~bkrobust.core.conventions.UNREACHED` if neither is defined.
    """
    a = descendant_lower_bound(cpdag, g0, z, x)
    b = backdoor_lower_bound(cpdag, g0, z, x, y)
    candidates = [v for v in (a, b) if v != UNREACHED]
    if not candidates:
        return UNREACHED
    return min(candidates)


# --------------------------------------------------------------------------------
# Deliverable 3: a tightened upper bound by local search from guided_upper_bound
# --------------------------------------------------------------------------------


def improved_upper_bound(
    cpdag: MPDAG,
    g0: MPDAG,
    fails: Callable[[MPDAG], bool],
    z: frozenset[Node],
    x: Node,
    *,
    max_tries: int = 64,
) -> tuple[int, str | None]:
    """Try to shrink :func:`~bkrobust.search.exact.guided_upper_bound`'s witness.

    ``guided_upper_bound`` already searches retraction-subset sizes in
    increasing order and returns the first that fails, so within its own
    search it is already minimal *unless its ``max_tries`` budget was
    exhausted first* -- with many retractable edges the combinatorial subset
    count can run past the default budget before reaching (or exhausting) a
    given size. This function does two things beyond that:

    1. Re-runs the same construction with a materially larger budget, so size
       classes that were previously cut off get fully explored.
    2. Local "drop": given the smallest failing subset found, tries removing
       one retraction at a time (``k`` candidates, not another combinatorial
       search) to see whether a strict subset of the SAME witness also fails
       -- catching the case where the witness itself contains a redundant
       retraction that increasing-size search wouldn't have tried (increasing
       search only tries fresh combinations at each size, never revisits a
       specific found witness's own subsets beyond what step 1 already
       covers). Recurses while a drop still fails.

    Returns:
        ``(k, witness_edge_string)`` with ``k`` the smallest retraction-subset
        size found to fail, or ``(UNREACHED, None)`` if none was found within
        budget. Still an **upper bound**, not a certified radius: absence of a
        smaller witness within budget is not proof none exists.
    """
    retr = retractable_edges(cpdag, g0)
    from bkrobust.demo.meek import meek_closure

    def test(subset: tuple[Edge, ...]) -> MPDAG | None:
        g = g0
        for a, b in subset:
            if not g.is_directed_edge(a, b):
                return None
            g = g.unoriented(a, b)
        h = meek_closure(g)
        if h is None:
            return None
        return h if fails(h) else None

    best_subset: tuple[Edge, ...] | None = None
    tries = 0
    for k in range(1, len(retr) + 1):
        for subset in itertools.combinations(retr, k):
            tries += 1
            if tries > max_tries:
                break
            h = test(subset)
            if h is not None:
                best_subset = subset
                break
        if best_subset is not None or tries > max_tries:
            break

    if best_subset is None:
        return UNREACHED, None

    # Local "drop": remove one retraction at a time from the found witness.
    improved = True
    while improved and len(best_subset) > 1:
        improved = False
        for i in range(len(best_subset)):
            smaller = best_subset[:i] + best_subset[i + 1 :]
            h = test(smaller)
            if h is not None:
                best_subset = smaller
                improved = True
                break

    g = g0
    for a, b in best_subset:
        g = g.unoriented(a, b)
    h = meek_closure(g)
    witness = h.edge_string() if h is not None else None
    return len(best_subset), witness


# --------------------------------------------------------------------------------
# Deliverable 2: exhaustive admissibility + tightness + certification study
# --------------------------------------------------------------------------------


@dataclass
class StudyTotals:
    """Totals for one admissibility/tightness sweep."""

    n_cpdags: int = 0
    n_cases: int = 0
    n_finite_radius: int = 0
    n_l_defined: int = 0
    n_l_admissible: int = 0
    n_l_violations: int = 0
    n_l_tight: int = 0
    n_l_desc_defined: int = 0
    n_l_back_defined: int = 0
    n_l_back_only: int = 0  # backdoor defined, descendant UNREACHED: the closed gap
    n_u_valid: int = 0
    n_lu_equal_and_exact: int = 0
    violations: list[dict[str, Any]] = field(default_factory=list)


def _row_for_instance(
    cpdag: MPDAG,
    space: Space,
    dists: dict[MPDAG, int],
    g0: MPDAG,
    x: Node,
    y: Node,
    z: frozenset[Node],
) -> dict[str, Any] | None:
    """One study row for a single ``(cpdag, g0, x, y, z)`` instance, or ``None`` if degenerate.

    Degenerate here means ``Z`` is not valid at ``g0`` itself (gated out, same
    convention as :mod:`bkrobust.core.conventions`), which should not occur
    for ``Z = optimal_adjustment_set_mpdag(g0, x, y)`` but is checked rather
    than assumed.
    """
    if not is_valid(z, g0, x, y):
        return None

    def fails(g: MPDAG, _z: frozenset[Node] = z, _x: Node = x, _y: Node = y) -> bool:
        return not is_valid(_z, g, _x, _y)

    r, _witness = radius(space, dists, fails)
    la = descendant_lower_bound(cpdag, g0, z, x)
    lb = backdoor_lower_bound(cpdag, g0, z, x, y)
    lc_candidates = [v for v in (la, lb) if v != UNREACHED]
    lc = min(lc_candidates) if lc_candidates else UNREACHED
    u_res = guided_upper_bound(cpdag, g0, fails, z, x)
    u = u_res.radius

    if lc == UNREACHED:
        admissible: bool | None = None  # no claim made
    elif r == UNREACHED:
        admissible = False  # L claims a finite bound where nothing in the space fails at all
    else:
        admissible = lc <= r

    if u == UNREACHED:
        u_valid: bool | None = None
    elif r == UNREACHED:
        u_valid = False
    else:
        u_valid = u >= r

    return {
        "cpdag": cpdag.edge_string(),
        "g0": g0.edge_string(),
        "x": x,
        "y": y,
        "z": ",".join(sorted(z)),
        "r": r,
        "L_descendant": la,
        "L_backdoor": lb,
        "L_combined": lc,
        "U": u,
        "admissible": "" if admissible is None else str(admissible),
        "tight": str(lc != UNREACHED and r != UNREACHED and lc == r),
        "u_valid": "" if u_valid is None else str(u_valid),
        "exact_cert": str(lc != UNREACHED and u != UNREACHED and r != UNREACHED and lc == u == r),
    }


def run_admissibility_study(
    n: int = 4,
    out_dir: str | Path = "results/axisb2/bounds",
    *,
    resume: bool = False,
) -> StudyTotals:
    """Exhaustive admissibility/tightness/certification sweep over all n-node CPDAGs.

    For every CPDAG on ``n`` nodes, every element of its **corrected** space
    (:func:`bkrobust.search.space_fixed.build_space_correct` -- see that
    module for why the corrected builder, not
    :func:`bkrobust.core.spacelib.build_space`, is required) as ``G0``, every
    ordered ``(X, Y)``, and ``Z`` the optimal adjustment set at ``G0`` where it
    is identified and valid: computes the exact BFS radius, the descendant and
    back-door lower bounds and their combination, and the guided upper bound.
    Writes one row per instance plus a totals summary, incrementally.

    Because this uses the corrected space (more elements per CPDAG than the
    space the 792-instance, 83.3%-certification baseline in
    ``results/search/bounds_n4.json`` was computed against -- see that file's
    originating commit, which predates the space correction), the instance
    count here will not equal 792 even at ``n=4``; the certification *rate* is
    the number to compare, with the denominator shift stated plainly rather
    than glossed over.

    Args:
        n: Node count.
        out_dir: Output directory. Gets ``rows.csv``, ``totals.json`` and
            ``manifest.json``.
        resume: Passed through to :class:`~bkrobust.core.resultsio.ResultWriter`.

    Returns:
        The :class:`StudyTotals` for the sweep (also written to
        ``totals.json``).
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    writer = ResultWriter(out / "rows.csv", resume=resume)
    totals = StudyTotals()

    cpdags = all_cpdags(n)
    totals.n_cpdags = len(cpdags)

    for cpdag in cpdags:
        space = build_space_correct(cpdag)
        for g0 in space.elements:
            dists = distances_from(space, g0)
            for x, y in itertools.permutations(cpdag.nodes, 2):
                o = optimal_adjustment_set_mpdag(g0, x, y)
                if o is None:
                    continue
                z = frozenset(o)
                row = _row_for_instance(cpdag, space, dists, g0, x, y, z)
                if row is None:
                    continue
                writer.write(row)
                totals.n_cases += 1

                r = row["r"]
                la, lb, lc, u = (
                    row["L_descendant"],
                    row["L_backdoor"],
                    row["L_combined"],
                    row["U"],
                )
                if r != UNREACHED:
                    totals.n_finite_radius += 1
                if la != UNREACHED:
                    totals.n_l_desc_defined += 1
                if lb != UNREACHED:
                    totals.n_l_back_defined += 1
                if lb != UNREACHED and la == UNREACHED:
                    totals.n_l_back_only += 1
                if lc != UNREACHED:
                    totals.n_l_defined += 1
                    if row["admissible"] == "True":
                        totals.n_l_admissible += 1
                    elif row["admissible"] == "False":
                        totals.n_l_violations += 1
                        totals.violations.append(row)
                    if row["tight"] == "True":
                        totals.n_l_tight += 1
                if u != UNREACHED and row["u_valid"] != "False":
                    totals.n_u_valid += 1
                if row["exact_cert"] == "True":
                    totals.n_lu_equal_and_exact += 1

    writer.close()

    write_manifest(
        out,
        seed=0,
        grid={"n": n},
        extra={
            "purpose": "Axis B, session 2: back-door lower bound admissibility study",
            "space": "build_space_correct (see space_fixed.py)",
            "baseline": {
                "source": "results/search/bounds_n4.json (commit b315704)",
                "space": "uncorrected build_space -- predates space_fixed.py",
                "n_cases": 792,
                "LU_equal_and_exact": 660,
                "certification_rate": 660 / 792,
            },
            "totals": totals_as_dict(totals),
        },
    )
    (out / "totals.json").write_text(_totals_json(totals))
    return totals


def totals_as_dict(t: StudyTotals) -> dict[str, Any]:
    """Flatten :class:`StudyTotals` for JSON, dropping the (possibly large) violation list."""
    d = {
        "n_cpdags": t.n_cpdags,
        "n_cases": t.n_cases,
        "n_finite_radius": t.n_finite_radius,
        "n_l_defined": t.n_l_defined,
        "n_l_admissible": t.n_l_admissible,
        "n_l_violations": t.n_l_violations,
        "n_l_tight": t.n_l_tight,
        "n_l_desc_defined": t.n_l_desc_defined,
        "n_l_back_defined": t.n_l_back_defined,
        "n_l_back_only_closed_gap": t.n_l_back_only,
        "n_u_valid": t.n_u_valid,
        "n_lu_equal_and_exact": t.n_lu_equal_and_exact,
        "certification_rate": (
            t.n_lu_equal_and_exact / t.n_finite_radius if t.n_finite_radius else None
        ),
        "baseline_certification_rate_uncorrected_space": 660 / 792,
        "n_violations_recorded": len(t.violations),
    }
    return d


def _totals_json(t: StudyTotals) -> str:
    import json

    payload = totals_as_dict(t)
    payload["violations"] = t.violations[:50]  # cap: a bug report, not a full dump
    return json.dumps(payload, indent=2, default=str)
