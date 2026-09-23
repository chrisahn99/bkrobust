r"""The generalised adjustment criterion lifted to an MPDAG, in polynomial time.

Why
---
:mod:`bkrobust.gac.dag_level` decides the GAC on a fully oriented DAG. The
question this repository actually asks is the MPDAG one -- "is ``Z`` a valid
adjustment set for ``(x, y)`` **in every DAG the partially oriented graph
represents**?" -- and the honest way to answer it is

.. code-block:: python

    all(is_gac_valid_dag(d, x, y, z) for d in extensions(g))

which is exponential in the number of undirected edges and is the cost that
dominates every sweep in the project. This module answers the same question
graphically, from ``g`` alone, in polynomial time. The enumeration above is not
a fallback here; it is the **semantics**, and the acceptance test for this
module is exact agreement with it (``tests/gac/test_gac.py``,
``results/axisa2/gac_agreement.json``).

Relation to :mod:`bkrobust.mpdag_criterion`
-------------------------------------------
Session 4 built the same machinery for the repository's *back-door* oracle. The
two criteria share their shape --

(a) ``g`` is amenable relative to ``(x, y)``;
(b) ``z`` misses a forbidden set;
(c) every proper definite-status non-causal path from ``x`` to ``y`` is blocked
    by ``z``

-- and (a) and (c) are **identical**, reused verbatim from
:mod:`bkrobust.mpdag_criterion`. The single difference is (b):

* back-door:  ``bd_forb(x, y, G) = possde(x, G) u {x, y}``;
* GAC:        ``forb(x, y, G) = possde(cn(x, y, G), G) u {x, y}``.

Everything below is about computing that one set correctly and cheaply.

The forbidden set, and why it is *not* ``possde`` of the textbook ``cn``
------------------------------------------------------------------------
What the semantics needs is the union over extensions,

    ``FORB = union_{D in [G]} forb(x, y, D)``
          ``= {x, y} u {v : exists D, exists c in cn(x, y, D), v in de(c, D)}``,

because condition (b) has to reproduce ``for all D: z ∩ forb(x, y, D) = {}``.
Two candidate formulas were implemented and measured against that union before
the one below was kept.

*Rejected: nodes on possibly causal paths.* The literature's ``cn(x, y, G)`` --
nodes on proper possibly causal paths -- is what
:func:`bkrobust.mpdag_criterion.criterion.causal_nodes` computes, and that
function is both exponential and (per its own docstring) over-inclusive on
MPDAGs that carry background knowledge. It is not usable here.

*Rejected: concatenating reachability.* The DAG formula
``cn = de(x) & (an(y) u {y})`` lifts naively to
``{c != x : c in possde(x, G) and y in possde(c, G)}``. That is polynomial but
**wrong**, because the two memberships can require *incompatible* extensions.
Minimal witness, found by the sweep at n=3::

    V0 -> V1,  V0 - V2        (x, y) = (V0, V1),  z = {V2}

``V2 in possde(V0)`` needs ``V0 -> V2``; ``V1 in possde(V2)`` needs the path
``V2 -> V0 -> V1``, i.e. ``V2 -> V0``. No single extension has both, and indeed
no extension has ``V2`` on a causal path from ``V0`` to ``V1`` at all -- yet the
naive formula forbids ``V2``, rejecting a set that is valid in both extensions.
This formula disagreed with the semantics on 6 of the 552 cases at n=3.

*Kept: anchor the causal node at a child of ``x``.* Suppose ``v in FORB``, so
some extension ``D`` has ``c in cn(x, y, D)`` with ``v in de(c, D) u {c}``. Let
``w`` be the **first node after ``x``** on the directed path ``x ~> c`` in ``D``
(so ``w = c`` when ``c`` is already a child). Then ``w`` is a child of ``x`` in
``D``, and since ``w`` is an ancestor of ``c`` it inherits both of ``c``'s
properties: ``y in de(w, D) u {w}`` and ``v in de(w, D) u {w}``. So the
quantifier over ``c`` can be replaced by a quantifier over **children of ``x``**
without losing anything. Moreover, **on an amenable graph** every such ``w``
is a child of ``x`` in ``G`` itself: amenability says exactly that no causal
path from ``x`` to ``y`` in any extension starts with an edge ``g`` leaves
undirected. Hence

    ``FORB(x, y, G) = {x, y} u union { possde(w, G) : w in ch(x, G),
                                       y in possde(w, G) }``

with ``possde`` the inclusive possible-descendant set of
:func:`bkrobust.mpdag_criterion.paths.possible_descendants` (the *unshielded*
one, which is the one that equals ``union_D de(., D)``; see that function for
why the naive walk rule is wrong here).

This is an over-approximation in principle -- ``y in possde(w)`` and
``v in possde(w)`` are still two separate existential statements, and nothing
above forces one extension to witness both -- and that residual gap is
**measured rather than assumed**: zero disagreements with the enumerated
semantics over the exhaustive n<=3 and n<=4 scope (74,568 cases) and over a
2.4M-case sampled n=5 stress pass. Should a counterexample ever be found, it
would be a case where ``z`` is reported invalid though it is valid in every
extension -- the conservative direction.

Cost
----
Polynomial, and this is a hot-path requirement, not a nicety (session 4's first
criterion was exhaustively verified *and* exponential; see
:mod:`bkrobust.mpdag_criterion.criterion`). Per call: one memoised Meek closure;
at most ``deg(x) + 1`` further closures for amenability; at most ``deg(x)``
memoised ``O(V*E)`` possible-descendant searches for the forbidden set; and one
``O(V*E)`` edge-state search for condition (c) plus one memoised
possible-descendant search per collider it meets. Nothing on the decision path
enumerates a path or an extension. Measured cold-cache timings on dense
Erdos-Renyi instances are in ``results/axisa2/gac_agreement.json``.

Empty extension sets
--------------------
Same policy, and same limits, as :mod:`bkrobust.mpdag_criterion.criterion`: a
directed cycle or a Meek-closure FAIL is a sound proof that ``[g]`` is empty and
is reported as ``"no_extensions"`` (False), matching the enumerating semantics;
detection is *not* complete, and callers are assumed to pass genuine MPDAGs
(every CPDAG and every element of
:func:`bkrobust.search.space_fixed.build_space_correct` qualifies).

Written to run on Python 3.9 as well as the repository's 3.11 target: builtin
generics and ``X | Y`` unions appear only in annotations, and neither
``zip(strict=)`` nor ``itertools.pairwise`` (3.10+) is used.
"""

from __future__ import annotations

from collections.abc import Iterable

from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import meek_closure
from bkrobust.mpdag_criterion.criterion import clear_cache as _clear_criterion_cache
from bkrobust.mpdag_criterion.criterion import is_amenable
from bkrobust.mpdag_criterion.paths import (
    open_definite_status_non_causal_path,
    possible_descendants,
)

Node = str

#: Reason strings returned by :func:`why_invalid_gac`, in the priority order in
#: which they are tested. The first condition that fails names the verdict, so
#: the reasons are mutually exclusive by construction. Deliberately the same
#: vocabulary as :data:`bkrobust.mpdag_criterion.criterion.REASONS`, since the
#: two criteria fail in the same five ways and only condition (b) differs.
GAC_REASONS: tuple[str, ...] = (
    "no_extensions",
    "degenerate_query",
    "not_amenable",
    "z_hits_forbidden",
    "open_noncausal_path",
)

#: Memoised Meek closures, keyed on the MPDAG (which hashes on a canonical key,
#: so the cache changes cost and nothing else). ``None`` records a closure FAIL,
#: which is a sound proof that ``[g]`` is empty.
_CLOSURE_CACHE: dict[MPDAG, MPDAG | None] = {}

#: Entries retained before :func:`_closed` drops the whole table. A plain
#: high-water mark rather than an LRU, matching
#: :mod:`bkrobust.mpdag_criterion.paths`: a perturbation search walks very many
#: distinct graphs, and an unbounded table would be a slow memory leak.
_CLOSURE_CACHE_LIMIT = 4096


def clear_cache() -> None:
    """Drop every memo this decision depends on. Only for tests and memory management.

    That is this module's Meek closures plus the caches held by
    :mod:`bkrobust.mpdag_criterion` (its own closures, the adjacency tables and
    the possible-descendant sets). All of them are pure functions of the graph,
    so clearing them changes cost and nothing else.
    """
    _CLOSURE_CACHE.clear()
    _clear_criterion_cache()


def _closed(g: MPDAG) -> MPDAG | None:
    """Return ``meek_closure(g)``, memoised; ``None`` on a closure FAIL.

    Args:
        g: The MPDAG.

    Returns:
        The Meek-closed graph, or ``None`` -- which is a sound proof that
        ``[g]`` is empty.
    """
    if g in _CLOSURE_CACHE:
        return _CLOSURE_CACHE[g]
    if len(_CLOSURE_CACHE) >= _CLOSURE_CACHE_LIMIT:
        _CLOSURE_CACHE.clear()
    closed = meek_closure(g)
    _CLOSURE_CACHE[g] = closed
    return closed


def gac_causal_children(g: MPDAG, x: Node, y: Node) -> set[Node]:
    """Children of ``x`` in ``g`` that can start a causal path to ``y``.

    The anchors of :func:`gac_forbidden_set`: the nodes ``w`` with ``x -> w``
    **directed in** ``g`` and ``y in possde(w, g)``, i.e. those for which some
    DAG extension has ``x -> w ~> y``. The module docstring argues why
    restricting the causal-node quantifier to children of ``x`` loses nothing on
    an amenable graph, and why the more obvious "all of ``possde(x)``" version
    is wrong.

    Args:
        g: The MPDAG, assumed Meek-closed.
        x: The treatment node.
        y: The outcome node.

    Returns:
        The qualifying children, possibly empty.
    """
    return {w for w in sorted(g.children(x)) if y in possible_descendants(g, w)}


def gac_forbidden_set(g: MPDAG, x: Node, y: Node) -> set[Node]:
    """``forb(x, y, g)``: the GAC's forbidden set, lifted to an MPDAG.

    ``{x, y} u union of possde(w, g) over w in`` :func:`gac_causal_children`.
    Intended to equal ``union over D in [g] of forb(x, y, D)``, which is what
    condition (b) must reproduce; the derivation, the two rejected formulas and
    the measured evidence are in the module docstring.

    Contrast with :func:`bkrobust.mpdag_criterion.criterion.backdoor_forbidden_set`,
    which is ``possde(x, g) u {x, y}``: this set is the smaller one whenever
    ``x`` has possible descendants that cannot reach ``y``, and that difference
    is the whole behavioural gap between the GAC and back-door at MPDAG level.

    Cost: at most ``deg(x)`` memoised ``O(V*E)`` searches.

    Args:
        g: The MPDAG, assumed Meek-closed.
        x: The treatment node.
        y: The outcome node.

    Returns:
        The forbidden set, always containing ``x`` and ``y``.
    """
    out: set[Node] = {x, y}
    for w in sorted(gac_causal_children(g, x, y)):
        out |= possible_descendants(g, w)
    return out


def why_invalid_gac(g: MPDAG, x: Node, y: Node, z: Iterable[Node]) -> str:
    """Why ``z`` fails the GAC for ``(x, y)`` in ``g``, or ``""`` if it does not.

    Conditions are tested in the fixed priority order of :data:`GAC_REASONS`, so
    exactly one reason is ever returned and the reasons are mutually exclusive:

    1. ``"no_extensions"`` -- ``g`` is provably contradictory (a directed cycle,
       or Meek's rules force a conflict), so ``[g]`` is empty and nothing can be
       certified. Detection is sound but not complete; see the module docstring.
    2. ``"degenerate_query"`` -- ``x == y``. There is no adjustment question to
       ask, and the enumerated semantics rejects it too.
    3. ``"not_amenable"`` -- some DAG extension has a causal path from ``x`` to
       ``y`` starting with an edge ``g`` leaves undirected, so ``g`` does not
       determine which edges at ``x`` are causal. **No** ``z`` is valid then.
    4. ``"z_hits_forbidden"`` -- ``z`` meets :func:`gac_forbidden_set`: it
       contains ``x``, ``y``, or a node that is on or below a causal path from
       ``x`` to ``y`` in some extension.
    5. ``"open_noncausal_path"`` -- some proper definite-status non-causal path
       from ``x`` to ``y`` is left open by ``z``.

    The order matters for interpretation, not for the verdict: 3 and 4 are
    properties of ``(g, x, y)`` and of ``z`` respectively and a case can fail
    both. Reporting amenability first means "this graph admits no valid set at
    all" is never mis-attributed to a bad choice of ``z``. Note that 3 must be
    tested before 4 for a second reason as well: the argument that
    :func:`gac_forbidden_set` may anchor on children of ``x`` in ``g`` assumes
    amenability.

    Args:
        g: The MPDAG, assumed to have a non-empty ``[g]`` (see the module
            docstring). It need not be Meek-closed; it is closed internally.
        x: The treatment node.
        y: The outcome node.
        z: The candidate adjustment set.

    Returns:
        One of :data:`GAC_REASONS`, or ``""`` when ``z`` is valid.
    """
    closed = _closed(g)
    if closed is None:
        return "no_extensions"
    if x == y:
        return "degenerate_query"
    if not is_amenable(closed, x, y):
        return "not_amenable"
    if set(z) & gac_forbidden_set(closed, x, y):
        return "z_hits_forbidden"
    if open_definite_status_non_causal_path(closed, x, y, z) is not None:
        return "open_noncausal_path"
    return ""


def is_gac_valid_mpdag(g: MPDAG, x: Node, y: Node, z: Iterable[Node]) -> bool:
    """Whether ``z`` satisfies the GAC for ``(x, y)`` in *every* DAG extension of ``g``.

    Decided graphically on ``g``, in polynomial time, without enumerating
    ``[g]``. The predicate this implements is exactly

    .. code-block:: python

        bool(extensions(g)) and all(is_gac_valid_dag(d, x, y, z) for d in extensions(g))

    including the empty-``[g]`` convention (a graph representing no model
    certifies nothing, so the answer is False), subject to the detection limits
    in the module docstring.

    Args:
        g: The MPDAG, assumed to have a non-empty ``[g]``.
        x: The treatment node.
        y: The outcome node.
        z: The candidate adjustment set.

    Returns:
        True iff ``z`` is a valid adjustment set for ``(x, y)`` in every DAG
        extension of ``g``.
    """
    return why_invalid_gac(g, x, y, z) == ""


__all__ = [
    "GAC_REASONS",
    "clear_cache",
    "gac_causal_children",
    "gac_forbidden_set",
    "is_gac_valid_mpdag",
    "why_invalid_gac",
]
