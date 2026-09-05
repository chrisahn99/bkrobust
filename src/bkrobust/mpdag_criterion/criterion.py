r"""The adjustment criterion stated directly on an MPDAG -- no extension enumeration.

Why
---
Everything else in this repository decides "is ``Z`` a valid adjustment set for
``(X, Y)`` in the MPDAG ``G``?" by *enumerating* ``[G]``, the DAG extensions of
``G``, and checking each one (:func:`bkrobust.core.oracle.is_valid`). That is
exponential in the number of undirected edges and is the dominant cost of every
sweep in the project. The same question has a purely graphical answer -- the
generalised adjustment criterion of Perkovic, Textor, Kalisch and Maathuis, read
on ``G`` itself. This module implements that answer and is validated against the
enumerating oracle rather than trusted.

On cost, precisely
------------------
The decision is **polynomial in the size of the graph**. Every step of
:func:`why_invalid` is:

* the Meek closure of ``g``, memoised, and one further closure per undirected
  neighbour of ``x`` for amenability -- at most ``deg(x) + 1`` closures;
* :func:`backdoor_forbidden_set`, which is one ``O(V*E)`` edge-state search
  (:func:`~bkrobust.mpdag_criterion.paths.unshielded_reachable`);
* condition (c), one ``O(V*E)`` edge-state search
  (:func:`~bkrobust.mpdag_criterion.paths.open_definite_status_non_causal_path`),
  plus at most one memoised ``O(V*E)`` possible-descendant search per collider
  it meets, so ``O(V^2 * E)`` in the worst case.

That was not true of the first version of this module, and the difference is the
whole point of the module existing. Condition (c) used to be decided by walking
every simple path from ``x`` to ``y`` -- ``Theta((n-2)!)`` of them on a dense
graph. Measured on this repository's own dense Erdos-Renyi instances, that cost
up to 48 seconds for a *single* query at ``n = 12`` and blew a 60-second cap at
``n = 14``, on instances the enumerating oracle it replaces finished comfortably.
Condition (c) is a reachability question, not an enumeration question, for
exactly the reason d-separation is; see
:func:`~bkrobust.mpdag_criterion.paths.open_definite_status_non_causal_path` for
the state space and the argument.

Two exported functions remain exponential and are documented as such:
:func:`causal_nodes` and :func:`forbidden_set`. Neither is on the decision path
-- they are reporting aids naming the textbook objects, and they carry warnings.

Which predicate, exactly
------------------------
There are two different notions in play, and this module deliberately targets
the repository's one:

* the **generalised adjustment criterion** (GAC) characterises *complete*
  adjustment-set validity: ``sum_z P(y | x, z) P(z) = P(y | do(x))``;
* the repository's oracle uses ``is_valid_adjustment_set_dag``, which is
  **Pearl's back-door criterion**, and back-door is sufficient but *not*
  necessary for adjustment.

They genuinely differ. In ``X -> Y``, ``X -> W`` with ``W`` unconnected to
``Y``, the set ``Z = {W}`` is a valid adjustment set (``W`` is independent of
``Y`` given ``X``, so the sum telescopes) but violates the back-door criterion
because ``W`` is a descendant of ``X``. So a literal GAC implementation would
*disagree* with :func:`bkrobust.core.oracle.is_valid`, and this module would
fail its own acceptance test.

The reconciliation is exact and is what is implemented here. In a DAG ``D``:

* back-door valid  <=>  ``x, y not in Z``, ``Z & de(x, D) = {}``, and ``Z``
  blocks every back-door path;
* given ``Z & de(x, D) = {}``, "blocks every back-door path" is *equivalent* to
  "blocks every non-causal path": a non-causal path leaving ``x`` forwards has a
  first backward edge, and the node just before it is a collider lying in
  ``de(x, D)``, so with ``Z & de(x, D) = {}`` that collider blocks the path
  outright;
* and ``forb(x, y, D) \subseteq de(x, D) u {x}``.

So back-door validity in ``D`` is exactly GAC validity in ``D`` *plus*
``Z & (de(x, D) u {x, y}) = {}``. Quantifying over ``D in [G]`` and using
``possde(x, G) = union_{D in [G]} de(x, D)``, the criterion this module decides
is the GAC with the forbidden set enlarged from ``forb(x, y, G)`` to

    ``bd_forb(x, y, G) = possde(x, G) u {x, y}``.

Conditions (a) amenability and (c) "every proper definite-status non-causal path
is blocked" are untouched -- they are the GAC's own. Both forbidden sets are
exported: :func:`forbidden_set` is the textbook ``forb``, kept because it is the
object the literature names but carrying a substantial caveat of its own, and
:func:`backdoor_forbidden_set` is the one the oracle's back-door semantics
require and the one the decision actually uses.

Two definitions that had to be repaired
---------------------------------------
Both were found by the differential sweep, not by reading, and both are cases
where a definition that is exact on a CPDAG is *not* exact on an MPDAG carrying
background knowledge -- a shielded path can be possibly causal and yet be
realisable in no extension at all:

* ``possde`` must range over **unshielded** possibly causal paths, or the
  identity ``possde(x, G) = union_{D} de(x, D)`` used just above is false. See
  :func:`~bkrobust.mpdag_criterion.paths.possible_descendants`.
* amenability must **not** use that unshielded reduction, because short-circuiting
  a path can replace the undirected first edge that is the whole question. It is
  decided by a Meek closure per undirected neighbour of ``x`` instead. See
  :func:`is_amenable`.

Getting either one wrong produced disagreements with the oracle; getting both
right produced none.

Empty DAG-extension sets
------------------------
:func:`bkrobust.core.oracle.is_valid` returns ``False`` when ``[G]`` is empty --
a graph representing no model certifies nothing. A graphical criterion has no
natural way to notice this, because the criterion's conditions are all "for
every path" statements that a contradictory graph satisfies vacuously. The
policy here, stated so callers can rely on it:

* **Assumption.** ``g`` is a genuine MPDAG: Meek-closed, with ``[g]`` non-empty.
  Every element of :func:`bkrobust.search.space_fixed.build_space_correct` and
  every CPDAG satisfies this, which covers every use in this repository.
* **Sound partial detection.** Two cheap, *sound* (never wrong when they fire)
  emptiness proofs are still run: a directed cycle in ``g``, and a Meek-closure
  FAIL. Both imply ``[g] = {}`` because Meek's rules only force orientations
  that hold in every extension. When either fires, ``False`` /
  ``"no_extensions"`` is returned, matching the oracle.
* **Not complete.** A contradictory graph that survives both checks is answered
  as though it had extensions. The smallest witness is the chordless undirected
  4-cycle ``V0 - V1 - V2 - V3 - V0``: every orientation of it either closes a
  directed cycle or creates a new unshielded collider, so ``[G]`` is empty, yet
  it is acyclic and Meek-closed and this module reports ``z = {}`` as valid for
  ``(V0, V1)`` where the oracle reports False. That graph is not an MPDAG (it is
  not even a CPDAG -- chain components of an essential graph are chordal), so it
  is outside the stated assumption. Note the converse trap: chordality of the
  undirected subgraph is **not** a sound emptiness proof for MPDAGs and is
  deliberately not used, because background knowledge can leave a chordless
  undirected subgraph on a graph with many extensions -- see
  :mod:`bkrobust.search.space_fixed`, whose whole subject is that mistake.
  Detecting emptiness in general is the extension-existence problem this module
  exists to avoid.

As a side benefit of computing the closure, the criterion is evaluated on
``meek_closure(g)`` rather than ``g``. For a Meek-closed ``g`` these are the same
graph; for a non-closed one, ``[meek_closure(g)] = [g]`` whenever ``[g]`` is
non-empty, and the closed graph is the one on which "definite status" has its
intended meaning.

Written to run on Python 3.9 as well as the repository's 3.11 target: builtin
generics and ``X | Y`` unions appear only in annotations, and ``zip`` is never
called with ``strict=`` (3.10+).
"""

from __future__ import annotations

from collections.abc import Iterable

from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations, meek_closure
from bkrobust.mpdag_criterion.paths import (
    Node,
    Path,
    clear_index_cache,
    is_blocked,
    is_definite_status_path,
    is_non_causal,
    open_definite_status_non_causal_path,
    possible_descendants,
    possibly_causal_paths,
    simple_paths,
)

#: Reason strings returned by :func:`why_invalid`, in the priority order in
#: which they are tested. The first condition that fails names the verdict, so
#: the reasons are mutually exclusive by construction.
REASONS: tuple[str, ...] = (
    "no_extensions",
    "degenerate_query",
    "not_amenable",
    "z_hits_forbidden",
    "open_noncausal_path",
)

#: Memoised Meek closures, keyed on the MPDAG (which hashes on a canonical key,
#: so the cache changes cost and nothing else). ``None`` records a closure FAIL.
_CLOSURE_CACHE: dict[MPDAG, MPDAG | None] = {}


def clear_cache() -> None:
    """Drop every memo this package keeps. Only for tests and memory management.

    That is the Meek closures held here and the adjacency tables and
    possible-descendant sets held by
    :mod:`bkrobust.mpdag_criterion.paths`. All of them are pure functions of the
    graph, so clearing them changes cost and nothing else.
    """
    _CLOSURE_CACHE.clear()
    clear_index_cache()


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
    closed = meek_closure(g)
    _CLOSURE_CACHE[g] = closed
    return closed


# --------------------------------------------------------------------------------
# Amenability
# --------------------------------------------------------------------------------


def is_amenable(g: MPDAG, x: Node, y: Node) -> bool:
    """Whether ``g`` is amenable relative to ``(x, y)``.

    Amenability asks that ``g`` already know which edges at ``x`` are causal: no
    proper possibly causal path from ``x`` to ``y`` may begin with an undirected
    edge ``x - v``. Read semantically -- and this is the reading that has to be
    decided, since the criterion must agree with a for-all over extensions -- it
    says: **in no DAG extension does a causal path from ``x`` to ``y`` start
    with an edge that ``g`` leaves undirected.**

    Why it is a precondition rather than a nicety: if some extension has
    ``x -> v -> ... -> y`` causal with ``x - v`` undirected in ``g``, then the
    extension orienting ``v -> x`` instead turns that same path into a back-door
    path ``x <- v -> ... -> y`` whose interior nodes all lie in ``possde(x, g)``
    and are therefore barred from ``Z`` by the forbidden-set condition. No ``Z``
    can then be valid, so the whole query is refused up front.

    How it is decided
    -----------------
    For each undirected neighbour ``v`` of ``x``: impose ``x -> v`` and
    Meek-close. A FAIL means no extension orients that edge that way, so ``v``
    cannot start a causal path. Otherwise ask whether ``y`` is a possible
    descendant of ``v`` **in that closed graph** -- which, by the
    :func:`~bkrobust.mpdag_criterion.paths.possible_descendants` identity, holds
    iff some extension with ``x -> v`` also has ``v ~> y``. The sub-path from
    ``v`` to ``y`` cannot revisit ``x``, since ``x -> v`` and the extension is
    acyclic, so the two pieces really do concatenate into a path.

    Two simpler formulations were tried and are both wrong; they are recorded
    here because the mistake is easy to repeat:

    * "no possibly causal path from ``x`` to ``y`` starts undirected" -- too
      strict. It counts shielded paths that no extension can realise, and
      disagreed with the semantics on 702 of the 19,056 ``(graph, x, y)``
      queries in the n<=4 sweep scope.
    * "no *unshielded* possibly causal path starts undirected" -- too lax, and
      worse (1,380 disagreements). Unshielded reduction is valid for
      ``possde``, where only the endpoint matters, but **not** here, where the
      identity of the *first edge* is the whole question: short-circuiting
      ``x -> v -> w -> y`` to ``x -> w -> y`` can replace an undirected first
      edge with a directed one and destroy the witness. Minimal witness
      ``V0 -> V1, V0 - V2, V1 - V2`` with ``(x, y) = (V0, V1)``: the extension
      ``V0 -> V2 -> V1`` is causal and starts undirected, yet ``<V0, V2, V1>``
      is shielded.

    The formulation implemented here matched the semantics on all 19,056
    queries.

    Args:
        g: The MPDAG.
        x: The treatment node.
        y: The outcome node.

    Returns:
        True iff no DAG extension of ``g`` has a causal path from ``x`` to ``y``
        whose first edge is undirected in ``g``.
    """
    for v in sorted(g.neighbors(x)):
        forced = apply_orientations(g, [(x, v)])
        if forced is None:
            continue
        if y in possible_descendants(forced, v):
            return False
    return True


# --------------------------------------------------------------------------------
# Causal nodes and forbidden sets
# --------------------------------------------------------------------------------


def causal_nodes(g: MPDAG, x: Node, y: Node) -> set[Node]:
    """``cn(x, y, g)``: nodes on proper possibly causal paths from ``x`` to ``y``, minus ``x``.

    .. warning::
       **Exponential in the number of vertices, and deliberately left that way.**
       It enumerates possibly causal paths, of which a dense graph has
       factorially many. Nothing in :func:`is_valid_mpdag` or
       :func:`why_invalid` calls it -- it is a reporting aid, exported because
       it names an object the literature names. The reachability shortcut that
       would make it polynomial is *not* exact on an MPDAG (see below), so
       there is no drop-in replacement to switch to; do not call it at scale.

    Computed by enumerating the possibly causal paths rather than as
    ``possde(x) & (possan(y) u {y})``: in an MPDAG the concatenation of a
    possibly causal path into ``v`` with one out of ``v`` need not be a *simple*
    path, so the reachability shortcut that is exact on a DAG is not exact here.

    Over-inclusive on knowledge-carrying MPDAGs
    -------------------------------------------
    This is the definition as stated in the literature, and it is **not** the
    same as "nodes lying on a causal path from ``x`` to ``y`` in some DAG
    extension". A *shielded* possibly causal path need not be realisable in any
    extension -- the phenomenon documented at length under
    :func:`~bkrobust.mpdag_criterion.paths.possible_descendants` -- and this
    function counts its nodes anyway. In ``V0 -> V1, V0 -> V2, V0 - V3, V1 - V3``
    it reports ``cn(V1, V0) = {V3, V0}`` although no extension has any causal
    path from ``V1`` to ``V0``.

    The obvious repair (restrict to unshielded paths) is **wrong** for a
    different reason: short-circuiting a shielded path deletes interior nodes
    that genuinely do lie on a causal path in some extension. An exact version
    would need a joint-realisability test per candidate node, which is not
    attempted here. Nothing in :func:`is_valid_mpdag` depends on this function;
    it and :func:`forbidden_set` are reporting aids, carrying the caveat.

    Args:
        g: The MPDAG.
        x: The treatment node.
        y: The outcome node.

    Returns:
        The causal nodes, excluding ``x`` and including ``y`` when a possibly
        causal path exists at all.
    """
    out: set[Node] = set()
    for path in possibly_causal_paths(g, x, y):
        out.update(path[1:])
    return out


def forbidden_set(g: MPDAG, x: Node, y: Node) -> set[Node]:
    r"""``forb(x, y, g)``: possible descendants of the causal nodes, plus ``x``.

    .. warning::
       **Exponential in the number of vertices, and deliberately left that way.**
       It is built on :func:`causal_nodes` and inherits that function's cost as
       well as its caveat. Not on the decision path: :func:`is_valid_mpdag`
       uses the polynomial :func:`backdoor_forbidden_set`.

    This is the textbook forbidden set of the generalised adjustment criterion,
    implemented exactly as stated. It is **not** the set used by
    :func:`is_valid_mpdag`, which must match the repository's back-door oracle
    and therefore uses :func:`backdoor_forbidden_set`; see the module docstring.
    It is exported because it is the object the literature names.

    Caveat, and it is not a small one
    ---------------------------------
    Because it is built on :func:`causal_nodes`, this set inherits that
    function's over-inclusiveness on MPDAGs carrying background knowledge, and
    is therefore **not** contained in :func:`backdoor_forbidden_set` -- the
    containment that holds on a DAG (``forb(x, y, D) \subseteq de(x, D) u {x}``)
    and that the module docstring's derivation uses. It failed on 1,146 of the
    19,056 ``(graph, x, y)`` queries in the n<=4 sweep scope, 240 of them with an
    amenable graph, so amenability does not rescue it. Witness::

        V0 -> V1,  V0 -> V2,  V0 - V3,  V1 - V3      (x, y) = (V1, V0)

    ``cn`` reports ``{V3, V0}`` off the shielded path ``<V1, V3, V0>``, so
    ``forb`` sweeps in all four nodes, while ``bd_forb = {V0, V1, V3}``. No
    extension has a causal path from ``V1`` to ``V0`` at all. Quote this set in a
    write-up only with that caveat attached; the derivation in the module
    docstring is carried out at the DAG level, where the containment does hold,
    and does not depend on this function.

    Args:
        g: The MPDAG.
        x: The treatment node.
        y: The outcome node.

    Returns:
        ``possde(cn(x, y, g), g) u {x}``.
    """
    cn = causal_nodes(g, x, y)
    if not cn:
        return {x}
    return possible_descendants(g, cn) | {x}


def backdoor_forbidden_set(g: MPDAG, x: Node, y: Node) -> set[Node]:
    r"""The forbidden set under the repository's *back-door* semantics.

    ``possde(x, g) u {x, y}``. Pearl's back-door criterion bars every descendant
    of the treatment, not merely the descendants of the causal nodes, so at the
    MPDAG level it bars every *possible* descendant -- a node that is a
    descendant of ``x`` in even one DAG extension already invalidates ``Z``
    there, and validity is a for-all over extensions. ``y`` is included because
    the outcome may never be adjusted for.

    This *is* the set :func:`is_valid_mpdag` tests against. Note that it is not
    in general a superset of :func:`forbidden_set` on an MPDAG with background
    knowledge, even though the corresponding DAG-level containment holds -- see
    that function's caveat.

    Args:
        g: The MPDAG.
        x: The treatment node.
        y: The outcome node.

    Returns:
        ``possde(x, g) u {x, y}``.
    """
    return possible_descendants(g, x) | {x, y}


# --------------------------------------------------------------------------------
# Condition (c): open non-causal paths
# --------------------------------------------------------------------------------


def definite_status_non_causal_paths(g: MPDAG, x: Node, y: Node) -> list[Path]:
    """Every proper definite-status non-causal path from ``x`` to ``y``.

    .. warning::
       **Exponential in the number of vertices, and not used by the decision.**
       Condition (c) is decided by
       :func:`~bkrobust.mpdag_criterion.paths.open_definite_status_non_causal_path`,
       which searches edge states instead of paths. This function is kept as
       that search's slow reference and for reports that need the paths
       themselves; do not call it at scale.

    These are the paths condition (c) of the criterion quantifies over. Paths
    through a node of indefinite status are excluded on purpose: such a node is
    a collider in some extensions and a non-collider in others, and the
    criterion's completeness argument shows that whenever such a path is open in
    some extension, a *definite-status* non-causal path is open too -- so
    restricting to definite status loses nothing and is what makes the criterion
    decidable on ``g`` alone.

    Args:
        g: The MPDAG.
        x: The treatment node.
        y: The outcome node.

    Returns:
        The qualifying paths, deterministically ordered.
    """
    return [
        p for p in simple_paths(g, x, y) if is_non_causal(g, p) and is_definite_status_path(g, p)
    ]


def open_non_causal_path(g: MPDAG, x: Node, y: Node, z: Iterable[Node]) -> Path | None:
    """The first proper definite-status non-causal path from ``x`` to ``y`` left open by ``z``.

    .. warning::
       **Exponential in the number of vertices, and not used by the decision.**
       :func:`why_invalid` calls
       :func:`~bkrobust.mpdag_criterion.paths.open_definite_status_non_causal_path`
       instead, which answers the same question by state search and returns a
       shortest witness rather than the enumeration-order-first one. This
       function is that search's reference: the two are asserted to agree
       (as predicates) on the whole exhaustive ``n <= 4`` scope in
       ``tests/criterion/test_criterion.py``. The *witness* the two return may
       differ, since "first in a depth-first enumeration" and "shortest" are
       different choices; only the ``None`` / not-``None`` verdict is shared.

    Args:
        g: The MPDAG.
        x: The treatment node.
        y: The outcome node.
        z: The candidate adjustment set.

    Returns:
        A witness path, or ``None`` if every such path is blocked.
    """
    z_set = set(z)
    for path in definite_status_non_causal_paths(g, x, y):
        if not is_blocked(g, path, z_set):
            return path
    return None


# --------------------------------------------------------------------------------
# The criterion
# --------------------------------------------------------------------------------


def why_invalid(g: MPDAG, x: Node, y: Node, z: frozenset[Node]) -> str:
    """Why ``z`` fails as an adjustment set for ``(x, y)`` in ``g``, or ``""`` if it does not.

    Conditions are tested in a fixed priority order, so exactly one reason is
    ever returned and the reasons are mutually exclusive:

    1. ``"no_extensions"`` -- ``g`` is provably contradictory (directed cycle,
       or Meek's rules force a conflict), so ``[g]`` is empty and nothing can be
       certified. See the module docstring on the limits of this detection.
    2. ``"degenerate_query"`` -- ``x == y``. The oracle rejects this (a node is
       never d-separated from itself); there is no adjustment question to ask.
    3. ``"not_amenable"`` -- some proper possibly causal path from ``x`` to
       ``y`` starts with an undirected edge, so ``g`` does not determine which
       edges at ``x`` are causal. *No* ``z`` is valid in this case.
    4. ``"z_hits_forbidden"`` -- ``z`` meets
       :func:`backdoor_forbidden_set`: it contains ``x``, ``y``, or a possible
       descendant of ``x``.
    5. ``"open_noncausal_path"`` -- some proper definite-status non-causal path
       from ``x`` to ``y`` is left open by ``z``.

    The order matters for interpretation, not for the verdict: 3 and 4 are
    properties of ``(g, x, y)`` and of ``z`` respectively, and a case can fail
    both. Reporting amenability first means "this graph admits no valid set at
    all" is never mis-attributed to a bad choice of ``z``.

    Args:
        g: The MPDAG, assumed Meek-closed with non-empty ``[g]`` (see the module
            docstring).
        x: The treatment node.
        y: The outcome node.
        z: The candidate adjustment set.

    Returns:
        One of :data:`REASONS`, or ``""`` when ``z`` is valid.
    """
    closed = _closed(g)
    if closed is None:
        return "no_extensions"
    if x == y:
        return "degenerate_query"
    graph = closed
    if not is_amenable(graph, x, y):
        return "not_amenable"
    if set(z) & backdoor_forbidden_set(graph, x, y):
        return "z_hits_forbidden"
    if open_definite_status_non_causal_path(graph, x, y, z) is not None:
        return "open_noncausal_path"
    return ""


def is_valid_mpdag(g: MPDAG, x: Node, y: Node, z: frozenset[Node]) -> bool:
    """Whether ``z`` is a valid adjustment set for ``(x, y)`` in the MPDAG ``g``.

    Decided graphically on ``g``, without enumerating ``[g]``. Agrees with
    :func:`bkrobust.core.oracle.is_valid` -- verified exhaustively over every
    CPDAG on 3 and 4 nodes with an undirected edge, every element of every
    corrected space over those CPDAGs, every ordered ``(x, y)`` pair and every
    candidate ``z`` (see ``tests/criterion/test_criterion.py`` and
    ``results/axisb4/criterion_agreement.json``).

    The criterion, spelled out: ``z`` is valid iff

    (a) ``g`` is amenable relative to ``(x, y)`` (:func:`is_amenable`);
    (b) ``z`` is disjoint from :func:`backdoor_forbidden_set`; and
    (c) every proper definite-status non-causal path from ``x`` to ``y`` is
        blocked by ``z``.

    (b) is the back-door strengthening of the textbook ``forb`` condition; the
    module docstring derives why that -- and only that -- is the difference
    between this predicate and the generalised adjustment criterion.

    Args:
        g: The MPDAG, assumed Meek-closed with non-empty ``[g]``.
        x: The treatment node.
        y: The outcome node.
        z: The candidate adjustment set.

    Returns:
        True iff ``z`` is a valid back-door adjustment set for ``(x, y)`` in
        every DAG extension of ``g``.
    """
    return why_invalid(g, x, y, z) == ""


__all__ = [
    "REASONS",
    "backdoor_forbidden_set",
    "causal_nodes",
    "clear_cache",
    "definite_status_non_causal_paths",
    "forbidden_set",
    "is_amenable",
    "is_valid_mpdag",
    "open_non_causal_path",
    "why_invalid",
]
