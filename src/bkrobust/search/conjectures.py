"""The two structural conjectures, and the machinery to test them exhaustively.

Conjecture 1 -- failure is upward-closed
---------------------------------------

    If ``Z`` fails in ``G`` and ``[G] subset-or-equal [G']``, then ``Z`` fails
    in ``G'``.

**This is a theorem, not a conjecture, and the proof is one line.** Validity is
a for-all over represented DAGs: ``is_valid(Z, G)`` iff every ``D in [G]`` has
``Z`` valid. So ``Z`` failing in ``G`` means some ``D in [G]`` has ``Z``
invalid. Since ``[G] subset-or-equal [G']``, that same ``D`` lies in ``[G']``,
witnessing failure there. QED.

There is exactly one caveat, and it is an artefact of a convention rather than
of the mathematics. This implementation defines ``is_valid`` to return ``False``
when ``[G]`` is **empty** -- a graph representing no model cannot certify
anything. Under that convention an empty-extension graph "fails" vacuously,
while ``[] subset-or-equal [G']`` holds for every ``G'``, so the implication
would break for any valid ``G'``. It does not bite here because
``is_valid_mpdag`` requires at least one extension and ``enumerate_space``
filters on it, so **no element of the space has empty extensions**. That is
checked empirically rather than assumed (:func:`check_no_empty_extensions`).

Consequences, which the search exploits: the failing set is an **up-set** in the
model-inclusion order; its boundary is an antichain of *minimal* failures; and
validity can be propagated by order comparison instead of recomputed per
element.

Conjecture 2 -- retraction-optimal witnesses
--------------------------------------------

    The nearest failure is always reachable from ``G0`` by a monotone **upward**
    path -- pure retractions -- so

        ``r_val`` = the minimum number of retractions from ``G0`` inducing failure.

This one is genuinely open and is *not* implied by Conjecture 1. Upward-closure
says failures persist as you move up; it does **not** say the nearest failure
lies above ``G0``. A failure incomparable to ``G0`` could sit at distance 2 while
the nearest failure in the up-set sits at distance 3, and the conjecture would
be false. Whether that configuration is realisable is an empirical question, and
:func:`conjecture2_counterexample_search` looks for it exhaustively.

If it holds, the search collapses from roughly ``3^m`` (every orientation state
of ``m`` undirected edges) to roughly ``2^k`` (every subset of ``k``
knowledge-oriented edges to retract).
"""

from __future__ import annotations

from collections import deque
from collections.abc import Callable
from dataclasses import dataclass

from bkrobust.core.conventions import UNREACHED
from bkrobust.core.spacelib import Space
from bkrobust.demo.graph import MPDAG


def check_no_empty_extensions(space: Space) -> bool:
    """Whether every element of ``space`` represents at least one DAG.

    The precondition that makes Conjecture 1 hold under this implementation's
    ``is_valid`` convention. Checked, not assumed.
    """
    return all(len(space.reps[g]) > 0 for g in space.elements)


def up_neighbours(space: Space, g: MPDAG) -> list[MPDAG]:
    """Neighbours strictly ABOVE ``g`` in model inclusion -- i.e. retractions.

    Moving up means representing more DAGs: less information, a retracted claim.
    Returned in a deterministic order.
    """
    reps_g = space.reps[g]
    out = [n for n in space.neighbours[g] if reps_g < space.reps[n]]
    return sorted(out, key=lambda h: h.edge_string())


def up_distances(space: Space, source: MPDAG) -> dict[MPDAG, int]:
    """BFS from ``source`` following only upward (retraction) covering edges.

    Args:
        space: The enumerated space.
        source: Usually ``G0``.

    Returns:
        Retraction-distance per element of the up-set of ``source``. Elements
        not above ``source`` are omitted.
    """
    dist = {source: 0}
    queue = deque([source])
    while queue:
        cur = queue.popleft()
        for nxt in up_neighbours(space, cur):
            if nxt not in dist:
                dist[nxt] = dist[cur] + 1
                queue.append(nxt)
    return dist


@dataclass(frozen=True)
class Conjecture1Result:
    """Outcome of checking upward-closure on one space and one ``Z``.

    Attributes:
        holds: Whether every comparable failing pair respected the implication.
        n_pairs_checked: Ordered pairs ``(G, G')`` with ``[G] <= [G']`` tested.
        n_failing: Elements where ``Z`` failed.
        violations: ``(G, G')`` pairs where ``G`` failed but ``G'`` did not.
        empty_extension_elements: Elements with ``[G]`` empty, which would break
            the implication under this implementation's convention.
    """

    holds: bool
    n_pairs_checked: int
    n_failing: int
    violations: tuple[tuple[str, str], ...]
    empty_extension_elements: tuple[str, ...]


def check_conjecture1(space: Space, fails: Callable[[MPDAG], bool]) -> Conjecture1Result:
    """Test upward-closure of failure over ALL comparable pairs of ``space``.

    Args:
        space: The enumerated space.
        fails: Predicate that is ``True`` where the property fails.

    Returns:
        A :class:`Conjecture1Result`.
    """
    failing = {g: fails(g) for g in space.elements}
    violations: list[tuple[str, str]] = []
    checked = 0
    for g in space.elements:
        for h in space.elements:
            if g is h:
                continue
            if space.reps[g] <= space.reps[h]:
                checked += 1
                if failing[g] and not failing[h]:
                    violations.append((g.edge_string(), h.edge_string()))
    empties = tuple(g.edge_string() for g in space.elements if not space.reps[g])
    return Conjecture1Result(
        holds=not violations,
        n_pairs_checked=checked,
        n_failing=sum(failing.values()),
        violations=tuple(violations),
        empty_extension_elements=empties,
    )


@dataclass(frozen=True)
class Conjecture2Result:
    """Outcome of comparing the full radius against the retraction-only radius.

    Attributes:
        r_full: Radius by BFS on the whole neighbour graph (the ground truth).
        r_up: Radius by BFS restricted to upward (retraction) moves.
        holds: Whether the two agree.
        witness_full: A nearest failure under the full search.
        witness_up: A nearest failure under the retraction-only search.
        n_space: Space size, for reporting enumeration scope.
    """

    r_full: int
    r_up: int
    holds: bool
    witness_full: str | None
    witness_up: str | None
    n_space: int


def check_conjecture2(
    space: Space,
    g0: MPDAG,
    fails: Callable[[MPDAG], bool],
    full_dists: dict[MPDAG, int],
) -> Conjecture2Result:
    """Compare the true radius with the minimum number of retractions to failure.

    Args:
        space: The enumerated space.
        g0: The analyst's graph.
        fails: Predicate that is ``True`` where the property fails.
        full_dists: BFS distances from ``g0`` on the full neighbour graph.

    Returns:
        A :class:`Conjecture2Result`. ``holds`` is ``True`` when the two radii
        agree, including when both are ``UNREACHED``.
    """
    best_full, wit_full = UNREACHED, None
    for g in space.elements:
        d = full_dists.get(g)
        if d is None or (best_full != UNREACHED and d >= best_full):
            continue
        if fails(g):
            best_full, wit_full = d, g

    ups = up_distances(space, g0)
    best_up, wit_up = UNREACHED, None
    for g in space.elements:
        d = ups.get(g)
        if d is None or (best_up != UNREACHED and d >= best_up):
            continue
        if fails(g):
            best_up, wit_up = d, g

    return Conjecture2Result(
        r_full=best_full,
        r_up=best_up,
        holds=(best_full == best_up),
        witness_full=None if wit_full is None else wit_full.edge_string(),
        witness_up=None if wit_up is None else wit_up.edge_string(),
        n_space=len(space.elements),
    )
