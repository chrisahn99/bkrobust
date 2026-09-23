"""Space construction, the model-inclusion order, and BFS distance. Frozen.

The space is

    G_C = { G a valid MPDAG : [G] subset-or-equal [C] }

ordered by **model inclusion**. The inclusion condition is load-bearing and was
a real bug once: filtering by self-validity alone admits graphs from a
*different* Markov equivalence class, because a fully oriented DAG is vacuously
Meek-closed and extends itself. Do not reintroduce that.

The covering relation is computed from model inclusion by enumerating
represented DAGs -- never from a conjectured characterisation -- and the
neighbour graph's edges are exactly the covering pairs. Distance is BFS hop
count on that graph.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Callable
from dataclasses import dataclass

from bkrobust.core.conventions import UNREACHED
from bkrobust.core.oracle import extensions
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.space import covering_pairs, enumerate_space, neighbour_graph


@dataclass(frozen=True)
class Space:
    """A fully enumerated perturbation space with its order and neighbour graph.

    Attributes:
        cpdag: The CPDAG the space refines.
        elements: Every element, in the enumerator's deterministic order.
        reps: ``[G]`` for each element.
        covers: Covering pairs as ``(lower, upper)`` with ``lower < upper``.
        neighbours: Symmetric adjacency induced by ``covers``.
    """

    cpdag: MPDAG
    elements: tuple[MPDAG, ...]
    reps: dict[MPDAG, frozenset[MPDAG]]
    covers: frozenset[tuple[MPDAG, MPDAG]]
    neighbours: dict[MPDAG, set[MPDAG]]

    def __len__(self) -> int:
        """Number of elements."""
        return len(self.elements)


def build_space(cpdag: MPDAG) -> Space:
    """Enumerate the space, its represented-DAG sets, covers and neighbour graph.

    Args:
        cpdag: The CPDAG to refine.

    Returns:
        A :class:`Space`.
    """
    elements = tuple(enumerate_space(cpdag))
    reps = {g: frozenset(extensions(g)) for g in elements}
    covers = frozenset(covering_pairs(list(elements), reps))
    nbrs = neighbour_graph(list(elements), covers)
    return Space(cpdag=cpdag, elements=elements, reps=reps, covers=covers, neighbours=nbrs)


def leq(space: Space, g1: MPDAG, g2: MPDAG) -> bool:
    """Model inclusion: ``[g1] subset-or-equal [g2]``, i.e. ``g1`` is at least as informative."""
    return space.reps[g1] <= space.reps[g2]


def distances_from(space: Space, source: MPDAG) -> dict[MPDAG, int]:
    """BFS hop count on the neighbour graph.

    Args:
        space: The enumerated space.
        source: Usually ``G0``.

    Returns:
        Distance per reachable element. Unreachable elements are **omitted**,
        never recorded as infinity or as a large number.

    Raises:
        KeyError: If ``source`` is not an element of ``space``.
    """
    if source not in space.neighbours:
        raise KeyError("source is not an element of the space")
    dist = {source: 0}
    queue = deque([source])
    while queue:
        cur = queue.popleft()
        for nxt in sorted(space.neighbours[cur], key=lambda g: g.edge_string()):
            if nxt not in dist:
                dist[nxt] = dist[cur] + 1
                queue.append(nxt)
    return dist


def radius(
    space: Space,
    dists: dict[MPDAG, int],
    fails: Callable[[MPDAG], bool],
) -> tuple[int, MPDAG | None]:
    """The normative radius: distance to the nearest element where ``fails`` holds.

    See :mod:`bkrobust.core.conventions`. Shells ``0 .. r-1`` are certified
    clean; shell ``r`` contains the returned witness.

    Args:
        space: The enumerated space.
        dists: Distances from ``G0``, as returned by :func:`distances_from`.
        fails: Predicate that is ``True`` where the property FAILS.

    Returns:
        ``(r, witness)``. ``r`` is :data:`~bkrobust.core.conventions.UNREACHED`
        with witness ``None`` when nothing reachable fails. ``r == 0`` means the
        property already fails at the source, which is degenerate.
    """
    best, witness = UNREACHED, None
    for g in space.elements:
        d = dists.get(g)
        if d is None:
            continue
        if best != UNREACHED and d >= best:
            continue
        if fails(g):
            best, witness = d, g
    return best, witness


def shell_sizes(dists: dict[MPDAG, int]) -> dict[int, int]:
    """Number of elements at each shell index."""
    out: dict[int, int] = {}
    for d in dists.values():
        out[d] = out.get(d, 0) + 1
    return dict(sorted(out.items()))
