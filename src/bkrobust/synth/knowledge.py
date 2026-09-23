"""Simulating an analyst's background knowledge, rather than hand-writing it.

``draw_k_true`` plays the role of an honest analyst: a random subset of the
CPDAG's undirected edges, oriented exactly as the (unknown to the analyst,
known to us) true DAG has them. The corruption operators then play the role
of an analyst who is imperfect in one of three ways:

* ``omit``: knows less than they claim -- some true claims are simply dropped.
* ``flip``: is confidently wrong about some claims -- their direction is
  reversed.
* ``compound``: both at once (``omit`` then ``flip`` on the remainder).
* ``tiered``: the realistic case. Rather than a bag of independent edge
  claims, the analyst asserts a *temporal ordering* over variables (baseline
  covariates before treatment before outcome, say) and every CPDAG edge that
  crosses tiers gets oriented at once, in the tier direction. A corrupted
  tier assignment misfiles some nodes into the wrong tier before those
  orientations are read off, so a single misplaced node can flip several
  assertions simultaneously -- this is what makes tiered knowledge behave so
  differently from the independent-edge corruptions above.

Call signatures. ``omit``, ``flip`` and ``compound`` share
``(k, rng, rate) -> list[Edge]``: they corrupt an existing knowledge list and
need nothing else. ``tiered`` is generative rather than corruptive -- it does
not take an existing ``k`` at all, since tiered knowledge is derived fresh
from ``(dag, cpdag)`` each time -- so its signature is
``(dag, cpdag, rng, n_tiers, corruption_rate) -> list[Edge]``. :data:`CORRUPTIONS`
holds both kinds under one registry; callers (the runner) dispatch on the
corruption name to know which call shape to use.

Determinism: every draw comes from the passed ``rng``. Iteration is always
over ``sorted(...)`` wherever the order could affect an RNG draw or a
downstream accumulation, since raw ``set``/``dict`` iteration order depends on
``PYTHONHASHSEED``.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Callable

import numpy as np

from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations

Edge = tuple[str, str]


def _check_rate(rate: float, name: str = "rate") -> None:
    if not (0.0 <= rate <= 1.0):
        raise ValueError(f"{name} must be in [0, 1], got {rate}")


# --- the analyst's true beliefs ----------------------------------------------


def draw_k_true(
    dag: MPDAG,
    cpdag: MPDAG,
    rng: np.random.Generator,
    knows_fraction: float,
) -> list[Edge]:
    """A random subset of ``cpdag``'s undirected edges, oriented as ``dag`` truly has them.

    Args:
        dag: The ground-truth DAG. Must be a consistent extension of
            ``cpdag`` (every one of ``cpdag``'s undirected edges must be
            oriented, one way or the other, in ``dag``).
        cpdag: The CPDAG background knowledge will be imposed on.
        rng: Sole source of randomness.
        knows_fraction: Fraction (``[0, 1]``) of ``cpdag``'s undirected edges
            the analyst claims to know. The count is
            ``round(knows_fraction * n_undirected)``, drawn without
            replacement.

    Returns:
        ``(tail, head)`` orientation claims, all true of ``dag``, sorted.

    Raises:
        ValueError: If ``knows_fraction`` is outside ``[0, 1]``, or if some
            selected edge is undirected in ``cpdag`` but not oriented in
            ``dag`` (i.e. ``dag`` does not actually extend ``cpdag``).
    """
    _check_rate(knows_fraction, "knows_fraction")
    undirected = sorted(cpdag.undirected_edges)
    n = len(undirected)
    if n == 0:
        return []
    k_count = round(knows_fraction * n)
    chosen = sorted(int(i) for i in rng.choice(n, size=k_count, replace=False))

    asserted: list[Edge] = []
    for i in chosen:
        a, b = undirected[i]
        if dag.is_directed_edge(a, b):
            asserted.append((a, b))
        elif dag.is_directed_edge(b, a):
            asserted.append((b, a))
        else:
            raise ValueError(
                f"dag does not orient cpdag's undirected edge {(a, b)}; "
                "dag is not a consistent extension of cpdag"
            )
    return sorted(asserted)


# --- corruption operators (act on an existing knowledge list) ---------------


def omit(k: list[Edge], rng: np.random.Generator, rate: float) -> list[Edge]:
    """Drop a ``rate`` fraction of ``k``'s claims, at random, without replacement.

    Args:
        k: The knowledge list to corrupt.
        rng: Sole source of randomness.
        rate: Fraction (``[0, 1]``) of claims to drop.

    Returns:
        A new list with ``round(rate * len(k))`` fewer entries.
    """
    _check_rate(rate)
    k_sorted = sorted(k)
    n = len(k_sorted)
    if n == 0:
        return []
    n_drop = round(rate * n)
    drop = set(int(i) for i in rng.choice(n, size=n_drop, replace=False))
    return [e for idx, e in enumerate(k_sorted) if idx not in drop]


def flip(k: list[Edge], rng: np.random.Generator, rate: float) -> list[Edge]:
    """Reverse a ``rate`` fraction of ``k``'s claims, at random, without replacement.

    Args:
        k: The knowledge list to corrupt.
        rng: Sole source of randomness.
        rate: Fraction (``[0, 1]``) of claims to reverse.

    Returns:
        A new list, same length as ``k``, with ``round(rate * len(k))``
        entries swapped ``(tail, head) -> (head, tail)`` and the rest
        unchanged.
    """
    _check_rate(rate)
    k_sorted = sorted(k)
    n = len(k_sorted)
    if n == 0:
        return []
    n_flip = round(rate * n)
    flip_idx = set(int(i) for i in rng.choice(n, size=n_flip, replace=False))
    return [(b, a) if idx in flip_idx else (a, b) for idx, (a, b) in enumerate(k_sorted)]


def compound(k: list[Edge], rng: np.random.Generator, rate: float) -> list[Edge]:
    """Both corruptions at once: :func:`omit` at ``rate``, then :func:`flip` the remainder.

    Args:
        k: The knowledge list to corrupt.
        rng: Sole source of randomness, threaded through both stages in order.
        rate: Fraction (``[0, 1]``) applied at each stage.

    Returns:
        A new list, shorter than ``k`` (from the omission) and with some of
        the survivors reversed (from the flip).
    """
    _check_rate(rate)
    return flip(omit(k, rng, rate), rng, rate)


# --- tiered / temporal knowledge ---------------------------------------------


def _topological_depths(dag: MPDAG) -> dict[str, int]:
    """Longest-path depth of every node (0 for a source), via Kahn's algorithm.

    Depth is well-defined independent of processing order within a layer: a
    node is only relaxed by :func:`dag.children` after every one of its own
    parents has already been processed (guaranteed by the indegree-zero
    admission condition), so the final depth does not depend on the order
    same-layer nodes are popped in -- only on the DAG's structure.

    Raises:
        ValueError: If ``dag`` is not acyclic (should not happen for a DAG
            produced by this package's generators, but checked rather than
            assumed).
    """
    nodes = sorted(dag.nodes)
    remaining_indegree = {v: len(dag.parents(v)) for v in nodes}
    depth = {v: 0 for v in nodes}
    frontier: deque[str] = deque(sorted(v for v in nodes if remaining_indegree[v] == 0))
    processed = 0
    while frontier:
        v = frontier.popleft()
        processed += 1
        for c in sorted(dag.children(v)):
            if depth[v] + 1 > depth[c]:
                depth[c] = depth[v] + 1
            remaining_indegree[c] -= 1
            if remaining_indegree[c] == 0:
                frontier.append(c)
    if processed != len(nodes):
        raise ValueError("dag is not acyclic: topological sort could not process every node")
    return depth


def tiered(
    dag: MPDAG,
    cpdag: MPDAG,
    rng: np.random.Generator,
    n_tiers: int,
    corruption_rate: float = 0.0,
) -> list[Edge]:
    """Tiered/temporal knowledge: partition nodes into ordered tiers, assert every crossing edge.

    Nodes are ordered by longest-path depth in ``dag`` (ties broken by node
    name) and split into ``n_tiers`` contiguous, near-equal groups in that
    order -- so, uncorrupted, every tier boundary respects ``dag``'s true
    causal order and the resulting assertions are all true. With
    ``corruption_rate > 0``, that many nodes (as a fraction of the total) are
    then reassigned to a uniformly random *different* tier before reading off
    assertions -- modelling an analyst who gets a variable's temporal position
    wrong, which can silently reverse or introduce several claims at once
    (every cross-tier edge touching that node), not just one.

    Every CPDAG undirected edge whose endpoints land in different tiers is
    asserted, oriented from the earlier tier to the later one; same-tier
    edges are left alone (the tiering says nothing about their order).

    Args:
        dag: The ground-truth DAG, used only to compute the depth-based tier
            order (its edges are not asserted directly).
        cpdag: The CPDAG whose undirected edges are the candidates for
            assertion.
        rng: Sole source of randomness (used only when ``corruption_rate > 0``).
        n_tiers: Number of tiers, >= 1.
        corruption_rate: Fraction (``[0, 1]``) of nodes misplaced into a
            random other tier before assertions are read off. ``0.0`` (the
            default) gives the true, uncorrupted tiering.

    Returns:
        ``(tail, head)`` orientation claims for every cross-tier CPDAG edge,
        sorted.
    """
    if n_tiers < 1:
        raise ValueError(f"n_tiers must be >= 1, got {n_tiers}")
    _check_rate(corruption_rate, "corruption_rate")

    depth = _topological_depths(dag)
    nodes_by_depth = sorted(dag.nodes, key=lambda v: (depth[v], v))
    n = len(nodes_by_depth)
    tier_of = {v: (idx * n_tiers) // n for idx, v in enumerate(nodes_by_depth)} if n else {}

    if corruption_rate > 0.0 and n > 0:
        nodes_sorted = sorted(dag.nodes)
        n_corrupt = round(corruption_rate * n)
        corrupt_positions = sorted(int(i) for i in rng.choice(n, size=n_corrupt, replace=False))
        for pos in corrupt_positions:
            v = nodes_sorted[pos]
            candidates = [t for t in range(n_tiers) if t != tier_of[v]]
            if candidates:
                tier_of[v] = candidates[int(rng.integers(0, len(candidates)))]

    asserted: list[Edge] = []
    for a, b in sorted(cpdag.undirected_edges):
        if tier_of[a] < tier_of[b]:
            asserted.append((a, b))
        elif tier_of[b] < tier_of[a]:
            asserted.append((b, a))
    return sorted(asserted)


# --- consistency check --------------------------------------------------


def check_consistent(cpdag: MPDAG, k: list[Edge]) -> bool:
    """Whether ``k`` is consistent with ``cpdag``: ``apply_orientations`` does not FAIL.

    Args:
        cpdag: The CPDAG to impose ``k`` on.
        k: The orientation claims.

    Returns:
        True iff :func:`bkrobust.demo.meek.apply_orientations` returns a
        graph rather than ``None``.
    """
    return apply_orientations(cpdag, k) is not None


CORRUPTIONS: dict[str, Callable[..., list[Edge]]] = {
    "omit": omit,
    "flip": flip,
    "compound": compound,
    "tiered": tiered,
}
