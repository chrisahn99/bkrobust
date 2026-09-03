"""Random DAG generators for Axis A instances.

Every generator has signature ``(n, rng, **params) -> MPDAG`` and returns a
*fully oriented* DAG (no undirected edges). Node names are zero-padded ``V``
labels (``V0``, ``V01`` ...) chosen so that string-sorted order equals
numeric order, since :class:`~bkrobust.demo.graph.MPDAG` sorts nodes as
strings everywhere.

Determinism: every draw comes from the passed ``rng``
(:class:`numpy.random.Generator`); nothing touches the global NumPy or stdlib
RNG. Every place that could be sensitive to Python's per-process hash
randomisation (set/dict iteration order) is avoided by construction: loops
range over integer indices, not over hashed containers.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

from bkrobust.demo.graph import MPDAG

Edge = tuple[str, str]


def _node_names(n: int) -> list[str]:
    """``n`` zero-padded node labels, ``V0000 ... V{n-1}``, string-sort == numeric sort."""
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")
    width = max(1, len(str(n - 1)))
    return [f"V{i:0{width}d}" for i in range(n)]


def _random_topological_order(n: int, rng: np.random.Generator) -> np.ndarray:
    """A random permutation of ``0..n-1``.

    ``order[k]`` is the node-index placed at topological position ``k``.
    """
    return rng.permutation(n)


def erdos_renyi_dag(n: int, rng: np.random.Generator, edge_prob: float) -> MPDAG:
    """A random DAG: draw a random topological order, keep each forward pair independently.

    Args:
        n: Number of nodes.
        rng: Sole source of randomness.
        edge_prob: Probability, independent per ordered pair, that the earlier
            node in the drawn topological order gets an edge into the later one.

    Returns:
        A fully oriented, acyclic :class:`MPDAG` on ``n`` nodes.
    """
    if not (0.0 <= edge_prob <= 1.0):
        raise ValueError(f"edge_prob must be in [0, 1], got {edge_prob}")
    nodes = _node_names(n)
    order = _random_topological_order(n, rng)
    node_at = [nodes[i] for i in order]
    directed: list[Edge] = []
    for i in range(n):
        for j in range(i + 1, n):
            if rng.uniform() < edge_prob:
                directed.append((node_at[i], node_at[j]))
    return MPDAG(nodes, directed=directed)


def scale_free_dag(n: int, rng: np.random.Generator, m_attach: int) -> MPDAG:
    """A Barabasi-Albert preferential-attachment DAG, oriented along a random topological order.

    Builds the classical BA skeleton by hand (never delegating randomness to
    networkx, which does not accept a ``numpy.random.Generator`` uniformly
    across versions): the first new node attaches to all ``m_attach`` seed
    nodes; every subsequent node attaches to ``m_attach`` distinct existing
    nodes drawn with probability proportional to current degree, via the
    standard "repeated node list" trick. The skeleton is then oriented edge
    by edge according to a random topological order, giving a DAG whose
    skeleton is the BA graph and whose CPDAG typically leaves much of it
    undirected (chains and forks are Markov-equivalent to many orderings).

    Args:
        n: Number of nodes.
        rng: Sole source of randomness.
        m_attach: Edges each new node attaches with. Must satisfy
            ``1 <= m_attach < n``.

    Returns:
        A fully oriented, acyclic :class:`MPDAG` on ``n`` nodes with exactly
        ``(n - m_attach) * m_attach`` skeleton edges.
    """
    if not (1 <= m_attach < n):
        raise ValueError(f"m_attach must satisfy 1 <= m_attach < n, got {m_attach} (n={n})")
    nodes = _node_names(n)

    edges_idx: set[tuple[int, int]] = set()
    repeated: list[int] = []

    first_source = m_attach
    for t in range(m_attach):
        edges_idx.add((t, first_source))
        repeated.append(t)
        repeated.append(first_source)

    for source in range(m_attach + 1, n):
        chosen: set[int] = set()
        while len(chosen) < m_attach:
            idx = int(rng.integers(0, len(repeated)))
            chosen.add(repeated[idx])
        for t in sorted(chosen):
            edges_idx.add((t, source))
            repeated.append(t)
            repeated.append(source)

    order = _random_topological_order(n, rng)
    pos = [0] * n
    for k in range(n):
        pos[order[k]] = k

    directed: list[Edge] = []
    for a, b in sorted(edges_idx):
        if pos[a] < pos[b]:
            directed.append((nodes[a], nodes[b]))
        else:
            directed.append((nodes[b], nodes[a]))
    return MPDAG(nodes, directed=directed)


def block_dag(
    n: int,
    rng: np.random.Generator,
    n_blocks: int,
    p_within: float,
    p_between: float,
) -> MPDAG:
    """A DAG with community structure: nodes partitioned into blocks, denser within than between.

    A single random topological order governs the whole graph (so it is
    always acyclic regardless of block membership); within that order, each
    forward pair is connected with probability ``p_within`` if the two nodes
    share a block, else ``p_between``.

    Args:
        n: Number of nodes.
        rng: Sole source of randomness.
        n_blocks: Number of blocks. Nodes are assigned to blocks
            round-robin (``node i`` in block ``i % n_blocks``) so block sizes
            differ by at most one.
        p_within: Edge probability for a same-block forward pair.
        p_between: Edge probability for a cross-block forward pair.

    Returns:
        A fully oriented, acyclic :class:`MPDAG` on ``n`` nodes.
    """
    if n_blocks < 1:
        raise ValueError(f"n_blocks must be >= 1, got {n_blocks}")
    if not (0.0 <= p_within <= 1.0) or not (0.0 <= p_between <= 1.0):
        raise ValueError("p_within and p_between must be in [0, 1]")
    nodes = _node_names(n)
    block_of = [i % n_blocks for i in range(n)]
    order = _random_topological_order(n, rng)
    node_at = [nodes[i] for i in order]
    block_at = [block_of[i] for i in order]

    directed: list[Edge] = []
    for i in range(n):
        for j in range(i + 1, n):
            p = p_within if block_at[i] == block_at[j] else p_between
            if rng.uniform() < p:
                directed.append((node_at[i], node_at[j]))
    return MPDAG(nodes, directed=directed)


def decoupled_backdoor_dag(n: int, rng: np.random.Generator, coupling: float) -> MPDAG:
    r"""The key designed family: two back-door routes through disjoint confounder groups.

    Builds treatment ``X`` and outcome ``Y`` confounded through **two**
    separate chains ``C1_0 -> C1_1 -> ... -> C1_{g1-1}`` and
    ``C2_0 -> ... -> C2_{g2-1}`` (each of length >= 2, splitting the
    ``n - 3`` remaining nodes as evenly as possible; requires ``n >= 7`` --
    two chains of length >= 2, plus ``X``, ``Y``, and the ``W`` node below).
    The last node of each chain directly confounds ``X`` and ``Y``, and a
    dedicated node ``W`` is Y's fourth parent:

        C1_last -> X,  C1_last -> Y,  C2_last -> X,  C2_last -> Y,
        X -> Y,        W -> Y

    ``W`` has no other edge at all -- it exists purely to force ``X -> Y``.

    Revision history, because the first version of this docstring described a
    design that turned out to be broken and it matters why: the original
    construction omitted ``W`` and relied on ``X - Y`` staying *undirected* in
    the CPDAG for the family's non-degeneracy (a reachable ``Y -> X`` graph
    elsewhere in the space, invalidating the fixed optimal set). That is a
    real mechanism, but it has a fatal cost: with ``X - Y`` undirected,
    ``all_valid_adjustment_sets_mpdag(cpdag, "X", "Y")`` is **empty**, because
    some CPDAG extensions have ``Y -> X`` and others ``X -> Y``, and no single
    adjustment set is valid for both an "X causes Y" world and a "Y causes X"
    world at once. Concretely: this made the degeneracy gate reject every
    single instance from this generator with
    ``reason="no_valid_adjustment_set"`` (checked against the true DAG) or, at
    the CPDAG level relevant to identification at ``G0``, zero valid sets
    outright -- the family was unusable. Fixed here the same way the prior
    worked example avoids the analogous problem for its own treatment/outcome
    edge (``Age -> CVD`` there): give ``Y`` a parent, ``W``, that is adjacent
    to nothing else, so ``(X, W)`` is an unshielded pair at ``Y`` and
    ``X -> Y`` is compelled. Verified directly in
    ``tests/synth/test_generators.py``: ``("X", "Y")`` is always in
    ``cpdag.directed_edges``, and
    ``len(all_valid_adjustment_sets_mpdag(cpdag, "X", "Y")) >= 2`` at every
    tested coupling.

    A second consequence of forcing full identification at the bare-CPDAG
    level, found while fixing the first bug and worth stating plainly rather
    than glossing over: it also forces ``C1_last -> X``, ``C2_last -> X``,
    ``C1_last -> Y`` and ``C2_last -> Y`` to be compelled **unconditionally**,
    for the same reason (any residual ambiguity in a confounder's edge into
    ``X`` lets some extension make it a descendant of ``X`` instead of an
    ancestor, which breaks CPDAG-level identification exactly as the ``X - Y``
    ambiguity did). Since the Henckel-Perkovic-Maathuis optimal set is a pure
    function of that compelled backbone (``O = pa(Y) \\ {X}`` here, because
    ``Y`` is ``X``'s only causal descendant), ``O`` is therefore *structurally
    invariant across the entire perturbation space* for this family, no
    matter what background knowledge does to the chain interiors: every
    element of the space keeps the same ``O``, so ``r_val``, ``r_opt`` and
    ``r_eps`` all come out ``UNREACHED`` and every instance is rejected by the
    degeneracy gate's sanity check
    (``reason="no_atomic_perturbation_changes_validity"``). This is not a
    residual bug to chase further; it is the other side of the same coin as
    the fix above, confirmed by exhausting the plausible alternative
    placements of the remaining ambiguity (see the session notes) -- CPDAG-
    level identification and a perturbable ``(X, Y)`` radius are in tension
    for any "two disjoint DIRECT confounder groups" construction of this
    shape. What the family *does* still deliver, and what the tests below
    check, is the CPDAG-structural claim it was built for: two backdoor
    routes whose blocking nodes provably sit in different chordal components
    at ``coupling=0`` and the same component at ``coupling=1``, with each
    chain's interior still genuinely perturbable (background knowledge can
    still be wrong about it) even though that perturbability does not reach
    ``(X, Y)``'s own radius.

    Coupling mechanism, unchanged in spirit from before but now sourced away
    from the compelled backbone: cross edges run from ``C1_{g1-2}`` (chain
    1's second-to-last node -- never ``C1_last`` itself, so ``C1_last``'s
    edges to ``X``/``Y`` are untouched) into a tail-first suffix of chain 2,
    ``round(coupling * g2)`` of them. At ``coupling=0`` the two chains share
    no edge and end up in different undirected components; at ``coupling=1``
    every edge is added and ``C1_{g1-2}`` (hence, transitively via the chain
    edge, ``C1_last``) becomes adjacent to every node of chain 2 (hence
    ``C2_last``), merging the two chains' components into one. Sourcing every
    cross edge from the single fixed node ``C1_{g1-2}`` (rather than a
    complete bipartite join) avoids the same fresh-collider cascade described
    in the previous revision of this docstring: whenever a chain-2 node
    ``C2_j`` gets the cross edge, so does its chain predecessor ``C2_{j-1}``
    (added one step earlier in the tail-first order), so the two are always
    adjacent and no new v-structure fires at ``C2_j``.

    Args:
        n: Number of nodes; must be >= 7 (two confounder chains of length
            >= 2, plus X, Y and W).
        rng: Sole source of randomness (unused: the construction is
            deterministic given ``(n, coupling)``; kept in the signature for
            interface uniformity with the other generators and the
            :data:`GENERATORS` registry).
        coupling: In ``[0, 1]``. Fraction of chain 2, counting from its tail,
            that gets a direct edge from ``C1_{g1-2}``; see above.

    Returns:
        A fully oriented, acyclic :class:`MPDAG` on ``n`` nodes.
    """
    if n < 7:
        raise ValueError(f"decoupled_backdoor_dag requires n >= 7, got {n}")
    if not (0.0 <= coupling <= 1.0):
        raise ValueError(f"coupling must be in [0, 1], got {coupling}")
    del rng  # reserved; the construction is otherwise deterministic given (n, coupling)

    remaining = n - 3
    g1 = -(-remaining // 2)  # ceil
    g2 = remaining - g1
    if g2 < 2:
        raise ValueError(f"n={n} does not split into two chains of length >= 2")

    c1 = [f"C1_{i}" for i in range(g1)]
    c2 = [f"C2_{j}" for j in range(g2)]
    nodes = ["X", "Y", "W", *c1, *c2]

    directed: list[Edge] = []
    for i in range(g1 - 1):
        directed.append((c1[i], c1[i + 1]))
    for j in range(g2 - 1):
        directed.append((c2[j], c2[j + 1]))
    directed.append((c1[-1], "X"))
    directed.append((c1[-1], "Y"))
    directed.append((c2[-1], "X"))
    directed.append((c2[-1], "Y"))
    directed.append(("X", "Y"))
    directed.append(("W", "Y"))

    cross_source = c1[g1 - 2]
    n_cross = round(coupling * g2)
    for k in range(n_cross):
        j = g2 - 1 - k
        directed.append((cross_source, c2[j]))

    return MPDAG(nodes, directed=directed)


GENERATORS: dict[str, Callable[..., MPDAG]] = {
    "erdos_renyi": erdos_renyi_dag,
    "scale_free": scale_free_dag,
    "block": block_dag,
    "decoupled_backdoor": decoupled_backdoor_dag,
}
