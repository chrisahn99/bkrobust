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
    ``n - 2`` remaining nodes as evenly as possible; requires ``n >= 6``).
    The last node of each chain directly confounds ``X`` and ``Y``:

        C1_last -> X,  C1_last -> Y,  C2_last -> X,  C2_last -> Y,  X -> Y

    Two structural facts, both load-bearing and verified directly in
    ``tests/synth/test_generators.py``:

    1. **The X-Y edge is always ambiguous.** Because ``{C1_last, C2_last}``
       are each adjacent to both ``X`` and ``Y``, neither pair is unshielded
       at ``Y``, so ``X -> Y`` is never marked by a v-structure and Meek's
       rules never force it either -- the CPDAG leaves it undirected. This is
       what keeps every instance from this family non-degenerate: the space
       always contains a reachable graph with ``Y -> X`` instead, at which
       point the optimal set (fixed from ``G0``) is no longer valid, so
       ``r_val`` and ``r_opt`` are finite whenever ``K_assumed`` orients
       ``X - Y``.
    2. **Coupling controls whether C1_last and C2_last are adjacent.** At
       ``coupling = 0.0`` they are not: {C1_last, C2_last} are then each
       other's only non-adjacent co-parent of ``X`` (and of ``Y``), so those
       four confounding edges get marked by a v-structure and are compelled
       -- frozen for the whole space. The two chains, having no edges to one
       another, end up as genuinely separate undirected (chordal) components
       of the CPDAG; ``C1_last``'s blocking role and ``C2_last``'s sit in
       different components (three components total: chain 1, chain 2, and
       the ``{X, Y}`` pair -- see point 1). Raising ``coupling`` adds cross
       edges, **all sourced from C1_last**, into a suffix of chain 2 taken
       tail-first: ``C1_last -> C2_last``, then ``C1_last -> C2_{g2-2}``,
       and so on. ``n_cross = round(coupling * g2)`` of these are added.
       Fanning every cross edge out of the *same single node* (rather than a
       complete bipartite join between the two chains) is deliberate and
       empirically load-bearing: a complete join creates fresh unshielded
       colliders deep in chain 2 whenever a chain has length >= 3 (two
       non-adjacent chain interior nodes both feeding the same chain-2 node)
       and those get compelled and cascade, silently undoing the very merge
       coupling is supposed to produce. The single-source fan-out never has
       that problem, because by the time a chain-2 node has C1_last as a
       parent, so does its own chain predecessor (added one step earlier in
       the tail-first order) -- so the two are always adjacent and no new
       v-structure fires. Empirically (see the generator test) this makes
       the merge monotonic in ``n_cross`` for every chain-length combination
       tried: 0 cross edges give three components, ``1..g2-1`` give two
       (``{X, Y}`` joins chain 1 through the critical edge; chain 2 stays
       separate), and ``g2`` (all of chain 2, i.e. ``coupling = 1.0``) give
       one -- the regime the prior worked example lives in.

    Args:
        n: Number of nodes; must be >= 6 (two confounder chains of length
            >= 2, plus X and Y).
        rng: Sole source of randomness (unused: the construction is
            deterministic given ``(n, coupling)``; kept in the signature for
            interface uniformity with the other generators and the
            :data:`GENERATORS` registry).
        coupling: In ``[0, 1]``. Fraction of chain 2, counting from its tail,
            that gets a direct edge from ``C1_last``; see above.

    Returns:
        A fully oriented, acyclic :class:`MPDAG` on ``n`` nodes.
    """
    if n < 6:
        raise ValueError(f"decoupled_backdoor_dag requires n >= 6, got {n}")
    if not (0.0 <= coupling <= 1.0):
        raise ValueError(f"coupling must be in [0, 1], got {coupling}")
    del rng  # reserved; the construction is otherwise deterministic given (n, coupling)

    remaining = n - 2
    g1 = -(-remaining // 2)  # ceil
    g2 = remaining - g1
    if g2 < 2:
        raise ValueError(f"n={n} does not split into two chains of length >= 2")

    c1 = [f"C1_{i}" for i in range(g1)]
    c2 = [f"C2_{j}" for j in range(g2)]
    nodes = ["X", "Y", *c1, *c2]

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

    n_cross = round(coupling * g2)
    for k in range(n_cross):
        j = g2 - 1 - k
        directed.append((c1[-1], c2[j]))

    return MPDAG(nodes, directed=directed)


GENERATORS: dict[str, Callable[..., MPDAG]] = {
    "erdos_renyi": erdos_renyi_dag,
    "scale_free": scale_free_dag,
    "block": block_dag,
    "decoupled_backdoor": decoupled_backdoor_dag,
}
