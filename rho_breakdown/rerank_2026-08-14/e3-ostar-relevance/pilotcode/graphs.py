"""
Graph machinery: DAG sampling, CPDAG, Meek R1-R4, MPDAG, background knowledge
consistency (Perkovic et al. UAI'17 Algorithm 1), brute-force MEC enumeration.

Representation: amat G (p x p int8).
  G[i,j]=1, G[j,i]=0  ->  i -> j   (directed)
  G[i,j]=1, G[j,i]=1  ->  i -- j   (undirected)
  G[i,j]=0, G[j,i]=0  ->  no edge
"""
from itertools import combinations, product

import numpy as np


# ---------------------------------------------------------------- basic queries
def is_directed(G, i, j):
    return G[i, j] == 1 and G[j, i] == 0


def is_undirected(G, i, j):
    return G[i, j] == 1 and G[j, i] == 1


def adjacent(G, i, j):
    return G[i, j] == 1 or G[j, i] == 1


def skeleton(G):
    return ((G + G.T) > 0).astype(np.int8)


def directed_edges(G):
    p = G.shape[0]
    return [(i, j) for i in range(p) for j in range(p) if is_directed(G, i, j)]


def undirected_edges(G):
    p = G.shape[0]
    return [(i, j) for i in range(p) for j in range(i + 1, p) if is_undirected(G, i, j)]


def has_directed_cycle(G):
    """Cycle among the *directed* edges only."""
    p = G.shape[0]
    A = np.zeros((p, p), dtype=np.int8)
    for i, j in directed_edges(G):
        A[i, j] = 1
    # DFS-based cycle detection
    color = np.zeros(p, dtype=np.int8)  # 0 white, 1 grey, 2 black

    def visit(u):
        color[u] = 1
        for v in np.flatnonzero(A[u]):
            if color[v] == 1:
                return True
            if color[v] == 0 and visit(v):
                return True
        color[u] = 2
        return False

    for u in range(p):
        if color[u] == 0 and visit(u):
            return True
    return False


def v_structures(G):
    """Unshielded colliders a -> b <- c with a,c non-adjacent. Only *directed* edges count."""
    p = G.shape[0]
    out = set()
    for b in range(p):
        pars = [a for a in range(p) if is_directed(G, a, b)]
        for a, c in combinations(sorted(pars), 2):
            if not adjacent(G, a, c):
                out.add((a, b, c))
    return out


# ---------------------------------------------------------------- DAG sampling
def random_dag(p, expected_degree, rng):
    """Erdos-Renyi over a random topological order. Returns amat (all directed)."""
    order = rng.permutation(p)
    prob = min(1.0, expected_degree / max(p - 1, 1))
    G = np.zeros((p, p), dtype=np.int8)
    for a in range(p):
        for b in range(a + 1, p):
            if rng.random() < prob:
                G[order[a], order[b]] = 1
    return G


def topological_order(G):
    """For a fully-directed graph."""
    p = G.shape[0]
    indeg = G.sum(axis=0).astype(int).copy()
    order, stack = [], [i for i in range(p) if indeg[i] == 0]
    stack.sort()
    while stack:
        u = stack.pop(0)
        order.append(u)
        for v in np.flatnonzero(G[u]):
            indeg[v] -= 1
            if indeg[v] == 0:
                stack.append(int(v))
        stack.sort()
    assert len(order) == p, "not a DAG"
    return order


# ---------------------------------------------------------------- Meek rules
def meek_closure(G):
    """Apply Meek's R1-R4 to fixpoint, in place on a copy. Returns new amat."""
    G = G.copy()
    p = G.shape[0]
    changed = True
    while changed:
        changed = False
        for a, b in [(a, b) for a in range(p) for b in range(p) if is_undirected(G, a, b)]:
            if not is_undirected(G, a, b):
                continue
            orient = False
            # R1: c -> a, a -- b, c not adj b  =>  a -> b
            for c in range(p):
                if c in (a, b):
                    continue
                if is_directed(G, c, a) and not adjacent(G, c, b):
                    orient = True
                    break
            # R2: a -> c -> b, a -- b  =>  a -> b
            if not orient:
                for c in range(p):
                    if c in (a, b):
                        continue
                    if is_directed(G, a, c) and is_directed(G, c, b):
                        orient = True
                        break
            # R3: a--c, a--d, c->b, d->b, c,d non-adjacent, a--b  =>  a -> b
            if not orient:
                cand = [c for c in range(p) if c not in (a, b)
                        and is_undirected(G, a, c) and is_directed(G, c, b)]
                for c, d in combinations(cand, 2):
                    if not adjacent(G, c, d):
                        orient = True
                        break
            # R4: a--d, d->c, c->b, a--c(or adjacent), b,d non-adjacent, a--b => a -> b
            if not orient:
                for d in range(p):
                    if d in (a, b) or not is_undirected(G, a, d) or adjacent(G, b, d):
                        continue
                    for c in range(p):
                        if c in (a, b, d):
                            continue
                        if is_directed(G, d, c) and is_directed(G, c, b) and adjacent(G, a, c):
                            orient = True
                            break
                    if orient:
                        break
            if orient:
                G[b, a] = 0
                changed = True
    return G


def dag_to_cpdag(D):
    """Skeleton + v-structures + Meek closure."""
    p = D.shape[0]
    G = skeleton(D)
    for (a, b, c) in v_structures(D):
        G[b, a] = 0
        G[b, c] = 0
    return meek_closure(G)


# ---------------------------------------------------------------- background knowledge
class MeekFail(Exception):
    pass


def apply_background_knowledge(C, K):
    """
    Perkovic et al. (UAI'17) Algorithm 1 / Meek (1995): orient each edge in K
    (list of (i,j) meaning i->j), closing under R1-R4 after each. Raise MeekFail
    if K is not consistent with C.

    Returns the MPDAG.
    """
    G = C.copy()
    for (i, j) in K:
        if is_directed(G, j, i):
            raise MeekFail(f"conflict: {i}->{j} but graph has {j}->{i}")
        if G[i, j] == 0 and G[j, i] == 0:
            raise MeekFail(f"{i}-{j} not an edge of the CPDAG")
        G[j, i] = 0
        G = meek_closure(G)
    if has_directed_cycle(G):
        raise MeekFail("directed cycle")
    if v_structures(G) != v_structures(C):
        raise MeekFail("new v-structure")
    return G


def is_consistent(C, K):
    try:
        apply_background_knowledge(C, K)
        return True
    except MeekFail:
        return False


# ---------------------------------------------------------------- brute force (ground truth for tests)
def consistent_dag_extensions(G, ref_vstructs=None, limit=None):
    """
    All DAGs obtained by orienting the undirected edges of G with no directed
    cycle and no new unshielded collider (relative to ref_vstructs, default =
    v-structures of G).
    """
    if ref_vstructs is None:
        ref_vstructs = v_structures(G)
    U = undirected_edges(G)
    out = []
    for bits in product([0, 1], repeat=len(U)):
        D = G.copy()
        for (i, j), b in zip(U, bits):
            if b:
                D[j, i] = 0
            else:
                D[i, j] = 0
        if has_directed_cycle(D):
            continue
        if v_structures(D) != ref_vstructs:
            continue
        out.append(D)
        if limit is not None and len(out) >= limit:
            break
    return out


def common_orientation_graph(dags, skel):
    """Edges oriented identically in every DAG -> directed; else undirected."""
    p = skel.shape[0]
    G = skel.copy()
    for i in range(p):
        for j in range(i + 1, p):
            if skel[i, j] == 0:
                continue
            fwd = all(D[i, j] == 1 for D in dags)
            bwd = all(D[j, i] == 1 for D in dags)
            if fwd:
                G[j, i] = 0
            elif bwd:
                G[i, j] = 0
    return G


def dag_agrees_with(D, G):
    """Is DAG D a consistent extension of (M)PDAG G? Same skeleton, all directed
    edges of G present in D with the same orientation, same v-structures."""
    if not np.array_equal(skeleton(D), skeleton(G)):
        return False
    for (i, j) in directed_edges(G):
        if not is_directed(D, i, j):
            return False
    return True
