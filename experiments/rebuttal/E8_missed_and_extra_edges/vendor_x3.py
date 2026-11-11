"""Vendored, verbatim, from the second author's working tree (X3 skeleton-perturbation study).

Source repository: the pilot tree
  commit 8486bb3509decf17adf490c078607b225d8df1b0
  (last commit touching rho_breakdown/x3_skeleton/code: 59066137b86022639c0eb5a4b068c1bb726f4a97)
Source files, in order: rho_breakdown/x3_skeleton/code/{graphs,adjust,scm,x1_ops}.py

Copied so that the compute machine needs only the submitted repository. Function bodies are
unchanged. The single edit is marked "E8 EDIT" (networkx < 3.3 has d_separated, not
is_d_separator). Representation: amat G[i,j]=1,G[j,i]=0 is i->j; both 1 is i--j.
"""
from itertools import combinations

import networkx as nx
import numpy as np

# E8 EDIT: networkx renamed d_separated -> is_d_separator in 3.3; keep both working.
_is_d_sep = getattr(nx, "is_d_separator", None) or nx.d_separated


# ==================== from graphs.py ====================

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

# ==================== from adjust.py ====================

def possibly_directed_neighbors(G, v):
    """w such that v -> w or v -- w."""
    p = G.shape[0]
    return [w for w in range(p) if w != v and G[v, w] == 1]

def poss_de(G, S):
    """Possible descendants of set S (including S): reachable by possibly-directed paths."""
    seen = set(S)
    stack = list(S)
    while stack:
        v = stack.pop()
        for w in possibly_directed_neighbors(G, v):
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return seen

PATH_CAP = 200_000

EXPAND_CAP = 5_000_000

class PathBlowup(Exception):
    pass

def possibly_causal_paths(G, x, y):
    """All simple possibly-causal paths x -> ... -> y (no edge v_{i+1} -> v_i).

    Raises PathBlowup past PATH_CAP complete paths or EXPAND_CAP expansions.
    """
    p = G.shape[0]
    paths = []
    nexp = [0]

    def dfs(path, visited):
        nexp[0] += 1
        if nexp[0] > EXPAND_CAP:
            raise PathBlowup("expansions")
        v = path[-1]
        if v == y:
            paths.append(list(path))
            if len(paths) > PATH_CAP:
                raise PathBlowup("paths")
            return
        for w in possibly_directed_neighbors(G, v):
            if w in visited:
                continue
            visited.add(w)
            path.append(w)
            dfs(path, visited)
            path.pop()
            visited.remove(w)

    dfs([x], {x})
    return paths

def causal_nodes(G, x, y):
    """cn(X,Y,G): nodes on proper possibly-causal paths from x to y, excluding x."""
    cn = set()
    for path in possibly_causal_paths(G, x, y):
        cn.update(path[1:])
    return cn

def forb(G, x, y):
    cn = causal_nodes(G, x, y)
    if not cn:
        return {x}
    return poss_de(G, cn) | {x}

def parents_of_set(G, S):
    p = G.shape[0]
    return {a for a in range(p) for b in S if is_directed(G, a, b)}

def is_amenable(G, x, y):
    """Every proper possibly-causal path from x to y starts with x -> v1."""
    paths = possibly_causal_paths(G, x, y)
    if not paths:
        return False  # no causal path; excluded upstream
    return all(is_directed(G, path[0], path[1]) for path in paths)

def optimal_adjustment_set(G, x, y):
    """O*(x,y,G). Returns None if G is not amenable rel. (x,y) (effect not
    identified by adjustment)."""
    if not is_amenable(G, x, y):
        return None
    cn = causal_nodes(G, x, y)
    O = parents_of_set(G, cn) - forb(G, x, y)
    return frozenset(O)

def _nx_from_dag(D):
    g = nx.DiGraph()
    g.add_nodes_from(range(D.shape[0]))
    for i in range(D.shape[0]):
        for j in range(D.shape[0]):
            if D[i, j] == 1:
                g.add_edge(i, j)
    return g

def is_valid_adjustment_set(D, x, y, Z):
    """Generalized adjustment criterion in the DAG D."""
    Z = set(Z)
    if Z & forb(D, x, y):
        return False
    cn = causal_nodes(D, x, y)
    Dp = D.copy()
    for w in cn:
        if is_directed(Dp, x, w):
            Dp[x, w] = 0  # proper back-door graph: cut first edges on causal paths
    g = _nx_from_dag(Dp)
    return _is_d_sep(  # E8 EDIT
        g, {x}, {y}, Z)

def total_effect_linear(A, x, y):
    """A[i,j] = structural coefficient of i in the equation for j.
    Total effect of x on y = ((I-A)^-1)[x,y]."""
    p = A.shape[0]
    return float(np.linalg.inv(np.eye(p) - A)[x, y])

def cov_linear(A, omega):
    """Sigma for x = A^T x + e, Var(e)=diag(omega)."""
    p = A.shape[0]
    M = np.linalg.inv(np.eye(p) - A.T)
    return M @ np.diag(omega) @ M.T

def ols_coefficient(Sigma, x, y, Z):
    """Population OLS coefficient of x in the regression of y on (x, Z)."""
    S = [x] + sorted(Z)
    Sxx = Sigma[np.ix_(S, S)]
    Sxy = Sigma[np.ix_(S, [y])]
    beta = np.linalg.solve(Sxx, Sxy)
    return float(beta[0, 0])

# ==================== from scm.py ====================

N_REF = 100_000

def _draw_weights(D, rng, lo=0.5, hi=2.0):
    W = np.zeros(D.shape)
    idx = np.argwhere(D == 1)
    signs = rng.choice([-1.0, 1.0], size=len(idx))
    mags = rng.uniform(lo, hi, size=len(idx))
    for (i, j), s, m in zip(idx, signs, mags):
        W[i, j] = s * m
    return W

def make_linear_iscm(D, rng, noise_scale=(0.5, 1.5)):
    """
    Returns dict with the FROZEN linear SEM: x = A^T x + c + e.
    A[i,j] = w_ij / sd_j ; e_j ~ N(0, (sigma_j/sd_j)^2).
    Intercepts are irrelevant for effects/covariances and are dropped.
    """
    p = D.shape[0]
    W = _draw_weights(D, rng)
    sigma = rng.uniform(*noise_scale, size=p)
    order = topological_order(D)

    # reference ancestral sample -> standardization constants (as in iSCM)
    Xref = np.zeros((N_REF, p))
    mu = np.zeros(p)
    sd = np.ones(p)
    for j in order:
        pa = np.flatnonzero(D[:, j])
        raw = Xref[:, pa] @ W[pa, j] + rng.normal(0, sigma[j], size=N_REF)
        mu[j] = raw.mean()
        sd[j] = raw.std()
        Xref[:, j] = (raw - mu[j]) / sd[j]

    A = W / sd[None, :]          # effective coefficients
    omega = (sigma / sd) ** 2    # effective noise variances
    return dict(kind="linear_iscm", D=D, W=W, sigma=sigma, mu=mu, sd=sd,
                A=A, omega=omega, order=order)

# ==================== from x1_ops.py ====================

UNREACH = 1 << 20

def hop_dist_from(Gsk, srcs, p):
    """BFS on a skeleton.  Sentinel UNREACH for unreachable vertices."""
    dist = np.full(p, UNREACH, dtype=int)
    frontier = list(srcs)
    for v in frontier:
        dist[v] = 0
    while frontier:
        nxt = []
        for v in frontier:
            for w in np.flatnonzero(Gsk[v]):
                if dist[w] > dist[v] + 1:
                    dist[w] = dist[v] + 1
                    nxt.append(int(w))
        frontier = nxt
    return dist

def stmt_hop(dist, a, b):
    return int(min(dist[a], dist[b]))
