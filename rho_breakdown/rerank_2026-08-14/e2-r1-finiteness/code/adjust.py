r"""
Adjustment-set machinery for MPDAGs.

O*(X,Y,G) = pa(cn(X,Y,G), G) \ forb(X,Y,G)
  Henckel, Perkovic & Maathuis (arXiv:1907.02435), building on Perkovic et al.'s
  complete graphical characterisation of adjustment in MPDAGs.

Validity in a DAG: generalized adjustment criterion (Perkovic et al. 2018):
  Z valid for (X,Y) iff Z n forb(X,Y,D) = {} and Z d-separates X,Y in the
  proper back-door graph.

Single-node X, single-node Y throughout (that is all this pilot needs).
"""
import networkx as nx
import numpy as np

from graphs import adjacent, is_directed, is_undirected


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


def possibly_causal_paths(G, x, y):
    """All simple possibly-causal paths x -> ... -> y (no edge v_{i+1} -> v_i)."""
    p = G.shape[0]
    paths = []

    def dfs(path, visited):
        v = path[-1]
        if v == y:
            paths.append(list(path))
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


# ----------------------------------------------------------- validity in the true DAG
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
    return nx.is_d_separator(g, {x}, {y}, Z)


# ----------------------------------------------------------- linear-SEM estimands
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
