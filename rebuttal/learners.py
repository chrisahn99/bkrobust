"""E30-E34: CPDAG learners on Gaussian data, all returning an MPDAG.

In-house (numpy/scipy, extending rebuttal/discovery.py, which is not changed):
  pc(alpha)          PC-stable with the sepset collider rule; with rule="sepset"
                     it reproduces discovery.pc exactly (checked in learners_check.py).
  pc(rule="cons")    conservative collider rule (Ramsey, Zhang & Spirtes 2006):
                     a-b-c unshielded is a collider iff b is in no separating
                     set of a, c found among subsets of adj(a) or adj(c), a
                     non-collider iff in all, ambiguous otherwise.
  pc(rule="maj")     majority rule (Colombo & Maathuis 2014): collider iff b is
                     in fewer than half of those sets, ambiguous at exactly half.
                     Ambiguous triples are left unoriented. Unlike Tetrad's CPC,
                     Meek's rules may still propagate through them, so that the
                     output is a closed PDAG the radius code accepts.
  hc()               greedy hill climbing over DAGs (add / delete / reverse,
                     best improvement, from the empty graph) with the Gaussian
                     BIC (lambda = 0.5), then the CPDAG of the learned DAG.
causal-learn 0.1.4.8 (installed with the user's approval, 11 November 2026):
  ges, boss, grasp   all with local_score_BIC_from_cov, lambda = 0.5 (standard BIC).

Every output is converted to an MPDAG over the network's node names.
"""

from __future__ import annotations

import warnings
from itertools import combinations
from math import log, sqrt

import numpy as np

import discovery as D

from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import meek_closure


# ---------------------------------------------------------------------------
# PC family


def _skeleton(t, p, alpha, max_depth):
    """Exactly the adjacency search of discovery.pc."""
    adj = {i: set(range(p)) - {i} for i in range(p)}
    sep, minp = {}, {}
    depth = 0
    while depth <= max_depth:
        snapshot = {i: set(a) for i, a in adj.items()}
        any_test = False
        for i in range(p):
            for j in sorted(snapshot[i]):
                if j not in adj[i]:
                    continue
                cand = sorted(snapshot[i] - {j})
                if len(cand) < depth:
                    continue
                any_test = True
                for S in combinations(cand, depth):
                    pv = t.pval(i, j, list(S))
                    if pv > alpha:
                        adj[i].discard(j)
                        adj[j].discard(i)
                        sep[frozenset((i, j))] = set(S)
                        minp[frozenset((i, j))] = pv
                        break
        if not any_test:
            break
        depth += 1
    return adj, sep, minp


def _all_sepsets(t, adj, a, c, alpha, max_depth):
    out = []
    for base in (adj[a] - {c}, adj[c] - {a}):
        base = sorted(base)
        for k in range(0, min(max_depth, len(base)) + 1):
            for S in combinations(base, k):
                if t.pval(a, c, list(S)) > alpha:
                    out.append(frozenset(S))
    return set(out)


def pc(nodes, X, alpha=0.01, max_depth=3, rule="sepset"):
    p = len(nodes)
    if max_depth is None:
        max_depth = p
    t = D.FisherZ(X)
    adj, sep, minp = _skeleton(t, p, alpha, max_depth)
    vs, n_amb = [], 0
    for b in range(p):
        for a, c in combinations(sorted(adj[b]), 2):
            if c in adj[a]:
                continue
            if rule == "sepset":
                coll = b not in sep.get(frozenset((a, c)), set())
            else:
                S = _all_sepsets(t, adj, a, c, alpha, max_depth)
                if not S:  # fall back on the skeleton's sepset
                    coll = b not in sep.get(frozenset((a, c)), set())
                else:
                    share = sum(b in s for s in S) / len(S)
                    if rule == "cons":
                        coll = share == 0.0
                        amb = 0.0 < share < 1.0
                    elif rule == "maj":
                        coll = share < 0.5
                        amb = share == 0.5
                    else:
                        raise ValueError(rule)
                    n_amb += amb
            if coll:
                vs.append((minp.get(frozenset((a, c)), 1.0), a, b, c))
    vs.sort(reverse=True)
    g, info = _orient(nodes, adj, vs)
    info["n_ambiguous"] = n_amb
    return g, info


def _orient(nodes, adj, vs):
    """discovery.pc's v-structure and closure stage, unchanged."""
    p = len(nodes)

    def pdag(vlist):
        head = {}
        for _, a, b, c in vlist:
            for u in (a, c):
                head.setdefault(frozenset((u, b)), set()).add(b)
        directed, undirected = [], []
        for i in range(p):
            for j in adj[i]:
                if i < j:
                    hs = head.get(frozenset((i, j)), set())
                    if hs == {j}:
                        directed.append((nodes[i], nodes[j]))
                    elif hs == {i}:
                        directed.append((nodes[j], nodes[i]))
                    else:
                        undirected.append((nodes[i], nodes[j]))
        g = MPDAG(nodes, directed=directed, undirected=undirected)
        return g if g.is_acyclic() else None

    used, dropped, heads = [], 0, set()
    for v in vs:
        _, a, b, c = v
        if (b, a) in heads or (b, c) in heads or pdag(used + [v]) is None:
            dropped += 1
            continue
        used.append(v)
        heads |= {(a, b), (c, b)}
    g = None
    while g is None:
        base = pdag(used)
        g = meek_closure(base) if base is not None else None
        if g is None:
            if not used:
                break
            used.pop()
            dropped += 1
    if g is None:
        g = MPDAG(nodes, undirected=[(nodes[i], nodes[j]) for i in range(p) for j in adj[i] if i < j])
    return g, {"n_vstructures": len(vs), "vstructures_dropped": dropped,
               "n_edges": sum(len(a) for a in adj.values()) // 2}


# ---------------------------------------------------------------------------
# hill climbing (DAG search) with the Gaussian BIC


def _local_bic(C, n, j, pa, lam=0.5):
    s = C[j, j]
    if pa:
        pa = list(pa)
        b = np.linalg.solve(C[np.ix_(pa, pa)], C[pa, j])
        s = s - C[j, pa] @ b
    s = max(s, 1e-12)
    return -0.5 * n * (1 + log(s)) - lam * (len(pa) + 1) * log(n)


def _reach(children, src, dst, skip=None):
    """Is there a directed path src ~> dst, optionally ignoring edge skip?"""
    stack, seen = [src], {src}
    while stack:
        v = stack.pop()
        for c in children[v]:
            if skip is not None and (v, c) == skip:
                continue
            if c == dst:
                return True
            if c not in seen:
                seen.add(c)
                stack.append(c)
    return False


def hc(nodes, X, lam=0.5, max_iter=10000):
    n, p = X.shape
    C = np.cov(X, rowvar=False)
    pa = [set() for _ in range(p)]
    ch = [set() for _ in range(p)]
    cur = [_local_bic(C, n, j, (), lam) for j in range(p)]
    add = np.full((p, p), -np.inf)   # add i -> j
    dele = np.full((p, p), -np.inf)  # delete i -> j

    def refresh(j):
        add[:, j] = -np.inf
        dele[:, j] = -np.inf
        for i in range(p):
            if i == j:
                continue
            if i in pa[j]:
                dele[i, j] = _local_bic(C, n, j, sorted(pa[j] - {i}), lam) - cur[j]
            elif j not in pa[i]:
                add[i, j] = _local_bic(C, n, j, sorted(pa[j] | {i}), lam) - cur[j]

    for j in range(p):
        refresh(j)
    it = 0
    while it < max_iter:
        it += 1
        # reversal i -> j  ==  delete at j, add j -> i at i
        rev = np.full((p, p), -np.inf)
        for j in range(p):
            for i in pa[j]:
                rev[i, j] = dele[i, j] + (_local_bic(C, n, i, sorted(pa[i] | {j}), lam) - cur[i])
        best = None
        cands = []
        for kind, M in (("add", add), ("del", dele), ("rev", rev)):
            flat = np.argsort(M, axis=None)[::-1]
            for f in flat[:200]:
                i, j = divmod(int(f), p)
                if M[i, j] <= 1e-9:
                    break
                cands.append((M[i, j], kind, i, j))
        cands.sort(reverse=True)
        for gain, kind, i, j in cands:
            if kind == "add" and _reach(ch, j, i):
                continue
            if kind == "rev" and _reach(ch, i, j, skip=(i, j)):
                continue
            best = (kind, i, j)
            break
        if best is None:
            break
        kind, i, j = best
        if kind == "add":
            pa[j].add(i); ch[i].add(j)
            cur[j] = _local_bic(C, n, j, sorted(pa[j]), lam)
            refresh(j)
        elif kind == "del":
            pa[j].discard(i); ch[i].discard(j)
            cur[j] = _local_bic(C, n, j, sorted(pa[j]), lam)
            refresh(j)
        else:
            pa[j].discard(i); ch[i].discard(j)
            pa[i].add(j); ch[j].add(i)
            cur[j] = _local_bic(C, n, j, sorted(pa[j]), lam)
            cur[i] = _local_bic(C, n, i, sorted(pa[i]), lam)
            refresh(j)
            refresh(i)
        # moves into other nodes whose legality changed are re-checked lazily
    dag = MPDAG(nodes, directed=[(nodes[i], nodes[j]) for j in range(p) for i in pa[j]])
    return dag_to_cpdag(dag), {"n_edges": sum(len(s) for s in pa), "iterations": it}


# ---------------------------------------------------------------------------
# causal-learn


def _from_general_graph(G, nodes):
    g = G.graph
    p = len(nodes)
    directed, undirected = [], []
    for i in range(p):
        for j in range(i + 1, p):
            a, b = g[i, j], g[j, i]
            if a == 0 and b == 0:
                continue
            if a == -1 and b == 1:
                directed.append((nodes[i], nodes[j]))
            elif a == 1 and b == -1:
                directed.append((nodes[j], nodes[i]))
            elif a == -1 and b == -1:
                undirected.append((nodes[i], nodes[j]))
            else:
                raise ValueError(f"unexpected endpoints {a},{b} on {nodes[i]},{nodes[j]}")
    return MPDAG(nodes, directed=directed, undirected=undirected)


def causal_learn(name, nodes, X, lam=0.5):
    warnings.filterwarnings("ignore")
    params = {"lambda_value": lam}
    if name == "ges":
        from causallearn.search.ScoreBased.GES import ges
        G = ges(X, score_func="local_score_BIC_from_cov", parameters=params)["G"]
    elif name == "boss":
        from causallearn.search.PermutationBased.BOSS import boss
        G = boss(X, score_func="local_score_BIC_from_cov", parameters=params, verbose=False)
    elif name == "grasp":
        from causallearn.search.PermutationBased.GRaSP import grasp
        G = grasp(X, score_func="local_score_BIC_from_cov", parameters=params, depth=3,
                  verbose=False)
    else:
        raise ValueError(name)
    g = _from_general_graph(G, nodes)
    return g, {"n_edges": len(g.skeleton())}


# ---------------------------------------------------------------------------
# registry


def alpha_n(n):
    return 1.0 / sqrt(n)


LEARNERS = {
    # E5's learner exactly (alpha 0.01, conditioning sets of size <= 3)
    "pc01d3": lambda nodes, X: pc(nodes, X, alpha=0.01, max_depth=3),
    # the rest have no depth cap (PC-stable runs until no pair has enough neighbours)
    "pc01": lambda nodes, X: pc(nodes, X, alpha=0.01, max_depth=None),
    "pc001": lambda nodes, X: pc(nodes, X, alpha=0.001, max_depth=None),
    "pcan": lambda nodes, X: pc(nodes, X, alpha=alpha_n(X.shape[0]), max_depth=None),
    "cpc": lambda nodes, X: pc(nodes, X, alpha=alpha_n(X.shape[0]), max_depth=None, rule="cons"),
    "mpc": lambda nodes, X: pc(nodes, X, alpha=alpha_n(X.shape[0]), max_depth=None, rule="maj"),
    "hc": lambda nodes, X: hc(nodes, X),
    "ges": lambda nodes, X: causal_learn("ges", nodes, X),
    "boss": lambda nodes, X: causal_learn("boss", nodes, X),
    "grasp": lambda nodes, X: causal_learn("grasp", nodes, X),
}


def learn(name, nodes, X):
    return LEARNERS[name](nodes, X)


# ---------------------------------------------------------------------------
# helpers for the analysis


def dor_tarsi(g: MPDAG):
    """One consistent DAG extension of a closed PDAG (Dor & Tarsi 1992), or None."""
    und = {v: set() for v in g.nodes}
    for a, b in g.undirected_edges:
        und[a].add(b); und[b].add(a)
    out_d = {v: set() for v in g.nodes}
    in_d = {v: set() for v in g.nodes}
    for a, b in g.directed_edges:
        out_d[a].add(b); in_d[b].add(a)
    alive = set(g.nodes)
    result = set(g.directed_edges)
    while alive:
        found = None
        for x in sorted(alive):
            if out_d[x] & alive:
                continue
            nb = und[x] & alive
            adjx = nb | (in_d[x] & alive)
            if all(all(w == y or g.has_edge(y, w) for w in adjx) for y in nb):
                found = (x, nb)
                break
        if found is None:
            return None
        x, nb = found
        for y in nb:
            result.add((y, x))
        alive.discard(x)
    return MPDAG(g.nodes, directed=result)
