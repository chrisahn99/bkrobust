"""Sampling from the fitted Gaussian networks and an in-house PC-stable.

Sampling: the four bnlearn Gaussian networks shipped in
results/axisa3/networks/example_models (ecoli70, magic-irri, magic-niab,
arth150) carry fitted linear-Gaussian CPDs (intercept, coefficients,
variance); data are drawn ancestrally from them.

PC-stable (Colombo & Maathuis 2014): Fisher-z partial-correlation tests at
level alpha, order-independent adjacency search up to --max-depth, unshielded
triples oriented as v-structures when the middle node is not in the
separating set (conflicting orientations are left undirected), then Meek
closure with the submitted ``bkrobust.demo.meek.meek_closure``. If the closure
FAILs, v-structures are dropped one at a time, weakest first, until it
succeeds. (Implemented as greedy insertion, strongest first, each kept only
if the closure succeeds.) No package beyond numpy/scipy is used.
"""

from __future__ import annotations

import json
from itertools import combinations
from math import atanh, sqrt

import numpy as np
from scipy.stats import norm

from common import ROOT

from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import meek_closure

NETDIR = ROOT / "results/axisa3/networks/example_models"
GAUSSIAN = ("ecoli70", "magic-irri", "magic-niab", "arth150")


def load_gbn(name):
    d = json.loads((NETDIR / f"{name}.json").read_text())
    return d["nodes"], d["cpds"]


def _topo(nodes, cpds):
    order, done = [], set()
    while len(order) < len(nodes):
        for v in nodes:
            if v not in done and all(p in done for p in cpds[v]["parents"]):
                order.append(v)
                done.add(v)
    return order


def sample(name, n, seed):
    nodes, cpds = load_gbn(name)
    rng = np.random.default_rng(seed)
    idx = {v: i for i, v in enumerate(nodes)}
    X = np.zeros((n, len(nodes)))
    for v in _topo(nodes, cpds):
        c = cpds[v]
        val = np.full(n, c["coefficients"]["(Intercept)"][0])
        for p in c["parents"]:
            val += c["coefficients"][p][0] * X[:, idx[p]]
        val += rng.normal(0, sqrt(c["variance"][0]), n)
        X[:, idx[v]] = val
    return nodes, X


def true_dag(name):
    nodes, cpds = load_gbn(name)
    return MPDAG(nodes, directed=[(p, v) for v in nodes for p in cpds[v]["parents"]])


class FisherZ:
    def __init__(self, X):
        self.n = X.shape[0]
        self.C = np.corrcoef(X, rowvar=False)

    def pval(self, i, j, S):
        idx = [i, j, *S]
        sub = self.C[np.ix_(idx, idx)]
        try:
            P = np.linalg.inv(sub)
        except np.linalg.LinAlgError:
            P = np.linalg.pinv(sub)
        r = -P[0, 1] / sqrt(abs(P[0, 0] * P[1, 1]))
        r = min(0.9999999, max(-0.9999999, r))
        z = sqrt(max(self.n - len(S) - 3, 1)) * abs(atanh(r))
        return 2 * (1 - norm.cdf(z))


def pc(nodes, X, alpha=0.01, max_depth=3):
    """Returns (cpdag MPDAG, info)."""
    p = len(nodes)
    t = FisherZ(X)
    adj = {i: set(range(p)) - {i} for i in range(p)}
    sep = {}
    minp = {}
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
    # v-structures
    arrows = {}
    vs = []
    for b in range(p):
        for a, c in combinations(sorted(adj[b]), 2):
            if c in adj[a]:
                continue
            if b not in sep.get(frozenset((a, c)), set()):
                strength = minp.get(frozenset((a, c)), 1.0)
                vs.append((strength, a, b, c))
    vs.sort(reverse=True)  # strongest (largest separating p-value) first

    def build(vlist):
        head = {}
        for _, a, b, c in vlist:
            for u in (a, c):
                head.setdefault(frozenset((u, b)), set()).add(b)
                arrows[(u, b)] = True
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
        if not g.is_acyclic():
            return None
        return meek_closure(g)

    # Greedy, strongest first: keep a v-structure only if the graph stays
    # acyclic and Meek-closable with it.
    used, dropped = [], 0
    g = build(used)
    for v in vs:
        h = build(used + [v])
        if h is None:
            dropped += 1
            continue
        used.append(v)
        g = h
    if g is None:
        g = MPDAG(nodes, undirected=[(nodes[i], nodes[j]) for i in range(p) for j in adj[i] if i < j])
    return g, {"n_vstructures": len(vs), "vstructures_dropped": dropped,
               "n_edges": sum(len(a) for a in adj.values()) // 2}


def skeleton_errors(learned, dag):
    ls = {frozenset(e) for e in learned.skeleton()}
    ts = {frozenset(e) for e in dag.skeleton()}
    return len(ts - ls), len(ls - ts)
