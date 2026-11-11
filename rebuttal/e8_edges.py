"""E8 (population level): which skeleton errors can break a committed set.

Random Erdos-Renyi DAGs (submitted generator), n in {10, 15, 20, 30}, expected
degree 2-4, 200 graphs per cell. The analyst's graph A is a DAG; Z = O(X, Y, A)
for every pair with a non-empty causal effect and valid Z. The truth T differs
from A by one edge:
  extra   T = A minus one edge (the analyst asserts an adjacency that is absent):
          the claim to test is that Z stays valid in T, always.
  missed  T = A plus one edge consistent with A's topological order (the
          analyst missed an adjacency): Z is checked in T, by the hop distance
          of the new edge from {X, Y} u Z in A.
Validity is the submitted GAC on DAGs. Population level: no data, no test.
"""

from __future__ import annotations

import json
import multiprocessing as mp
from collections import Counter, defaultdict

import numpy as np

from common import OUT, is_gac_valid_mpdag

from bkrobust.demo.graph import MPDAG
from bkrobust.demo.evaluate import optimal_adjustment_set_dag
from bkrobust.synth.generators import erdos_renyi_dag


def hop(dag, sources, edge):
    import networkx as nx
    g = nx.Graph()
    g.add_nodes_from(dag.nodes)
    g.add_edges_from(dag.skeleton())
    best = 99
    for s in sources:
        lens = nx.single_source_shortest_path_length(g, s)
        for v in edge:
            best = min(best, lens.get(v, 99))
    return best


def topo(dag):
    indeg = {v: len(dag.parents(v)) for v in dag.nodes}
    order, st = [], sorted(v for v, d in indeg.items() if d == 0)
    while st:
        v = st.pop(0)
        order.append(v)
        for c in sorted(dag.children(v)):
            indeg[c] -= 1
            if indeg[c] == 0:
                st.append(c)
    return order


def work(spec):
    n, deg, seed = spec
    rng = np.random.default_rng(seed)
    A = erdos_renyi_dag(n, rng, deg / (n - 1))
    order = topo(A)
    pos = {v: i for i, v in enumerate(order)}
    edges = sorted(A.directed_edges)
    absent = [(a, b) for a in order for b in order if pos[a] < pos[b] and not A.has_edge(a, b)]
    res = {"extra": Counter(), "missed": defaultdict(Counter)}
    queries = []
    for x in A.nodes:
        for y in sorted(A.descendants(x)):
            z = frozenset(optimal_adjustment_set_dag(A, x, y))
            if is_gac_valid_mpdag(A, x, y, z):
                queries.append((x, y, z))
    if not queries:
        return res
    qi = rng.permutation(len(queries))[:10]
    for i in qi:
        x, y, z = queries[i]
        if edges:
            e = edges[int(rng.integers(len(edges)))]
            T = MPDAG(A.nodes, directed=[f for f in edges if f != e])
            res["extra"][is_gac_valid_mpdag(T, x, y, z)] += 1
        if absent:
            e = absent[int(rng.integers(len(absent)))]
            T = MPDAG(A.nodes, directed=edges + [e])
            h = hop(A, {x, y} | set(z), e)
            res["missed"][min(h, 3)][is_gac_valid_mpdag(T, x, y, z)] += 1
    return res


if __name__ == "__main__":
    specs = [(n, d, s) for n in (10, 15, 20, 30) for d in (2, 3, 4) for s in range(200)]
    extra = Counter()
    missed = defaultdict(Counter)
    with mp.Pool(2) as P:
        for r in P.imap_unordered(work, specs, chunksize=20):
            extra.update(r["extra"])
            for h, c in r["missed"].items():
                missed[h].update(c)
    out = {"extra_edge": {"valid": extra[True], "invalid": extra[False]},
           "missed_edge_by_hop": {str(h): {"valid": c[True], "invalid": c[False],
                                           "share_invalid": c[False] / max(1, c[True] + c[False])}
                                  for h, c in sorted(missed.items())}}
    p = OUT / "e8" / "edges.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
