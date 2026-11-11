"""Validation of rebuttal/learners.py (run before E30; output in results/rebuttal/e30/checks.json).

1. pc(rule="sepset", alpha=0.01, depth 3) equals discovery.pc on 12 samples
   (4 networks x n in {200, 5000} x seeds {0, 1}, partly) -- the E5 learner.
2. Oracle d-separation tests: pc with each collider rule returns the true CPDAG,
   at max depth 3 and at the maximum degree.
3. Score-based learners at n = 10^6 (seed 0): distance to the true CPDAG.
4. dor_tarsi returns a consistent extension of every true CPDAG.
"""

from __future__ import annotations

import json
import sys
import time

import networkx as nx

import discovery as D
import learners as L
from common import OUT

from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import is_consistent_extension


class Oracle:
    def __init__(self, dag, nodes):
        self.G = nx.DiGraph()
        self.G.add_nodes_from(nodes)
        self.G.add_edges_from(dag.directed_edges)
        self.nodes = nodes

    def pval(self, i, j, S):
        return 1.0 if nx.is_d_separator(self.G, {self.nodes[i]}, {self.nodes[j]},
                                        {self.nodes[k] for k in S}) else 0.0


def cmp(g, cp):
    return {"equal": g == cp,
            "skel_missing": len(set(map(frozenset, cp.skeleton())) - set(map(frozenset, g.skeleton()))),
            "skel_extra": len(set(map(frozenset, g.skeleton())) - set(map(frozenset, cp.skeleton()))),
            "dir_diff": len(set(g.directed_edges) ^ set(cp.directed_edges)),
            "und_diff": len(set(g.undirected_edges) ^ set(cp.undirected_edges))}


def main(parts):
    res = {}
    if "1" in parts:
        r = []
        for net in D.GAUSSIAN:
            for n in (200, 5000):
                for s in (0, 1):
                    nodes, X = D.sample(net, n, s)
                    a, _ = D.pc(nodes, X, alpha=0.01)
                    b, _ = L.pc(nodes, X, alpha=0.01)
                    r.append({"net": net, "n": n, "seed": s, "equal": a == b})
        res["pc_reproduces_discovery_pc"] = r
        print("1", all(x["equal"] for x in r), flush=True)
    if "2" in parts:
        r = []
        orig = D.FisherZ
        for net in D.GAUSSIAN:
            dag = D.true_dag(net)
            cp = dag_to_cpdag(dag)
            nodes, X = D.sample(net, 20, 0)
            deg = max(len(dag.adjacent(v)) for v in dag.nodes)
            D.FisherZ = lambda X_, dag=dag, nodes=nodes: Oracle(dag, nodes)
            for rule in ("sepset", "cons", "maj"):
                for depth in (3, deg):
                    g, info = L.pc(nodes, X, alpha=0.5, max_depth=depth, rule=rule)
                    r.append({"net": net, "rule": rule, "max_depth": depth, **cmp(g, cp),
                              "n_ambiguous": info["n_ambiguous"]})
                    print(r[-1], flush=True)
            D.FisherZ = orig
        res["oracle_pc"] = r
    if "3" in parts:
        r = []
        for net in D.GAUSSIAN:
            dag = D.true_dag(net)
            cp = dag_to_cpdag(dag)
            nodes, X = D.sample(net, 1_000_000, 0)
            for name in ("hc", "ges", "boss", "grasp", "pcan", "cpc", "mpc"):
                t = time.perf_counter()
                g, _ = L.learn(name, nodes, X)
                r.append({"net": net, "learner": name, "secs_laptop": round(time.perf_counter() - t, 1),
                          **cmp(g, cp)})
                print(r[-1], flush=True)
        res["large_n"] = r
    if "4" in parts:
        r = []
        for net in D.GAUSSIAN:
            cp = dag_to_cpdag(D.true_dag(net))
            e = L.dor_tarsi(cp)
            r.append({"net": net, "consistent": e is not None and is_consistent_extension(e, cp)})
        res["dor_tarsi"] = r
        print("4", r, flush=True)
    path = OUT / "e30" / f"checks_{''.join(parts)}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(sys.argv[1:] or ["1", "2", "3", "4"])
