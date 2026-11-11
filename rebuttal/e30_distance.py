"""E30 exploratory (not registered): where is the nearest class error?

For every committed set of the TRUE arm, the hop distance, in the true skeleton, from
{x, y} u Z to the nearest pair on which C_hat and the true CPDAG differ (adjacency or
mark): 0 if an endpoint lies in {x, y} u Z. Output: results/rebuttal/e30/distance_rows.csv.
"""

from __future__ import annotations

import json
from collections import deque

import pandas as pd

import discovery as D
from common import OUT

from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG


def marks(g):
    m = {}
    for a, b in g.directed_edges:
        m[frozenset((a, b))] = (a, b)
    for a, b in g.undirected_edges:
        m[frozenset((a, b))] = "u"
    return m


def bfs(dag, src):
    dist = {v: 0 for v in src}
    q = deque(src)
    while q:
        v = q.popleft()
        for w in dag.adjacent(v):
            if w not in dist:
                dist[w] = dist[v] + 1
                q.append(w)
    return dist


def main():
    truth = {net: D.true_dag(net) for net in D.GAUSSIAN}
    tm = {net: marks(dag_to_cpdag(d)) for net, d in truth.items()}
    out = []
    for line in (OUT / "e30" / "rows.jsonl").open():
        r = json.loads(line)
        if "rows" not in r or r["n"] == "oracle":
            continue
        dag = truth[r["network"]]
        cp = MPDAG(dag.nodes, directed=[tuple(e) for e in r["cpdag"]["directed"]],
                   undirected=[tuple(e) for e in r["cpdag"]["undirected"]])
        lm = marks(cp)
        t = tm[r["network"]]
        errs = [k for k in set(lm) | set(t) if lm.get(k) != t.get(k)]
        for x in r["rows"]:
            if x["arm"] != "TRUE" or x["status"] != "ok":
                continue
            S = {x["x"], x["y"]} | (set(x["z"].split(",")) if x["z"] else set())
            dist = bfs(dag, S)
            d = min((min(dist.get(a, 99), dist.get(b, 99)) for a, b in map(tuple, errs)), default=99)
            out.append({"network": r["network"], "n": r["n"], "seed": r["seed"], "learner": r["learner"],
                        "x": x["x"], "y": x["y"], "valid_truth": x["valid_truth"],
                        "class_ok": x.get("class_ok"), "d_err": d, "n_err": len(errs)})
    d = pd.DataFrame(out)
    d.to_csv(OUT / "e30" / "distance_rows.csv", index=False)
    d["inv"] = ~d.valid_truth.astype(bool)
    d["bin"] = d.d_err.clip(upper=3)
    print(d.groupby("bin").agg(rows=("inv", "size"), invalid=("inv", "mean")).round(3))


if __name__ == "__main__":
    main()
