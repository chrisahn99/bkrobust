"""E12: the radius on the five published networks above 500 nodes.

link (724 nodes), munin (1,041), munin2, munin3, munin4 (1,003-1,041) were
left out of the paper's per-query records. For each, the 20 frozen
CPDAG-only frame queries (results/frame/frame.jsonl) under the eight
real-name elicited conditions (results/elicit/knowledge.json, the paper's
build_state): network size, undirected edges, largest chain component; the
committed set Z = O(G0) (closed form, E13 b; the paper's enumeration-based
O is intractable on components this size), its validity at G0, the CPDAG
check, the radius by the submitted hybrid (GAC, SAT rung limit 120 s, whole
query killed at 900 s), validity in the data-generating DAG, and time per call.
"""

from __future__ import annotations

import csv
import json
import multiprocessing as mp
import signal
import time

from common import OUT, ROOT, UNREACHED, Timer, hybrid, is_gac_valid_mpdag, k_of
from oset import optimal_set_closed_form

from bkrobust.demo.graph import undirected_components
from bkrobust.robustness import llm_survival as ls
from bkrobust.robustness import real_survival as rs

NETS = ("link", "munin", "munin2", "munin3", "munin4")
REAL = ["D_LLM", "D_LLM_INSTR", "D_LLM_SRC", "D_LLM_72B", "D_LLM_72B_INSTR", "D_LLM_32B",
        "D_LLM_SRC_70B", "D_LLM_GEMMA_27B"]


class Kill(Exception):
    pass


def _a(*_):
    raise Kill


def work(net):
    signal.signal(signal.SIGALRM, _a)
    with Timer() as tl:
        dag, cpdag = rs.load_network(net)
    comps = [len(c) for c in undirected_components(cpdag)]
    info = {"network": net, "nodes": len(cpdag.nodes), "edges": len(dag.directed_edges),
            "undirected": len(cpdag.undirected_edges), "largest_component": max(comps or [0]),
            "load_s": round(tl.s, 2)}
    frame = [json.loads(line) for line in open(ROOT / "results/frame/frame.jsonl")]
    q = [(r["x"], r["y"]) for r in frame if r["network"] == net]
    bundle, _ = ls.load_bundle()
    rows = []
    for cond in REAL:
        k, g0, why, shared = ls.build_state(bundle, cond, net, dag, cpdag)
        base = dict(info, condition=cond, n_k=len(k), g0_status=why,
                    k_g0=(len(k_of(cpdag, g0)) if g0 is not None else None))
        for x, y in q:
            row = dict(base, x=x, y=y)
            if g0 is None:
                rows.append(dict(row, status=f"g0_{why}"))
                continue
            with Timer() as to:
                o, st = optimal_set_closed_form(g0, x, y)
            row["o_s"] = round(to.s, 4)
            if st != "ok":
                rows.append(dict(row, status=st))
                continue
            if not is_gac_valid_mpdag(g0, x, y, o):
                rows.append(dict(row, status="z_invalid_at_g0"))
                continue
            row.update(status="ok", z_size=len(o))
            signal.alarm(900)
            t = time.perf_counter()
            try:
                h = hybrid(cpdag, g0, x, y, o, time_limit_s=120)
                signal.alarm(0)
                row.update(r=("inf" if h.radius == UNREACHED else h.radius), method=h.method,
                           exact=h.exact, top_s=round(h.top_seconds, 4),
                           total_s=round(time.perf_counter() - t, 4))
            except Kill:
                row.update(r="", method="killed", exact=False, total_s=900)
            row["valid_truth"] = is_gac_valid_mpdag(dag, x, y, o)
            rows.append(row)
    print(net, info, flush=True)
    return rows


if __name__ == "__main__":
    res = []
    with mp.Pool(2) as P:
        for r in P.imap_unordered(work, NETS):
            res.extend(r)
    keys = sorted({k for r in res for k in r})
    p = OUT / "e12" / "large_networks.csv"
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(res)
    print("wrote", p)
