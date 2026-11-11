"""E2: running time of the final implementation (GAC hybrid) for n up to 400.

Instances: Erdos-Renyi DAGs from the submitted generator
(bkrobust.synth.generators.erdos_renyi_dag), edge probability d/(n-1) for an
expected degree d, n in {8,...,400}; plus the dense grid of App. D for n <= 24.
G0 = Meek(CPDAG, K) with K the recovering knowledge (full) or its first half
(half), as bkrobust.gac.sweep.dense_instance builds them. Queries: pairs in a
seeded random order, the first --queries whose closed-form optimal set is
non-empty and GAC-valid at G0 (App. D rule, applied to several pairs).

Each query: one GAC validity call at G0, the CPDAG check, the depth-3 search,
the SAT ladder (per-rung limit), each timed; whole query killed at --kill-s and
recorded as censored. Laptop timings (Apple M5), one worker per process.
"""

from __future__ import annotations

import argparse
import csv
import multiprocessing as mp
import signal
import time

import numpy as np

from common import OUT, UNREACHED, Timer, hybrid, is_gac_valid_mpdag, k_of
from oset import optimal_set_closed_form

from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover
from bkrobust.demo.meek import apply_orientations
from bkrobust.synth.generators import erdos_renyi_dag


class Kill(Exception):
    pass


def _alarm(*_):
    raise Kill


def comp_stats(cpdag):
    import networkx as nx
    g = nx.Graph()
    g.add_edges_from(cpdag.undirected_edges)
    comps = [len(c) for c in nx.connected_components(g)]
    return max(comps) if comps else 0


def instance(n, p, seed, half):
    rng = np.random.default_rng(seed)
    dag = erdos_renyi_dag(n, rng, p)
    cpdag = dag_to_cpdag(dag)
    k = sorted(knowledge_to_recover(dag, cpdag))
    if half:
        k = k[: len(k) // 2]
    g0 = apply_orientations(cpdag, k)
    return dag, cpdag, g0, k


def work(spec):
    n, p, seed, half, nq, kill_s, rung_s, family = spec
    signal.signal(signal.SIGALRM, _alarm)
    out = []
    with Timer() as tb:
        dag, cpdag, g0, k = instance(n, p, seed, half)
    base = {"family": family, "n": n, "p": round(p, 4), "seed": seed,
            "knowledge": "half" if half else "full", "build_s": round(tb.s, 4),
            "n_edges": len(dag.directed_edges), "k_undirected": len(cpdag.undirected_edges),
            "k_g0": len(k_of(cpdag, g0)), "largest_component": comp_stats(cpdag)}
    rng = np.random.default_rng(10_000 + seed)
    nodes = list(g0.nodes)
    pairs = [(a, b) for a in nodes for b in nodes if a != b]
    order = rng.permutation(len(pairs))
    got = {"appD": 0, "finite": 0}
    tried = 0
    for idx in order:
        if min(got.values()) >= nq or tried > 4000:
            break
        tried += 1
        x, y = pairs[idx]
        o, st = optimal_set_closed_form(g0, x, y)
        if st != "ok" or not o or not is_gac_valid_mpdag(g0, x, y, o):
            continue
        # appD: the App. D rule (any valid non-empty O); finite: Z fails at the
        # CPDAG, so the radius is finite and the search (or SAT) must run.
        if got["appD"] < nq:
            cls = "appD"
        elif not is_gac_valid_mpdag(cpdag, x, y, o) and got["finite"] < nq:
            cls = "finite"
        else:
            continue
        got[cls] += 1
        row = dict(base, x=x, y=y, z_size=len(o), query_class=cls)
        with Timer() as tv:
            is_gac_valid_mpdag(g0, x, y, o)
        row["validity_call_s"] = round(tv.s, 6)
        signal.alarm(int(kill_s))
        t = time.perf_counter()
        try:
            h = hybrid(cpdag, g0, x, y, o, time_limit_s=rung_s)
            signal.alarm(0)
            row.update(r="inf" if h.radius == UNREACHED else h.radius, method=h.method,
                       exact=h.exact, top_s=round(h.top_seconds, 6),
                       search_s=round(h.search_seconds, 6), ladder_s=round(h.ladder_seconds, 6),
                       total_s=round(time.perf_counter() - t, 6), censored=False)
        except Kill:
            row.update(r="", method="killed", exact=False, censored=True,
                       total_s=round(time.perf_counter() - t, 3))
        out.append(row)
    if sum(got.values()) == 0:
        out.append(dict(base, x="", y="", method="no_query"))
    return out


def specs(a):
    s = []
    for n in (8, 10, 12, 14, 16, 18, 20, 24):
        for p in (0.3, 0.5, 0.7, 0.85):
            for seed in range(a.seeds):
                for half in (False, True):
                    s.append((n, p, seed, half, a.queries, a.kill_s, a.rung_s, "dense_appD"))
    for n in (16, 24, 32, 48, 64, 100, 150, 200, 300, 400):
        for d in (2, 3, 4):
            for seed in range(a.seeds):
                for half in (False, True):
                    s.append((n, d / (n - 1), seed, half, a.queries, a.kill_s, a.rung_s,
                              f"sparse_deg{d}"))
    return s


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--queries", type=int, default=3)
    ap.add_argument("--kill-s", type=float, default=600)
    ap.add_argument("--rung-s", type=float, default=120)
    ap.add_argument("--out", default="e2/timing.csv")
    a = ap.parse_args()
    sp = specs(a)
    # biggest first so the long tail starts early
    sp.sort(key=lambda s: -s[0] * (1 + s[1] * s[0]))
    print(len(sp), "instances", flush=True)
    rows = []
    path = OUT / a.out
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["family", "n", "p", "seed", "knowledge", "build_s", "n_edges", "k_undirected",
              "k_g0", "largest_component", "query_class", "x", "y", "z_size", "validity_call_s", "r",
              "method", "exact", "top_s", "search_s", "ladder_s", "total_s", "censored"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        with mp.Pool(a.workers, maxtasksperchild=20) as pool:
            for i, rs_ in enumerate(pool.imap_unordered(work, sp)):
                for r in rs_:
                    w.writerow({k: r.get(k, "") for k in fields})
                f.flush()
                if i % 20 == 0:
                    last = rs_[-1]
                    print(i, last["family"], last["n"], last.get("method"), last.get("total_s"),
                          flush=True)
    print("wrote", path)
