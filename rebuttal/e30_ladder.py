"""E30-E33: sample-size ladder x CPDAG learners, on the four fitted Gaussian networks.

One task = (network, n, seed, learner). Data are drawn with discovery.sample (the
E5 sampler; seeds 0-4 at n in {200, 1000, 5000} are E5's samples). The query
frame is E5's (e5_discovery.frame, 40 pairs per network) and the knowledge arms
are E5's (TRUE, LLM, NONE; e5_discovery.arms). Per task one JSON line holds the
graph-level metrics and one record per (arm, query).

Per committed set (status ok): Z = O(G0), validity in the true DAG (GAC), the
radius on the learned class (hybrid, 30 s per rung), and the local check:
  A = An_true({x, y} u Z) (inclusive).
  class_ok  : no true edge inside A is missing from C_hat, no directed edge of
              C_hat inside A reverses a true edge, and orienting every undirected
              C_hat edge inside A that is a true edge as in the truth Meek-closes.
              Then some D in [C_hat] contains the truth on A with its orientations.
  g0_ok     : the same with G0 in place of C_hat (class and knowledge both right on A).
  k_err     : |K_G0 minus dir(D)| for one such D (closure of the class_ok orientations
              plus the G0 closure orientations that still close, in sorted order,
              then a Dor-Tarsi extension).
  Theorem (Remark 8(i) localised, with Corollary 1(c)): class_ok and (r = inf or
  k_err < r) imply Z valid in the truth. A violation is a bug (consistency check).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import os
import time
import traceback

import numpy as np

import discovery as D
import learners as L
from common import OUT, UNREACHED, hybrid, is_gac_valid_mpdag, k_of
from e5_discovery import arms, frame
from oset import optimal_set_closed_form

from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations

NS = (100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000)
LEARNER_ORDER = ("pc01d3", "pc01", "pc001", "pcan", "cpc", "mpc", "hc", "ges", "boss", "grasp")
RUNG_S = 30
ODIR = OUT / "e30"


def vstructs(g):
    out = set()
    for b in g.nodes:
        pa = sorted(g.parents(b))
        for i, a in enumerate(pa):
            for c in pa[i + 1:]:
                if not g.has_edge(a, c):
                    out.add((a, b, c))
    return out


def graph_metrics(cp, true_cp):
    ls = {frozenset(e) for e in cp.skeleton()}
    ts = {frozenset(e) for e in true_cp.skeleton()}
    tp = len(ls & ts)
    shd = len(ls ^ ts)
    for e in ls & ts:
        a, b = tuple(e)
        mark = lambda g: ("d", a, b) if g.is_directed_edge(a, b) else (("d", b, a) if g.is_directed_edge(b, a) else "u")
        shd += mark(cp) != mark(true_cp)
    lv, tv = vstructs(cp), vstructs(true_cp)
    ext = L.dor_tarsi(cp)
    return {"skel_prec": tp / len(ls) if ls else None, "skel_rec": tp / len(ts),
            "missing": len(ts - ls), "extra": len(ls - ts),
            "vs_prec": len(lv & tv) / len(lv) if lv else None, "vs_rec": len(lv & tv) / len(tv) if tv else None,
            "shd_cpdag": shd, "n_undirected": len(cp.undirected_edges),
            "has_extension": ext is not None,
            "is_cpdag": ext is not None and dag_to_cpdag(ext) == cp}


def anc_set(dag, nodes):
    out = set(nodes)
    for v in nodes:
        out |= dag.ancestors(v)
    return out


def local_check(cp, g0, dag, x, y, z):
    A = anc_set(dag, {x, y} | set(z))
    tedges = [(a, b) for a, b in dag.directed_edges if a in A and b in A]
    res = {"A_size": len(A), "A_true_edges": len(tedges)}
    for name, g in (("class", cp), ("g0", g0)):
        miss = sum(not g.has_edge(a, b) for a, b in tedges)
        rev = sum(g.is_directed_edge(b, a) for a, b in tedges)
        orient = [(a, b) for a, b in tedges if g.is_undirected_edge(a, b)]
        gc = apply_orientations(g, orient) if (miss == 0 and rev == 0) else None
        res[f"{name}_missing_A"] = miss
        res[f"{name}_reversed_A"] = rev
        res[f"{name}_ok"] = gc is not None
        if name == "class":
            gclass = gc
    res["k_err"] = None
    if gclass is not None:
        g = gclass
        for a, b in sorted(k_of(cp, g0)):
            if g.is_undirected_edge(a, b):
                h = apply_orientations(g, [(a, b)])
                if h is not None:
                    g = h
        d = L.dor_tarsi(g)
        if d is not None:
            res["k_err"] = len(k_of(cp, g0) - set(d.directed_edges))
    return res


def hop_sets(dag, src, h):
    seen = set(src)
    front = set(src)
    for _ in range(h):
        front = {w for v in front for w in dag.adjacent(v)} - seen
        seen |= front
    return seen


def local_agree(cp, true_cp, region):
    for a in region:
        for b in region:
            if a < b:
                ma = (cp.has_edge(a, b), cp.is_directed_edge(a, b), cp.is_directed_edge(b, a))
                mb = (true_cp.has_edge(a, b), true_cp.is_directed_edge(a, b), true_cp.is_directed_edge(b, a))
                if ma != mb:
                    return False
    return True


def score(cp, dag, true_cp, net, with_local=True):
    out = []
    g0s = arms(net, cp, dag)
    for arm, g0 in g0s.items():
        for x, y in frame(net):
            row = {"arm": arm, "x": x, "y": y}
            o, st = optimal_set_closed_form(g0, x, y)
            if st == "ok" and not is_gac_valid_mpdag(g0, x, y, o):
                st = "z_invalid_at_g0"
            row["status"] = st
            if st == "ok":
                t = time.perf_counter()
                h = hybrid(cp, g0, x, y, o, time_limit_s=RUNG_S)
                row.update(z=",".join(sorted(o)), z_size=len(o),
                           r=("inf" if h.radius == UNREACHED else int(h.radius)),
                           r_exact=bool(h.exact), r_method=str(h.method),
                           r_secs=round(time.perf_counter() - t, 3),
                           n_k_g0=len(k_of(cp, g0)),
                           valid_truth=bool(is_gac_valid_mpdag(dag, x, y, o)))
                if with_local:
                    row.update(local_check(cp, g0, dag, x, y, o))
                    for hh in (1, 2):
                        row[f"agree_hop{hh}"] = local_agree(cp, true_cp, hop_sets(dag, {x, y}, hh))
            out.append(row)
    return out


def task_id(net, n, seed, learner):
    return f"{net}|{n}|{seed}|{learner}"


def run_task(task):
    net, n, seed, learner = task
    try:
        dag = D.true_dag(net)
        true_cp = dag_to_cpdag(dag)
        if n == "oracle":
            cp, info, secs = true_cp, {}, 0.0
        else:
            nodes, X = D.sample(net, n, seed)
            t = time.perf_counter()
            cp, info = L.learn(learner, nodes, X)
            secs = time.perf_counter() - t
        rec = {"id": task_id(net, n, seed, learner), "network": net, "n": n, "seed": seed,
               "learner": learner, "learn_secs_laptop": round(secs, 2),
               **{f"info_{k}": v for k, v in info.items()}, **graph_metrics(cp, true_cp)}
        t = time.perf_counter()
        rec["rows"] = score(cp, dag, true_cp, net)
        rec["score_secs_laptop"] = round(time.perf_counter() - t, 2)
        rec["cpdag"] = {"directed": sorted(cp.directed_edges), "undirected": sorted(cp.undirected_edges)}
        return rec
    except Exception as e:  # noqa: BLE001
        return {"id": task_id(net, n, seed, learner), "network": net, "n": n, "seed": seed,
                "learner": learner, "error": f"{type(e).__name__}: {e}", "tb": traceback.format_exc()}


def tasks(seeds, ns, learners, nets):
    out = [(net, "oracle", "oracle", "oracle") for net in nets]
    for net in nets:
        for n in ns:
            for s in range(seeds):
                for lr in learners:
                    out.append((net, n, s, lr))
    return out


def cost(t):
    net, n, s, lr = t
    big = net in ("arth150", "magic-irri")
    return (lr == "ges") * (100 if big else 5) + (lr in ("cpc", "mpc")) * 10 + 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--ns", default=",".join(map(str, NS)))
    ap.add_argument("--learners", default=",".join(LEARNER_ORDER))
    ap.add_argument("--nets", default=",".join(D.GAUSSIAN))
    ap.add_argument("--out", default="rows.jsonl")
    a = ap.parse_args()
    ODIR.mkdir(parents=True, exist_ok=True)
    path = ODIR / a.out
    done = set()
    if path.exists():
        for line in path.open():
            r = json.loads(line)
            if "error" not in r:
                done.add(r["id"])
    T = [t for t in tasks(a.seeds, [int(x) for x in a.ns.split(",")], a.learners.split(","),
                          a.nets.split(",")) if task_id(*t) not in done]
    T.sort(key=cost, reverse=True)
    print(len(T), "tasks to run;", len(done), "done", flush=True)
    t0 = time.time()
    with mp.Pool(a.workers, maxtasksperchild=20) as P, path.open("a") as f:
        for i, r in enumerate(P.imap_unordered(run_task, T)):
            f.write(json.dumps(r, default=str) + "\n")
            f.flush()
            msg = r.get("error", f"learn {r.get('learn_secs_laptop')}s score {r.get('score_secs_laptop')}s")
            print(i, r["id"], msg, f"elapsed {time.time() - t0:.0f}s", flush=True)
