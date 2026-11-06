"""E5, E6, E7: the radius when the CPDAG is learned from data.

Networks: the four fitted Gaussian networks (ecoli70, magic-irri, magic-niab,
arth150). Query frame (fixed, seeded, reads the truth like the Table 1 frame):
up to 40 ordered pairs per network with y a descendant of x in the true DAG
and x incident to an undirected edge of the oracle CPDAG.

Knowledge arms, applied to whichever CPDAG is used (claims on pairs that are
not undirected there are ignored; claims that would contradict are skipped,
in sorted order):
  TRUE  the true orientation of every undirected edge (truthful analyst);
  LLM   the D_LLM elicited claims of the paper's release;
  NONE  no knowledge.
For each (network, n, seed, arm, query): Z = O(G0) (closed form, valid at
G0), whether Z is valid in the true DAG (GAC), the radius on the learned
CPDAG, and skeleton errors near the query. The oracle-CPDAG pipeline is the
reference (seed "oracle").

E6 (--bootstrap B): for each learned graph at the given n, resample rows B
times, rerun PC, rebuild G0 with the same claims, and record the share of
replicates in which the original Z is still valid at the replicate's G0
("stability"), plus the share where the replicate's radius is >= 2.

E7 (--hide): one variable hidden before PC, at hop 1 of {x, y} or at hop >= 2,
paired with the run without it (same data seed); scored at the full true DAG.
"""

from __future__ import annotations

import argparse
import csv
import json
import multiprocessing as mp

import numpy as np

import discovery as D
from common import OUT, RELEASE, UNREACHED, hybrid, is_gac_valid_mpdag
from oset import optimal_set_closed_form

from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations

NQ = 40


def frame(net):
    dag = D.true_dag(net)
    cp = dag_to_cpdag(dag)
    und_nodes = {v for e in cp.undirected_edges for v in e}
    pairs = sorted((x, y) for x in sorted(und_nodes) for y in sorted(dag.descendants(x))
                   if y != x)
    rng = np.random.default_rng(abs(hash(net)) % 1000 + 17)
    idx = rng.permutation(len(pairs))[:NQ]
    return [pairs[i] for i in sorted(idx)]


def llm_claims(net):
    ks = json.loads((RELEASE / "knowledge_states.json").read_text())
    return [tuple(e) for e in ks["D_LLM"]["networks"].get(net, {}).get("k", [])]


def apply_greedy(cpdag, claims):
    g = cpdag
    for a, b in sorted(claims):
        if not g.is_undirected_edge(a, b):
            continue
        h = apply_orientations(g, [(a, b)])
        if h is not None:
            g = h
    return g


def arms(net, cpdag, dag):
    true_k = [(a, b) if dag.is_directed_edge(a, b) else (b, a)
              for a, b in cpdag.undirected_edges if dag.has_edge(a, b)]
    return {"TRUE": apply_greedy(cpdag, true_k), "LLM": apply_greedy(cpdag, llm_claims(net)),
            "NONE": cpdag}


def near_errors(learned, dag, nodes):
    ls = {frozenset(e) for e in learned.skeleton()}
    ts = {frozenset(e) for e in dag.skeleton()}
    bad = (ts - ls) | (ls - ts)
    return sum(1 for e in bad if e & nodes)


def score(net, cpdag, dag, g0s, q, rung_s=60):
    out = []
    for arm, g0 in g0s.items():
        for x, y in q:
            row = {"arm": arm, "x": x, "y": y}
            if x not in g0.nodes or y not in g0.nodes:
                out.append(dict(row, status="hidden"))
                continue
            o, st = optimal_set_closed_form(g0, x, y)
            if st != "ok":
                out.append(dict(row, status=st))
                continue
            if not is_gac_valid_mpdag(g0, x, y, o):
                out.append(dict(row, status="z_invalid_at_g0"))
                continue
            h = hybrid(cpdag, g0, x, y, o, time_limit_s=rung_s)
            out.append(dict(row, status="ok", z=",".join(sorted(o)), z_size=len(o),
                            r=("inf" if h.radius == UNREACHED else h.radius),
                            r_exact=h.exact, r_method=h.method,
                            valid_truth=is_gac_valid_mpdag(dag, x, y, o),
                            near_skel_err=near_errors(cpdag, dag, {x, y} | set(o))))
    return out


def run_learned(task):
    net, n, seed, boot, alpha = task
    dag = D.true_dag(net)
    q = frame(net)
    nodes, X = D.sample(net, n, seed)
    cp, info = D.pc(nodes, X, alpha=alpha)
    g0s = arms(net, cp, dag)
    miss, extra = D.skeleton_errors(cp, dag)
    rows = score(net, cp, dag, g0s, q)
    base = {"network": net, "n": n, "seed": seed, "missing": miss, "extra": extra,
            "reversed": sum(1 for a, b in cp.directed_edges if dag.is_directed_edge(b, a)),
            **{f"pc_{k}": v for k, v in info.items()}}
    rows = [dict(base, **r) for r in rows]
    if boot:
        rng = np.random.default_rng(99 + seed)
        reps = []
        for _ in range(boot):
            Xb = X[rng.integers(0, n, n)]
            cpb, _ = D.pc(nodes, Xb, alpha=alpha)
            reps.append((cpb, arms(net, cpb, dag)))
        for r in rows:
            if r.get("status") != "ok":
                continue
            z = frozenset(r["z"].split(",")) if r["z"] else frozenset()
            stay = r2 = 0
            for cpb, g0b in reps:
                g = g0b[r["arm"]]
                v = is_gac_valid_mpdag(g, r["x"], r["y"], z)
                stay += v
                if v:
                    hb = hybrid(cpb, g, r["x"], r["y"], z, time_limit_s=30)
                    r2 += hb.radius == UNREACHED or hb.radius >= 2
            r["boot_stability"] = stay / boot
            r["boot_share_r_ge2"] = r2 / boot
            r["boot_B"] = boot
    return rows


def run_oracle(net):
    dag = D.true_dag(net)
    cp = dag_to_cpdag(dag)
    rows = score(net, cp, dag, arms(net, cp, dag), frame(net))
    return [dict({"network": net, "n": "oracle", "seed": "oracle"}, **r) for r in rows]


def run_hidden(task):
    """E7: hide one variable at hop 1 or >= 2 of a query, rerun PC on the rest."""
    net, n, seed, alpha = task
    dag = D.true_dag(net)
    nodes, X = D.sample(net, n, seed)
    out = []
    rng = np.random.default_rng(500 + seed)
    for x, y in frame(net)[:15]:
        near = (dag.adjacent(x) | dag.adjacent(y)) - {x, y}
        far = set(dag.nodes) - near - {x, y}
        for where, pool in (("hop1", sorted(near)), ("hop2plus", sorted(far))):
            if not pool:
                continue
            h = pool[int(rng.integers(len(pool)))]
            keep = [i for i, v in enumerate(nodes) if v != h]
            sub_nodes = [nodes[i] for i in keep]
            cp, _ = D.pc(sub_nodes, X[:, keep], alpha=alpha)
            # truthful claims on the observed part, from the full DAG's orientation
            tk = [(a, b) if dag.is_directed_edge(a, b) else (b, a)
                  for a, b in cp.undirected_edges if dag.has_edge(a, b)]
            g0 = apply_greedy(cp, tk)
            o, st = optimal_set_closed_form(g0, x, y)
            row = {"network": net, "n": n, "seed": seed, "x": x, "y": y, "where": where,
                   "hidden": h, "hidden_is_confounder": bool(
                       (h in dag.ancestors(x)) and (h in dag.ancestors(y))), "status": st}
            if st == "ok" and is_gac_valid_mpdag(g0, x, y, o):
                hr = hybrid(cp, g0, x, y, o, time_limit_s=30)
                row.update(z=",".join(sorted(o)),
                           r=("inf" if hr.radius == UNREACHED else hr.radius),
                           valid_truth=is_gac_valid_mpdag(dag, x, y, o))
            elif st == "ok":
                row["status"] = "z_invalid_at_g0"
            out.append(row)
    return out


def write(rows, name):
    path = OUT / "e5" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = sorted({k for r in rows for k in r})
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print("wrote", path, len(rows), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["learned", "hidden"], default="learned")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--boot", type=int, default=0)
    ap.add_argument("--boot-n", type=int, default=5000)
    ap.add_argument("--alpha", type=float, default=0.01)
    a = ap.parse_args()
    if a.mode == "learned":
        tasks = [(net, n, s, a.boot if (a.boot and n == a.boot_n and s < 2) else 0, a.alpha)
                 for net in D.GAUSSIAN for n in (200, 1000, 5000) for s in range(a.seeds)]
        rows = []
        with mp.Pool(a.workers) as P:
            for r in P.imap_unordered(run_oracle, D.GAUSSIAN):
                rows.extend(r)
            for i, r in enumerate(P.imap_unordered(run_learned, tasks)):
                rows.extend(r)
                print(i, r[0]["network"] if r else "", r[0]["n"] if r else "", flush=True)
        write(rows, "learned_rows.csv" if not a.boot else "learned_rows_boot.csv")
    else:
        tasks = [(net, 5000, s, a.alpha) for net in D.GAUSSIAN for s in range(a.seeds)]
        rows = []
        with mp.Pool(a.workers) as P:
            for r in P.imap_unordered(run_hidden, tasks):
                rows.extend(r)
        write(rows, "hidden_rows.csv")
