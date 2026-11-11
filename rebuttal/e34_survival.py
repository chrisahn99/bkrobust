"""E34: does the radius still rank survival under knowledge corruption on learned classes?

Reads the learned CPDAGs stored by E30 (results/rebuttal/e30/rows.jsonl); no discovery
is rerun. Units: committed sets (status ok) of the TRUE (registered) and LLM arms, for
learners pcan, ges, boss at n in {500, 5000, 50000}, seeds 0-4, plus the true CPDAG.

For each unit the claims K are the orientations the arm actually applied (sorted
greedy, as e5_discovery.apply_greedy). The paper's flip arm
(bkrobust.robustness.real_survival): at each depth of depth_grid(|K|), DRAWS draws of
flip_draw (family = unit id), Meek closure of the corrupted claims on C_hat, survival =
share of non-contradictory draws in which Z stays GAC-valid; survival AUC =
auc_over_grid over the depth grid. One JSON line per unit.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import time

import discovery as D
from common import OUT, is_gac_valid_mpdag
from e5_discovery import llm_claims

from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.robustness import real_survival as rs

ODIR = OUT / "e34"
LEARNERS = ("pcan", "ges", "boss")
NS = (500, 5000, 50000)
SEEDS = range(5)


def applied(cpdag, claims):
    g, kept = cpdag, []
    for a, b in sorted(claims):
        if not g.is_undirected_edge(a, b):
            continue
        h = apply_orientations(g, [(a, b)])
        if h is not None:
            g, kept = h, kept + [(a, b)]
    return g, kept


def arm_claims(net, cp, dag, arm):
    if arm == "TRUE":
        return [(a, b) if dag.is_directed_edge(a, b) else (b, a)
                for a, b in cp.undirected_edges if dag.has_edge(a, b)]
    return llm_claims(net)


def unit_task(u):
    t0 = time.perf_counter()
    cp = MPDAG(D.true_dag(u["network"]).nodes, directed=[tuple(e) for e in u["cpdag"]["directed"]],
               undirected=[tuple(e) for e in u["cpdag"]["undirected"]])
    dag = D.true_dag(u["network"])
    g0, k = applied(cp, arm_claims(u["network"], cp, dag, u["arm"]))
    z = frozenset(u["z"].split(",")) if u["z"] else frozenset()
    assert is_gac_valid_mpdag(g0, u["x"], u["y"], z), "Z not valid at the rebuilt G0"
    out = {key: u[key] for key in ("graph_id", "network", "n", "seed", "learner", "arm", "x", "y",
                                   "z", "r", "valid_truth")}
    out["n_k"] = len(k)
    if not k:
        out["status"] = "no_claims"
        return out
    fam = f"e34|{u['graph_id']}|{u['arm']}|{u['x']}|{u['y']}"
    cells, cache = {}, {}
    grid = rs.depth_grid(len(k))
    for d in grid:
        nc = ns = 0
        for rep in range(u["draws"]):
            kc, _ = rs.flip_draw(k, d, rep, family=fam)
            key = tuple(sorted(kc))
            if key not in cache:
                g = apply_orientations(cp, list(key))
                cache[key] = None if g is None else bool(is_gac_valid_mpdag(g, u["x"], u["y"], z))
            v = cache[key]
            if v is None:
                nc += 1
            else:
                ns += v
        cells[d] = rs.cell_statistics(u["draws"], nc, ns)
    out.update(status="ok", grid=list(grid), S={str(d): c["S"] for d, c in cells.items()},
               contra={str(d): c["contradiction_rate"] for d, c in cells.items()},
               auc=rs.auc_over_grid(cells, grid), secs_laptop=round(time.perf_counter() - t0, 2))
    return out


def units(draws, max_per_graph):
    src = OUT / "e30" / "rows.jsonl"
    out = []
    for line in src.open():
        r = json.loads(line)
        if "error" in r:
            continue
        keep = r["n"] == "oracle" or (r["learner"] in LEARNERS and r["n"] in NS and r["seed"] in SEEDS)
        if not keep:
            continue
        for arm in ("TRUE", "LLM"):
            rows = [x for x in r["rows"] if x["arm"] == arm and x["status"] == "ok"]
            for x in rows[:max_per_graph]:
                out.append({"graph_id": r["id"], "network": r["network"], "n": r["n"], "seed": r["seed"],
                            "learner": r["learner"], "arm": arm, "x": x["x"], "y": x["y"], "z": x["z"],
                            "r": x["r"], "valid_truth": x["valid_truth"], "cpdag": r["cpdag"],
                            "draws": draws})
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--draws", type=int, default=200)
    ap.add_argument("--max-per-graph", type=int, default=40)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", default="units.jsonl")
    a = ap.parse_args()
    U = units(a.draws, a.max_per_graph)
    ODIR.mkdir(parents=True, exist_ok=True)
    path = ODIR / a.out
    done = set()
    if path.exists():
        for line in path.open():
            r = json.loads(line)
            done.add((r["graph_id"], r["arm"], r["x"], r["y"]))
    U = [u for u in U if (u["graph_id"], u["arm"], u["x"], u["y"]) not in done]
    if a.limit:
        U = U[: a.limit]
    print(len(U), "units to run", flush=True)
    t0 = time.time()
    with mp.Pool(a.workers, maxtasksperchild=50) as P, path.open("a") as f:
        for i, r in enumerate(P.imap_unordered(unit_task, U, chunksize=2)):
            f.write(json.dumps(r, default=str) + "\n")
            f.flush()
            if i % 50 == 0:
                print(i, r["graph_id"], r["arm"], r.get("secs_laptop"), f"{time.time() - t0:.0f}s", flush=True)
