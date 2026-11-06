"""E3: depth budget b and SAT mode.

Pool: dense Erdos-Renyi instances (submitted generator), n in 10..24, edge
probability 0.5..0.85, full or half recovering knowledge, seeds 0..S. Queries
with a finite radius (Z valid and non-empty at G0, invalid at the CPDAG) are
classified by the default hybrid; all deep ones (r >= 4, i.e. answered by the
SAT stage at b = 3) are kept plus an equal-size seeded sample of shallow ones.

Arms: b in {0, 1, 2, 3, 4, 6} x SAT mode in {ladder, minimise}.
  ladder   = the submitted E1 encoding, rungs k = b+1, b+2, ... (Algorithm 1);
  minimise = the same encoding, one CP-SAT call minimising the number of
             retractions (no rungs).
Every arm must return the same radius (any difference is a bug).
"""

from __future__ import annotations

import argparse
import csv
import multiprocessing as mp
import time

import numpy as np
from ortools.sat.python import cp_model

from common import OUT, UNREACHED, hybrid, is_gac_valid_mpdag
from oset import optimal_set_closed_form

from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover
from bkrobust.demo.meek import apply_orientations
from bkrobust.sat.e1 import _build
from bkrobust.search.exact_fast import radius_local_up_fast
from bkrobust.synth.generators import erdos_renyi_dag

BUDGETS = (0, 1, 2, 3, 4, 6)
RUNG_S = 120.0


def sat_ladder(cpdag, g0, x, y, z, start):
    model, _ov, retracted, _d = _build(cpdag, g0, x, y, z, criterion="gac")
    idx = [r.Index() for r in retracted]
    for k in range(start, len(retracted) + 1):
        sub = model.Clone()
        if idx:
            sub.Add(sum(sub.GetBoolVarFromProtoIndex(i) for i in idx) <= k)
        s = cp_model.CpSolver()
        s.parameters.num_workers = 1
        s.parameters.random_seed = 0
        s.parameters.max_time_in_seconds = RUNG_S
        st = s.Solve(sub)
        if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            return k, True
        if st != cp_model.INFEASIBLE:
            return None, False
    return UNREACHED, True


def sat_minimise(cpdag, g0, x, y, z, start):
    model, _ov, retracted, _d = _build(cpdag, g0, x, y, z, criterion="gac")
    if start > 0 and retracted:
        model.Add(sum(retracted) >= start)
    model.Minimize(sum(retracted))
    s = cp_model.CpSolver()
    s.parameters.num_workers = 1
    s.parameters.random_seed = 0
    s.parameters.max_time_in_seconds = RUNG_S * 4
    st = s.Solve(model)
    if st == cp_model.OPTIMAL:
        return int(round(s.ObjectiveValue())), True
    if st == cp_model.INFEASIBLE:
        return UNREACHED, True
    return None, False


def arm(cpdag, g0, x, y, z, b, mode):
    """Algorithm 1 with depth budget b and the given SAT mode. Returns (r, exact, s, leg)."""
    t = time.perf_counter()
    if is_gac_valid_mpdag(cpdag, x, y, z):
        return UNREACHED, True, time.perf_counter() - t, "top"

    def fails(g):
        return not is_gac_valid_mpdag(g, x, y, z)

    if b > 0:
        f = radius_local_up_fast(cpdag, g0, fails, max_depth=b)
        if f.radius != UNREACHED:
            return f.radius, True, time.perf_counter() - t, "search"
        if f.exact:
            return UNREACHED, True, time.perf_counter() - t, "search"
    fn = sat_ladder if mode == "ladder" else sat_minimise
    r, ex = fn(cpdag, g0, x, y, z, b + 1 if b > 0 else 0)
    return r, ex, time.perf_counter() - t, "sat"


def pool_instance(spec):
    n, p, seed, half = spec
    rng = np.random.default_rng(seed)
    dag = erdos_renyi_dag(n, rng, p)
    cpdag = dag_to_cpdag(dag)
    k = sorted(knowledge_to_recover(dag, cpdag))
    if half:
        k = k[: len(k) // 2]
    g0 = apply_orientations(cpdag, k)
    q = []
    r2 = np.random.default_rng(777 + seed)
    nodes = list(g0.nodes)
    pairs = [(a, c) for a in nodes for c in nodes if a != c]
    for i in r2.permutation(len(pairs))[:60]:
        x, y = pairs[i]
        o, st = optimal_set_closed_form(g0, x, y)
        if st != "ok" or not o or not is_gac_valid_mpdag(g0, x, y, o):
            continue
        if is_gac_valid_mpdag(cpdag, x, y, o):
            continue
        h = hybrid(cpdag, g0, x, y, o, time_limit_s=RUNG_S)
        q.append({"n": n, "p": p, "seed": seed, "half": half, "x": x, "y": y,
                  "z": sorted(o), "r_default": h.radius, "method_default": h.method,
                  "k_g0": len(set(g0.directed_edges) - set(cpdag.directed_edges))})
    return q


def run_arms(q):
    rng = np.random.default_rng(q["seed"])
    dag = erdos_renyi_dag(q["n"], rng, q["p"])
    cpdag = dag_to_cpdag(dag)
    k = sorted(knowledge_to_recover(dag, cpdag))
    if q["half"]:
        k = k[: len(k) // 2]
    g0 = apply_orientations(cpdag, k)
    z = frozenset(q["z"])
    out = []
    for b in BUDGETS:
        for mode in ("ladder", "minimise"):
            r, ex, s, leg = arm(cpdag, g0, q["x"], q["y"], z, b, mode)
            out.append({**{kk: q[kk] for kk in ("n", "p", "seed", "half", "x", "y", "k_g0",
                                                 "r_default")},
                        "b": b, "mode": mode, "r": r, "exact": ex, "seconds": round(s, 6),
                        "leg": leg, "agrees": r == q["r_default"]})
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--seeds", type=int, default=12)
    ap.add_argument("--max-deep", type=int, default=80)
    a = ap.parse_args()
    specs = [(n, p, s, h) for n in (10, 12, 14, 16, 18, 20, 24) for p in (0.5, 0.7, 0.85)
             for s in range(a.seeds) for h in (False, True)]
    pool = []
    with mp.Pool(a.workers) as P:
        for qs in P.imap_unordered(pool_instance, specs):
            pool.extend(qs)
    deep = [q for q in pool if q["r_default"] >= 4]
    shallow = [q for q in pool if 1 <= q["r_default"] <= 3]
    rng = np.random.default_rng(0)
    deep = [deep[i] for i in rng.permutation(len(deep))[: a.max_deep]]
    shallow = [shallow[i] for i in rng.permutation(len(shallow))[: max(len(deep), 40)]]
    print(f"pool {len(pool)} finite queries; deep {len(deep)}, shallow {len(shallow)}", flush=True)
    res = []
    with mp.Pool(a.workers) as P:
        for o in P.imap_unordered(run_arms, deep + shallow):
            res.extend(o)
    path = OUT / "e3" / "budget_sat.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(res[0]))
        w.writeheader()
        w.writerows(res)
    from collections import Counter
    print("radius distribution of pool:", Counter(q["r_default"] for q in pool))
    print("wrote", path)
