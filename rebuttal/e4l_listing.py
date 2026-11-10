"""E4L: listing the Markov equivalence class against the hybrid, at scale.

Preregistered in rebuttal/PREREGISTRATION.md (section E4L). Complements José's
E4 (branch ze/rebuttal), whose listing, pruning and counting rules are followed
so the two runs are comparable: depth-first listing of [C^] with Meek closure,
branch-and-bound on reversed K_G0 orientations, class size by root-picking
(He, Jia & Yu 2015).

Each unit is one query (C^, G0, X, Y, Z = O(G0)). Rows are appended to
results/rebuttal/e4l/rows.jsonl as they finish; a restart skips done units.

    python rebuttal/e4l_listing.py --smoke
    python rebuttal/e4l_listing.py --workers 8            # all populations
    python rebuttal/e4l_listing.py --long --workers 4     # long-cap subsample
"""

from __future__ import annotations

import argparse
import csv
import json
import multiprocessing as mp
import signal
import time

from common import OUT, UNREACHED, hybrid, is_gac_valid_mpdag, k_of
from oset import optimal_set_closed_form

from bkrobust.demo.graph import MPDAG, undirected_components
from bkrobust.demo.meek import apply_orientations
from bkrobust.gac.dag_level import is_gac_valid_dag

LIST_FULL_MAX = 10**6      # list to completion when |[C^]| <= this
LIST_CAP_S = 300.0         # cap for a full listing
LIST_PROBE_S = 20.0        # throughput probe on larger classes
BNB_CAP_S = 120.0
COUNT_CAP_S = 60.0
UNIT_KILL_S = 900          # whole unit (hybrid included)
ODIR = OUT / "e4l"


class _Timeout(Exception):
    pass


class _Kill(Exception):
    pass


def _alarm(*_):
    raise _Kill


# ------------------------------------------------------------------ class size


def _count_uccg(nodes, edges, memo, deadline):
    if not edges:
        return 1
    key = (nodes, edges)
    if key in memo:
        return memo[key]
    if time.perf_counter() > deadline:
        raise _Timeout
    g = MPDAG(nodes, [], edges)
    total = 0
    for v in sorted(nodes):
        h = apply_orientations(g, [(v, u) for u in sorted(g.adjacent(v))])
        prod = 1
        for comp in undirected_components(h):
            ce = frozenset(e for e in h.undirected_edges if e[0] in comp and e[1] in comp)
            prod *= _count_uccg(frozenset(comp), ce, memo, deadline)
        total += prod
    memo[key] = total
    return total


def class_size(cpdag, limit_s):
    deadline = time.perf_counter() + limit_s
    memo = {}
    size = 1
    try:
        for comp in undirected_components(cpdag):
            ce = frozenset(e for e in cpdag.undirected_edges if e[0] in comp and e[1] in comp)
            size *= _count_uccg(frozenset(comp), ce, memo, deadline)
    except _Timeout:
        return None
    return size


# ------------------------------------------------------------------ listing


def iter_extensions(g):
    if not g.undirected_edges:
        yield g
        return
    a, b = min(g.undirected_edges)
    for o in ((a, b), (b, a)):
        h = apply_orientations(g, [o])
        if h is None:
            raise RuntimeError(f"dead branch {o}")
        yield from iter_extensions(h)


def leg_list(cpdag, x, y, z, k0, limit_s):
    t0 = time.perf_counter()
    best, n, stop = None, 0, None
    for d in iter_extensions(cpdag):
        n += 1
        if not is_gac_valid_dag(d, x, y, z):
            c = len(k0 - d.directed_edges)
            best = c if best is None else min(best, c)
        if time.perf_counter() - t0 > limit_s:
            stop = "time"
            break
    sec = time.perf_counter() - t0
    return {"radius": None if best is None else best, "n_dags": n, "complete": stop is None,
            "seconds": round(sec, 4), "dags_per_s": n / sec if sec > 0 else None}


def leg_bnb(cpdag, x, y, z, k0, limit_s):
    k_rev = frozenset((b, a) for a, b in k0)
    st = {"best": None, "leaves": 0, "nodes": 0, "stop": None}
    t0 = time.perf_counter()

    def rec(g):
        if st["stop"]:
            return
        st["nodes"] += 1
        if time.perf_counter() - t0 > limit_s:
            st["stop"] = "time"
            return
        lb = len(k_rev & g.directed_edges)
        if st["best"] is not None and lb >= st["best"]:
            return
        if not g.undirected_edges:
            st["leaves"] += 1
            if not is_gac_valid_dag(g, x, y, z):
                st["best"] = lb
            return
        a, b = min(g.undirected_edges)
        for o in (((b, a), (a, b)) if (b, a) in k0 else ((a, b), (b, a))):
            rec(apply_orientations(g, [o]))

    rec(cpdag)
    return {"radius": st["best"], "leaves": st["leaves"], "tree_nodes": st["nodes"],
            "complete": st["stop"] is None, "seconds": round(time.perf_counter() - t0, 4)}


# ------------------------------------------------------------------ populations


def units_e2():
    for i, r in enumerate(csv.DictReader(open(OUT / "e2" / "timing.csv"))):
        if r["method"] in ("no_query", "") or not r["x"]:
            continue
        yield {"uid": f"e2:{i}", "pop": "P-E2", "n": int(r["n"]), "p": float(r["p"]),
               "seed": int(r["seed"]), "half": r["knowledge"] == "half", "family": r["family"],
               "x": r["x"], "y": r["y"]}


def units_rel():
    from release_rows import knowledge, rows
    for r in rows(set(knowledge())):
        if r["status"] == "ok":
            yield {"uid": f"rel:{r['arm']}:{r['network']}:{r['x']}:{r['y']}", "pop": "P-REL",
                   "arm": r["arm"], "network": r["network"], "x": r["x"], "y": r["y"]}


def units_big():
    for r in csv.DictReader(open(OUT / "e12" / "large_networks.csv")):
        if r["status"] == "ok":
            yield {"uid": f"big:{r['condition']}:{r['network']}:{r['x']}:{r['y']}", "pop": "P-BIG",
                   "condition": r["condition"], "network": r["network"], "x": r["x"], "y": r["y"]}


def units_syn():
    from e15_aid import FILES
    for which, p in FILES.items():
        for r in csv.DictReader(open(p)):
            yield {"uid": f"syn:{r['instance_id']}", "pop": "P-SYN", "which": which, "row": r}


def build(u):
    """(cpdag, g0, x, y) for a unit."""
    if u["pop"] == "P-E2":
        import e2_timing as e2
        _dag, cpdag, g0, _k = e2.instance(u["n"], u["p"], u["seed"], u["half"])
        return cpdag, g0, u["x"], u["y"]
    if u["pop"] == "P-REL":
        from release_rows import state
        _dag, cpdag, _k, g0, _ = state(u["arm"], u["network"])
        return cpdag, g0, u["x"], u["y"]
    if u["pop"] == "P-BIG":
        from bkrobust.robustness import llm_survival as ls
        from bkrobust.robustness import real_survival as rs
        dag, cpdag = rs.load_network(u["network"])
        bundle, _ = ls.load_bundle()
        k, _g0, _why, _ = ls.build_state(bundle, u["condition"], u["network"], dag, cpdag)
        return cpdag, apply_orientations(cpdag, list(k)), u["x"], u["y"]
    from e15_aid import rebuild
    inst, _ = rebuild(u["row"], u["which"])
    return inst["cpdag"], inst["g0"], inst["x"], inst["y"]


def run_unit(args):
    u, long_cap = args
    signal.signal(signal.SIGALRM, _alarm)
    out = {k: v for k, v in u.items() if k != "row"}
    t_unit = time.perf_counter()
    signal.alarm(int(long_cap + 600) if long_cap else UNIT_KILL_S)
    try:
        cpdag, g0, x, y = build(u)
        z, st = optimal_set_closed_form(g0, x, y)
        if st != "ok" or not is_gac_valid_mpdag(g0, x, y, z):
            signal.alarm(0)
            return dict(out, status=f"no_O:{st}")
        k0 = k_of(cpdag, g0)
        out.update(n_nodes=len(cpdag.nodes), k_undirected=len(cpdag.undirected_edges),
                   k_g0=len(k0), z_size=len(z),
                   largest_component=max((len(c) for c in undirected_components(cpdag)), default=1))
        t = time.perf_counter()
        size = class_size(cpdag, COUNT_CAP_S)
        out["class_size"] = str(size) if size is not None else None
        out["class_size_log10"] = (len(str(size)) - 1 + float("0." + str(size)[:6]) if size else None)
        out["count_s"] = round(time.perf_counter() - t, 4)
        t = time.perf_counter()
        h = hybrid(cpdag, g0, x, y, z, time_limit_s=120)
        out.update(hyb_s=round(time.perf_counter() - t, 6),
                   hyb_r=None if h.radius == UNREACHED else h.radius,
                   hyb_method=h.method, hyb_exact=h.exact)
        if long_cap:
            lim = long_cap
        else:
            lim = LIST_CAP_S if (size is not None and size <= LIST_FULL_MAX) else LIST_PROBE_S
        out["list"] = leg_list(cpdag, x, y, z, k0, lim)
        out["bnb"] = leg_bnb(cpdag, x, y, z, k0, BNB_CAP_S) if not long_cap else None
        signal.alarm(0)
        lst = out["list"]
        if lst["dags_per_s"] and size:
            out["list_projected_s"] = size / lst["dags_per_s"]
        # consistency
        bad = []
        for leg in ("list", "bnb"):
            L = out.get(leg)
            if L and L["complete"] and h.exact and L["radius"] != out["hyb_r"]:
                bad.append(f"{leg}_radius")
        if lst["complete"] and size is not None and lst["n_dags"] != size:
            bad.append("list_count")
        if (not lst["complete"]) and lst["radius"] is not None and out["hyb_r"] is not None \
                and lst["radius"] < out["hyb_r"]:
            bad.append("list_upper_bound")
        out["consistency_fail"] = bad
        out["status"] = "ok"
    except _Kill:
        out["status"] = "killed"
    except Exception as e:  # noqa: BLE001
        signal.alarm(0)
        out["status"] = f"error:{type(e).__name__}:{e}"[:300]
    out["unit_s"] = round(time.perf_counter() - t_unit, 3)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--long", action="store_true")
    ap.add_argument("--pops", default="P-SYN,P-REL,P-BIG,P-E2")
    a = ap.parse_args()
    ODIR.mkdir(parents=True, exist_ok=True)
    if a.long:
        rows = [json.loads(line) for line in open(ODIR / "rows.jsonl")]
        cand = [r for r in rows if r.get("status") == "ok" and r.get("class_size_log10")
                and 5 <= r["class_size_log10"] < 8 and not r["list"]["complete"]]
        import numpy as np
        rng = np.random.default_rng(0)
        pick = []
        for dec in (5, 6, 7):
            c = [r for r in cand if int(r["class_size_log10"]) == dec]
            pick += [c[i] for i in rng.permutation(len(c))[:10]]
        allu = {u["uid"]: u for g in (units_e2, units_rel, units_big, units_syn) for u in g()}
        jobs = [(allu[r["uid"]], 1800.0) for r in pick]
        path = ODIR / "rows_long.jsonl"
    else:
        gens = {"P-E2": units_e2, "P-REL": units_rel, "P-BIG": units_big, "P-SYN": units_syn}
        jobs = []
        for p in a.pops.split(","):
            us = list(gens[p]())
            if a.smoke:
                us = us[:: max(1, len(us) // 6)][:6]
            jobs += [(u, None) for u in us]
        path = ODIR / ("rows_smoke.jsonl" if a.smoke else "rows.jsonl")
    done = set()
    if path.exists():
        done = {json.loads(line)["uid"] for line in open(path)}
    jobs = [j for j in jobs if j[0]["uid"] not in done]
    # big populations first so the long tail starts early
    order = {"P-BIG": 0, "P-E2": 1, "P-REL": 2, "P-SYN": 3}
    jobs.sort(key=lambda j: (order[j[0]["pop"]], -j[0].get("n", 0)))
    print(len(jobs), "units to run;", len(done), "already done", flush=True)
    t0 = time.time()
    with open(path, "a") as f, mp.Pool(a.workers, maxtasksperchild=50) as P:
        for i, o in enumerate(P.imap_unordered(run_unit, jobs)):
            f.write(json.dumps(o, default=str) + "\n")
            f.flush()
            if i % 50 == 0:
                print(f"{i} {time.time() - t0:.0f}s {o['pop']} {o.get('status')} "
                      f"size={o.get('class_size_log10')} unit_s={o.get('unit_s')}", flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
