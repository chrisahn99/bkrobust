"""E3: depth budget b x SAT mode (ladder / minimise) on the 144 high-k instances.

Run from the repository root:

    .venv/bin/python experiments/rebuttal/E3_depth_budget_and_sat_mode/run.py --smoke
    .venv/bin/python experiments/rebuttal/E3_depth_budget_and_sat_mode/run.py

One JSON line per instance (all 12 arms plus a reference call to the submitted
``breakdown_radius``), flushed as written; a rerun skips instances already in the
file. Exit 2 if any instance has arms that disagree on the radius or a SAT
witness that does not fail on replay.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import multiprocessing as mp
import platform
import subprocess
import sys
import threading
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402

from bkrobust.core import oracle as _oracle  # noqa: E402
from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover  # noqa: E402
from bkrobust.demo.graph import MPDAG, undirected_components  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.gac import mpdag_level as _gac  # noqa: E402
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag  # noqa: E402
from bkrobust.hybrid import DEFAULT_SEARCH_BUDGET, breakdown_radius  # noqa: E402
from bkrobust.sat.e1 import radius_e1  # noqa: E402
from bkrobust.sat.verify import to_mpdag  # noqa: E402
from bkrobust.search.exact import SearchStats, retractable_edges  # noqa: E402
from bkrobust.search.exact_fast import radius_local_up_fast  # noqa: E402
from bkrobust.synth.generators import erdos_renyi_dag  # noqa: E402
from sat_minimise import radius_minimise  # noqa: E402

CRITERION = "gac"
BUDGETS = (0, 1, 2, 3, 4, 6)
MODES = ("ladder", "minimise")
ARMS = [(b, m) for b in BUDGETS for m in MODES]
SOLVER_SEED = 0


class ArmCap(Exception):
    """The search leg of one arm ran past its wall-clock cap."""


def repo_commit() -> dict:
    def git(*a):
        return subprocess.run(["git", "-C", str(ROOT), *a], capture_output=True, text=True, timeout=30)

    try:
        head = git("rev-parse", "HEAD").stdout.strip() or "unknown"
        dirty = bool(git("status", "--porcelain", "--", "src").stdout.strip())
    except Exception:  # no git on the node: record it rather than fail
        head, dirty = "unknown", None
    return {"commit": head, "src_dirty": dirty}


def _unpadded(dag: MPDAG) -> MPDAG:
    # highk.jsonl names nodes V0..V17; the current generator pads (V00..V17).
    # Relabelling changes no edge, only the names the recorded (x, y) refer to.
    ren = {v: f"V{int(v[1:])}" for v in dag.nodes}
    return MPDAG([ren[v] for v in dag.nodes], directed=[(ren[a], ren[b]) for a, b in dag.directed_edges])


def build_instance(row: dict):
    rng = np.random.default_rng(int(row["seed"]))
    dag = _unpadded(erdos_renyi_dag(int(row["n"]), rng, float(row["edge_prob"])))
    cpdag = dag_to_cpdag(dag)
    g0 = apply_orientations(cpdag, sorted(knowledge_to_recover(dag, cpdag)))
    if g0 is None:
        raise RuntimeError("knowledge did not close")
    comps = undirected_components(cpdag)
    largest = max(comps, key=len) if comps else frozenset()
    c = len(largest)
    e_in = sum(1 for a, b in cpdag.undirected_edges if a in largest and b in largest)
    stats = {
        "m_undirected": len(cpdag.undirected_edges),
        "k_g0": len(retractable_edges(cpdag, g0)),
        "largest_component": c,
        "largest_component_density": round(e_in / (c * (c - 1) / 2), 4) if c > 1 else 0.0,
    }
    return cpdag, g0, stats


def _good_query(g0, x, y):
    if x not in g0.nodes or y not in g0.nodes:
        return None
    o = optimal_adjustment_set_mpdag(g0, x, y)
    if not o:
        return None
    z = frozenset(o)
    return z if is_gac_valid_mpdag(g0, x, y, z) else None


def choose_query(g0, row):
    z = _good_query(g0, row["x"], row["y"])
    if z is not None:
        return row["x"], row["y"], z, "recorded"
    for x, y in itertools.permutations(sorted(g0.nodes), 2):
        z = _good_query(g0, x, y)
        if z is not None:
            return x, y, z, "fallback"
    return None


def _clear_caches():
    # Every arm starts cold, so no arm is timed on closures another arm paid for.
    _gac.clear_cache()
    _oracle.clear_cache()


def _check_witness(cpdag, g0, x, y, z, witness, radius):
    w = to_mpdag(cpdag, frozenset(tuple(e) for e in witness))
    fails = not is_gac_valid_mpdag(w, x, y, z)
    kept = set(w.directed_edges)
    retractions = sum(1 for e in retractable_edges(cpdag, g0) if e not in kept)
    return {"witness_fails": fails, "witness_retractions": retractions,
            "witness_ok": bool(fails and retractions == radius)}


def run_arm(cpdag, g0, x, y, z, b, mode, *, sat_limit, cap):
    _clear_caches()
    deadline = [float("inf")]  # armed after the top-state check: the cap bounds the search leg only

    def fails(g):
        if time.perf_counter() > deadline[0]:
            raise ArmCap
        return not is_gac_valid_mpdag(g, x, y, z)

    rec = {"b": b, "mode": mode, "radius": None, "exact": True, "capped": False,
           "method": None, "top_s": 0.0, "search_s": 0.0, "sat_s": 0.0, "sat_build_s": 0.0}
    t0 = time.perf_counter()
    top_valid = not fails(cpdag)
    rec["top_s"] = time.perf_counter() - t0
    if top_valid:
        rec.update(radius=UNREACHED, method="top_state")
    else:
        st = SearchStats()
        t0 = time.perf_counter()
        deadline[0] = t0 + cap
        try:
            found = radius_local_up_fast(cpdag, g0, fails, max_depth=b, stats=st)
        except ArmCap:
            found = None
        rec["search_s"] = time.perf_counter() - t0
        rec["search_stats"] = st.as_dict()
        if found is None:
            rec.update(capped=True, exact=False, method="capped")
        elif found.radius != UNREACHED:
            rec.update(radius=found.radius, method="degenerate" if found.radius == 0 else "local_up_fast")
        elif found.exact:
            rec.update(radius=UNREACHED, method="local_up_fast")
        else:
            t0 = time.perf_counter()
            if mode == "ladder":
                lad = radius_e1(cpdag, g0, x, y, z, time_limit_s=sat_limit, seed=SOLVER_SEED, criterion=CRITERION)
                rec.update(radius=lad.radius, exact=lad.exact, method="e1_ladder", sat_build_s=lad.build_seconds,
                           rungs=len(lad.rung_seconds), rung_status={str(k): v for k, v in lad.rung_status.items()},
                           conflicts=lad.conflicts, branches=lad.branches)
                witness = lad.witness_orientations
            else:
                mn = radius_minimise(cpdag, g0, x, y, z, time_limit_s=sat_limit, seed=SOLVER_SEED, criterion=CRITERION)
                rec.update(radius=mn["radius"], exact=mn["exact"], method="e1_minimise", sat_build_s=mn["build_seconds"],
                           sat_status=mn["status"], conflicts=mn["conflicts"], branches=mn["branches"])
                witness = mn["witness_orientations"]
            rec["sat_s"] = time.perf_counter() - t0
            if rec["exact"] and rec["radius"] not in (UNREACHED, -2):
                rec.update(_check_witness(cpdag, g0, x, y, z, witness, rec["radius"]))
    rec["total_s"] = rec["top_s"] + rec["search_s"] + rec["sat_s"]
    return rec


def process_unit(task):
    idx, row, cfg = task
    out = {"unit": idx, "criterion": CRITERION, **cfg["commit"], "smoke": cfg["smoke"],
           "instance": {k: row[k] for k in ("n", "edge_prob", "seed", "x", "y", "z_size")},
           "hunt_e1_radius_backdoor": row.get("e1", {}).get("radius")}
    t_unit = time.perf_counter()
    try:
        cpdag, g0, stats = build_instance(row)
        q = choose_query(g0, row)
    except Exception as exc:  # recorded, never dropped
        out.update(error=f"build: {exc!r}", consistent=None)
        return out
    out["rebuilt"] = stats
    out["rebuilt_matches_file"] = {k: stats[k] == row.get(k) for k in stats}
    if q is None:
        out.update(error="no valid query", consistent=None)
        return out
    x, y, z, src = q
    out["query"] = {"x": x, "y": y, "z": sorted(z), "z_size": len(z), "source": src,
                    "z_size_matches_file": len(z) == row.get("z_size")}

    arms = cfg["arms"]
    rot = idx % len(arms)  # rotate arm order so no arm always runs first or last
    order = arms[rot:] + arms[:rot]
    recs = []
    for pos, (b, mode) in enumerate(order):
        r = run_arm(cpdag, g0, x, y, z, b, mode, sat_limit=cfg["sat_limit"], cap=cfg["cap"])
        r["order"] = pos
        recs.append(r)
    recs.sort(key=lambda r: (r["b"], MODES.index(r["mode"])))
    out["arms"] = recs

    failures = []
    if cfg["reference"]:
        _clear_caches()
        t0 = time.perf_counter()
        ref = breakdown_radius(cpdag, None, x, y, z, g0=g0, search_budget=DEFAULT_SEARCH_BUDGET,
                               criterion=CRITERION, time_limit_s=cfg["sat_limit"])
        out["reference"] = {"b": DEFAULT_SEARCH_BUDGET, "radius": ref.radius, "exact": ref.exact,
                            "method": ref.method, "oracle": ref.oracle, "total_s": ref.total_seconds,
                            "wall_s": time.perf_counter() - t0}
        if ref.exact:
            recs_for_check = recs + [{"b": "ref", "mode": "library", "radius": ref.radius, "exact": True}]
        else:
            recs_for_check = recs
    else:
        recs_for_check = recs

    radii = {r["radius"] for r in recs_for_check if r["exact"] and not r.get("capped")}
    if len(radii) > 1:
        failures.append("radius disagreement: " + ", ".join(
            f"b={r['b']}/{r['mode']}={r['radius']}" for r in recs_for_check if r["exact"]))
    for r in recs:
        if "witness_ok" in r and not r["witness_ok"]:
            failures.append(f"witness b={r['b']}/{r['mode']}: fails={r['witness_fails']} "
                            f"retractions={r['witness_retractions']} radius={r['radius']}")
    out["radius"] = radii.pop() if len(radii) == 1 else None
    out["n_capped"] = sum(r["capped"] for r in recs)
    out["n_inexact"] = sum((not r["exact"]) and not r["capped"] for r in recs)
    out["failures"] = failures
    out["consistent"] = not failures

    if cfg["smoke"] and cfg.get("cap_probe_unit") == idx:
        # Exercise the cap path: a zero cap must be reported as capped, not as a radius.
        p = run_arm(cpdag, g0, x, y, z, max(BUDGETS), "ladder", sat_limit=cfg["sat_limit"], cap=0.0)
        out["cap_probe"] = {"capped": p["capped"], "method": p["method"], "radius": p["radius"]}
        if p["method"] != "top_state" and not p["capped"]:
            out["failures"].append("cap probe was not capped")
            out["consistent"] = False
    out["unit_s"] = time.perf_counter() - t_unit
    return out


def smoke_units(rows):
    # Smallest instances first; one with no failure, one at radius 1 and one beyond,
    # judged from the file's (back-door) radius, so the top-state, search and SAT
    # paths are all reached at some b.
    small = [(i, r) for i, r in enumerate(rows) if r["n"] == 8]
    pick = []
    for want in (lambda v: v == -1, lambda v: v == 1, lambda v: v >= 2):
        for i, r in small:
            if want(r["e1"]["radius"]) and i not in pick:
                pick.append(i)
                break
    return pick


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--instances", default=str(ROOT / "results/axisb3/highk.jsonl"),
                    help="high-k instance file (default: results/axisb3/highk.jsonl of this repository)")
    ap.add_argument("--out", default=str(ROOT / "results/rebuttal/E3"))
    ap.add_argument("--smoke", action="store_true", help="3 small instances, all arms, plus a cap probe")
    ap.add_argument("--workers", type=int, default=1, help="1 = serial, as registered (timings)")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--sat-time-limit", type=float, default=300.0, help="CP-SAT seconds per call (per rung)")
    ap.add_argument("--arm-cap", type=float, default=600.0, help="wall-clock cap on one arm's search leg")
    ap.add_argument("--no-reference", action="store_true", help="skip the library breakdown_radius check")
    args = ap.parse_args()
    if not 0 <= args.shard < args.n_shards:
        ap.error("need 0 <= --shard < --n-shards")

    inst_path = Path(args.instances)
    raw = inst_path.read_bytes()
    rows = [json.loads(line) for line in raw.decode().splitlines() if line.strip()]
    unit_ids = smoke_units(rows) if args.smoke else list(range(len(rows)))
    unit_ids = [i for i in unit_ids if i % args.n_shards == args.shard]

    out_dir = Path(args.out) / ("smoke" if args.smoke else "full")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"rows_shard{args.shard}of{args.n_shards}.jsonl"

    import ortools

    commit = repo_commit()
    cfg_record = {
        "experiment": "E3", "criterion": CRITERION, "budgets": list(BUDGETS), "modes": list(MODES),
        "top_state_check": True, "sat_time_limit_s": args.sat_time_limit, "arm_cap_s": args.arm_cap,
        "solver": {"num_workers": 1, "random_seed": SOLVER_SEED}, "reference_call": not args.no_reference,
        "reference_budget": DEFAULT_SEARCH_BUDGET, "workers": args.workers, "smoke": args.smoke,
        "instances": str(inst_path), "instances_sha256": hashlib.sha256(raw).hexdigest(),
        "n_instances_in_file": len(rows), **commit, "ortools": ortools.__version__,
        "python": platform.python_version(), "host": platform.node(),
    }
    (out_dir / f"config_shard{args.shard}of{args.n_shards}.json").write_text(json.dumps(cfg_record, indent=1))

    done = set()
    if out_file.exists():
        for line in out_file.read_text().splitlines():
            try:
                done.add(json.loads(line)["unit"])
            except (json.JSONDecodeError, KeyError):
                pass  # a line cut by a kill; that unit is simply redone
    todo = [i for i in unit_ids if i not in done]
    print(f"[E3] {len(unit_ids)} units in shard, {len(done)} already written, {len(todo)} to run -> {out_file}",
          flush=True)

    cfg = {"arms": ARMS, "sat_limit": args.sat_time_limit, "cap": args.arm_cap, "reference": not args.no_reference,
           "smoke": args.smoke, "commit": commit,
           "cap_probe_unit": unit_ids[-1] if (args.smoke and unit_ids) else None}
    tasks = [(i, rows[i], cfg) for i in todo]

    state = {"n": 0, "t0": time.time(), "stop": False}

    def heartbeat():
        while not state["stop"]:
            time.sleep(30)
            if not state["stop"]:
                print(f"[E3] alive: {state['n']}/{len(todo)} written, {time.time() - state['t0']:.0f}s elapsed",
                      flush=True)

    threading.Thread(target=heartbeat, daemon=True).start()
    with out_file.open("a") as fh:
        if args.workers > 1:
            pool = mp.Pool(args.workers)
            results = pool.imap_unordered(process_unit, tasks)
        else:
            pool = None
            results = map(process_unit, tasks)
        for res in results:
            fh.write(json.dumps(res) + "\n")
            fh.flush()
            state["n"] += 1
            flag = "OK" if res.get("consistent") else ("ERR" if res.get("consistent") is None else "INCONSISTENT")
            print(f"[E3] unit {res['unit']} r={res.get('radius')} {flag} {res.get('unit_s', 0):.2f}s "
                  f"({state['n']}/{len(todo)})", flush=True)
        if pool is not None:
            pool.close()
            pool.join()
    state["stop"] = True

    bad = []
    for line in out_file.read_text().splitlines():
        try:
            r = json.loads(line)
        except json.JSONDecodeError:
            continue
        if r.get("consistent") is False:
            bad.append((r["unit"], r["failures"]))
    for u, f in bad:
        print(f"[E3] CONSISTENCY FAILURE unit {u}: {f}", flush=True)
    sys.exit(2 if bad else 0)


if __name__ == "__main__":
    main()
