"""E2: running time of the submitted hybrid under the GAC, by graph size (n = 8 to 400).

Each unit is one (family, parameter, n, seed) instance. Every timed stage runs in a fresh
child process of this same script (``--leg``), so that caches are cold, a stage can be killed
at the cap, and its memory can be limited. The parent only schedules, checks consistency and
writes one JSON line per unit. See PREREG.md for the protocol.

    .venv/bin/python experiments/rebuttal/E2_running_time_by_graph_size/run.py --smoke
    .venv/bin/python experiments/rebuttal/E2_running_time_by_graph_size/run.py --out results/rebuttal/E2
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import platform
import socket
import statistics
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SRC = ROOT / "src"

CRITERION = "gac"
SEARCH_BUDGET = 3  # bkrobust.hybrid.DEFAULT_SEARCH_BUDGET, the submitted default
GATE_MAX_ATTEMPTS = 100
COLD_REPEATS = 5

APPD_N = (8, 10, 12, 14, 16, 18, 20, 24)
APPD_P = (0.3, 0.5, 0.7, 0.85)
SPARSE_N = (16, 24, 32, 48, 64, 100, 150, 200, 300, 400)
SPARSE_DEG = (2, 4)
SEEDS = tuple(range(8))

SKIP_MIN_ATTEMPTS = 4
SKIP_CENSORED_SHARE = 0.75
CENSORED = ("timeout", "mem", "over_cap", "skipped")


# --------------------------------------------------------------------------- the grid


def unit_grid(smoke: bool) -> list[dict]:
    """Units sorted by n first, so the per-series skip rule sees smaller n before larger n."""
    units = []
    if smoke:
        cells = [("appD", 0.5, 8, 0), ("appD", 0.85, 10, 1), ("sparse", 2, 16, 0)]
    else:
        cells = [("appD", p, n, s) for n in APPD_N for p in APPD_P for s in SEEDS]
        cells += [("sparse", d, n, s) for n in SPARSE_N for d in SPARSE_DEG for s in SEEDS]
    for fam, param, n, seed in cells:
        edge_prob = float(param) if fam == "appD" else min(1.0, param / (n - 1))
        units.append({
            "uid": f"{fam}_{param:g}_n{n}_s{seed}",
            "family": fam, "param": param, "n": n, "seed": seed, "edge_prob": edge_prob,
        })
    units.sort(key=lambda u: (u["n"], u["family"], u["param"], u["seed"]))
    return units


# --------------------------------------------------------------------------- child legs
# Everything below up to the parent section runs inside a child process.


def _graphs(p: dict):
    from bkrobust.demo.graph import MPDAG

    cpdag = MPDAG(p["nodes"], directed=[tuple(e) for e in p["cpdag_directed"]],
                  undirected=[tuple(e) for e in p["cpdag_undirected"]])
    g0 = MPDAG(p["nodes"], directed=[tuple(e) for e in p["g0_directed"]])
    return cpdag, g0, p["x"], p["y"], frozenset(p["z"])


def _fast_query(dag, is_valid):
    """First sorted ordered pair whose optimal set (HPM 2022) in the DAG is non-empty and valid.

    Same definition as bkrobust.demo.evaluate.optimal_adjustment_set_dag, with descendants and
    ancestors computed once per node: at n = 400 the pair loop visits up to 160k pairs.
    """
    nodes = sorted(dag.nodes)
    desc = {v: dag.descendants(v) for v in nodes}
    anc = {v: dag.ancestors(v) for v in nodes}
    for x, y in itertools.permutations(nodes, 2):
        if y not in desc[x]:
            continue
        cn = desc[x] & (anc[y] | {y})
        o = set()
        for c in cn:
            o |= dag.parents(c)
        o -= cn | {x}
        if o and is_valid(x, y, frozenset(o)):
            return x, y, frozenset(o)
    return None


def leg_prep(p: dict) -> dict:
    """Draw the instance (with the k >= 1 gate), choose the query, time the cheap stages."""
    import networkx as nx
    import numpy as np

    from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover
    from bkrobust.demo.meek import apply_orientations
    from bkrobust.gac.mpdag_level import clear_cache, is_gac_valid_mpdag, why_invalid_gac
    from bkrobust.search.exact import retractable_edges
    from bkrobust.synth.generators import erdos_renyi_dag

    n, seed, edge_prob = p["n"], p["seed"], p["edge_prob"]
    for attempt in range(GATE_MAX_ATTEMPTS + 1):
        consistency: list[str] = []
        eff_seed = seed + 10_000 * attempt
        t0 = time.perf_counter()
        dag = erdos_renyi_dag(n, np.random.default_rng(eff_seed), edge_prob)
        t1 = time.perf_counter()
        cpdag = dag_to_cpdag(dag)
        cpdag_s = time.perf_counter() - t1
        dag_s = t1 - t0
        if not cpdag.undirected_edges:
            continue
        g0 = dag
        g0_s = None
        if n <= 24:
            # Appendix D's construction, kept where it is cheap, to check it yields the DAG.
            t2 = time.perf_counter()
            built = apply_orientations(cpdag, sorted(knowledge_to_recover(dag, cpdag)))
            g0_s = time.perf_counter() - t2
            if built != dag:
                consistency.append("appD_g0_differs_from_dag")
        t3 = time.perf_counter()
        q = _fast_query(g0, lambda a, b, c: is_gac_valid_mpdag(g0, a, b, c))
        query_s = time.perf_counter() - t3
        if q is None:
            continue
        break
    else:
        return {"status": "gate_exhausted"}

    x, y, z = q
    if n <= 24:
        from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag

        if optimal_adjustment_set_mpdag(g0, x, y) != set(z):
            consistency.append("fast_optimal_set_differs")

    def cold(g):
        times, verdict = [], None
        for _ in range(COLD_REPEATS):
            clear_cache()
            t = time.perf_counter()
            verdict = is_gac_valid_mpdag(g, x, y, z)
            times.append(time.perf_counter() - t)
        return verdict, times

    top_valid, top_t = cold(cpdag)
    g0_valid, g0_t = cold(g0)
    clear_cache()
    if not g0_valid:
        consistency.append("query_invalid_at_g0")

    und = nx.Graph()
    und.add_nodes_from(cpdag.nodes)
    und.add_edges_from(cpdag.undirected_edges)
    comps = [len(c) for c in nx.connected_components(und) if len(c) > 1]
    return {
        "status": "ok", "attempt": attempt, "effective_seed": eff_seed,
        "n_edges": len(dag.directed_edges),
        "k_undirected": len(cpdag.undirected_edges),
        "k_g0": len(retractable_edges(cpdag, g0)),
        "largest_component": max(comps) if comps else 0,
        "x": x, "y": y, "z": sorted(z),
        "dag_s": dag_s, "cpdag_s": cpdag_s, "appD_g0_s": g0_s, "query_s": query_s,
        "top_valid": top_valid, "top_reason": why_invalid_gac(cpdag, x, y, z),
        "top_s": top_t, "top_median_s": statistics.median(top_t),
        "g0_valid": g0_valid, "g0_call_s": g0_t, "g0_call_median_s": statistics.median(g0_t),
        "consistency": consistency,
        "graph": {
            "nodes": list(cpdag.nodes),
            "cpdag_directed": sorted(cpdag.directed_edges),
            "cpdag_undirected": sorted(cpdag.undirected_edges),
            "g0_directed": sorted(g0.directed_edges),
            "x": x, "y": y, "z": sorted(z),
        },
    }


def _search(p: dict, max_depth) -> dict:
    from bkrobust.gac.mpdag_level import clear_cache, is_gac_valid_mpdag
    from bkrobust.search.exact import SearchStats
    from bkrobust.search.exact_fast import radius_local_up_fast

    cpdag, g0, x, y, z = _graphs(p)
    clear_cache()
    st = SearchStats()
    t = time.perf_counter()
    res = radius_local_up_fast(cpdag, g0, lambda g: not is_gac_valid_mpdag(g, x, y, z),
                               max_depth=max_depth, stats=st)
    return {"seconds": time.perf_counter() - t, "radius": res.radius, "exact": res.exact,
            "method": res.method, **st.as_dict()}


def leg_search(p: dict) -> dict:
    return _search(p, SEARCH_BUDGET)


def leg_search_full(p: dict) -> dict:
    return _search(p, None)


def leg_sat(p: dict) -> dict:
    """SAT build always; the ladder alone (no search first) when the top check failed.

    The build line is printed and flushed before the ladder starts, so a ladder killed at the
    cap still leaves its build time in the partial output.
    """
    from bkrobust.sat.e1 import _build, radius_e1

    cpdag, g0, x, y, z = _graphs(p)
    t = time.perf_counter()
    model, _ov, retracted, _d = _build(cpdag, g0, x, y, z, criterion=CRITERION)
    build_s = time.perf_counter() - t
    proto = model.Proto()
    build = {"build_seconds": build_s, "n_variables": len(proto.variables),
             "n_constraints": len(proto.constraints), "n_retractable": len(retracted)}
    print("PARTIAL " + json.dumps(build), flush=True)
    del model, proto
    if not p["ladder"]:
        return {"mode": "build_only", **build}
    lad = radius_e1(cpdag, g0, x, y, z, time_limit_s=p["cap_s"], seed=0, criterion=CRITERION)
    return {"mode": "ladder", **build, "radius": lad.radius, "exact": lad.exact,
            "ladder_build_seconds": lad.build_seconds, "ladder_total_seconds": lad.total_seconds,
            "rung_seconds": {str(k): v for k, v in lad.rung_seconds.items()},
            "rung_status": {str(k): v for k, v in lad.rung_status.items()}}


def _hybrid(p: dict, budget: int) -> dict:
    from bkrobust.gac.mpdag_level import clear_cache
    from bkrobust.hybrid import breakdown_radius

    cpdag, g0, x, y, z = _graphs(p)
    clear_cache()
    t = time.perf_counter()
    r = breakdown_radius(cpdag, None, x, y, z, g0=g0, search_budget=budget,
                         criterion=CRITERION, time_limit_s=p["cap_s"])
    wall = time.perf_counter() - t
    return {"seconds": wall, "radius": r.radius, "method": r.method, "oracle": r.oracle,
            "exact": r.exact, "top_seconds": r.top_seconds, "search_seconds": r.search_seconds,
            "ladder_seconds": r.ladder_seconds, "total_seconds": r.total_seconds,
            **r.stats.as_dict()}


def leg_hybrid(p: dict) -> dict:
    return _hybrid(p, SEARCH_BUDGET)


def leg_hybrid_b0(p: dict) -> dict:
    # Budget 0 sends every finite-radius instance to the SAT stage: exercises that branch.
    return _hybrid(p, 0)


LEGS = {"prep": leg_prep, "search": leg_search, "search_full": leg_search_full,
        "sat": leg_sat, "hybrid": leg_hybrid, "hybrid_b0": leg_hybrid_b0}


def child_main(leg: str) -> None:
    # Set here rather than via preexec_fn, which is unsafe while the parent runs threads.
    gb = float(os.environ.get("E2_MEM_GB", "0"))
    if gb > 0:
        try:
            import resource

            resource.setrlimit(resource.RLIMIT_AS, (int(gb * 2**30), int(gb * 2**30)))
        except (ImportError, ValueError, OSError):
            pass  # macOS may refuse RLIMIT_AS; the cap then bounds time only
    sys.path.insert(0, str(SRC))
    payload = json.loads(sys.stdin.read())
    out = LEGS[leg](payload)
    print("RESULT " + json.dumps(out), flush=True)


# --------------------------------------------------------------------------- parent


def run_leg(leg: str, payload: dict, cap_s: float, args) -> dict:
    """Run one leg in a fresh child; never raises. Censored legs keep any partial output."""
    cmd = [sys.executable, str(Path(__file__).resolve()), "--leg", leg]
    env = dict(os.environ, PYTHONPATH=str(SRC), PYTHONHASHSEED="0", E2_MEM_GB=str(args.mem_gb))
    t = time.perf_counter()
    status, out, err = "ok", "", ""
    try:
        cp = subprocess.run(cmd, input=json.dumps(payload), capture_output=True, text=True,
                            timeout=cap_s + args.startup_slack, env=env)
        out, err = cp.stdout, cp.stderr
        if cp.returncode != 0:
            oom = "MemoryError" in err or "bad_alloc" in err or cp.returncode in (-9, -6, 137)
            status = "mem" if oom else "error"
    except subprocess.TimeoutExpired as e:
        status = "timeout"
        out = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
    rec: dict = {"status": status, "wall_s": time.perf_counter() - t}
    for line in out.splitlines():
        if line.startswith("PARTIAL "):
            rec["partial"] = json.loads(line[8:])
        elif line.startswith("RESULT "):
            rec.update(json.loads(line[7:]))
    measured = rec.get("seconds", rec.get("ladder_total_seconds", rec.get("build_seconds", 0.0)))
    if rec["status"] == "ok" and measured > cap_s:
        rec["status"] = "over_cap"
    if status == "error":
        rec["stderr_tail"] = err.strip().splitlines()[-5:]
    return rec


class State:
    """Rows written so far, shared by the worker threads, for the skip rule and progress."""

    def __init__(self, rows: list[dict], total: int):
        self.lock = threading.Lock()
        self.rows = rows
        self.total = total
        self.done = 0
        self.running: dict[str, str] = {}
        self.failures = 0
        self.started = time.time()

    def should_skip(self, unit: dict, leg: str) -> bool:
        with self.lock:
            by_n: dict[int, list[str]] = {}
            for r in self.rows:
                if (r.get("family"), r.get("param")) != (unit["family"], unit["param"]):
                    continue
                if r.get("n", 0) >= unit["n"] or leg not in r.get("legs", {}):
                    continue
                by_n.setdefault(r["n"], []).append(r["legs"][leg]["status"])
            for sts in by_n.values():
                att = [s for s in sts if s not in ("not_needed", "skipped")]
                cens = [s for s in att if s in CENSORED]
                if len(att) >= SKIP_MIN_ATTEMPTS and len(cens) / len(att) >= SKIP_CENSORED_SHARE:
                    return True
            return False


def check_consistency(row: dict) -> list[str]:
    from bkrobust.core.conventions import UNREACHED

    bad = list(row["prep"].get("consistency", []))
    L = row["legs"]
    ok = {k: v for k, v in L.items() if v.get("status") == "ok"}
    top_valid = row["prep"]["top_valid"]
    h, s, sat, full, b0 = (ok.get(k) for k in ("hybrid", "search", "sat", "search_full", "hybrid_b0"))
    if h:
        if top_valid and not (h["method"] == "top_state" and h["radius"] == UNREACHED):
            bad.append("top_valid_but_hybrid_not_top_state")
        if not top_valid and h["method"] == "top_state":
            bad.append("top_invalid_but_hybrid_top_state")
        if h["oracle"] != "gac_criterion":
            bad.append("hybrid_oracle_not_gac")
    if h and s and s["radius"] != UNREACHED and h["radius"] != s["radius"]:
        bad.append("hybrid_vs_bounded_search")
    if h and s and s["exact"] and s["radius"] == UNREACHED and h["radius"] != UNREACHED:
        bad.append("hybrid_vs_exhausted_search")
    finite = [("bounded_search", s["radius"])] if s and s["radius"] != UNREACHED else []
    exact = []
    for name, leg in (("hybrid", h), ("sat_ladder", sat), ("search_full", full), ("hybrid_b0", b0)):
        if leg and leg.get("exact") and "radius" in leg:
            exact.append((name, leg["radius"]))
    if not top_valid and any(r == UNREACHED for _, r in exact):
        bad.append("top_fails_but_exact_leg_unreached")
    vals = {r for _, r in exact + finite}
    if len(vals) > 1:
        bad.append("radius_disagreement:" + ",".join(f"{a}={b}" for a, b in exact + finite))
    if h and h["method"] == "e1_ladder" and h["exact"] and h["radius"] != UNREACHED \
            and h["radius"] <= SEARCH_BUDGET:
        bad.append("ladder_radius_within_search_budget")
    return bad


def process_unit(unit: dict, args, state: State, meta: dict) -> dict:
    tag = unit["uid"]
    row = {**unit, **meta, "legs": {}}

    def leg(name, payload, cap):
        with state.lock:
            state.running[tag] = name
        if name != "prep" and state.should_skip(unit, name):
            row["legs"][name] = {"status": "skipped"}
            return row["legs"][name]
        rec = run_leg(name, payload, cap, args)
        row["legs"][name] = rec
        return rec

    prep = leg("prep", unit, args.prep_cap)
    graph = prep.pop("graph", None)
    row["prep"] = row["legs"].pop("prep")
    row["status"] = prep["status"]  # ok, gate_exhausted, or the child's failure (timeout, mem, error)
    if prep["status"] != "ok":
        return row
    payload = {**graph, "cap_s": args.cap}
    top_valid = prep["top_valid"]

    if top_valid:
        row["legs"]["search"] = {"status": "not_needed"}
    else:
        leg("search", payload, args.cap)
    leg("sat", {**payload, "ladder": not top_valid}, args.cap)
    h = leg("hybrid", payload, args.cap)
    if args.with_b0:
        leg("hybrid_b0", payload, args.cap)
    sat = row["legs"]["sat"]
    # The unbounded search is the independent check of any radius the SAT stage produced.
    needs_full = h.get("method") == "e1_ladder" or args.with_b0 or (
        sat.get("mode") == "ladder" and sat.get("radius", -1) > SEARCH_BUDGET)
    if needs_full and not top_valid:
        leg("search_full", payload, args.cap)
    else:
        row["legs"]["search_full"] = {"status": "not_needed"}
    row["consistency_failures"] = check_consistency(row)
    return row


def repo_commit(override: str | None) -> tuple[str, bool | None]:
    if override:
        return override, None
    try:
        c = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True,
                           text=True, timeout=10).stdout.strip()
        d = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "src"],
                           capture_output=True, text=True, timeout=10).stdout.strip()
        return c or "unknown", bool(d)
    except (OSError, subprocess.SubprocessError):
        return "unknown", None


def load_rows(out: Path) -> list[dict]:
    rows = []
    for f in sorted(out.glob("rows_*.jsonl")):
        for line in f.read_text().splitlines():
            if line.strip():
                try:
                    rows.append(json.loads(line))
                except json.JSONDecodeError:
                    pass  # a line cut by a kill; the unit is simply rerun
    return rows


def heartbeat(state: State, stop: threading.Event):
    while not stop.wait(30):
        with state.lock:
            el = time.time() - state.started
            print(f"[{time.strftime('%H:%M:%S')}] {state.done}/{state.total} units this run, "
                  f"{el / 60:.1f} min, consistency failures {state.failures}, "
                  f"running {dict(state.running)}", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--smoke", action="store_true", help="3 tiny units, cap 20 s, every branch")
    ap.add_argument("--out", default=None, help="output directory (default results/rebuttal/E2)")
    ap.add_argument("--workers", type=int, default=1, help="concurrent units (serial by default)")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--cap", type=float, default=120.0, help="seconds per leg (Appendix D)")
    ap.add_argument("--prep-cap", type=float, default=900.0, help="seconds for instance + query")
    ap.add_argument("--mem-gb", type=float, default=4.0, help="address-space limit per leg, 0 = none")
    ap.add_argument("--startup-slack", type=float, default=20.0, help="import time allowed per child")
    ap.add_argument("--with-b0", action="store_true",
                    help="also run the hybrid with search budget 0 (always on in --smoke)")
    ap.add_argument("--retry-errors", action="store_true", help="rerun units whose status is error")
    ap.add_argument("--commit", default=None, help="repository commit if .git is absent")
    ap.add_argument("--leg", default=None, help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args.leg:
        child_main(args.leg)
        return 0
    sys.path.insert(0, str(SRC))
    if args.smoke:
        args.cap, args.prep_cap, args.with_b0 = min(args.cap, 20.0), 60.0, True
    out = Path(args.out or ("results/rebuttal/E2_smoke" if args.smoke else "results/rebuttal/E2"))
    if not out.is_absolute():
        out = ROOT / out
    out.mkdir(parents=True, exist_ok=True)

    commit, dirty = repo_commit(args.commit)
    meta = {"criterion": CRITERION, "search_budget": SEARCH_BUDGET, "cap_s": args.cap,
            "mem_gb": args.mem_gb, "commit": commit, "src_dirty": dirty,
            "python": platform.python_version(), "host": socket.gethostname(),
            "smoke": args.smoke}
    config = {**meta, "workers": args.workers, "n_shards": args.n_shards,
              "appD_n": APPD_N, "appD_p": APPD_P, "sparse_n": SPARSE_N, "sparse_deg": SPARSE_DEG,
              "seeds": SEEDS, "gate_max_attempts": GATE_MAX_ATTEMPTS, "cold_repeats": COLD_REPEATS,
              "skip_min_attempts": SKIP_MIN_ATTEMPTS, "skip_censored_share": SKIP_CENSORED_SHARE}
    (out / f"run_config_shard{args.shard}of{args.n_shards}.json").write_text(json.dumps(config, indent=1))

    units = [u for i, u in enumerate(unit_grid(args.smoke)) if i % args.n_shards == args.shard]
    rows = load_rows(out)
    last = {r["uid"]: r for r in rows}
    todo = [u for u in units if u["uid"] not in last
            or (args.retry_errors and last[u["uid"]].get("status") == "error")]
    print(f"E2: {len(units)} units in shard {args.shard}/{args.n_shards}, {len(todo)} to run, "
          f"commit {commit}, criterion {CRITERION}, cap {args.cap:g} s, out {out}", flush=True)

    state = State(rows, len(todo))
    stop = threading.Event()
    threading.Thread(target=heartbeat, args=(state, stop), daemon=True).start()
    fout = open(out / f"rows_shard{args.shard}of{args.n_shards}.jsonl", "a")

    def work(u):
        try:
            row = process_unit(u, args, state, meta)
        except Exception as e:  # noqa: BLE001 -- one broken unit must not stop the sweep
            row = {**u, **meta, "status": "error", "error": repr(e)}
        with state.lock:
            fout.write(json.dumps(row) + "\n")
            fout.flush()
            state.rows.append(row)
            state.done += 1
            state.running.pop(u["uid"], None)
            if row.get("consistency_failures"):
                state.failures += 1
                print(f"CONSISTENCY FAILURE {u['uid']}: {row['consistency_failures']}", flush=True)
            h = row.get("legs", {}).get("hybrid", {})
            print(f"done {u['uid']}: {row['status']} k={row.get('prep', {}).get('k_undirected')} "
                  f"hybrid {h.get('status')} {h.get('method')} r={h.get('radius')} "
                  f"{h.get('seconds', float('nan')):.3g}s", flush=True)

    if args.workers <= 1:
        for u in todo:
            work(u)
    else:
        with ThreadPoolExecutor(max_workers=args.workers) as ex:
            list(ex.map(work, todo))
    stop.set()
    fout.close()
    print(f"E2 finished: {state.done} units, {state.failures} with consistency failures", flush=True)
    return 1 if state.failures else 0


if __name__ == "__main__":
    sys.exit(main())
