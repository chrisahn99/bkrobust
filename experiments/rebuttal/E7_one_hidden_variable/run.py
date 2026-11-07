"""E7: one variable of the true DAG hidden before PC runs, near or far from the query.

For each of the four Gaussian networks of the learned-CPDAG panel and each of the
panel's n = 5000 samples (seeds 0-4), PC runs once on the full sample and once per
hidden variable h with column h dropped. On every learned graph, every frame query
(x, y) of that network gets its committed adjustment set under two knowledge arms,
the set's validity in the full true DAG (GAC), and, for the registered arm, its GAC
radius on the learned graph. One JSON line per unit (network, seed, h), flushed as it
goes; a rerun skips units already written. analyse.py pairs the hidden units with the
nothing-hidden unit of the same network and seed.

    python experiments/rebuttal/E7_one_hidden_variable/run.py --smoke
    python experiments/rebuttal/E7_one_hidden_variable/run.py --workers 6
    python experiments/rebuttal/E7_one_hidden_variable/run.py --check-samples

Samples. Each unit regenerates its panel sample with the panel's ``sample_for`` and keeps
it when its hash equals the one stored in structure.csv. A few panel samples do not
regenerate bit for bit on another platform (one standard-normal draw in 220,000 differs
by one ulp; see DEVIATIONS.md); for those the unit loads the panel's own sample from
panel_samples/ (written by make_panel_samples.py on the machine that made the panel) and
keeps it only when that file's hash equals the stored one. Every unit records which of
the two it used; a cell that matches neither stops the run before any PC runs (exit 3).

Exit codes: 0 done, 2 bad arguments, missing inputs or a config clash with an
existing output, 3 a consistency check failed (see PREREG.md).
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import importlib.metadata
import json
import multiprocessing as mp
import os
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

# one PC per worker; BLAS threads would oversubscribe the cores
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
for _p in (ROOT / "src", HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import numpy as np  # noqa: E402

from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.meek import apply_orientations, is_consistent_extension, meek_closure  # noqa: E402
from bkrobust.gac.dag_level import is_gac_valid_dag  # noqa: E402
from bkrobust.gac.mpdag_level import clear_cache, is_gac_valid_mpdag  # noqa: E402
from bkrobust.hybrid import DEFAULT_SEARCH_BUDGET, breakdown_radius  # noqa: E402
from bkrobust.sat.e1 import radius_e1  # noqa: E402
from bkrobust.search.exact import SearchStats  # noqa: E402
from bkrobust.search.exact_fast import radius_local_up_fast  # noqa: E402
from vendor_pest import (  # noqa: E402
    ALPHA,
    committed_adjustment_set,
    data_sha,
    from_endpoint_matrix,
    graph_from_json,
    is_valid_adjustment_set_gac_dag,
    learn_graph,
    load_gaussian,
    sample_for,
    sha,
    stable_seed,
)

CRITERION = "gac"
NETWORKS = ("ecoli70", "arth150", "magic-niab", "magic-irri")
SMOKE_NETWORKS = ("ecoli70", "magic-niab")
PANEL_N = 5000
SEEDS = (0, 1, 2, 3, 4)
HIDE_SEED = 20261106
PC_WALL_S = 900
ARMS = ("A_ORACLE", "D_DEGEN")
RADIUS_ARMS = ("A_ORACLE",)
NONE = "-"  # the h of the nothing-hidden unit
PANEL_SAMPLES = HERE / "panel_samples"


def panel_sample_path(net: str, seed: int) -> Path:
    return PANEL_SAMPLES / f"{net}_n{PANEL_N}_seed{seed}.npy"


def panel_sample(dag, weights, noise, net: str, seed: int, expected: str):
    """The panel's sample for (net, seed): regenerated, or the vendored copy when regeneration differs.

    Returns ``(data, source, regenerated_sha)``; ``source`` is ``regenerated``, ``vendored``
    or ``unresolved`` (then ``data`` is the regenerated sample and its hash is not ``expected``).
    """
    data = sample_for(dag, weights, noise, net, PANEL_N, seed)
    regen = data_sha(data)
    if regen == expected:
        return data, "regenerated", regen
    path = panel_sample_path(net, seed)
    if path.exists():
        stored = np.load(path, allow_pickle=False)
        if (stored.dtype == np.float64 and stored.shape == data.shape
                and data_sha(np.ascontiguousarray(stored)) == expected):
            return np.ascontiguousarray(stored), "vendored", regen
    return data, "unresolved", regen


class _Wall(Exception):
    pass


def _fire(signum, frame):
    raise _Wall


def repo_commit() -> str:
    try:
        out = subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True, timeout=10
        )
        return out.stdout.strip() or "unknown"
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"


def file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def skeleton_distances(dag: MPDAG, sources: set[str]) -> dict[str, int]:
    """Hop distance from ``sources`` in the true skeleton; unreachable nodes are absent."""
    dist = {s: 0 for s in sources}
    frontier = sorted(sources)
    while frontier:
        nxt = []
        for v in frontier:
            for u in sorted(dag.adjacent(v)):
                if u not in dist:
                    dist[u] = dist[v] + 1
                    nxt.append(u)
        frontier = nxt
    return dist


def hidden_set(dag: MPDAG, net: str, queries: list[dict]) -> dict[str, list[int]]:
    """The network's hidden variables, each with the query indices that drew it."""
    drawn: dict[str, list[int]] = collections.defaultdict(list)
    for q in queries:
        x, y = q["x"], q["y"]
        dist = skeleton_distances(dag, {x, y})
        near = sorted(v for v in dag.nodes if dist.get(v) == 1)
        far = sorted(v for v in dag.nodes if v not in (x, y) and dist.get(v, 10**9) >= 2)
        rng = np.random.default_rng([HIDE_SEED, stable_seed(net), q["i"]])
        # draw from both pools in a fixed order, so one empty pool never shifts the other's draw
        for pool in (near, far):
            if pool:
                drawn[pool[int(rng.integers(len(pool)))]].append(q["i"])
    return dict(sorted(drawn.items()))


# --- the radius ----------------------------------------------------------------------


def radius_dispatch(cpdag: MPDAG, g0: MPDAG, x: str, y: str, z: frozenset, ladder: bool, rung_s: float) -> dict:
    """The submitted hybrid's steps, with the SAT ladder switchable off for large graphs.

    Same order and predicates as ``bkrobust.hybrid.breakdown_radius`` with
    ``criterion="gac"``: top state, bounded upward search, then the E1 ladder. The
    panel sweep skipped the ladder above 60 nodes, which the submitted function
    cannot do, so this is the one place the dispatch is re-written.
    """

    def fails(g: MPDAG) -> bool:
        return not is_gac_valid_mpdag(g, x, y, z)

    if not fails(cpdag):
        return {"radius": UNREACHED, "exact": True, "method": "top_state", "lb": None}
    found = radius_local_up_fast(cpdag, g0, fails, max_depth=DEFAULT_SEARCH_BUDGET, stats=SearchStats())
    if found.radius != UNREACHED:
        method = "degenerate" if found.radius == 0 else "local_up_fast"
        return {"radius": found.radius, "exact": True, "method": method, "lb": None}
    if found.exact:
        return {"radius": UNREACHED, "exact": True, "method": "local_up_fast", "lb": None}
    if not ladder:
        return {"radius": None, "exact": False, "method": "search_lower_bound", "lb": DEFAULT_SEARCH_BUDGET + 1}
    lad = radius_e1(cpdag, g0, x, y, z, time_limit_s=rung_s, criterion=CRITERION)
    if lad.exact:
        return {"radius": lad.radius, "exact": True, "method": "e1_ladder", "lb": None}
    unsat = [k for k, s in lad.rung_status.items() if s == "UNSAT"]
    lb = max(DEFAULT_SEARCH_BUDGET + 1, max(unsat, default=-1) + 1)
    return {"radius": None, "exact": False, "method": "e1_ladder_timeout", "lb": lb}


def crosscheck(mine: dict, cpdag: MPDAG, g0: MPDAG, x: str, y: str, z: frozenset, rung_s: float) -> str | None:
    """Compare the dispatch with the submitted ``breakdown_radius``; a message on disagreement."""
    hr = breakdown_radius(cpdag, None, x, y, z, g0=g0, criterion=CRITERION, time_limit_s=rung_s)
    if hr.exact and mine["exact"]:
        if hr.radius != mine["radius"]:
            return f"breakdown_radius {hr.radius} != dispatch {mine['radius']}"
    elif hr.exact and mine["lb"] is not None and hr.radius != UNREACHED and hr.radius < mine["lb"]:
        return f"breakdown_radius {hr.radius} below dispatch lower bound {mine['lb']}"
    elif mine["exact"] and not hr.exact and mine["radius"] != UNREACHED:
        lb = DEFAULT_SEARCH_BUDGET + 1
        if mine["radius"] < lb:
            return f"dispatch {mine['radius']} below breakdown_radius lower bound {lb}"
    return None


def timed_radius(cpdag, g0, x, y, z, ladder, rung_s, wall_s) -> dict:
    signal.signal(signal.SIGALRM, _fire)
    signal.setitimer(signal.ITIMER_REAL, wall_s)
    t0 = time.perf_counter()
    try:
        res = radius_dispatch(cpdag, g0, x, y, z, ladder, rung_s)
        res["status"] = "ok"
    except _Wall:
        res = {"radius": None, "exact": False, "method": "wall_timeout", "lb": None, "status": "wall_timeout"}
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    res["seconds"] = round(time.perf_counter() - t0, 3)
    return res


# --- one unit: one PC run and every query on its graph ---------------------------------


def unit_task(task: dict) -> dict:
    net, seed, h = task["network"], task["seed"], task["h"]
    dag, weights, noise, _ = load_gaussian(Path(task["models"]), net)
    data, source, regen = panel_sample(dag, weights, noise, net, seed, task["sample_sha256"])
    rec: dict = {
        "unit": f"{net}|{seed}|{h}",
        "network": net,
        "seed": seed,
        "h": h,
        "criterion": CRITERION,
        "repo_commit": task["commit"],
        "causal_learn_version": task["cl_version"],
        "alpha": ALPHA,
        "n": PANEL_N,
        "sample_sha256": data_sha(data),
        "sample_source": source,
        "sample_sha256_regenerated": regen,
        "problems": [],
    }
    if rec["sample_sha256"] != task["sample_sha256"]:
        rec["problems"].append(f"sample hash {rec['sample_sha256']} != panel {task['sample_sha256']}")
    observed = [v for v in dag.nodes if v != h]
    cols = [i for i, v in enumerate(dag.nodes) if v != h]
    signal.signal(signal.SIGALRM, _fire)
    signal.setitimer(signal.ITIMER_REAL, task["pc_wall"])
    t0 = time.perf_counter()
    try:
        matrix = learn_graph(data[:, cols], "pc", ALPHA)
    except _Wall:
        rec.update(status="timeout", seconds=round(time.perf_counter() - t0, 1), rows=[])
        return rec
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    rec["seconds"] = round(time.perf_counter() - t0, 1)
    g, counts = from_endpoint_matrix(matrix, observed)
    closure = meek_closure(g)
    # the panel used the raw graph where the closure failed; kept, and flagged
    used = closure if closure is not None else g
    true_skel = {e for e in dag.skeleton() if h not in e}
    rec.update(
        status="ok",
        closure_ok=int(closure is not None),
        n_conflict_edges=counts["conflict"],
        graph_sha256=sha(used),
        skeleton_missing_vs_truth=len(true_skel - used.skeleton()),
        skeleton_extra_vs_truth=len(used.skeleton() - true_skel),
        graph={"directed": sorted(used.directed_edges), "undirected": sorted(used.undirected_edges)},
    )
    if h == NONE and rec["graph_sha256"] != task["stored_graph_sha256"]:
        rec["problems"].append(f"learned graph {rec['graph_sha256']} != stored panel graph {task['stored_graph_sha256']}")

    knowledge = {
        "A_ORACLE": [
            (a, b) if dag.is_directed_edge(a, b) else (b, a)
            for (a, b) in sorted(used.undirected_edges)
            if dag.has_edge(a, b)
        ],
        "D_DEGEN": [],
    }
    states = {}
    for arm, k in knowledge.items():
        try:
            states[arm] = (apply_orientations(used, k), "")
        except Exception as e:  # a raw, non-closed graph can break the closure code
            states[arm] = (None, type(e).__name__)
    rows = []
    for q in task["queries"]:
        x, y = q["x"], q["y"]
        if h in (x, y):
            continue
        dist = q["dist"].get(h) if h != NONE else None
        for arm in ARMS:
            g0, err = states[arm]
            row = {
                "q": q["i"], "x": x, "y": y, "arm": arm, "n_k": len(knowledge[arm]),
                "dist": dist, "h_common_ancestor": int(h != NONE and h in q["anc_x"] and h in q["anc_y"]),
            }
            if g0 is None:
                row.update(status="knowledge_inconsistent" if not err else "state_error", err=err, committed=0, invalid_truth=0)
                rows.append(row)
                continue
            try:
                cs = committed_adjustment_set(g0, x, y)
            except Exception as e:  # counted as no committed set, never dropped
                row.update(status="error", err=type(e).__name__, committed=0, invalid_truth=0)
                rows.append(row)
                continue
            row.update(status="ok", verdict=cs.verdict, committed=int(cs.z is not None))
            if cs.z is None:
                row["invalid_truth"] = 0
                rows.append(row)
                continue
            z = cs.z
            v_sub = is_gac_valid_dag(dag, x, y, z)
            v_nx = is_valid_adjustment_set_gac_dag(dag, x, y, z)
            if v_sub != v_nx:
                rec["problems"].append(f"{x}->{y} {arm}: GAC at truth submitted {v_sub} != networkx {v_nx}")
            if not is_gac_valid_mpdag(g0, x, y, z):
                rec["problems"].append(f"{x}->{y} {arm}: committed set not GAC-valid on its own graph")
            row.update(z=sorted(z), z_size=len(z), invalid_truth=int(not v_sub))
            if h == NONE:
                ext = is_consistent_extension(dag, g0)
                row["truth_extends_g0"] = int(ext)
                if ext and not v_sub:
                    rec["problems"].append(f"{x}->{y} {arm}: truth extends G0 but Z invalid at truth")
            if arm in RADIUS_ARMS and task["radius"]:
                rr = timed_radius(used, g0, x, y, z, task["ladder"], task["rung_s"], task["radius_wall"])
                row.update({f"radius_{k}": v for k, v in rr.items()})
                if rr["exact"] and rr["radius"] == 0:
                    rec["problems"].append(f"{x}->{y} {arm}: radius 0 for a set valid on G0")
                # deterministic subsample, so a shard split never changes which rows are checked
                tick = int(hashlib.sha256(f"{rec['unit']}|{q['i']}".encode()).hexdigest()[:8], 16)
                if task["ladder"] and rr["status"] == "ok" and tick % task["crosscheck_every"] == 0:
                    # the submitted function has no wall of its own; a timeout here is not a disagreement
                    signal.signal(signal.SIGALRM, _fire)
                    signal.setitimer(signal.ITIMER_REAL, task["radius_wall"])
                    try:
                        msg = crosscheck(rr, used, g0, x, y, z, task["rung_s"])
                        row["crosschecked"] = 1
                    except _Wall:
                        msg, row["crosschecked"] = None, "wall_timeout"
                    finally:
                        signal.setitimer(signal.ITIMER_REAL, 0)
                    if msg:
                        rec["problems"].append(f"{x}->{y} {arm}: {msg}")
            rows.append(row)
    rec["rows"] = rows
    # the closure caches are keyed by graph; every unit's graphs are new
    clear_cache()
    return rec


# --- driver -----------------------------------------------------------------------------


def load_inputs(pest_dir: Path, models: Path, networks, seeds, n_queries):
    """Queries, hidden sets, panel sample hashes and stored graph hashes, per network."""
    frame = [json.loads(line) for line in (pest_dir / "frame.jsonl").open()]
    sample_sha = {
        (r["network"], int(r["seed"])): r["sample_sha256"]
        for r in csv.DictReader((pest_dir / "structure.csv").open())
        if r["method"] == "pc" and int(r["n"]) == PANEL_N and r["status"] == "ok"
    }
    stored = {}
    for line in (pest_dir / "graphs.jsonl").open():
        e = json.loads(line)
        if e["method"] == "pc" and e["n"] == PANEL_N:
            stored[(e["network"], e["seed"])] = sha(graph_from_json(e))
    out = {}
    for net in networks:
        dag = load_gaussian(models, net)[0]
        rows = [r for r in frame if r["network"] == net]
        queries = []
        for i, r in enumerate(rows[:n_queries]):
            x, y = r["x"], r["y"]
            queries.append({
                "i": i, "x": x, "y": y,
                "dist": skeleton_distances(dag, {x, y}),
                "anc_x": sorted(dag.ancestors(x)), "anc_y": sorted(dag.ancestors(y)),
            })
        hidden = hidden_set(dag, net, queries)
        out[net] = {
            "n_nodes": len(dag.nodes),
            "queries": queries,
            "hidden": hidden,
            "sample_sha": {s: sample_sha.get((net, s)) for s in seeds},
            "stored_graph": {s: stored.get((net, s)) for s in seeds},
        }
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--smoke", action="store_true", help="two small networks, seed 0, two queries each")
    ap.add_argument("--out", default=None, help="output folder (default results/rebuttal/E7)")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--networks", default=",".join(NETWORKS))
    ap.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    ap.add_argument("--pc-wall", type=int, default=PC_WALL_S, help="seconds per PC run before a timeout")
    ap.add_argument("--radius-wall", type=float, default=120.0, help="seconds per radius row")
    ap.add_argument("--rung-s", type=float, default=20.0, help="per-rung limit of the SAT ladder")
    ap.add_argument("--ladder-max-nodes", type=int, default=60, help="no SAT ladder above this many nodes")
    ap.add_argument("--crosscheck-every", type=int, default=10)
    ap.add_argument("--no-radius", action="store_true", help="validity only, no radius")
    ap.add_argument("--check-samples", action="store_true",
                    help="only check that every (network, seed) panel sample of the grid is reproduced, then exit")
    ap.add_argument(
        "--pest-dir",
        default=str(ROOT / "experiments" / "rebuttal" / "_data" / "pest"),
        help="folder with frame.jsonl, structure.csv and graphs.jsonl of the learned-CPDAG panel",
    )
    ap.add_argument(
        "--models",
        default=str(ROOT / "results" / "axisa3" / "networks" / "example_models"),
        help="folder with the four Gaussian network JSON files",
    )
    args = ap.parse_args(argv)

    networks = tuple(s for s in args.networks.split(",") if s)
    seeds = tuple(int(s) for s in args.seeds.split(",") if s)
    n_queries, ladder_max, every = 10**9, args.ladder_max_nodes, args.crosscheck_every
    if args.smoke:
        # 46 and 44 nodes: a cut at 45 sends ecoli70 down the no-ladder path and magic-niab up the ladder
        networks, seeds, n_queries, ladder_max, every = SMOKE_NETWORKS, (0,), 2, 45, 1
    workers = args.workers or (1 if args.smoke else 6)
    out = Path(args.out) if args.out else ROOT / "results" / "rebuttal" / ("E7_smoke" if args.smoke else "E7")
    if (not (0 <= args.shard < args.n_shards) or workers < 1 or every < 1
            or set(networks) - set(NETWORKS) or set(seeds) - set(range(10))):
        print("bad --shard/--n-shards/--workers/--networks/--seeds/--crosscheck-every", file=sys.stderr)
        return 2
    pest_dir, models = Path(args.pest_dir), Path(args.models)
    for f in ("frame.jsonl", "structure.csv", "graphs.jsonl"):
        if not (pest_dir / f).exists():
            print(f"missing {pest_dir / f}", file=sys.stderr)
            return 2
    try:
        cl_version = importlib.metadata.version("causal-learn")
    except importlib.metadata.PackageNotFoundError:
        print("causal-learn is not installed (see requirements.txt)", file=sys.stderr)
        return 2
    inputs = load_inputs(pest_dir, models, networks, seeds, n_queries)
    for net, e in inputs.items():
        for s in seeds:
            if e["sample_sha"][s] is None or e["stored_graph"][s] is None:
                print(f"{net} seed {s}: no PC n={PANEL_N} row in structure.csv/graphs.jsonl", file=sys.stderr)
                return 2
    # PREREG check (1) before any PC runs: every cell's sample is the panel's, by hash
    sources: dict[str, dict[str, str]] = {}
    for net in networks:
        dag, weights, noise, _ = load_gaussian(models, net)
        for s in seeds:
            expected = inputs[net]["sample_sha"][s]
            _, src, regen = panel_sample(dag, weights, noise, net, s, expected)
            sources.setdefault(net, {})[str(s)] = src
            print(f"sample {net} seed {s}: {src} (regenerated {regen}, panel {expected})", flush=True)
    counts = collections.Counter(src for d in sources.values() for src in d.values())
    print(f"samples: {sum(counts.values())} cells, {counts['regenerated']} regenerated, "
          f"{counts['vendored']} vendored, {counts['unresolved']} unresolved", flush=True)
    if counts["unresolved"]:
        bad = [f"{n} seed {s}" for n, d in sources.items() for s, v in d.items() if v == "unresolved"]
        print("CONSISTENCY FAILURE: panel sample not reproduced for " + ", ".join(bad), file=sys.stderr)
        return 3
    if args.check_samples:
        return 0
    out.mkdir(parents=True, exist_ok=True)
    commit = repo_commit()
    config = {
        "experiment": "E7",
        "criterion": CRITERION,
        "method": "pc",
        "indep_test": "fisherz",
        "alpha": ALPHA,
        "panel_n": PANEL_N,
        "seeds": list(seeds),
        "networks": list(networks),
        "arms": list(ARMS),
        "radius_arms": [] if args.no_radius else list(RADIUS_ARMS),
        "hide_seed": HIDE_SEED,
        "hidden": {net: list(e["hidden"]) for net, e in inputs.items()},
        "queries_per_network": {net: len(e["queries"]) for net, e in inputs.items()},
        "pc_wall_s": args.pc_wall,
        "radius_wall_s": args.radius_wall,
        "rung_s": args.rung_s,
        "search_budget": DEFAULT_SEARCH_BUDGET,
        "ladder_max_nodes": ladder_max,
        "crosscheck_every": every,
        "smoke": bool(args.smoke),
        "causal_learn_version": cl_version,
        "numpy_version": np.__version__,
        "inputs_sha256": {f: file_sha(pest_dir / f) for f in ("frame.jsonl", "structure.csv", "graphs.jsonl")},
        "models_sha256": {net: file_sha(models / f"{net}.json") for net in networks},
        "vendor_sha256": file_sha(HERE / "vendor_pest.py"),
        "sample_sources": sources,
        "panel_samples_sha256": {
            f"{n}|{s}": file_sha(panel_sample_path(n, int(s)))
            for n, d in sources.items() for s, v in d.items() if v == "vendored"
        },
        "run_sha256": file_sha(Path(__file__)),
        # a later edit of the preregistration then shows up as a config clash
        "prereg_sha256": file_sha(HERE / "PREREG.md"),
    }
    cfg_path = out / "config.json"
    if cfg_path.exists():
        old = json.loads(cfg_path.read_text())
        if {k: v for k, v in old.items() if k != "repo_commit"} != config:
            print(f"{cfg_path} holds a different configuration; use another --out", file=sys.stderr)
            return 2
    cfg_path.write_text(json.dumps({**config, "repo_commit": commit}, indent=1))

    # nothing-hidden units first, so a pipeline that does not reproduce the panel stops early;
    # within that, the largest network first, for balance
    size = {net: inputs[net]["n_nodes"] for net in networks}
    units = [(net, s, NONE) for s in seeds for net in networks]
    units += sorted(
        ((net, s, h) for net in networks for s in seeds for h in inputs[net]["hidden"]),
        key=lambda u: (-size[u[0]], u[0], u[1], u[2]),
    )
    units = [u for i, u in enumerate(units) if i % args.n_shards == args.shard]
    units_path = out / f"units_{args.shard:03d}_of_{args.n_shards:03d}.jsonl"
    done = set()
    if units_path.exists():
        for line in units_path.read_text().splitlines():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue  # a line cut by a kill; the unit is simply redone
            if r.get("problems"):
                print(f"CONSISTENCY FAILURE already recorded in {units_path}: {r['unit']}: {r['problems'][0]}", file=sys.stderr)
                return 3
            done.add(r["unit"])
    todo = [u for u in units if f"{u[0]}|{u[1]}|{u[2]}" not in done]
    print(f"shard {args.shard}/{args.n_shards}: {len(units)} units, {len(done)} already written, "
          f"{len(todo)} to run on {workers} workers", flush=True)
    tasks = [
        {
            "network": net, "seed": s, "h": h, "models": str(models),
            "queries": inputs[net]["queries"],
            "sample_sha256": inputs[net]["sample_sha"][s],
            "stored_graph_sha256": inputs[net]["stored_graph"][s],
            "commit": commit, "cl_version": cl_version, "pc_wall": args.pc_wall,
            "radius": not args.no_radius, "ladder": size[net] <= ladder_max, "rung_s": args.rung_s,
            "radius_wall": args.radius_wall, "crosscheck_every": every,
        }
        for net, s, h in todo
    ]

    t0 = time.perf_counter()
    state = {"done": 0, "last": "-"}
    stop = threading.Event()

    def heartbeat():
        while not stop.wait(30):
            print(f"  .. {state['done']}/{len(tasks)} units, last {state['last']}, {time.perf_counter() - t0:.0f} s", flush=True)

    threading.Thread(target=heartbeat, daemon=True).start()
    failures: list[str] = []
    with units_path.open("a") as fh, mp.get_context("spawn").Pool(workers) as pool:
        for rec in pool.imap_unordered(unit_task, tasks):
            # the line carries every row, so the unit is written whole or not at all
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            state["done"] += 1
            state["last"] = f"{rec['unit']} {rec['status']} {rec.get('seconds')} s"
            print(f"  [{state['done']}/{len(tasks)}] {state['last']} ({time.perf_counter() - t0:.0f} s)", flush=True)
            if rec["problems"]:
                failures += [f"{rec['unit']}: {p}" for p in rec["problems"]]
                pool.terminate()
                break
    stop.set()
    if failures:
        print("CONSISTENCY FAILURE:", file=sys.stderr)
        for f in failures[:20]:
            print("  " + f, file=sys.stderr)
        return 3
    print(f"shard {args.shard}: finished in {time.perf_counter() - t0:.0f} s -> {units_path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
