"""E10: the breakdown radius across every valid adjustment set of a query, criterion gac.

Two parts (see PREREG.md):

* ``networks`` -- each certifiable released query (``status = ok``) of the elicited arms and
  D_RAND, rebuilt from the release and checked against its hashes. Every set valid in G0 is
  listed (exhaustively when at most 12 nodes lie outside the forbidden set, else at most 20 by a
  seeded walk), and each gets a radius and its validity at the true DAG.
* ``census`` -- every CPDAG on 4 and 5 nodes, every knowledge state G0, every ordered pair; every
  valid set gets the submitted hybrid radius and the exact BFS radius on the enumerated space.

One JSON line per unit (a released row, or one CPDAG of the census), flushed as it goes and
skipped on resume. Exit 1 on a consistency failure, 2 if some units raised.

    .venv/bin/python experiments/rebuttal/E10_radius_across_valid_sets/run.py --smoke
    .venv/bin/python experiments/rebuttal/E10_radius_across_valid_sets/run.py --workers 1
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import multiprocessing as mp
import subprocess
import sys
import threading
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
for p in (ROOT / "src", HERE):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import numpy as np  # noqa: E402

from bkrobust.benchmarks.describe import parse_file, to_mpdag  # noqa: E402
from bkrobust.core import oracle as core_oracle  # noqa: E402
from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.core.spacelib import distances_from, radius as space_radius  # noqa: E402
from bkrobust.demo.evaluate import (  # noqa: E402
    asymptotic_variance,
    optimal_adjustment_set_mpdag,
    random_sem,
)
from bkrobust.demo.example import dag_to_cpdag  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.gac.dag_level import is_gac_valid_dag  # noqa: E402
from bkrobust.gac.mpdag_level import clear_cache, gac_forbidden_set, is_gac_valid_mpdag  # noqa: E402
from bkrobust.hybrid import breakdown_radius  # noqa: E402
from bkrobust.search.conjecture_study import all_cpdags  # noqa: E402
from bkrobust.search.exact import retractable_edges  # noqa: E402
from bkrobust.search.space_fixed import build_space_correct  # noqa: E402
from vendor_ledger import (  # noqa: E402
    RELEASE_SEED,
    committed_adjustment_set,
    is_valid_adjustment_set_gac_dag,
    k_sha,
    possible_ancestors,
    rebuild_knowledge,
    sha,
)
from vendor_search import bfs_radius  # noqa: E402

CRITERION = "gac"
DATA = ROOT / "experiments" / "rebuttal" / "_data"
ELICITED_ARMS = (
    "D_LLM", "D_LLM_INSTR", "D_LLM_SRC", "D_SCRAMBLED", "D_LLM_72B", "D_LLM_72B_INSTR",
    "D_SCRAMBLED_72B", "D_LLM_32B", "D_LLM_SRC_70B", "D_LLM_GEMMA_27B",
)
ARMS = (*ELICITED_ARMS, "D_RAND")
# Released sweep settings (ledger_sweep.py): exhaustive search depth when few edges are
# retractable, else the hybrid's default 3; the ladder only up to 60 nodes.
EXHAUSTIVE_MAX_RETRACTABLE = 12
DEFAULT_DEPTH = 3
LADDER_MAX_NODES = 60
CENSUS_SIZES = ((4, 6), (5, 4))  # (nodes, max undirected edges), as in frontier_totals.json
SAMPLE_SALT = 20261106
N_SEMS = 5
SMOKE_ROWS = (("asia", "D_LLM", None), ("mediator", "D_RAND", None),
              ("Sebastiani_2005", None, "e1_ladder"), ("magic-irri", "D_LLM", None))


def zkey(z) -> str:
    return "|".join(sorted(z))


def unit_seed(key: str) -> int:
    return int(hashlib.sha256(f"{SAMPLE_SALT}|{key}".encode()).hexdigest()[:12], 16)


def file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repo_commit() -> tuple[str, bool | None]:
    try:
        head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--", "src"],
                               capture_output=True, text=True, check=True).stdout.strip()
        return head, bool(dirty)
    except (OSError, subprocess.CalledProcessError):
        return "unknown", None


def clear_all() -> None:
    clear_cache()  # gac and criterion memo tables: pure functions of the graph, cost only
    core_oracle.clear_cache()


# ------------------------------------------------------------------ radius of one set


def set_radius(cpdag, g0, x, y, z, n_nodes: int, n_ret: int, cfg: dict) -> dict:
    """Radius of one set under gac, with the method, exactness and a lower bound if censored."""
    depth = n_ret if n_ret <= EXHAUSTIVE_MAX_RETRACTABLE else DEFAULT_DEPTH

    def fails(g):
        return not is_gac_valid_mpdag(g, x, y, z)

    t0 = time.perf_counter()
    if n_nodes <= LADDER_MAX_NODES:
        hb = breakdown_radius(cpdag, None, x, y, z, g0=g0, search_budget=depth,
                              criterion=CRITERION, time_limit_s=cfg["ladder_time_limit"])
        rec = {"r": hb.radius if hb.exact else None, "exact": hb.exact, "method": hb.method,
               "lb": None if hb.exact else depth + 1, "budget": depth}
        if hb.method == "e1_ladder":
            # The ladder shares no search code with the BFS: run both, they must agree.
            bfs = bfs_radius(cpdag, g0, fails, max_depth=None, max_closures=cfg["closure_budget"],
                             time_limit_s=cfg["set_time_limit"])
            rec["bfs_check"] = {"r": bfs["radius"], "exact": bfs["exact"], "lb": bfs["lb"],
                                "stop": bfs["stop"]}
            if not hb.exact and bfs["exact"]:
                rec.update({"r": bfs["radius"], "exact": True, "lb": None, "method": "e1_ladder+bfs"})
    else:
        if not fails(cpdag):
            rec = {"r": UNREACHED, "exact": True, "method": "top_state", "lb": None, "budget": None}
        else:
            bfs = bfs_radius(cpdag, g0, fails, max_depth=None, max_closures=cfg["closure_budget"],
                             time_limit_s=cfg["set_time_limit"])
            rec = {"r": bfs["radius"] if bfs["exact"] else None, "exact": bfs["exact"],
                   "method": "bfs_budget", "lb": bfs["lb"], "budget": cfg["closure_budget"],
                   "stop": bfs["stop"]}
    rec["seconds"] = round(time.perf_counter() - t0, 4)
    return rec


# ------------------------------------------------------------------ valid sets of one query


def exhaustive_sets(g0, x, y, pool: list[str], cap: int, rng) -> tuple[list, int, bool]:
    found = []
    for size in range(len(pool) + 1):
        for combo in itertools.combinations(pool, size):
            z = frozenset(combo)
            if is_gac_valid_mpdag(g0, x, y, z):
                found.append(z)
    total = len(found)
    if total <= cap:
        return found, total, False
    keep = sorted(rng.choice(total, size=cap, replace=False).tolist())
    return [found[i] for i in keep], total, True


def sampled_sets(g0, x, y, z_o, forb, pool, rng, cfg) -> tuple[dict, int]:
    """O, the canonical set, and valid neighbours found by a seeded walk; at most max_sets."""
    found: dict[frozenset, str] = {z_o: "optimal"}
    checks = 0
    t_end = time.perf_counter() + cfg["sample_time_limit"]

    def valid(z) -> bool:
        nonlocal checks
        checks += 1
        return is_gac_valid_mpdag(g0, x, y, z)

    canon = frozenset((possible_ancestors(g0, x) | possible_ancestors(g0, y)) - forb - {x, y})
    if canon not in found and valid(canon):
        found[canon] = "canonical"
    anchor = {x, y} | set(z_o) | set(canon)
    near = sorted(v for v in pool if any(g0.has_edge(v, a) for a in anchor)) or list(pool)
    if canon in found and canon:
        # One greedy shrink from the canonical set; every intermediate set is valid.
        cur = set(canon)
        for v in rng.permutation(sorted(canon)).tolist():
            if len(found) >= cfg["max_sets"] or checks >= cfg["sample_checks"] or time.perf_counter() > t_end:
                break
            cand = frozenset(cur - {v})
            if valid(cand):
                cur = set(cand)
                found.setdefault(cand, "shrink")
    while len(found) < cfg["max_sets"] and checks < cfg["sample_checks"] and time.perf_counter() < t_end:
        base = list(found)[int(rng.integers(len(found)))]
        if base and rng.random() < 0.5:
            cand = frozenset(base - {sorted(base)[int(rng.integers(len(base)))]})
            origin = "drop"
        else:
            extra = [v for v in near if v not in base]
            if not extra:
                checks += 1  # counts as an attempt, so the loop always ends
                continue
            cand = frozenset(base | {extra[int(rng.integers(len(extra)))]})
            origin = "add"
        if cand in found:
            checks += 1
            continue
        if valid(cand):
            found[cand] = origin
    return found, checks


# ------------------------------------------------------------------ networks part

_NETS: dict[str, tuple] = {}
_KNOW: dict[str, dict] = {}


def knowledge(path: str) -> dict:
    # Read once per process, so the task tuples sent to workers stay small.
    if path not in _KNOW:
        _KNOW[path] = json.loads(Path(path).read_text())
    return _KNOW[path]


def load_network(models: Path, name: str):
    """DAG, CPDAG and file name, by ledger_sweep.load's rule; one network cached per process."""
    if name not in _NETS:
        path = next(p for p in sorted(models.iterdir()) if p.name.startswith(name))
        dag = to_mpdag(parse_file(path, file_sha(path)))
        _NETS.clear()
        _NETS[name] = (dag, dag_to_cpdag(dag), path.name)
    return _NETS[name]


def net_key(r: dict) -> str:
    return "net|" + "|".join((r["network"], r["arm"], r["rep"], r["x"], r["y"]))


def more_robust(a: int, b: int) -> bool:
    if a == UNREACHED:
        return b != UNREACHED
    if b == UNREACHED:
        return False
    return a > b


def variance_cost(dag, x, y, z_o, z_best, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    ratios = []
    for _ in range(N_SEMS):
        sem = random_sem(dag, rng)
        try:
            ratios.append(asymptotic_variance(sem, x, y, z_best) / asymptotic_variance(sem, x, y, z_o))
        except Exception:  # a singular covariance gives no ratio, not a crash
            ratios.append(None)
    return {"z_best": sorted(z_best), "avar_ratio_best_over_o": ratios}


def compute_network(row: dict, cfg: dict) -> dict:
    dag, cpdag, model_file = load_network(Path(cfg["models"]), row["network"])
    x, y, arm, rep = row["x"], row["y"], row["arm"], int(row["rep"])
    res: dict = {"model_file": model_file, "n_nodes": len(cpdag.nodes)}
    mism: list[str] = []

    def expect(name, got, want):
        res.setdefault("rebuild", {})[name] = got
        if str(got) != str(want):
            mism.append(f"{name}: got {got}, released {want}")

    expect("chat_source_sha256", sha(cpdag), row["chat_source_sha256"])
    expect("generating_dag_sha256", sha(dag), row["generating_dag_sha256"])
    k = rebuild_knowledge(row["network"], arm, rep, dag, cpdag,
                          knowledge(cfg["knowledge"]), cfg["seed"])
    g0 = apply_orientations(cpdag, k) if k is not None else None
    if g0 is None:
        return {**res, "status": "rebuild_mismatch", "hard_failures": ["knowledge or G0 not rebuilt"]}
    expect("k_sha256", k_sha(k), row["k_sha256"])
    expect("n_k", len(k), row["n_k"])
    expect("g0_undirected", len(g0.undirected_edges), row["g0_undirected"])
    cs = committed_adjustment_set(g0, x, y)
    expect("committed_verdict", cs.verdict, row["committed_verdict"])
    if cs.z is None:
        return {**res, "status": "rebuild_mismatch", "hard_failures": mism + ["no committed set"]}
    z_o = cs.z
    expect("z_size", len(z_o), row["z_size"])
    if mism:
        return {**res, "status": "rebuild_mismatch", "hard_failures": mism}

    hard: list[str] = []
    rng = np.random.default_rng(unit_seed(net_key(row)))
    forb = gac_forbidden_set(g0, x, y)
    pool = sorted(set(cpdag.nodes) - forb - {x, y})
    n_ret = len(retractable_edges(cpdag, g0))
    t0 = time.perf_counter()
    if len(pool) <= cfg["exhaustive_max_candidates"]:
        sets, n_total, truncated = exhaustive_sets(g0, x, y, pool, cfg["max_sets_exhaustive"], rng)
        if z_o not in sets:
            if n_total <= cfg["max_sets_exhaustive"]:
                hard.append("committed set missing from the exhaustive list")
            sets = [z_o, *sets]
        origins = {z: ("optimal" if z == z_o else "exhaustive") for z in sets}
        mode, checks = "exhaustive", 2 ** len(pool)
    else:
        origins, checks = sampled_sets(g0, x, y, z_o, forb, pool, rng, cfg)
        sets, n_total, truncated, mode = list(origins), None, None, "sampled"
    enum_s = round(time.perf_counter() - t0, 3)
    sets = [z_o] + sorted((z for z in sets if z != z_o), key=lambda z: (len(z), sorted(z)))

    out_sets = []
    for z in sets:
        r = set_radius(cpdag, g0, x, y, z, len(cpdag.nodes), n_ret, cfg)
        t_sub = is_gac_valid_dag(dag, x, y, z)
        t_nx = is_valid_adjustment_set_gac_dag(dag, x, y, z)
        if t_sub != t_nx:
            hard.append(f"truth oracles disagree on {zkey(z)}")
        bc = r.get("bfs_check")
        if bc and r["exact"] and bc["exact"] and bc["r"] != r["r"] and r["method"] == "e1_ladder":
            hard.append(f"ladder {r['r']} and search {bc['r']} disagree on {zkey(z)}")
        if bc and bc["exact"] is False and r["exact"] and r["r"] != UNREACHED and bc["lb"] > r["r"]:
            hard.append(f"search certifies r >= {bc['lb']} but ladder says {r['r']} on {zkey(z)}")
        out_sets.append({"z": sorted(z), "size": len(z), "origin": origins.get(z, "exhaustive"),
                         "valid_at_truth": t_sub, **r})
        if len(out_sets) % 50 == 0:
            clear_all()

    o = out_sets[0]
    if int(o["valid_at_truth"]) != int(row["z_valid_at_truth"]):
        hard.append(f"O valid_at_truth {o['valid_at_truth']} vs released {row['z_valid_at_truth']}")
    released_exact = row["hop_exact"] == "1" and row["r_hop"] not in ("",)
    if released_exact and o["exact"] and int(row["r_hop"]) != o["r"]:
        hard.append(f"O radius {o['r']} vs released r_hop {row['r_hop']}")

    res.update({
        "status": "ok", "z_o": sorted(z_o), "n_retractable": n_ret, "n_k": len(k),
        "pool_size": len(pool), "mode": mode, "n_valid_total": n_total, "truncated": truncated,
        "validity_checks": checks, "enum_seconds": enum_s, "n_sets": len(out_sets), "sets": out_sets,
        "released_r_hop": row["r_hop"], "released_hop_exact": row["hop_exact"],
        "hard_failures": hard,
    })
    # Variance cost of the most robust set when some set is certified more robust than O.
    if o["exact"]:
        better = [s for s in out_sets[1:] if (s["exact"] and more_robust(s["r"], o["r"]))
                  or (not s["exact"] and o["r"] != UNREACHED and s["lb"] is not None and s["lb"] > o["r"])]
        if better:
            best = max(better, key=lambda s: (s["exact"] and s["r"] == UNREACHED,
                                              s["r"] if s["exact"] else s["lb"], -s["size"]))
            res["variance_cost"] = variance_cost(dag, x, y, z_o, frozenset(best["z"]),
                                                 unit_seed("sem|" + net_key(row)))
            res["variance_cost"]["best_valid_at_truth"] = best["valid_at_truth"]
    return res


def run_network_unit(row: dict, cfg: dict) -> dict:
    rec = {"key": net_key(row), "part": "networks",
           **{f: row[f] for f in ("network", "tier", "x", "y", "arm", "rep")},
           "criterion": CRITERION, "commit": cfg["commit"], "smoke": cfg["smoke"]}
    try:
        rec.update(compute_network(row, cfg))
    except Exception:  # reported, and retried on resume
        rec["status"] = "error"
        rec["error"] = traceback.format_exc()[-2000:]
    return rec


# ------------------------------------------------------------------ census part

_CPDAGS: dict[tuple[int, int], list] = {}


def census_cpdags(n: int, max_und: int) -> list:
    if (n, max_und) not in _CPDAGS:
        _CPDAGS[(n, max_und)] = [c for c in all_cpdags(n) if 1 <= len(c.undirected_edges) <= max_und]
    return _CPDAGS[(n, max_und)]


def compute_census(n: int, max_und: int, idx: int, cfg: dict) -> dict:
    cpdag = census_cpdags(n, max_und)[idx]
    space = build_space_correct(cpdag)
    hard: list[str] = []
    if cpdag not in space.neighbours:
        return {"status": "ok", "cpdag": cpdag.edge_string(), "space_size": len(space),
                "instances": [], "hard_failures": ["CPDAG missing from its own space"]}
    nodes = list(cpdag.nodes)
    instances, n_pairs, n_lt2, n_cens = [], 0, 0, 0
    for gi, g0 in enumerate(space.elements):
        dists = distances_from(space, g0)
        for x, y in itertools.permutations(nodes, 2):
            n_pairs += 1
            others = [v for v in nodes if v not in (x, y)]
            valid = [frozenset(c) for s in range(len(others) + 1)
                     for c in itertools.combinations(others, s)
                     if is_gac_valid_mpdag(g0, x, y, frozenset(c))]
            if len(valid) < 2:
                n_lt2 += 1
                continue
            cs = committed_adjustment_set(g0, x, y)
            o_enum = optimal_adjustment_set_mpdag(g0, x, y)
            sets = []
            for z in valid:
                hb = breakdown_radius(cpdag, None, x, y, z, g0=g0, criterion=CRITERION,
                                      time_limit_s=cfg["ladder_time_limit"])

                def fails(g, _z=z, _x=x, _y=y):
                    return not is_gac_valid_mpdag(g, _x, _y, _z)

                rb, _ = space_radius(space, dists, fails)
                if not hb.exact:
                    n_cens += 1
                elif hb.radius != rb:
                    hard.append(f"g0 {gi} {x}->{y} Z={zkey(z)}: hybrid {hb.radius} vs space BFS {rb}")
                # [Z, hybrid radius or None, method, space-BFS radius]
                sets.append([zkey(z), hb.radius if hb.exact else None, hb.method, rb])
            instances.append({
                "g0": gi, "x": x, "y": y, "verdict": cs.verdict,
                "o": None if cs.z is None else zkey(cs.z),
                "o_enum": None if o_enum is None else zkey(o_enum),
                "sets": sets,
            })
    return {"status": "ok", "cpdag": cpdag.edge_string(), "n_undirected": len(cpdag.undirected_edges),
            "space_size": len(space), "n_pairs": n_pairs, "n_fewer_than_2_sets": n_lt2,
            "n_censored": n_cens, "instances": instances, "hard_failures": hard}


def run_census_unit(task: tuple, cfg: dict) -> dict:
    n, max_und, idx = task
    rec = {"key": f"census|{n}|{max_und}|{idx}", "part": "census", "n": n, "max_undirected": max_und,
           "cpdag_index": idx, "criterion": CRITERION, "commit": cfg["commit"], "smoke": cfg["smoke"],
           "space": "build_space_correct"}
    try:
        rec.update(compute_census(n, max_und, idx, cfg))
    except Exception:
        rec["status"] = "error"
        rec["error"] = traceback.format_exc()[-2000:]
    clear_all()
    return rec


def run_task(task: tuple) -> dict:
    kind, payload, cfg = task
    t0 = time.perf_counter()
    rec = run_census_unit(payload, cfg) if kind == "census" else run_network_unit(payload, cfg)
    clear_all()
    rec["seconds"] = round(time.perf_counter() - t0, 3)
    return rec


# ------------------------------------------------------------------ driver


def select_network_rows(records: Path, arms: list[str], smoke: bool) -> list[dict]:
    rows = [r for r in csv.DictReader(records.open()) if r["status"] == "ok" and r["arm"] in arms]
    rows.sort(key=net_key)
    if not smoke:
        return rows
    picked = []
    for net, arm, method in SMOKE_ROWS:
        hit = next((r for r in rows if r["network"] == net and (arm is None or r["arm"] == arm)
                    and (method is None or r["hop_method"] == method)), None)
        if hit is not None:
            picked.append(hit)
    return picked


def census_tasks(sizes, smoke: bool) -> list[tuple]:
    if smoke:
        # every 3-node CPDAG, and the first 4-node one with two undirected edges
        first4 = next(i for i, c in enumerate(census_cpdags(4, 6)) if len(c.undirected_edges) == 2)
        return [(3, 3, i) for i in range(len(census_cpdags(3, 3)))] + [(4, 6, first4)]
    return [(n, m, i) for n, m in sizes for i in range(len(census_cpdags(n, m)))]


def done_keys(out: Path) -> set[str]:
    keys = set()
    for f in sorted(out.glob("rows_*.jsonl")):
        for line in f.open():
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue  # a line cut by a kill; that unit is rerun
            if rec.get("status") != "error":
                keys.add(rec["key"])
    return keys


def task_key(kind: str, payload) -> str:
    return f"census|{payload[0]}|{payload[1]}|{payload[2]}" if kind == "census" else net_key(payload)


def parse_sizes(text: str) -> tuple:
    return tuple(tuple(int(v) for v in part.split(":")) for part in text.split(",") if part)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--part", choices=("all", "census", "networks"), default="all")
    ap.add_argument("--records", type=Path, default=DATA / "per_query_records.csv")
    ap.add_argument("--knowledge", type=Path, default=DATA / "knowledge_states.json")
    ap.add_argument("--models", type=Path, default=ROOT / "results/axisa3/networks/example_models")
    ap.add_argument("--out", type=Path, default=None, help="default results/rebuttal/E10 (E10_smoke with --smoke)")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--arms", default=",".join(ARMS))
    ap.add_argument("--census-sizes", default=",".join(f"{n}:{m}" for n, m in CENSUS_SIZES))
    ap.add_argument("--seed", type=int, default=RELEASE_SEED, help="D_RAND knowledge seed (released)")
    ap.add_argument("--max-sets", type=int, default=20, help="sets per query on large pools")
    ap.add_argument("--max-sets-exhaustive", type=int, default=512)
    ap.add_argument("--exhaustive-max-candidates", type=int, default=12)
    ap.add_argument("--sample-checks", type=int, default=200, help="validity checks per sampled query")
    ap.add_argument("--sample-time-limit", type=float, default=120.0)
    ap.add_argument("--closure-budget", type=int, default=4000)
    ap.add_argument("--set-time-limit", type=float, default=20.0)
    ap.add_argument("--ladder-time-limit", type=float, default=60.0)
    args = ap.parse_args(argv)
    if args.smoke:
        # tiny budgets: every code path runs once, the numbers mean nothing
        args.max_sets, args.max_sets_exhaustive, args.sample_checks = 3, 64, 30
        args.sample_time_limit, args.closure_budget = 5.0, 300
        args.set_time_limit, args.ladder_time_limit = 1.0, 3.0
    out = args.out or ROOT / ("results/rebuttal/E10_smoke" if args.smoke else "results/rebuttal/E10")
    out.mkdir(parents=True, exist_ok=True)
    if not 0 <= args.shard < args.n_shards:
        ap.error("need 0 <= --shard < --n-shards")

    commit, src_dirty = repo_commit()
    arms = [a for a in args.arms.split(",") if a]
    sizes = parse_sizes(args.census_sizes)
    tasks_all: list[tuple] = []
    cfg = {
        "experiment": "E10", "criterion": CRITERION, "commit": commit, "src_dirty": src_dirty,
        "smoke": args.smoke, "part": args.part, "arms": arms, "seed": args.seed,
        "census_sizes": [list(s) for s in sizes], "max_sets": args.max_sets,
        "max_sets_exhaustive": args.max_sets_exhaustive,
        "exhaustive_max_candidates": args.exhaustive_max_candidates,
        "sample_checks": args.sample_checks, "sample_time_limit": args.sample_time_limit,
        "closure_budget": args.closure_budget, "set_time_limit": args.set_time_limit,
        "ladder_time_limit": args.ladder_time_limit, "ladder_max_nodes": LADDER_MAX_NODES,
        "exhaustive_max_retractable": EXHAUSTIVE_MAX_RETRACTABLE, "models": str(args.models),
        "shard": args.shard, "n_shards": args.n_shards, "workers": args.workers,
    }
    if args.part in ("all", "networks"):
        cfg["records_sha256"] = file_sha(args.records)
        cfg["knowledge_sha256"] = file_sha(args.knowledge)
        cfg["knowledge"] = str(args.knowledge)
        rows = select_network_rows(args.records, arms, args.smoke)
        cfg["n_selected_networks"] = len(rows)
        tasks_all += [("networks", r) for r in rows]
    if args.part in ("all", "census"):
        ct = census_tasks(sizes, args.smoke)
        cfg["n_selected_census"] = len(ct)
        tasks_all += [("census", t) for t in ct]
    tasks_all = [t for i, t in enumerate(tasks_all) if i % args.n_shards == args.shard]
    suffix = "" if args.n_shards == 1 else f".shard{args.shard}of{args.n_shards}"
    (out / f"config{suffix}.json").write_text(
        json.dumps(cfg, indent=1, sort_keys=True))

    have = done_keys(out)
    todo = [(kind, p, cfg) for kind, p in tasks_all if task_key(kind, p) not in have]
    print(f"[E10] {len(tasks_all)} units in this shard, {len(tasks_all) - len(todo)} already done, "
          f"{len(todo)} to run; commit {commit[:12]}, src dirty {src_dirty}", flush=True)

    files = {kind: (out / f"rows_{kind}{suffix}.jsonl").open("a") for kind in ("networks", "census")}
    state = {"done": 0, "hard": 0, "err": 0, "t0": time.perf_counter(), "stop": False}

    def heartbeat():
        while not state["stop"]:
            time.sleep(30)
            if not state["stop"]:
                print(f"[E10] {state['done']}/{len(todo)} done, hard {state['hard']}, errors {state['err']}, "
                      f"{time.perf_counter() - state['t0']:.0f} s", flush=True)

    threading.Thread(target=heartbeat, daemon=True).start()

    def consume(rec: dict):
        fh = files[rec["part"]]
        fh.write(json.dumps(rec, sort_keys=True) + "\n")
        fh.flush()
        state["done"] += 1
        if rec.get("status") == "error":
            state["err"] += 1
            print(f"[E10] ERROR {rec['key']}: {rec['error'].strip().splitlines()[-1]}", flush=True)
        if rec.get("hard_failures"):
            state["hard"] += 1
            print(f"[E10] CONSISTENCY {rec['key']}: {rec['hard_failures'][:3]}", flush=True)
        if args.smoke:
            print(f"[E10] {rec['key']} {rec.get('status')} {rec.get('seconds')} s", flush=True)

    try:
        if args.workers <= 1:
            for t in todo:
                consume(run_task(t))
        else:
            with mp.get_context("spawn").Pool(args.workers) as pool:
                for rec in pool.imap_unordered(run_task, todo, chunksize=1):
                    consume(rec)
    finally:
        state["stop"] = True
        for fh in files.values():
            fh.close()
    print(f"[E10] finished {state['done']} units in {time.perf_counter() - state['t0']:.0f} s; "
          f"hard {state['hard']}, errors {state['err']}", flush=True)
    if state["hard"]:
        print("[E10] CONSISTENCY FAILURE", flush=True)
        return 1
    return 2 if state["err"] else 0


if __name__ == "__main__":
    sys.exit(main())
