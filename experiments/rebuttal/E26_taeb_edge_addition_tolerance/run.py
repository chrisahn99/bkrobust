"""E26: an edge-addition tolerance in the sense of Taeb, Guo and Henckel, Section 6 (see PREREG.md).

For every released query with a committed set (status ok in per_query_records.csv, all arms),
the knowledge state G0 is rebuilt from the release and checked against its hashes, and one DAG D0
of G0's class is fixed (a deterministic consistent extension). The tolerance of the committed set
Z is the least number of directed edges that, added to D0, make Z fail the generalized adjustment
criterion (GAC) for (x, y). Two move sets:

* ``order`` (registered): an added edge a -> b must respect one fixed causal order of D0, the
  topological order of D0 with ties broken by node name (Taeb et al. restrict their supergraphs by
  a temporal order; this is the total-order version of that restriction).
* ``acyclic`` (secondary): any edge that keeps the graph acyclic. Computed only where the order
  tolerance exceeds 1, since every order move is an acyclic move.

Search: iterative deepening to depth 3, exact by construction (validity depends only on the
subgraph induced by the ancestors of {x, y} and Z, so the first edge of a minimal breaking set has
its head there, and each next edge has its head in the ancestors of the graph built so far). A
stopped search gives a lower bound. Every breaking set is re-tested by two independent GAC
implementations from the submitted and released code; a sample of non-breaking moves too.

One JSON line per query, flushed as it goes and skipped on resume. Exit 1 on a consistency
failure, 2 if some units raised.

    .venv/bin/python experiments/rebuttal/E26_taeb_edge_addition_tolerance/run.py --smoke
    .venv/bin/python experiments/rebuttal/E26_taeb_edge_addition_tolerance/run.py --workers 10
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import heapq
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
from bkrobust.demo.example import dag_to_cpdag  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.gac.dag_level import is_gac_valid_dag  # noqa: E402
from vendor_ledger import (  # noqa: E402
    RELEASE_SEED,
    committed_adjustment_set,
    is_valid_adjustment_set_gac_dag,
    k_sha,
    rebuild_knowledge,
    sha,
)

EXPERIMENT = "E26"
CRITERION = "gac"
DATA = ROOT / "experiments" / "rebuttal" / "_data"
MAX_DEPTH = 3
N_RANDOM_EXTENSIONS = 2
XCHECK_NONBREAKING = 10
SALT = 20261107
SMOKE_ROWS = (("asia", "D_LLM"), ("mediator", "D_RAND"), ("Acid_1996", "D_DEGEN"),
              ("alarm", "A_ORACLE"), ("diabetes", "D_LLM_72B"), ("andes", "D_RAND"))
SMOKE_PER_NETWORK = 2


def unit_seed(key: str) -> int:
    return int(hashlib.sha256(f"{SALT}|{key}".encode()).hexdigest()[:12], 16)


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


# ------------------------------------------------------------------ a mutable DAG and the GAC


class DG:
    """Adjacency sets of a DAG; edges are added and removed in place during the search."""

    def __init__(self, nodes, edges):
        self.nodes = list(nodes)
        self.ch = {v: set() for v in self.nodes}
        self.pa = {v: set() for v in self.nodes}
        for a, b in edges:
            self.ch[a].add(b)
            self.pa[b].add(a)

    def add(self, a, b):
        self.ch[a].add(b)
        self.pa[b].add(a)

    def remove(self, a, b):
        self.ch[a].discard(b)
        self.pa[b].discard(a)

    def edges(self):
        return [(a, b) for a in self.nodes for b in self.ch[a]]


def reach(adj: dict, start) -> set:
    seen = set(start)
    stack = list(start)
    while stack:
        v = stack.pop()
        for w in adj[v]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return seen


def gac_valid(g: DG, x, y, z: frozenset) -> bool:
    """GAC in a DAG (Perkovic et al. 2018): forbidden set, then d-separation in the proper
    back-door graph, by a Bayes-ball search. Mirrors vendor_ledger.is_valid_adjustment_set_gac_dag."""
    if x in z or y in z:
        return False
    de_x = reach(g.ch, [x]) - {x}
    causal = de_x & reach(g.pa, [y])
    forb = reach(g.ch, causal) | {x}
    if z & forb:
        return False
    an_z = reach(g.pa, z)
    # Bayes ball from x; the edges x -> c, c causal, are removed.
    seen = set()
    stack = [(x, "up")]
    while stack:
        v, d = stack.pop()
        if (v, d) in seen:
            continue
        seen.add((v, d))
        if v == y:
            return False
        if v in z:
            if d == "down" and v in an_z:
                for p in g.pa[v]:
                    if not (p == x and v in causal):
                        stack.append((p, "up"))
            continue
        if d == "up":
            for p in g.pa[v]:
                if not (p == x and v in causal):
                    stack.append((p, "up"))
            for c in g.ch[v]:
                if not (v == x and c in causal):
                    stack.append((c, "down"))
        else:
            for c in g.ch[v]:
                if not (v == x and c in causal):
                    stack.append((c, "down"))
            if v in an_z:
                for p in g.pa[v]:
                    if not (p == x and v in causal):
                        stack.append((p, "up"))
    return True


def topo_order(g: DG) -> dict:
    indeg = {v: len(g.pa[v]) for v in g.nodes}
    heap = [v for v in g.nodes if indeg[v] == 0]
    heapq.heapify(heap)
    order = []
    while heap:
        v = heapq.heappop(heap)
        order.append(v)
        for c in sorted(g.ch[v]):
            indeg[c] -= 1
            if indeg[c] == 0:
                heapq.heappush(heap, c)
    if len(order) != len(g.nodes):
        raise ValueError("graph has a cycle")
    return {v: i for i, v in enumerate(order)}


# ------------------------------------------------------------------ one DAG of G0's class


def extension(g0: MPDAG, rng=None) -> MPDAG:
    """A consistent DAG extension of an MPDAG: orient one undirected edge, Meek-close, repeat.
    In a Meek-closed graph either orientation of any undirected edge is consistent (Meek 1995),
    so the loop never fails; the second orientation is a guard. rng=None is the deterministic
    extension (smallest edge first, its sorted direction)."""
    g = g0
    while g.undirected_edges:
        und = sorted(g.undirected_edges)
        if rng is None:
            a, b = und[0]
        else:
            a, b = und[int(rng.integers(len(und)))]
            if rng.integers(2):
                a, b = b, a
        nxt = apply_orientations(g, [(a, b)])
        if nxt is None:
            nxt = apply_orientations(g, [(b, a)])
        if nxt is None:
            raise ValueError(f"no consistent orientation of {a} - {b}")
        g = nxt
    return g


# ------------------------------------------------------------------ the tolerance search


class Stop(Exception):
    pass


def candidates(g: DG, x, y, z: frozenset, mode: str, pos: dict) -> list:
    """Moves whose head is an ancestor of {x, y} and Z in g, ordered so likely breakers come first."""
    target = {x, y} | set(z)
    an_t = reach(g.pa, target)
    rank = {y: 0, x: 1}
    out = []
    for b in an_t:
        if mode == "order":
            tails = [a for a in g.nodes if pos[a] < pos[b] and a not in g.pa[b]]
        else:
            de_b = reach(g.ch, [b])
            tails = [a for a in g.nodes if a not in de_b and a not in g.pa[b]]
        hr = rank.get(b, 2 if b in z else 3)
        for a in tails:
            out.append((hr, 0 if a in an_t else 1, a, b))
    out.sort()
    return [(a, b) for _, _, a, b in out]


def tolerance(d0: MPDAG, x, y, z: frozenset, mode: str, max_depth: int, time_limit: float,
              pos: dict | None = None) -> dict:
    g = DG(d0.nodes, d0.directed_edges)
    if pos is None:
        pos = topo_order(g)
    t_end = time.perf_counter() + time_limit
    stats = {"checks": 0, "nonbreaking_sample": []}

    def dfs(k: int, added: list):
        for a, b in candidates(g, x, y, z, mode, pos):
            if time.perf_counter() > t_end:
                raise Stop
            g.add(a, b)
            added.append((a, b))
            try:
                if k == 1:
                    stats["checks"] += 1
                    if not gac_valid(g, x, y, z):
                        return list(added)
                    if len(stats["nonbreaking_sample"]) < XCHECK_NONBREAKING and len(added) == 1:
                        stats["nonbreaking_sample"].append((a, b))
                else:
                    hit = dfs(k - 1, added)
                    if hit:
                        return hit
            finally:
                added.pop()
                g.remove(a, b)
        return None

    t0 = time.perf_counter()
    rec = {"mode": mode, "max_depth": max_depth}
    done_level = 0
    try:
        for k in range(1, max_depth + 1):
            hit = dfs(k, [])
            if hit:
                rec.update(t=k, exact=True, lb=None, witness=hit)
                break
            done_level = k
        else:
            rec.update(t=None, exact=False, lb=max_depth + 1, witness=None, stop="depth")
    except Stop:
        rec.update(t=None, exact=False, lb=done_level + 1, witness=None, stop="time")
    rec["checks"] = stats["checks"]
    rec["seconds"] = round(time.perf_counter() - t0, 4)
    rec["_sample"] = stats["nonbreaking_sample"]
    return rec


def with_edges(d0: MPDAG, extra) -> MPDAG:
    return MPDAG(d0.nodes, set(d0.directed_edges) | {tuple(e) for e in extra}, [])


# ------------------------------------------------------------------ one released query

_NETS: dict[str, tuple] = {}
_KNOW: dict[str, dict] = {}


def knowledge(path: str) -> dict:
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


def row_key(r: dict) -> str:
    return "q|" + "|".join((r["network"], r["arm"], r["rep"], r["x"], r["y"]))


def compute(row: dict, cfg: dict) -> dict:
    dag, cpdag, model_file = load_network(Path(cfg["models"]), row["network"])
    x, y, arm, rep = row["x"], row["y"], row["arm"], int(row["rep"])
    res: dict = {"model_file": model_file, "n_nodes": len(cpdag.nodes)}
    mism: list[str] = []

    def expect(name, got, want):
        if str(got) != str(want):
            mism.append(f"{name}: got {got}, released {want}")

    expect("chat_source_sha256", sha(cpdag), row["chat_source_sha256"])
    expect("generating_dag_sha256", sha(dag), row["generating_dag_sha256"])
    k = rebuild_knowledge(row["network"], arm, rep, dag, cpdag, knowledge(cfg["knowledge"]), cfg["seed"])
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
    z = frozenset(cs.z)
    expect("z_size", len(z), row["z_size"])
    t_truth = is_gac_valid_dag(dag, x, y, z)
    expect("z_valid_at_truth", int(t_truth), row["z_valid_at_truth"])
    if mism:
        return {**res, "status": "rebuild_mismatch", "hard_failures": mism}

    hard: list[str] = []
    d0 = extension(g0)
    if not d0.is_dag():
        hard.append("D0 is not a DAG")
    if not (is_gac_valid_dag(d0, x, y, z) and is_valid_adjustment_set_gac_dag(d0, x, y, z)):
        hard.append("Z not valid in D0 although valid in G0")
    gd0 = DG(d0.nodes, d0.directed_edges)
    if not gac_valid(gd0, x, y, z):
        hard.append("own GAC says Z invalid in D0")
    pos = topo_order(gd0)
    y_desc_x = y in reach(gd0.ch, [x])

    def check(rec: dict, label: str):
        if rec.get("witness"):
            g1 = with_edges(d0, rec["witness"])
            if not g1.is_dag():
                hard.append(f"{label}: witness graph has a cycle")
            if is_gac_valid_dag(g1, x, y, z) or is_valid_adjustment_set_gac_dag(g1, x, y, z):
                hard.append(f"{label}: witness {rec['witness']} not confirmed by the reference oracles")
            if label == "order" and any(pos[a] > pos[b] for a, b in rec["witness"]):
                hard.append("order: witness breaks the causal order")
        for a, b in rec.pop("_sample", []):
            g1 = with_edges(d0, [(a, b)])
            if not is_valid_adjustment_set_gac_dag(g1, x, y, z):
                hard.append(f"{label}: own GAC valid but reference invalid after {a}->{b}")
        if rec.get("witness"):
            rec["witness"] = [list(e) for e in rec["witness"]]

    t_ord = tolerance(d0, x, y, z, "order", cfg["max_depth"], cfg["time_limit"], pos)
    check(t_ord, "order")
    if t_ord["t"] == 1:
        t_acy = {"mode": "acyclic", "t": 1, "exact": True, "lb": None, "implied_by_order": True}
    else:
        t_acy = tolerance(d0, x, y, z, "acyclic", cfg["max_depth"], cfg["time_limit_acyclic"])
        check(t_acy, "acyclic")

    # Sensitivity to the DAG chosen in G0's class: depth-1 order tolerance on seeded extensions.
    ext = []
    rng = np.random.default_rng(unit_seed(row_key(row)))
    for i in range(cfg["n_random_extensions"]):
        di = extension(g0, rng)
        ti = tolerance(di, x, y, z, "order", 1, cfg["time_limit_ext"])
        ti.pop("_sample", None)
        ext.append({"i": i, "same_as_d0": di == d0, "t1": ti["t"] == 1, "exact": ti["exact"],
                    "decided": bool(ti["exact"] or ti.get("stop") == "depth"),
                    "seconds": ti["seconds"]})

    res.update({
        "status": "ok", "z": sorted(z), "z_size": len(z), "n_k": len(k),
        "g0_undirected": len(g0.undirected_edges), "z_valid_at_truth": int(t_truth),
        "y_descendant_of_x_in_d0": y_desc_x, "d0_sha": sha(d0),
        "released_r_hop": row["r_hop"], "released_hop_exact": row["hop_exact"],
        "released_hop_lower_bound": row["hop_lower_bound"], "released_hop_method": row["hop_method"],
        "tol_order": t_ord, "tol_acyclic": t_acy, "extensions": ext, "hard_failures": hard,
    })
    return res


def run_unit(row: dict, cfg: dict) -> dict:
    t0 = time.perf_counter()
    rec = {"key": row_key(row), **{f: row[f] for f in ("network", "tier", "x", "y", "arm", "rep")},
           "experiment": EXPERIMENT, "criterion": CRITERION, "commit": cfg["commit"], "smoke": cfg["smoke"]}
    try:
        rec.update(compute(row, cfg))
    except Exception:  # reported, and retried on resume
        rec["status"] = "error"
        rec["error"] = traceback.format_exc()[-2000:]
    rec["seconds"] = round(time.perf_counter() - t0, 3)
    return rec


def run_task(task):
    row, cfg = task
    return run_unit(row, cfg)


# ------------------------------------------------------------------ driver


def select_rows(records: Path, smoke: bool) -> list[dict]:
    rows = [r for r in csv.DictReader(records.open()) if r["status"] == "ok"]
    # network-major order keeps one network per worker cache for long stretches
    rows.sort(key=row_key)
    if not smoke:
        return rows
    picked = []
    for net, arm in SMOKE_ROWS:
        hits = [r for r in rows if r["network"] == net and r["arm"] == arm]
        picked += hits[:SMOKE_PER_NETWORK]
    return picked


def done_keys(out: Path) -> set[str]:
    keys = set()
    for f in sorted(out.glob("rows*.jsonl")):
        for line in f.open():
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue  # a line cut by a kill; that unit is rerun
            if rec.get("status") != "error":
                keys.add(rec["key"])
    return keys


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--records", type=Path, default=DATA / "per_query_records.csv")
    ap.add_argument("--knowledge", type=Path, default=DATA / "knowledge_states.json")
    ap.add_argument("--models", type=Path, default=ROOT / "results/axisa3/networks/example_models")
    ap.add_argument("--out", type=Path, default=None, help="default results/rebuttal/E26 (E26_smoke with --smoke)")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--seed", type=int, default=RELEASE_SEED, help="D_RAND knowledge seed (released)")
    ap.add_argument("--max-depth", type=int, default=MAX_DEPTH)
    ap.add_argument("--time-limit", type=float, default=120.0, help="seconds per query, order tolerance")
    ap.add_argument("--time-limit-acyclic", type=float, default=60.0)
    ap.add_argument("--time-limit-ext", type=float, default=30.0, help="seconds per random extension")
    ap.add_argument("--n-random-extensions", type=int, default=N_RANDOM_EXTENSIONS)
    args = ap.parse_args(argv)
    if args.smoke:
        args.time_limit, args.time_limit_acyclic, args.time_limit_ext = 20.0, 10.0, 5.0
    out = args.out or ROOT / ("results/rebuttal/E26_smoke" if args.smoke else "results/rebuttal/E26")
    out.mkdir(parents=True, exist_ok=True)
    if not 0 <= args.shard < args.n_shards:
        ap.error("need 0 <= --shard < --n-shards")

    commit, src_dirty = repo_commit()
    cfg = {
        "experiment": EXPERIMENT, "criterion": CRITERION, "commit": commit, "src_dirty": src_dirty,
        "smoke": args.smoke, "seed": args.seed, "max_depth": args.max_depth,
        "time_limit": args.time_limit, "time_limit_acyclic": args.time_limit_acyclic,
        "time_limit_ext": args.time_limit_ext, "n_random_extensions": args.n_random_extensions,
        "models": str(args.models), "knowledge": str(args.knowledge),
        "records_sha256": file_sha(args.records), "knowledge_sha256": file_sha(args.knowledge),
        "shard": args.shard, "n_shards": args.n_shards, "workers": args.workers,
        "extension": "deterministic: smallest undirected edge, sorted direction, Meek closure, repeat",
        "order": "topological order of D0, ties broken by node name",
    }
    rows = select_rows(args.records, args.smoke)
    cfg["n_selected"] = len(rows)
    rows = [r for i, r in enumerate(rows) if i % args.n_shards == args.shard]
    suffix = "" if args.n_shards == 1 else f".shard{args.shard}of{args.n_shards}"
    (out / f"config{suffix}.json").write_text(json.dumps(cfg, indent=1, sort_keys=True))

    have = done_keys(out)
    todo = [(r, cfg) for r in rows if row_key(r) not in have]
    print(f"[E26] {len(rows)} units in this shard, {len(rows) - len(todo)} already done, {len(todo)} to run; "
          f"commit {commit[:12]}, src dirty {src_dirty}", flush=True)
    fh = (out / f"rows{suffix}.jsonl").open("a")
    state = {"done": 0, "hard": 0, "err": 0, "t0": time.perf_counter(), "stop": False}

    def heartbeat():
        while not state["stop"]:
            time.sleep(60)
            if not state["stop"]:
                print(f"[E26] {state['done']}/{len(todo)} done, hard {state['hard']}, errors {state['err']}, "
                      f"{time.perf_counter() - state['t0']:.0f} s", flush=True)

    threading.Thread(target=heartbeat, daemon=True).start()

    def consume(rec: dict):
        fh.write(json.dumps(rec, sort_keys=True) + "\n")
        fh.flush()
        state["done"] += 1
        if rec.get("status") == "error":
            state["err"] += 1
            print(f"[E26] ERROR {rec['key']}: {rec['error'].strip().splitlines()[-1]}", flush=True)
        if rec.get("hard_failures"):
            state["hard"] += 1
            print(f"[E26] CONSISTENCY {rec['key']}: {rec['hard_failures'][:3]}", flush=True)
        if args.smoke:
            to, ta = rec.get("tol_order", {}), rec.get("tol_acyclic", {})
            print(f"[E26] {rec['key']} {rec.get('status')} n={rec.get('n_nodes')} r={rec.get('released_r_hop')} "
                  f"t_order={to.get('t')} lb={to.get('lb')} checks={to.get('checks')} "
                  f"t_acyclic={ta.get('t')} {rec.get('seconds')} s", flush=True)

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
        fh.close()
    print(f"[E26] finished {state['done']} units in {time.perf_counter() - state['t0']:.0f} s; "
          f"hard {state['hard']}, errors {state['err']}", flush=True)
    if state["hard"]:
        print("[E26] CONSISTENCY FAILURE", flush=True)
        return 1
    return 2 if state["err"] else 0


if __name__ == "__main__":
    sys.exit(main())
