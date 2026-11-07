#!/usr/bin/env python3
"""E24: the real far-away discovery errors transplanted into the true graph.

For every learned CPDAG of the P-EST panel and every frame query (x, y) that is certifiable under
the repaired oracle arm of E5, build hybrid DAGs that carry the truth on pairs within d hops of
{x, y} and a consistent extension of the learned graph elsewhere, then score the committed set of
each hybrid's CPDAG (plus truthful repaired knowledge) in the true DAG. See PREREG.md.

    .venv/bin/python experiments/rebuttal/E24_e24_far_error_transplant/run.py --smoke
    .venv/bin/python experiments/rebuttal/E24_e24_far_error_transplant/run.py --workers 2

Exit 0 when every consistency check holds, 2 when one fails (rows are still written), 1 on error.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import sys
import threading
import time
from collections import deque
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
for _p in (REPO / "src", HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from bkrobust.benchmarks.describe import parse_file, to_mpdag  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.meek import apply_orientations, is_consistent_extension  # noqa: E402
from vendor_e5 import (  # noqa: E402
    canon,
    committed_adjustment_set,
    file_sha,
    graph_from_json,
    greedy_consistent,
    is_valid_adjustment_set_gac_dag,
    repo_commit,
    sha,
    truthful,
)

import hashlib  # noqa: E402

EXPERIMENT = "E24_e24_far_error_transplant"
NETWORKS = ("ecoli70", "arth150", "magic-niab", "magic-irri")
LEVELS = ("none", 1, 2, 3, "all")
CRITERION = "gac"
SMOKE_GRAPHS = {
    "magic-niab": [("pc", 5000, 0), ("ges", 1000, 0)],
    "magic-irri": [("pc", 5000, 0)],
}
FAR = 10**6  # distance of a node not connected to {x, y} in the true skeleton


def k_hash(k: list) -> str:
    """E5 run.py::k_hash."""
    return hashlib.sha256(json.dumps(sorted(k)).encode()).hexdigest()[:16]


# --- graph operations ---------------------------------------------------------------


def dor_tarsi(g: MPDAG) -> MPDAG | None:
    """A consistent DAG extension of a partially directed graph (Dor and Tarsi 1992), or None.

    Repeatedly removes the first node, in sorted order, that has no outgoing directed edge and
    whose every undirected neighbour is adjacent to all of its other adjacent nodes, orienting
    its undirected edges into it. Deterministic. The result keeps the skeleton and every
    directed edge, and creates no unshielded collider that the input did not have.
    """
    par = {v: set() for v in g.nodes}
    ch = {v: set() for v in g.nodes}
    nb = {v: set() for v in g.nodes}
    for a, b in g.directed_edges:
        ch[a].add(b)
        par[b].add(a)
    for a, b in g.undirected_edges:
        nb[a].add(b)
        nb[b].add(a)
    alive = set(g.nodes)
    out = set(g.directed_edges)
    while alive:
        pick = None
        for v in sorted(alive):
            if ch[v]:
                continue
            adj_v = par[v] | nb[v]
            if all((adj_v - {u}) <= (par[u] | ch[u] | nb[u] | {u}) for u in nb[v]):
                pick = v
                break
        if pick is None:
            return None
        v = pick
        for u in nb[v]:
            out.add((u, v))
            nb[u].discard(v)
        for p in par[v]:
            ch[p].discard(v)
        alive.discard(v)
        par[v], nb[v] = set(), set()
    d = MPDAG(g.nodes, out, [])
    return d if d.is_acyclic() else None


def hop_distances(dag: MPDAG, x: str, y: str) -> dict:
    """Hops from {x, y} in the true skeleton (multi-source BFS); unreachable nodes get FAR."""
    adj = {v: set() for v in dag.nodes}
    for a, b in dag.skeleton():
        adj[a].add(b)
        adj[b].add(a)
    dist = {x: 0, y: 0}
    q = deque([x, y])
    while q:
        u = q.popleft()
        for w in adj[u]:
            if w not in dist:
                dist[w] = dist[u] + 1
                q.append(w)
    return {v: dist.get(v, FAR) for v in dag.nodes}


def hybrid(dag: MPDAG, ext: MPDAG, inside: set) -> MPDAG | None:
    """Truth on pairs with both ends inside, the extension elsewhere; None when cyclic."""
    edges = {(a, b) for (a, b) in dag.directed_edges if a in inside and b in inside}
    edges |= {(a, b) for (a, b) in ext.directed_edges if not (a in inside and b in inside)}
    h = MPDAG(dag.nodes, edges, [])
    return h if h.is_acyclic() else None


def error_profile(dag: MPDAG, h: MPDAG, dist: dict) -> dict:
    """Pairs on which h differs from the truth, and how far the nearest one is from {x, y}."""
    ts, hs = dag.skeleton(), h.skeleton()
    missing, extra = ts - hs, hs - ts
    reversed_ = {canon(e) for e in h.directed_edges if (e[1], e[0]) in dag.directed_edges}
    wrong = missing | extra | reversed_

    def hop(e):
        return min(dist[e[0]], dist[e[1]])

    hops = [hop(e) for e in wrong]
    return {"err_missing": len(missing), "err_extra": len(extra), "err_reversed": len(reversed_),
            "err_total": len(wrong), "err_min_hop": min(hops) if hops else None,
            "err_within_1hop": sum(1 for v in hops if v <= 1)}


def repair_incremental(g: MPDAG, claims: list) -> tuple[MPDAG | None, list, int]:
    """E5's greedy repair, applied to the running graph instead of from scratch each time."""
    cur, kept, dropped = g, [], 0
    for e in claims:
        r = apply_orientations(cur, [e])
        if r is None:
            dropped += 1
        else:
            cur, kept = r, [*kept, e]
    return cur, kept, dropped


# --- inputs -------------------------------------------------------------------------


def load_inputs(cfg: dict) -> dict:
    data_dir, models_dir = Path(cfg["data_dir"]), Path(cfg["models_dir"])
    frame = [json.loads(line) for line in (data_dir / "frame.jsonl").open() if line.strip()]
    dags = {}
    for net in cfg["networks"]:
        path = models_dir / f"{net}.json"
        dags[net] = to_mpdag(parse_file(path, file_sha(path)))
    graphs: dict[str, dict] = {}
    for net in cfg["networks"]:
        graphs[f"{net}|true_cpdag|true|-"] = {"net": net, "method": "true_cpdag", "n": None, "seed": None,
                                              "kind": "true_cpdag", "g": dag_to_cpdag(dags[net])}
    for name in ("graphs.jsonl", "graphs_retry.jsonl"):
        p = data_dir / name
        if not p.exists():
            continue
        for line in p.open():
            if not line.strip():
                continue
            e = json.loads(line)
            if e["network"] not in cfg["networks"]:
                continue
            gid = f"{e['network']}|{e['method']}|{e['n']}|{e['seed']}"
            if gid in graphs:
                raise SystemExit(f"duplicate learned graph {gid}")
            graphs[gid] = {"net": e["network"], "method": e["method"], "n": e["n"], "seed": e["seed"],
                           "kind": "learned", "g": graph_from_json(e)}
    e5 = {}
    e5_dir = Path(cfg["e5_dir"])
    for p in sorted(e5_dir.glob("rows.*.jsonl")) if e5_dir.exists() else []:
        for line in p.open():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            if r.get("phase") == "B" and r.get("arm") == "A_ORACLE_R":
                e5[(r["gid"], r["x"], r["y"])] = r
    return {"frame": frame, "dags": dags, "graphs": graphs, "e5": e5}


# --- one graph ----------------------------------------------------------------------


def base_rec(cfg: dict, info: dict, gid: str, x: str, y: str) -> dict:
    return {"gid": gid, "network": info["net"], "method": info["method"], "n": info["n"],
            "seed": info["seed"], "graph_kind": info["kind"], "x": x, "y": y, "arm": "A_ORACLE_R",
            "criterion": CRITERION, "repo_commit": cfg["repo_commit"]}


def e5_check(inp: dict, gid: str, x: str, y: str, status: str, valid, kh: str) -> list[str]:
    st = inp["e5"].get((gid, x, y))
    if st is None:
        return [] if not inp["e5"] else ["e5_row_missing"]
    out = []
    if st.get("status") != status:
        out.append(f"e5_status:{status}!={st.get('status')}")
    if status == "ok" and int(st.get("z_valid_at_truth")) != int(valid):
        out.append(f"e5_valid:{valid}!={st.get('z_valid_at_truth')}")
    if st.get("k_sha256") != kh:
        out.append(f"e5_k_sha:{kh}!={st.get('k_sha256')}")
    return out


def run_block(block: dict) -> list[dict]:
    cfg, inp = _CFG, _INP
    gid = block["gid"]
    info = inp["graphs"][gid]
    dag, g = inp["dags"][info["net"]], info["g"]
    t0 = time.perf_counter()
    oracle = truthful(dag, sorted(g.undirected_edges))
    kept, dropped = greedy_consistent(g, oracle)
    g0 = apply_orientations(g, kept)
    kh = k_hash(kept)
    rows = []

    if info["kind"] == "true_cpdag":  # the true-class anchor
        for x, y in block["pairs"]:
            rec = base_rec(cfg, info, gid, x, y)
            rec.update(level="anchor", unit=f"anchor|{gid}|{x}|{y}", k_sha256=kh, n_k=len(kept))
            fails = []
            cs = committed_adjustment_set(g0, x, y)
            rec["status"] = "ok" if cs.z is not None else cs.verdict
            if cs.z is not None:
                rec["z_size"] = len(cs.z)
                rec["z_valid_at_truth"] = int(is_valid_adjustment_set_gac_dag(dag, x, y, cs.z))
                if not rec["z_valid_at_truth"]:
                    fails.append("anchor_set_invalid_at_truth")
            fails += e5_check(inp, gid, x, y, rec["status"], rec.get("z_valid_at_truth"), kh)
            rec["consistency_fail"] = fails
            rows.append(rec)
        return rows

    # one consistent extension of G0, shared by every query of this graph
    ext = dor_tarsi(g0) if g0 is not None else None
    ext_ok = ext is not None and is_consistent_extension(ext, g0)
    ext_fail = [] if ext is None or ext_ok else ["extension_not_consistent"]
    first = True
    for x, y in block["pairs"]:
        rec = base_rec(cfg, info, gid, x, y)
        cs0 = committed_adjustment_set(g0, x, y) if g0 is not None else None
        st0 = "knowledge_inconsistent" if cs0 is None else ("ok" if cs0.z is not None else cs0.verdict)
        v0 = int(is_valid_adjustment_set_gac_dag(dag, x, y, cs0.z)) if st0 == "ok" else None
        rec.update(level="g0", unit=f"g0|{gid}|{x}|{y}", status=st0, z_valid_at_truth=v0, k_sha256=kh,
                   n_k=len(kept), n_oracle_dropped=dropped, has_extension=int(ext_ok),
                   consistency_fail=e5_check(inp, gid, x, y, st0, v0, kh) + ext_fail)
        rows.append(rec)
        if st0 != "ok":
            continue
        dist = hop_distances(dag, x, y)
        for level in LEVELS:
            t1 = time.perf_counter()
            r = base_rec(cfg, info, gid, x, y)
            r.update(level=str(level), unit=f"h|{gid}|{x}|{y}|{level}", g0_valid_at_truth=v0)
            fails: list[str] = []
            if not ext_ok:
                r["hstatus"] = "no_extension"
                r["consistency_fail"] = ext_fail
                rows.append(r)
                continue
            if level == "none":
                inside: set = set()
            elif level == "all":
                inside = set(dag.nodes)
            else:
                inside = {v for v, dv in dist.items() if dv <= level}
            r["n_inside"] = len(inside)
            h = hybrid(dag, ext, inside)
            if h is None:
                r["hstatus"] = "cyclic"
                r["consistency_fail"] = fails
                rows.append(r)
                continue
            # check 3: truth inside N_d
            for a, b in dag.directed_edges:
                if a in inside and b in inside and not h.is_directed_edge(a, b):
                    fails.append("hybrid_differs_from_truth_inside")
                    break
            for a, b in h.directed_edges:
                if a in inside and b in inside and not dag.is_directed_edge(a, b):
                    fails.append("hybrid_extra_inside")
                    break
            if level == "all" and h != dag:
                fails.append("all_level_hybrid_not_truth")
            r.update(error_profile(dag, h, dist))
            r["hybrid_sha"] = sha(h)
            hcp = dag_to_cpdag(h)
            r["hcp_undirected"] = len(hcp.undirected_edges)
            oh = truthful(dag, sorted(hcp.undirected_edges))
            g0h, kept_h, drop_h = repair_incremental(hcp, oh)
            if first:  # check 6: the incremental repair equals E5's repair
                kept_e5, _ = greedy_consistent(hcp, oh)
                if kept_e5 != kept_h or apply_orientations(hcp, kept_e5) != g0h:
                    fails.append("incremental_repair_differs_from_e5")
                r["repair_crosschecked"] = 1
            r.update(n_k=len(kept_h), n_k_dropped=drop_h, g0h_undirected=len(g0h.undirected_edges))
            cs = committed_adjustment_set(g0h, x, y)
            r["committed_verdict"] = cs.verdict
            if cs.z is None:
                r["hstatus"] = cs.verdict
            else:
                r["hstatus"] = "ok"
                r["z_size"] = len(cs.z)
                ok = is_valid_adjustment_set_gac_dag(dag, x, y, cs.z)
                r["z_valid_at_truth"] = int(ok)
                r["truth_extends_g0h"] = int(is_consistent_extension(dag, g0h))
                if r["truth_extends_g0h"] and not ok:
                    fails.append("truth_extends_g0h_but_set_invalid")
                if level == "all" and not ok:
                    fails.append("all_level_set_invalid")
            r["seconds"] = round(time.perf_counter() - t1, 3)
            r["consistency_fail"] = fails
            rows.append(r)
        first = False
    if rows:
        rows[0]["block_seconds"] = round(time.perf_counter() - t0, 3)
    return rows


_INP: dict = {}
_CFG: dict = {}


def init_worker(cfg: dict) -> None:
    _CFG.update(cfg)
    _INP.update(load_inputs(cfg))


def make_blocks(inp: dict, args) -> list[dict]:
    pairs_by_net: dict[str, list] = {}
    for r in inp["frame"]:
        pairs_by_net.setdefault(r["network"], []).append((r["x"], r["y"]))
    blocks = []
    for gid, info in sorted(inp["graphs"].items()):
        if args.smoke and info["kind"] == "learned" and \
                (info["method"], info["n"], info["seed"]) not in SMOKE_GRAPHS.get(info["net"], []):
            continue
        blocks.append({"gid": gid, "pairs": pairs_by_net.get(info["net"], [])})
    size = {n: len(inp["dags"][n].nodes) for n in inp["dags"]}
    mine = [b for i, b in enumerate(blocks) if i % args.n_shards == args.shard]
    mine.sort(key=lambda b: (-size[inp["graphs"][b["gid"]]["net"]], b["gid"]))
    return mine


def block_units(inp: dict, b: dict) -> list[str]:
    """The g0/anchor units of a block; hybrid units exist only for certifiable pairs."""
    kind = inp["graphs"][b["gid"]]["kind"]
    pre = "anchor" if kind == "true_cpdag" else "g0"
    return [f"{pre}|{b['gid']}|{x}|{y}" for x, y in b["pairs"]]


def done_units(out: Path) -> tuple[set, list]:
    done, fails = set(), []
    for p in sorted(out.glob("rows.*.jsonl")):
        for line in p.open():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            done.add(r["unit"])
            if r.get("consistency_fail"):
                fails.append((r["unit"], r["consistency_fail"]))
    return done, fails


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--smoke", action="store_true", help="two networks, three learned graphs, the anchor")
    ap.add_argument("--out", default="results/rebuttal/E24")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--data-dir", default="experiments/rebuttal/_data/pest")
    ap.add_argument("--models-dir", default="results/axisa3/networks/example_models")
    ap.add_argument("--e5-dir", default="results/rebuttal/E5", help="E5 rows, for consistency check 1")
    ap.add_argument("--networks", default=",".join(NETWORKS))
    args = ap.parse_args(argv)
    if args.workers > 6:
        ap.error("at most 6 workers")
    if not 0 <= args.shard < args.n_shards:
        ap.error("--shard must be in [0, --n-shards)")

    def rel(p):
        p = Path(p)
        return p if p.is_absolute() else REPO / p

    out = rel(args.out)
    if args.smoke and args.out == "results/rebuttal/E24":
        out = REPO / "results/rebuttal/E24_smoke"
    networks = list(SMOKE_GRAPHS) if args.smoke else [n for n in args.networks.split(",") if n]
    out.mkdir(parents=True, exist_ok=True)
    data_dir = rel(args.data_dir)
    e5_dir = rel(args.e5_dir)
    e5_files = {p.name: file_sha(p) for p in sorted(e5_dir.glob("rows.*.jsonl"))} if e5_dir.exists() else {}
    cfg = {
        "experiment": EXPERIMENT, "criterion": CRITERION, "repo_commit": repo_commit(REPO),
        "run_sha256": file_sha(Path(__file__)), "vendor_sha256": file_sha(HERE / "vendor_e5.py"),
        "smoke": bool(args.smoke), "networks": networks, "levels": [str(v) for v in LEVELS],
        "arm": "A_ORACLE_R", "extension": "dor_tarsi_sorted",
        "data_dir": str(data_dir), "models_dir": str(rel(args.models_dir)), "e5_dir": str(e5_dir),
        "input_sha256": {p.name.replace(".", "_"): file_sha(p) for p in sorted(data_dir.iterdir()) if p.is_file()},
        "e5_rows_sha256": e5_files, "shard": args.shard, "n_shards": args.n_shards, "workers": args.workers,
    }
    tag = f"shard{args.shard}of{args.n_shards}"
    (out / f"config.{tag}.json").write_text(json.dumps(cfg, indent=1))
    inp = load_inputs(cfg)
    if not inp["e5"]:
        print("warning: no E5 rows found; consistency check 1 is skipped", flush=True)
    blocks = make_blocks(inp, args)
    done, _ = done_units(out)
    todo = []
    for b in blocks:
        left = [p for p, u in zip(b["pairs"], block_units(inp, b)) if u not in done]
        if left:
            todo.append({**b, "pairs": left})
    n_pairs = sum(len(b["pairs"]) for b in todo)
    print(f"E24 {tag}: {len(blocks)} graphs, {len(todo)} to run ({n_pairs} queries), {args.workers} workers, "
          f"criterion {CRITERION}, commit {cfg['repo_commit'][:12]}", flush=True)

    t0 = time.perf_counter()
    written = 0
    stop = threading.Event()

    def heartbeat():
        while not stop.wait(60):
            print(f"  [alive] {written} rows, {time.perf_counter() - t0:.0f} s", flush=True)

    threading.Thread(target=heartbeat, daemon=True).start()
    fh = (out / f"rows.{tag}.jsonl").open("a")

    def write(rows):
        # g0/anchor rows last, so a crash mid-block leaves the block to be redone, never half-marked
        nonlocal written
        for r in sorted(rows, key=lambda r: r["level"] in ("g0", "anchor")):
            fh.write(json.dumps(r, default=str) + "\n")
        fh.flush()
        written += len(rows)

    if args.workers <= 1:
        init_worker(cfg)
        for i, b in enumerate(todo, 1):
            write(run_block(b))
            print(f"  {i}/{len(todo)} graphs, {written} rows, {time.perf_counter() - t0:.0f} s", flush=True)
    else:
        with mp.get_context("spawn").Pool(args.workers, initializer=init_worker, initargs=(cfg,)) as pool:
            for i, rows in enumerate(pool.imap_unordered(run_block, todo), 1):
                write(rows)
                print(f"  {i}/{len(todo)} graphs, {written} rows, {time.perf_counter() - t0:.0f} s", flush=True)
    fh.close()
    stop.set()

    _, fails = done_units(out)
    print(f"done: {written} rows written in {time.perf_counter() - t0:.0f} s; "
          f"{len(fails)} rows with a consistency failure", flush=True)
    if fails:
        (out / f"consistency_failures.{tag}.json").write_text(json.dumps(fails, indent=1))
        for u, f in fails[:20]:
            print(f"  FAIL {u}: {f}", flush=True)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
