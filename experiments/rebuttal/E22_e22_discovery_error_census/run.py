#!/usr/bin/env python3
"""E22: which discovery error breaks each committed set, a census on the stored learned graphs.

Reads the certifiable A_ORACLE_R rows of E5 on learned graphs, rebuilds G0 and the committed set Z with
E5's functions (and checks them against the stored row), diffs each learned CPDAG against the true CPDAG,
names the failing adjustment condition at the truth, and tests the hypothesis of statement S3 restricted to
the true ancestral closure of {x, y} and Z. See PREREG.md.

    .venv/bin/python experiments/rebuttal/E22_e22_discovery_error_census/run.py --smoke
    .venv/bin/python experiments/rebuttal/E22_e22_discovery_error_census/run.py --workers 2

Exit 0 when every consistency check holds, 2 when one fails (a reconstruction mismatch stops the run at
once), 1 on error.
"""

from __future__ import annotations

import argparse
import collections
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
for _p in (REPO / "src", HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import networkx as nx  # noqa: E402

from bkrobust.benchmarks.describe import parse_file, to_mpdag  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.meek import apply_orientations, is_consistent_extension  # noqa: E402
from bkrobust.mpdag_criterion.paths import possible_descendants  # noqa: E402
from vendor_bkrobust import (  # noqa: E402
    _d_separated,
    committed_adjustment_set,
    is_valid_adjustment_set_gac_dag,
)
from vendor_e5 import canon, file_sha, k_hash, load_learned, oracle_r_knowledge, repo_commit  # noqa: E402

EXPERIMENT = "E22_e22_discovery_error_census"
CRITERION = "gac"
ARM = "A_ORACLE_R"
NETWORKS = ("ecoli70", "arth150", "magic-niab", "magic-irri")
# tie order for the nearest error of a row; also the order of the census table
ERROR_TYPES = ("query_direction_wrong", "missing_adjacency", "reversed_compelled", "false_v_structure",
               "missed_v_structure", "extra_adjacency", "under_oriented", "over_oriented")
SMOKE_PER_NET = {"ecoli70": 4, "magic-niab": 4, "arth150": 4}  # rows from the first graph of each


# --- graph helpers -------------------------------------------------------------------


def v_structures(g: MPDAG) -> set[tuple[str, str, str]]:
    """(a, b, c) with a -> b <- c, a < c, a and c non-adjacent."""
    out = set()
    for b in g.nodes:
        ps = sorted(g.parents(b))
        for i, a in enumerate(ps):
            for c in ps[i + 1:]:
                if not g.has_edge(a, c):
                    out.add((a, b, c))
    return out


def digraph(dag: MPDAG) -> nx.DiGraph:
    d = nx.DiGraph()
    d.add_nodes_from(dag.nodes)
    d.add_edges_from(dag.directed_edges)
    return d


def hop_from(dag: MPDAG, sources) -> dict[str, int]:
    """Shortest-path distance in the true skeleton from the nearest source."""
    adj = collections.defaultdict(set)
    for a, b in dag.skeleton():
        adj[a].add(b)
        adj[b].add(a)
    dist = {s: 0 for s in sources}
    frontier = list(sources)
    while frontier:
        nxt = []
        for u in frontier:
            for v in adj[u]:
                if v not in dist:
                    dist[v] = dist[u] + 1
                    nxt.append(v)
        frontier = nxt
    return dist


def graph_errors(learned: MPDAG, true_cp: MPDAG) -> list[dict]:
    """Net difference between the learned CPDAG and the true CPDAG, one dict per error (no row fields)."""
    errs = []
    ls, ts = learned.skeleton(), true_cp.skeleton()
    for e in sorted(ts - ls):
        errs.append({"type": "missing_adjacency", "nodes": list(e), "detail": f"{e[0]}-{e[1]}"})
    for e in sorted(ls - ts):
        errs.append({"type": "extra_adjacency", "nodes": list(e), "detail": f"{e[0]}-{e[1]}"})
    for a, b in sorted(ls & ts):
        l_ab, l_ba, l_u = learned.is_directed_edge(a, b), learned.is_directed_edge(b, a), learned.is_undirected_edge(a, b)
        t_ab, t_ba, t_u = true_cp.is_directed_edge(a, b), true_cp.is_directed_edge(b, a), true_cp.is_undirected_edge(a, b)
        if (t_ab and l_ba) or (t_ba and l_ab):
            tt, th = (a, b) if t_ab else (b, a)
            errs.append({"type": "reversed_compelled", "nodes": [a, b], "detail": f"true {tt}->{th}, learned {th}->{tt}"})
        elif (t_ab or t_ba) and l_u:
            tt, th = (a, b) if t_ab else (b, a)
            errs.append({"type": "under_oriented", "nodes": [a, b], "detail": f"true {tt}->{th}, learned {a}-{b}"})
        elif t_u and (l_ab or l_ba):
            lt, lh = (a, b) if l_ab else (b, a)
            errs.append({"type": "over_oriented", "nodes": [a, b], "detail": f"true {a}-{b}, learned {lt}->{lh}"})
    lv, tv = v_structures(learned), v_structures(true_cp)
    for t in sorted(lv - tv):
        errs.append({"type": "false_v_structure", "nodes": list(t), "detail": f"{t[0]}->{t[1]}<-{t[2]}"})
    for t in sorted(tv - lv):
        errs.append({"type": "missed_v_structure", "nodes": list(t), "detail": f"{t[0]}->{t[1]}<-{t[2]}"})
    return errs


def failing_conditions(d: nx.DiGraph, x: str, y: str, z: frozenset) -> dict:
    """The three GAC conditions at the truth, in the order of PREREG.md."""
    desc_x = nx.descendants(d, x)
    y_desc = y in desc_x
    causal = desc_x & (nx.ancestors(d, y) | {y})
    forbidden = {x}
    for c in causal:
        forbidden |= nx.descendants(d, c) | {c}
    forb_in_z = sorted(set(z) & forbidden)
    proper = d.copy()
    proper.remove_edges_from([(x, c) for c in list(d.successors(x)) if c in causal])
    open_path = not _d_separated(proper, {x}, {y}, set(z))
    first = None
    if not forb_in_z and not open_path:
        first = None
    elif not y_desc:
        first = "y_not_true_descendant"
    elif forb_in_z:
        first = "forbidden_in_z"
    else:
        first = "open_noncausal_path"
    return {"y_true_descendant": int(y_desc), "forbidden_in_z": forb_in_z, "open_noncausal_path": int(open_path),
            "invalid_by_conditions": int(bool(forb_in_z) or open_path), "first_failing": first}


def one_extension(g: MPDAG) -> MPDAG | None:
    """A DAG in [g]: orient the first open edge, close, repeat."""
    cur = g
    while cur.undirected_edges:
        a, b = sorted(cur.undirected_edges)[0]
        cur = apply_orientations(cur, [(a, b)])
        if cur is None:
            return None
    return cur


def s3_local(dag: MPDAG, g0: MPDAG, A: set, x: str, y: str, z: frozenset, build_extension: bool) -> dict:
    """The hypothesis of S3 restricted to A, comparing G0 with the true DAG; plus a built witness D."""
    true_in_a = [(a, b) for (a, b) in dag.directed_edges if a in A and b in A]
    missing = [f"{a}->{b}" for a, b in true_in_a if not g0.has_edge(a, b)]
    reversed_ = [f"{a}->{b}" for a, b in true_in_a if g0.is_directed_edge(b, a)]
    extra = [f"{a}-{b}" for a, b in sorted(g0.skeleton()) if a in A and b in A and not dag.has_edge(a, b)]
    to_orient = [(a, b) for a, b in true_in_a if g0.is_undirected_edge(a, b)]
    out = {"s3_missing_in_A": missing, "s3_reversed_in_A": reversed_, "s3_n_extra_in_A": len(extra),
           "s3_n_to_orient": len(to_orient), "s3_closure_ok": None, "s3_hypothesis": 0,
           "s3_extension_built": None, "s3_extension_contains_truth_on_A": None, "s3_z_valid_in_extension": None}
    if missing or reversed_:
        return out
    g1 = apply_orientations(g0, to_orient)
    out["s3_closure_ok"] = int(g1 is not None)
    if g1 is None:
        return out
    out["s3_hypothesis"] = 1
    if build_extension:
        d = one_extension(g1)
        ok = d is not None and is_consistent_extension(d, g1)
        out["s3_extension_built"] = int(ok)
        if ok:
            out["s3_extension_contains_truth_on_A"] = int(all(d.is_directed_edge(a, b) for a, b in true_in_a))
            out["s3_z_valid_in_extension"] = int(is_valid_adjustment_set_gac_dag(d, x, y, z))
    return out


# --- workers -------------------------------------------------------------------------

_W: dict = {}


def init_worker(cfg: dict) -> None:
    _W["cfg"] = cfg
    models = Path(cfg["models_dir"])
    _W["dags"], _W["true_cp"] = {}, {}
    for net in cfg["networks"]:
        p = models / f"{net}.json"
        dag = to_mpdag(parse_file(p, file_sha(p)))
        _W["dags"][net] = dag
        _W["true_cp"][net] = dag_to_cpdag(dag)
    _W["graphs"] = load_learned(Path(cfg["data_dir"]), tuple(cfg["networks"]))


def run_block(block: dict) -> list[dict]:
    cfg = _W["cfg"]
    info = _W["graphs"][block["gid"]]
    net, g = info["net"], info["g"]
    dag, true_cp = _W["dags"][net], _W["true_cp"][net]
    d = digraph(dag)
    k, n_dropped = oracle_r_knowledge(dag, g)
    g0 = apply_orientations(g, k)
    base_errs = graph_errors(g, true_cp)
    rows = []
    for st in block["rows"]:
        x, y = st["x"], st["y"]
        t0 = time.perf_counter()
        rec = {"unit": st["unit"], "gid": block["gid"], "network": net, "method": info["method"], "n": info["n"],
               "seed": info["seed"], "x": x, "y": y, "arm": ARM, "criterion": CRITERION,
               "repo_commit": cfg["repo_commit"], "n_oracle_dropped": n_dropped, "n_k": len(k)}
        fails = []
        if g0 is None:
            rec["consistency_fail"] = ["g0_rebuild_failed"]
            rows.append(rec)
            continue
        cs = committed_adjustment_set(g0, x, y)
        z = cs.z
        got = {"k_sha256": k_hash(k), "n_k": len(k), "committed_verdict": cs.verdict,
               "z_size": None if z is None else len(z), "y_true_descendant": int(y in dag.descendants(x))}
        for f, v in got.items():
            if v != st[f]:
                fails.append(f"reconstruct:{f}:{v!r}!={st[f]!r}")
        if z is None:
            rec["consistency_fail"] = fails or ["reconstruct:no_committed_set"]
            rows.append(rec)
            continue
        z_valid = int(is_valid_adjustment_set_gac_dag(dag, x, y, z))
        if z_valid != st["z_valid_at_truth"]:
            fails.append(f"reconstruct:z_valid_at_truth:{z_valid}!={st['z_valid_at_truth']}")
        rec["reconstruct_ok"] = int(not fails)
        rec.update(z=sorted(z), z_size=len(z), committed_verdict=cs.verdict, z_valid_at_truth=z_valid)
        cond = failing_conditions(d, x, y, z)
        rec.update(cond)
        if cond["invalid_by_conditions"] != 1 - z_valid:
            fails.append("conditions_disagree_with_gac")
        # the true ancestral closure of {x, y} and Z, inclusive
        A = {x, y} | set(z)
        for v in list(A):
            A |= nx.ancestors(d, v)
        rec["size_A"] = len(A)
        hop = hop_from(dag, [x, y])
        errs = []
        for e in base_errs:
            hs = [hop[v] for v in e["nodes"] if v in hop]
            errs.append({"type": e["type"], "nodes": e["nodes"], "detail": e["detail"],
                         "hop": min(hs) if hs else None, "in_A": int(all(v in A for v in e["nodes"]))})
        y_possdesc = y in possible_descendants(g0, x)
        rec["y_possible_descendant_in_g0"] = int(y_possdesc)
        if y_possdesc and not cond["y_true_descendant"]:
            errs.append({"type": "query_direction_wrong", "nodes": [x, y], "detail": f"{y} not a descendant of {x}",
                         "hop": 0, "in_A": 1})
        rec["errors"] = errs
        rec["n_errors"] = {t: sum(e["type"] == t for e in errs) for t in ERROR_TYPES}
        rec["n_errors_within2"] = {t: sum(e["type"] == t and e["hop"] is not None and e["hop"] <= 2 for e in errs)
                                   for t in ERROR_TYPES}
        rec["n_errors_in_A"] = {t: sum(e["type"] == t and e["in_A"] for e in errs) for t in ERROR_TYPES}
        reach = [e for e in errs if e["hop"] is not None]
        if reach:
            mh = min(e["hop"] for e in reach)
            rec["nearest_hop"] = mh
            rec["nearest_type"] = next(t for t in ERROR_TYPES if any(e["type"] == t and e["hop"] == mh for e in reach))
        else:
            rec["nearest_hop"], rec["nearest_type"] = None, "none"
        rec["any_error_within2"] = int(any(e["hop"] is not None and e["hop"] <= 2 for e in errs))
        in_a = [e for e in errs if e["in_A"]]
        rec["census_only_extra_in_A"] = int(all(e["type"] == "extra_adjacency" for e in in_a))
        rec["census_n_in_A"] = len(in_a)
        s3 = s3_local(dag, g0, A, x, y, z, cfg["build_extension"])
        rec.update(s3)
        if s3["s3_hypothesis"] and cfg["build_extension"]:
            if s3["s3_extension_built"] != 1:
                fails.append("s3_extension_not_built")
            elif s3["s3_extension_contains_truth_on_A"] != 1:
                fails.append("s3_extension_misses_truth_on_A")
            elif s3["s3_z_valid_in_extension"] != 1:
                fails.append("s3_z_invalid_in_extension_of_g0")
        rec["seconds"] = round(time.perf_counter() - t0, 3)
        rec["consistency_fail"] = fails
        rows.append(rec)
    return rows


# --- planning ------------------------------------------------------------------------


def select_rows(e5_rows: Path, networks) -> dict[str, list[dict]]:
    keep = ("unit", "gid", "x", "y", "k_sha256", "n_k", "committed_verdict", "z_size", "y_true_descendant",
            "z_valid_at_truth")
    by_gid: dict[str, list[dict]] = collections.defaultdict(list)
    seen = set()
    for line in e5_rows.open():
        r = json.loads(line)
        if (r["phase"] == "B" and r["graph_kind"] == "learned" and r["arm"] == ARM and r["status"] == "ok"
                and r["network"] in networks and r["unit"] not in seen):
            seen.add(r["unit"])
            by_gid[r["gid"]].append({f: r.get(f) for f in keep})
    return dict(by_gid)


def done_units(out: Path) -> set:
    done = set()
    for p in sorted(out.glob("rows.*.jsonl")):
        for line in p.open():
            try:
                done.add(json.loads(line)["unit"])
            except (json.JSONDecodeError, KeyError):
                continue
    return done


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--smoke", action="store_true", help="four rows from one graph of three networks")
    ap.add_argument("--out", default="results/rebuttal/E22")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--e5-rows", default="results/rebuttal/E5/rows.shard0of1.jsonl")
    ap.add_argument("--data-dir", default="experiments/rebuttal/_data/pest")
    ap.add_argument("--models-dir", default="results/axisa3/networks/example_models")
    ap.add_argument("--no-extension", action="store_true", help="skip building the witness DAG D")
    args = ap.parse_args(argv)
    if args.workers > 6:
        ap.error("at most 6 workers")
    if not 0 <= args.shard < args.n_shards:
        ap.error("--shard must be in [0, --n-shards)")

    def rel(p):
        p = Path(p)
        return p if p.is_absolute() else REPO / p

    e5_rows, data_dir, models_dir, out = rel(args.e5_rows), rel(args.data_dir), rel(args.models_dir), rel(args.out)
    if args.smoke and args.out == "results/rebuttal/E22":
        out = REPO / "results/rebuttal/E22_smoke"
    out.mkdir(parents=True, exist_ok=True)

    by_gid = select_rows(e5_rows, NETWORKS)
    n_all = sum(len(v) for v in by_gid.values())
    n_inv = sum(1 for v in by_gid.values() for r in v if r["z_valid_at_truth"] == 0)
    if args.smoke:
        picked = {}
        for net, m in SMOKE_PER_NET.items():
            gid = next(g for g in sorted(by_gid) if g.startswith(net + "|")
                       and {r["z_valid_at_truth"] for r in by_gid[g]} == {0, 1})
            rows = sorted(by_gid[gid], key=lambda r: r["unit"])
            inv = [r for r in rows if r["z_valid_at_truth"] == 0][: m // 2]
            val = [r for r in rows if r["z_valid_at_truth"] == 1][: m - len(inv)]
            picked[gid] = inv + val
        by_gid = picked

    inputs = {p.name.replace(".", "_"): file_sha(p) for p in sorted(data_dir.iterdir())
              if p.is_file() and p.name in ("graphs.jsonl", "graphs_retry.jsonl")}
    inputs["e5_rows"] = file_sha(e5_rows)
    for net in NETWORKS:
        inputs[f"model_{net.replace('-', '_')}"] = file_sha(models_dir / f"{net}.json")
    cfg = {
        "experiment": EXPERIMENT, "criterion": CRITERION, "arm": ARM, "repo_commit": repo_commit(REPO),
        "smoke": bool(args.smoke), "networks": list(NETWORKS), "e5_rows": str(e5_rows),
        "n_rows_population": n_all, "n_invalid_population": n_inv, "n_graphs_population": len(select_rows(e5_rows, NETWORKS)),
        "build_extension": not args.no_extension, "data_dir": str(data_dir), "models_dir": str(models_dir),
        "input_sha256": inputs, "run_sha256": file_sha(Path(__file__)),
        "vendor_bkrobust_sha256": file_sha(HERE / "vendor_bkrobust.py"), "vendor_e5_sha256": file_sha(HERE / "vendor_e5.py"),
        "shard": args.shard, "n_shards": args.n_shards, "workers": args.workers,
    }
    tag = f"shard{args.shard}of{args.n_shards}"
    (out / f"config.{tag}.json").write_text(json.dumps(cfg, indent=1))

    done = done_units(out)
    gids = sorted(by_gid)
    mine = [g for i, g in enumerate(gids) if i % args.n_shards == args.shard]
    todo = []
    for g in mine:
        left = [r for r in by_gid[g] if r["unit"] not in done]
        if left:
            todo.append({"gid": g, "rows": left})
    n_rows = sum(len(b["rows"]) for b in todo)
    print(f"E22 {tag}: population {n_all} rows ({n_inv} invalid) on {cfg['n_graphs_population']} graphs; "
          f"{len(todo)} graphs to run ({n_rows} rows), {args.workers} workers, criterion {CRITERION}, "
          f"commit {cfg['repo_commit'][:12]}", flush=True)

    t0 = time.perf_counter()
    written, fails = 0, []
    fh = (out / f"rows.{tag}.jsonl").open("a")
    stopped = False

    def write(rows) -> bool:
        nonlocal written
        bad = False
        for r in rows:
            fh.write(json.dumps(r, default=str) + "\n")
            if r.get("consistency_fail"):
                fails.append((r["unit"], r["consistency_fail"]))
                bad = bad or any(f.startswith("reconstruct") or f == "g0_rebuild_failed" for f in r["consistency_fail"])
        fh.flush()
        written += len(rows)
        return bad

    if args.workers <= 1:
        init_worker(cfg)
        for i, b in enumerate(todo, 1):
            if write(run_block(b)):
                stopped = True
                break
            print(f"  {i}/{len(todo)} graphs, {written}/{n_rows} rows, {time.perf_counter() - t0:.0f} s", flush=True)
    else:
        with mp.get_context("spawn").Pool(args.workers, initializer=init_worker, initargs=(cfg,)) as pool:
            for i, rows in enumerate(pool.imap_unordered(run_block, todo), 1):
                if write(rows):
                    stopped = True
                    pool.terminate()
                    break
                print(f"  {i}/{len(todo)} graphs, {written}/{n_rows} rows, {time.perf_counter() - t0:.0f} s", flush=True)
    fh.close()
    print(f"done: {written} rows in {time.perf_counter() - t0:.0f} s; {len(fails)} rows with a consistency failure"
          + ("; STOPPED on a reconstruction mismatch" if stopped else ""), flush=True)
    if fails:
        (out / f"consistency_failures.{tag}.json").write_text(json.dumps(fails, indent=1))
        for u, f in fails[:20]:
            print(f"  FAIL {u}: {f}", flush=True)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
