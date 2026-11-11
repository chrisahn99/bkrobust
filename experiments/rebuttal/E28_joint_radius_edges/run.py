#!/usr/bin/env python3
"""E28: joint radius over orientation retractions and edge additions near the query (see PREREG.md).

For every scored E22 row (learned PC/GES graph, A_ORACLE_R knowledge K, committed set Z valid in G0),
rebuilds G (the learned CPDAG), K and Z with E5's functions and checks them against the stored row. Then:

* r_orient: the orientation-only radius of Z, the ledger's ``breaking_depths`` (retractions of K only,
  depth 3, 300-subset budget), as E5 and the paper compute it;
* r_joint: the fewest moves after which Z fails the GAC in the resulting state, where a move is either
  the retraction of one orientation of K or the addition of one directed edge a -> b between two nodes
  that are non-adjacent in G and both within two hops of {x, y} u Z in the skeleton of G. The state
  after retractions R and additions A is meek_closure(G + A), then ``apply_orientations`` with K minus R
  (the radius's own primitive). An inconsistent state (closure FAIL) is not a break.
  Depth 1 is exhaustive and lists every breaking move (the witnesses); depth 2 is searched only when
  depth 1 finds no break, in the order retraction pairs, retraction + addition, addition pairs, and
  stops at the first break or at its time limit;
* coverage: for invalid rows, whether a true edge that G0 misses inside An_true({x, y} u Z) (E22's
  ``s3_missing_in_A``) is a depth-1 addition witness, in either direction and in its true direction.

One JSON line per row, flushed as it goes; a restarted shard skips rows already written.
Exit 0 when every consistency check holds, 2 when one fails, 1 on error.

    .venv/bin/python experiments/rebuttal/E28_joint_radius_edges/run.py --smoke
    .venv/bin/python experiments/rebuttal/E28_joint_radius_edges/run.py --shard 0 --n-shards 6
"""

from __future__ import annotations

import argparse
import collections
import itertools
import json
import os
import resource
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

from bkrobust.benchmarks.describe import parse_file, to_mpdag  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.meek import apply_orientations, meek_closure  # noqa: E402
from bkrobust.gac.mpdag_level import clear_cache, is_gac_valid_mpdag, why_invalid_gac  # noqa: E402
from vendor_bkrobust import breaking_depths, committed_adjustment_set, is_valid_adjustment_set_gac_dag  # noqa: E402
from vendor_e5 import file_sha, k_hash, load_learned, oracle_r_knowledge, repo_commit  # noqa: E402

EXPERIMENT = "E28_joint_radius_edges"
CRITERION = "gac"
ARM = "A_ORACLE_R"
NETWORKS = ("ecoli70", "arth150", "magic-niab", "magic-irri")
HOPS = 2
MAX_JOINT_DEPTH = 2
DEPTH1_LIMIT_S = 1800.0
DEPTH2_LIMIT_S = 600.0
N_STATE_CHECKS = 3          # per row: breaking addition states checked for skeleton and edge direction
SMOKE = {"ecoli70": 2, "magic-niab": 2, "arth150": 2}  # rows from the first mixed graph of each network


# --- states ----------------------------------------------------------------------------


def hop_set(g: MPDAG, sources, hops: int) -> set:
    adj = collections.defaultdict(set)
    for a, b in g.skeleton():
        adj[a].add(b)
        adj[b].add(a)
    seen = set(sources)
    frontier = set(sources)
    for _ in range(hops):
        frontier = {w for v in frontier for w in adj[v]} - seen
        seen |= frontier
    return seen


def state(g: MPDAG, k_kept: list, adds: tuple) -> MPDAG | None:
    """meek_closure(G + adds), then the kept orientations with the radius's apply_orientations."""
    cur = g
    if adds:
        for a, b in adds:
            if cur.has_edge(a, b):
                return None
        cur = MPDAG(g.nodes, directed=set(g.directed_edges) | set(adds), undirected=g.undirected_edges)
        cur = meek_closure(cur)
        if cur is None:
            return None
    if not k_kept:
        return cur
    return apply_orientations(cur, k_kept)


def breaks(st: MPDAG | None, x, y, z) -> bool | None:
    """None for an inconsistent state, else whether Z fails the GAC in it."""
    if st is None:
        return None
    return not is_gac_valid_mpdag(st, x, y, z)


# --- one row -----------------------------------------------------------------------------

_W: dict = {}


def init_worker(cfg: dict) -> None:
    _W["cfg"] = cfg
    models = Path(cfg["models_dir"])
    _W["dags"] = {}
    for net in cfg["networks"]:
        p = models / f"{net}.json"
        _W["dags"][net] = to_mpdag(parse_file(p, file_sha(p)))
    _W["graphs"] = load_learned(Path(cfg["data_dir"]), tuple(cfg["networks"]))


def run_row(st: dict) -> dict:
    cfg = _W["cfg"]
    clear_cache()
    info = _W["graphs"][st["gid"]]
    net, g = info["net"], info["g"]
    dag = _W["dags"][net]
    x, y = st["x"], st["y"]
    t0 = time.perf_counter()
    rec = {"unit": st["unit"], "gid": st["gid"], "network": net, "method": info["method"], "n": info["n"],
           "seed": info["seed"], "x": x, "y": y, "arm": ARM, "criterion": CRITERION,
           "repo_commit": cfg["repo_commit"], "z_valid_at_truth": st["z_valid_at_truth"],
           "first_failing": st.get("first_failing")}
    fails = []
    k, _ = oracle_r_knowledge(dag, g)
    g0 = apply_orientations(g, k)
    if g0 is None:
        rec["consistency_fail"] = ["g0_rebuild_failed"]
        return rec
    cs = committed_adjustment_set(g0, x, y)
    if cs.z is None or sorted(cs.z) != sorted(st["z"]):
        fails.append(f"reconstruct:z:{None if cs.z is None else sorted(cs.z)}!={st['z']}")
    if len(k) != st["n_k"]:
        fails.append(f"reconstruct:n_k:{len(k)}!={st['n_k']}")
    if cs.z is None:
        rec["consistency_fail"] = fails
        return rec
    z = cs.z
    if not is_gac_valid_mpdag(g0, x, y, z):
        fails.append("z_not_valid_in_g0")
    if int(is_valid_adjustment_set_gac_dag(dag, x, y, z)) != st["z_valid_at_truth"]:
        fails.append("z_valid_at_truth_differs")
    rec.update(z=sorted(z), z_size=len(z), n_k=len(k), k_sha256=k_hash(k), reconstruct_ok=int(not fails))

    # orientation-only radius, the ledger function
    r_claim, status, r_id, witness = breaking_depths(g, k, x, y, z)
    rec.update(r_orient=r_claim, r_orient_status=status, r_orient_id=r_id,
               r_orient_witness=[f"{a}->{b}" for a, b in witness])

    # candidate additions
    region = hop_set(g, {x, y} | set(z), HOPS)
    nodes = sorted(region)
    pairs = [(a, b) for a, b in itertools.combinations(nodes, 2) if not g.has_edge(a, b)]
    adds = [(a, b) for a, b in pairs] + [(b, a) for a, b in pairs]
    rec.update(region_size=len(region), n_candidate_pairs=len(pairs))

    # depth 1: every single move
    t1 = time.perf_counter()
    w_ret, w_add, n_incons, n_tried, equiv_done = [], [], 0, 0, 0
    reasons = collections.Counter()
    complete1 = True
    for e in k:
        kept = [f for f in k if f != e]
        b = breaks(state(g, kept, ()), x, y, z)
        n_tried += 1
        if b is None:
            n_incons += 1
        elif b:
            w_ret.append(f"{e[0]}->{e[1]}")
    for a, b_ in adds:
        if time.perf_counter() - t1 > cfg["depth1_limit_s"]:
            complete1 = False
            break
        s1 = state(g, k, ((a, b_),))
        n_tried += 1
        if s1 is None:
            n_incons += 1
            continue
        if breaks(s1, x, y, z):
            w_add.append(f"{a}->{b_}")
            if equiv_done < N_STATE_CHECKS:
                equiv_done += 1
                if not s1.is_directed_edge(a, b_):
                    fails.append(f"added_edge_not_directed_as_added:{a}->{b_}")
                if s1.skeleton() != g.skeleton() | {tuple(sorted((a, b_)))}:
                    fails.append(f"skeleton_not_g_plus_added:{a}->{b_}")
                reasons[why_invalid_gac(s1, x, y, z)] += 1
    # the depth-1 retractions are the radius's own depth-1 subsets: they must agree with r_orient
    if len(k) <= 300 and (r_claim == 1) != bool(w_ret):
        fails.append(f"depth1_retractions_disagree_with_r_orient:{r_claim}:{len(w_ret)}")
    w_pairs = sorted({"-".join(sorted(e.split("->"))) for e in w_add})
    rec.update(depth1_complete=int(complete1), depth1_seconds=round(time.perf_counter() - t1, 2),
               depth1_n_tried=n_tried, depth1_n_inconsistent=n_incons, witnesses_retract=w_ret,
               witnesses_add=w_add, witness_pairs=w_pairs, n_witness_pairs=len(w_pairs),
               n_witness_retract=len(w_ret), first_add_witness_reasons=dict(reasons))

    # joint radius
    r_joint, r_status, r_witness = None, None, []
    if w_ret or w_add:
        r_joint, r_status = 1, "exact"
        r_witness = (w_ret or w_add)[:1]
    elif not complete1:
        r_status = "censored_at_1"
    else:
        t2 = time.perf_counter()
        ret_moves = [("r", e) for e in k]
        add_moves = [("a", e) for e in adds]
        r_status = "gt_2"

        def combos():
            yield from itertools.combinations(ret_moves, 2)
            yield from itertools.product(ret_moves, add_moves)
            yield from itertools.combinations(add_moves, 2)

        n2 = 0
        for m1, m2 in combos():
            if time.perf_counter() - t2 > cfg["depth2_limit_s"]:
                r_status = "censored_at_2"
                break
            rs = {m[1] for m in (m1, m2) if m[0] == "r"}
            ad = tuple(m[1] for m in (m1, m2) if m[0] == "a")
            if len(ad) == 2 and {ad[0][0], ad[0][1]} == {ad[1][0], ad[1][1]}:
                continue
            n2 += 1
            if breaks(state(g, [f for f in k if f not in rs], ad), x, y, z):
                r_joint, r_status = 2, "exact"
                r_witness = [("retract " if m[0] == "r" else "add ") + f"{m[1][0]}->{m[1][1]}" for m in (m1, m2)]
                break
        rec.update(depth2_seconds=round(time.perf_counter() - t2, 2), depth2_n_tried=n2)
    rec.update(r_joint=r_joint, r_joint_status=r_status, r_joint_witness=r_witness)

    # coverage of the true missing edges inside An_true({x, y} u Z), from E22
    miss = [tuple(e.split("->")) for e in st.get("s3_missing_in_A") or []]
    add_set = set(w_add)
    pair_set = set(w_pairs)
    rec["true_missing_in_A"] = [f"{a}->{b}" for a, b in miss]
    rec["true_missing_in_region"] = [f"{a}->{b}" for a, b in miss if a in region and b in region]
    rec["true_missing_witness_pair"] = [f"{a}->{b}" for a, b in miss if "-".join(sorted((a, b))) in pair_set]
    rec["true_missing_witness_oriented"] = [f"{a}->{b}" for a, b in miss if f"{a}->{b}" in add_set]
    rec["covered_pair"] = int(bool(rec["true_missing_witness_pair"]))
    rec["covered_oriented"] = int(bool(rec["true_missing_witness_oriented"]))
    rec["n_reversed_in_A"] = len(st.get("s3_reversed_in_A") or [])
    # secondary: the fewest of the row's true missing edges in A, added in their true direction, that break Z
    # (subsets up to size 3, or all of them when there are more), and whether adding all of them breaks Z
    r_true, r_true_w, all_b = None, [], None
    if miss:
        all_b = breaks(state(g, k, tuple(miss)), x, y, z)
        for size in range(1, min(3, len(miss)) + 1):
            for sub in itertools.combinations(miss, size):
                if breaks(state(g, k, sub), x, y, z):
                    r_true, r_true_w = size, [f"{a}->{b}" for a, b in sub]
                    break
            if r_true is not None:
                break
    rec.update(r_true_missing=r_true, r_true_missing_witness=r_true_w,
               all_true_missing_added_breaks=None if all_b is None else int(all_b))
    rec["seconds"] = round(time.perf_counter() - t0, 2)
    rec["consistency_fail"] = fails
    return rec


# --- planning ----------------------------------------------------------------------------


def select_rows(e22_rows: Path) -> list[dict]:
    keep = ("unit", "gid", "x", "y", "z", "n_k", "z_valid_at_truth", "first_failing", "s3_missing_in_A",
            "s3_reversed_in_A", "network")
    out, seen = [], set()
    for line in e22_rows.open():
        r = json.loads(line)
        if r["unit"] in seen:
            continue
        if r.get("reconstruct_ok") != 1 or r.get("excluded") or r.get("consistency_fail") or r.get("z") is None:
            continue
        if r.get("network") not in NETWORKS:
            continue
        seen.add(r["unit"])
        out.append({f: r.get(f) for f in keep})
    return sorted(out, key=lambda r: r["unit"])


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
    ap.add_argument("--smoke", action="store_true", help="two rows (one invalid, one valid) from one graph of three networks")
    ap.add_argument("--out", default="results/rebuttal/E28")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--e22-rows", default="results/rebuttal/E22/rows.shard0of1.jsonl")
    ap.add_argument("--data-dir", default="experiments/rebuttal/_data/pest")
    ap.add_argument("--models-dir", default="results/axisa3/networks/example_models")
    ap.add_argument("--mem-gb", type=float, default=2.0, help="address-space cap of this process")
    ap.add_argument("--depth1-limit", type=float, default=DEPTH1_LIMIT_S)
    ap.add_argument("--depth2-limit", type=float, default=DEPTH2_LIMIT_S)
    args = ap.parse_args(argv)
    if not 0 <= args.shard < args.n_shards:
        ap.error("--shard must be in [0, --n-shards)")
    if args.mem_gb > 0:
        cap = int(args.mem_gb * 2**30)
        resource.setrlimit(resource.RLIMIT_AS, (cap, cap))

    def rel(p):
        p = Path(p)
        return p if p.is_absolute() else REPO / p

    e22_rows, data_dir, models_dir, out = rel(args.e22_rows), rel(args.data_dir), rel(args.models_dir), rel(args.out)
    if args.smoke and args.out == "results/rebuttal/E28":
        out = REPO / "results/rebuttal/E28_smoke"
    out.mkdir(parents=True, exist_ok=True)

    rows = select_rows(e22_rows)
    n_all = len(rows)
    n_inv = sum(r["z_valid_at_truth"] == 0 for r in rows)
    n_inv_missing = sum(r["z_valid_at_truth"] == 0 and bool(r["s3_missing_in_A"]) for r in rows)
    if args.smoke:
        picked = []
        for net, m in SMOKE.items():
            by = collections.defaultdict(list)
            for r in rows:
                if r["network"] == net:
                    by[r["gid"]].append(r)
            gid = next(gd for gd in sorted(by) if {r["z_valid_at_truth"] for r in by[gd]} == {0, 1})
            inv = [r for r in by[gid] if r["z_valid_at_truth"] == 0][: m // 2]
            val = [r for r in by[gid] if r["z_valid_at_truth"] == 1][: m - len(inv)]
            picked += inv + val
        rows = picked

    inputs = {"e22_rows": file_sha(e22_rows)}
    for name in ("graphs.jsonl", "graphs_retry.jsonl"):
        if (data_dir / name).exists():
            inputs[name.replace(".", "_")] = file_sha(data_dir / name)
    for net in NETWORKS:
        inputs[f"model_{net.replace('-', '_')}"] = file_sha(models_dir / f"{net}.json")
    cfg = {
        "experiment": EXPERIMENT, "criterion": CRITERION, "arm": ARM, "repo_commit": repo_commit(REPO),
        "smoke": bool(args.smoke), "networks": list(NETWORKS), "hops": HOPS, "max_joint_depth": MAX_JOINT_DEPTH,
        "depth1_limit_s": args.depth1_limit, "depth2_limit_s": args.depth2_limit,
        "n_rows_population": n_all, "n_invalid_population": n_inv, "n_invalid_missing_in_A": n_inv_missing,
        "data_dir": str(data_dir), "models_dir": str(models_dir), "input_sha256": inputs,
        "run_sha256": file_sha(Path(__file__)), "prereg_sha256": file_sha(HERE / "PREREG.md"),
        "vendor_bkrobust_sha256": file_sha(HERE / "vendor_bkrobust.py"),
        "vendor_e5_sha256": file_sha(HERE / "vendor_e5.py"),
        "shard": args.shard, "n_shards": args.n_shards, "mem_gb": args.mem_gb,
    }
    tag = f"shard{args.shard}of{args.n_shards}"
    (out / f"config.{tag}.json").write_text(json.dumps(cfg, indent=1))

    done = done_units(out)
    mine = [r for i, r in enumerate(rows) if i % args.n_shards == args.shard and r["unit"] not in done]
    mine.sort(key=lambda r: (r["gid"], r["unit"]))
    print(f"E28 {tag}: population {n_all} rows ({n_inv} invalid, {n_inv_missing} with a true edge missing in A); "
          f"{len(mine)} rows to run, criterion {CRITERION}, commit {cfg['repo_commit'][:12]}", flush=True)

    init_worker(cfg)
    t0 = time.perf_counter()
    fails = []
    with (out / f"rows.{tag}.jsonl").open("a") as fh:
        for i, st in enumerate(mine, 1):
            try:
                rec = run_row(st)
            except MemoryError:
                clear_cache()
                rec = {"unit": st["unit"], "gid": st["gid"], "error": "MemoryError", "consistency_fail": []}
                with (out / f"errors.{tag}.jsonl").open("a") as eh:
                    eh.write(json.dumps({"unit": st["unit"], "error": "MemoryError"}) + "\n")
                print(f"  {st['unit']}: MemoryError, not written", flush=True)
                continue
            fh.write(json.dumps(rec, default=str) + "\n")
            fh.flush()
            if rec.get("consistency_fail"):
                fails.append((rec["unit"], rec["consistency_fail"]))
            print(f"  {i}/{len(mine)} {st['unit']} inv={1 - st['z_valid_at_truth']} r_orient={rec.get('r_orient')} "
                  f"r_joint={rec.get('r_joint')}/{rec.get('r_joint_status')} pairs={rec.get('n_candidate_pairs')} "
                  f"wp={rec.get('n_witness_pairs')} cov={rec.get('covered_pair')} {rec.get('seconds')} s "
                  f"[{time.perf_counter() - t0:.0f} s]", flush=True)
    print(f"done: {len(mine)} rows in {time.perf_counter() - t0:.0f} s; {len(fails)} with a consistency failure", flush=True)
    if fails:
        (out / f"consistency_failures.{tag}.json").write_text(json.dumps(fails, indent=1))
        for u, f in fails[:20]:
            print(f"  FAIL {u}: {f}", flush=True)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
