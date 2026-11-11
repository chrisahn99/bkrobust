"""E17: worked applied case, and validity of the committed set at the true graph, by radius.

Three kinds of unit, one JSON line each:
  worked  every released D_LLM query of the worked-case networks, recomputed and compared
          with the released row (radius, claim radius, |Z|, validity at truth);
  table   the matched-control frame (5 queries per network): D_LLM, and every stored
          D_RANDOM_MATCHED draw, each with Z, GAC radius and validity at truth;
  draw    each stored random draw re-drawn from its seed, to prove the states are the ones
          the matched control used.

Run from the repository root:
    .venv/bin/python experiments/rebuttal/E17_worked_case_and_validity_by_radius/run.py --smoke
Exit 0 = every consistency check held; 2 = at least one failed (rows are still written).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
for p in (ROOT / "src", ROOT / "experiments", HERE):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import numpy as np  # noqa: E402

import final_table as ft  # noqa: E402  (imported, never modified)
import matched_control as mc  # noqa: E402  (imported, never modified)
import vendor_witness as vw  # noqa: E402
from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.demo.evaluate import (  # noqa: E402
    is_valid_adjustment_set_dag,
    optimal_adjustment_set_mpdag,
)
from bkrobust.demo.graph import canon  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.gac.dag_level import is_gac_valid_dag  # noqa: E402
from bkrobust.hybrid import breakdown_radius  # noqa: E402
from bkrobust.mpdag_criterion.criterion import is_amenable  # noqa: E402

CRITERION = "gac"
RECORDS_SHA256 = "13541300e559b7aae464d35121b1f5f8a9678209dac99063ba1e31768412452e"
SMOKE_NETWORKS = ("Acid_1996", "Polzer_2012")

# Filled in the parent before the pool forks, so workers inherit the parsed graphs.
NETS: dict = {}
CFG: dict = {}
RELEASE: dict = {}


def sha256(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def git(*cmd: str) -> str:
    try:
        return subprocess.run(["git", "-C", str(ROOT), *cmd], capture_output=True,
                              text=True, timeout=30).stdout.strip()
    except Exception:
        return ""


def load_release(path: Path) -> dict:
    import csv
    out = {}
    with path.open(newline="") as fh:
        for r in csv.DictReader(fh):
            if r["arm"] == "D_LLM" and r["rep"] == "0":
                out[(r["network"], r["x"], r["y"])] = r
    return out


def build_units(args, knowledge, queries, draws, stored_rval) -> list[dict]:
    """Deterministic unit list; the order fixes the sharding."""
    units = []
    for net in args.worked_networks:
        for (n, x, y), rel in sorted(RELEASE.items()):
            if n == net and net in NETS and net in knowledge:
                units.append({"kind": "worked", "arm": "D_LLM", "network": net, "seed_idx": -1,
                              "x": x, "y": y, "k": knowledge[net], "release": rel})
    universe = sorted(set(knowledge) & set(queries) & set(NETS) - set(args.skip))
    if args.smoke:
        universe = [n for n in universe if n in SMOKE_NETWORKS]
    for net in universe:
        for x, y in queries[net]:
            units.append({"kind": "table", "arm": "D_LLM", "network": net, "seed_idx": -1,
                          "x": x, "y": y, "k": knowledge[net]})
    for d in draws:
        net = d["network"]
        if net not in universe or (args.smoke and d["seed_idx"] >= 2):
            continue
        units.append({"kind": "draw", "arm": "D_RANDOM_MATCHED", "network": net,
                      "seed_idx": d["seed_idx"], "x": "", "y": "", "k": d["k"], "stored": d,
                      "target_k": len(knowledge[net])})
        if d["exhausted"]:
            continue  # no state exists for this draw; counted from the draw units
        for x, y in queries[net]:
            units.append({"kind": "table", "arm": "D_RANDOM_MATCHED", "network": net,
                          "seed_idx": d["seed_idx"], "x": x, "y": y, "k": d["k"],
                          "stored_r_val": stored_rval.get((net, d["seed_idx"], x, y))})
    for u in units:
        u["key"] = f"{u['kind']}|{u['arm']}|{u['network']}|{u['seed_idx']}|{u['x']}|{u['y']}"
    return units


def check_draw(u: dict) -> dict:
    cpdag = NETS[u["network"]]["cpdag"]
    rng = np.random.default_rng(mc._draw_seed(CFG["seed"], u["network"], u["seed_idx"]))
    k, n_rej, n_att = mc.draw_random_k(rng, sorted(cpdag.undirected_edges), u["target_k"],
                                       cpdag, CFG["max_attempts"])
    stored = u["stored"]
    same = ((k is None) == stored["exhausted"]
            and (k is None or [list(e) for e in k] == [list(e) for e in stored["k"]])
            and n_att == stored["n_attempts"] and u["target_k"] == stored["target_k"])
    return {"redrawn_k": k and [list(e) for e in k], "n_attempts": n_att,
            "exhausted": stored["exhausted"], "checks": [] if same else ["draw_not_reproduced"]}


def compute_query(u: dict) -> dict:
    net, x, y = u["network"], u["x"], u["y"]
    dag, cpdag = NETS[net]["dag"], NETS[net]["cpdag"]
    k = [tuple(e) for e in u["k"]]
    out: dict = {"status": None, "checks": [], "notes": []}
    g0 = apply_orientations(cpdag, k)
    if g0 is None:
        out["status"] = "inconsistent"
        return out
    truth = dag.directed_edges
    closure = [(a, b) for (a, b) in g0.directed_edges if canon(a, b) in cpdag.undirected_edges]
    out["k_g0"] = len(closure)
    out["n_wrong_closure"] = sum(1 for (a, b) in closure if (b, a) in truth)
    out["n_wrong_claims"] = sum(1 for (a, b) in k if (b, a) in truth)
    out["wrong_claims"] = [[a, b] for (a, b) in k if (b, a) in truth]
    out["amenable"] = bool(is_amenable(g0, x, y))
    if not out["amenable"]:
        out["status"] = "not_amenable"
        return out
    zset = optimal_adjustment_set_mpdag(g0, x, y)
    if zset is None:
        out["status"] = "no_z"
        return out
    z = frozenset(zset)
    out["z"] = sorted(z)
    out["z_size"] = len(z)
    hy = breakdown_radius(cpdag, None, x, y, z, g0=g0, search_budget=CFG["search_budget"],
                          criterion=CRITERION, time_limit_s=CFG["ladder_time_limit"])
    out.update(radius=hy.radius, radius_method=hy.method, radius_exact=hy.exact,
               radius_witness=hy.witness)
    # A second oracle that shares no code with the criterion, where enumeration is affordable.
    out["enum_checked"] = len(cpdag.undirected_edges) <= CFG["enum_max_undirected"]
    if out["enum_checked"]:
        he = breakdown_radius(cpdag, None, x, y, z, g0=g0, search_budget=CFG["search_budget"],
                              use_criterion=False, criterion=CRITERION,
                              time_limit_s=CFG["ladder_time_limit"])
        out["radius_enum"] = he.radius
        if hy.exact and he.exact and he.radius != hy.radius:
            out["checks"].append("criterion_vs_enumeration")
    depth, subset, n_sets, cost = vw.witness(cpdag, k, x, y, z)
    out.update(claim_radius=depth, claim_witness=vw.render(subset) if subset else None,
               claim_n_sets=n_sets, claim_cost=cost or None)
    if subset:
        out["claim_witness_wrong_at_truth"] = [((b, a) in truth) for (a, b) in subset]
    out["valid_gac_truth"] = bool(is_gac_valid_dag(dag, x, y, z))
    out["valid_backdoor_truth"] = bool(is_valid_adjustment_set_dag(dag, x, y, z))
    out["status"] = "ok"
    # Theorem checks: the true DAG lies in the CPDAG's class, and in [G0] when no closure
    # orientation is wrong, so a set valid at the top (or at G0) cannot fail at truth.
    if hy.exact and hy.radius == UNREACHED and not out["valid_gac_truth"]:
        out["checks"].append("unreachable_but_invalid_at_truth")
    if hy.exact and hy.radius == UNREACHED and depth is not None:
        out["checks"].append("unreachable_but_claim_retraction_breaks")
    if hy.radius != 0 and not out["valid_gac_truth"] and out["n_wrong_closure"] == 0:
        out["checks"].append("invalid_at_truth_without_wrong_closure")
    return out


def compare_release(row: dict, rel: dict) -> None:
    """Fail on a radius, |Z|, claim radius or validity that differs from the released row."""
    rs = rel["status"]
    row["release"] = {c: rel[c] for c in ("status", "z_size", "r_hop", "hop_exact", "hop_method",
                                           "r_claim", "r_claim_status", "witness",
                                           "z_valid_at_truth", "z_backdoor_at_truth",
                                           "n_wrong_claims", "k_sha256")}
    if CFG["k_sha16"].get(row["network"]) != rel["k_sha256"]:
        row["checks"].append("release_k_differs")
    if rs == "y_not_possible_descendant" or row["status"] == "timeout":
        return
    if (rs == "ok") != (row["status"] == "ok"):
        row["checks"].append(f"release_status_differs:{rs}/{row['status']}")
        return
    if rs != "ok":
        return
    if int(rel["z_size"]) != row["z_size"]:
        row["checks"].append("release_z_size_differs")
    if rel["hop_exact"] == "1" and row["radius_exact"] and int(rel["r_hop"]) != row["radius"]:
        row["checks"].append("release_radius_differs")
    st = rel["r_claim_status"]
    if st == "exact" and int(rel["r_claim"]) != row["claim_radius"]:
        row["checks"].append("release_claim_radius_differs")
    if st in ("unreached", "gt_3", "no_retractable_edges") and row["claim_radius"] is not None:
        row["checks"].append("release_claim_radius_differs")
    if rel["witness"] and row["claim_witness"] != rel["witness"]:
        row["notes"].append("witness_differs_from_release")  # several sets can tie at one depth
    if int(rel["z_valid_at_truth"]) != int(row["valid_gac_truth"]):
        row["checks"].append("release_validity_differs")


def run_unit(u: dict) -> dict:
    row = {"key": u["key"], "kind": u["kind"], "arm": u["arm"], "network": u["network"],
           "seed_idx": u["seed_idx"], "x": u["x"], "y": u["y"], "k": [list(e) for e in u["k"] or []],
           "len_k": len(u["k"] or []), "criterion": CRITERION, "commit": CFG["commit"],
           "smoke": CFG["smoke"]}
    t0 = time.perf_counter()
    try:
        fn = check_draw if u["kind"] == "draw" else compute_query
        row.update(ft.run_with_timeout(CFG["time_limit"], fn, u))
    except ft.InstanceTimeout as exc:
        row.update(status="timeout", error=str(exc), checks=[], notes=[])
    except Exception as exc:
        row.update(status="error", error=f"{type(exc).__name__}: {exc}", checks=[], notes=[])
    if u["kind"] == "draw":
        row.setdefault("status", "ok")
    if "stored_r_val" in u:
        row["stored_r_val_backdoor_era"] = u["stored_r_val"]
    if u["kind"] == "worked":
        compare_release(row, u["release"])
    row["seconds"] = round(time.perf_counter() - t0, 4)
    return row


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smoke", action="store_true", help="two small networks, two random seeds")
    ap.add_argument("--out", default="results/rebuttal/E17")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--release-records", default="experiments/rebuttal/_data/per_query_records.csv")
    ap.add_argument("--knowledge", default="results/elicit/knowledge.json")
    ap.add_argument("--draws", default="results/matched_control/arms/D_RANDOM_MATCHED/draws.jsonl")
    ap.add_argument("--stored-instances",
                    default="results/matched_control/arms/D_RANDOM_MATCHED/instances.jsonl")
    ap.add_argument("--worked-networks", default="Acid_1996,Polzer_2012")
    ap.add_argument("--queries-per-network", type=int, default=mc.DEFAULT_QPN)
    ap.add_argument("--seed", type=int, default=mc.DEFAULT_SEED)
    ap.add_argument("--max-attempts", type=int, default=200)
    ap.add_argument("--skip", default="pathfinder")
    ap.add_argument("--time-limit", type=float, default=120.0, help="seconds per unit")
    ap.add_argument("--ladder-time-limit", type=float, default=60.0, help="seconds per SAT rung")
    ap.add_argument("--search-budget", type=int, default=3)
    ap.add_argument("--enum-max-undirected", type=int, default=8)
    ap.add_argument("--claim-depth", type=int, default=3)
    ap.add_argument("--claim-budget", type=int, default=400)
    args = ap.parse_args(argv)
    args.worked_networks = [n for n in args.worked_networks.split(",") if n]
    args.skip = [n for n in args.skip.split(",") if n]
    if not 0 <= args.shard < args.n_shards:
        ap.error("need 0 <= --shard < --n-shards")
    rp = lambda s: Path(s) if Path(s).is_absolute() else ROOT / s  # noqa: E731

    # The released search went to depth 3; the vendored default is 2.
    vw.MAX_DEPTH, vw.SUBSET_BUDGET = args.claim_depth, args.claim_budget

    records, kpath, dpath = rp(args.release_records), rp(args.knowledge), rp(args.draws)
    for p in (records, kpath, dpath):
        if not p.exists():
            print(f"missing input {p}", file=sys.stderr)
            return 1
    RELEASE.update(load_release(records))
    kn_full = json.loads(kpath.read_text())["D_LLM"]["networks"]
    knowledge = {n: [tuple(e) for e in v["k"]] for n, v in kn_full.items()}
    queries = ft.load_queries(args.queries_per_network)
    draws = [json.loads(line) for line in dpath.open() if line.strip()]
    stored_rval = {}
    sp = rp(args.stored_instances)
    if sp.exists():
        for line in sp.open():
            r = json.loads(line)
            stored_rval[(r["network"], r["seed_idx"], r["x"], r["y"])] = r["r_val"]

    needed = set(knowledge) - set(args.skip)
    if args.smoke:
        needed &= set(SMOKE_NETWORKS)
    t0 = time.perf_counter()
    NETS.update(ft.load_networks(names=needed, skip=set(args.skip), max_nodes=0))
    print(f"[E17] parsed {len(NETS)} networks in {time.perf_counter() - t0:.1f}s", flush=True)

    CFG.update(seed=args.seed, max_attempts=args.max_attempts, time_limit=args.time_limit,
               ladder_time_limit=args.ladder_time_limit, search_budget=args.search_budget,
               enum_max_undirected=args.enum_max_undirected, smoke=args.smoke,
               commit=git("rev-parse", "HEAD") or "unknown",
               k_sha16={n: v.get("k_sha256", "")[:16] for n, v in kn_full.items()})

    units = build_units(args, knowledge, queries, draws, stored_rval)
    mine = [u for i, u in enumerate(units) if i % args.n_shards == args.shard]
    out = rp(args.out)
    out.mkdir(parents=True, exist_ok=True)
    tag = f"shard{args.shard}of{args.n_shards}"
    rows_path = out / f"rows.{tag}.jsonl"
    config = {
        "experiment": "E17", "criterion": CRITERION, "commit": CFG["commit"],
        "src_dirty": bool(git("status", "--porcelain", "--", "src")), "smoke": args.smoke,
        "seed": args.seed, "queries_per_network": args.queries_per_network,
        "max_attempts": args.max_attempts, "time_limit": args.time_limit,
        "ladder_time_limit": args.ladder_time_limit, "search_budget": args.search_budget,
        "enum_max_undirected": args.enum_max_undirected, "claim_depth": args.claim_depth,
        "claim_budget": args.claim_budget, "worked_networks": args.worked_networks,
        "skip": args.skip, "records_sha256": sha256(records),
        "records_sha256_expected": RECORDS_SHA256, "knowledge_sha256": sha256(kpath),
        "draws_sha256": sha256(dpath), "frame_sha256": sha256(ft.FRAME_PATH),
        "n_units_total": len(units), "n_units_shard": len(mine), "shard": args.shard,
        "n_shards": args.n_shards, "workers": args.workers, "python": sys.version.split()[0],
        "networks_parsed": sorted(NETS),
    }
    (out / f"config.{tag}.json").write_text(json.dumps(config, indent=2))

    done = set()
    if rows_path.exists():
        for line in rows_path.open():
            try:
                done.add(json.loads(line)["key"])
            except (json.JSONDecodeError, KeyError):
                pass  # a line cut by a crash; its unit is simply run again
    todo = [u for u in mine if u["key"] not in done]
    print(f"[E17] {tag}: {len(mine)} units, {len(done)} already written, {len(todo)} to run, "
          f"workers={args.workers}", flush=True)

    t_start = last = time.perf_counter()
    n_done = 0
    with rows_path.open("a") as fh:
        def write(row: dict) -> None:
            nonlocal n_done
            fh.write(json.dumps(row) + "\n")
            fh.flush()
            n_done += 1
            if row.get("checks"):
                print(f"[E17] CHECK FAILED {row['key']}: {row['checks']}", flush=True)

        def progress() -> None:
            nonlocal last
            if time.perf_counter() - last >= 60 or n_done == len(todo):
                last = time.perf_counter()
                print(f"[E17] {n_done}/{len(todo)} done, {last - t_start:.0f}s", flush=True)

        if args.workers <= 1:
            for u in todo:
                write(run_unit(u))
                progress()
        else:
            # fork: workers inherit the parsed graphs instead of re-parsing a thousand-node net each.
            with mp.get_context("fork").Pool(args.workers) as pool:
                it = pool.imap_unordered(run_unit, todo, chunksize=1)
                while n_done < len(todo):
                    try:
                        write(it.next(timeout=30))
                    except mp.TimeoutError:
                        pass
                    progress()

    rows = [json.loads(line) for line in rows_path.open() if line.strip()]
    failed = [r for r in rows if r.get("checks")]
    status = {}
    for r in rows:
        status[r.get("status")] = status.get(r.get("status"), 0) + 1
    print(f"[E17] {tag} finished: {len(rows)} rows, status {status}, "
          f"{len(failed)} with failed checks", flush=True)
    if config["records_sha256"] != RECORDS_SHA256:
        print("[E17] WARNING: per_query_records.csv is not the released file", flush=True)
    return 2 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
