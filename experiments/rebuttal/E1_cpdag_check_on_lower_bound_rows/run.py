"""E1: the validity check at the CPDAG on the 372 released lower-bound rows.

For each released row with ``hop_method = local_up_fast_budget``, rebuild its
state exactly as ``experiments/ledger_sweep.py`` (second author's working tree)
built it, verify the rebuild against the row's hashes and counts, then, under the
GAC and under the back-door criterion:

1. the top-state check, ``Z`` valid at the CPDAG (Algorithm 1, line 2), and the
   submitted ``breakdown_radius`` wherever that check answers;
2. the search with the released settings (reproduces ``hop_lower_bound``);
3. on rows the check does not settle, the search with a larger budget.

One JSON line per released row. See PREREG.md for the registered statistic.

    python experiments/rebuttal/E1_cpdag_check_on_lower_bound_rows/run.py --smoke
    python experiments/rebuttal/E1_cpdag_check_on_lower_bound_rows/run.py --workers 8
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import multiprocessing as mp
import subprocess
import sys
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
from bkrobust.benchmarks.measure import component_of  # noqa: E402
from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.gac.mpdag_level import clear_cache, is_gac_valid_mpdag  # noqa: E402
from bkrobust.hybrid import breakdown_radius  # noqa: E402
from bkrobust.mpdag_criterion import is_valid_mpdag  # noqa: E402
from bkrobust.search.exact import retractable_edges  # noqa: E402
from search import bfs_radius  # noqa: E402
from vendor_optimal import committed_adjustment_set  # noqa: E402

DATA = ROOT / "experiments" / "rebuttal" / "_data"
CRITERIA = ("gac", "backdoor")
PRIMARY = "gac"
# Released sweep settings (experiments/ledger_sweep.py): exhaustive depth when the
# retractable set is small, else the hybrid's default depth 3; 4000 closures.
EXHAUSTIVE_MAX_RETRACTABLE = 12
RELEASE_DEPTH = 3
RELEASE_SEED = 20260913
# Released cache limits keep these at 4096 entries, but the criterion's closure
# cache is unbounded; a long search on a 400-node graph would grow it without end.
CACHE_CLEAR_EVERY = 2000
SMOKE_NETWORK = "magic-irri"  # the smallest of the seven networks (64 nodes)
SMOKE_ARMS = ("D_LLM", "D_RAND", "A_TRUE", "A_ORACLE")  # one per knowledge-rebuild path


def unit_key(r: dict) -> str:
    return "|".join((r["network"], r["x"], r["y"], r["arm"], r["rep"]))


def sha16(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()[:16]


def file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repo_commit() -> tuple[str, bool | None]:
    try:
        head = subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True, check=True
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "-C", str(ROOT), "status", "--porcelain", "--", "src"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        return head, bool(dirty)
    except (OSError, subprocess.CalledProcessError):
        return "unknown", None


# ------------------------------------------------------------------ state rebuild

_NETS: dict[str, tuple] = {}


def load_network(models: Path, name: str):
    """DAG, CPDAG and node count of a network, cached per worker process."""
    if name not in _NETS:
        # Same file rule as ledger_sweep.load: first file, sorted, whose name starts with it.
        path = next(p for p in sorted(models.iterdir()) if p.name.startswith(name))
        parsed = parse_file(path, file_sha(path))
        dag = to_mpdag(parsed)
        _NETS.clear()  # one network at a time keeps memory flat across 8 workers
        _NETS[name] = (dag, dag_to_cpdag(dag), path.name)
    return _NETS[name]


def chance_level(cpdag, base: list, rng: np.random.Generator, retries: int = 20):
    """Copied from ledger_sweep.chance_level: same edges, uniform Meek-consistent orientations."""
    for _ in range(retries):
        drawn = [(a, b) if rng.integers(2) else (b, a) for (a, b) in base]
        if apply_orientations(cpdag, drawn) is not None:
            return drawn
    return None


def rebuild_knowledge(row: dict, dag, cpdag, know: dict, seed: int) -> list | None:
    """The arm's knowledge on this network, as ledger_sweep.main builds it."""
    net, arm = row["network"], row["arm"]
    base_llm = [tuple(e) for e in (know["D_LLM"]["networks"].get(net) or {}).get("k", [])]
    if arm == "D_RAND":
        net_seed = int(hashlib.sha256(net.encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng([seed, net_seed, int(row["rep"])])
        return chance_level(cpdag, base_llm, rng)
    if arm == "A_TRUE":
        return [(a, b) if dag.is_directed_edge(a, b) else (b, a) for (a, b) in base_llm]
    if arm == "A_ORACLE":
        return sorted(knowledge_to_recover(dag, cpdag))
    if arm == "D_DEGEN":
        return []
    entry = know[arm]["networks"].get(net) or {}
    return [tuple(e) for e in entry.get("k", [])]


def make_predicate(criterion: str, x: str, y: str, z: frozenset):
    """``valid(g)`` under one criterion, clearing the memo tables now and then."""
    calls = [0]
    base = is_gac_valid_mpdag if criterion == "gac" else is_valid_mpdag

    def valid(g) -> bool:
        calls[0] += 1
        if calls[0] % CACHE_CLEAR_EVERY == 0:
            clear_cache()  # changes cost only; both modules' caches are pure functions of g
        return base(g, x, y, z)

    return valid


# ------------------------------------------------------------------ one unit


def run_unit(task: tuple) -> dict:
    row, cfg = task
    t0 = time.perf_counter()
    rec = {
        "key": unit_key(row),
        **{f: row[f] for f in ("network", "x", "y", "arm", "rep")},
        "criterion": PRIMARY,
        "criteria": list(CRITERIA),
        "commit": cfg["commit"],
        "src_dirty": cfg["src_dirty"],
        "smoke": cfg["smoke"],
        "released": {
            f: row[f]
            for f in (
                "hop_method",
                "hop_lower_bound",
                "r_hop",
                "r_claim",
                "r_claim_status",
                "n_k",
                "z_size",
                "component_size",
                "severity",
            )
        },
    }
    try:
        rec.update(compute(row, cfg))
    except Exception:  # an unrun unit is reported, not swallowed, and retried on resume
        rec["status"] = "error"
        rec["error"] = traceback.format_exc()[-2000:]
    clear_cache()
    rec["seconds"] = round(time.perf_counter() - t0, 3)
    return rec


def compute(row: dict, cfg: dict) -> dict:
    know = cfg["know"]
    dag, cpdag, model_file = load_network(Path(cfg["models"]), row["network"])
    x, y = row["x"], row["y"]
    res: dict = {"model_file": model_file, "n_nodes": len(cpdag.nodes)}
    mism: list[str] = []

    def expect(name, got, want):
        res.setdefault("rebuild", {})[name] = {"got": got, "released": want}
        if str(got) != str(want):
            mism.append(name)

    expect("chat_source_sha256", sha16(cpdag.edge_string()), row["chat_source_sha256"])
    expect("generating_dag_sha256", sha16(dag.edge_string()), row["generating_dag_sha256"])
    k = rebuild_knowledge(row, dag, cpdag, know, cfg["seed"])
    if k is None:
        mism.append("knowledge_not_rebuilt")
        return {**res, "status": "rebuild_mismatch", "mismatch": mism}
    expect("k_sha256", sha16(json.dumps(sorted(k))), row["k_sha256"])
    expect("n_k", len(k), row["n_k"])
    g0 = apply_orientations(cpdag, k)
    if g0 is None:
        mism.append("g0_inconsistent")
        return {**res, "status": "rebuild_mismatch", "mismatch": mism}
    expect("g0_undirected", len(g0.undirected_edges), row["g0_undirected"])
    k_g0 = sum(
        1 for (a, b) in g0.directed_edges if tuple(sorted((a, b))) in cpdag.undirected_edges
    )
    expect("k_g0", k_g0, row["k_g0"])
    cs = committed_adjustment_set(g0, x, y)
    expect("committed_verdict", cs.verdict, row["committed_verdict"])
    if cs.z is None:
        return {**res, "status": "rebuild_mismatch", "mismatch": mism + ["no_committed_set"]}
    z = cs.z
    expect("z_size", len(z), row["z_size"])
    comp = component_of(cpdag, x)
    expect("component_size", len(comp) if comp else 0, row["component_size"])
    if mism:
        return {**res, "status": "rebuild_mismatch", "mismatch": mism}

    n_ret = len(retractable_edges(cpdag, g0))
    release_depth = n_ret if n_ret <= EXHAUSTIVE_MAX_RETRACTABLE else RELEASE_DEPTH
    res.update({"status": "ok", "z": sorted(z), "n_retractable": n_ret})
    released_lb = int(row["hop_lower_bound"])
    released_failure = row["r_claim_status"] == "exact" or row["r_hop"] not in ("", "-1", str(UNREACHED))
    hard: list[str] = []

    for crit in CRITERIA:
        valid = make_predicate(crit, x, y, z)
        r: dict = {"criterion": crit, "valid_at_g0": valid(g0)}
        t = time.perf_counter()
        r["valid_at_cpdag"] = valid(cpdag)
        r["top_seconds"] = round(time.perf_counter() - t, 6)
        if r["valid_at_cpdag"]:
            # Only here: past the top-state check the submitted function would run the
            # E1 ladder, whose model build does not return on graphs of this size.
            hb = breakdown_radius(cpdag, None, x, y, z, g0=g0, criterion=crit)
            r["hybrid"] = {"method": hb.method, "radius": hb.radius, "exact": hb.exact,
                           "oracle": hb.oracle, "seconds": round(hb.total_seconds, 6)}
            if hb.method != "top_state" or hb.radius != UNREACHED:
                hard.append(f"{crit}: breakdown_radius gave {hb.method}/{hb.radius} on a CPDAG-valid set")
        else:
            r["hybrid"] = None

        def fails(g, valid=valid):
            return not valid(g)

        rel = bfs_radius(cpdag, g0, fails, max_depth=release_depth,
                         max_closures=cfg["release_closures"], time_limit_s=cfg["release_time_limit"])
        r["release_search"] = rel
        r["reproduces_released_lb"] = rel["lb_release_formula"] == released_lb
        r["released_lb_overstated"] = rel["lb"] is not None and rel["lb_release_formula"] == released_lb \
            and rel["lb"] < released_lb
        large = None
        if not r["valid_at_cpdag"] or cfg["search_all"]:
            large = bfs_radius(cpdag, g0, fails, max_depth=None,
                               max_closures=cfg["closure_budget"], time_limit_s=cfg["search_time_limit"])
        r["large_search"] = large

        # Outcome on the unsettled rows, against the released bound.
        if r["valid_at_cpdag"]:
            r["outcome"] = "settled_infinite"
        elif large["radius"] != UNREACHED:
            r["outcome"] = "radius_found"
        elif large["exact"]:
            r["outcome"] = "exhausted_without_failure"
        elif large["lb"] > released_lb:
            r["outcome"] = "bound_raised"
        else:
            r["outcome"] = "unchanged"
        r["best_lb"] = None if large is None or large["lb"] is None else max(large["lb"], rel["lb"] or 0)

        # Hard consistency checks: each is a contradiction, never a matter of budget.
        if crit == "gac" and not r["valid_at_g0"]:
            hard.append("gac: committed set invalid at G0, though it is read off G0")
        found = [s for s in (rel, large) if s is not None and s["radius"] != UNREACHED]
        if r["valid_at_cpdag"] and found:
            hard.append(f"{crit}: valid at the CPDAG yet a search found a failure at {found[0]['radius']}")
        if r["valid_at_cpdag"] and crit == PRIMARY and released_failure:
            hard.append("gac: valid at the CPDAG yet the released row records a failure")
        if not r["valid_at_cpdag"] and any(s["stop"] == "exhausted" for s in (rel, large) if s):
            hard.append(f"{crit}: fails at the CPDAG yet an exhaustive search found no failure")
        if large is not None and large["radius"] != UNREACHED and rel["lb"] is not None \
                and large["radius"] < rel["lb"]:
            hard.append(f"{crit}: larger search found {large['radius']} below its own certified {rel['lb']}")
        res[crit] = r

    g, b = res["gac"], res["backdoor"]
    res["criteria_disagree_at_cpdag"] = g["valid_at_cpdag"] != b["valid_at_cpdag"]
    # A back-door-valid set is GAC-valid (the back-door forbidden set contains the GAC one).
    for where in ("valid_at_cpdag", "valid_at_g0"):
        if b[where] and not g[where]:
            hard.append(f"back-door valid but GAC invalid ({where})")
    # A GAC failure is a back-door failure, so r_backdoor <= r_gac when both are found.
    rg = g["large_search"]["radius"] if g["large_search"] else UNREACHED
    rb = b["large_search"]["radius"] if b["large_search"] else UNREACHED
    if rg != UNREACHED and rb != UNREACHED and rg < rb:
        hard.append(f"GAC radius {rg} below back-door radius {rb}")
    res["hard_failures"] = hard
    return res


# ------------------------------------------------------------------ driver


def select_rows(records: Path, smoke: bool) -> list[dict]:
    rows = [r for r in csv.DictReader(records.open()) if r["hop_method"] == "local_up_fast_budget"]
    rows.sort(key=unit_key)
    if smoke:
        picked = []
        for arm in SMOKE_ARMS:
            hit = next((r for r in rows if r["network"] == SMOKE_NETWORK and r["arm"] == arm), None)
            if hit is not None:
                picked.append(hit)
        return picked
    return rows


def done_keys(out: Path) -> set[str]:
    keys = set()
    for f in sorted(out.glob("rows*.jsonl")):
        for line in f.open():
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue  # a line cut by a kill; the unit is rerun
            if rec.get("status") != "error":
                keys.add(rec["key"])
    return keys


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--records", type=Path, default=DATA / "per_query_records.csv")
    ap.add_argument("--knowledge", type=Path, default=DATA / "knowledge_states.json")
    ap.add_argument("--models", type=Path, default=ROOT / "results/axisa3/networks/example_models")
    ap.add_argument("--out", type=Path, default=ROOT / "results/rebuttal/E1")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--smoke", action="store_true", help="4 rows of magic-irri, tiny budgets, < 60 s")
    ap.add_argument("--seed", type=int, default=RELEASE_SEED, help="D_RAND seed of the released sweep")
    ap.add_argument("--release-closures", type=int, default=4000)
    ap.add_argument("--release-time-limit", type=float, default=900.0,
                    help="safety cap only; the released search had none")
    ap.add_argument("--closure-budget", type=int, default=40000)
    ap.add_argument("--search-time-limit", type=float, default=300.0)
    ap.add_argument("--search-all", action="store_true",
                    help="also run the larger search on rows the CPDAG check settles (cross-check)")
    args = ap.parse_args(argv)
    if not 0 <= args.shard < args.n_shards:
        ap.error("--shard must be in [0, --n-shards)")
    if args.smoke:
        # Smoke budgets exercise every leg (both searches, both stop reasons) quickly.
        args.release_closures, args.closure_budget = 150, 300
        args.release_time_limit = args.search_time_limit = 2.0
        args.search_all = True
        if args.out == ROOT / "results/rebuttal/E1":
            args.out = ROOT / "results/rebuttal/E1_smoke"

    commit, dirty = repo_commit()
    rows = select_rows(args.records, args.smoke)
    n_selected = len(rows)
    rows = [r for i, r in enumerate(rows) if i % args.n_shards == args.shard]
    args.out.mkdir(parents=True, exist_ok=True)
    tag = "" if args.n_shards == 1 else f".shard{args.shard}of{args.n_shards}"
    config = {
        "experiment": "E1",
        "criterion": PRIMARY,
        "criteria": list(CRITERIA),
        "commit": commit,
        "src_dirty": dirty,
        "smoke": args.smoke,
        "population": "hop_method == local_up_fast_budget",
        "n_selected": n_selected,
        "n_shard_rows": len(rows),
        "shard": args.shard,
        "n_shards": args.n_shards,
        "seed": args.seed,
        "release_closures": args.release_closures,
        "release_depth_rule": f"n_retractable if <= {EXHAUSTIVE_MAX_RETRACTABLE} else {RELEASE_DEPTH}",
        "closure_budget": args.closure_budget,
        "search_time_limit": args.search_time_limit,
        "search_all": args.search_all,
        "workers": args.workers,
        "records_sha256": file_sha(args.records),
        "knowledge_sha256": file_sha(args.knowledge),
    }
    (args.out / f"config{tag}.json").write_text(json.dumps(config, indent=1))

    skip = done_keys(args.out)
    todo = [r for r in rows if unit_key(r) not in skip]
    print(f"[E1] {len(rows)} rows in shard {args.shard}/{args.n_shards} "
          f"({n_selected} selected), {len(rows) - len(todo)} already done, {len(todo)} to run; "
          f"commit {commit[:12]}{' (src dirty)' if dirty else ''}", flush=True)
    cfg = {
        "know": json.loads(args.knowledge.read_text()),
        "models": str(args.models),
        "seed": args.seed,
        "commit": commit,
        "src_dirty": dirty,
        "smoke": args.smoke,
        "release_closures": args.release_closures,
        "release_time_limit": args.release_time_limit,
        "closure_budget": args.closure_budget,
        "search_time_limit": args.search_time_limit,
        "search_all": args.search_all,
    }
    # Grouping by network lets a worker reuse its loaded CPDAG on consecutive units.
    todo.sort(key=lambda r: (r["network"], unit_key(r)))
    tasks = [(r, cfg) for r in todo]
    n_hard = n_err = n_settled = n_ok = 0
    t0 = last = time.perf_counter()
    with (args.out / f"rows{tag}.jsonl").open("a") as fh:
        # A pool even for one worker, so the parent can print a heartbeat while a
        # long unit runs.
        pool = mp.get_context("spawn").Pool(max(1, args.workers))
        results = pool.imap_unordered(run_unit, tasks, chunksize=1)
        done = 0
        while done < len(tasks):
            try:
                rec = results.next(timeout=60)
            except mp.TimeoutError:
                print(f"[E1] {done}/{len(tasks)} done, {time.perf_counter() - t0:.0f} s, waiting", flush=True)
                continue
            done += 1
            fh.write(json.dumps(rec, default=str) + "\n")
            fh.flush()
            n_err += rec.get("status") == "error"
            n_ok += rec.get("status") == "ok"
            n_hard += bool(rec.get("hard_failures"))
            n_settled += bool(rec.get("gac", {}).get("valid_at_cpdag"))
            if rec.get("hard_failures"):
                print(f"[E1] HARD {rec['key']}: {rec['hard_failures']}", flush=True)
            if rec.get("status") == "error":
                print(f"[E1] ERROR {rec['key']}: {rec['error'].strip().splitlines()[-1]}", flush=True)
            now = time.perf_counter()
            if now - last >= 30 or done == len(tasks) or args.smoke:
                last = now
                print(f"[E1] {done}/{len(tasks)} done, ok {n_ok}, gac-settled {n_settled}, "
                      f"hard {n_hard}, errors {n_err}, {now - t0:.0f} s", flush=True)
        pool.close()
        pool.join()
    if n_hard:
        print(f"[E1] CONSISTENCY FAILURE on {n_hard} rows", flush=True)
        return 1
    if n_err:
        print(f"[E1] {n_err} rows raised; rerun to retry them", flush=True)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
