"""E12: the radius and the audit at truth on the five networks above 500 nodes.

Per unit (frame row x knowledge condition) and in this order: G0 and the committed set, the
GAC check at the CPDAG, the budgeted upward search where that check fails, the submitted
``breakdown_radius`` under a wall clock and a memory cap, the claim radius; only then the true
DAG, to score. See PREREG.md.

    python experiments/rebuttal/E12_networks_above_500_nodes/run.py --smoke
    python experiments/rebuttal/E12_networks_above_500_nodes/run.py --workers 6
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import signal
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

import networkx as nx  # noqa: E402

if not hasattr(nx, "is_d_separator"):  # networkx < 3.3
    nx.is_d_separator = nx.d_separated

from bkrobust.benchmarks.describe import parse_file, to_mpdag  # noqa: E402
from bkrobust.benchmarks.measure import component_of, separation  # noqa: E402
from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag  # noqa: E402
from bkrobust.hybrid import DEFAULT_SEARCH_BUDGET, breakdown_radius  # noqa: E402
from bkrobust.mpdag_criterion import is_valid_mpdag  # noqa: E402
from bkrobust.search.exact import retractable_edges  # noqa: E402
from search import budgeted_up_search  # noqa: E402
from vendor_ledger import (  # noqa: E402
    NO_RETRACTABLE_EDGES,
    RUNGS,
    audit_row,
    breaking_depths,
    chance_level,
    claim_lower_bound,
)
from vendor_optimal import committed_adjustment_set  # noqa: E402

import numpy as np  # noqa: E402

NETWORKS = ("link", "munin", "munin2", "munin3", "munin4")
SMOKE_NETWORKS = ("alarm", "link")
DEPLOYMENT_ARMS = (
    "D_LLM", "D_LLM_INSTR", "D_LLM_SRC", "D_SCRAMBLED", "D_LLM_72B", "D_LLM_72B_INSTR",
    "D_SCRAMBLED_72B", "D_LLM_32B", "D_LLM_SRC_70B", "D_LLM_GEMMA_27B",
)  # the ledger's list minus D_LLM_PLUS_ANC, which lives in a second file not released
EXHAUSTIVE_MAX_RETRACTABLE = 12  # the ledger's depth rule


class Deadline(BaseException):
    """BaseException so that no library `except Exception` can swallow the alarm."""


def _alarm(signum, frame):
    raise Deadline()


def with_deadline(seconds: float, fn, *a, **kw):
    """Run fn under SIGALRM; Python code (the ladder's model build) is interruptible."""
    old = signal.signal(signal.SIGALRM, _alarm)
    signal.setitimer(signal.ITIMER_REAL, max(seconds, 0.01))
    try:
        return fn(*a, **kw)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old)


def sha_file(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def wsha(s: str | None) -> str | None:
    # a witness is a whole 1000-node edge string; keep its hash
    return None if s is None else hashlib.sha256(s.encode()).hexdigest()[:16]


def sha_graph(g) -> str:
    return hashlib.sha256(g.edge_string().encode()).hexdigest()[:16]


def git(*cmd: str) -> str:
    try:
        return subprocess.check_output(["git", "-C", str(ROOT), *cmd], text=True,
                                       stderr=subprocess.DEVNULL).strip()
    except Exception:
        return "unknown"


# ---------------------------------------------------------------- per-process graph cache
_GRAPHS: dict = {}
_G0: dict = {}
SET: dict = {}  # settings, set in every process


def load(net: str):
    if net not in _GRAPHS:
        models = Path(SET["models"])
        path = next(p for p in sorted(models.iterdir()) if p.name.split(".")[0] == net)
        parsed = parse_file(path, sha_file(path))
        dag = to_mpdag(parsed)
        _GRAPHS[net] = (dag, dag_to_cpdag(dag), dict(parsed.name_map))
    return _GRAPHS[net]


def init_worker(settings: dict) -> None:
    SET.update(settings)
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    if settings["mem_limit_gb"] > 0:
        import resource
        lim = int(settings["mem_limit_gb"] * 2**30)
        # RLIMIT_DATA counts heap and anonymous maps, not the virtual reservations of
        # imported libraries, so a runaway ladder build raises MemoryError instead of
        # taking the 15 GB machine down
        try:
            resource.setrlimit(resource.RLIMIT_DATA, (lim, lim))
        except (ValueError, OSError) as e:
            print(f"  memory cap not set: {e}", flush=True)


def fails_fn(x, y, z):
    if SET["criterion"] == "gac":
        return lambda g: not is_gac_valid_mpdag(g, x, y, z)
    return lambda g: not is_valid_mpdag(g, x, y, z)


# ---------------------------------------------------------------- tasks
def oracle_task(net: str) -> dict:
    dag, cpdag, _ = load(net)
    t0 = time.perf_counter()
    try:
        k = with_deadline(SET["oracle_limit"], knowledge_to_recover, dag, cpdag)
        return {"network": net, "k": sorted(k), "seconds": time.perf_counter() - t0}
    except Deadline:
        return {"network": net, "k": None, "seconds": time.perf_counter() - t0}


def hop(cpdag, g0, x, y, z, rec: dict, cons: list) -> int | None:
    """The hop radius by the registered dispatch; returns the certificate for the B6 bucket."""
    fails = fails_fn(x, y, z)
    t0 = time.perf_counter()
    top_valid = not fails(cpdag)
    rec["top_valid"] = int(top_valid)
    rec["top_seconds"] = round(time.perf_counter() - t0, 6)
    n_ret = len(retractable_edges(cpdag, g0))
    depth = n_ret if n_ret <= EXHAUSTIVE_MAX_RETRACTABLE else DEFAULT_SEARCH_BUDGET
    rec["n_retractable"], rec["search_depth_budget"] = n_ret, depth

    srch = None
    if not top_valid:
        s = budgeted_up_search(cpdag, g0, fails, max_depth=depth,
                               closure_budget=SET["closure_budget"],
                               time_limit_s=SET["search_limit"])
        srch = s
        rec["search"] = {"stop": s.stop, "radius": s.radius, "exact": int(s.exact),
                         "lower_bound": s.lower_bound, "depth_reached": s.depth_reached,
                         "witness_sha": wsha(s.witness),
                         "seconds": round(s.seconds, 4), **s.stats.as_dict()}
        if s.stop == "degenerate":
            cons.append("search: committed set fails at G0")
        if s.stop == "exhausted":
            cons.append("search exhausted the up-set without a failure although Z fails at the CPDAG")
    else:
        rec["search"] = None

    limit = SET["hybrid_limit"]
    if srch is not None and srch.exact:
        limit = max(limit, 2 * srch.seconds + 30)
    h: dict = {"limit_s": round(limit, 1)}
    t0 = time.perf_counter()
    try:
        res = with_deadline(limit, breakdown_radius, cpdag, None, x, y, z, g0=g0,
                            search_budget=depth, criterion=SET["criterion"],
                            time_limit_s=limit)
        h.update({"outcome": "completed", "radius": res.radius, "method": res.method,
                  "exact": int(res.exact), "witness_sha": wsha(res.witness),
                  "top_seconds": round(res.top_seconds, 6),
                  "search_seconds": round(res.search_seconds, 4),
                  "ladder_seconds": round(res.ladder_seconds, 4)})
    except Deadline:
        h["outcome"] = "timeout"
    except MemoryError:
        h["outcome"] = "memory"
    except Exception as e:  # recorded, never swallowed silently
        h["outcome"] = f"error:{type(e).__name__}:{str(e)[:200]}"
    h["wall_seconds"] = round(time.perf_counter() - t0, 4)
    rec["hybrid"] = h

    done = h["outcome"] == "completed" and h["exact"]
    if done:
        if top_valid != (h["method"] == "top_state"):
            cons.append(f"CPDAG check top_valid={top_valid} but hybrid method={h['method']}")
        if top_valid and h["radius"] != UNREACHED:
            cons.append("top valid but hybrid radius finite")
        if not top_valid and h["radius"] == UNREACHED:
            cons.append("Z fails at the CPDAG but hybrid returned an infinite radius")
        if h["method"] == "degenerate":
            cons.append("hybrid: committed set fails at G0")
        if srch is not None and srch.exact and srch.stop == "found":
            if srch.radius != h["radius"]:
                cons.append(f"search radius {srch.radius} != hybrid radius {h['radius']}")
            elif h["method"] == "local_up_fast" and wsha(srch.witness) != h["witness_sha"]:
                cons.append("search and hybrid found different witnesses at the same depth")

    if top_valid:
        r, exact, method, lb = UNREACHED, True, "top_state", None
    elif done:
        r, exact, method, lb = h["radius"], True, "hybrid_" + h["method"], None
    elif srch is not None and srch.exact:
        r, exact, method, lb = srch.radius, True, "e12_search_" + srch.stop, None
    else:
        r, exact, method = UNREACHED, False, "censored"
        lb = (srch.lower_bound if srch is not None and srch.lower_bound else 1)
    rec.update({"r_hop": r, "hop_exact": int(exact), "hop_method": method,
                "hop_lower_bound": lb, "censored": int(not exact)})
    return (r if r != UNREACHED else None) if exact else lb


def single_breaks(cpdag, k, x, y, z) -> int:
    """phi_1's numerator, as in ledger_sweep.instrument_row."""
    n = 0
    for i in range(len(k)):
        gi = apply_orientations(cpdag, [e for j, e in enumerate(k) if j != i])
        if gi is not None and not is_gac_valid_mpdag(gi, x, y, z):
            n += 1
    return n


def unit_task(u: dict) -> dict:
    t_unit = time.perf_counter()
    rec = {k: u[k] for k in ("key", "network", "x", "y", "arm", "rep", "true_dag_on_path",
                             "k_sha256", "n_k")}
    rec.update({"criterion": SET["criterion"], "commit": SET["commit"], "tier": "T1",
                "smoke": SET["smoke"]})
    cons: list = []
    try:
        dag, cpdag, names = load(u["network"])
        rec.update({"nodes": len(dag.nodes), "cpdag_undirected": len(cpdag.undirected_edges),
                    "chat_source_sha256": sha_graph(cpdag),
                    "generating_dag_sha256": sha_graph(dag)})
        if u.get("status"):  # decided upstream (oracle knowledge timed out)
            rec["status"] = u["status"]
            return rec
        k = [tuple(e) for e in u["k"]]
        gkey = (u["network"], u["k_sha256"])
        if gkey not in _G0:
            _G0.clear()  # one G0 at a time; a 1000-node MPDAG per condition adds up
            _G0[gkey] = apply_orientations(cpdag, k)
        g0 = _G0[gkey]
        x, y = u["x"], u["y"]
        if g0 is None:
            rec["status"] = "knowledge_inconsistent"
            return rec
        # ---- observables only, from here to the freeze
        rec["g0_undirected"] = len(g0.undirected_edges)
        rec["k_g0"] = sum(1 for (a, b) in g0.directed_edges
                          if tuple(sorted((a, b))) in cpdag.undirected_edges)
        cs = committed_adjustment_set(g0, x, y)
        rec["committed_verdict"] = cs.verdict
        if cs.z is None:
            rec["status"] = cs.verdict
            return rec
        z = cs.z
        rec["status"], rec["z_size"] = "ok", len(z)
        sep, sep_status = separation(cpdag, x, z)
        comp = component_of(cpdag, x)
        rec.update({"separation": sep, "separation_status": sep_status,
                    "component_size": len(comp) if comp else 0})
        if not retractable_edges(cpdag, g0):
            # nothing to retract: the ledger's convention, no call needed
            rec.update({"top_valid": None, "search": None, "hybrid": None, "r_hop": UNREACHED,
                        "hop_exact": 1, "hop_method": NO_RETRACTABLE_EDGES,
                        "hop_lower_bound": None, "censored": 0})
            hop_cert = None
        else:
            hop_cert = hop(cpdag, g0, x, y, z, rec, cons)

        t0 = time.perf_counter()
        try:
            r_claim, r_status, r_id, witness = with_deadline(
                SET["claim_limit"], breaking_depths, cpdag, k, x, y, z)
            left = SET["claim_limit"] - (time.perf_counter() - t0)
            n_break = with_deadline(left, single_breaks, cpdag, k, x, y, z)
            rec.update({"r_claim": r_claim, "r_claim_status": r_status, "r_id": r_id,
                        "phi_1": round(n_break / len(k), 4) if k else None,
                        "witness": " and ".join(f"{names.get(a, a)}->{names.get(b, b)}"
                                                for a, b in witness)})
            claim_cert = claim_lower_bound(r_claim, r_status)
        except Deadline:
            rec.update({"r_claim": None, "r_claim_status": "censored_time", "r_id": None,
                        "phi_1": None, "witness": ""})
            claim_cert = 1  # Z is valid at G0, nothing more is certified
        rec["claim_seconds"] = round(time.perf_counter() - t0, 4)
        cert = {"B2_free_rule": sep if sep_status == "measured" else 1, "B4_kg0": rec["k_g0"],
                "B5_k": len(k), "B6_hop": hop_cert, "B7_claim": claim_cert}
        rec["frozen_at"] = round(time.perf_counter() - t_unit, 4)
        # ---- the truth, revealed only now
        rec.update(audit_row(dag, cpdag, k, x, y, z, cert))
    except MemoryError:
        rec["status"] = "error"
        rec["error"] = "MemoryError outside the hybrid call"
    except Exception:
        rec["status"] = "error"
        rec["error"] = traceback.format_exc()[-1500:]
    finally:
        rec["consistency"] = cons
        rec["unit_seconds"] = round(time.perf_counter() - t_unit, 4)
    return rec


# ---------------------------------------------------------------- plan
def unit_key(net, x, y, arm, rep) -> str:
    return f"{net}|{x}|{y}|{arm}|{rep}"


def in_shard(key: str, shard: int, n_shards: int) -> bool:
    return int(hashlib.sha256(key.encode()).hexdigest()[:8], 16) % n_shards == shard


def build_units(args, know, frame, done: set, pool_map) -> tuple[list, dict]:
    present = [a for a in DEPLOYMENT_ARMS if a in know]
    arms = [a for a in args.arms.split(",") if a] or [*present, "D_RAND", "D_DEGEN",
                                                       "A_TRUE", "A_ORACLE"]
    nets = args.networks
    rows = {n: [r for r in frame if r["network"] == n] for n in nets}
    if args.smoke:
        rows = {n: v[:2] for n, v in rows.items()}
    # A_ORACLE knowledge is expensive on 1000 nodes: build it only where this shard needs it
    want_oracle = [n for n in nets if "A_ORACLE" in arms and any(
        in_shard(k, args.shard, args.n_shards) and k not in done
        for k in (unit_key(n, r["x"], r["y"], "A_ORACLE", 0) for r in rows[n]))]
    oracle: dict = {}
    for n in list(want_oracle):
        cache = Path(args.out) / f"oracle_knowledge.{n}.json"
        if cache.exists():
            oracle[n] = json.loads(cache.read_text())
            want_oracle.remove(n)
    for res in pool_map(oracle_task, want_oracle):
        oracle[res["network"]] = res
        (Path(args.out) / f"oracle_knowledge.{res['network']}.json").write_text(json.dumps(res))
        print(f"  A_ORACLE knowledge {res['network']}: "
              f"{'timeout' if res['k'] is None else len(res['k'])} in {res['seconds']:.0f} s",
              flush=True)

    units = []
    for n in nets:
        dag, cpdag, _ = load(n)
        plans: list = []
        base_llm: list = []
        for a in present:
            kk = [tuple(e) for e in (know[a]["networks"].get(n) or {}).get("k", [])]
            if a == "D_LLM":
                base_llm = kk
            plans.append((a, 0, kk, "via_Chat_only", None))
        if "D_RAND" in arms and base_llm:
            for rep in range(args.rand_reps):
                net_seed = int(hashlib.sha256(n.encode()).hexdigest()[:8], 16)
                rng = np.random.default_rng([args.seed, net_seed, rep])
                kk = chance_level(cpdag, base_llm, rng)
                if kk is not None:
                    plans.append(("D_RAND", rep, kk, "via_Chat_only", None))
        if "D_DEGEN" in arms:
            plans.append(("D_DEGEN", 0, [], "via_Chat_only", None))
        if "A_TRUE" in arms and base_llm:
            truthful = [(a, b) if dag.is_directed_edge(a, b) else (b, a) for (a, b) in base_llm]
            plans.append(("A_TRUE", 0, truthful, "via_K", None))
        if "A_ORACLE" in arms and n in oracle:
            kt = oracle[n]["k"]
            plans.append(("A_ORACLE", 0, [tuple(e) for e in kt] if kt else [], "via_K",
                          None if kt is not None else "oracle_knowledge_timeout"))
        for arm, rep, k, prov, status in plans:
            if arm not in arms:
                continue
            k_sha = hashlib.sha256(json.dumps(sorted(k)).encode()).hexdigest()[:16]
            for r in rows[n]:
                key = unit_key(n, r["x"], r["y"], arm, rep)
                if key in done or not in_shard(key, args.shard, args.n_shards):
                    continue
                units.append({"key": key, "network": n, "x": r["x"], "y": r["y"], "arm": arm,
                              "rep": rep, "true_dag_on_path": prov, "k": [list(e) for e in k],
                              "k_sha256": k_sha, "n_k": len(k), "status": status})
    # longest-first is unknowable; group by condition so each worker reuses its G0
    units.sort(key=lambda u: (u["network"], u["arm"], u["rep"], u["x"], u["y"]))
    return units, {"arms": arms, "rows_per_network": {n: len(v) for n, v in rows.items()},
                   "oracle_seconds": {n: round(v.get("seconds", 0), 1) for n, v in oracle.items()}}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--smoke", action="store_true", help="2 rows of alarm and link, tiny budgets")
    ap.add_argument("--out", default="", help="default results/rebuttal/E12 (smoke: .../E12/smoke)")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--knowledge", default=str(ROOT / "experiments/rebuttal/_data/knowledge_states.json"))
    ap.add_argument("--frame", default=str(ROOT / "results/frame/frame.jsonl"))
    ap.add_argument("--models", default=str(ROOT / "results/axisa3/networks/example_models"))
    ap.add_argument("--networks", default=",".join(NETWORKS))
    ap.add_argument("--arms", default="", help="comma list; empty = all released conditions")
    ap.add_argument("--criterion", default="gac", choices=["gac", "backdoor"])
    ap.add_argument("--rand-reps", type=int, default=3)
    ap.add_argument("--seed", type=int, default=20260913)
    ap.add_argument("--closure-budget", type=int, default=40000)
    ap.add_argument("--search-limit", type=float, default=120.0)
    ap.add_argument("--hybrid-limit", type=float, default=60.0)
    ap.add_argument("--claim-limit", type=float, default=120.0)
    ap.add_argument("--oracle-limit", type=float, default=900.0)
    ap.add_argument("--mem-limit-gb", type=float, default=2.0, help="per worker; 0 disables")
    args = ap.parse_args(argv)
    if args.smoke:
        args.networks = ",".join(SMOKE_NETWORKS)
        args.rand_reps = 1
        args.closure_budget, args.search_limit = 30, 5.0
        args.hybrid_limit, args.claim_limit, args.oracle_limit = 5.0, 5.0, 10.0
    args.networks = [n for n in args.networks.split(",") if n]
    args.out = args.out or str(ROOT / ("results/rebuttal/E12/smoke" if args.smoke
                                       else "results/rebuttal/E12"))
    if not 0 <= args.shard < args.n_shards:
        ap.error("need 0 <= shard < n-shards")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    suffix = "" if args.n_shards == 1 else f".shard{args.shard}of{args.n_shards}"

    knowledge_p, frame_p = Path(args.knowledge), Path(args.frame)
    settings = {
        "smoke": args.smoke, "criterion": args.criterion, "commit": git("rev-parse", "HEAD"),
        "src_dirty": bool(git("status", "--porcelain", "src")), "models": args.models,
        "closure_budget": args.closure_budget, "search_limit": args.search_limit,
        "hybrid_limit": args.hybrid_limit, "claim_limit": args.claim_limit,
        "oracle_limit": args.oracle_limit, "mem_limit_gb": args.mem_limit_gb,
    }
    SET.update(settings)
    know = json.loads(knowledge_p.read_text())
    frame = [json.loads(line) for line in frame_p.open()]

    done: set = set()
    for f in out.glob("rows*.jsonl"):
        for line in f.open():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue  # a line cut by a kill; the unit is redone
            if r.get("status") != "error":
                done.add(r["key"])

    workers = 1 if args.smoke else max(1, min(args.workers, 8))
    for n in args.networks:
        load(n)  # before the fork, so workers share the parsed graphs
    pool = None
    if workers > 1:
        ctx = mp.get_context("fork" if sys.platform.startswith("linux") else "spawn")
        pool = ctx.Pool(workers, initializer=init_worker, initargs=(settings,),
                        maxtasksperchild=50)  # bounds memory left behind by an interrupted call
    else:
        init_worker(settings)

    def pool_map(fn, items):
        return pool.map(fn, items, chunksize=1) if pool else [fn(i) for i in items]

    units, plan = build_units(args, know, frame, done, pool_map)
    config = {**settings, "networks": args.networks, "seed": args.seed,
              "rand_reps": args.rand_reps, "shard": args.shard, "n_shards": args.n_shards,
              "workers": workers, "knowledge_sha256": sha_file(knowledge_p),
              "frame_sha256": sha_file(frame_p), "units_todo": len(units),
              "units_already_done": len(done), "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
              "python": sys.version.split()[0], "networkx": nx.__version__, **plan}
    (out / f"config{suffix}.json").write_text(json.dumps(config, indent=1))
    print(f"E12: {len(units)} units to run, {len(done)} already done, {workers} workers, "
          f"criterion {args.criterion}, commit {settings['commit'][:12]}", flush=True)

    n_cons = n_done = 0
    paths: dict = {}
    t0 = last = time.perf_counter()
    with (out / f"rows{suffix}.jsonl").open("a") as fh:
        it = pool.imap_unordered(unit_task, units, chunksize=1) if pool else map(unit_task, units)
        while True:
            try:
                rec = it.next(timeout=60) if pool else next(it)
            except mp.TimeoutError:
                rec = None
            except StopIteration:
                break
            if rec is not None:
                fh.write(json.dumps(rec, default=str) + "\n")
                fh.flush()
                n_done += 1
                n_cons += bool(rec.get("consistency"))
                tag = rec.get("hop_method") or rec.get("status")
                paths[tag] = paths.get(tag, 0) + 1
                if rec.get("consistency"):
                    print(f"  CONSISTENCY {rec['key']}: {rec['consistency']}", flush=True)
            now = time.perf_counter()
            if now - last >= 60 or (rec is not None and not pool) or n_done == len(units):
                last = now
                print(f"  [{time.strftime('%H:%M:%S')}] {n_done}/{len(units)} units, "
                      f"{now - t0:.0f} s, consistency failures {n_cons}, paths {paths}",
                      flush=True)
            if n_done == len(units):
                break
    if pool:
        pool.close()
        pool.join()
    print(f"E12 finished: {n_done} units in {time.perf_counter() - t0:.0f} s; "
          f"{n_cons} consistency failures", flush=True)
    return 3 if n_cons else 0


if __name__ == "__main__":
    sys.exit(main())
