"""E27: the amenability identity r_val == r_amen on populations not yet scored (see PREREG.md).

Population A: every loadable network not among the 33 of the CPDAG-only frame, queries by the
frame_build.py rule (census; expected empty). Population B: the frozen CPDAG-only frame x the ten
elicited knowledge states, states built as E11 builds them.

One work item is (condition, network). Output: one JSON line per unit in
<out>/units.shard-<I>-of-<K>.jsonl; a rerun skips every unit already present.
Exit codes: 0 ok; 1 an item raised; 2 consistency failure; 3 input hash mismatch.

    .venv/bin/python experiments/rebuttal/E27_amenability_identity_new_populations/run.py --smoke
    .venv/bin/python experiments/rebuttal/E27_amenability_identity_new_populations/run.py --workers 6
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import resource
import signal
import subprocess
import sys
import threading
import time
import traceback
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "src"))

from bkrobust.benchmarks import describe as _describe  # noqa: E402
from bkrobust.benchmarks.acquire import load_cached  # noqa: E402
from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag  # noqa: E402
from bkrobust.demo.graph import undirected_components  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag  # noqa: E402
from bkrobust.hybrid import breakdown_radius  # noqa: E402
from bkrobust.mpdag_criterion import is_amenable  # noqa: E402
from bkrobust.mpdag_criterion.paths import possible_descendants  # noqa: E402
from bkrobust.robustness import llm_survival as ls  # noqa: E402
from bkrobust.robustness import real_survival as rs  # noqa: E402
from bkrobust.search.exact_fast import radius_local_up_fast  # noqa: E402

EXPERIMENT = "E27"
CRITERION = "gac"
CPDAG_FRAME = ("results/frame/frame.jsonl",
               "6f05f044b13117d731aa345d4d549eebfa0586d7da11ceff953502003c7843fc")
TABLE1_FRAME = ("results/axis_robustness_real/frame.jsonl",
                "73881798319206be0a657c849431ad35c4cbddc5c4c3eecc4790d18c84981f38")
BUNDLE_SHA = "45119393225fee57ed1dd4c002130f805d6a655c2ee7648006470918c2dfc7bb"
FRAME_CAP = 20          # frame_build.py --cap default, the frozen frame's cap
MIN_COMPONENT = 2       # frame_build.MIN_COMPONENT
ENUM_MAX_K = 12         # enumeration check when G0 commits at most this many orientations
CALL_WALL_S = 600       # hard wall clock per radius call
HYBRID_LIMIT_S = rs.RADIUS_TIME_LIMIT_S  # 300 s, as in E11

SMOKE_CONDITIONS = ["D_LLM"]
SMOKE_NETWORKS = ["Acid_1996", "asia", "sachs"]
SMOKE_UNITS_PER_NETWORK = 3


class Wall(Exception):
    pass


def _alarm(*_):
    raise Wall


def timed(fn, *a, **kw) -> tuple[Any, float, str]:
    """(value, seconds, status) under the per-call wall clock; status ok|wall|memory."""
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(CALL_WALL_S)
    t = time.perf_counter()
    try:
        v = fn(*a, **kw)
        st = "ok"
    except Wall:
        v, st = None, "wall"
    except MemoryError:
        v, st = None, "memory"
    finally:
        signal.alarm(0)
    return v, round(time.perf_counter() - t, 6), st


def sha_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repo_commit() -> tuple[str, bool]:
    try:
        head = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(REPO), "status", "--porcelain",
                                "--untracked-files=no", "--", "src"],
                               capture_output=True, text=True, check=True).stdout.strip()
        return head, bool(dirty)
    except Exception:
        return "unknown", True


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


# --- per-process state ------------------------------------------------------------

_NET: dict[str, Any] = {}
_BUNDLE: dict[str, Any] = {}


def init_worker(mem_gb: float) -> None:
    if mem_gb > 0:
        cap = int(mem_gb * 1024 ** 3)
        try:
            resource.setrlimit(resource.RLIMIT_DATA, (cap, cap))
        except (ValueError, OSError):
            pass


def load_network(name: str, networks_dir: str):
    if name not in _NET:
        acq = load_cached(Path(networks_dir))
        recs = {r.name: r for r in acq.files}
        fname = rs.network_filename(name, recs)
        parsed = _describe.parse_file(recs[fname].path, recs[fname].sha256)
        dag = _describe.to_mpdag(parsed)
        _NET[name] = (dag, dag_to_cpdag(dag))
    return _NET[name]


def get_bundle(path: str) -> dict[str, Any]:
    if path not in _BUNDLE:
        _BUNDLE[path] = ls.load_bundle(Path(path))[0]
    return _BUNDLE[path]


# --- the two radii -------------------------------------------------------------------


def k_of(cpdag, g) -> frozenset:
    return frozenset(set(g.directed_edges) - set(cpdag.directed_edges))


def r_amen(cpdag, g0, x, y) -> tuple[int, bool]:
    """First author's e10b definition: UNREACHED if amenable at the CPDAG, else the upward search."""
    if is_amenable(cpdag, x, y):
        return UNREACHED, True
    f = radius_local_up_fast(cpdag, g0, lambda g: not is_amenable(g, x, y), max_depth=None)
    return f.radius, bool(f.exact)


def r_val(cpdag, g0, x, y, z) -> tuple[Any, bool, str]:
    out = breakdown_radius(cpdag, None, x, y, frozenset(z), g0=g0, criterion=CRITERION,
                           time_limit_s=HYBRID_LIMIT_S)
    return out.radius, bool(out.exact), out.method


def enum_radii(cpdag, g0, x, y, z) -> dict[str, Any]:
    """Both radii by plain enumeration of closed subsets of K_G0 (no shared search code)."""
    kg0 = sorted(k_of(cpdag, g0))
    inf = float("inf")
    best_amen, best_val = inf, inf
    n_closed = 0
    for keep in range(len(kg0), -1, -1):
        drop = len(kg0) - keep
        if drop >= best_amen and (z is None or drop >= best_val):
            break
        for s in combinations(kg0, keep):
            g = apply_orientations(cpdag, list(s))
            if g is None or k_of(cpdag, g) != frozenset(s):
                continue
            n_closed += 1
            if drop < best_amen and not is_amenable(g, x, y):
                best_amen = drop
            if z is not None and drop < best_val and not is_gac_valid_mpdag(g, x, y, frozenset(z)):
                best_val = drop
    back = lambda r: UNREACHED if r == inf else int(r)  # noqa: E731
    return {"enum_r_amen": back(best_amen), "enum_r_val": back(best_val) if z is not None else None,
            "enum_k_g0": len(kg0), "enum_n_closed": n_closed}


def as_num(r) -> float:
    """UNREACHED (-1) read as infinity, for order comparisons."""
    return float("inf") if r == UNREACHED else float(r)


# --- population A: census ---------------------------------------------------------------


def frame_rule_pairs(cpdag) -> tuple[list[tuple[str, str]], int]:
    """frame_build.py's CPDAG-only rule, unchanged: (sampled pairs, all frame pairs)."""
    comps = [c for c in undirected_components(cpdag) if len(c) >= MIN_COMPONENT]
    cands: set[str] = set()
    for comp in comps:
        cands |= set(comp)
        for node in comp:
            cands |= {nb for nb in cpdag.neighbors(node) if nb not in comp}
    frame = sorted((x, y) for x in sorted(cands) for y in sorted(possible_descendants(cpdag, x) - {x}))
    step = max(1, len(frame) // FRAME_CAP) if len(frame) > FRAME_CAP else 1
    return frame[::step][:FRAME_CAP], len(frame)


def census(networks_dir: Path, frame_networks: set[str]) -> dict[str, Any]:
    acq = load_cached(networks_dir)
    rows = []
    for rec in sorted(acq.files, key=lambda r: r.name):
        parsed = _describe.parse_file(rec.path, rec.sha256)
        if parsed.name in frame_networks:
            continue
        row = {"file": rec.name, "network": parsed.name, "sha256": rec.sha256}
        if parsed.bidirected:
            row.update(status="not_a_dag", frame_pairs=0, sampled_pairs=0)
        else:
            cpdag = dag_to_cpdag(_describe.to_mpdag(parsed))
            sampled, n_frame = frame_rule_pairs(cpdag)
            row.update(status="ok", nodes=len(cpdag.nodes),
                       undirected_edges=len(cpdag.undirected_edges),
                       frame_pairs=n_frame, sampled_pairs=len(sampled),
                       pairs=[list(p) for p in sampled],
                       knowledge_available=False)
        rows.append(row)
    return {"n_files_manifest": len(acq.files), "n_networks_outside_33": len(rows),
            "n_sampled_pairs": sum(r["sampled_pairs"] for r in rows),
            "n_scored_units": 0,
            "note": "no elicited knowledge exists for these networks; any sampled pair would be "
                    "recorded here and reported, not scored with invented knowledge",
            "networks": rows}


# --- population B: scoring ---------------------------------------------------------------


def unit_key(cond: str, u: dict[str, Any]) -> str:
    return f"B|{cond}|{u['network']}|{u['x']}|{u['y']}"


def score_item(item: dict[str, Any]) -> dict[str, Any]:
    t0 = time.perf_counter()
    cond, net, cfg = item["condition"], item["network"], item["cfg"]
    bundle = get_bundle(cfg["elicit"])
    common = {"experiment": EXPERIMENT, "criterion": CRITERION, "commit": cfg["commit"],
              "src_dirty": cfg["src_dirty"], "condition": cond, "smoke": cfg["smoke"],
              "naming": bundle[cond].get("naming")}
    rows, failures = [], []
    if not ls.has_network(bundle, cond, net):
        for u in item["units"]:
            rows.append({**common, **u, "key": unit_key(cond, u), "status": "network_absent_from_bundle"})
        return {"rows": rows, "failures": failures, "elapsed_s": time.perf_counter() - t0}
    dag, cpdag = load_network(net, cfg["networks_dir"])
    k, g0, g0_reason, shared = ls.build_state(bundle, cond, net, dag, cpdag)
    for u in item["units"]:
        x, y = u["x"], u["y"]
        row = {**common, **u, "key": unit_key(cond, u), "n_k": shared.get("n_k"),
               "k_g0": shared.get("k_g0"), "g0_status": g0_reason}
        if g0 is None:
            row["status"] = g0_reason
            rows.append(row)
            continue
        ra, ta, sa = timed(r_amen, cpdag, g0, x, y)
        row.update(r_amen=ra[0] if ra else None, r_amen_exact=ra[1] if ra else False,
                   t_amen_s=ta, r_amen_call=sa, amenable_at_g0=is_amenable(g0, x, y),
                   amenable_at_cpdag=is_amenable(cpdag, x, y))
        z, reason = rs.commit_z_star(g0, x, y)
        if z is None:
            row["status"] = reason
        else:
            row.update(status="ok", z_star=sorted(z), z_size=len(z))
            rv, tv, sv = timed(r_val, cpdag, g0, x, y, z)
            row.update(r_val=rv[0] if rv else None, r_val_exact=rv[1] if rv else False,
                       r_val_method=rv[2] if rv else None, t_val_s=tv, r_val_call=sv)
            both = row["r_val_exact"] and row["r_amen_exact"]
            row["both_exact"] = bool(both)
            row["equal"] = (row["r_val"] == row["r_amen"]) if both else None
            if both and row["r_val"] == 0:
                failures.append({"kind": "r_val_zero_on_valid_z", "key": row["key"]})
            if both and as_num(row["r_val"]) > as_num(row["r_amen"]):
                # every valid set needs amenability, so r_val <= r_amen must hold
                failures.append({"kind": "r_val_above_r_amen", "key": row["key"],
                                 "r_val": row["r_val"], "r_amen": row["r_amen"]})
        if row.get("k_g0") is not None and row["k_g0"] <= ENUM_MAX_K:
            e, te, se = timed(enum_radii, cpdag, g0, x, y, z)
            row.update(e or {}, t_enum_s=te, enum_call=se)
            if e and se == "ok":
                if row.get("r_amen_exact") and e["enum_r_amen"] != row["r_amen"]:
                    failures.append({"kind": "enum_vs_search_r_amen", "key": row["key"],
                                     "enum": e["enum_r_amen"], "search": row["r_amen"]})
                if z is not None and row.get("r_val_exact") and e["enum_r_val"] != row["r_val"]:
                    failures.append({"kind": "enum_vs_hybrid_r_val", "key": row["key"],
                                     "enum": e["enum_r_val"], "hybrid": row["r_val"]})
        rows.append(row)
    return {"rows": rows, "failures": failures, "elapsed_s": time.perf_counter() - t0}


def safe_score(item: dict[str, Any]) -> dict[str, Any]:
    try:
        return score_item(item)
    except Exception:
        return {"error": traceback.format_exc(), "item": f"{item['condition']}|{item['network']}"}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--mem-gb", type=float, default=2.0, help="data-memory cap per worker")
    ap.add_argument("--elicit", default=str(REPO / "results/elicit/knowledge.json"))
    ap.add_argument("--networks-dir", default=str(REPO / "results/axisa3/networks"))
    args = ap.parse_args(argv)
    if args.workers > 10:
        raise SystemExit("at most 10 workers")
    out = Path(args.out or REPO / ("results/rebuttal/E27_smoke" if args.smoke else "results/rebuttal/E27"))
    out.mkdir(parents=True, exist_ok=True)

    cp_path, t1_path = REPO / CPDAG_FRAME[0], REPO / TABLE1_FRAME[0]
    hashes = {"cpdag_only": sha_file(cp_path), "table1": sha_file(t1_path), "bundle": sha_file(Path(args.elicit))}
    bad = [n for n, want in (("cpdag_only", CPDAG_FRAME[1]), ("table1", TABLE1_FRAME[1]),
                             ("bundle", BUNDLE_SHA)) if hashes[n] != want]
    if bad:
        print(f"input hash mismatch: {bad}", file=sys.stderr)
        return 3

    frame = read_jsonl(cp_path)
    frame_nets = {r["network"] for r in frame}
    table1_pairs = {(r["network"], r["x"], r["y"]) for r in read_jsonl(t1_path)}
    tag = f"shard-{args.shard}-of-{args.n_shards}"

    if args.shard == 0:
        cen = census(Path(args.networks_dir), frame_nets)
        (out / "census_population_A.json").write_text(json.dumps(cen, indent=1))
        print(f"population A: {cen['n_networks_outside_33']} networks outside the 33, "
              f"{cen['n_sampled_pairs']} pairs by the CPDAG-only rule: "
              + ", ".join(f"{r['network']}={r['status']}/{r['sampled_pairs']}" for r in cen["networks"]),
              flush=True)

    bundle = get_bundle(args.elicit)
    conds = ls.panel_conditions(bundle, naming=None)
    nets = sorted(frame_nets)
    if args.smoke:
        conds = [c for c in SMOKE_CONDITIONS if c in conds]
        nets = [n for n in SMOKE_NETWORKS if n in nets]
    commit, src_dirty = repo_commit()
    cfg = {"elicit": args.elicit, "networks_dir": args.networks_dir, "commit": commit,
           "src_dirty": src_dirty, "smoke": args.smoke}

    done = {r["key"] for f in sorted(out.glob("units*.jsonl")) for r in read_jsonl(f)}
    items = []
    for net in nets:
        nu = [{"network": r["network"], "x": r["x"], "y": r["y"],
               "component_size": r["component_size"],
               "in_table1_frame": (r["network"], r["x"], r["y"]) in table1_pairs}
              for r in frame if r["network"] == net]
        if args.smoke:
            nu = nu[:SMOKE_UNITS_PER_NETWORK]
        for cond in conds:
            items.append({"condition": cond, "network": net, "units": nu})
    items = [it for i, it in enumerate(items) if i % args.n_shards == args.shard]
    pending = []
    for it in items:
        todo = [u for u in it["units"] if unit_key(it["condition"], u) not in done]
        if todo:
            pending.append({**it, "units": todo, "cfg": cfg})
    n_units = sum(len(it["units"]) for it in pending)
    config = {"experiment": EXPERIMENT, "criterion": CRITERION, "smoke": args.smoke,
              "population_B": "cpdag_only_frame x all bundle conditions", "conditions": conds,
              "n_conditions": len(conds), "n_frame_pairs": len(frame), "input_sha256": hashes,
              "hash_ok": not bad, "commit": commit, "src_dirty": src_dirty,
              "call_wall_s": CALL_WALL_S, "hybrid_limit_s": HYBRID_LIMIT_S, "enum_max_k": ENUM_MAX_K,
              "mem_gb_per_worker": args.mem_gb, "workers": args.workers,
              "n_items_pending": len(pending), "n_units_pending": n_units,
              "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    (out / f"config.{tag}.json").write_text(json.dumps(config, indent=1))
    print(f"E27 {tag}: {len(pending)} items, {n_units} units, {len(conds)} conditions, commit={commit[:10]}",
          flush=True)

    state = {"done_items": 0, "units": 0, "fail": 0, "err": 0}
    t_start = time.time()
    stop = threading.Event()

    def heartbeat() -> None:
        while not stop.wait(60.0):
            print(f"[{time.time() - t_start:8.0f}s] items {state['done_items']}/{len(pending)} "
                  f"units {state['units']}/{n_units} failures {state['fail']} errors {state['err']}", flush=True)

    threading.Thread(target=heartbeat, daemon=True).start()
    units_fh = (out / f"units.{tag}.jsonl").open("a", encoding="utf-8")
    fail_fh = (out / f"consistency_failures.{tag}.jsonl").open("a", encoding="utf-8")
    err_fh = (out / f"errors.{tag}.jsonl").open("a", encoding="utf-8")

    def consume(res: dict[str, Any]) -> None:
        state["done_items"] += 1
        if "error" in res:
            state["err"] += 1
            err_fh.write(json.dumps(res) + "\n")
            err_fh.flush()
            print(f"ERROR in {res['item']}:\n{res['error']}", flush=True)
            return
        for f in res["failures"]:
            state["fail"] += 1
            fail_fh.write(json.dumps(f, default=str) + "\n")
            fail_fh.flush()
        for r in res["rows"]:
            units_fh.write(json.dumps(r, default=str) + "\n")
            state["units"] += 1
        units_fh.flush()
        rows = res["rows"]
        ok = [r for r in rows if r.get("status") == "ok"]
        neq = sum(1 for r in ok if r.get("equal") is False)
        print(f"[{time.time() - t_start:8.0f}s] {rows[0]['condition'] if rows else '?'}|"
              f"{rows[0]['network'] if rows else '?'}: {len(rows)} units, ok={len(ok)}, unequal={neq}, "
              f"{res['elapsed_s']:.1f}s", flush=True)

    if args.workers <= 1:
        init_worker(args.mem_gb)
        for it in pending:
            consume(safe_score(it))
    else:
        with mp.get_context("spawn").Pool(args.workers, initializer=init_worker,
                                          initargs=(args.mem_gb,), maxtasksperchild=20) as pool:
            for res in pool.imap_unordered(safe_score, pending):
                consume(res)
    stop.set()
    for fh in (units_fh, fail_fh, err_fh):
        fh.close()
    print(json.dumps({"tag": tag, **state, "elapsed_s": round(time.time() - t_start, 1)}), flush=True)
    if state["fail"]:
        print(f"CONSISTENCY FAILURE: {state['fail']}", file=sys.stderr)
        return 2
    return 1 if state["err"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
