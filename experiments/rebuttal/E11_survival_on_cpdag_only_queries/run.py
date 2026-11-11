"""E11: elicited-knowledge survival on the CPDAG-only query frame (see PREREG.md).

One work item is (frame, condition, network). Inside an item the corruption draws
are shared across its queries exactly as run_llm_survival.run_shard shares them,
and the seeds come from the committed family llm|condition|network, so a query
scored here gets the same draws it got (or would have got) in the committed panel.

Output: one JSON line per unit (frame, condition, network, x, y) in
<out>/units.shard-<I>-of-<K>.jsonl, written when its item finishes. A rerun skips
every unit already present in any units*.jsonl of <out>.

Exit codes: 0 ok; 1 an item raised; 2 consistency failure; 3 input hash mismatch.

    .venv/bin/python experiments/rebuttal/E11_survival_on_cpdag_only_queries/run.py --smoke
    .venv/bin/python experiments/rebuttal/E11_survival_on_cpdag_only_queries/run.py --workers 8
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import statistics
import subprocess
import sys
import threading
import time
import traceback
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
from bkrobust.hybrid import breakdown_radius  # noqa: E402
from bkrobust.robustness import llm_survival as ls  # noqa: E402
from bkrobust.robustness import real_survival as rs  # noqa: E402
from bkrobust.robustness.llm_analyse import ENDPOINT, _contradiction_columns  # noqa: E402
from bkrobust.robustness.real_analyse import CENSORED_STATUS, _build_unit  # noqa: E402

EXPERIMENT = "E11"
CRITERION = "gac"  # registered; "backdoor" is computed beside it
FRAMES = {
    # name: (default path relative to the repository root, frozen sha256)
    "cpdag_only": ("results/frame/frame.jsonl",
                   "6f05f044b13117d731aa345d4d549eebfa0586d7da11ceff953502003c7843fc"),
    "table1": ("results/axis_robustness_real/frame.jsonl",
               "73881798319206be0a657c849431ad35c4cbddc5c4c3eecc4790d18c84981f38"),
}
BUNDLE_SHA = "45119393225fee57ed1dd4c002130f805d6a655c2ee7648006470918c2dfc7bb"
SHARD_WALL_CAP_S = 18000.0  # same safety net as run_llm_survival; a hit censors the tail

SMOKE_CONDITIONS = ["D_LLM", "D_LLM_32B"]
SMOKE_NETWORKS = ["Acid_1996", "asia", "sachs", "pathfinder"]
SMOKE_UNITS_PER_NETWORK = 3
SMOKE_N_DRAWS = 25


def sha_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repo_commit() -> tuple[str, bool]:
    """HEAD of the repository and whether src/ differs from it (the code under test)."""
    try:
        head = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(REPO), "status", "--porcelain",
                                "--untracked-files=no", "--", "src"],
                               capture_output=True, text=True, check=True).stdout.strip()
        return head, bool(dirty)
    except Exception:
        return "unknown", True


# --- frames -------------------------------------------------------------------


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def frame_units(name: str, rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Adapt either frame to one record per distinct (network, x, y).

    The Table 1 frame repeats a pair once per coverage level and carries row ids
    and committed columns; the CPDAG-only frame has one row per pair and none of
    them, so its row id is derived from the pair itself.
    """
    out: list[dict[str, Any]] = []
    if name == "table1":
        for net in sorted({r["network"] for r in rows}):
            for r in ls.distinct_pairs(rows, net):
                out.append({
                    "network": net, "x": r["x"], "y": r["y"],
                    "frame_row_id": r["row_id"],
                    "component_size": r["component_size_recomputed"],
                    "separation": r["separation_recomputed"],
                    "separation_status": r["separation_status_recomputed"],
                    "radius_committed_programmatic_k": r["radius_recomputed"],
                    "in_committed_admissible": True,
                })
        return out
    seen: set[tuple[str, str, str]] = set()
    for r in rows:
        key = (r["network"], r["x"], r["y"])
        if key in seen:
            raise SystemExit(f"cpdag_only frame repeats {key}")
        seen.add(key)
        out.append({
            "network": r["network"], "x": r["x"], "y": r["y"],
            "frame_row_id": hashlib.sha256(f"cpdag_only|{r['network']}|{r['x']}|{r['y']}"
                                           .encode()).hexdigest()[:16],
            "component_size": r["component_size"],
            "separation": None, "separation_status": "not_computed_e11",
            "radius_committed_programmatic_k": None,
            "in_committed_admissible": r["in_committed_admissible"],
            "inclusion_probability": r["inclusion_probability"],
        })
    return out


def unit_key(frame: str, cond: str, u: dict[str, Any]) -> str:
    return f"{frame}|{cond}|{u['network']}|{u['x']}|{u['y']}"


# --- per-process caches ---------------------------------------------------------

_NET: dict[str, Any] = {}
_BUNDLE: dict[str, Any] = {}


def load_network(name: str, networks_dir: str):
    """rs.load_network with an explicit cache directory, so the run does not depend on cwd."""
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


# --- scoring --------------------------------------------------------------------


def radius_cols(cpdag, g0, x, y, z, criterion: str, time_limit_s: float) -> dict[str, Any]:
    """rs.radius_columns with the criterion passed explicitly (it uses the default)."""
    out = breakdown_radius(cpdag, None, x, y, frozenset(z), g0=g0,
                           criterion=criterion, time_limit_s=time_limit_s)
    row = {"dispatch_leg": out.method, "oracle": out.oracle, "exact": out.exact,
           "assumes": out.assumes, "r_total_seconds": round(out.total_seconds, 4)}
    if not out.exact:
        row.update(radius=None, r_status="censored_wall_cap")
    elif out.radius == UNREACHED:
        row.update(radius=UNREACHED, r_status="unreached")
    else:
        row.update(radius=out.radius, r_status="ok")
    return row


def score_item(item: dict[str, Any]) -> dict[str, Any]:
    """Score every pending unit of one (frame, condition, network)."""
    t0 = time.perf_counter()
    frame, cond, net, cfg = item["frame"], item["condition"], item["network"], item["cfg"]
    bundle = get_bundle(cfg["elicit"])
    common = {"experiment": EXPERIMENT, "criterion": CRITERION, "commit": cfg["commit"],
              "src_dirty": cfg["src_dirty"], "frame": frame, "condition": cond,
              "n_draws_per_grid_point": cfg["n_draws"], "smoke": cfg["smoke"]}
    rows: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []

    if not ls.has_network(bundle, cond, net):
        # pathfinder was never elicited: written so it stays in the denominator
        for u in item["units"]:
            rows.append({**common, **u, "key": unit_key(frame, cond, u), "arm": "llm",
                         "status": "network_absent_from_bundle", "radius": None,
                         ENDPOINT: None, "naming": bundle[cond].get("naming")})
        return {"rows": rows, "failures": failures, "elapsed_s": time.perf_counter() - t0}

    dag, cpdag = load_network(net, cfg["networks_dir"])
    k, g0, g0_reason, shared = ls.build_state(bundle, cond, net, dag, cpdag)
    family = ls.panel_family(cond, net)
    largest = max((len(c) for c in undirected_components(cpdag)), default=0)

    irows: list[dict[str, Any]] = []
    scored: list[tuple[dict[str, Any], frozenset]] = []
    for u in item["units"]:
        irow = {**u, **shared, "arm": "llm", "gate": rs.GATE_NAME,
                "shard_id": f"e11__{frame}__{cond}__{net}", "largest_component_size": largest}
        if g0 is None:
            irow["status"] = g0_reason
        else:
            z, reason = rs.commit_z_star(g0, u["x"], u["y"])
            if z is None:
                irow["status"] = reason
            else:
                irow.update(status="ok", z_star=sorted(z), z_size=len(z))
                irow.update(radius_cols(cpdag, g0, u["x"], u["y"], z, "gac", cfg["radius_time_limit_s"]))
                bd = radius_cols(cpdag, g0, u["x"], u["y"], z, "backdoor", cfg["radius_time_limit_s"])
                irow.update({f"bd_{c}": v for c, v in bd.items()})
                if irow["radius"] == 0:
                    # Z* is GAC-valid at G0 by construction, so a GAC radius of 0 is a bug
                    failures.append({"kind": "gac_radius_zero_on_valid_z", "key": unit_key(frame, cond, u)})
                scored.append((irow, z))
        irows.append(irow)

    # Differential check: the enumeration oracle shares no code with the criterion.
    n_x = 0
    for irow, z in scored:
        if n_x >= cfg["xcheck_per_item"] or len(cpdag.undirected_edges) > cfg["xcheck_max_undirected"]:
            break
        if irow.get("r_status") != "ok":
            continue
        n_x += 1
        enum = breakdown_radius(cpdag, None, irow["x"], irow["y"], z, g0=g0, use_criterion=False,
                                criterion="gac", time_limit_s=cfg["radius_time_limit_s"])
        irow["xcheck_enum"] = {"radius": enum.radius, "exact": enum.exact, "method": enum.method}
        if enum.exact and enum.radius != irow["radius"]:
            failures.append({"kind": "criterion_vs_enumeration", "key": unit_key(frame, cond, irow),
                             "criterion": irow["radius"], "enumeration": enum.radius})

    cells: list[list[dict[str, Any]]] = [[] for _ in scored]
    curves: list[list[list[int]]] = [[] for _ in scored]
    n_censored = 0
    if scored:
        cache = rs.ClosureCache(cpdag)
        grid = rs.depth_grid(len(k))
        for gi, d in enumerate(grid):
            if time.perf_counter() - t0 > cfg["wall_cap_s"]:
                for gp in grid[gi:]:
                    for c in cells:
                        c.append({"grid_point": gp, "status": CENSORED_STATUS, "S": None, "n_eval": None})
                n_censored = len(grid) - gi
                break
            n_contra = 0
            survived = [0] * len(scored)
            for rep in range(cfg["n_draws"]):
                k_cor, _seed = rs.flip_draw(k, d, rep, family=family)
                key, g = cache.graph(k_cor)
                if g is None:
                    n_contra += 1
                    continue
                for i, (irow, z) in enumerate(scored):
                    if cache.survived(key, g, irow["x"], irow["y"], z):
                        survived[i] += 1
            for i in range(len(scored)):
                cells[i].append({"grid_point": d, "status": "ok",
                                 **rs.cell_statistics(cfg["n_draws"], n_contra, survived[i])})
                curves[i].append([d, cfg["n_draws"] - n_contra, survived[i]])

    by_id = {id(irow): i for i, (irow, _z) in enumerate(scored)}
    for irow in irows:
        i = by_id.get(id(irow))
        c = cells[i] if i is not None else []
        unit = _build_unit(irow, c, ENDPOINT)
        unit.update(_contradiction_columns(irow, c))
        keep = ("r_status", "oracle", "exact", "r_total_seconds", "z_star", "z_size",
                "bd_radius", "bd_r_status", "bd_oracle", "bd_exact", "bd_dispatch_leg",
                "xcheck_enum", "model", "family", "naming", "n_asked", "assert_rate",
                "k_accuracy", "g0_status", "g0_undirected_edges", "in_committed_admissible",
                "inclusion_probability", "radius_committed_programmatic_k")
        unit.update({c2: irow.get(c2) for c2 in keep})
        unit.update(common)
        unit["key"] = unit_key(frame, cond, irow)
        unit["in_table1_frame"] = (net, irow["x"], irow["y"]) in item["table1_pairs"]
        unit["curve_d_neval_nsurvived"] = curves[i] if i is not None else []
        unit["n_grid_points_censored_item"] = n_censored
        rows.append(unit)
    return {"rows": rows, "failures": failures, "elapsed_s": time.perf_counter() - t0}


def safe_score(item: dict[str, Any]) -> dict[str, Any]:
    try:
        return score_item(item)
    except Exception:
        return {"error": traceback.format_exc(), "item": f"{item['frame']}|{item['condition']}|{item['network']}"}


# --- driver ---------------------------------------------------------------------


def signature(r: dict[str, Any]) -> tuple:
    """What must agree when one (condition, network, x, y) is scored in both frames."""
    return (r.get("status"), r.get("radius"), r.get("bd_radius"), tuple(r.get("z_star") or ()),
            r.get(ENDPOINT), r.get("AUC_frac_contra_as_fail"))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smoke", action="store_true", help="tiny deterministic subset, < 60 s")
    ap.add_argument("--out", default=None, help="default results/rebuttal/E11 (E11_smoke with --smoke)")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--frames", default="cpdag_only,table1")
    ap.add_argument("--cpdag-frame", default=str(REPO / FRAMES["cpdag_only"][0]))
    ap.add_argument("--table1-frame", default=str(REPO / FRAMES["table1"][0]))
    ap.add_argument("--elicit", default=str(REPO / "results/elicit/knowledge.json"))
    ap.add_argument("--networks-dir", default=str(REPO / "results/axisa3/networks"))
    ap.add_argument("--include-scrambled", action="store_true", help="also the two scrambled-name controls")
    ap.add_argument("--n-draws", type=int, default=rs.N_DRAWS_DEFAULT)
    ap.add_argument("--radius-time-limit-s", type=float, default=rs.RADIUS_TIME_LIMIT_S)
    ap.add_argument("--wall-cap-s", type=float, default=SHARD_WALL_CAP_S)
    ap.add_argument("--xcheck-per-item", type=int, default=1)
    ap.add_argument("--xcheck-max-undirected", type=int, default=8)
    ap.add_argument("--allow-hash-change", action="store_true")
    args = ap.parse_args(argv)
    if args.workers > 8:
        raise SystemExit("at most 8 workers")
    if not 0 <= args.shard < args.n_shards:
        raise SystemExit("need 0 <= --shard < --n-shards")
    n_draws = SMOKE_N_DRAWS if args.smoke else args.n_draws
    out = Path(args.out or REPO / ("results/rebuttal/E11_smoke" if args.smoke else "results/rebuttal/E11"))
    out.mkdir(parents=True, exist_ok=True)

    paths = {"cpdag_only": Path(args.cpdag_frame), "table1": Path(args.table1_frame)}
    hashes = {n: sha_file(p) for n, p in paths.items()}
    bundle_sha = sha_file(Path(args.elicit))
    bad = [n for n in FRAMES if hashes[n] != FRAMES[n][1]] + (["bundle"] if bundle_sha != BUNDLE_SHA else [])
    if bad and not args.allow_hash_change:
        print(f"input hash mismatch: {bad}", file=sys.stderr)
        return 3

    frames = [f for f in args.frames.split(",") if f]
    units_by_frame = {n: frame_units(n, read_jsonl(paths[n])) for n in FRAMES}
    table1_pairs = {(u["network"], u["x"], u["y"]) for u in units_by_frame["table1"]}
    overlap = table1_pairs & {(u["network"], u["x"], u["y"]) for u in units_by_frame["cpdag_only"]}
    bundle = get_bundle(args.elicit)
    conds = ls.panel_conditions(bundle, naming=None if args.include_scrambled else ls.REAL_NAMING)
    if args.smoke:
        conds = [c for c in SMOKE_CONDITIONS if c in conds]
    commit, src_dirty = repo_commit()
    cfg = {"elicit": args.elicit, "networks_dir": args.networks_dir, "n_draws": n_draws,
           "wall_cap_s": args.wall_cap_s, "radius_time_limit_s": args.radius_time_limit_s,
           "xcheck_per_item": args.xcheck_per_item, "xcheck_max_undirected": args.xcheck_max_undirected,
           "commit": commit, "src_dirty": src_dirty, "smoke": args.smoke}

    done: dict[str, dict[str, Any]] = {}
    for f in sorted(out.glob("units*.jsonl")):
        for r in read_jsonl(f):
            done[r["key"]] = r

    items: list[dict[str, Any]] = []
    for fr in frames:
        units = units_by_frame[fr]
        nets = sorted({u["network"] for u in units})
        if args.smoke:
            nets = [n for n in SMOKE_NETWORKS if n in nets]
        for net in nets:
            nu = [u for u in units if u["network"] == net]
            if args.smoke:
                # pairs in both frames first, so the cross-frame check runs
                nu = sorted(nu, key=lambda u: ((u["network"], u["x"], u["y"]) not in overlap,
                                               u["x"], u["y"]))[:SMOKE_UNITS_PER_NETWORK]
            for cond in conds:
                items.append({"frame": fr, "condition": cond, "network": net, "units": nu})
    items = [it for i, it in enumerate(items) if i % args.n_shards == args.shard]
    pending = []
    for it in items:
        todo = [u for u in it["units"] if unit_key(it["frame"], it["condition"], u) not in done]
        if todo:
            pending.append({**it, "units": todo, "cfg": cfg, "table1_pairs": table1_pairs})
    n_units = sum(len(it["units"]) for it in pending)

    tag = f"shard-{args.shard}-of-{args.n_shards}"
    config = {"experiment": EXPERIMENT, "criterion": CRITERION, "criteria": ["gac", "backdoor"],
              "smoke": args.smoke, "n_draws": n_draws, "frames": frames, "frame_sha256": hashes,
              "bundle_sha256": bundle_sha, "hash_ok": not bad, "conditions": conds,
              "n_conditions": len(conds), "naming": "all" if args.include_scrambled else "real",
              "commit": commit, "src_dirty": src_dirty, "n_items_this_shard": len(items),
              "n_items_pending": len(pending), "n_units_pending": n_units,
              "n_cpdag_only_units_in_frame": len(units_by_frame["cpdag_only"]),
              "workers": args.workers, "wall_cap_s": args.wall_cap_s,
              "radius_time_limit_s": args.radius_time_limit_s,
              "xcheck_per_item": args.xcheck_per_item, "xcheck_max_undirected": args.xcheck_max_undirected,
              "seed_family": "llm|<condition>|<network> (llm_survival.panel_family)",
              "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    (out / f"config.{tag}.json").write_text(json.dumps(config, indent=1))
    print(f"E11 {tag}: {len(pending)} items pending of {len(items)}, {n_units} units, "
          f"{len(conds)} conditions, frames={frames}, n_draws={n_draws}, commit={commit[:10]}", flush=True)

    # signatures for the cross-frame check, from rows already on disk
    sig: dict[str, tuple] = {}
    for r in done.values():
        if r.get("status") == "ok":
            sig[r["key"].split("|", 1)[1]] = signature(r)

    state = {"done_items": 0, "units": 0, "fail": 0, "err": 0}
    t_start = time.time()
    stop = threading.Event()

    def heartbeat() -> None:
        while not stop.wait(30.0):
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
            if r.get("status") == "ok":
                pk = r["key"].split("|", 1)[1]
                s = signature(r)
                if pk in sig and sig[pk] != s:
                    state["fail"] += 1
                    fail_fh.write(json.dumps({"kind": "cross_frame_mismatch", "key": r["key"],
                                              "this": s, "other": sig[pk]}, default=str) + "\n")
                    fail_fh.flush()
                sig.setdefault(pk, s)
            units_fh.write(json.dumps(r, default=str) + "\n")
            units_fh.flush()
            state["units"] += 1
        rows = res["rows"]
        meds = [r[ENDPOINT] for r in rows if r.get(ENDPOINT) is not None]
        print(f"[{time.time() - t_start:8.0f}s] done {rows[0]['key'].rsplit('|', 2)[0] if rows else '?'}: "
              f"{len(rows)} units, ok={sum(r.get('status') == 'ok' for r in rows)}, "
              f"AUC_median={statistics.median(meds) if meds else None}, {res['elapsed_s']:.1f}s", flush=True)

    if args.workers <= 1:
        for it in pending:
            consume(safe_score(it))
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            for res in pool.imap_unordered(safe_score, pending):
                consume(res)
    stop.set()
    for fh in (units_fh, fail_fh, err_fh):
        fh.close()
    print(json.dumps({"tag": tag, **state, "elapsed_s": round(time.time() - t_start, 1)}), flush=True)
    if state["fail"]:
        print(f"CONSISTENCY FAILURE: {state['fail']} (see consistency_failures.{tag}.jsonl)", file=sys.stderr)
        return 2
    return 1 if state["err"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
