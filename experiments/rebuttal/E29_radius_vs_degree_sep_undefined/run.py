"""E29: the radius against the undirected degree of X where separation is undefined (see PREREG.md).

Population: the 1,033 real-name units of the submitted Table 1 cohort that have a GAC radius, a survival
AUC and no separation value (E16's undefined population, selected by E16's rule). For each unit the state
is rebuilt as the GAC recompute rebuilt it (`real_survival.load_network`, `llm_survival.build_state`,
committed set read from `results/radius_gac_recompute/`, checked against `real_survival.commit_z_star`),
and the run records: the undirected degree of X and the size of X's chain component in the CPDAG; the
analyst's orientations at X (K_G0 edges incident to X); the radius by an upward search that returns its
witness state, checked against the submitted `radius_local_up_fast` and the stored GAC radius; every
failing upper cover of G0 (the distance-1 witnesses) with the orientations it retracts; and the failure
mechanism of each witness (amenability lost, forbidden node in Z, or a non-causal path left open).

One work item is (condition, network). Output: one JSON line per unit in
<out>/units.shard-<I>-of-<K>.jsonl; a rerun skips every unit already written.
Exit codes: 0 ok; 1 an item raised; 2 consistency failure; 3 input hash mismatch.

    PYTHONPATH=src .venv/bin/python experiments/rebuttal/E29_radius_vs_degree_sep_undefined/run.py --smoke
    PYTHONPATH=src .venv/bin/python experiments/rebuttal/E29_radius_vs_degree_sep_undefined/run.py --workers 2
"""

from __future__ import annotations

import argparse
import csv
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
from collections import deque
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "src"))

from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.demo.graph import undirected_components  # noqa: E402
from bkrobust.gac.mpdag_level import gac_forbidden_set, is_gac_valid_mpdag  # noqa: E402
from bkrobust.mpdag_criterion import is_amenable  # noqa: E402
from bkrobust.robustness import llm_survival as ls  # noqa: E402
from bkrobust.robustness import real_survival as rs  # noqa: E402
from bkrobust.search.exact_fast import local_up_covers_fast, radius_local_up_fast  # noqa: E402

EXPERIMENT = "E29"
CRITERION = "gac"
ENDPOINT = "AUC_frac"
INPUTS = {
    "units_csv": ("results/axis_robustness_llm/analysis_units.csv",
                  "58de266a17e2dfd4acddfba3d76cc3b1afd5d6981809d667dbbc1caa18307658"),
    "gac_csv": ("results/radius_gac_recompute/axis_robustness_llm.csv",
                "5b3e2191736d5c8bf3b7c31e8ecc6c718781594f6432c394d9a2160f48e373ef"),
    "gac_large_csv": ("results/radius_gac_recompute/axis_robustness_llm_large.csv",
                      "817addcf08e1a0160f4d1a960e503fea7c872a81204689c77a6e6b65e6f6ebd7"),
    "bundle": ("results/elicit/knowledge.json",
               "45119393225fee57ed1dd4c002130f805d6a655c2ee7648006470918c2dfc7bb"),
}
EXPECTED_N = 1033
CALL_WALL_S = 600
SMOKE_ITEMS = 4


class Wall(Exception):
    pass


def _alarm(*_):
    raise Wall


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


def _num(v: str | None) -> float | None:
    if v is None or v.strip() == "":
        return None
    try:
        return float(v)
    except ValueError:
        return None


# --- population: E16's rule, reimplemented on the same committed files ------------------


def select_population() -> tuple[list[dict[str, Any]], list[str]]:
    """Real-name submitted-cohort units with an exact GAC radius equal to the stored one, an AUC and
    no separation value (E16 run.py: build_cohorts, scorable, separation None)."""
    gac: dict[tuple[str, str], dict[str, str]] = {}
    for key in ("gac_csv", "gac_large_csv"):
        with (REPO / INPUTS[key][0]).open(encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                if r["gac_radius"] == "":
                    continue
                gac[(r["shard_id"], r["frame_row_id"])] = r
    units, fails = [], []
    with (REPO / INPUTS["units_csv"][0]).open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("naming") != "real":
                continue
            stored, auc, sep = _num(r["radius"]), _num(r[ENDPOINT]), _num(r["separation"])
            if stored is None or stored == UNREACHED or auc is None:
                continue
            g = gac.get((r["shard_id"], r["frame_row_id"]))
            if g is None or g["gac_exact"] != "True" or float(g["gac_radius"]) != stored:
                fails.append(f"GAC radius missing, inexact or different for {r['shard_id']}/{r['frame_row_id']}")
                continue
            if sep is not None:
                continue
            units.append({
                "network": r["network"], "x": r["x"], "y": r["y"], "condition": r["condition"],
                "shard_id": r["shard_id"], "frame_row_id": r["frame_row_id"],
                "radius_stored": int(stored), "AUC_frac": auc,
                "separation_status": r["separation_status"],
                "component_size_csv": _num(r["component_size"]),
                "k_g0_csv": _num(r["k_g0"]), "assumes": r["assumes"],
                "z_star_gac_csv": json.loads(g["z_star"]),
            })
    return units, fails


# --- per-process state ---------------------------------------------------------------------

_BUNDLE: dict[str, Any] = {}


def init_worker(mem_gb: float) -> None:
    if mem_gb > 0:
        cap = int(mem_gb * 1024 ** 3)
        try:
            resource.setrlimit(resource.RLIMIT_DATA, (cap, cap))
        except (ValueError, OSError):
            pass


def get_bundle() -> dict[str, Any]:
    if "b" not in _BUNDLE:
        _BUNDLE["b"] = ls.load_bundle(REPO / INPUTS["bundle"][0])[0]
    return _BUNDLE["b"]


# --- measurements ------------------------------------------------------------------------------


def edges_at(edges, v) -> list[tuple[str, str]]:
    return sorted(e for e in edges if v in e)


def mechanism(h, x, y, z) -> str:
    """Why z fails in h: amenability lost, a forbidden node in z, or (else) a non-causal path left open."""
    if not is_amenable(h, x, y):
        return "not_amenable"
    if set(z) & gac_forbidden_set(h, x, y):
        return "forbidden_node"
    return "open_path"


def radius_with_witness(cpdag, g0, fails):
    """Upward BFS in the order of `radius_local_up_fast`; returns (radius, exact, witness MPDAG or None)."""
    if not fails(cpdag):
        return UNREACHED, True, None
    if fails(g0):
        return 0, True, g0
    seen = {g0.edge_string()}
    frontier = deque([(g0, 0)])
    while frontier:
        cur, d = frontier.popleft()
        for h in local_up_covers_fast(cpdag, cur):
            key = h.edge_string()
            if key in seen:
                continue
            seen.add(key)
            if fails(h):
                return d + 1, True, h
            frontier.append((h, d + 1))
    return UNREACHED, True, None


def measure(cpdag, g0, x, y, z) -> dict[str, Any]:
    k_g0 = frozenset(g0.directed_edges) - frozenset(cpdag.directed_edges)
    k_x = edges_at(k_g0, x)
    comp = next((c for c in undirected_components(cpdag) if x in c), frozenset({x}))
    out: dict[str, Any] = {
        "deg_u_x": len(cpdag.neighbors(x)),
        "deg_x": len(cpdag.neighbors(x)) + len(cpdag.parents(x)) + len(cpdag.children(x)),
        "deg_u_x_g0": len(g0.neighbors(x)),
        "comp_size": len(comp),
        "z_in_component": len(set(z) & set(comp)),
        "k_g0": len(k_g0),
        "n_k_x": len(k_x), "k_x": [list(e) for e in k_x],
        "n_k_x_out": sum(1 for a, _ in k_x if a == x), "n_k_x_in": sum(1 for _, b in k_x if b == x),
        "amenable_at_g0": is_amenable(g0, x, y), "amenable_at_cpdag": is_amenable(cpdag, x, y),
    }

    def fails(g):
        return not is_gac_valid_mpdag(g, x, y, z)

    r, exact, w = radius_with_witness(cpdag, g0, fails)
    out.update(radius=r, radius_exact=exact)
    ref = radius_local_up_fast(cpdag, g0, fails, max_depth=None) if fails(cpdag) else None
    out["radius_ref"] = UNREACHED if ref is None else ref.radius
    out["witness_same_as_ref"] = (w is None and ref is None) or (
        w is not None and ref is not None and w.edge_string() == ref.witness)
    if w is not None and r > 0:
        retracted = sorted(frozenset(g0.directed_edges) - frozenset(w.directed_edges))
        out.update(
            witness_retracted=[list(e) for e in retracted],
            witness_n_retracted=len(retracted),
            witness_n_retracted_at_x=sum(1 for e in retracted if x in e),
            witness_mechanism=mechanism(w, x, y, z),
        )
    # every distance-1 witness: failing upper covers of g0
    covers = local_up_covers_fast(cpdag, g0)
    d1 = []
    for h in covers:
        if fails(h):
            ret = sorted(frozenset(g0.directed_edges) - frozenset(h.directed_edges))
            d1.append({"retracted": [list(e) for e in ret], "n_retracted": len(ret),
                       "n_retracted_at_x": sum(1 for e in ret if x in e),
                       "mechanism": mechanism(h, x, y, z)})
    out["n_covers"] = len(covers)
    out["d1_witnesses"] = d1
    out["n_d1_witnesses"] = len(d1)
    kx = {tuple(e) for e in k_x}
    # registered quantity: X carries exactly one analyst orientation, and a distance-1 witness retracts it
    out["only_kx_retracted"] = bool(len(kx) == 1 and any(kx <= {tuple(e) for e in w_["retracted"]} for w_ in d1))
    out["some_d1_retracts_at_x"] = any(w_["n_retracted_at_x"] > 0 for w_ in d1)
    return out


def unit_key(u: dict[str, Any]) -> str:
    return f"{u['shard_id']}|{u['frame_row_id']}"


def score_item(item: dict[str, Any]) -> dict[str, Any]:
    t0 = time.perf_counter()
    cond, net, cfg = item["condition"], item["network"], item["cfg"]
    bundle = get_bundle()
    dag, cpdag = rs.load_network(net)
    k, g0, reason, shared = ls.build_state(bundle, cond, net, dag, cpdag)
    rows, failures = [], []
    common = {"experiment": EXPERIMENT, "criterion": CRITERION, "commit": cfg["commit"],
              "src_dirty": cfg["src_dirty"], "smoke": cfg["smoke"], "g0_status": reason,
              "k_g0_state": shared.get("k_g0")}
    for u in item["units"]:
        row = {**common, **u, "key": unit_key(u)}
        x, y = u["x"], u["y"]
        if g0 is None:
            row["status"] = f"g0_{reason}"
            failures.append({"kind": "g0_missing", "key": row["key"], "reason": reason})
            rows.append(row)
            continue
        z = frozenset(u["z_star_gac_csv"])
        z2, zr = rs.commit_z_star(g0, x, y)
        row["z_matches_commit"] = z2 is not None and z2 == z
        if not row["z_matches_commit"]:
            failures.append({"kind": "z_differs_from_commit_z_star", "key": row["key"],
                             "csv": sorted(z), "commit": sorted(z2) if z2 is not None else zr})
        signal.signal(signal.SIGALRM, _alarm)
        signal.alarm(CALL_WALL_S)
        t = time.perf_counter()
        try:
            row.update(measure(cpdag, g0, x, y, z))
            row["status"] = "ok"
        except Wall:
            row["status"] = "wall"
        except MemoryError:
            row["status"] = "memory"
        finally:
            signal.alarm(0)
        row["t_s"] = round(time.perf_counter() - t, 4)
        if row["status"] == "ok":
            if row["radius"] != u["radius_stored"]:
                failures.append({"kind": "radius_differs_from_stored", "key": row["key"],
                                 "got": row["radius"], "stored": u["radius_stored"]})
            if row["radius_ref"] != row["radius"] or not row["witness_same_as_ref"]:
                failures.append({"kind": "radius_or_witness_differs_from_radius_local_up_fast",
                                 "key": row["key"], "got": row["radius"], "ref": row["radius_ref"]})
            if (row["radius"] == 1) != (row["n_d1_witnesses"] > 0):
                failures.append({"kind": "radius_one_vs_d1_witnesses", "key": row["key"]})
            if row["k_g0"] != u["k_g0_csv"]:
                failures.append({"kind": "k_g0_differs_from_csv", "key": row["key"],
                                 "got": row["k_g0"], "csv": u["k_g0_csv"]})
        rows.append(row)
    return {"rows": rows, "failures": failures, "elapsed_s": time.perf_counter() - t0}


def safe_score(item: dict[str, Any]) -> dict[str, Any]:
    try:
        return score_item(item)
    except Exception:
        return {"error": traceback.format_exc(), "item": f"{item['condition']}|{item['network']}"}


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    out = []
    for line in path.read_text().splitlines():
        try:
            out.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smoke", action="store_true", help=f"first {SMOKE_ITEMS} (condition, network) items")
    ap.add_argument("--out", default=None)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--mem-gb", type=float, default=2.0, help="data-memory cap per worker")
    args = ap.parse_args(argv)
    if args.workers > 6:
        raise SystemExit("at most 6 workers")
    if not 0 <= args.shard < args.n_shards:
        ap.error("need 0 <= --shard < --n-shards")
    out = Path(args.out or REPO / ("results/rebuttal/E29_smoke" if args.smoke else "results/rebuttal/E29"))
    out.mkdir(parents=True, exist_ok=True)

    hashes = {k: sha_file(REPO / p) for k, (p, _) in INPUTS.items()}
    bad = [k for k, (_, want) in INPUTS.items() if hashes[k] != want]
    if bad:
        print(f"input hash mismatch: {bad}", file=sys.stderr)
        return 3
    units, sel_fails = select_population()
    if sel_fails or len(units) != EXPECTED_N:
        for f in sel_fails[:20]:
            print("CONSISTENCY FAILURE:", f, flush=True)
        print(f"population has {len(units)} units, expected {EXPECTED_N}; stopping", file=sys.stderr)
        return 2
    by_item: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for u in units:
        by_item.setdefault((u["condition"], u["network"]), []).append(u)
    keys = sorted(by_item)
    if args.smoke:
        keys = keys[:SMOKE_ITEMS]
    keys = [k for i, k in enumerate(keys) if i % args.n_shards == args.shard]

    commit, src_dirty = repo_commit()
    cfg = {"commit": commit, "src_dirty": src_dirty, "smoke": args.smoke}
    done = {r["key"] for f in sorted(out.glob("units*.jsonl")) for r in read_jsonl(f) if r.get("status") == "ok"}
    pending = []
    for k in keys:
        todo = [u for u in by_item[k] if unit_key(u) not in done]
        if todo:
            pending.append({"condition": k[0], "network": k[1], "units": todo, "cfg": cfg})
    n_units = sum(len(it["units"]) for it in pending)
    tag = f"shard-{args.shard}-of-{args.n_shards}"
    config = {"experiment": EXPERIMENT, "criterion": CRITERION, "endpoint": ENDPOINT, "smoke": args.smoke,
              "population": "E16 separation-undefined units, submitted cohort", "n_population": len(units),
              "n_items_total": len(by_item), "input_sha256": hashes, "hash_ok": not bad,
              "commit": commit, "src_dirty": src_dirty, "call_wall_s": CALL_WALL_S,
              "mem_gb_per_worker": args.mem_gb, "workers": args.workers, "shard": args.shard,
              "n_shards": args.n_shards, "n_items_pending": len(pending), "n_units_pending": n_units,
              "n_units_already_done": len(done),
              "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    (out / f"config.{tag}.json").write_text(json.dumps(config, indent=1))
    print(f"E29 {tag}: {len(pending)} items, {n_units} units pending ({len(done)} done), commit={commit[:10]}",
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
        print(f"[{time.time() - t_start:8.0f}s] {rows[0]['condition'] if rows else '?'}|"
              f"{rows[0]['network'] if rows else '?'}: {len(rows)} units, ok={len(ok)}, "
              f"failures={len(res['failures'])}, {res['elapsed_s']:.1f}s", flush=True)

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
