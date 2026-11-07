#!/usr/bin/env python3
"""E23: PC and GES agreement as a check for discovery error, on the stored learned graphs.

For each certifiable A_ORACLE_R row of a learned graph (E5's arm), the committed set Z is checked
in the state built from the other method's graph learned on the same sample, with the same knowledge
rule. One JSON line per row. See PREREG.md.

    .venv/bin/python experiments/rebuttal/E23_e23_pc_ges_disagreement_check/run.py --smoke
    .venv/bin/python experiments/rebuttal/E23_e23_pc_ges_disagreement_check/run.py --workers 2

Exit 0 when every consistency check holds, 2 when one fails (rows are still written), 1 on error.
"""

from __future__ import annotations

import argparse
import hashlib
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

from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag  # noqa: E402
from vendor_bkrobust import committed_adjustment_set, is_valid_adjustment_set_gac_dag, sha  # noqa: E402
from vendor_e5 import NETWORKS, canon, file_sha, greedy_consistent, load_inputs, plans_transfer, repo_commit  # noqa: E402

CRITERION = "gac"
ARM = "A_ORACLE_R"
N_GRID = (200, 1000, 5000)
OTHER = {"pc": "ges", "ges": "pc"}
SMOKE_NETWORKS = ("ecoli70", "magic-niab")
SMOKE_SEEDS = (0,)
SMOKE_N = (200, 5000)


def k_hash(k: list) -> str:
    """As in E5 run.py."""
    return hashlib.sha256(json.dumps(sorted(k)).encode()).hexdigest()[:16]


_INP: dict = {}
_CFG: dict = {}
_PLAN: dict = {}


def init_worker(cfg: dict) -> None:
    _CFG.update(cfg)
    _INP.update(load_inputs(Path(cfg["data_dir"]), Path(cfg["models_dir"]), tuple(cfg["networks"])))


def oracle_state(gid: str):
    """E5's A_ORACLE_R knowledge on this graph and its closed state (cached per graph)."""
    if gid not in _PLAN:
        plan = [p for p in plans_transfer(_INP, gid, 0) if p[0] == ARM]
        assert len(plan) == 1, gid
        k = [tuple(e) for e in plan[0][2]]
        _PLAN[gid] = (k, apply_orientations(_INP["graphs"][gid]["g"], k))
    return _PLAN[gid]


def run_triple(task: dict) -> list[dict]:
    """Every certifiable row of the PC and GES graphs of one (network, n, seed)."""
    net, n, seed = task["net"], task["n"], task["seed"]
    dag = _INP["dags"][net]
    gids = {m: f"{net}|{m}|{n}|{seed}" for m in ("pc", "ges")}
    present = {m: g for m, g in gids.items() if g in _INP["graphs"]}
    pairs = [(r["x"], r["y"]) for r in _INP["frame"] if r["network"] == net]
    out = []
    for m, gid in sorted(present.items()):
        k, g0 = oracle_state(gid)
        partner = present.get(OTHER[m])
        for x, y in pairs:
            unit = f"{gid}|{x}|{y}"
            if unit in task["done"]:
                continue
            t0 = time.perf_counter()
            if g0 is None:
                continue  # knowledge_inconsistent: not certifiable, not a row of the population
            cs = committed_adjustment_set(g0, x, y)
            if cs.z is None:
                continue
            z = cs.z
            rec = {"unit": unit, "gid": gid, "network": net, "n": n, "seed": seed, "method": m,
                   "direction": f"{m}->{OTHER[m]}", "x": x, "y": y, "arm": ARM,
                   "e5_unit": f"B|{gid}|{ARM}|0|{x}|{y}", "k_sha256": k_hash(k), "n_k": len(k),
                   "z": sorted(z), "z_size": len(z), "committed_verdict": cs.verdict,
                   "z_valid_at_truth": int(is_valid_adjustment_set_gac_dag(dag, x, y, z)),
                   "criterion": CRITERION, "repo_commit": _CFG["repo_commit"], "partner": partner}
            fails = []
            if not is_gac_valid_mpdag(g0, x, y, z):
                fails.append("z_not_valid_in_own_state")
            if partner is not None:
                k2, g2 = oracle_state(partner)
                rec["partner_n_k"] = len(k2)
                rec["partner_state_ok"] = int(g2 is not None)
                rec["same_graph"] = int(sha(_INP["graphs"][gid]["g"]) == sha(_INP["graphs"][partner]["g"]))
                if g2 is None:
                    rec.update(z_valid_other=0, z_same_other=0, other_verdict="knowledge_inconsistent")
                else:
                    cs2 = committed_adjustment_set(g2, x, y)
                    rec["z_valid_other"] = int(is_gac_valid_mpdag(g2, x, y, z))
                    rec["z_same_other"] = int(cs2.z is not None and cs2.z == z)
                    rec["other_verdict"] = cs2.verdict
                rec["flag"] = 1 - rec["z_valid_other"]
                rec["flag_or_diff"] = int(rec["flag"] or not rec["z_same_other"])
                # secondary: the source knowledge itself, moved to the partner graph
                gp = _INP["graphs"][partner]["g"]
                moved, dropped = greedy_consistent(gp, [e for e in k if canon(e) in gp.undirected_edges])
                g3 = apply_orientations(gp, moved)
                rec["moved_n_k"], rec["moved_dropped"] = len(moved), dropped
                rec["z_valid_moved"] = int(g3 is not None and is_gac_valid_mpdag(g3, x, y, z))
                rec["flag_moved"] = 1 - rec["z_valid_moved"]
                if rec["same_graph"] and rec["flag"]:
                    fails.append("identical_graphs_but_flagged")
            rec["seconds"] = round(time.perf_counter() - t0, 3)
            rec["consistency_fail"] = fails
            out.append(rec)
    return out


def done_units(out: Path) -> dict:
    rows = {}
    for p in sorted(out.glob("rows.*.jsonl")):
        for line in p.open():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue  # a line cut by a crash; recomputed
            rows[r["unit"]] = r
    return rows


def e5_check(rows: dict, e5_dir: Path, gids: set) -> tuple[list, int]:
    """Consistency check 1 against E5's stored rows, restricted to the graphs this run covered."""
    files = sorted(e5_dir.glob("rows.*.jsonl"))
    if not files:
        return [], 0
    e5 = {}
    for p in files:
        for line in p.open():
            r = json.loads(line)
            if r["phase"] == "B" and r["arm"] == ARM and r["gid"] in gids and r["status"] == "ok":
                e5[r["unit"]] = r
    mine = {r["e5_unit"]: r for r in rows.values() if r["gid"] in gids}
    fails = [(u, "in_e5_not_rebuilt") for u in sorted(set(e5) - set(mine))]
    fails += [(u, "rebuilt_not_in_e5") for u in sorted(set(mine) - set(e5))]
    for u in sorted(set(e5) & set(mine)):
        for f in ("k_sha256", "z_size", "z_valid_at_truth", "committed_verdict"):
            if e5[u].get(f) != mine[u].get(f):
                fails.append((u, f"{f}:{mine[u].get(f)!r}!={e5[u].get(f)!r}"))
    return fails, len(e5)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--smoke", action="store_true", help="two networks, seed 0, n 200 and 5000, all 20 queries each")
    ap.add_argument("--out", default="results/rebuttal/E23")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--data-dir", default="experiments/rebuttal/_data/pest")
    ap.add_argument("--models-dir", default="results/axisa3/networks/example_models")
    ap.add_argument("--e5-dir", default="results/rebuttal/E5", help="E5 rows for the rebuild check")
    args = ap.parse_args(argv)
    if args.workers > 6:
        ap.error("at most 6 workers")
    if not 0 <= args.shard < args.n_shards:
        ap.error("--shard must be in [0, --n-shards)")

    def rel(p):
        p = Path(p)
        return p if p.is_absolute() else REPO / p

    data_dir, models_dir, out, e5_dir = rel(args.data_dir), rel(args.models_dir), rel(args.out), rel(args.e5_dir)
    if args.smoke and args.out == "results/rebuttal/E23":
        out = REPO / "results/rebuttal/E23_smoke"
    networks = SMOKE_NETWORKS if args.smoke else NETWORKS
    out.mkdir(parents=True, exist_ok=True)
    e5_files = sorted(e5_dir.glob("rows.*.jsonl"))
    cfg = {
        "experiment": "E23_e23_pc_ges_disagreement_check", "criterion": CRITERION, "arm": ARM,
        "repo_commit": repo_commit(), "run_sha256": file_sha(Path(__file__)),
        "vendor_sha256": {p.name: file_sha(p) for p in sorted(HERE.glob("vendor_*.py"))},
        "smoke": bool(args.smoke), "networks": list(networks),
        "n_grid": list(SMOKE_N if args.smoke else N_GRID), "methods": ["pc", "ges"],
        "data_dir": str(data_dir), "models_dir": str(models_dir), "e5_dir": str(e5_dir),
        "input_sha256": {p.name.replace(".", "_"): file_sha(p) for p in sorted(data_dir.iterdir()) if p.is_file()},
        "e5_rows_sha256": {p.name.replace(".", "_"): file_sha(p) for p in e5_files},
        "shard": args.shard, "n_shards": args.n_shards, "workers": args.workers,
    }
    tag = f"shard{args.shard}of{args.n_shards}"
    (out / f"config.{tag}.json").write_text(json.dumps(cfg, indent=1))

    init_worker(cfg)
    triples = sorted({(i["net"], i["n"], i["seed"]) for i in _INP["graphs"].values()
                      if i["kind"] == "learned" and i["n"] in cfg["n_grid"]})
    if args.smoke:
        triples = [t for t in triples if t[2] in SMOKE_SEEDS]
    triples = [t for i, t in enumerate(triples) if i % args.n_shards == args.shard]
    done = set(done_units(out))
    tasks = [{"net": a, "n": b, "seed": c, "done": done} for a, b, c in triples]
    print(f"E23 {tag}: {len(tasks)} (network, n, seed) triples, {len(done)} rows already written, "
          f"{args.workers} workers, criterion {CRITERION}, commit {cfg['repo_commit'][:12]}", flush=True)

    t0 = time.perf_counter()
    written = 0
    with (out / f"rows.{tag}.jsonl").open("a") as fh:
        def write(rows):
            nonlocal written
            for r in rows:
                fh.write(json.dumps(r) + "\n")
            fh.flush()
            written += len(rows)

        if args.workers <= 1:
            for i, t in enumerate(tasks, 1):
                write(run_triple(t))
                print(f"  {i}/{len(tasks)} triples, {written} rows, {time.perf_counter() - t0:.0f} s", flush=True)
        else:
            with mp.get_context("spawn").Pool(args.workers, initializer=init_worker, initargs=(cfg,)) as pool:
                for i, rows in enumerate(pool.imap_unordered(run_triple, tasks), 1):
                    write(rows)
                    print(f"  {i}/{len(tasks)} triples, {written} rows, {time.perf_counter() - t0:.0f} s", flush=True)

    rows = done_units(out)
    fails = [(u, r["consistency_fail"]) for u, r in rows.items() if r["consistency_fail"]]
    covered = {f"{a}|{m}|{b}|{c}" for a, b, c in triples for m in ("pc", "ges")}
    if args.n_shards == 1:
        e5_fails, n_e5 = e5_check(rows, e5_dir, covered)
        print(f"E5 rebuild check: {n_e5} E5 rows, {len(e5_fails)} disagreements" if e5_files
              else "E5 rebuild check: E5 rows absent, NOT checked")
    else:
        e5_fails, n_e5 = [], 0
        print("E5 rebuild check: left to analyse.py on a sharded run")
    fails += e5_fails
    (out / f"e5_check.{tag}.json").write_text(json.dumps({"e5_rows": n_e5, "fails": e5_fails,
                                                          "e5_present": bool(e5_files)}, indent=1))
    print(f"done: {written} rows written in {time.perf_counter() - t0:.0f} s, {len(rows)} in total; "
          f"{len(fails)} consistency failures", flush=True)
    if fails:
        (out / f"consistency_failures.{tag}.json").write_text(json.dumps(fails, indent=1))
        for u, f in fails[:20]:
            print(f"  FAIL {u}: {f}", flush=True)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
