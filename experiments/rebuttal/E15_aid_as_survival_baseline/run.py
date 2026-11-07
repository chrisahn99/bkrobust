"""E15: adjustment identification distance (AID) beside the breakdown radius.

Regenerates every stored survival instance from its seed, checks it reproduces
the stored x, Z*, radius and separation, and writes per instance the GAC radius
and the AID from G0 to the true DAG (oracle) and to the CPDAG (analyst), ours and
gadjid's. Survival (AUC_frac) is carried over from the stored table, not redrawn.

    PYTHONPATH=src .venv/bin/python experiments/rebuttal/E15_aid_as_survival_baseline/run.py --smoke
    PYTHONPATH=src .venv/bin/python experiments/rebuttal/E15_aid_as_survival_baseline/run.py --workers 8

Exit codes: 0 ok, 3 consistency failure (rows are still written, flagged), 4 gadjid missing
on a full run without --allow-no-gadjid.
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
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(HERE))

from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.hybrid import breakdown_radius  # noqa: E402
from bkrobust.robustness import survival as sv  # noqa: E402
from bkrobust.robustness import unconfounded_generator as ug  # noqa: E402
from bkrobust.synth.component_generator import (  # noqa: E402
    component_containing,
    distances_within_component,
)

import aid  # noqa: E402

CRITERION = "gac"
UNTARGETED_CSV = "results/axis_robustness_unconfounded_v2/survival_instances.csv"
TARGETED_CSV = "results/axis_robustness_p6/survival_instances.csv"
# Stored columns carried into each row so analyse.py needs nothing else.
CARRY = ["arm", "component_size", "n_tiers", "coverage", "base_wrongness", "seed",
         "separation", "separation_status", "n_k", "shd_truth", "shd_cpdag", "k_g0",
         "r_val", "AUC_frac", "AUC_abs", "s0_ok"]

_GADJID: dict = {}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repo_commit() -> tuple[str, bool]:
    try:
        sha = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                             capture_output=True, text=True, check=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "-C", str(REPO), "status", "--porcelain", "--", "src"],
                                    capture_output=True, text=True).stdout.strip())
        return sha, dirty
    except Exception:  # noqa: BLE001
        return "unknown", True


def load_units(untargeted: Path, targeted: Path, smoke: bool) -> list[dict]:
    units = []
    for dataset, path in (("untargeted", untargeted), ("targeted", targeted)):
        with open(path) as f:
            rows = list(csv.DictReader(f))
        if smoke:
            # Two rows per arm from the smallest component size: deterministic and cheap,
            # and it touches all three builders.
            picked, seen = [], {}
            for r in rows:
                if r["component_size"] == "6" and seen.get(r["arm"], 0) < 2:
                    picked.append(r)
                    seen[r["arm"]] = seen.get(r["arm"], 0) + 1
            rows = picked
        for r in rows:
            units.append({"dataset": dataset, "unit_id": f"{dataset}|{r['instance_id']}", "row": r})
    return units


def rebuild(row: dict):
    arm = row["arm"]
    c, seed = int(row["component_size"]), int(row["seed"])
    if arm == "tiered_unconfounded":
        return ug.build_unconfounded_tiered_instance(c, seed, int(row["n_tiers"]))
    if arm == "flip":
        return sv.build_flip_instance(c, int(row["separation"]), float(row["coverage"]),
                                      float(row["base_wrongness"]), seed)
    if arm == "tiered":
        return sv.build_tiered_instance(c, int(row["separation"]), seed, int(row["n_tiers"]))
    return None, f"unknown_arm:{arm}"


def measured_separation(cpdag, x, z_star):
    """BFS distance from x to the nearest Z* member inside x's undirected component
    of the CPDAG; the rule unconfounded_generator._measure_separation applies."""
    comp = component_containing(cpdag, x)
    if comp is None:
        return None
    dist = distances_within_component(cpdag, comp, x)
    reach = [dist[v] for v in z_star if v in comp and v in dist]
    return min(reach) if reach else None


def _worker_init(gadjid_cfg: dict) -> None:
    _GADJID.update(gadjid_cfg)


def process(unit: dict, time_limit_s: float) -> dict:
    row = unit["row"]
    out = {"unit_id": unit["unit_id"], "dataset": unit["dataset"],
           "instance_id": row["instance_id"], "criterion": CRITERION}
    out.update({f"stored_{k}": row.get(k, "") for k in CARRY})
    fails: list[str] = []
    t0 = time.perf_counter()
    try:
        inst, reason = rebuild(row)
    except Exception as exc:  # noqa: BLE001
        inst, reason = None, f"exception:{exc!r}"
    if inst is None:
        out.update(status="rebuild_failed", reason=reason, consistency_failures=["rebuild_failed"])
        return out
    cpdag, g0, dag, x, y, z = inst["cpdag"], inst["g0"], inst["dag"], inst["x"], inst["y"], inst["z_star"]

    z_stored = frozenset(s for s in row["z_star"].split(",") if s)
    if x != row["x"]:
        fails.append("x_mismatch")
    if z != z_stored:
        fails.append("z_mismatch")
    sep_meas = measured_separation(cpdag, x, z)
    stored_sep = int(row["separation"]) if row["separation"] not in ("", None) else None
    out.update(x=x, y=y, z_star=sorted(z), n_nodes=len(dag.nodes),
               separation_rebuilt=inst.get("separation"), separation_measured=sep_meas)
    # On the untargeted arm the stored separation was measured this way; on the
    # targeted arms it is the generator's requested value, recorded only.
    if unit["dataset"] == "untargeted" and inst.get("separation") != stored_sep:
        fails.append("separation_mismatch")
    out["separation_measured_equals_stored"] = sep_meas == stored_sep

    res = breakdown_radius(cpdag, None, x, y, z, g0=g0, criterion=CRITERION, time_limit_s=time_limit_s)
    out.update(radius_gac=res.radius, radius_exact=res.exact, radius_method=res.method,
               radius_oracle=res.oracle, radius_seconds=res.total_seconds)
    stored_r = int(row["r_val"]) if row["r_val"] not in ("", None) else None
    if not res.exact:
        fails.append("radius_not_exact")
    elif res.radius != stored_r:
        fails.append("radius_mismatch")

    try:
        answers = aid.guess_answers(g0)
        out["guess_answers"] = aid.answer_summary(answers)
        for ref_name, ref in (("oracle", dag), ("cpdag", cpdag)):
            for k, v in aid.aid_scopes(answers, ref, x, y).items():
                out[f"aid_{ref_name}_{k}"] = v
        # Both say "Z* fails in some member of the class": the pair AID to the
        # CPDAG and a finite radius.
        if out["aid_cpdag_pair"] != int(res.radius != UNREACHED):
            fails.append("pair_aid_vs_radius")
        out["aid_status"] = "ok"
    except Exception as exc:  # noqa: BLE001
        out["aid_status"] = f"error:{exc!r}"[:300]
        fails.append("aid_error")

    if _GADJID.get("available"):
        sym = _GADJID["symmetric_undirected"]
        for ref_name, ref in (("oracle", dag), ("cpdag", cpdag)):
            g = aid.gadjid_aids(ref, g0, sym)
            out.update({f"gadjid_{ref_name}_{k}": v for k, v in g.items()})
            if "oset_aid_count" in g and out.get("aid_status") == "ok":
                out[f"gadjid_{ref_name}_oset_agrees"] = g["oset_aid_count"] == out[f"aid_{ref_name}_graph_count"]
    out["status"] = "ok" if not fails else "consistency_failure"
    out["consistency_failures"] = fails
    out["seconds"] = time.perf_counter() - t0
    return out


def _task(args):
    unit, time_limit_s = args
    return process(unit, time_limit_s)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smoke", action="store_true", help="two units per arm, component size 6")
    ap.add_argument("--out", default=None, help="default results/rebuttal/E15 (E15_smoke with --smoke)")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--untargeted-csv", default=UNTARGETED_CSV, help="relative to the repository root")
    ap.add_argument("--targeted-csv", default=TARGETED_CSV, help="relative to the repository root")
    ap.add_argument("--time-limit", type=float, default=300.0, help="per-rung limit of the radius ladder")
    ap.add_argument("--allow-no-gadjid", action="store_true",
                    help="run a full pass without gadjid (its columns are then absent)")
    args = ap.parse_args(argv)

    out_dir = Path(args.out or ("results/rebuttal/E15_smoke" if args.smoke else "results/rebuttal/E15"))
    if not out_dir.is_absolute():
        out_dir = REPO / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = [Path(p) if Path(p).is_absolute() else REPO / p for p in (args.untargeted_csv, args.targeted_csv)]

    gadjid_cfg = aid.gadjid_probe()
    if not gadjid_cfg.get("available") and not (args.smoke or args.allow_no_gadjid):
        print(f"gadjid unavailable ({gadjid_cfg.get('error')}); install it or pass --allow-no-gadjid",
              file=sys.stderr)
        return 4

    commit, dirty = repo_commit()
    units = load_units(paths[0], paths[1], args.smoke)
    n_units_total = len(units)
    units = [u for i, u in enumerate(units) if i % args.n_shards == args.shard]
    out_file = out_dir / f"units_shard{args.shard:03d}of{args.n_shards:03d}.jsonl"
    done, prior_fail = set(), set()
    for f in out_dir.glob("units_shard*.jsonl"):
        for line in f.read_text().splitlines():
            if line.strip():
                rec = json.loads(line)
                done.add(rec["unit_id"])
                if rec.get("status") != "ok":
                    prior_fail.add(rec["unit_id"])
    todo = [u for u in units if u["unit_id"] not in done]
    # A resumed run must still fail on rows an earlier pass flagged.
    n_prior_fail = sum(1 for u in units if u["unit_id"] in prior_fail)

    config = {
        "experiment": "E15", "criterion": CRITERION, "repo_commit": commit, "src_dirty": dirty,
        "smoke": args.smoke, "shard": args.shard, "n_shards": args.n_shards, "workers": args.workers,
        "untargeted_csv": str(paths[0]), "untargeted_sha256": sha256(paths[0]),
        "targeted_csv": str(paths[1]), "targeted_sha256": sha256(paths[1]),
        "n_units_total": n_units_total, "n_units_shard": len(units), "n_todo": len(todo),
        "gadjid": gadjid_cfg, "gadjid_available": bool(gadjid_cfg.get("available")),
        "time_limit_s": args.time_limit, "python": sys.version.split()[0],
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    (out_dir / f"config_shard{args.shard:03d}of{args.n_shards:03d}.json").write_text(json.dumps(config, indent=2))
    print(f"E15: {len(units)} units in shard {args.shard}/{args.n_shards}, {len(todo)} to do, "
          f"gadjid={gadjid_cfg.get('available')}, commit {commit[:10]}", flush=True)

    n_fail, n_done, t_start, t_print = 0, 0, time.time(), time.time()
    tasks = [(u, args.time_limit) for u in todo]
    with open(out_file, "a") as fh:
        def write(rec):
            nonlocal n_fail, n_done
            rec.update(repo_commit=commit, criterion=CRITERION)
            fh.write(json.dumps(rec, default=str) + "\n")
            fh.flush()
            n_done += 1
            if rec["status"] != "ok":
                n_fail += 1
                print(f"  CONSISTENCY {rec['unit_id']}: {rec.get('consistency_failures')} {rec.get('reason', '')}",
                      flush=True)

        if args.workers <= 1:
            _worker_init(gadjid_cfg)
            for t in tasks:
                write(_task(t))
                if time.time() - t_print > 30:
                    t_print = time.time()
                    print(f"  {n_done}/{len(todo)} done, {n_fail} failures, {time.time() - t_start:.0f}s", flush=True)
        else:
            with mp.get_context("spawn").Pool(args.workers, initializer=_worker_init,
                                              initargs=(gadjid_cfg,)) as pool:
                it = pool.imap_unordered(_task, tasks, chunksize=1)
                while n_done < len(tasks):
                    try:
                        write(it.next(timeout=30))
                    except mp.TimeoutError:
                        pass
                    if time.time() - t_print > 30:
                        t_print = time.time()
                        print(f"  {n_done}/{len(todo)} done, {n_fail} failures, "
                              f"{time.time() - t_start:.0f}s", flush=True)

    print(f"E15: wrote {n_done} rows to {out_file}, {n_fail} consistency failures, "
          f"{time.time() - t_start:.0f}s", flush=True)
    if n_prior_fail:
        print(f"E15: {n_prior_fail} consistency failures in rows written by an earlier pass", flush=True)
    return 3 if (n_fail or n_prior_fail) else 0


if __name__ == "__main__":
    sys.exit(main())
