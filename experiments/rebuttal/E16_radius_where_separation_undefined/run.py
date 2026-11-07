"""E16: does the radius rank survival on the Table 1 units where separation is undefined?

Pure re-analysis of committed per-unit tables (no solver is called): the radius used is the GAC
radius of results/radius_gac_recompute/, checked equal to the stored radius on every scorable unit.
One JSON line per analysis cell (cohort x statistic). See PREREG.md.

    PYTHONPATH=src .venv/bin/python experiments/rebuttal/E16_radius_where_separation_undefined/run.py --smoke
    PYTHONPATH=src .venv/bin/python experiments/rebuttal/E16_radius_where_separation_undefined/run.py --workers 4
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import multiprocessing as mp
import subprocess
import sys
import threading
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(REPO_ROOT / "experiments"))

from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.robustness.paired_resample import _kendall_tau_b, paired_comparison  # noqa: E402
from bkrobust.robustness.real_analyse import tau_for_stratum  # noqa: E402

import llm_paired_v2 as base  # noqa: E402  (panel definitions and the loop reference of Table 1)

CRITERION = "gac"
ENDPOINT = "AUC_frac"
SEED = 0
FULL_N_BOOT = 10_000
SMOKE_N_BOOT = 200
COHORTS = ("submitted", "corrected")
CELLS = (
    "describe_undefined",
    "pooled_tau_radius_undefined",      # registered (submitted cohort)
    "within_tau_radius_undefined",
    "pooled_tau_radius_all",
    "pooled_tau_radius_defined",
    "pooled_tau_separation_defined",
    "pooled_paired_radius_minus_separation_defined",
    "within_tau_radius_all",
    "within_paired_radius_minus_separation_matched",  # secondary (corrected cohort)
)

# Published values these cells must reproduce (results/axis_robustness_llm_v2/*.json).
# Points are hard checks (same data, same code); intervals are soft (tolerance 0.02).
REFERENCE = {
    ("submitted", "pooled_tau_radius_all"): {"point": 0.31906985313150255, "n": 2422, "ci": (0.122, 0.485)},
    ("submitted", "pooled_tau_separation_defined"): {"point": 0.2837657748918398, "n": 1389,
                                                      "ci": (-0.028961397095442657, 0.49215872363025825)},
    ("submitted", "pooled_paired_radius_minus_separation_defined"): {
        "point": 0.0386086177695657, "n": 1389, "ci": (-0.18151365039418357, 0.2838412156272436)},
    ("submitted", "within_tau_radius_all"): {"point": 0.4191447599843434, "n": 1973,
                                              "ci": (0.05496887195304333, 0.6788482254279423)},
    ("submitted", "within_paired_radius_minus_separation_matched"): {
        "point": 0.18707401519160471, "n": 1213, "ci": (0.06633205125485178, 0.8948111786895783),
        "tau_a": 0.678428355109646, "tau_b": 0.49135433991804134},
    ("submitted", "pooled_tau_radius_undefined"): {"n": 1033},
    ("corrected", "pooled_tau_radius_all"): {"point": 0.32400745110548795, "n": 2516},
    ("corrected", "within_tau_radius_all"): {"point": 0.42938618097202336, "n": 2067,
                                              "ci": (0.07128736910518067, 0.6867286421502002)},
    ("corrected", "within_paired_radius_minus_separation_matched"): {
        "n": 1307, "tau_a": 0.6834986035430098, "tau_b": 0.4537640473531176},
    ("corrected", "pooled_tau_radius_undefined"): {"n": 1033},
}
POINT_TOL = 1e-6
CI_TOL = 0.02

# The 94 empty-valid rows carry a newer wording of the same `assumes` label; tau_for_stratum
# refuses mixed labels, so normalise exactly as experiments/table1_corrected_cohort.py does.
_COMMITTED_ASSUMES = "Conjecture 2 (hence Anti-Exchange Case B, verified not proved)"
_CURRENT_ASSUMES = "Conjecture 2 (proved: Anti-Exchange Case B, THEOREMS.md section 4)"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_info() -> tuple[str, bool | None]:
    try:
        head = subprocess.run(["git", "-C", str(REPO_ROOT), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(REPO_ROOT), "status", "--porcelain", "--", "src"],
                               capture_output=True, text=True, check=True).stdout.strip() != ""
        return head, dirty
    except (OSError, subprocess.CalledProcessError):
        return "unknown", None


def load_gac(gac_dir: Path) -> dict[tuple[str, str], dict]:
    out = {}
    for name in ("axis_robustness_llm.csv", "axis_robustness_llm_large.csv",
                 "axis_robustness_llm_v2_empty_valid.csv"):
        with (gac_dir / name).open(encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                if r["gac_radius"] == "":
                    continue  # not recomputed in that file; the large file covers those units
                out[(r["shard_id"], r["frame_row_id"])] = r
    return out


def build_cohorts(args) -> tuple[dict[str, list[dict]], list[str]]:
    """Real-name rows of both cohorts with `radius` replaced by the GAC radius.

    Returns the cohorts and the list of consistency failures found while joining."""
    old = base.load_units(Path(args.units_csv))
    new = [u for u in base.load_units(Path(args.empty_valid_csv)) if u.get("naming") == "real"]
    for u in new:
        if u.get("assumes") == _CURRENT_ASSUMES:
            u["assumes"] = _COMMITTED_ASSUMES
    gac = load_gac(Path(args.gac_dir))
    failures = []
    cohorts = {"submitted": base.real_naming(old), "corrected": base.real_naming(old) + new}
    seen = set()
    for rows in cohorts.values():
        for u in rows:
            if id(u) in seen:
                continue
            seen.add(id(u))
            stored = u["radius"]
            u["radius_stored_backdoor"] = stored
            g = gac.get((u["shard_id"], u["frame_row_id"]))
            u["radius"] = None if g is None else float(g["gac_radius"])
            scorable = stored is not None and stored != UNREACHED and u.get(ENDPOINT) is not None
            if not scorable:
                continue
            if g is None:
                failures.append(f"no GAC radius for scorable unit {u['shard_id']}/{u['frame_row_id']}")
            elif g["gac_exact"] != "True":
                failures.append(f"GAC radius not exact for {u['shard_id']}/{u['frame_row_id']}")
            elif float(g["gac_radius"]) != stored:
                failures.append(f"GAC radius {g['gac_radius']} != stored {stored} "
                                f"for {u['shard_id']}/{u['frame_row_id']}")
            sep_measured = u.get("separation_status") == "measured"
            if sep_measured != (u.get("separation") is not None):
                failures.append(f"separation value and status disagree for {u['shard_id']}/{u['frame_row_id']}")
    if len(new) != 94:
        failures.append(f"expected 94 real-name empty-valid units, got {len(new)}")
    return cohorts, failures


def scorable(rows: list[dict]) -> list[dict]:
    return [u for u in rows if u.get("radius") is not None and u["radius"] != UNREACHED
            and u.get(ENDPOINT) is not None]


def radius_varying_states(rows: list[dict]) -> dict[tuple, list[dict]]:
    """Table 1's within-state population: scorable rows grouped by (network, condition),
    keeping states where the radius varies (llm_paired_v2.task4_within_state)."""
    by_state = defaultdict(list)
    for u in scorable(rows):
        by_state[(u["network"], u["condition"])].append(u)
    return {k: v for k, v in by_state.items() if len({r["radius"] for r in v}) > 1}


def restrict(states: dict[tuple, list[dict]], keep) -> dict[tuple, list[dict]]:
    out = {}
    for k, rows in states.items():
        sub = [r for r in rows if keep(r)]
        if len(sub) >= 2:
            out[k] = sub
    return out


# ------------------------------------------------------------------ within-state tau_b, vectorised
# A resample's stratified tau_b is a sum of fixed per-state pair counts, so summing each drawn
# network's counts gives exactly what llm_paired_v2._within_state_bootstrap computes pair by pair.

def pair_counts(xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
    i, j = np.triu_indices(len(xs), 1)
    dx, dy = xs[i] - xs[j], ys[i] - ys[j]
    keep = ~((dx == 0) & (dy == 0))
    dx, dy = dx[keep], dy[keep]
    tx = np.sum(dx == 0)
    ty = np.sum((dx != 0) & (dy == 0))
    both = (dx != 0) & (dy != 0)
    c = np.sum(both & (dx * dy > 0))
    d = np.sum(both & (dx * dy < 0))
    return np.array([c, d, tx, ty], dtype=np.int64)


def tau_from_counts(cnt: np.ndarray) -> float | None:
    c, d, tx, ty = (int(v) for v in cnt)
    denom = math.sqrt((c + d + tx) * (c + d + ty))
    return (c - d) / denom if denom > 0 else None


def network_count_matrix(states, col):
    networks = sorted({k[0] for k in states})
    m = np.zeros((len(networks), 4), dtype=np.int64)
    pos = {n: i for i, n in enumerate(networks)}
    for k, rows in states.items():
        xs = np.array([float(r[col]) for r in rows])
        ys = np.array([float(r[ENDPOINT]) for r in rows])
        m[pos[k[0]]] += pair_counts(xs, ys)
    return networks, m


def within_bootstrap(states, cols: list[str], n_boot: int, seed: int) -> np.ndarray:
    """(n_boot, len(cols)) resampled within-state tau_b, nan where undefined; one shared draw
    per resample, with the same rng calls as llm_paired_v2._within_state_bootstrap."""
    mats = []
    networks = None
    for col in cols:
        networks, m = network_count_matrix(states, col)
        mats.append(m)
    n_nets = len(networks)
    rng = np.random.default_rng(seed)
    out = np.full((n_boot, len(cols)), np.nan)
    for b in range(n_boot):
        idx = rng.integers(0, n_nets, size=n_nets)
        for c, m in enumerate(mats):
            t = tau_from_counts(m[idx].sum(axis=0))
            if t is not None:
                out[b, c] = t
    return out


def loop_reference_draws(states, col: str, n_draws: int, seed: int) -> list[float | None]:
    """First draws of the bootstrap, recomputed pair by pair with the submitted
    llm_paired_v2._stratified_tau_b (the code behind the published within-state numbers)."""
    keys = list(states)
    networks = sorted({k[0] for k in keys})
    by_net = defaultdict(list)
    for k in keys:
        by_net[k[0]].append(k)
    net_arr = np.array(networks, dtype=object)
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_draws):
        drawn = net_arr[rng.integers(0, len(networks), size=len(networks))]
        resampled = {(rep, k): states[k] for rep, net in enumerate(drawn) for k in by_net[net]}
        vals.append(base._stratified_tau_b(resampled, col, ENDPOINT)["tau_b"])
    return vals


def percentile_ci(vals: np.ndarray) -> tuple[float | None, float | None, int]:
    valid = vals[~np.isnan(vals)]
    if len(valid) < 100:  # the rule of llm_paired_v2._within_state_bootstrap
        return None, None, int(len(vals) - len(valid))
    lo, hi = np.percentile(valid, [2.5, 97.5])
    return float(lo), float(hi), int(len(vals) - len(valid))


def within_cell(states, cols, n_boot, seed, n_check):
    """Within-state tau_b for one column, or the paired lead cols[0] - cols[1]."""
    fails = []
    points = []
    if not states:
        return {"n_units": 0, "n_states": 0, "n_networks": 0, "point": None,
                "ci_lo_2p5": None, "ci_hi_97p5": None, "status": "no state"}, fails
    for col in cols:
        ref = base._stratified_tau_b(states, col, ENDPOINT)["tau_b"]
        _, m = network_count_matrix(states, col)
        mine = tau_from_counts(m.sum(axis=0))
        points.append(mine)
        if (ref is None) != (mine is None) or (ref is not None and abs(ref - mine) > 1e-12):
            fails.append(f"within-state tau_b of {col}: vectorised {mine} != loop {ref}")
    boot = within_bootstrap(states, cols, n_boot, seed)
    for c, col in enumerate(cols):
        ref = loop_reference_draws(states, col, n_check, seed)
        for b, rv in enumerate(ref):
            mv = None if np.isnan(boot[b, c]) else float(boot[b, c])
            if (rv is None) != (mv is None) or (rv is not None and abs(rv - mv) > 1e-12):
                fails.append(f"bootstrap draw {b} of {col}: vectorised {mv} != loop {rv}")
    res = {
        "n_units": sum(len(v) for v in states.values()),
        "n_states": len(states),
        "n_networks": len({k[0] for k in states}),
        "n_states_with_variation": {col: sum(1 for v in states.values() if len({r[col] for r in v}) > 1)
                                    for col in cols},
        "bootstrap_crosscheck_draws": n_check,
    }
    if len(cols) == 1:
        lo, hi, nan = percentile_ci(boot[:, 0])
        res.update(point=points[0], ci_lo_2p5=lo, ci_hi_97p5=hi, n_boot_nan=nan)
    else:
        delta = boot[:, 0] - boot[:, 1]
        lo, hi, nan = percentile_ci(delta)
        valid = delta[~np.isnan(delta)]
        res.update(tau_a=points[0], tau_b=points[1],
                   point=None if None in points else points[0] - points[1],
                   ci_lo_2p5=lo, ci_hi_97p5=hi, n_boot_nan=nan,
                   frac_delta_gt_0=float(np.mean(valid > 0)) if len(valid) else None,
                   n_equal_radius_separation=sum(1 for v in states.values() for r in v
                                                 if r[cols[0]] == r[cols[1]]))
        # the separate CIs, from the same draws, so the reader sees both rankings
        for c, col in enumerate(cols):
            lo_c, hi_c, _ = percentile_ci(boot[:, c])
            res[f"ci_{col}"] = [lo_c, hi_c]
    return res, fails


def pooled_cross_check(rows, col) -> list[str]:
    """tau_b three ways on the same rows: scipy via tau_for_stratum, paired_resample's helper,
    and this file's pair counting as a single stratum."""
    clean = [r for r in rows if r.get(col) is not None and r.get(ENDPOINT) is not None
             and not (col == "radius" and r[col] == UNREACHED)]
    if len(clean) < 2:
        return []
    x = np.array([float(r[col]) for r in clean])
    y = np.array([float(r[ENDPOINT]) for r in clean])
    a = _kendall_tau_b(x, y)
    b = tau_from_counts(pair_counts(x, y)) if len(clean) <= 6000 else a
    fails = []
    if (a is None) != (b is None) or (a is not None and abs(a - b) > 1e-9):
        fails.append(f"pooled tau_b of {col}: kendalltau {a} != pair counts {b}")
    return fails


# ------------------------------------------------------------------ cells

_COHORTS: dict[str, list[dict]] = {}


def _init(args_dict):
    global _COHORTS
    args = argparse.Namespace(**args_dict)
    _COHORTS, _ = build_cohorts(args)


def run_cell(task):
    cohort, cell, n_boot, n_check = task
    t0 = time.perf_counter()
    rows = _COHORTS[cohort]
    sc = scorable(rows)
    undefined = [u for u in sc if u.get("separation") is None]
    defined = [u for u in sc if u.get("separation") is not None]
    fails: list[str] = []
    if cell == "describe_undefined":
        rad = defaultdict(int)
        for u in undefined:
            rad[int(u["radius"])] += 1
        res = {
            "n_scorable": len(sc), "n_separation_undefined": len(undefined),
            "n_separation_defined": len(defined),
            "networks_undefined": sorted({u["network"] for u in undefined}),
            "n_states_undefined": len({(u["network"], u["condition"]) for u in undefined}),
            "radius_counts_undefined": dict(sorted(rad.items())),
            "separation_status_undefined": sorted({u.get("separation_status") for u in undefined}),
            "auc_mean_undefined": float(np.mean([u[ENDPOINT] for u in undefined])) if undefined else None,
            "n_added_units_undefined": sum(1 for u in undefined if "empty_valid" in u["shard_id"]),
        }
        if len(undefined) + len(defined) != len(sc):
            fails.append("undefined + defined != scorable")
    elif cell.startswith("pooled_tau_"):
        col = "separation" if "separation" in cell else "radius"
        subset = {"undefined": undefined, "defined": defined, "all": sc}[cell.rsplit("_", 1)[1]]
        res = tau_for_stratum(subset, col, ENDPOINT, n_boot=n_boot, seed=SEED)
        res["point"] = res["tau_b"]
        fails += pooled_cross_check(subset, col)
    elif cell == "pooled_paired_radius_minus_separation_defined":
        res = paired_comparison(defined, "radius", "separation", ENDPOINT, "network", n_boot=n_boot, seed=SEED)
        res["point"] = res["delta_point"]
        res.pop("loo_by_cluster", None)
    elif cell == "within_tau_radius_undefined":
        states = restrict(radius_varying_states(rows), lambda r: r.get("separation") is None)
        res, fails = within_cell(states, ["radius"], n_boot, SEED, n_check)
    elif cell == "within_tau_radius_all":
        res, fails = within_cell(radius_varying_states(rows), ["radius"], n_boot, SEED, n_check)
    elif cell == "within_paired_radius_minus_separation_matched":
        states = restrict(radius_varying_states(rows), lambda r: r.get("separation") is not None)
        res, fails = within_cell(states, ["radius", "separation"], n_boot, SEED, n_check)
    else:
        raise ValueError(cell)
    n_key = "n_units" if "n_units" in res else ("n" if "n" in res else "n_separation_undefined")
    if cell == "describe_undefined":
        n_key = "n_separation_undefined"
    ref = REFERENCE.get((cohort, cell if cell != "describe_undefined" else "pooled_tau_radius_undefined"))
    repro = {}
    if ref:
        if "n" in ref and res.get(n_key) != ref["n"]:
            fails.append(f"n = {res.get(n_key)}, published {ref['n']}")
        # in the paired within-state cell tau_a / tau_b are the radius and separation taus
        for k in ("point", "tau_a", "tau_b"):
            if k not in ref or cell == "describe_undefined":
                continue
            got = res.get(k)
            if got is None or abs(got - ref[k]) > POINT_TOL:
                fails.append(f"{k} = {got}, published {ref[k]}")
        if "ci" in ref and n_boot == FULL_N_BOOT:
            lo, hi = res.get("ci_lo_2p5"), res.get("ci_hi_97p5")
            repro["ci_published"] = list(ref["ci"])
            repro["ci_ok"] = (lo is not None and hi is not None and abs(lo - ref["ci"][0]) <= CI_TOL
                              and abs(hi - ref["ci"][1]) <= CI_TOL)
    return {"cohort": cohort, "cell": cell, "result": res, "reproduction": repro,
            "consistency_failures": fails, "seconds": round(time.perf_counter() - t0, 2)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--smoke", action="store_true", help="200 resamples, all cells; under a minute")
    ap.add_argument("--out", default=None, help="default results/rebuttal/E16 (smoke: results/rebuttal/E16_smoke)")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--n-boot", type=int, default=None, help="default 10000 (smoke 200)")
    ap.add_argument("--units-csv", default=str(REPO_ROOT / "results/axis_robustness_llm/analysis_units.csv"))
    ap.add_argument("--empty-valid-csv",
                    default=str(REPO_ROOT / "results/axis_robustness_llm_v2/empty_valid_units.csv"))
    ap.add_argument("--gac-dir", default=str(REPO_ROOT / "results/radius_gac_recompute"))
    args = ap.parse_args(argv)
    if not 0 <= args.shard < args.n_shards:
        ap.error("need 0 <= --shard < --n-shards")
    n_boot = args.n_boot or (SMOKE_N_BOOT if args.smoke else FULL_N_BOOT)
    n_check = 3 if args.smoke else 5
    out = Path(args.out or REPO_ROOT / ("results/rebuttal/E16_smoke" if args.smoke else "results/rebuttal/E16"))
    out.mkdir(parents=True, exist_ok=True)
    suffix = "" if args.n_shards == 1 else f".shard{args.shard}of{args.n_shards}"
    rows_path = out / f"rows{suffix}.jsonl"

    commit, src_dirty = git_info()
    t0 = time.perf_counter()
    cohorts, join_fails = build_cohorts(args)
    if join_fails:
        for f in join_fails[:20]:
            print("CONSISTENCY FAILURE:", f, flush=True)
        print(f"{len(join_fails)} join failures; stopping before any statistic", flush=True)
        return 2
    print(f"loaded: submitted {len(scorable(cohorts['submitted']))} scorable, "
          f"corrected {len(scorable(cohorts['corrected']))} scorable; GAC radius equals stored "
          f"on every scorable unit ({time.perf_counter() - t0:.1f}s)", flush=True)

    inputs = {"units_csv": args.units_csv, "empty_valid_csv": args.empty_valid_csv}
    config = {
        "experiment": "E16", "criterion": CRITERION, "endpoint": ENDPOINT, "seed": SEED,
        "n_boot": n_boot, "smoke": args.smoke, "commit": commit, "src_dirty": src_dirty,
        "shard": args.shard, "n_shards": args.n_shards, "cells": list(CELLS), "cohorts": list(COHORTS),
        "sha256": {k: sha256(Path(v)) for k, v in inputs.items()},
        "gac_sha256": {p.name: sha256(p) for p in sorted(Path(args.gac_dir).glob("axis_robustness_llm*.csv"))},
        "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    (out / f"config{suffix}.json").write_text(json.dumps(config, indent=2))

    done = set()
    if rows_path.exists():
        for line in rows_path.read_text().splitlines():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            if not r.get("consistency_failures") and r.get("n_boot") == n_boot:
                done.add(r["key"])
    all_tasks = [(c, cell) for c in COHORTS for cell in CELLS]
    tasks = [(c, cell, n_boot, n_check) for i, (c, cell) in enumerate(all_tasks)
             if i % args.n_shards == args.shard and f"{c}/{cell}" not in done]
    print(f"{len(tasks)} cells to run in shard {args.shard}/{args.n_shards} "
          f"({len(done)} already written), n_boot={n_boot}, workers={args.workers}", flush=True)

    state = {"done": 0}
    stop = threading.Event()

    def heartbeat():
        while not stop.wait(30):
            print(f"[{time.perf_counter() - t0:7.0f}s] {state['done']}/{len(tasks)} cells done", flush=True)

    threading.Thread(target=heartbeat, daemon=True).start()
    n_fail = 0
    init_args = {k: getattr(args, k) for k in ("units_csv", "empty_valid_csv", "gac_dir")}
    with rows_path.open("a", encoding="utf-8") as fh:
        if args.workers > 1:
            pool = mp.get_context("spawn").Pool(min(args.workers, 8), initializer=_init, initargs=(init_args,))
            results = pool.imap_unordered(run_cell, tasks)
        else:
            global _COHORTS
            _COHORTS = cohorts
            pool, results = None, map(run_cell, tasks)
        for r in results:
            r.update(key=f"{r['cohort']}/{r['cell']}", experiment="E16", criterion=CRITERION,
                     commit=commit, src_dirty=src_dirty, n_boot=n_boot, seed=SEED, smoke=args.smoke)
            fh.write(json.dumps(r, default=float) + "\n")
            fh.flush()
            state["done"] += 1
            res = r["result"]
            print(f"[{time.perf_counter() - t0:7.0f}s] {r['key']}: point={res.get('point')} "
                  f"CI=[{res.get('ci_lo_2p5')}, {res.get('ci_hi_97p5')}] ({r['seconds']}s)", flush=True)
            if r["reproduction"].get("ci_ok") is False:
                print(f"  WARNING interval differs from published {r['reproduction']['ci_published']}", flush=True)
            for f in r["consistency_failures"]:
                n_fail += 1
                print(f"  CONSISTENCY FAILURE: {f}", flush=True)
        if pool is not None:
            pool.close()
            pool.join()
    stop.set()
    print(f"finished {state['done']} cells in {time.perf_counter() - t0:.0f}s, "
          f"{n_fail} consistency failures", flush=True)
    return 2 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
