"""E18: Section 7 mean-bias profiles under standardised models, ten coefficient draws.

For each of the 49 queries of Table 6, three arms (naive, posthoc, internal) times ten
coefficient draws. One JSON line per (query, draw, arm). See PREREG.md.

The graph-only part of a query (the slice of states agreeing with G0 outside X's
component, their distances, and the possible parent sets of X in each) does not depend
on the SEM, so it is computed once per query and reused by the 30 units; only the
regressions are redone per unit. That is the "cached enumeration" of the plan.

Usage, from the repository root:
    PYTHONPATH=src .venv/bin/python experiments/rebuttal/E18_section7_standardised_models/run.py --smoke
    PYTHONPATH=src .venv/bin/python experiments/rebuttal/E18_section7_standardised_models/run.py --workers 6
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
for p in (ROOT / "src", ROOT / "experiments", HERE):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag, random_sem  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.epsilon.bias import (  # noqa: E402
    ZERO_TOL,
    make_context,
    make_context_from_covariance,
    possible_parent_sets_enumerated,
)
from bkrobust.hybrid import breakdown_radius  # noqa: E402
from bkrobust.search.space_fixed import knowledge_of  # noqa: E402

import final_table as ft  # noqa: E402
import mean_table_relevant_v2 as rel  # noqa: E402
import mean_table_xcomp_extend_v2 as xc  # noqa: E402
import vendor_scm  # noqa: E402

CRITERION = "gac"
ARMS = ("naive", "posthoc", "internal")
EPSILONS = [0.02, 0.10, 0.25, 0.40, 0.60, 0.90]
BASE_SEED = 20260917  # the seed of the stored Section 7 run
WEIGHT_RANGE = (0.3, 1.5)  # random_sem defaults, recorded so a changed default is visible
NOISE_RANGE = (0.5, 1.5)
SMOKE_QUERIES = [
    ("Acid_1996", "x1", "x10"),
    ("asia", "bronc", "dysp"),
    ("barley", "dg25", "s2225"),
    ("Schipf_2010", "A", "U"),
    ("magic-niab", "G1217", "FT"),  # an extended row: exercises the stored b_top check
]
# Tolerances for "posthoc equals naive". Relative 1e-9, plus 1e-8 on the relative scale
# because bias_at clamps |bias| < 1e-9 to zero in absolute units, and rescaling moves
# what sits near that clamp.
REL_TOL, ABS_TOL = 1e-9, 1e-8


class Timeout(Exception):
    pass


def _alarm(_signum: int, _frame: Any) -> None:
    raise Timeout()


# --- population -------------------------------------------------------------


def load_population(data_dir: Path) -> list[dict[str, Any]]:
    """The 49 stored rows of Table 6: 39 x_component rows, 10 extended rows."""
    rows = []
    for line in (data_dir / "x_component" / "instances.jsonl").open():
        r = json.loads(line)
        if r.get("status") == "ok":
            r["_source"] = "x_component"
            rows.append(r)
    for line in (data_dir / "extended" / "instances.jsonl").open():
        r = json.loads(line)
        if r.get("status") == "ok" and r.get("has_finite_r_val") and r.get("has_nonzero_bias"):
            r["_source"] = "extended"
            rows.append(r)
    return sorted(rows, key=lambda r: (r["network"], r["x"], r["y"]))


def draw_seed(net: str, x: str, y: str, draw: int) -> int:
    # Draw 0 is the stored seed exactly; later draws salt the network name so that no
    # draw can coincide with another query's stored seed.
    return ft.derive_seed(BASE_SEED, net if draw == 0 else f"{net}|draw{draw}", x, y)


def git_commit() -> str:
    try:
        out = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=10)
        dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--", "src"],
                               capture_output=True, text=True, timeout=10).stdout.strip()
        return out.stdout.strip() + ("+dirty-src" if dirty else "")
    except Exception:
        return "unknown"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# --- worker state ------------------------------------------------------------

_G: dict[str, Any] = {}


def init_worker(networks: list[str], condition: str) -> None:
    signal.signal(signal.SIGALRM, _alarm)
    _G["knowledge"] = ft.load_knowledge(condition)
    _G["parsed"] = ft.load_networks(names=set(networks), skip=set(), max_nodes=0)


# --- the SEM-free part of a query -------------------------------------------


def prepare_query(net: str, x: str, y: str) -> dict[str, Any]:
    cpdag, dag = _G["parsed"][net]["cpdag"], _G["parsed"][net]["dag"]
    g0 = apply_orientations(cpdag, list(_G["knowledge"][net]))
    if g0 is None:
        raise RuntimeError("g0 build failed")
    z = optimal_adjustment_set_mpdag(g0, x, y)
    if z is None:
        raise RuntimeError("degenerate query")
    z = frozenset(z)

    # Same construction as xc.build_record, kept step for step so draw 0 is bit-identical.
    x_nodes, x_edges = xc.x_component_edges(cpdag, x)
    outside_fixed = xc.outside_fixed_orientations(cpdag, g0, x_nodes)
    slice_elems = ({g0.edge_string(): g0} if not x_edges
                   else xc.enumerate_slice(cpdag, x_edges, outside_fixed))
    outside_resolution = xc.resolve_outside_via_dag(cpdag, dag, x_nodes, outside_fixed)
    elems = list(slice_elems.values())
    dists = {g: xc.slice_distance(cpdag, g, g0) for g in elems}
    psets = {}
    for g in elems:
        g_eval = apply_orientations(cpdag, outside_resolution + list(knowledge_of(cpdag, g)))
        assert g_eval is not None, "resolved graph must stay consistent"
        psets[g] = possible_parent_sets_enumerated(g_eval, x)
    top = apply_orientations(cpdag, xc.resolve_outside_via_dag(cpdag, dag, x_nodes, []))
    ps_top = possible_parent_sets_enumerated(top, x)

    rad = breakdown_radius(cpdag, None, x, y, z, g0=g0, criterion=CRITERION)
    return {"cpdag": cpdag, "dag": dag, "g0": g0, "z": z, "x_nodes": x_nodes,
            "k_x": len(x_edges), "elems": elems, "dists": dists, "psets": psets,
            "ps_top": ps_top, "radius": rad}


def worst_rel(ctx: Any, parent_sets: frozenset, taus: dict, scale: float) -> float:
    """bias_at(...).worst / scale, from cached parent sets (same arithmetic as bias_at)."""
    if not parent_sets:
        raise ValueError("no DAG extension; B is undefined")
    vals = []
    for p in parent_sets:
        if p not in taus:
            taus[p] = ctx.total_effect_under(p)
        vals.append(abs(ctx.theta_z - taus[p]))
    w = max(vals)
    return (0.0 if w < ZERO_TOL else float(w)) / scale


def profile_fields(ctx: Any, q: dict[str, Any]) -> dict[str, Any]:
    scale = abs(ctx.theta_z)
    taus: dict = {}
    b = {g: worst_rel(ctx, q["psets"][g], taus, scale) for g in q["elems"]}
    prof = rel.profile_from(b, q["elems"], q["dists"])
    mus = [pp["mu"] for pp in prof]
    return {
        "max_distance_rel": max(pp["d"] for pp in prof) if prof else 0,
        "n_reachable_rel": len(q["elems"]),
        "mu_peak_rel": max(mus) if mus else 0.0,
        "mu_peak_rel_at": prof[mus.index(max(mus))]["d"] if mus else None,
        "mu_monotone_rel": all(mus[i] >= mus[i - 1] - 1e-12 for i in range(1, len(mus))),
        "r_mu_rel": {str(e): rel.r_onset(prof, e) for e in EPSILONS},
        "r_mu_rel_persists": {str(e): rel.persists(prof, e) for e in EPSILONS},
        "b_top": worst_rel(ctx, q["ps_top"], taus, scale),
        "profile_rel": prof,
    }


def arm_contexts(q: dict[str, Any], x: str, y: str, seed: int) -> dict[str, tuple[Any, Any]]:
    """(BiasContext, the SEM's covariance and nodes) per arm, all from one coefficient draw."""
    dag, cpdag, z = q["dag"], q["cpdag"], q["z"]
    sem = random_sem(dag, np.random.default_rng(seed), WEIGHT_RANGE, NOISE_RANGE)
    sigma = sem.covariance()
    nodes = tuple(dag.nodes)
    out = {"naive": (make_context(sem, cpdag, x, y, z), (sigma, nodes))}

    sd = np.sqrt(np.diag(sigma))
    sigma_std = sigma / np.outer(sd, sd)
    i = {n: k for k, n in enumerate(nodes)}
    tau_std = sem.true_total_effect(x, y) * sd[i[x]] / sd[i[y]]
    out["posthoc"] = (make_context_from_covariance(sigma_std, nodes, cpdag, x, y, z,
                                                   tau_true=tau_std), (sigma_std, nodes))

    sem_int, _ = vendor_scm.iscm_population(sem)
    out["internal"] = (make_context(sem_int, cpdag, x, y, z), (sem_int.covariance(), nodes))
    return out


def close(a: float, b: float) -> bool:
    return abs(a - b) <= REL_TOL * max(1.0, abs(a), abs(b)) + ABS_TOL


def same_as_naive(mine: dict[str, Any], naive: dict[str, Any]) -> list[str]:
    bad = []
    if mine["r_mu_rel"] != naive["r_mu_rel"]:
        bad.append(f"r_mu {mine['r_mu_rel']} != {naive['r_mu_rel']}")
    if mine["r_mu_rel_persists"] != naive["r_mu_rel_persists"]:
        bad.append("r_mu_persists differ")
    if not close(mine["b_top"], naive["b_top"]):
        bad.append(f"b_top {mine['b_top']} != {naive['b_top']}")
    for a, b in zip(mine["profile_rel"], naive["profile_rel"], strict=True):
        if a["d"] != b["d"] or not close(a["mu"], b["mu"]) or not close(a["sphere_max"], b["sphere_max"]):
            bad.append(f"profile differs at d={a['d']}")
            break
    return bad


# --- one query, all its units ------------------------------------------------


def run_query(task: dict[str, Any]) -> dict[str, Any]:
    stored, todo, opts = task["stored"], set(map(tuple, task["todo"])), task["opts"]
    net, x, y = stored["network"], stored["x"], stored["y"]
    base = {"experiment": "E18", "criterion": CRITERION, "commit": opts["commit"],
            "network": net, "x": x, "y": y, "r_val_stored": stored["r_val"],
            "source": stored["_source"]}
    rows, failures = [], []
    t0 = time.perf_counter()
    signal.alarm(opts["budget"])
    try:
        q = prepare_query(net, x, y)
        signal.alarm(0)
    except Timeout:
        for d, arm in sorted(todo):
            rows.append({**base, "draw": d, "arm": arm, "status": "timeout",
                         "budget": opts["budget"]})
        return {"rows": rows, "failures": [], "key": (net, x, y)}
    except Exception as exc:  # every stored query built before; failing now is a consistency failure
        for d, arm in sorted(todo):
            rows.append({**base, "draw": d, "arm": arm, "status": "error", "error": repr(exc)})
        return {"rows": rows, "failures": [f"{net} {x}->{y}: {exc!r}"], "key": (net, x, y)}
    finally:
        signal.alarm(0)
    prep_s = time.perf_counter() - t0
    rad = q["radius"]
    base.update({"radius_gac": rad.radius, "radius_exact": rad.exact, "radius_method": rad.method,
                 "z": sorted(q["z"]), "k_x": q["k_x"], "slice_size": len(q["elems"])})
    if rad.radius != stored["r_val"] or not rad.exact:
        failures.append(f"{net} {x}->{y}: GAC radius {rad.radius} (exact={rad.exact}) "
                        f"!= stored r_val {stored['r_val']}")

    printed = opts["printed"].get(f"{net}|{x}|{y}")
    for d in range(opts["draws"]):
        seed = draw_seed(net, x, y, d)
        ctxs = arm_contexts(q, x, y, seed)
        rec: dict[str, dict[str, Any]] = {}
        for arm in ARMS:  # naive first: the posthoc check needs it, even when naive is done
            t1 = time.perf_counter()
            ctx, (sigma, nodes) = ctxs[arm]
            row = {**base, "draw": d, "arm": arm, "seed": seed}
            if abs(ctx.theta_z) < 1e-12:
                row["status"] = "zero_estimate"
                rec[arm] = row
            else:
                row.update({"status": "ok", "theta_z": ctx.theta_z, "tau_true": ctx.tau_true,
                            "scale_std": ctx.scale_std})
                row.update(profile_fields(ctx, q))
                row["sortability"] = vendor_scm.sortability(sigma, nodes, q["dag"])
                checks: dict[str, Any] = {}
                if arm == "naive" and d == 0:
                    probs = xc.compare_records(row, stored, EPSILONS)
                    if stored.get("b_top") is not None and not close(row["b_top"], stored["b_top"]):
                        probs.append(f"b_top {row['b_top']} != stored {stored['b_top']}")
                    if printed is not None and abs(row["b_top"] - printed) > 5e-4 + 1e-9:
                        probs.append(f"b_top {row['b_top']:.6f} != printed Table 6 {printed}")
                    checks["repro_stored"] = not probs
                    checks["repro_problems"] = probs
                    if len(q["elems"]) <= opts["direct_check_max"]:
                        # The original code path, uncached, on the same context.
                        direct = xc.build_record(net, x, y, stored["r_val"], None, q["cpdag"],
                                                 q["dag"], q["g0"], ctx, abs(ctx.theta_z),
                                                 EPSILONS)
                        dprobs = xc.compare_records(row, direct, EPSILONS)
                        if not close(row["b_top"], direct["b_top"]):
                            dprobs.append("b_top differs from build_record")
                        checks["direct_match"] = not dprobs
                        probs = probs + [f"direct: {p}" for p in dprobs]
                    for p in probs:
                        failures.append(f"{net} {x}->{y} draw 0 naive: {p}")
                if arm == "posthoc" and rec["naive"].get("status") == "ok":
                    diff = same_as_naive(row, rec["naive"])
                    checks["posthoc_match"] = not diff
                    for p in diff:
                        failures.append(f"{net} {x}->{y} draw {d} posthoc vs naive: {p}")
                if arm == "internal":
                    ok = bool(np.allclose(np.diag(sigma), 1.0, rtol=0, atol=1e-9))
                    checks["unit_variance"] = ok
                    if not ok:
                        failures.append(f"{net} {x}->{y} draw {d} internal: variances not 1")
                row["checks"] = checks
                row["seconds"] = round(time.perf_counter() - t1, 4)
                rec[arm] = row
            if (d, arm) in todo:
                row["prepare_seconds"] = round(prep_s, 2)
                rows.append(row)
    return {"rows": rows, "failures": failures, "key": (net, x, y)}


# --- driver --------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smoke", action="store_true", help="5 queries, 2 draws, all arms (< 60 s)")
    ap.add_argument("--out", default=None, help="default results/rebuttal/E18 (smoke: .../E18/smoke)")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--draws", type=int, default=10)
    ap.add_argument("--data-dir", default="results/mean_table_relevant_v2",
                    help="stored Section 7 run (x_component/ and extended/), relative to the repo root")
    ap.add_argument("--condition", default="D_LLM")
    ap.add_argument("--budget", type=int, default=3600, help="seconds for one query's enumeration")
    ap.add_argument("--direct-check-max", type=int, default=200,
                    help="re-run the uncached build_record on draw 0 when the slice has at most this many states")
    args = ap.parse_args()
    if args.smoke:
        args.draws = min(args.draws, 2)
    if not 0 <= args.shard < args.n_shards:
        ap.error("need 0 <= --shard < --n-shards")
    out_dir = Path(args.out or ("results/rebuttal/E18/smoke" if args.smoke else "results/rebuttal/E18"))
    out_dir = out_dir if out_dir.is_absolute() else ROOT / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    data_dir = Path(args.data_dir)
    data_dir = data_dir if data_dir.is_absolute() else ROOT / data_dir

    population = load_population(data_dir)
    if len(population) != 49:
        print(f"[E18] population has {len(population)} rows, expected 49", flush=True)
        return 3
    if args.smoke:
        population = [r for r in population if (r["network"], r["x"], r["y"]) in set(SMOKE_QUERIES)]
    mine = [r for i, r in enumerate(population) if i % args.n_shards == args.shard]

    printed = {f"{r['network']}|{r['x']}|{r['y']}": r["b_rel_chat"]
               for r in json.loads((HERE / "table6_printed.json").read_text())}
    tag = f"s{args.shard}of{args.n_shards}"
    units_path = out_dir / f"units.{tag}.jsonl"
    fail_path = out_dir / f"failures.{tag}.jsonl"
    commit = git_commit()
    config = {
        "experiment": "E18", "criterion": CRITERION, "smoke": args.smoke, "arms": list(ARMS),
        "draws": args.draws, "epsilons": EPSILONS, "seed": BASE_SEED,
        "weight_range": list(WEIGHT_RANGE), "noise_range": list(NOISE_RANGE),
        "condition": args.condition, "n_population": 49, "n_selected": len(population),
        "n_this_shard": len(mine), "shard": args.shard, "n_shards": args.n_shards,
        "budget": args.budget, "direct_check_max": args.direct_check_max, "commit": commit,
        "x_component_sha256": sha256(data_dir / "x_component" / "instances.jsonl"),
        "extended_sha256": sha256(data_dir / "extended" / "instances.jsonl"),
        "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    (out_dir / f"config.{tag}.json").write_text(json.dumps(config, indent=2) + "\n")

    done: set[tuple] = set()
    if units_path.exists():
        for line in units_path.open():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue  # a line cut by a kill; its unit is recomputed
            done.add((r["network"], r["x"], r["y"], r["draw"], r["arm"]))
    opts = {"commit": commit, "budget": args.budget, "draws": args.draws,
            "direct_check_max": args.direct_check_max, "printed": printed}
    tasks = []
    for r in mine:
        todo = [(d, a) for d in range(args.draws) for a in ARMS
                if (r["network"], r["x"], r["y"], d, a) not in done]
        if todo:
            tasks.append({"stored": r, "todo": todo, "opts": opts})
    n_units = len(mine) * args.draws * len(ARMS)
    print(f"[E18] {len(mine)} queries in shard {tag}, {n_units} units, {n_units - sum(len(t['todo']) for t in tasks)} "
          f"already written, {len(tasks)} queries to run, workers={args.workers}, commit={commit}", flush=True)

    state = {"queries": 0, "units": 0, "fail": 0, "t0": time.time()}
    stop = threading.Event()

    def heartbeat() -> None:
        while not stop.wait(60):
            print(f"[E18] progress: {state['queries']}/{len(tasks)} queries, {state['units']} units written, "
                  f"{state['fail']} failures, {time.time() - state['t0']:.0f}s", flush=True)

    threading.Thread(target=heartbeat, daemon=True).start()
    networks = sorted({t["stored"]["network"] for t in tasks})
    all_failures: list[str] = []
    with units_path.open("a") as fu, fail_path.open("a") as ff:
        def take(res: dict[str, Any]) -> None:
            for row in res["rows"]:
                fu.write(json.dumps(row) + "\n")
            fu.flush()
            for f in res["failures"]:
                ff.write(json.dumps({"failure": f, "commit": commit}) + "\n")
                print(f"[E18] CONSISTENCY FAILURE: {f}", flush=True)
            ff.flush()
            all_failures.extend(res["failures"])
            state["queries"] += 1
            state["units"] += len(res["rows"])
            state["fail"] += len(res["failures"])
            net, x, y = res["key"]
            print(f"[E18] done {net} {x}->{y}: {len(res['rows'])} units "
                  f"({state['queries']}/{len(tasks)}, {time.time() - state['t0']:.0f}s)", flush=True)

        if tasks and args.workers <= 1:
            init_worker(networks, args.condition)
            for t in tasks:
                take(run_query(t))
        elif tasks:
            # spawn, not fork: the heartbeat thread is already running in this process.
            ctx = mp.get_context("spawn")
            # Largest slices first, so the slow queries do not start last.
            tasks.sort(key=lambda t: -t["stored"].get("n_reachable_rel", 0))
            with ctx.Pool(args.workers, initializer=init_worker,
                          initargs=(networks, args.condition)) as pool:
                for res in pool.imap_unordered(run_query, tasks):
                    take(res)
    stop.set()
    print(f"[E18] finished: {state['units']} units written, {len(all_failures)} consistency failures, "
          f"{time.time() - state['t0']:.0f}s -> {units_path}", flush=True)
    return 2 if all_failures else 0


if __name__ == "__main__":
    sys.exit(main())
