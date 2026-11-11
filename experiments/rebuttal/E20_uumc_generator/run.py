"""E20: Section 7 mean-bias profiles under naive, iSCM and UUMC coefficients, ten draws each.

One unit = (query, arm, draw); one JSON line per unit. The graph-only work (the slice of states
agreeing with G0 outside X's component, their distances and their possible parent sets of X) is
done once per query and reused by every draw, because B(G) depends on the coefficients only
through the regressions tau(pa) and theta_Z.

    PYTHONPATH=src .venv/bin/python experiments/rebuttal/E20_uumc_generator/run.py --smoke
    PYTHONPATH=src .venv/bin/python experiments/rebuttal/E20_uumc_generator/run.py --workers 8

Exit 0: all consistency checks passed. Exit 2: at least one failed (rows are still written).
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))
sys.path.insert(0, str(HERE))

from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.epsilon.bias import (  # noqa: E402
    ZERO_TOL, bias_at, make_context, make_context_from_covariance,
    possible_parent_sets_enumerated, possible_parent_sets_semilocal)
from bkrobust.hybrid import breakdown_radius  # noqa: E402
from bkrobust.search.space_fixed import knowledge_of  # noqa: E402

import final_table as ft  # noqa: E402
import mean_table_relevant_v2 as rel  # noqa: E402
import mean_table_xcomp_extend_v2 as xc  # noqa: E402

import generators as gen  # noqa: E402

EPSILONS = [0.02, 0.10, 0.25, 0.40, 0.60, 0.90]
ARMS = ("naive", "iscm", "uumc")
CRITERION = "gac"
# Small slices covering both sources (x_component, extended) and two of the Table 2 queries.
SMOKE_QUERIES = [("asia", "asia", "dysp"), ("magic-niab", "G1217", "FT"),
                 ("Schipf_2010", "A", "U"), ("barley", "dg25", "s2225")]
TOL = 1e-9


class Timeout(Exception):
    pass


def _alarm(_s, _f):
    raise Timeout()


def repo_commit() -> str:
    try:
        return subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True,
                              text=True, timeout=10).stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def src_dirty() -> bool | None:
    """Whether src/ differs from the commit; a run on a patched library is not this experiment."""
    try:
        out = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "src"],
                             capture_output=True, text=True, timeout=10)
        return bool(out.stdout.strip()) if out.returncode == 0 else None
    except Exception:
        return None


def load_queries(data_dir: Path) -> list[dict[str, Any]]:
    """The 49 rows of Table 6, with the committed values draw 0 must reproduce."""
    fullspace = {(r["network"], r["x"], r["y"]): r for r in map(
        json.loads, (data_dir / "mean_table_fullspace/instances.jsonl").open())}
    out = []
    for line in (data_dir / "mean_table_relevant_v2/x_component/instances.jsonl").open():
        r = json.loads(line)
        if r.get("status") == "ok":
            r["b_top_stored"] = fullspace[(r["network"], r["x"], r["y"])]["beta_top"]
            r["source"] = "x_component"
            out.append(r)
    for line in (data_dir / "mean_table_relevant_v2/extended/instances.jsonl").open():
        r = json.loads(line)
        if r.get("status") == "ok" and r.get("has_finite_r_val") and r.get("has_nonzero_bias"):
            r["b_top_stored"] = r["b_top"]
            r["source"] = "extended"
            out.append(r)
    return sorted(out, key=lambda r: (r["network"], r["x"], r["y"]))


def draw_sem(arm: str, draw: int, seed: int, dag, net: str, x: str, y: str):
    base = ft.derive_seed(seed + draw, net, x, y)
    if arm == "naive":
        return gen.naive_sem(dag, np.random.default_rng(base)), base
    if arm == "iscm":
        # Paired with naive draw k, so the only difference is the standardisation.
        return gen.iscm_from(gen.naive_sem(dag, np.random.default_rng(base))), base
    if arm == "uumc":
        return gen.uumc_sem(dag, np.random.default_rng(np.random.SeedSequence([base, 1]))), base
    raise ValueError(arm)


def b_from_parent_sets(theta: float, tau: dict, psets) -> float:
    """Same arithmetic as bias_at: max |theta - tau| over the parent sets, zeroed below ZERO_TOL."""
    worst = max(abs(theta - tau[p]) for p in psets)
    return 0.0 if worst < ZERO_TOL else float(worst)


def build_cache(q: dict, cpdag, dag, g0, budget: int) -> dict[str, Any]:
    """Graph-only part, once per query; mirrors mean_table_xcomp_extend_v2.build_record."""
    x = q["x"]
    x_nodes, x_edges = xc.x_component_edges(cpdag, x)
    outside_fixed = xc.outside_fixed_orientations(cpdag, g0, x_nodes)
    signal.alarm(budget)
    try:
        if not x_edges:
            slice_elems = {g0.edge_string(): g0}
        else:
            slice_elems = xc.enumerate_slice(cpdag, x_edges, outside_fixed)
        outside_res = xc.resolve_outside_via_dag(cpdag, dag, x_nodes, outside_fixed)
        states = []
        semilocal_mismatch = 0
        for key in sorted(slice_elems):
            g = slice_elems[key]
            g_eval = apply_orientations(cpdag, outside_res + list(knowledge_of(cpdag, g)))
            assert g_eval is not None
            ps = possible_parent_sets_enumerated(g_eval, x)
            if ps != possible_parent_sets_semilocal(cpdag, g_eval, x):
                semilocal_mismatch += 1
            states.append({"g": g, "g_eval": g_eval, "d": xc.slice_distance(cpdag, g, g0),
                           "psets": sorted(ps, key=sorted)})
        top = apply_orientations(cpdag, xc.resolve_outside_via_dag(cpdag, dag, x_nodes, []))
        top_ps = sorted(possible_parent_sets_enumerated(top, x), key=sorted)
        signal.alarm(0)
    except Timeout:
        signal.alarm(0)
        return {"timeout": True}
    anc, n_lengths = gen.path_structure(dag)
    return {"timeout": False, "k_x": len(x_edges), "states": states, "top": top, "top_ps": top_ps,
            "semilocal_mismatch": semilocal_mismatch, "anc": anc, "n_lengths": n_lengths}


def rel_profile(ctx, cache) -> tuple[list[dict], float, list[float]]:
    """mu profile, B_rel(Chat) and per-state B_rel under one context."""
    scale = abs(ctx.theta_z)
    tau: dict = {}
    for st in cache["states"]:
        for p in st["psets"]:
            if p not in tau:
                tau[p] = ctx.total_effect_under(p)
    for p in cache["top_ps"]:
        if p not in tau:
            tau[p] = ctx.total_effect_under(p)
    bs = [b_from_parent_sets(ctx.theta_z, tau, st["psets"]) / scale for st in cache["states"]]
    b_top = b_from_parent_sets(ctx.theta_z, tau, cache["top_ps"]) / scale
    b_map = {st["g"]: b for st, b in zip(cache["states"], bs, strict=True)}
    dists = {st["g"]: st["d"] for st in cache["states"]}
    prof = rel.profile_from(b_map, [st["g"] for st in cache["states"]], dists)
    return prof, b_top, bs


def compare_stored(prof, b_top, theta, q) -> list[str]:
    bad = []
    if abs(theta - q["theta_z"]) > TOL:
        bad.append(f"theta_z {theta} != stored {q['theta_z']}")
    if abs(b_top - q["b_top_stored"]) > TOL:
        bad.append(f"B(Chat) {b_top} != stored {q['b_top_stored']}")
    mine = {p["d"]: p for p in prof}
    theirs = {p["d"]: p for p in q["profile_rel"]}
    if set(mine) != set(theirs):
        bad.append("profile distances differ")
    else:
        for d in mine:
            for k in ("mu", "beta", "sphere_max"):
                if abs(mine[d][k] - theirs[d][k]) > TOL:
                    bad.append(f"d={d} {k} {mine[d][k]} != {theirs[d][k]}")
            if mine[d]["n_sphere"] != theirs[d]["n_sphere"]:
                bad.append(f"d={d} n_sphere differs")
    for e in EPSILONS:
        if rel.r_onset(prof, e) != q["r_mu_rel"][str(e)]:
            bad.append(f"r_mu({e}) {rel.r_onset(prof, e)} != {q['r_mu_rel'][str(e)]}")
    return bad


def process_query(job: dict) -> list[dict]:
    """All missing units of one query. Runs in a worker process."""
    signal.signal(signal.SIGALRM, _alarm)
    q, todo, args = job["q"], job["todo"], job["args"]
    net, x, y = q["network"], q["x"], q["y"]
    base = {"network": net, "x": x, "y": y, "criterion": CRITERION, "commit": job["commit"],
            "r_val": q["r_val"], "source": q["source"]}
    t0 = time.perf_counter()
    knowledge = ft.load_knowledge(args["condition"])
    parsed = ft.load_networks(names={net}, skip=set(), max_nodes=0)
    if net not in parsed or net not in knowledge:
        return [{**base, "arm": a, "draw": k, "status": "unavailable"} for a, k in todo]
    cpdag, dag = parsed[net]["cpdag"], parsed[net]["dag"]
    g0 = apply_orientations(cpdag, list(knowledge[net]))
    z = frozenset(optimal_adjustment_set_mpdag(g0, x, y) or ())
    graph_fail: list[str] = []
    if sorted(z) != sorted(q["z"]):
        graph_fail.append(f"Z {sorted(z)} != stored {q['z']}")

    # The radius is graph-only; recompute it under the GAC as a second route to the stored value.
    r_val_gac = None
    signal.alarm(args["radius_budget"])
    try:
        r_val_gac = breakdown_radius(cpdag, None, x, y, z, g0=g0, criterion=CRITERION).radius
        signal.alarm(0)
    except Timeout:
        signal.alarm(0)
    if r_val_gac is not None and str(r_val_gac) != str(q["r_val"]):
        graph_fail.append(f"GAC radius {r_val_gac} != stored r_val {q['r_val']}")
    r_check = int(q["r_val"])

    cache = build_cache(q, cpdag, dag, g0, args["budget"])
    if cache["timeout"]:
        return [{**base, "arm": a, "draw": k, "status": "timeout", "budget": args["budget"]}
                for a, k in todo]
    if cache["semilocal_mismatch"]:
        graph_fail.append(f"{cache['semilocal_mismatch']} states: enumerated != semilocal parent sets")
    cache_s = time.perf_counter() - t0
    # States checked against the library's own bias_at, chosen by position so the choice is fixed.
    n = len(cache["states"])
    probe = sorted({0, n // 2, n - 1})

    rows = []
    for arm, k in todo:
        fails = list(graph_fail)
        sem, base_seed = draw_sem(arm, k, args["seed"], dag, net, x, y)
        sigma = sem.covariance()
        if arm in ("iscm", "uumc") and np.max(np.abs(np.diag(sigma) - 1.0)) > 1e-9:
            fails.append("standardised arm has a variance != 1")
        ctx = make_context(sem, cpdag, x, y, z)
        row = {**base, "arm": arm, "draw": k, "seed": base_seed, "r_val_gac": r_val_gac,
               "k_x": cache["k_x"], "slice_size": n, "theta_z": ctx.theta_z,
               "tau_true": ctx.tau_true, **gen.sortabilities(sigma, cache["anc"], cache["n_lengths"])}
        if abs(ctx.theta_z) < 1e-12:
            rows.append({**row, "status": "zero_estimate", "consistency_failures": fails})
            continue
        prof, b_top, bs = rel_profile(ctx, cache)
        for i in probe:
            direct = bias_at(ctx, cache["states"][i]["g_eval"], method="enumerate").worst / abs(ctx.theta_z)
            if abs(direct - bs[i]) > TOL:
                fails.append(f"cached B {bs[i]} != bias_at {direct} at state {i}")
        below = [b for st, b in zip(cache["states"], bs, strict=True) if st["d"] < r_check]
        if any(b > 0 for b in below):
            fails.append("nonzero B below r_val (Theorem 1)")
        if arm == "naive" and k == 0:
            fails += compare_stored(prof, b_top, ctx.theta_z, q)
            # Post hoc standardisation rescales theta and every tau by sd_x/sd_y: B_rel must not move.
            sd = np.sqrt(np.diag(sigma))
            ctx_std = make_context_from_covariance(sigma / np.outer(sd, sd), ctx.nodes, cpdag, x, y, z)
            prof_std, b_top_std, _ = rel_profile(ctx_std, cache)
            if abs(b_top_std - b_top) > TOL or any(abs(a["mu"] - b["mu"]) > TOL
                                                    for a, b in zip(prof, prof_std, strict=True)):
                fails.append("post hoc standardisation changed the relative profile")
        mus = [p["mu"] for p in prof]
        row.update({
            "status": "ok",
            "D": max(p["d"] for p in prof),
            "b_top_rel": b_top,
            "mu_peak_rel": max(mus),
            "mu_monotone_rel": all(mus[i] >= mus[i - 1] - 1e-12 for i in range(1, len(mus))),
            "r_mu_rel": {str(e): rel.r_onset(prof, e) for e in EPSILONS},
            "r_mu_rel_persists": {str(e): rel.persists(prof, e) for e in EPSILONS},
            "profile_rel": [{"d": p["d"], "n_sphere": p["n_sphere"], "mu": p["mu"],
                             "sphere_max": p["sphere_max"]} for p in prof],
            "consistency_failures": fails,
            "cache_seconds": round(cache_s, 2),
        })
        rows.append(row)
    return rows


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--smoke", action="store_true", help="4 small queries, draws 0-1, all arms")
    p.add_argument("--out", default=None, help="default results/rebuttal/E20 (E20_smoke with --smoke)")
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--shard", type=int, default=0)
    p.add_argument("--n-shards", type=int, default=1)
    p.add_argument("--draws", type=int, default=10)
    p.add_argument("--arms", default=",".join(ARMS))
    p.add_argument("--seed", type=int, default=20260917)
    p.add_argument("--condition", default="D_LLM")
    p.add_argument("--budget", type=int, default=3600, help="seconds for one query's state cache")
    p.add_argument("--radius-budget", type=int, default=300, help="seconds for the GAC radius")
    p.add_argument("--data-dir", default="results",
                   help="folder holding mean_table_fullspace/ and mean_table_relevant_v2/ (repo-relative)")
    args = p.parse_args()

    arms = [a for a in args.arms.split(",") if a]
    if any(a not in ARMS for a in arms):
        p.error(f"--arms must be among {ARMS}")
    draws = 2 if args.smoke else args.draws
    out = ROOT / (args.out or ("results/rebuttal/E20_smoke" if args.smoke else "results/rebuttal/E20"))
    out.mkdir(parents=True, exist_ok=True)
    data_dir = Path(args.data_dir)
    data_dir = data_dir if data_dir.is_absolute() else ROOT / data_dir

    queries = load_queries(data_dir)
    if args.smoke:
        want = set(SMOKE_QUERIES)
        queries = [q for q in queries if (q["network"], q["x"], q["y"]) in want]
    queries = [q for i, q in enumerate(queries) if i % args.n_shards == args.shard]
    commit = repo_commit()
    tag = f"shard{args.shard}of{args.n_shards}"
    units_path = out / f"units.{tag}.jsonl"
    (out / f"config.{tag}.json").write_text(json.dumps({
        "experiment": "E20", "arms": arms, "draws": draws, "seed": args.seed,
        "criterion": CRITERION, "condition": args.condition, "epsilons": EPSILONS,
        "smoke": args.smoke, "n_queries": len(queries), "commit": commit,
        "src_dirty": src_dirty(), "n_shards": args.n_shards,
        "naive_generator": "bkrobust.demo.evaluate.random_sem defaults |w|~U(0.3,1.5) noise var~U(0.5,1.5)",
        "iscm_generator": "population internal standardisation of the naive draw",
        "uumc_generator": "Herman et al. 2025 Algorithm 1, population form",
        "adjustment_set": "optimal_adjustment_set_mpdag(G0)", "budget": args.budget,
    }, indent=2) + "\n")

    done = set()
    n_prior_fail = 0
    if units_path.exists():
        for line in units_path.open():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue  # a line cut by a kill; its unit is redone
            done.add((r["network"], r["x"], r["y"], r["arm"], r["draw"]))
            # A failure written by an earlier, interrupted invocation still fails this run.
            n_prior_fail += bool(r.get("consistency_failures"))
    jobs = []
    for q in queries:
        todo = [(a, k) for a in arms for k in range(draws)
                if (q["network"], q["x"], q["y"], a, k) not in done]
        if todo:
            jobs.append({"q": q, "todo": todo, "commit": commit,
                         "args": {"condition": args.condition, "seed": args.seed,
                                  "budget": args.budget, "radius_budget": args.radius_budget}})
    n_units = sum(len(j["todo"]) for j in jobs)
    print(f"[E20] {len(queries)} queries in {tag}, {len(done)} units already written, "
          f"{n_units} to run, {args.workers} workers, commit {commit[:12]}", flush=True)

    t0 = last = time.time()
    n_written, n_fail = 0, n_prior_fail
    with units_path.open("a") as f, mp.Pool(max(1, args.workers), maxtasksperchild=1) as pool:
        it = pool.imap_unordered(process_query, jobs)
        remaining = len(jobs)
        while remaining:
            try:
                rows = it.next(timeout=30)
            except mp.TimeoutError:
                if time.time() - last >= 60:
                    print(f"[E20] {n_written}/{n_units} units, {remaining} queries running or queued, "
                          f"{time.time() - t0:.0f}s", flush=True)
                    last = time.time()
                continue
            remaining -= 1
            for r in rows:
                f.write(json.dumps(r) + "\n")
                n_written += 1
                if r.get("consistency_failures"):
                    n_fail += 1
                    print(f"[E20] CONSISTENCY {r['network']} {r['x']}->{r['y']} {r['arm']} "
                          f"draw {r['draw']}: {r['consistency_failures']}", flush=True)
            f.flush()
            r0 = rows[0]
            print(f"[E20] done {r0['network']} {r0['x']}->{r0['y']} ({len(rows)} units, "
                  f"status {sorted({r['status'] for r in rows})}); {n_written}/{n_units}, "
                  f"{time.time() - t0:.0f}s", flush=True)
            last = time.time()
    print(f"[E20] finished: {n_written} units written, {n_fail} with consistency failures "
          f"(including {n_prior_fail} from earlier invocations), "
          f"{time.time() - t0:.0f}s", flush=True)
    return 2 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
