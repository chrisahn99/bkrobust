"""E14: does the Su-Henckel coefficient-agreement test flag committed sets invalid at truth?

For each distinct released state (network, x, y, knowledge) with a committed set: rebuild G0 as
the released sweep did and check it against the row; list sets GAC-valid in G0 (the committed
set first, at most --max-sets); decide each set's validity in the true DAG; then, for each
replicate, draw a non-Gaussian linear SEM on the true DAG, sample n rows for each n, and run the
test across the sets. One JSON line per state. See PREREG.md.

    python experiments/rebuttal/E14_data_test_for_invalid_sets/run.py --smoke
    python experiments/rebuttal/E14_data_test_for_invalid_sets/run.py --workers 8
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

from bkrobust.demo.evaluate import is_valid_adjustment_set_dag  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.gac.dag_level import is_gac_valid_dag  # noqa: E402
from bkrobust.gac.mpdag_level import clear_cache, gac_forbidden_set, is_gac_valid_mpdag  # noqa: E402
from coef_test import agreement_test, coefficient_and_influence, pairwise_bonferroni  # noqa: E402
from sem_sampler import (ancestral_closure, check_covariance, draw_sem, population_coefficient,  # noqa: E402
                         sample)
from vendor_ledger import RELEASE_SEED, knowledge_for_row, load  # noqa: E402
from vendor_optimal import committed_adjustment_set, possible_ancestors  # noqa: E402

DATA = ROOT / "experiments" / "rebuttal" / "_data"
CRITERION = "gac"
SEED = 20261106
N_VALUES = (1000, 10000)
REPS = 50
ALPHA = 0.05
MAX_SETS = 20
GREEDY_ORDERS = 3
MAX_SET_FRACTION = 0.1  # a set with more than n/10 members is dropped at that n
POP_TOL = 1e-6  # population coefficients are unit-scale (standardised SEM)
COV_TOL = 1e-8
SMOKE_NETWORKS = ("asia", "sachs")
SMOKE_PER_CLASS = 4
SMOKE_REPS = 3


def sha16(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()[:16]


def file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def state_key(r: dict) -> str:
    return "|".join((r["network"], r["x"], r["y"], r["k_sha256"]))


def repo_commit() -> tuple[str, bool | None]:
    try:
        head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--", "src"],
                               capture_output=True, text=True, check=True).stdout.strip()
        return head, bool(dirty)
    except (OSError, subprocess.CalledProcessError):
        return "unknown", None


def radius_stratum(row: dict) -> str:
    r = int(row["r_hop"]) if row["r_hop"] not in ("", None) else -1
    if r < 0:
        return "no_failure_found"
    return "3+" if r >= 3 else str(r)


# ------------------------------------------------------------------ candidate sets


def greedy_minimal(g0, x, y, start: frozenset, rng: np.random.Generator) -> frozenset:
    """Drop members of a valid set one at a time in a random order while it stays valid."""
    cur = set(start)
    for v in rng.permutation(sorted(start)):
        trial = frozenset(cur - {v})
        if is_gac_valid_mpdag(g0, x, y, trial):
            cur.discard(v)
    return frozenset(cur)


def candidate_sets(g0, x, y, z: frozenset, seed: list, max_sets: int) -> tuple[list, dict]:
    """Sets GAC-valid in G0, committed set first, then a fixed proposal order, stopping at the cap.

    Proposals are built lazily: on the large networks a greedy pass costs |canonical| validity
    calls, so nothing past the cap is computed.
    """
    forb = gac_forbidden_set(g0, x, y)
    canon = frozenset((possible_ancestors(g0, x) | possible_ancestors(g0, y)) - forb - {x, y})
    canon_valid = is_gac_valid_mpdag(g0, x, y, canon)

    def proposals():
        yield "committed", z
        yield "canonical", canon
        yield "parents_x", frozenset(g0.parents(x)) - {x, y}
        if canon_valid:
            for o in range(GREEDY_ORDERS):
                yield f"minimal_{o}", greedy_minimal(g0, x, y, canon, np.random.default_rng(seed + [o]))
        for v in sorted(z):
            yield f"drop_{v}", z - {v}
        for v in sorted(canon - z):
            yield f"add_{v}", z | {v}

    sets, origin, seen = [], [], set()
    for name, s in proposals():
        if s in seen:
            continue
        seen.add(s)
        if is_gac_valid_mpdag(g0, x, y, s):
            sets.append(s)
            origin.append(name)
            if len(sets) == max_sets:
                break
    return sets, {"origin": origin, "n_proposals_tested": len(seen), "canonical_size": len(canon),
                  "canonical_valid": canon_valid, "capped": len(sets) == max_sets}


# ------------------------------------------------------------------ one unit

_NETS: dict = {}


def network(models: str, name: str):
    if name not in _NETS:
        # One network per worker at a time keeps memory flat. The criterion's caches are kept
        # within a network: states of one arm share G0, and its possible-ancestor sets are the
        # expensive part on the large networks.
        _NETS.clear()
        clear_cache()
        _NETS[name] = load(Path(models), name)
    return _NETS[name]


def run_unit(task: tuple) -> dict:
    row, cfg = task
    t0 = time.perf_counter()
    rec = {"key": state_key(row), "network": row["network"], "x": row["x"], "y": row["y"],
           "k_sha256": row["k_sha256"], "arms": row["_arms"], "criterion": CRITERION,
           "commit": cfg["commit"], "src_dirty": cfg["src_dirty"], "smoke": cfg["smoke"],
           "radius_stratum": radius_stratum(row),
           "released": {f: row[f] for f in ("arm", "rep", "r_hop", "hop_method", "hop_lower_bound",
                                             "r_claim", "r_claim_status", "z_size", "z_valid_at_truth",
                                             "z_backdoor_at_truth", "severity")}}
    try:
        rec.update(compute(row, cfg))
    except Exception:  # reported, and retried on resume
        rec["status"] = "error"
        rec["error"] = traceback.format_exc()[-2000:]
    rec["seconds"] = round(time.perf_counter() - t0, 3)
    return rec


def compute(row: dict, cfg: dict) -> dict:
    dag, cpdag, model_file = network(cfg["models"], row["network"])
    x, y = row["x"], row["y"]
    res: dict = {"model_file": model_file, "n_nodes": len(dag.nodes)}
    mism = []
    if sha16(cpdag.edge_string()) != row["chat_source_sha256"]:
        mism.append("chat_source_sha256")
    if sha16(dag.edge_string()) != row["generating_dag_sha256"]:
        mism.append("generating_dag_sha256")
    k = knowledge_for_row(row, dag, cpdag, cfg["know"], cfg["release_seed"])
    if k is None or sha16(json.dumps(sorted(k))) != row["k_sha256"] or len(k) != int(row["n_k"]):
        mism.append("knowledge")
    g0 = apply_orientations(cpdag, k or [])
    cs = committed_adjustment_set(g0, x, y) if g0 is not None else None
    if g0 is None or cs is None or cs.z is None:
        mism.append("no_committed_set")
    elif cs.verdict != row["committed_verdict"] or len(cs.z) != int(row["z_size"]):
        mism.append("committed_set")
    if mism:
        return {**res, "status": "rebuild_mismatch", "mismatch": mism}
    z = cs.z
    hard: list[str] = []
    uh = int(hashlib.sha256(state_key(row).encode()).hexdigest()[:8], 16)
    sets, meta = candidate_sets(g0, x, y, z, [cfg["seed"], uh], cfg["max_sets"])
    res.update({"z": sorted(z), "sets": [sorted(s) for s in sets], **meta})
    if not sets or sets[0] != z:
        hard.append("committed set not GAC-valid in G0")
    for s in sets:
        if not is_gac_valid_mpdag(g0, x, y, s):
            hard.append("candidate set not GAC-valid in G0")
    truth = [is_gac_valid_dag(dag, x, y, s) for s in sets]
    if any(t != is_gac_valid_mpdag(dag, x, y, s) for t, s in zip(truth, sets)):
        hard.append("DAG-level and MPDAG-level GAC disagree on the true DAG")
    res["valid_at_truth"] = truth
    res["z_valid_at_truth"] = bool(truth[0]) if truth else None
    res["z_backdoor_at_truth"] = bool(is_valid_adjustment_set_dag(dag, x, y, z))
    res["matches_released_truth"] = res["z_valid_at_truth"] == (row["z_valid_at_truth"] == "1")
    res["all_valid_at_truth"] = bool(all(truth))
    res["hard_failures"] = hard
    if len(sets) < 2:
        res["status"] = "fewer_than_two_sets"
        return res

    keep = ancestral_closure(dag, {x, y}.union(*sets))
    res["n_simulated_nodes"] = len(keep)
    reps = []
    for r in range(cfg["reps"]):
        model = draw_sem(dag, keep, np.random.default_rng([cfg["seed"], uh, r]))
        gap = check_covariance(model)
        tau = model["sem"].true_total_effect(x, y)
        pop = [population_coefficient(model, x, y, s) for s in sets]
        if gap > COV_TOL:
            hard.append(f"rep {r}: covariance gap {gap:.2e}")
        for s_ok, b in zip(truth, pop):
            if s_ok and abs(b - tau) > POP_TOL:
                hard.append(f"rep {r}: a truth-valid set has population coefficient {b:.6g} != effect {tau:.6g}")
                break
        spread = max(abs(b - pop[0]) for b in pop[1:])
        rep = {"rep": r, "tau": tau, "beta_z_pop": pop[0], "bias_z_pop": pop[0] - tau,
               "pop_spread": spread, "pop_detectable": bool(spread > POP_TOL), "by_n": {}}
        for n in cfg["n_values"]:
            usable = [i for i, s in enumerate(sets) if len(s) <= MAX_SET_FRACTION * n]
            if not usable or usable[0] != 0 or len(usable) < 2:
                rep["by_n"][str(n)] = {"skipped": "fewer than two usable sets", "n_sets": len(usable)}
                continue
            data = sample(model, n, np.random.default_rng([cfg["seed"], uh, r, n]))
            col = model["idx"]
            xv, yv = data[:, col[x]], data[:, col[y]]
            betas, psis = [], []
            for i in usable:
                zc = [col[v] for v in sorted(sets[i])]
                b, psi = coefficient_and_influence(xv, yv, data[:, zc])
                betas.append(b)
                psis.append(psi)
            betas_a, psi_m = np.array(betas), np.column_stack(psis)
            tol = n ** -0.5
            w = agreement_test(betas_a, psi_m, tol, cfg["alpha"])
            pw = pairwise_bonferroni(betas_a, psi_m, tol, cfg["alpha"])
            rep["by_n"][str(n)] = {"n_sets": len(usable), "beta_z": betas[0], "stat": w["stat"],
                                   "df": w["df"], "p": w["p"], "reject": w["reject"],
                                   "pairwise_m": pw["m"], "pairwise_min_p": pw["min_p"],
                                   "pairwise_reject": pw["reject"]}
        reps.append(rep)
    res["reps"] = reps
    summary = {}
    for n in cfg["n_values"]:
        rows = [rp["by_n"][str(n)] for rp in reps if "reject" in rp["by_n"][str(n)]]
        summary[str(n)] = {"n_reps": len(rows),
                           "reject_rate": float(np.mean([q["reject"] for q in rows])) if rows else None,
                           "pairwise_reject_rate": float(np.mean([q["pairwise_reject"] for q in rows])) if rows else None,
                           "mean_df": float(np.mean([q["df"] for q in rows])) if rows else None}
    res["summary"] = summary
    res["pop_detectable_rate"] = float(np.mean([rp["pop_detectable"] for rp in reps]))
    res["status"] = "ok"
    return res


# ------------------------------------------------------------------ driver


def select_states(records: Path, smoke: bool) -> list[dict]:
    """Released ok rows, one per distinct state, the first in sorted (arm, rep) order kept."""
    rows = [r for r in csv.DictReader(records.open()) if r["status"] == "ok"]
    rows.sort(key=lambda r: (state_key(r), r["arm"], r["rep"]))
    states: dict[str, dict] = {}
    for r in rows:
        k = state_key(r)
        if k not in states:
            states[k] = {**r, "_arms": []}
        states[k]["_arms"].append(f"{r['arm']}:{r['rep']}")
    out = sorted(states.values(), key=state_key)
    if smoke:
        small = [r for r in out if r["network"] in SMOKE_NETWORKS]
        bad = [r for r in small if r["z_valid_at_truth"] == "0"][:SMOKE_PER_CLASS]
        good = [r for r in small if r["z_valid_at_truth"] == "1"][:SMOKE_PER_CLASS]
        out = bad + good
    return out


def done_keys(out: Path) -> set[str]:
    keys = set()
    for f in sorted(out.glob("units*.jsonl")):
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
    ap.add_argument("--out", type=Path, default=ROOT / "results/rebuttal/E14")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--smoke", action="store_true", help="8 states of asia and sachs, 3 replicates, < 60 s")
    ap.add_argument("--reps", type=int, default=REPS)
    ap.add_argument("--n-values", type=int, nargs="+", default=list(N_VALUES))
    ap.add_argument("--alpha", type=float, default=ALPHA)
    ap.add_argument("--max-sets", type=int, default=MAX_SETS)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--release-seed", type=int, default=RELEASE_SEED, help="D_RAND seed of the released sweep")
    args = ap.parse_args(argv)
    if not 0 <= args.shard < args.n_shards:
        ap.error("--shard must be in [0, --n-shards)")
    if args.smoke:
        args.reps = SMOKE_REPS
        if args.out == ROOT / "results/rebuttal/E14":
            args.out = ROOT / "results/rebuttal/E14_smoke"

    commit, dirty = repo_commit()
    states = select_states(args.records, args.smoke)
    n_selected = len(states)
    states = [r for i, r in enumerate(states) if i % args.n_shards == args.shard]
    args.out.mkdir(parents=True, exist_ok=True)
    tag = "" if args.n_shards == 1 else f".shard{args.shard}of{args.n_shards}"
    config = {"experiment": "E14", "criterion": CRITERION, "commit": commit, "src_dirty": dirty,
              "smoke": args.smoke, "population": "status == ok, distinct (network, x, y, k_sha256)",
              "n_selected": n_selected, "n_shard_units": len(states), "shard": args.shard,
              "n_shards": args.n_shards, "reps": args.reps, "n_values": args.n_values,
              "alpha": args.alpha, "max_sets": args.max_sets, "greedy_orders": GREEDY_ORDERS,
              "max_set_fraction": MAX_SET_FRACTION, "rank_rule": "eigenvalue > n^-0.5 * var(psi_Z)",
              "noise": "uniform|laplace|exponential per node, unit-variance SEM",
              "seed": args.seed, "release_seed": args.release_seed, "workers": args.workers,
              "records_sha256": file_sha(args.records), "knowledge_sha256": file_sha(args.knowledge)}
    (args.out / f"config{tag}.json").write_text(json.dumps(config, indent=1))

    skip = done_keys(args.out)
    todo = [r for r in states if state_key(r) not in skip]
    print(f"[E14] {len(states)} states in shard {args.shard}/{args.n_shards} ({n_selected} selected), "
          f"{len(states) - len(todo)} already done, {len(todo)} to run; commit {commit[:12]}"
          f"{' (src dirty)' if dirty else ''}", flush=True)
    cfg = {"know": json.loads(args.knowledge.read_text()), "models": str(args.models),
           "seed": args.seed, "release_seed": args.release_seed, "commit": commit,
           "src_dirty": dirty, "smoke": args.smoke, "reps": args.reps,
           "n_values": args.n_values, "alpha": args.alpha, "max_sets": args.max_sets}
    todo.sort(key=lambda r: (r["network"], state_key(r)))  # a worker reuses its loaded network
    tasks = [(r, cfg) for r in todo]
    n_hard = n_err = n_ok = 0
    t0 = last = time.perf_counter()
    with (args.out / f"units{tag}.jsonl").open("a") as fh:
        pool = mp.get_context("spawn").Pool(max(1, args.workers))
        results = pool.imap_unordered(run_unit, tasks, chunksize=1)
        done = 0
        while done < len(tasks):
            try:
                rec = results.next(timeout=60)
            except mp.TimeoutError:
                print(f"[E14] {done}/{len(tasks)} done, {time.perf_counter() - t0:.0f} s, waiting", flush=True)
                continue
            done += 1
            fh.write(json.dumps(rec, default=str) + "\n")
            fh.flush()
            n_err += rec.get("status") == "error"
            n_ok += rec.get("status") == "ok"
            if rec.get("hard_failures"):
                n_hard += 1
                print(f"[E14] HARD {rec['key']}: {rec['hard_failures'][:3]}", flush=True)
            if rec.get("status") == "error":
                print(f"[E14] ERROR {rec['key']}: {rec['error'].strip().splitlines()[-1]}", flush=True)
            now = time.perf_counter()
            if now - last >= 30 or done == len(tasks) or args.smoke:
                last = now
                print(f"[E14] {done}/{len(tasks)} done, ok {n_ok}, hard {n_hard}, errors {n_err}, "
                      f"{now - t0:.0f} s", flush=True)
        pool.close()
        pool.join()
    if n_hard:
        print(f"[E14] CONSISTENCY FAILURE on {n_hard} states", flush=True)
        return 1
    if n_err:
        print(f"[E14] {n_err} states raised; rerun to retry them", flush=True)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
