"""E19: the App. E.4 bias-radius study (900 instances) under three SEM parameterisations.

Each unit is one stored instance of results/epsilon/study/results.csv (status ok),
rebuilt from its index exactly as experiments/run_r_epsilon.py built it; the same
coefficient draw is then evaluated as naive, posthoc-standardised and internally
standardised (see sem_variants.py and PREREG.md). Validity uses the generalised
adjustment criterion, the default of certify -> breakdown_radius.

    .venv/bin/python experiments/rebuttal/E19_bias_radius_internal_standardisation/run.py --smoke
    .venv/bin/python experiments/rebuttal/E19_bias_radius_internal_standardisation/run.py --workers 8

Exit status: 0 clean, 3 if any unit failed a consistency check (the rows are still written).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
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
for p in (ROOT, ROOT / "src", HERE):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.core.oracle import is_valid  # noqa: E402
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag, random_sem  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.epsilon.bias import ZERO_TOL, knowledge_of, make_context  # noqa: E402
from bkrobust.epsilon.certify import certify  # noqa: E402
from bkrobust.epsilon.profile import r_epsilon, r_epsilon_from_profile, top_bias  # noqa: E402
from bkrobust.synth.generators import GENERATORS  # noqa: E402
from bkrobust.synth.knowledge import CORRUPTIONS, check_consistent, draw_k_true  # noqa: E402
from bkrobust.synth.runner import _pick_treatment_outcome, gate  # noqa: E402
from experiments.run_r_epsilon import (  # noqa: E402 - the study's own grid, seed and caps
    CAP_UNDIRECTED,
    CELLS,
    EPS_RELATIVE,
    MAX_X_UNDIRECTED_DEGREE,
    R0_EPS_ABSOLUTE,
    ROOT_SEED,
)
from sem_variants import ARMS, UNIT_VAR_TOL, make_arm, max_unit_var_error, sortabilities  # noqa: E402

CRITERION = "gac"
EXPERIMENT = "E19"
REL_TOL = 1e-9
EPS_KEYS = [f"{e:g}" for e in EPS_RELATIVE]
# DEVIATIONS.md (2026-11-06): a relative bias equal to a tolerance (to 1e-9) is a tie. The staircase is
# extended one step past the largest tolerance so that both conventions can be read off it.
TIE_TOL = 1e-9
EPS_EXTENSION = max(EPS_RELATIVE) * (1 + 1e-6)


def _first_depth(stair: list[tuple[int, float]], test: Any) -> int:
    for d, b in stair:
        if test(b):
            return d
    return UNREACHED


def _truncate(stair: list[tuple[int, float]], test: Any) -> list[tuple[int, float]]:
    out = []
    for d, b in stair:
        out.append((d, b))
        if test(b):
            break
    return out


def tie_flags(stair: list[tuple[int, float]], scale: float) -> list[str]:
    """Tolerances some staircase value equals to TIE_TOL (relative), e.g. ['1']."""
    return [k for e, k in zip(EPS_RELATIVE, EPS_KEYS, strict=True)
            if any(abs(b - e * scale) <= TIE_TOL * e * scale for _, b in stair)]


class _Timeout(Exception):  # noqa: N818
    pass


def _alarm(signum: int, frame: Any) -> None:
    raise _Timeout()


def _clean(v: Any) -> Any:
    # strict JSON: NaN/inf become null so any reader can parse the file
    if isinstance(v, float) and not math.isfinite(v):
        return None
    if isinstance(v, (np.floating, np.integer, np.bool_)):
        return _clean(v.item())
    if isinstance(v, dict):
        return {str(k): _clean(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_clean(x) for x in v]
    return v


def repo_state() -> dict[str, Any]:
    def git(*a: str) -> str:
        try:
            return subprocess.run(["git", "-C", str(ROOT), *a], capture_output=True, text=True,
                                  check=True).stdout.strip()
        except Exception:
            return "nogit"
    # dirtiness of the code this run executes only; other rebuttal folders are always untracked
    dirty = git("status", "--porcelain", "--", "src", "experiments/run_r_epsilon.py")
    return {"repo_commit": git("rev-parse", "HEAD"), "repo_code_dirty": bool(dirty and dirty != "nogit")}


# --- baseline ------------------------------------------------------------------


def load_baseline(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if r["status"] == "ok"]
    out = []
    for r in rows:
        out.append({
            "instance_id": r["instance_id"],
            "index": int(r["instance_id"].split("_")[1]),
            "generator": r["generator"],
            "n": int(r["n"]),
            "knows_fraction": float(r["knows_fraction"]),
            "corruption_name": r["corruption_name"],
            "x": r["x"], "y": r["y"], "z": r["z"],
            "n_k0": int(r["n_k0"]),
            "r_val": int(r["r_val"]),
            "r0_d": int(r["r0_d"]),
            "theta_z": float(r["theta_z"]),
            "staircase": json.loads(r["staircase_json"]),
            "staircase_shape": r["staircase_shape"],
            "r_eps": {k: int(r[f"r_eps_{k}"]) for k in EPS_KEYS},
            "timing_total_seconds": float(r["timing_total_seconds"]),
        })
    return sorted(out, key=lambda b: b["index"])


def smoke_subset(base: list[dict[str, Any]], k: int = 6) -> list[dict[str, Any]]:
    # one fast instance per code path that matters: step, ramp, flat, ceiling reached / not
    fast = [b for b in base if b["timing_total_seconds"] < 0.5]
    tests = [
        lambda b: b["staircase_shape"] == "step",
        lambda b: b["staircase_shape"] == "ramp",
        lambda b: b["staircase_shape"] == "flat",
        lambda b: b["r_eps"]["1"] == UNREACHED,
        lambda b: b["r_eps"]["1"] != UNREACHED,
        lambda b: b["r_eps"]["0.01"] > b["r_val"],
    ]
    picked: list[dict[str, Any]] = []
    for t in tests:
        hit = next((b for b in fast if t(b) and b not in picked), None)
        if hit is not None:
            picked.append(hit)
    return picked[:k]


# --- one unit ------------------------------------------------------------------


def rebuild(index: int) -> tuple[dict[str, Any] | None, str | None]:
    """Lines 208-247 of experiments/run_r_epsilon.py, same RNG consumption order."""
    (gen, params), n, knows, (corr, rate) = CELLS[index % len(CELLS)]
    rng = np.random.default_rng(ROOT_SEED + index)
    if gen == "decoupled_backdoor" and n < 7:
        return None, "generator_incompatible_n"
    dag = GENERATORS[gen](n, rng, **params)
    cpdag = dag_to_cpdag(dag)
    x, y = _pick_treatment_outcome(gen, dag, rng)
    if len(cpdag.undirected_edges) > CAP_UNDIRECTED:
        return None, "cpdag_too_large"
    if len(cpdag.neighbors(x)) > MAX_X_UNDIRECTED_DEGREE:
        return None, "x_undirected_degree_too_high"
    ok, reason = gate(dag, cpdag, x, y)
    if not ok:
        return None, reason
    k_true = draw_k_true(dag, cpdag, rng, knows)
    k_assumed = list(k_true) if corr == "none" else CORRUPTIONS[corr](k_true, rng, rate)
    if not check_consistent(cpdag, k_assumed):
        return None, "k_assumed_inconsistent"
    g0 = apply_orientations(cpdag, k_assumed)
    if g0 is None:
        return None, "g0_none"
    z = optimal_adjustment_set_mpdag(g0, x, y)
    if z is None:
        return None, "z_not_identified"
    zf = frozenset(z)
    if not is_valid(zf, g0, x, y):
        return None, "z_invalid_at_g0"
    sem0 = random_sem(dag, rng)  # last RNG use, as in the study
    return {"gen": gen, "params": params, "n": n, "knows": knows, "corr": corr, "rate": rate,
            "dag": dag, "cpdag": cpdag, "x": x, "y": y, "k_assumed": k_assumed, "g0": g0,
            "z": zf, "sem0": sem0}, None


def shape_of(seq: list[float], r_val_int: int, n_k0: int) -> str:
    # copied from run_r_epsilon.py, including its never-true UNREACHED test, so the
    # classification is the one App. E.4 printed
    increases = sum(1 for i in range(1, len(seq)) if seq[i] > seq[i - 1] + 1e-12)
    shape = "flat" if len(seq) < 2 or increases == 0 else ("step" if increases == 1 else "ramp")
    if r_val_int == UNREACHED or n_k0 == 0:
        shape = "n/a"
    return shape


def run_arm(inst: dict[str, Any], arm: str, timeout_s: int) -> dict[str, Any]:
    out: dict[str, Any] = {"status": "ok", "notes": ""}
    t0 = time.perf_counter()
    old = signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(timeout_s)
    try:
        sem = make_arm(inst["sem0"], arm)
        out["max_unit_var_error"] = max_unit_var_error(sem)
        out.update(sortabilities(sem))
        cpdag, g0, x, y, zf = inst["cpdag"], inst["g0"], inst["x"], inst["y"], inst["z"]
        n_k0 = len(knowledge_of(cpdag, g0))
        cert = certify(cpdag, list(inst["k_assumed"]), x, y, sem=sem, z=zf, g0=g0,
                       epsilons=list(EPS_RELATIVE) + [EPS_EXTENSION], units="relative",
                       search_budget=max(3, n_k0))
        ctx = make_context(sem, cpdag, x, y, zf)
        r_val_int = cert.r_val if cert.r_val != UNREACHED else 0
        scale = abs(ctx.theta_z)  # certify's own conversion for units="relative"
        top = max(EPS_RELATIVE) * scale
        full = [(s.d, s.beta_up) for s in cert.profile]
        # the staircase the study stored: the library stops after the first step strictly above the largest tolerance
        stair = _truncate(full, lambda b: b > top)
        stair_ge = _truncate(full, lambda b: b >= top * (1 - TIE_TOL))
        stair_gt = _truncate(full, lambda b: b > top * (1 + TIE_TOL))
        seq = [b for _, b in stair]
        out.update({
            "r_val": cert.r_val, "r_val_exact": cert.status == "ok", "r_val_method": cert.method,
            "theta_z": cert.theta_z, "tau_true": ctx.tau_true,
            "r_eps": {f"{e:g}": cert.r_eps[e] for e in EPS_RELATIVE},
            "r_eps_ge": {k: _first_depth(full, lambda b, e=e: b >= e * scale * (1 - TIE_TOL))
                         for e, k in zip(EPS_RELATIVE, EPS_KEYS, strict=True)},
            "r_eps_gt": {k: _first_depth(full, lambda b, e=e: b > e * scale * (1 + TIE_TOL))
                         for e, k in zip(EPS_RELATIVE, EPS_KEYS, strict=True)},
            "staircase_shape_ge": shape_of([b for _, b in stair_ge], r_val_int, n_k0),
            "staircase_shape_gt": shape_of([b for _, b in stair_gt], r_val_int, n_k0),
            "threshold_ties": tie_flags(full, scale),
            "staircase_extended_len": len(full),
            "staircase": stair,
            "staircase_rel": [b / abs(cert.theta_z) for b in seq],
            "r0_d": r_epsilon_from_profile(cert.profile, R0_EPS_ABSOLUTE),
            "beta_top_rel": top_bias(ctx).worst / abs(cert.theta_z),
            "monotone_violation": any(seq[i] > seq[i + 1] + 1e-9 for i in range(len(seq) - 1)),
            "staircase_shape": shape_of(seq, r_val_int, n_k0),
        })
        # certificate audit, as in the study
        k_false = sum(1 for a, b in knowledge_of(cpdag, g0) if (a, b) not in inst["dag"].directed_edges)
        realised = abs(cert.theta_z - ctx.tau_true)
        beta_k = max([0.0] + [b for d, b in stair if d <= k_false])
        out.update({
            "k_false": k_false,
            "beta_up_at_kfalse_is_lower_bound": k_false > (stair[-1][0] if stair else -1),
            "certificate_holds": bool(realised <= beta_k + 1e-9),
            "conservativeness_ratio": realised / beta_k if beta_k > ZERO_TOL else float("nan"),
        })
        # two exact strategies must agree; a disagreement is a bug in the library
        eff = 1.00 * cert.eps_scale
        inc = r_epsilon(ctx, g0, eff, r_val=r_val_int, strategy="incremental").radius
        bis = r_epsilon(ctx, g0, eff, r_val=r_val_int, strategy="bisection").radius
        out.update({"radius_incremental_eps1": inc, "radius_bisection_eps1": bis})
    except _Timeout:
        out["status"] = "capped_timeout"
        out["notes"] = f"timeout {timeout_s}s"
    except Exception as exc:
        out["status"] = "run_error"
        out["notes"] = f"{type(exc).__name__}: {exc}"
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)
        out["seconds"] = time.perf_counter() - t0
    return out


def _close(a: float, b: float) -> bool:
    return abs(a - b) <= REL_TOL * max(1.0, abs(a), abs(b))


def _near_tie(rel_seq: list[float]) -> bool:
    return any(abs(v - e) <= 1e-7 * max(1.0, e) for v in rel_seq for e in EPS_RELATIVE)


TIE_SENSITIVE = ("naive vs baseline r_eps", "naive vs baseline shape", "naive vs baseline staircase differs",
                 "posthoc vs naive shape", "posthoc vs naive r_eps")


def checks(base: dict[str, Any], inst: dict[str, Any], arms: dict[str, dict[str, Any]]) -> tuple[list[str], list[str]]:
    """Returns (failures, tie_only): differences explained by a relative bias equal to a tolerance go to tie_only."""
    fails: list[str] = []
    ident = {"x": inst["x"], "y": inst["y"], "z": "|".join(sorted(inst["z"])),
             "generator": inst["gen"], "n": inst["n"],
             "n_k0": len(knowledge_of(inst["cpdag"], inst["g0"]))}
    for k, v in ident.items():
        if base[k] != v:
            fails.append(f"identity_{k}: baseline {base[k]!r} rebuilt {v!r}")
    ok = {a: r for a, r in arms.items() if r["status"] == "ok"}
    for a, r in ok.items():
        if r["radius_incremental_eps1"] != r["radius_bisection_eps1"]:
            fails.append(f"{a}: incremental {r['radius_incremental_eps1']} != bisection {r['radius_bisection_eps1']}")
        if r["monotone_violation"]:
            fails.append(f"{a}: staircase not monotone")
        if a != "naive" and r["max_unit_var_error"] > UNIT_VAR_TOL:
            fails.append(f"{a}: variance off unity by {r['max_unit_var_error']:.3g}")
    rvals = {r["r_val"] for r in ok.values()}
    if len(rvals) > 1:
        fails.append(f"r_val differs across arms: { {a: r['r_val'] for a, r in ok.items()} }")
    nv = ok.get("naive")
    if nv is not None:
        if nv["r_val"] != base["r_val"] or nv["r0_d"] != base["r0_d"]:
            fails.append(f"naive vs baseline r_val/r0: {nv['r_val']},{nv['r0_d']} vs {base['r_val']},{base['r0_d']}")
        if nv["r_eps"] != base["r_eps"]:
            fails.append(f"naive vs baseline r_eps: {nv['r_eps']} vs {base['r_eps']}")
        if nv["staircase_shape"] != base["staircase_shape"]:
            fails.append(f"naive vs baseline shape: {nv['staircase_shape']} vs {base['staircase_shape']}")
        if not _close(nv["theta_z"], base["theta_z"]):
            fails.append(f"naive vs baseline theta_z: {nv['theta_z']} vs {base['theta_z']}")
        bs = base["staircase"]
        if len(bs) != len(nv["staircase"]) or not all(
                d1 == d2 and _close(b1, b2) for (d1, b1), (d2, b2) in zip(bs, nv["staircase"], strict=False)):
            fails.append("naive vs baseline staircase differs")
    ph = ok.get("posthoc")
    if nv is not None and ph is not None:
        # relative bias is scale-free, so this arm must be the naive arm up to rounding
        if ph["r_eps"] != nv["r_eps"] and not _near_tie(nv["staircase_rel"]):
            fails.append(f"posthoc vs naive r_eps: {ph['r_eps']} vs {nv['r_eps']}")
        if ph["staircase_shape"] != nv["staircase_shape"]:
            fails.append(f"posthoc vs naive shape: {ph['staircase_shape']} vs {nv['staircase_shape']}")
        sp, sn = ph["staircase_rel"], nv["staircase_rel"]
        n_common = min(len(sp), len(sn))
        if not all(abs(a - b) <= 1e-7 * max(1.0, abs(b)) for a, b in zip(sp[:n_common], sn[:n_common], strict=True)):
            fails.append("posthoc vs naive relative staircase differs")
    # DEVIATIONS.md: a difference in r_eps, shape or staircase length is a tie, not a failure, when the naive
    # or the stored staircase holds a relative bias equal to a tolerance (to 1e-9) and nothing else differs
    ties = set(nv["threshold_ties"]) if nv is not None else set()
    if nv is not None and abs(base["theta_z"]) > 1e-12:
        ties |= set(tie_flags([(d, b) for d, b in base["staircase"]], abs(base["theta_z"])))
    tie_only: list[str] = []
    if ties and nv is not None:
        diff_keys = {k for k in EPS_KEYS if nv["r_eps"][k] != base["r_eps"][k]}
        if ph is not None:
            diff_keys |= {k for k in EPS_KEYS if ph["r_eps"][k] != nv["r_eps"][k]}
        top_key = EPS_KEYS[EPS_RELATIVE.index(max(EPS_RELATIVE))]
        # r_eps may differ only at tied tolerances; shape and staircase length only when the largest one is tied
        explained = diff_keys <= ties and (top_key in ties or not any(
            f.startswith(("naive vs baseline shape", "naive vs baseline staircase", "posthoc vs naive shape"))
            for f in fails))
        if explained:
            tie_only = [f for f in fails if f.startswith(TIE_SENSITIVE)]
            fails = [f for f in fails if not f.startswith(TIE_SENSITIVE)]
    return fails, tie_only


def process_unit(task: tuple[dict[str, Any], int, dict[str, Any]]) -> dict[str, Any]:
    base, timeout_s, repo = task
    t0 = time.perf_counter()
    row: dict[str, Any] = {"experiment": EXPERIMENT, "instance_id": base["instance_id"],
                           "index": base["index"], "criterion": CRITERION, **repo}
    inst, reason = rebuild(base["index"])
    if inst is None:
        row.update({"status_all_ok": False, "arms": {},
                    "consistency_failures": [f"rebuild rejected the stored instance: {reason}"]})
        return _clean(row)
    row.update({"generator": inst["gen"], "generator_params": inst["params"], "n": inst["n"],
                "knows_fraction": inst["knows"], "corruption_name": inst["corr"],
                "corruption_rate": inst["rate"], "x": inst["x"], "y": inst["y"],
                "z": "|".join(sorted(inst["z"])),
                "n_k0": len(knowledge_of(inst["cpdag"], inst["g0"]))})
    arms = {a: run_arm(inst, a, timeout_s) for a in ARMS}
    row["arms"] = arms
    row["status_all_ok"] = all(r["status"] == "ok" for r in arms.values())
    row["consistency_failures"], row["tie_only_differences"] = checks(base, inst, arms)
    row["seconds"] = time.perf_counter() - t0
    return _clean(row)


# --- driver ----------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--smoke", action="store_true", help="six fast instances, every code path")
    ap.add_argument("--out", default=str(ROOT / "results" / "rebuttal" / EXPERIMENT))
    ap.add_argument("--baseline-csv", default=str(ROOT / "results" / "epsilon" / "study" / "results.csv"))
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--timeout", type=int, default=60, help="seconds per arm per instance")
    ap.add_argument("--limit", type=int, default=0, help="first N units only (0 = all)")
    args = ap.parse_args()
    if not 0 <= args.shard < args.n_shards:
        ap.error("need 0 <= --shard < --n-shards")

    base = load_baseline(Path(args.baseline_csv))
    if args.smoke:
        base = smoke_subset(base)
        args.timeout = min(args.timeout, 20)
    out_dir = Path(args.out) / ("smoke" if args.smoke else "")
    base = [b for i, b in enumerate(base) if i % args.n_shards == args.shard]
    if args.limit:
        base = base[: args.limit]

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"units.shard{args.shard}of{args.n_shards}.jsonl"
    done: set[str] = set()
    if out_path.exists():
        with out_path.open() as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["instance_id"])
                except (json.JSONDecodeError, KeyError):
                    pass  # a line cut by a kill; the unit is simply redone
    todo = [b for b in base if b["instance_id"] not in done]
    repo = repo_state()
    # what the fidelity gate reads: one config per shard, written before any unit runs
    cfg = {"experiment": EXPERIMENT, "criterion": CRITERION, "arms": list(ARMS),
           "eps_relative": list(EPS_RELATIVE), "root_seed": ROOT_SEED, "n_cells": len(CELLS),
           "sem_generator": "bkrobust.demo.evaluate.random_sem defaults",
           "weight_range": [0.3, 1.5], "noise_range": [0.5, 1.5],
           "timeout_s": args.timeout, "smoke": args.smoke, "shard": args.shard, "n_shards": args.n_shards,
           "n_units_in_shard": len(base), "limit": args.limit,
           "baseline_csv": str(args.baseline_csv),
           "baseline_csv_sha256": hashlib.sha256(Path(args.baseline_csv).read_bytes()).hexdigest(),
           "python": sys.version.split()[0], "numpy": np.__version__, **repo}
    (out_dir / f"run_config.shard{args.shard}of{args.n_shards}.json").write_text(json.dumps(cfg, indent=2))
    print(f"[E19] {len(base)} units in shard {args.shard}/{args.n_shards}, {len(done)} already written, "
          f"{len(todo)} to do; workers={args.workers} timeout={args.timeout}s criterion={CRITERION} "
          f"commit={repo['repo_commit']} dirty={repo['repo_code_dirty']} -> {out_path}", flush=True)

    counters = {"written": 0, "failures": 0, "not_ok": 0}
    t_start = time.perf_counter()
    stop = threading.Event()

    def heartbeat() -> None:
        while not stop.wait(60):
            el = time.perf_counter() - t_start
            print(f"[E19 {el:7.0f}s] written {counters['written']}/{len(todo)} "
                  f"not_ok={counters['not_ok']} consistency_failures={counters['failures']}", flush=True)

    threading.Thread(target=heartbeat, daemon=True).start()
    tasks = [(b, args.timeout, repo) for b in todo]
    with out_path.open("a") as fh:
        def write(row: dict[str, Any]) -> None:
            fh.write(json.dumps(row) + "\n")
            fh.flush()
            counters["written"] += 1
            counters["not_ok"] += 0 if row.get("status_all_ok") else 1
            if row["consistency_failures"]:
                counters["failures"] += 1
                print(f"[E19] CONSISTENCY FAILURE {row['instance_id']}: {row['consistency_failures']}", flush=True)

        if args.workers <= 1:
            for t in tasks:
                write(process_unit(t))
        else:
            with mp.get_context("spawn").Pool(args.workers, maxtasksperchild=50) as pool:
                for row in pool.imap_unordered(process_unit, tasks, chunksize=1):
                    write(row)
    stop.set()
    # the exit status covers the whole file, so a resumed run cannot hide an earlier failure
    with out_path.open() as fh:
        rows = [json.loads(line) for line in fh if line.strip()]
    n_fail = sum(1 for r in rows if r.get("consistency_failures"))
    n_tie = sum(1 for r in rows if r.get("tie_only_differences"))
    print(f"[E19] rows whose only differences are threshold ties (counted, not failures): {n_tie}", flush=True)
    el = time.perf_counter() - t_start
    print(f"[E19] done in {el:.0f}s: written now {counters['written']}, rows in file {len(rows)}, "
          f"not_ok now {counters['not_ok']}, rows with consistency failures {n_fail}", flush=True)
    return 3 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
