"""Empirical validation study for the epsilon-bias radius (`docs/R_EPSILON_THEORY.md`).

Mirrors the Axis-A instance pipeline of `bkrobust.synth.runner.run_instance`, but
computes the *bias* certificate (`bkrobust.epsilon.certify.certify`) instead of the
old sampled `mean_abs_bias` radius, and audits it against the true DAG.

This script owns `experiments/run_r_epsilon.py` and everything under
`results/epsilon/study/`. It does not modify `src/bkrobust/epsilon/*`; any bug
found there is reported in the console summary at the end of the run, not fixed
here.

Usage::

    .venv/bin/python experiments/run_r_epsilon.py [--target 1200] [--time-budget 2100]

Writes:
    results/epsilon/study/results.csv       -- one row per accepted instance
    results/epsilon/study/manifest.json     -- seed, grid, counts, git sha, timings
    results/epsilon/study/analysis.json     -- analyses 1-6, machine readable
    results/epsilon/study/SUMMARY.md        -- analyses 1-6, human readable
    results/epsilon/study/examples.txt      -- describe(cert) for a handful of instances
    results/epsilon/study/hard_errors.jsonl -- any incremental/bisection disagreement
"""

from __future__ import annotations

import argparse
import itertools
import json
import signal
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.core.oracle import bias_stats, is_valid  # noqa: E402
from bkrobust.core.resultsio import ResultWriter, git_sha, write_manifest  # noqa: E402
from bkrobust.core.spacelib import distances_from  # noqa: E402
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag, random_sem  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.epsilon.bias import ZERO_TOL, knowledge_of, make_context  # noqa: E402
from bkrobust.epsilon.certify import certify  # noqa: E402
from bkrobust.epsilon.profile import (  # noqa: E402
    greedy_chain_bound,
    r_epsilon,
    r_epsilon_from_profile,
    retraction_shell,
    top_bias,
)
from bkrobust.search.space_fixed import build_space_correct  # noqa: E402
from bkrobust.synth.generators import GENERATORS  # noqa: E402
from bkrobust.synth.knowledge import CORRUPTIONS, check_consistent, draw_k_true  # noqa: E402
from bkrobust.synth.runner import (  # noqa: E402
    _pick_treatment_outcome,
    gate,
)
from bkrobust.synth.runner import (  # noqa: E402 - follows the deliberate sys.path setup block above
    _radius_eps as old_radius_eps,
)

OUT_DIR = ROOT / "results" / "epsilon" / "study"
ROOT_SEED = 20260916  # today's date, arbitrary but fixed

# --- the grid ---------------------------------------------------------------

CAP_UNDIRECTED = 9  # task cap is <=10; 9 matches the project's own established safe
# cutoff (bkrobust.synth.runner.MAX_UNDIRECTED_EDGES, whose docstring measures
# k=9 ~1s vs k=10 ~12s for the OLD full-space enumeration; the epsilon module's
# per-shell semi-local evaluator is cheaper per-state but the same combinatorial
# wall lands in the same place once |K_G0| and deg_undirected(X) both grow).
MAX_X_UNDIRECTED_DEGREE = 6  # bias_at costs 2**deg_undirected(X) per state visited;
# combined with up to 2**|K_G0| states in the largest shell this is the real
# driver of pathological per-instance cost, not |K_G0| alone -- gated separately.
N_VALUES = [6, 7, 8, 10, 12]
KNOWS_FRACTIONS = [0.4, 0.7, 1.0]
CORRUPTION_GRID: list[tuple[str, float]] = [
    ("none", 0.0),
    ("flip", 0.15),
    ("flip", 0.3),
    ("omit", 0.3),
]
GENERATOR_VARIANTS: list[tuple[str, dict[str, Any]]] = [
    ("erdos_renyi", {"edge_prob": 0.25}),
    ("erdos_renyi", {"edge_prob": 0.35}),
    ("erdos_renyi", {"edge_prob": 0.5}),
    ("scale_free", {"m_attach": 1}),
    ("scale_free", {"m_attach": 2}),
    ("block", {"n_blocks": 2, "p_within": 0.6, "p_between": 0.15}),
    ("decoupled_backdoor", {"coupling": 0.5}),
]

EPS_RELATIVE: tuple[float, ...] = (0.01, 0.02, 0.05, 0.10, 0.25, 0.50, 1.00)
R0_EPS_ABSOLUTE = 1e-12

STATES_SAVED_STRIDE = 8  # ~1/8 of accepted instances get the explicit shell-count audit
OLD_COMPARISON_CAP = 50  # instances (small enough for exact enumeration) for analysis 6
OLD_COMPARISON_MAX_UNDIRECTED = 5  # cap on |cpdag.undirected_edges| for build_space_correct
OLD_COMPARISON_MAX_SPACE = 150  # cap on |space| for the mean_abs_bias monotonicity check
N_EXAMPLES = 8  # describe(cert) written to examples.txt

INSTANCE_HARD_TIMEOUT_S = 8  # SIGALRM safety net per instance

CELLS: list[tuple[tuple[str, dict[str, Any]], int, float, tuple[str, float]]] = list(
    itertools.product(GENERATOR_VARIANTS, N_VALUES, KNOWS_FRACTIONS, CORRUPTION_GRID)
)


class _Timeout(Exception):  # noqa: N818 - a control-flow signal, not an error condition
    pass


def _alarm_handler(signum: int, frame: Any) -> None:
    raise _Timeout("instance exceeded hard wall-cap")


# --- row schema ---------------------------------------------------------------


def _row_template() -> dict[str, Any]:
    row: dict[str, Any] = {
        "instance_id": "",
        "generator": "",
        "generator_params_json": "",
        "n": 0,
        "seed": 0,
        "knows_fraction": 0.0,
        "corruption_name": "",
        "corruption_rate": 0.0,
        "n_undirected_cpdag": 0,
        "n_k0": 0,
        "x": "",
        "y": "",
        "z_size": 0,
        "z": "",
        "r_val": 0,
        "r_val_method": "",
        "r_val_seconds": 0.0,
        "r_val_exact": True,
        "theta_z": 0.0,
        "tau_true": 0.0,
        "beta_top": float("nan"),
        "r0_d": UNREACHED,
        "staircase_json": "[]",
        "traversal_limit_d": 0,
        "monotone_violation": False,
        "k_false": 0,
        "realised_error_effect": float("nan"),
        "beta_up_at_kfalse": float("nan"),
        "beta_up_at_kfalse_is_lower_bound": False,
        "certificate_holds": True,
        "conservativeness_ratio": float("nan"),
        "staircase_shape": "n/a",
        "staircase_n_distinct": 0,
        "staircase_total_rise": 0.0,
        "cost_incremental_seconds": float("nan"),
        "cost_incremental_states": 0,
        "cost_incremental_shells": 0,
        "cost_incremental_radius": UNREACHED,
        "cost_bisection_seconds": float("nan"),
        "cost_bisection_states": 0,
        "cost_bisection_shells": 0,
        "cost_bisection_radius": UNREACHED,
        "cost_disagreement": False,
        "cost_greedy_seconds": float("nan"),
        "cost_greedy_states": 0,
        "cost_greedy_radius": UNREACHED,
        "cost_greedy_tight": False,
        "states_saved_subsample": False,
        "states_saved_below_rval": -1,
        "old_reps_subsample": False,
        "old_r_eps_json": "",
        "old_new_disagree_count": -1,
        "old_new_disagree_direction": "",
        "old_mean_abs_bias_monotone": "",
        "status": "ok",
        "notes": "",
        "timing_total_seconds": 0.0,
    }
    for e in EPS_RELATIVE:
        row[f"r_eps_{e:g}"] = UNREACHED
    return row


# --- one instance -------------------------------------------------------------


def attempt_instance(
    index: int,
    cell: tuple[tuple[str, dict[str, Any]], int, float, tuple[str, float]],
    do_states_saved: bool,
    do_old_comparison: bool,
) -> tuple[str, dict[str, Any] | None, str | None]:
    """Build and certify one instance.

    Returns ``(status, row, reject_reason)``: ``status`` is ``"accepted"`` or
    ``"rejected"``; ``row`` is the CSV row (only when accepted); ``reject_reason``
    is the rejection code (only when rejected).
    """
    (generator_name, generator_params), n, knows_fraction, (corruption_name, corruption_rate) = cell
    seed = ROOT_SEED + index
    rng = np.random.default_rng(seed)

    if generator_name == "decoupled_backdoor" and n < 7:
        return "rejected", None, "generator_incompatible_n"

    dag = GENERATORS[generator_name](n, rng, **generator_params)
    cpdag = dag_to_cpdag(dag)
    x, y = _pick_treatment_outcome(generator_name, dag, rng)

    if len(cpdag.undirected_edges) > CAP_UNDIRECTED:
        return "rejected", None, "cpdag_too_large"
    if len(cpdag.neighbors(x)) > MAX_X_UNDIRECTED_DEGREE:
        return "rejected", None, "x_undirected_degree_too_high"

    accepted, reason = gate(dag, cpdag, x, y)
    if not accepted:
        return "rejected", None, reason

    k_true = draw_k_true(dag, cpdag, rng, knows_fraction)
    if corruption_name == "none":
        k_assumed = list(k_true)
        applied_rate = 0.0
    else:
        k_assumed = CORRUPTIONS[corruption_name](k_true, rng, corruption_rate)
        applied_rate = corruption_rate

    if not check_consistent(cpdag, k_assumed):
        return "rejected", None, "k_assumed_inconsistent"

    g0 = apply_orientations(cpdag, k_assumed)
    assert g0 is not None
    z = optimal_adjustment_set_mpdag(g0, x, y)
    if z is None:
        return "rejected", None, "z_not_identified"
    zf = frozenset(z)

    if not is_valid(zf, g0, x, y):
        return "rejected", None, "z_invalid_at_g0"

    # --- accepted: everything from here on is recorded, never silently dropped ---
    row = _row_template()
    instance_id = f"inst_{index:06d}"
    row.update(
        {
            "instance_id": instance_id,
            "generator": generator_name,
            "generator_params_json": json.dumps(generator_params, sort_keys=True),
            "n": n,
            "seed": seed,
            "knows_fraction": knows_fraction,
            "corruption_name": corruption_name,
            "corruption_rate": float(applied_rate),
            "n_undirected_cpdag": len(cpdag.undirected_edges),
            "x": x,
            "y": y,
            "z_size": len(zf),
            "z": "|".join(sorted(zf)),
        }
    )

    t_start = time.perf_counter()
    old_handler = signal.signal(signal.SIGALRM, _alarm_handler)
    signal.alarm(INSTANCE_HARD_TIMEOUT_S)
    try:
        k0_edges = knowledge_of(cpdag, g0)
        n_k0 = len(k0_edges)
        row["n_k0"] = n_k0
        search_budget = max(3, n_k0)

        sem = random_sem(dag, rng)

        cert = certify(
            cpdag,
            list(k_assumed),
            x,
            y,
            sem=sem,
            z=zf,
            g0=g0,
            epsilons=EPS_RELATIVE,
            units="relative",
            search_budget=search_budget,
        )

        row["r_val"] = cert.r_val
        row["r_val_method"] = cert.method
        row["r_val_seconds"] = cert.seconds_r_val
        row["r_val_exact"] = cert.status == "ok"
        row["theta_z"] = cert.theta_z
        r_val_int = cert.r_val if cert.r_val != UNREACHED else 0

        # theta_z ~ 0 makes "relative" undefined for the certificate proper, but
        # certify() would already have raised in that case (units="relative").
        # ctx rebuilt independently for the cost-benchmark calls below, which
        # need a BiasContext object rather than a Certificate.
        ctx = make_context(sem, cpdag, x, y, zf)
        row["tau_true"] = ctx.tau_true

        for e, radius in cert.r_eps.items():
            row[f"r_eps_{e:g}"] = radius

        staircase = [(s.d, s.beta_up) for s in cert.profile]
        row["staircase_json"] = json.dumps(staircase)
        row["traversal_limit_d"] = staircase[-1][0] if staircase else 0

        # Monotonicity sanity check (Theorem C.1) -- cheap, and a real bug would
        # be an important finding.
        beta_seq = [b for _, b in staircase]
        row["monotone_violation"] = any(
            beta_seq[i] > beta_seq[i + 1] + 1e-9 for i in range(len(beta_seq) - 1)
        )

        # beta_top: exact B(Chat), not read off a possibly-truncated traversal.
        top = top_bias(ctx)
        row["beta_top"] = top.worst

        # r_0: eps ~ 0 in absolute units, read straight off the same staircase
        # (Theorem C.4: one traversal answers every eps).
        row["r0_d"] = r_epsilon_from_profile(cert.profile, R0_EPS_ABSOLUTE)

        # --- certificate audit -------------------------------------------------
        true_dag_directed = dag.directed_edges
        k_false = sum(1 for a, b in k0_edges if (a, b) not in true_dag_directed)
        row["k_false"] = k_false
        realised_effect = abs(cert.theta_z - ctx.tau_true)
        row["realised_error_effect"] = realised_effect

        max_profiled_d = staircase[-1][0] if staircase else -1
        beta_up_kfalse = 0.0
        for d, b in staircase:
            if d <= k_false:
                beta_up_kfalse = max(beta_up_kfalse, b)
        row["beta_up_at_kfalse"] = beta_up_kfalse
        is_lower_bound = k_false > max_profiled_d
        row["beta_up_at_kfalse_is_lower_bound"] = is_lower_bound
        row["certificate_holds"] = bool(realised_effect <= beta_up_kfalse + 1e-9)
        row["conservativeness_ratio"] = (
            realised_effect / beta_up_kfalse if beta_up_kfalse > ZERO_TOL else float("nan")
        )

        # --- staircase shape -----------------------------------------------
        seq = beta_seq
        n_distinct = len(set(round(v, 12) for v in seq))
        increases = sum(1 for i in range(1, len(seq)) if seq[i] > seq[i - 1] + 1e-12)
        if len(seq) < 2 or increases == 0:
            shape = "flat"
        elif increases == 1:
            shape = "step"
        else:
            shape = "ramp"
        if r_val_int == UNREACHED or n_k0 == 0:
            shape = "n/a"
        row["staircase_shape"] = shape
        row["staircase_n_distinct"] = n_distinct
        row["staircase_total_rise"] = (seq[-1] - seq[0]) if seq else 0.0

        # --- cost: incremental vs bisection vs greedy, on the largest eps ---
        largest_eff = 1.00 * cert.eps_scale
        res_inc = r_epsilon(ctx, g0, largest_eff, r_val=r_val_int, strategy="incremental")
        res_bis = r_epsilon(ctx, g0, largest_eff, r_val=r_val_int, strategy="bisection")
        row["cost_incremental_seconds"] = res_inc.seconds
        row["cost_incremental_states"] = res_inc.states_evaluated
        row["cost_incremental_shells"] = res_inc.shells_evaluated
        row["cost_incremental_radius"] = res_inc.radius
        row["cost_bisection_seconds"] = res_bis.seconds
        row["cost_bisection_states"] = res_bis.states_evaluated
        row["cost_bisection_shells"] = res_bis.shells_evaluated
        row["cost_bisection_radius"] = res_bis.radius
        disagree = res_inc.radius != res_bis.radius
        row["cost_disagreement"] = disagree
        if disagree:
            (OUT_DIR / "hard_errors.jsonl").parent.mkdir(parents=True, exist_ok=True)
            with (OUT_DIR / "hard_errors.jsonl").open("a") as fh:
                fh.write(
                    json.dumps(
                        {
                            "instance_id": instance_id,
                            "kind": "incremental_bisection_disagreement",
                            "incremental_radius": res_inc.radius,
                            "bisection_radius": res_bis.radius,
                            "eps": largest_eff,
                            "g0": g0.edge_string(),
                            "cpdag": cpdag.edge_string(),
                        }
                    )
                    + "\n"
                )

        res_greedy = greedy_chain_bound(ctx, g0, largest_eff)
        row["cost_greedy_seconds"] = res_greedy.seconds
        row["cost_greedy_states"] = res_greedy.states_evaluated
        row["cost_greedy_radius"] = res_greedy.radius
        row["cost_greedy_tight"] = bool(res_greedy.radius == res_inc.radius)

        # --- states-not-evaluated saving (subsample) ------------------------
        if do_states_saved and r_val_int > 0:
            saved = sum(len(retraction_shell(cpdag, g0, d)) for d in range(r_val_int))
            row["states_saved_subsample"] = True
            row["states_saved_below_rval"] = saved
        elif do_states_saved:
            row["states_saved_subsample"] = True
            row["states_saved_below_rval"] = 0

        # --- comparison with the incumbent mean_abs_bias radius (subsample) -
        # Gated on |cpdag.undirected_edges| *before* doing any work: build_space_correct
        # is exponential in that count, so this must not run speculatively.
        if do_old_comparison and len(cpdag.undirected_edges) <= OLD_COMPARISON_MAX_UNDIRECTED:
            row["old_reps_subsample"] = True
            try:
                space = build_space_correct(cpdag)
                if g0 in space.neighbours:
                    dists = distances_from(space, g0)
                    abs_eps_list = [e * cert.eps_scale for e in EPS_RELATIVE]
                    old_rng = np.random.default_rng(seed + 9_000_000)
                    old_reps = old_radius_eps(
                        space, dists, zf, x, y, abs_eps_list, old_rng, n_bias_draws=25
                    )
                    row["old_r_eps_json"] = json.dumps(old_reps)
                    disagree_count = 0
                    directions = set()
                    for e_rel, e_abs in zip(EPS_RELATIVE, abs_eps_list, strict=True):
                        new_r = cert.r_eps.get(float(e_rel), UNREACHED)
                        old_r = old_reps.get(f"{e_abs:g}", UNREACHED)
                        if old_r != new_r:
                            disagree_count += 1
                            if old_r == UNREACHED:
                                directions.add("old_larger")
                            elif new_r == UNREACHED:
                                directions.add("new_larger")
                            elif old_r > new_r:
                                directions.add("old_larger")
                            else:
                                directions.add("new_larger")
                    row["old_new_disagree_count"] = disagree_count
                    row["old_new_disagree_direction"] = (
                        "equal"
                        if disagree_count == 0
                        else ("mixed" if len(directions) > 1 else next(iter(directions)))
                    )

                    # Is mean_abs_bias even monotone here? Check on covering
                    # pairs only: monotonicity there implies it on the whole
                    # finite poset by transitivity along chains.
                    if len(space.elements) <= OLD_COMPARISON_MAX_SPACE:
                        mono_rng = np.random.default_rng(seed + 9_100_000)
                        mean_bias_cache: dict[Any, float] = {}

                        def mean_bias_of(g: Any) -> float:
                            if g not in mean_bias_cache:
                                stats = bias_stats(zf, g, x, y, mono_rng, n_draws=20)
                                mean_bias_cache[g] = (
                                    stats.mean_abs_bias if stats.n_evaluations > 0 else 0.0
                                )
                            return mean_bias_cache[g]

                        violations = 0
                        n_pairs = 0
                        for lower, upper in space.covers:
                            n_pairs += 1
                            b_lo, b_hi = mean_bias_of(lower), mean_bias_of(upper)
                            if b_lo > b_hi + max(1e-6, 0.02 * max(b_lo, b_hi, 1e-9)):
                                violations += 1
                        row["old_mean_abs_bias_monotone"] = (
                            "monotone" if violations == 0 else f"violated_{violations}_{n_pairs}"
                        )
                    else:
                        row["old_mean_abs_bias_monotone"] = "skipped_large_space"
                else:
                    row["old_r_eps_json"] = json.dumps({"error": "g0_not_in_corrected_space"})
            except Exception as exc:  # never let the audit subsample kill an instance
                row["old_r_eps_json"] = json.dumps({"error": f"{type(exc).__name__}: {exc}"})

        row["status"] = "ok"
    except _Timeout:
        row["status"] = "capped_timeout"
        row["notes"] = f"hard wall-cap of {INSTANCE_HARD_TIMEOUT_S}s hit"
    except Exception as exc:  # an unexpected bug -- record, don't crash the sweep
        row["status"] = "run_error"
        row["notes"] = f"{type(exc).__name__}: {exc}"
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)
        row["timing_total_seconds"] = time.perf_counter() - t_start

    return "accepted", row, None


# --- main ---------------------------------------------------------------------


def run(target: int, time_budget_s: float) -> dict[str, Any]:
    """Sweep the grid, writing each accepted instance as it is produced."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUT_DIR / "results.csv"
    if (OUT_DIR / "hard_errors.jsonl").exists():
        (OUT_DIR / "hard_errors.jsonl").unlink()

    rejection_counts: dict[str, int] = {}
    status_counts: dict[str, int] = {}
    n_accepted = 0
    n_attempted = 0
    states_saved_count = 0
    old_comparison_count = 0
    examples: list[str] = []
    t0 = time.perf_counter()

    with ResultWriter(csv_path, resume=False) as writer:
        rep = 0
        while n_accepted < target and (time.perf_counter() - t0) < time_budget_s:
            for cell_idx, cell in enumerate(CELLS):
                if n_accepted >= target or (time.perf_counter() - t0) >= time_budget_s:
                    break
                index = rep * len(CELLS) + cell_idx
                n_attempted += 1

                do_states_saved = n_accepted % STATES_SAVED_STRIDE == 0 and states_saved_count < 250
                # cheap pre-check for the old-comparison eligibility happens
                # inside attempt_instance once |cpdag.undirected_edges| is known;
                # here we just cap the *count* so we don't try forever.
                do_old_comparison = old_comparison_count < OLD_COMPARISON_CAP

                status, row, reason = attempt_instance(
                    index, cell, do_states_saved, do_old_comparison
                )
                if status == "rejected":
                    rejection_counts[reason] = rejection_counts.get(reason, 0) + 1
                    continue

                assert row is not None
                if row["old_reps_subsample"]:
                    old_comparison_count += 1
                if do_states_saved:
                    states_saved_count += 1

                status_counts[row["status"]] = status_counts.get(row["status"], 0) + 1
                writer.write(row)
                n_accepted += 1
                if len(examples) < N_EXAMPLES and row["status"] == "ok":
                    examples.append(row["instance_id"])
                if n_accepted % 100 == 0:
                    elapsed = time.perf_counter() - t0
                    print(
                        f"[{elapsed:7.1f}s] accepted={n_accepted} attempted={n_attempted} "
                        f"(rejection rate "
                        f"{sum(rejection_counts.values()) / max(1, n_attempted):.2%})",
                        flush=True,
                    )
            rep += 1

    total_seconds = time.perf_counter() - t0
    return {
        "n_attempted": n_attempted,
        "n_accepted": n_accepted,
        "rejection_counts": dict(sorted(rejection_counts.items())),
        "status_counts": dict(sorted(status_counts.items())),
        "states_saved_count": states_saved_count,
        "old_comparison_count": old_comparison_count,
        "example_ids": examples,
        "total_seconds": total_seconds,
        "csv_path": str(csv_path),
    }


# --- post-hoc examples.txt (re-derives a handful of certificates for describe()) -----


def write_examples(example_ids: list[str], csv_path: Path) -> None:
    """Write a few full `describe(cert)` strings, so the prose form is inspectable."""
    import csv as csv_mod

    if not example_ids:
        return
    wanted = set(example_ids)
    lines: list[str] = []
    with csv_path.open(newline="") as fh:
        for row in csv_mod.DictReader(fh):
            if row["instance_id"] not in wanted:
                continue
            lines.append(f"=== {row['instance_id']} ===")
            lines.append(
                f"generator={row['generator']} n={row['n']} x={row['x']} y={row['y']} "
                f"z={row['z']} r_val={row['r_val']} theta_z={row['theta_z']}"
            )
            lines.append(f"staircase (d, beta_up effect units): {row['staircase_json']}")
            lines.append(
                "r_eps (relative): "
                + json.dumps({e: row.get(f"r_eps_{e:g}") for e in EPS_RELATIVE})
            )
            lines.append(
                f"k_false={row['k_false']} realised={row['realised_error_effect']} "
                f"beta_up(k_false)={row['beta_up_at_kfalse']} holds={row['certificate_holds']} "
                f"conservativeness={row['conservativeness_ratio']}"
            )
            lines.append("")
    (OUT_DIR / "examples.txt").write_text("\n".join(lines))


def _finite(series: pd.Series) -> pd.Series:
    """Drop UNREACHED (-1) sentinels and NaNs -- never average a status."""
    return series[(series != UNREACHED) & series.notna()]


def _describe(series: pd.Series) -> dict[str, float]:
    s = _finite(series.astype(float))
    if len(s) == 0:
        return {"n": 0}
    return {
        "n": len(s),
        "mean": float(s.mean()),
        "median": float(s.median()),
        "std": float(s.std()) if len(s) > 1 else 0.0,
        "min": float(s.min()),
        "max": float(s.max()),
        "p25": float(s.quantile(0.25)),
        "p75": float(s.quantile(0.75)),
    }


def write_analysis(csv_path: Path) -> dict[str, Any]:
    """Analyses 1-6 over the accepted-instance CSV. Writes analysis.json + SUMMARY.md."""
    df = pd.read_csv(csv_path)
    n = len(df)
    ok = df[df["status"] == "ok"]
    eps_cols = [(e, f"r_eps_{e:g}") for e in EPS_RELATIVE]

    analysis: dict[str, Any] = {"n_rows": n, "n_ok": len(ok)}

    # --- 1. distribution of r_val and each r_eps; fraction with r_eps > r_val ---
    a1: dict[str, Any] = {"r_val": _describe(ok["r_val"])}
    a1["r_val_unreached_fraction"] = float((ok["r_val"] == UNREACHED).mean()) if len(ok) else None
    per_eps = {}
    for e, col in eps_cols:
        sub = ok[col]
        unreached_frac = float((sub == UNREACHED).mean()) if len(ok) else None
        both_finite = ok[(ok["r_val"] != UNREACHED) & (sub != UNREACHED)]
        strictly_greater = both_finite[both_finite[col] > both_finite["r_val"]]
        # r_eps == UNREACHED while r_val is finite also "buys more" (infinitely so).
        r_eps_unreached_r_val_finite = ok[(ok["r_val"] != UNREACHED) & (sub == UNREACHED)]
        n_buys_more = len(strictly_greater) + len(r_eps_unreached_r_val_finite)
        n_comparable = len(ok[ok["r_val"] != UNREACHED])
        per_eps[e] = {
            "distribution": _describe(sub),
            "unreached_fraction": unreached_frac,
            "fraction_r_eps_gt_r_val": (n_buys_more / n_comparable) if n_comparable else None,
            "mean_extra_shells_when_both_finite": (
                float((strictly_greater[col] - strictly_greater["r_val"]).mean())
                if len(strictly_greater)
                else None
            ),
        }
    a1["per_eps"] = per_eps
    analysis["1_r_val_and_r_eps_distributions"] = a1

    # --- 2. r_0 == r_val (Remark 3's "generically equal" claim) ---
    comparable = ok.copy()
    equal_mask = comparable["r0_d"] == comparable["r_val"]
    n_comp = len(comparable)
    exceptions = comparable[~equal_mask]
    witnesses = exceptions.head(20)[
        ["instance_id", "r_val", "r0_d", "generator", "n", "n_k0", "staircase_shape"]
    ].to_dict(orient="records")
    analysis["2_r0_equals_r_val"] = {
        "n": n_comp,
        "fraction_equal": float(equal_mask.mean()) if n_comp else None,
        "n_exceptions": int((~equal_mask).sum()),
        "exception_witnesses": witnesses,
    }

    # --- 3. staircase shape distribution ---
    shape_counts = ok["staircase_shape"].value_counts(dropna=False).to_dict()
    non_na = ok[ok["staircase_shape"] != "n/a"]
    analysis["3_staircase_shape"] = {
        "counts": {str(k): int(v) for k, v in shape_counts.items()},
        "fraction_ramp_among_nontrivial": (
            float((non_na["staircase_shape"] == "ramp").mean()) if len(non_na) else None
        ),
        "fraction_step_among_nontrivial": (
            float((non_na["staircase_shape"] == "step").mean()) if len(non_na) else None
        ),
        "fraction_flat_among_nontrivial": (
            float((non_na["staircase_shape"] == "flat").mean()) if len(non_na) else None
        ),
        "mean_n_distinct": float(ok["staircase_n_distinct"].mean()) if len(ok) else None,
        "mean_total_rise": float(ok["staircase_total_rise"].mean()) if len(ok) else None,
    }

    # --- 4. certificate conservativeness ---
    #
    # The comparison is only meaningful where beta_up(k_false) was actually
    # COMPUTED. The traversal stops as soon as it crosses the largest threshold
    # on the grid, so when k_false lies beyond the last evaluated shell the
    # recorded beta_up(k_false) is a LOWER bound on the real one -- and a
    # realised error above a lower bound is not a violation of anything. Those
    # rows are separated out rather than counted as failures or silently
    # dropped: three of them were reported as certificate violations once, and
    # recomputing their full staircases showed the certificate holding in all
    # three, with equality (which is Theorem D-prime, tightness, showing up).
    comparable = ok[~ok["beta_up_at_kfalse_is_lower_bound"].astype(bool)]
    truncated = ok[ok["beta_up_at_kfalse_is_lower_bound"].astype(bool)]
    cons = comparable["conservativeness_ratio"]
    holds_violations = comparable[~comparable["certificate_holds"].astype(bool)]
    truncated_over = truncated[~truncated["certificate_holds"].astype(bool)]
    analysis["4_conservativeness"] = {
        "distribution": _describe(cons),
        "n_comparable": len(comparable),
        "n_truncated_excluded": len(truncated),
        "certificate_always_holds": bool(len(holds_violations) == 0),
        "n_violations": len(holds_violations),
        "violation_witnesses": holds_violations.head(10)[
            ["instance_id", "realised_error_effect", "beta_up_at_kfalse", "k_false"]
        ].to_dict(orient="records"),
        "n_exceeding_a_truncated_lower_bound": len(truncated_over),
        "truncated_note": (
            "Rows where k_false exceeded the last evaluated shell, so the recorded "
            "beta_up(k_false) is a lower bound and the comparison is not a test of "
            "Theorem D. Excluded from the violation count and from the ratio."
        ),
        "monotone_violations": int(ok["monotone_violation"].sum()) if len(ok) else 0,
    }

    # --- 5. cost: incremental vs bisection vs greedy, and the shell-skip saving ---
    cost_cols = [
        "n_k0",
        "cost_incremental_seconds",
        "cost_incremental_states",
        "cost_bisection_seconds",
        "cost_bisection_states",
        "cost_greedy_seconds",
        "cost_greedy_states",
        "cost_greedy_tight",
        "cost_disagreement",
    ]
    cost_df = ok[cost_cols].copy()
    saved = ok[ok["states_saved_subsample"]]
    analysis["5_cost"] = {
        "incremental_states": _describe(cost_df["cost_incremental_states"]),
        "bisection_states": _describe(cost_df["cost_bisection_states"]),
        "incremental_seconds": _describe(cost_df["cost_incremental_seconds"]),
        "bisection_seconds": _describe(cost_df["cost_bisection_seconds"]),
        "incremental_cheaper_states_fraction": float(
            (cost_df["cost_incremental_states"] < cost_df["cost_bisection_states"]).mean()
        )
        if len(cost_df)
        else None,
        "greedy_tight_fraction": float(cost_df["cost_greedy_tight"].mean())
        if len(cost_df)
        else None,
        "greedy_seconds": _describe(cost_df["cost_greedy_seconds"]),
        "n_hard_disagreements": int(cost_df["cost_disagreement"].sum()),
        "states_by_n_k0": {
            str(k): {
                "incremental": float(v["cost_incremental_states"].mean()),
                "bisection": float(v["cost_bisection_states"].mean()),
                "n": len(v),
            }
            for k, v in cost_df.groupby("n_k0")
        },
        "shell_skip_saving": {
            "n_subsample": len(saved),
            "states_saved_below_r_val": _describe(saved["states_saved_below_rval"]),
        },
    }

    # --- 6. comparison with the incumbent mean_abs_bias radius ---
    old_sub = ok[ok["old_reps_subsample"]]
    direction_counts = old_sub["old_new_disagree_direction"].value_counts(dropna=False).to_dict()
    monotone_counts = old_sub["old_mean_abs_bias_monotone"].value_counts(dropna=False).to_dict()
    analysis["6_incumbent_comparison"] = {
        "n_subsample": len(old_sub),
        "fraction_any_disagreement": (
            float((old_sub["old_new_disagree_count"] > 0).mean()) if len(old_sub) else None
        ),
        "mean_n_thresholds_disagreeing_of_7": (
            float(old_sub["old_new_disagree_count"].mean()) if len(old_sub) else None
        ),
        "direction_counts": {str(k): int(v) for k, v in direction_counts.items()},
        "old_mean_abs_bias_monotone_counts": {str(k): int(v) for k, v in monotone_counts.items()},
    }

    (OUT_DIR / "analysis.json").write_text(json.dumps(analysis, indent=2, default=str))
    write_summary_md(analysis)
    return analysis


def write_summary_md(a: dict[str, Any]) -> None:
    """Render the analyses as the human-readable summary."""
    lines: list[str] = ["# r_epsilon empirical validation -- summary", ""]
    lines.append(f"Accepted instances analysed: {a['n_ok']} / {a['n_rows']} rows in results.csv.")
    lines.append("")

    lines.append("## 1. r_val and r_eps distributions")
    a1 = a["1_r_val_and_r_eps_distributions"]
    lines.append(f"- r_val: {a1['r_val']} (UNREACHED fraction {a1['r_val_unreached_fraction']})")
    for e, d in a1["per_eps"].items():
        lines.append(
            f"- eps={e}: r_eps {d['distribution']}, UNREACHED fraction {d['unreached_fraction']}, "
            f"fraction(r_eps > r_val)={d['fraction_r_eps_gt_r_val']}, "
            f"mean extra shells (both finite)={d['mean_extra_shells_when_both_finite']}"
        )
    lines.append("")

    lines.append("## 2. r_0 == r_val (Remark 3)")
    a2 = a["2_r0_equals_r_val"]
    lines.append(
        f"- n={a2['n']}, fraction equal={a2['fraction_equal']}, exceptions={a2['n_exceptions']}"
    )
    if a2["exception_witnesses"]:
        lines.append("- example exceptions:")
        for w in a2["exception_witnesses"][:5]:
            lines.append(f"    {w}")
    lines.append("")

    lines.append("## 3. staircase shape")
    a3 = a["3_staircase_shape"]
    lines.append(f"- counts: {a3['counts']}")
    lines.append(
        f"- among non-trivial: ramp={a3['fraction_ramp_among_nontrivial']}, "
        f"step={a3['fraction_step_among_nontrivial']}, flat={a3['fraction_flat_among_nontrivial']}"
    )
    lines.append("")

    lines.append("## 4. certificate conservativeness")
    a4 = a["4_conservativeness"]
    lines.append(f"- realised/beta_up(k_false): {a4['distribution']}")
    lines.append(
        f"- comparable rows: {a4['n_comparable']}; "
        f"excluded because beta_up(k_false) was only a lower bound: "
        f"{a4['n_truncated_excluded']} (of which {a4['n_exceeding_a_truncated_lower_bound']} "
        "had a realised error above that lower bound, which is not a violation -- see the note)"
    )
    lines.append(
        f"- certificate always holds: {a4['certificate_always_holds']} "
        f"(violations: {a4['n_violations']}); monotonicity violations: {a4['monotone_violations']}"
    )
    lines.append("")

    lines.append("## 5. cost")
    a5 = a["5_cost"]
    lines.append(f"- incremental states: {a5['incremental_states']}")
    lines.append(f"- bisection states: {a5['bisection_states']}")
    lines.append(
        f"- incremental cheaper (states) fraction: {a5['incremental_cheaper_states_fraction']}"
    )
    lines.append(f"- greedy_chain_bound tight fraction: {a5['greedy_tight_fraction']}")
    lines.append(f"- hard incremental/bisection disagreements: {a5['n_hard_disagreements']}")
    lines.append(
        f"- states saved below r_val (subsample n={a5['shell_skip_saving']['n_subsample']}): "
        f"{a5['shell_skip_saving']['states_saved_below_r_val']}"
    )
    lines.append("")

    lines.append("## 6. comparison with the incumbent mean_abs_bias radius")
    a6 = a["6_incumbent_comparison"]
    lines.append(f"- subsample n={a6['n_subsample']}")
    lines.append(f"- fraction with any disagreement: {a6['fraction_any_disagreement']}")
    lines.append(
        f"- mean # thresholds disagreeing (of 7): {a6['mean_n_thresholds_disagreeing_of_7']}"
    )
    lines.append(f"- direction counts: {a6['direction_counts']}")
    lines.append(f"- old mean_abs_bias monotone counts: {a6['old_mean_abs_bias_monotone_counts']}")
    lines.append("")

    (OUT_DIR / "SUMMARY.md").write_text("\n".join(lines))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", type=int, default=1200)
    ap.add_argument("--time-budget", type=float, default=2100.0, help="seconds for the grid sweep")
    args = ap.parse_args()

    print(f"git sha: {git_sha(short=False)}", flush=True)
    print(f"python: {sys.version.split()[0]}", flush=True)
    print(f"grid cells: {len(CELLS)}, target accepted: {args.target}", flush=True)

    summary = run(args.target, args.time_budget)
    print(json.dumps(summary, indent=2), flush=True)

    write_examples(summary["example_ids"], OUT_DIR / "results.csv")

    write_manifest(
        OUT_DIR,
        seed=ROOT_SEED,
        grid={
            "n_cells": len(CELLS),
            "generators": [g for g, _ in GENERATOR_VARIANTS],
            "n_values": N_VALUES,
            "knows_fractions": KNOWS_FRACTIONS,
            "corruptions": CORRUPTION_GRID,
            "cap_undirected": CAP_UNDIRECTED,
            "eps_relative_grid": EPS_RELATIVE,
            "target_accepted": args.target,
            "time_budget_s": args.time_budget,
        },
        extra={
            "n_attempted": summary["n_attempted"],
            "n_accepted": summary["n_accepted"],
            "rejection_counts": summary["rejection_counts"],
            "status_counts": summary["status_counts"],
            "states_saved_count": summary["states_saved_count"],
            "old_comparison_count": summary["old_comparison_count"],
            "total_seconds": summary["total_seconds"],
        },
    )
    # The sweep and its analysis are one command: a results.csv on disk with no
    # analysis.json beside it is a half-finished run, and the two drifting apart
    # is exactly how a stale number reaches a paper.
    analysis = write_analysis(OUT_DIR / "results.csv")
    print(
        json.dumps(
            {k: v for k, v in analysis.items() if k != "per_instance"}, indent=2, default=str
        )[:4000],
        flush=True,
    )
    print("done.", flush=True)
