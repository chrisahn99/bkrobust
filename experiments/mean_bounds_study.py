r"""Cheap statistics that bracket ``mu(d)`` and the mean-based radius ``r_mean``.

``mu(d)`` (``src/bkrobust/epsilon/meanprofile.py``) is the subset-weighted mean
bias over the retraction shell at depth ``d`` -- it varies with ``d`` where
``beta_up`` does not, but exact enumeration of a shell costs ``C(|K_G0|, d)``
Meek closures, which is expensive for the deeper shells of a large instance.
This script asks: what CHEAP statistics bracket ``mu(d)``, and therefore
bracket ``r_mean(eps)`` (the first ``d`` with ``mu(d) > eps``), well enough to
be useful when exhaustive enumeration is not affordable?

It evaluates, against exhaustive ground truth on every shell of every
instance:

1. The contamination sandwich already in ``meanprofile.mean_bounds``
   (``p * min+ <= mu <= p * beta_up``, ``p`` = contaminated fraction). This
   needs full shell enumeration to get ``p``, ``min+`` and the shell's own
   ``beta_up`` -- it is evaluated here as a reference, not as a cheap bound.
   The UPPER sandwich bound gives the SMALLER radius estimate and the LOWER
   sandwich bound gives the LARGER one (a bound that is too small on ``mu``
   crosses ``eps`` later, i.e. at a larger ``d``); this script gets that
   direction right and reports the resulting bracket width in shells.

2. A bound that needs NO shell enumeration: the running ``beta_up(d)`` that
   ``certify``/the staircase already produces (the max over the *ball* of
   radius ``d``, not the shell) divided by ``C(|K_G0|, d)``. This is only
   PROVABLY a valid lower bound on ``mu(d)`` at ``d = r_val`` (where the ball
   and the shell coincide, since everything below ``r_val`` is exactly zero)
   and at ``d = |K_G0|`` (where the shell has a single state, equal to
   ``beta_top``, so ``mu`` is known exactly for free). For intermediate
   ``d`` the ball's running max can exceed the shell's own true max, so the
   "bound" can be unsound; this script tests validity shell by shell rather
   than assuming it.

3. Whether the contaminated fraction ``p(d)`` is itself close to ``d / n``
   (``n = |K_G0|``), and how well the resulting near-free estimator
   ``mu_hat(d) = (d/n) * beta_up(d)`` (using the certify staircase's
   ``beta_up``, not the shell's true max) predicts both ``mu(d)`` and the
   radius read off it.

4. Whether ``mu`` is monotone enough, instance by instance, that evaluating
   only the shells ``r_val`` and ``|K_G0|`` brackets every intermediate shell
   for free (no interpolation, just the two endpoints as bounds).

Every bias number here is conditional on the seeded linear-Gaussian SEM
attached to the ground-truth DAG (the corpus ships no data); this is stated
once here and must travel with any number quoted from the outputs.

Usage::

    PYTHONPATH=src .venv/bin/python experiments/mean_bounds_study.py [options]

Writes ``<out>/shell_bounds.jsonl`` (one record per instance/shell),
``<out>/radius_bounds.jsonl`` (one record per instance/eps), and
``<out>/summary.json`` (counts and fractions aggregated over instances --
never per-network medians, per house rules).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.demo.evaluate import (  # noqa: E402
    optimal_adjustment_set_mpdag,
    random_sem,
)
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.epsilon.bias import knowledge_of, make_context  # noqa: E402
from bkrobust.epsilon.meanprofile import (  # noqa: E402
    ShellStat,
    is_monotone,
    mean_bounds,
    r_mean,
    shell_stat,
)

import final_table as ft  # noqa: E402

#: Same 15 instances as ``mean_vs_max_shell.py`` -- exhaustive ground truth
#: for these already exists in ``results/mean_vs_max/shells.jsonl`` and
#: ``results/final_table_eps_gt1/instances.jsonl`` (for ``r_val``,
#: ``beta_top`` and the certify staircase). Spans ``|K_G0|`` 3..14, r_val
#: 1..11.
DEFAULT_INSTANCES = (
    "asia:bronc:dysp",
    "asia:asia:either",
    "mediator:I:Y",
    "Sebastiani_2005:ANXA2.5:ANXA2.11",
    "Didelez_2010:Age:HRT",
    "magic-niab:G1217:YLD",
    "water:CKNI_12_15:CBODN_12_45",
    "Schipf_2010:A:TT",
    "barley:dg25:s2225",
    "magic-irri:G3212:FT",
    "child:CO2:DuctFlow",
    "Kampen_2014:AFF:SAN",
    "andes:HORIZ53:SNode_118",
    "ecoli70:cspG:hupB",
    "paths:11:14",
)

#: Same relative-units epsilon grid as ``results/final_table_eps_gt1``, so
#: the mean-based radius here is directly comparable to the existing
#: max-based ``r_eps``.
EPS_GRID: tuple[float, ...] = (1.0, 2.0, 5.0, 10.0, 20.0)

TOL = 1e-9


def load_ground_truth(path: Path) -> dict[tuple[str, str, str], dict[str, Any]]:
    """Index ``final_table_eps_gt1/instances.jsonl`` by ``(network, x, y)``.

    Args:
        path: Path to the jsonl file.

    Returns:
        Mapping to the parsed record, which carries ``r_val``, ``beta_top``
        and the ``staircase`` (``{str(d): beta_up(d)}``, the running max over
        the *ball* of radius ``d``, not the shell).
    """
    out: dict[tuple[str, str, str], dict[str, Any]] = {}
    with path.open() as f:
        for line in f:
            r = json.loads(line)
            out[(r["network"], r["x"], r["y"])] = r
    return out


def study_instance(
    net: str,
    x: str,
    y: str,
    cpdag,
    dag,
    k,
    seed: int,
    max_subsets: int,
    gt: dict[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Enumerate every shell exhaustively and score the candidate bounds on it.

    Args:
        net: Network name.
        x: Treatment.
        y: Outcome.
        cpdag: The estimated CPDAG.
        dag: The ground-truth DAG, for the seeded SEM.
        k: The elicited knowledge.
        seed: Base seed; the per-instance seed is derived as in final_table.
        max_subsets: Refuse instances needing more than ``2**n > max_subsets``.
        gt: This instance's record from ``results/final_table_eps_gt1``,
            supplying ``r_val``, ``beta_top`` and the certify staircase.

    Returns:
        ``(meta, shell_records)``. ``shell_records`` is empty unless
        ``meta["status"] == "ok"``.
    """
    t0 = time.perf_counter()
    g0 = apply_orientations(cpdag, list(k))
    if g0 is None:
        return {"network": net, "x": x, "y": y, "status": "inconsistent"}, []
    z_set = optimal_adjustment_set_mpdag(g0, x, y)
    if z_set is None:
        return {"network": net, "x": x, "y": y, "status": "degenerate"}, []

    rng = np.random.default_rng(ft.derive_seed(seed, net, x, y))
    sem = random_sem(dag, rng)
    ctx = make_context(sem, cpdag, x, y, frozenset(z_set))
    scale = abs(ctx.theta_z)
    if scale < 1e-12:
        return {"network": net, "x": x, "y": y, "status": "zero_estimate"}, []

    k0 = knowledge_of(cpdag, g0)
    n = len(k0)
    total = 2**n
    if total > max_subsets:
        return (
            {"network": net, "x": x, "y": y, "status": "too_large", "n_knowledge": n},
            [],
        )

    r_val = int(gt["r_val"])
    if r_val <= 0:
        # House rule: r_val == 0 is degenerate (adjustment set already
        # invalid) and must be excluded from everything.
        return {"network": net, "x": x, "y": y, "status": "degenerate_r_val0"}, []
    beta_top = float(gt["beta_top"])
    staircase = {int(dd): float(v) for dd, v in gt["staircase"].items()}

    cache: dict[str, float] = {}
    stats: dict[int, ShellStat] = {}
    for d in range(n + 1):
        stats[d] = shell_stat(ctx, cpdag, g0, d, scale=scale, cache=cache)

    records: list[dict[str, Any]] = []
    for d in range(n + 1):
        stat = stats[d]
        n_subsets = math.comb(n, d)
        lower_sandwich, upper_sandwich = mean_bounds(stat)

        certify_beta_up = staircase.get(d)
        free_lower_bound = (
            certify_beta_up / n_subsets if certify_beta_up is not None else None
        )
        free_lower_valid = (
            None
            if free_lower_bound is None
            else bool(free_lower_bound <= stat.mean + TOL)
        )
        shell_true_eq_certify_ball = (
            None
            if certify_beta_up is None
            else bool(abs(stat.beta_up - certify_beta_up) <= TOL)
        )

        p_hat_dn = d / n if n > 0 else 0.0
        mu_hat_free = (
            p_hat_dn * certify_beta_up if certify_beta_up is not None else None
        )
        mu_hat_oracle_p = p_hat_dn * stat.beta_up

        records.append(
            {
                "network": net,
                "x": x,
                "y": y,
                "d": d,
                "r_val": r_val,
                "n_knowledge": n,
                "n_subsets_total": n_subsets,
                "mu_true": stat.mean,
                "beta_up_shell_true": stat.beta_up,
                "beta_up_certify_ball": certify_beta_up,
                "shell_true_eq_certify_ball": shell_true_eq_certify_ball,
                "min_nonzero_true": stat.min_nonzero,
                "frac_nonzero_true": stat.frac_nonzero,
                "p_hat_d_over_n": p_hat_dn,
                "p_hat_abs_err": abs(p_hat_dn - stat.frac_nonzero),
                "lower_sandwich": lower_sandwich,
                "upper_sandwich": upper_sandwich,
                "sandwich_upper_exact": bool(
                    abs(upper_sandwich - stat.mean) <= TOL
                ),
                "sandwich_violated": bool(
                    lower_sandwich > stat.mean + TOL
                    or upper_sandwich < stat.mean - TOL
                ),
                "free_lower_bound": free_lower_bound,
                "free_lower_bound_valid": free_lower_valid,
                "mu_hat_free": mu_hat_free,
                "mu_hat_free_abs_err": (
                    None if mu_hat_free is None else abs(mu_hat_free - stat.mean)
                ),
                "mu_hat_oracle_p": mu_hat_oracle_p,
                "mu_hat_oracle_p_abs_err": abs(mu_hat_oracle_p - stat.mean),
                "below_r_val_zero_ok": bool(d >= r_val or stat.mean == 0.0),
            }
        )

    meta = {
        "network": net,
        "x": x,
        "y": y,
        "status": "ok",
        "n_knowledge": n,
        "r_val": r_val,
        "beta_top": beta_top,
        "mu_at_n_matches_beta_top": bool(abs(stats[n].mean - beta_top) <= TOL),
        "beta_up_shell_at_rval_matches_certify": bool(
            abs(stats[r_val].beta_up - staircase.get(r_val, float("nan"))) <= TOL
        ),
        "is_monotone_from_rval": is_monotone(
            [stats[d] for d in range(r_val, n + 1)]
        ),
        "seconds": round(time.perf_counter() - t0, 3),
    }
    return meta, records


def radius_records(
    net: str, x: str, y: str, meta: dict[str, Any], shells: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    """Build the per-``eps`` radius comparison for one instance.

    Compares ``r_mean`` read off the true profile against the same read off
    the sandwich bounds, the near-free estimator, and a monotone two-point
    bracket using only the shells ``r_val`` and ``n``.

    Args:
        net: Network name.
        x: Treatment.
        y: Outcome.
        meta: This instance's meta record from :func:`study_instance`.
        shells: This instance's shell records from :func:`study_instance`,
            one per ``d`` in ``0..n``.

    Returns:
        One record per ``eps`` in :data:`EPS_GRID`.
    """
    n = meta["n_knowledge"]
    r_val = meta["r_val"]
    by_d = {r["d"]: r for r in shells}

    def series(key: str) -> list[tuple[int, float]]:
        # Only d >= r_val is meaningful (below is exactly zero, by
        # construction of shell_stat + Theorem A); r_mean() only needs
        # increasing d and tolerates the profile starting at r_val.
        return [(d, by_d[d][key]) for d in range(r_val, n + 1)]

    def r_of(vals: list[tuple[int, float]], eps: float) -> int:
        for d, v in vals:
            if v > eps:
                return d
        return UNREACHED

    mu_true = series("mu_true")
    lower = series("lower_sandwich")
    upper = series("upper_sandwich")
    mu_hat_free_series = [
        (d, by_d[d]["mu_hat_free"])
        for d in range(r_val, n + 1)
        if by_d[d]["mu_hat_free"] is not None
    ]

    # Monotone two-point bracket: only shells r_val and n are "evaluated";
    # every d strictly between them is bracketed, under monotonicity, by
    # mu(r_val) <= mu(d) <= mu(n). r_mean read off the constant-mu(r_val)
    # floor gives the latest possible radius consistent with that bracket,
    # and off the constant-mu(n) ceiling gives the earliest.
    mu_rval = by_d[r_val]["mu_true"]
    mu_n = by_d[n]["mu_true"]

    out = []
    for eps in EPS_GRID:
        r_true = r_of(mu_true, eps)
        r_lower = r_of(lower, eps)
        r_upper = r_of(upper, eps)
        r_free = r_of(mu_hat_free_series, eps) if mu_hat_free_series else UNREACHED

        def bracket_width(lo_r: int, hi_r: int) -> int | None:
            if lo_r == UNREACHED or hi_r == UNREACHED:
                return None
            return abs(hi_r - lo_r)

        # Monotone bracket on the radius itself: if mu(r_val) > eps the
        # radius is r_val outright; if mu(n) <= eps it is UNREACHED
        # outright; otherwise the true crossing lies somewhere in
        # (r_val, n], which is the bracket width reported.
        if mu_rval > eps:
            r_mono_lo = r_mono_hi = r_val
        elif mu_n <= eps:
            r_mono_lo = r_mono_hi = UNREACHED
        else:
            r_mono_lo, r_mono_hi = r_val + 1, n

        out.append(
            {
                "network": net,
                "x": x,
                "y": y,
                "eps": eps,
                "r_val": r_val,
                "n_knowledge": n,
                "r_true": r_true,
                "r_from_upper_sandwich": r_upper,
                "r_from_lower_sandwich": r_lower,
                "sandwich_bracket_width": bracket_width(r_upper, r_lower),
                "r_from_mu_hat_free": r_free,
                "mu_hat_free_err": (
                    None
                    if r_free == UNREACHED or r_true == UNREACHED
                    else r_free - r_true
                ),
                "r_monotone_bracket_lo": r_mono_lo,
                "r_monotone_bracket_hi": r_mono_hi,
                "monotone_bracket_width": bracket_width(r_mono_lo, r_mono_hi),
                "monotone_bracket_contains_r_true": (
                    None
                    if r_true == UNREACHED and r_mono_lo == UNREACHED
                    else bool(
                        (r_mono_lo != UNREACHED and r_true != UNREACHED
                         and r_mono_lo <= r_true <= r_mono_hi)
                        or (r_mono_lo == UNREACHED and r_true == UNREACHED)
                    )
                ),
            }
        )
    return out


def main() -> None:
    """Run the bounding-statistics study and write its three output files."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--condition", default="D_LLM")
    p.add_argument("--instances", default=",".join(DEFAULT_INSTANCES))
    p.add_argument("--seed", type=int, default=20260917)
    p.add_argument(
        "--max-subsets",
        type=int,
        default=2**15,
        help="skip an instance needing more than this many subsets",
    )
    p.add_argument(
        "--ground-truth",
        default="results/final_table_eps_gt1/instances.jsonl",
        help="source of r_val, beta_top and the certify staircase",
    )
    p.add_argument("--out", default="results/mean_bounds")
    args = p.parse_args()

    wanted = [tuple(s.split(":", 2)) for s in args.instances.split(",") if s]
    nets_wanted = {w[0] for w in wanted}
    knowledge = ft.load_knowledge(args.condition)
    parsed = ft.load_networks(names=nets_wanted, skip=set(), max_nodes=0)
    gt = load_ground_truth(ROOT / args.ground_truth)

    out_dir = ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    metas: list[dict[str, Any]] = []
    shell_records: list[dict[str, Any]] = []
    radius_out: list[dict[str, Any]] = []

    for net, x, y in wanted:
        if net not in parsed or net not in knowledge:
            print(f"[mean_bounds] {net}: unavailable", file=sys.stderr)
            continue
        key = (net, x, y)
        if key not in gt:
            print(f"[mean_bounds] {net} {x}->{y}: no ground truth row, skipping",
                  file=sys.stderr)
            continue
        meta, recs = study_instance(
            net, x, y, parsed[net]["cpdag"], parsed[net]["dag"],
            knowledge[net], args.seed, args.max_subsets, gt[key],
        )
        metas.append(meta)
        shell_records.extend(recs)
        print(
            f"[mean_bounds] {net} {x}->{y} status={meta['status']} "
            f"|K_G0|={meta.get('n_knowledge')} r_val={meta.get('r_val')} "
            f"seconds={meta.get('seconds')}",
            flush=True,
        )
        if meta["status"] == "ok":
            radius_out.extend(radius_records(net, x, y, meta, recs))

    with (out_dir / "shell_bounds.jsonl").open("w") as f:
        for r in shell_records:
            f.write(json.dumps(r) + "\n")
    with (out_dir / "radius_bounds.jsonl").open("w") as f:
        for r in radius_out:
            f.write(json.dumps(r) + "\n")

    ok_metas = [m for m in metas if m["status"] == "ok"]
    ok_shells = [r for r in shell_records if True]
    below_rval = [r for r in ok_shells if r["d"] < r["r_val"]]
    at_or_above = [r for r in ok_shells if r["d"] >= r["r_val"]]

    def frac(xs: list[bool]) -> float | None:
        return (sum(1 for v in xs if v) / len(xs)) if xs else None

    summary = {
        "args": vars(args),
        "note": (
            "All bias numbers are conditional on the seeded linear-Gaussian "
            "SEM attached to each instance's ground-truth DAG (the corpus "
            "ships no data). Condition D_LLM only. Scope: the "
            f"{len(ok_metas)} status=ok instances actually run here."
        ),
        "n_instances_requested": len(wanted),
        "n_instances_ok": len(ok_metas),
        "instance_statuses": {
            m["network"] + ":" + m["x"] + ":" + m["y"]: m["status"] for m in metas
        },
        "fact_a_mu_zero_below_rval": {
            "n_shells_checked": len(below_rval),
            "n_violations": sum(1 for r in below_rval if r["mu_true"] != 0.0),
        },
        "fact_b_sandwich_never_violated": {
            "n_shells_checked": len(ok_shells),
            "n_violations": sum(1 for r in ok_shells if r["sandwich_violated"]),
            "frac_upper_exact": frac([r["sandwich_upper_exact"] for r in ok_shells]),
            "frac_upper_exact_at_or_above_rval": frac(
                [r["sandwich_upper_exact"] for r in at_or_above]
            ),
        },
        "mu_at_n_equals_beta_top": {
            "n_instances_checked": len(ok_metas),
            "n_violations": sum(
                1 for m in ok_metas if not m["mu_at_n_matches_beta_top"]
            ),
        },
        "beta_up_shell_eq_certify_ball_at_rval": {
            "n_instances_checked": len(ok_metas),
            "n_violations": sum(
                1 for m in ok_metas
                if not m["beta_up_shell_at_rval_matches_certify"]
            ),
        },
        "shell_true_eq_certify_ball_above_rval": {
            "n_shells_checked": len([r for r in ok_shells if r["d"] > r["r_val"]]),
            "n_matching": sum(
                1
                for r in ok_shells
                if r["d"] > r["r_val"] and r["shell_true_eq_certify_ball"]
            ),
        },
        "free_lower_bound_tightness": (
            lambda ratios: {
                "n": len(ratios),
                "min_ratio_mu_over_bound": min(ratios) if ratios else None,
                "median_ratio_mu_over_bound": (
                    sorted(ratios)[len(ratios) // 2] if ratios else None
                ),
                "max_ratio_mu_over_bound": max(ratios) if ratios else None,
            }
        )(
            [
                r["mu_true"] / r["free_lower_bound"]
                for r in ok_shells
                if r["d"] >= r["r_val"]
                and r["free_lower_bound"]
                and r["free_lower_bound"] > 0
            ]
        ),
        "free_lower_bound_validity": {
            "at_r_val": {
                "n": len([r for r in ok_shells if r["d"] == r["r_val"]]),
                "frac_valid": frac(
                    [
                        r["free_lower_bound_valid"]
                        for r in ok_shells
                        if r["d"] == r["r_val"]
                        and r["free_lower_bound_valid"] is not None
                    ]
                ),
            },
            "above_r_val": {
                "n": len([r for r in ok_shells if r["d"] > r["r_val"]]),
                "frac_valid": frac(
                    [
                        r["free_lower_bound_valid"]
                        for r in ok_shells
                        if r["d"] > r["r_val"]
                        and r["free_lower_bound_valid"] is not None
                    ]
                ),
            },
        },
        "p_hat_d_over_n_hypothesis": {
            "n_shells_checked": len(at_or_above),
            "mean_abs_err": (
                sum(r["p_hat_abs_err"] for r in at_or_above) / len(at_or_above)
                if at_or_above else None
            ),
            "max_abs_err": (
                max((r["p_hat_abs_err"] for r in at_or_above), default=None)
            ),
            "n_within_0p05": sum(
                1 for r in at_or_above if r["p_hat_abs_err"] <= 0.05
            ),
            "n_exact": sum(
                1 for r in at_or_above if r["p_hat_abs_err"] <= TOL
            ),
        },
        "mu_hat_free_estimator": {
            "n_shells_checked": len(
                [r for r in at_or_above if r["mu_hat_free"] is not None]
            ),
            "mean_abs_err": (
                sum(
                    r["mu_hat_free_abs_err"]
                    for r in at_or_above
                    if r["mu_hat_free_abs_err"] is not None
                )
                / max(
                    1,
                    len(
                        [
                            r
                            for r in at_or_above
                            if r["mu_hat_free_abs_err"] is not None
                        ]
                    ),
                )
                if any(r["mu_hat_free_abs_err"] is not None for r in at_or_above)
                else None
            ),
        },
        "monotonicity": {
            "n_instances_checked": len(ok_metas),
            "n_monotone": sum(1 for m in ok_metas if m["is_monotone_from_rval"]),
        },
        "radius_brackets": {
            "n_records": len(radius_out),
            "n_both_finite_sandwich": len(
                [r for r in radius_out if r["sandwich_bracket_width"] is not None]
            ),
            "sandwich_bracket_widths": [
                r["sandwich_bracket_width"]
                for r in radius_out
                if r["sandwich_bracket_width"] is not None
            ],
            "n_both_finite_monotone": len(
                [r for r in radius_out if r["monotone_bracket_width"] is not None]
            ),
            "monotone_bracket_widths": [
                r["monotone_bracket_width"]
                for r in radius_out
                if r["monotone_bracket_width"] is not None
            ],
            "monotone_bracket_all_contain_r_true": all(
                r["monotone_bracket_contains_r_true"] is not False
                for r in radius_out
            ),
            "mu_hat_free_radius_errors": [
                r["mu_hat_free_err"]
                for r in radius_out
                if r["mu_hat_free_err"] is not None
            ],
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"[mean_bounds] wrote {out_dir}/shell_bounds.jsonl, "
          f"radius_bounds.jsonl and summary.json")


if __name__ == "__main__":
    main()
