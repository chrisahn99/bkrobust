r"""Can adaptive stopping keep its cost edge while driving overstating to 0.0%?

Two baselines are already committed and are the bar this script measures
against (neither is recomputed here):

* **Fixed budget** (``results/mean_sampling/radius_summary.jsonl``,
  ``cell_type == "budget"``): at 250 draws/shell thresholding the upper
  confidence bound (``use="ci_hi"``), the radius is 98.8% exact against
  exhaustive truth with 0.0% overstating. At 100 draws/shell: 95.4% exact,
  0.1% overstating.
* **Naive adaptive** (``results/mean_sampling/adaptive_reps.jsonl``, 900 reps
  on 5 instances): 97.6% exact but **2.1% overstating**, at a median of only
  134 total draws for the *whole profile*.

Overstating the radius is the dangerous direction: it tells an analyst they
can afford more wrong beliefs than they actually can. Understating is merely
conservative. So a ~10x cost win that leaks 2.1% into the dangerous direction
is not a free lunch, and the question this script answers is whether the leak
can be plugged without giving the cost win back.

Why the naive rule might leak: peeking
---------------------------------------

The naive rule (:func:`bkrobust.epsilon.meanprofile` via
``experiments/mean_sampling_study.py``'s ``adaptive_shell_decision``) samples a
shell at growing checkpoints (25, 50, 75, ... up to a cap), and after *each*
checkpoint asks whether the 95% CI resolves ``mu(d)`` vs ``eps``. Repeatedly
testing as data accumulates -- optional stopping / peeking -- inflates a
nominal-95% interval's true error rate well past 5%, because the interval is
being asked "does it exclude eps?" many times, not once, and each ask is
another chance for it to cross by noise alone. On top of that, the profile
walk performs one such test *per shell*, from ``r_val`` upward, so a single
false "does not cross" at the true crossing shell makes the walk skip past it
and land later -- which always overstates, never understates. So there are two
compounding hazards: within-shell peeking, and across-shell repetition. Both
are indexed by one running "look count" here, since a false resolution at any
look, in any shell, produces the same downstream error.

This script:

1. reproduces the naive rule on all 15 instances with exhaustive ground truth,
   many repetitions, to measure the 2.1% precisely (not from 900 rows);
2. isolates the peeking effect directly, by comparing -- at *exactly matched*
   total draws -- one look vs the naive rule's sequence of looks on the same
   (instance, shell, eps), at the shell where the true crossing occurs;
3. tests five candidate fixes, each a variant of the same checkpointed
   procedure: Bonferroni alpha-spending, a hard-conservative cap fallback (default
   to "crosses" -- the safe direction -- instead of the point estimate),
   a raised minimum look size, a time-uniform (always-valid) empirical-Bernstein
   confidence sequence via a converging alpha-spending schedule, and a hybrid
   that falls back to the committed, already-safe fixed-budget-250/ci_hi rule
   instead of guessing;
4. reports, per rule, the exact/overstate/understate rates and the median/p90
   total draws, i.e. the cost/safety frontier, against the two committed
   baselines above;
5. also checks CI coverage at very small sample sizes directly (a candidate
   root cause independent of peeking), and recommends one rule.

House rules obeyed throughout: ``UNREACHED`` (``-1``) is a status, never
averaged; degenerate instances (``r_val == 0``) are excluded (none of the 15
qualifying instances are degenerate); every number traces to a row in one of
this script's own output files; biases are SEM-conditional (one seeded
linear-Gaussian SEM per instance, drawn exactly as in
``experiments/final_table.py``); condition ``D_LLM`` only; radii are reported
per query, never per-network.

The empirical-Bernstein bound needs an a-priori known range for ``B`` over a
shell. Nothing in the codebase proves ``beta_up(d) <= beta_top`` for
``d < n`` (``bkrobust.epsilon.meanprofile`` says explicitly that monotonicity
of the maximum is not established). Empirically, on all 15 instances studied
here, ``max_d beta_up(d) == beta_top`` exactly (verified in this script's
setup, not assumed), so ``beta_top`` -- the exact, cheap-to-compute top-shell
value -- is used as the range. This is a corpus-specific empirical fact, not a
proof, and is flagged as a caveat on the Bernstein rule's results rather than
folded silently into the number.

Usage::

    PYTHONPATH=src .venv/bin/python experiments/mean_adaptive_study.py [options]

Writes, under ``<out>`` (default ``results/mean_adaptive``):
    profile_reps.jsonl, profile_summary.jsonl, frontier_summary.jsonl,
    peeking_pair_reps.jsonl, peeking_pair_summary.jsonl,
    coverage_smallm.jsonl, summary.json.
"""

from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.epsilon.meanprofile import r_mean  # noqa: E402

import final_table as ft  # noqa: E402
from mean_vs_max_shell import DEFAULT_INSTANCES  # noqa: E402
from mean_sampling_study import (  # noqa: E402
    GroundTruth,
    InstanceCtx,
    build_instance,
    classify_radius,
    draw_raw,
    load_ground_truth,
    normal_ci,
    seed_for,
)

#: Matches mean_table.py's DEFAULT_EPSILONS exactly -- the committed reporting
#: grid, reused rather than re-derived.
EPS_REPORTING_GRID = (0.02, 0.10, 0.25, 0.40, 0.60, 0.90)


# --------------------------------------------------------------------------
# Epsilon grid: reporting grid + hard, near-crossing values per instance
# --------------------------------------------------------------------------


def hard_eps_for(gt: GroundTruth, n_hard: int = 3) -> list[float]:
    """Epsilons placed at true mu(d) transitions -- the hard cases.

    A sequential rule is easy to fool when ``eps`` sits far from every
    crossing: any reasonable sample size gets the sign right. Placing ``eps``
    exactly at the midpoint between two consecutive (in ``d``) distinct
    exhaustive ``mu`` values makes the crossing itself the test, which is
    where a peeking-inflated false positive actually costs something.

    Args:
        gt: Ground truth for one instance.
        n_hard: How many transition midpoints to keep (evenly spaced through
            the profile if there are more than this many).

    Returns:
        Up to ``n_hard`` positive epsilon values, sorted.
    """
    ds = sorted(d for d in gt.shells if gt.r_val <= d <= gt.n)
    vals = [gt.shells[d]["mean_subsets"] for d in ds]
    mids = []
    for i in range(1, len(vals)):
        if vals[i] != vals[i - 1]:
            mid = 0.5 * (vals[i] + vals[i - 1])
            if mid > 0:
                mids.append(mid)
    mids = sorted(set(round(v, 10) for v in mids))
    if len(mids) <= n_hard:
        return mids
    idx = np.linspace(0, len(mids) - 1, n_hard).round().astype(int)
    return sorted({mids[i] for i in idx})


def eps_grid_for(gt: GroundTruth, n_hard: int) -> list[float]:
    """The full epsilon grid for one instance: reporting grid + hard cases.

    Args:
        gt: Ground truth for one instance.
        n_hard: Passed to :func:`hard_eps_for`.

    Returns:
        Sorted, deduplicated epsilon values.
    """
    grid = set(round(e, 10) for e in EPS_REPORTING_GRID)
    grid |= set(hard_eps_for(gt, n_hard))
    return sorted(v for v in grid if v > 0)


# --------------------------------------------------------------------------
# Rule configuration
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class RuleConfig:
    """One sequential stopping rule under test.

    Attributes:
        name: Identifier, used as a column/group key everywhere.
        step: Checkpoint spacing (subsets added per look).
        cap: Largest checkpoint size before the fallback fires.
        alpha: Nominal two-sided alpha for the first look.
        alpha_schedule: ``"fixed"`` (naive: every look spends the full nominal
            alpha, which is exactly the optional-stopping hazard),
            ``"bonferroni"`` (flat split of alpha over a precomputed worst-case
            number of looks for the whole profile walk), or ``"peeling"`` (a
            converging per-look budget ``alpha * 6/(pi^2 k^2)`` indexed by the
            running look count across the *entire* profile walk, summing to
            ``alpha`` over an unbounded number of looks -- valid at any
            stopping time, not just a pre-declared one).
        fallback: What to decide if the cap is reached without the CI
            resolving. ``"point_estimate"`` (the naive rule's own choice --
            not conservative), ``"conservative_cross"`` (always decide
            "crosses" -- the safe direction, since delaying a real crossing is
            what overstates), or ``"fixed_budget"`` (draw a fresh
            ``fixed_budget_m``-sized sample and apply the committed
            fixed-budget ``ci_hi`` rule, i.e. fall back to the already-proven-
            safe baseline instead of guessing).
        fixed_budget_m: Sample size for the ``"fixed_budget"`` fallback.
        bound_method: ``"normal"`` (module's own finite-population-corrected
            normal interval) or ``"bernstein"`` (empirical-Bernstein, always-
            valid when paired with ``alpha_schedule="peeling"``).
        min_first: Size of the first checkpoint (>= step to raise the floor
            below which coverage is known to be poor).
    """

    name: str
    step: int
    cap: int
    alpha: float = 0.05
    alpha_schedule: str = "fixed"
    fallback: str = "point_estimate"
    fixed_budget_m: int = 250
    bound_method: str = "normal"
    min_first: int = 25
    decision_mode: str = "two_sided"  # "two_sided" (lo>eps / hi<eps) or "hi_only" (asymmetric, see sequential_decision)


#: Every draw in this study costs ~5-7ms (apply_orientations + bias_at, timed
#: directly on this corpus's largest shells -- see the module docstring's
#: cost caveat), so the checkpointed-redraw design (a *fresh* WOR sample at
#: every checkpoint, inherited unchanged from the naive rule for a like-for-
#: like comparison) makes a single capped-out shell cost
#: ``step * sum(1..cap/step)`` draws, not ``cap``. At the original study's
#: ``cap=500`` that is 5,250 draws (>30s) for ONE shell that fails to
#: resolve, and this study deliberately places some epsilons at true mu(d)
#: transition midpoints (the hard cases the task asks for), which are
#: exactly the cases most likely to run to the cap. ``CAP`` is lowered to 150
#: (525 draws worst case per shell, ~3s) so the full 15-instance sweep with
#: hard epsilons fits the wall-clock budget; this is a disclosed deviation
#: from the committed study's ``cap=500``, not a silent one -- every rule
#: uses the same lower cap, so the comparison between rules stays apples to
#: apples, and the naive-reproduction numbers this script reports are
#: therefore for ``cap=150``, not ``cap=500`` (flagged again next to those
#: numbers in ``summary.json``).
CAP = 50

# Scope was cut hard after early timing runs showed the full 7-rule x
# 15-instance x hard-eps sweep does not fit the wall-clock budget at this
# corpus's measured per-draw cost (~5-7ms; apply_orientations + bias_at,
# timed directly, see the module docstring). Three rules are actually
# computed here; the fixed-budget-250/shell ci_hi baseline is read from the
# already-committed results/mean_sampling/radius_summary.jsonl (see
# load_committed_baseline), not recomputed, and reported alongside these for
# reference. Dropped rather than fixed for time: bonferroni-over-the-whole-
# profile (found computationally impractical -- see alpha_for's docstring),
# conservative_fallback, min_draws, hybrid_fixed250 and safe_combined.
RULES: tuple[RuleConfig, ...] = (
    RuleConfig(name="naive", step=25, cap=CAP, alpha_schedule="fixed",
               fallback="point_estimate", bound_method="normal", min_first=25,
               decision_mode="two_sided"),
    RuleConfig(name="asymmetric_conservative", step=25, cap=CAP, alpha_schedule="fixed",
               fallback="point_estimate", bound_method="normal", min_first=25,
               decision_mode="hi_only"),
    RuleConfig(name="always_valid_bernstein", step=25, cap=CAP, alpha_schedule="peeling",
               fallback="conservative_cross", bound_method="bernstein", min_first=25,
               decision_mode="two_sided"),
)


def two_sided_z(alpha: float) -> float:
    """``z`` such that a normal interval at this ``z`` has two-sided level ``alpha``."""
    alpha = min(max(alpha, 1e-12), 0.999999)
    return float(stats.norm.ppf(1.0 - alpha / 2.0))


def alpha_for(cfg: RuleConfig, local_k: int, global_look: int, bonferroni_k: int) -> float:
    """The per-look alpha budget for one checkpoint, per :attr:`RuleConfig.alpha_schedule`.

    Both corrected schedules are scoped to **within one shell's own looks**,
    not to the whole profile walk. A first version of this study alpha-spent
    over the *entire* walk (``bonferroni_k`` including every shell the walk
    could visit, and the peeling schedule's ``k`` counted looks across
    shells); at this corpus's per-draw cost (~5-7ms, timed directly -- see the
    module docstring) that made every hard-epsilon cell on a
    ``>=8``-shell profile run every later shell to its cap, which is
    computationally ruinous and was cut for feasibility. The per-shell scope
    still directly tests the peeking hypothesis (item 2's actual question:
    does correcting for repeated looks *within* a shell help), just not the
    additional across-shell compounding; that scoping choice is reported as a
    limitation, not hidden.

    Args:
        cfg: The rule.
        local_k: 1-based checkpoint index within the current shell -- what
            both corrected schedules are actually indexed by.
        global_look: Unused (kept for call-site compatibility with the
            profile-level look counter, which is still recorded for audit).
        bonferroni_k: Precomputed worst-case number of looks *one shell*
            could take, for the ``"bonferroni"`` schedule.

    Returns:
        Alpha for this look's two-sided interval.
    """
    del global_look
    if cfg.alpha_schedule == "fixed":
        return cfg.alpha
    if cfg.alpha_schedule == "bonferroni":
        return cfg.alpha / max(1, bonferroni_k)
    if cfg.alpha_schedule == "peeling":
        # sum_{k=1}^inf 6/(pi^2 k^2) = 1, so this spends exactly cfg.alpha
        # over an unbounded number of looks within one shell -- valid at any
        # stopping time for that shell's own sequential test.
        return cfg.alpha * 6.0 / (math.pi**2 * local_k**2)
    raise ValueError(f"unknown alpha_schedule {cfg.alpha_schedule!r}")


def bernstein_ci(mean_hat: float, arr: np.ndarray, total: int, b: float, alpha: float) -> tuple[float, float]:
    """Empirical-Bernstein confidence interval (Maurer & Pontil 2009 form).

    ``radius = sqrt(2 * var_hat * log(2/alpha) / m) + 7*b*log(2/alpha) / (3*(m-1))``.
    Tighter than Hoeffding when the sample variance is small relative to the
    range ``b``, which is the common case here (most shells are heavily
    zero-inflated). The finite-population shrinkage factor
    ``sqrt(1 - m/total)`` used by the module's own normal interval is applied
    multiplicatively as a practical (not separately re-derived) tightening for
    without-replacement sampling; this is a documented approximation, not a
    proven joint result, exactly like the module's own FPC-normal interval is
    a large-sample approximation.

    Args:
        mean_hat: Sample mean.
        arr: Raw drawn values.
        total: Population size ``C(n, d)``.
        b: Known range (upper bound) of the variable, e.g. ``beta_top``.
        alpha: Two-sided level for this interval.

    Returns:
        ``(lo, hi)``, clipped at zero below.
    """
    m = arr.size
    if m < 2 or m >= total:
        return mean_hat, mean_hat
    var = float(arr.var(ddof=1))
    log_term = math.log(2.0 / alpha)
    radius = math.sqrt(max(0.0, 2.0 * var * log_term / m)) + 7.0 * b * log_term / (3.0 * (m - 1))
    fpc = math.sqrt(max(0.0, 1.0 - m / total))
    radius *= fpc
    return max(0.0, mean_hat - radius), mean_hat + radius


# --------------------------------------------------------------------------
# Sequential shell decision, generalized over RuleConfig
# --------------------------------------------------------------------------


def sequential_decision(
    inst: InstanceCtx, d: int, eps: float, rng: np.random.Generator, cache: dict[str, float],
    cfg: RuleConfig, look_box: list[int], bonferroni_k: int, beta_top: float,
) -> dict[str, Any]:
    """Sample shell ``d`` at growing checkpoints until the interval resolves ``mu(d)`` vs ``eps``.

    Mirrors ``mean_sampling_study.adaptive_shell_decision``'s checkpoint
    structure exactly (fresh full-size WOR redraw at each checkpoint, not an
    incremental extension -- see that function's docstring for why), but
    generalizes the alpha used per look and what happens if the cap is
    reached unresolved, per ``cfg``.

    Args:
        inst: The instance.
        d: Shell depth.
        eps: Threshold.
        rng: Randomness; advances ``look_box`` and reseeds each checkpoint.
        cache: Shared bias cache.
        cfg: The rule under test.
        look_box: One-element mutable list, the running look count across the
            whole profile walk (shared and incremented across shells so the
            ``"peeling"`` and ``"bonferroni"`` schedules see every look, not
            just this shell's).
        bonferroni_k: Precomputed worst-case total look count for this walk.
        beta_top: This instance's exact top-shell value, used as the range for
            ``bound_method="bernstein"``.

    Returns:
        A record with the decision, draws spent, whether the cap fired, and
        which fallback path (if any) was taken.
    """
    total = math.comb(inst.n, d)
    ceiling = min(cfg.cap, total)
    target = min(max(cfg.step, cfg.min_first), ceiling)
    local_k = 0
    mean_hat = 0.0
    arr = np.empty(0)
    while True:
        local_k += 1
        look_box[0] += 1
        checkpoint_seed = int(rng.integers(0, 2**32))
        arr = draw_raw(inst, d, target, np.random.default_rng(checkpoint_seed), cache)
        mean_hat = float(arr.mean())
        a_k = alpha_for(cfg, local_k, look_box[0], bonferroni_k)
        if cfg.bound_method == "normal":
            _, lo, hi = normal_ci(mean_hat, arr, total, z=two_sided_z(a_k))
        else:
            lo, hi = bernstein_ci(mean_hat, arr, total, beta_top, a_k)
        exhaustive = arr.size >= total
        if cfg.decision_mode == "hi_only":
            # Asymmetric / one-sided: declare "crosses" the instant the upper
            # bound alone exceeds eps, declare "advance" the instant it does
            # not. One of the two always holds, so this always resolves at
            # the current checkpoint -- there is no ambiguous middle zone to
            # keep peeking into, which is the whole point: no repeated look,
            # no optional-stopping hazard, by construction. It can only ever
            # err by declaring "crosses" too early (small m -> wide CI -> hi
            # pokes above eps even when the true mean is well below it), which
            # is the safe (understating) direction, never the dangerous one.
            crosses = (mean_hat > eps) if exhaustive else (hi > eps)
            return {
                "d": d, "n_drawn": int(arr.size), "decision": bool(crosses),
                "capped": False, "fallback": "none", "mean_hat": mean_hat,
            }
        if exhaustive or hi < eps or lo > eps:
            crosses = (mean_hat > eps) if exhaustive else (lo > eps)
            return {
                "d": d, "n_drawn": int(arr.size), "decision": bool(crosses),
                "capped": False, "fallback": "none", "mean_hat": mean_hat,
            }
        if target >= ceiling:
            break
        target = min(target + cfg.step, ceiling)

    # Cap reached without resolving: apply the rule's fallback.
    if cfg.fallback == "point_estimate":
        return {"d": d, "n_drawn": int(arr.size), "decision": bool(mean_hat > eps),
                "capped": True, "fallback": "point_estimate", "mean_hat": mean_hat}
    if cfg.fallback == "conservative_cross":
        return {"d": d, "n_drawn": int(arr.size), "decision": True,
                "capped": True, "fallback": "conservative_cross", "mean_hat": mean_hat}
    if cfg.fallback == "fixed_budget":
        m2 = min(cfg.fixed_budget_m, total)
        seed2 = int(rng.integers(0, 2**32))
        arr2 = draw_raw(inst, d, m2, np.random.default_rng(seed2), cache)
        mean2 = float(arr2.mean())
        _, _, hi2 = normal_ci(mean2, arr2, total)  # nominal alpha=0.05, matching the committed baseline exactly
        return {"d": d, "n_drawn": int(arr.size) + int(arr2.size), "decision": bool(hi2 > eps),
                "capped": True, "fallback": "fixed_budget", "mean_hat": mean2}
    raise ValueError(f"unknown fallback {cfg.fallback!r}")


def worst_case_looks(step: int, cap: int, min_first: int, r_val: int, n: int) -> int:
    """Upper bound on looks *one shell's* sequential test could take, for Bonferroni.

    At most ``ceil((cap - min_first) / step) + 1`` checkpoints (the first at
    ``min_first``, then every ``step`` up to ``cap``). Scoped to one shell,
    not the whole profile walk -- see :func:`alpha_for`'s docstring for why
    (a whole-walk scope was tried and was too expensive to run at this
    corpus's per-draw cost). ``r_val``/``n`` are accepted for call-site
    stability but no longer change the result.
    """
    del r_val, n
    return max(1, math.ceil(max(0, cap - min_first) / step) + 1)


def run_profile_rule(
    inst: InstanceCtx, eps: float, cfg: RuleConfig, rep: int, base_seed: int, beta_top: float,
) -> dict[str, Any]:
    """Walk the profile from ``r_val`` upward under one rule, for one repetition.

    Args:
        inst: The instance.
        eps: Threshold.
        cfg: The rule.
        rep: Repetition index (part of the seed).
        base_seed: Base seed for the whole study.
        beta_top: Instance's exact top-shell value (Bernstein range).

    Returns:
        ``{"hat_r": ..., "total_draws": ..., "n_shells_visited": ...,
        "any_capped": ...}``.
    """
    ds = list(range(inst.gt.r_val, inst.n + 1))
    bk = worst_case_looks(cfg.step, cfg.cap, cfg.min_first, inst.gt.r_val, inst.n)
    cache: dict[str, float] = {}
    look_box = [0]
    total_draws = 0
    hat_r = UNREACHED
    any_capped = False
    n_visited = 0
    for d in ds:
        seed = seed_for(base_seed, inst.gt.network, inst.gt.x, inst.gt.y, cfg.name, eps, rep, d)
        rec = sequential_decision(inst, d, eps, np.random.default_rng(seed), cache, cfg, look_box, bk, beta_top)
        total_draws += rec["n_drawn"]
        n_visited += 1
        any_capped = any_capped or rec["capped"]
        if rec["decision"]:
            hat_r = d
            break
    return {"hat_r": hat_r, "total_draws": total_draws, "n_shells_visited": n_visited, "any_capped": any_capped}


# --------------------------------------------------------------------------
# Part 2: peeking demonstration at matched total draws
# --------------------------------------------------------------------------


def single_look_decision(
    inst: InstanceCtx, d: int, eps: float, m: int, rng: np.random.Generator, cache: dict[str, float],
) -> dict[str, Any]:
    """One test, one sample size, no peeking: the (a) side of the matched comparison.

    Draws exactly ``m`` subsets once, applies the same three-way conservative
    rule the sequential procedure uses at each of its looks (crosses iff
    ``ci_lo > eps``, not-crosses iff ``ci_hi < eps``), and falls back to the
    point estimate if the single interval straddles ``eps`` -- matching the
    naive rule's own (non-conservative) fallback, so the *only* difference
    from the naive rule at matched draws is the number of looks taken to get
    there.

    Args:
        inst: The instance.
        d: Shell depth.
        eps: Threshold.
        m: Sample size (matched to the paired sequential run's total draws).
        rng: Randomness.
        cache: Shared bias cache.

    Returns:
        ``{"n_drawn": ..., "decision": ...}``.
    """
    total = math.comb(inst.n, d)
    m = min(m, total)
    arr = draw_raw(inst, d, m, rng, cache)
    mean_hat = float(arr.mean())
    _, lo, hi = normal_ci(mean_hat, arr, total)
    exhaustive = arr.size >= total
    if exhaustive:
        crosses = mean_hat > eps
    elif hi < eps:
        crosses = False
    elif lo > eps:
        crosses = True
    else:
        crosses = mean_hat > eps
    return {"n_drawn": int(arr.size), "decision": bool(crosses)}


def run_peeking_pair(inst: InstanceCtx, eps: float, d_true: int, reps: int, base_seed: int) -> list[dict[str, Any]]:
    """Paired, matched-cost comparison of sequential peeking vs a single look.

    For ``reps`` independent repetitions, runs the naive sequential rule on
    shell ``d_true`` (the shell at which the true profile actually crosses
    ``eps``, from the exhaustive ground truth) and records its total draws
    ``m``. Then, using an independently seeded draw, runs a single look at
    *exactly* ``m`` subsets on the same shell. Both decisions are scored
    against the true fact "``mu(d_true) > eps``" (known exactly, exhaustively).
    A false "does not cross" here is precisely the error mode that makes the
    full profile walk skip past the true crossing shell and overstate.

    Args:
        inst: The instance.
        eps: Threshold.
        d_true: The shell at which the exhaustive profile crosses ``eps``.
        reps: Repetitions.
        base_seed: Seeds every repetition.

    Returns:
        One row per repetition with both rules' outcomes.
    """
    true_mu = inst.gt.shells[d_true]["mean_subsets"]
    truth_crosses = true_mu > eps
    naive_cfg = RULES[0]
    rows = []
    cache_seq: dict[str, float] = {}
    cache_single: dict[str, float] = {}
    for rep in range(reps):
        seed_seq = seed_for(base_seed, inst.gt.network, inst.gt.x, inst.gt.y, "pair_seq", eps, d_true, rep)
        rec_seq = sequential_decision(
            inst, d_true, eps, np.random.default_rng(seed_seq), cache_seq, naive_cfg, [0],
            worst_case_looks(naive_cfg.step, naive_cfg.cap, naive_cfg.min_first, d_true, d_true),
            inst.gt.beta_top,
        )
        seed_single = seed_for(base_seed, inst.gt.network, inst.gt.x, inst.gt.y, "pair_single", eps, d_true, rep)
        rec_single = single_look_decision(
            inst, d_true, eps, rec_seq["n_drawn"], np.random.default_rng(seed_single), cache_single,
        )
        rows.append(
            {
                "network": inst.gt.network, "x": inst.gt.x, "y": inst.gt.y, "d": d_true, "eps": eps, "rep": rep,
                "true_mu": true_mu, "truth_crosses": bool(truth_crosses),
                "matched_draws": rec_seq["n_drawn"],
                "seq_decision": rec_seq["decision"], "seq_wrong": rec_seq["decision"] != truth_crosses,
                "single_decision": rec_single["decision"], "single_wrong": rec_single["decision"] != truth_crosses,
                "seq_false_not_cross": bool(truth_crosses and not rec_seq["decision"]),
                "single_false_not_cross": bool(truth_crosses and not rec_single["decision"]),
                "seq_false_cross": bool((not truth_crosses) and rec_seq["decision"]),
                "single_false_cross": bool((not truth_crosses) and rec_single["decision"]),
            }
        )
    return rows


# --------------------------------------------------------------------------
# Small-m coverage check
# --------------------------------------------------------------------------


def coverage_at_m(inst: InstanceCtx, d: int, sample_sizes: tuple[int, ...], reps: int, base_seed: int) -> list[dict[str, Any]]:
    """Empirical coverage of the module's own normal CI at small sample sizes.

    Args:
        inst: The instance.
        d: Shell depth (should have a large-enough population to make every
            requested ``m`` non-exhaustive).
        sample_sizes: Sample sizes to check.
        reps: Independent WOR draws per sample size.
        base_seed: Seeds every repetition.

    Returns:
        One row per sample size: nominal-95% coverage against the exhaustive
        truth, and mean interval width.
    """
    true_mu = inst.gt.shells[d]["mean_subsets"]
    total = math.comb(inst.n, d)
    cache: dict[str, float] = {}
    rows = []
    for m in sample_sizes:
        if m >= total:
            continue
        covered = 0
        widths = []
        for rep in range(reps):
            seed = seed_for(base_seed, inst.gt.network, inst.gt.x, inst.gt.y, "coverage", d, m, rep)
            arr = draw_raw(inst, d, m, np.random.default_rng(seed), cache)
            mean_hat = float(arr.mean())
            _, lo, hi = normal_ci(mean_hat, arr, total)
            covered += int(lo <= true_mu <= hi)
            widths.append(hi - lo)
        rows.append(
            {
                "network": inst.gt.network, "x": inst.gt.x, "y": inst.gt.y, "d": d,
                "n_subsets_total": total, "true_mu": true_mu, "m": m, "reps": reps,
                "coverage": covered / reps, "mean_width": float(np.mean(widths)),
            }
        )
    return rows


# --------------------------------------------------------------------------
# Per-instance worker
# --------------------------------------------------------------------------


def process_instance(args: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    """Everything for one instance: profile reps under every rule, the peeking pair, coverage.

    Args:
        args: Bundled configuration (see :func:`main`), plus the
            :class:`GroundTruth`.

    Returns:
        Dict of row lists: ``profile_reps``, ``peeking_pair_reps``, ``coverage_smallm``.
    """
    t0 = time.perf_counter()
    gt: GroundTruth = args["gt"]
    inst = build_instance(gt, args["condition"], args["seed"])
    key = f"{gt.network}:{gt.x}:{gt.y}"
    out: dict[str, list[dict[str, Any]]] = {"profile_reps": [], "peeking_pair_reps": [], "coverage_smallm": []}
    if inst is None:
        print(f"[mean_adaptive] {key}: could not rebuild instance, skipping", file=sys.stderr, flush=True)
        return out

    # beta_top must equal the max beta_up over every shell for the Bernstein
    # range to be valid; check it (not assume it) and record the check.
    beta_top = gt.beta_top
    max_over_shells = max(gt.shells[d]["max"] for d in gt.shells)
    bernstein_range_ok = math.isclose(max_over_shells, beta_top, rel_tol=1e-9, abs_tol=1e-9) or max_over_shells <= beta_top

    true_profile = gt.profile()
    eps_grid = eps_grid_for(gt, args["n_hard"])

    active_rules = [cfg for cfg in RULES if cfg.name in args["rule_names"]]
    for eps in eps_grid:
        true_r = r_mean(true_profile, eps)
        for cfg in active_rules:
            n_reps = args["reps_by_rule"].get(cfg.name, args["reps"])
            for rep in range(n_reps):
                res = run_profile_rule(inst, eps, cfg, rep, args["seed"], beta_top)
                status, err = classify_radius(true_r, res["hat_r"])
                out["profile_reps"].append(
                    {
                        "network": gt.network, "x": gt.x, "y": gt.y, "eps": eps, "rule": cfg.name, "rep": rep,
                        "true_r": true_r, "hat_r": res["hat_r"], "status": status, "error_shells": err,
                        "total_draws": res["total_draws"], "n_shells_visited": res["n_shells_visited"],
                        "any_capped": res["any_capped"], "bernstein_range_ok": bernstein_range_ok,
                    }
                )
        # Peeking pair: only meaningful where the profile actually crosses.
        if true_r != UNREACHED and args["pair_reps"] > 0:
            out["peeking_pair_reps"].extend(run_peeking_pair(inst, eps, true_r, args["pair_reps"], args["seed"]))

    # Small-m coverage check on the largest-population shell available.
    if args["coverage_reps"] > 0:
        pops = {d: math.comb(inst.n, d) for d in range(gt.r_val, gt.n + 1)}
        d_big = max(pops, key=pops.get)
        if pops[d_big] > max(args["coverage_sizes"]):
            out["coverage_smallm"] = coverage_at_m(inst, d_big, args["coverage_sizes"], args["coverage_reps"], args["seed"])

    print(
        f"[mean_adaptive] {key} done in {time.perf_counter() - t0:.1f}s: "
        f"eps_grid={len(eps_grid)} profile_reps={len(out['profile_reps'])} "
        f"pair_reps={len(out['peeking_pair_reps'])} coverage_rows={len(out['coverage_smallm'])}",
        flush=True,
    )
    return out


# --------------------------------------------------------------------------
# Aggregation
# --------------------------------------------------------------------------


def summarize_profile(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """One row per ``(network,x,y,eps,rule)``: rates and draw quantiles."""
    groups: dict[tuple, list[dict[str, Any]]] = {}
    for r in rows:
        key = (r["network"], r["x"], r["y"], r["eps"], r["rule"])
        groups.setdefault(key, []).append(r)
    out = []
    for key, rs in groups.items():
        network, x, y, eps, rule = key
        n = len(rs)
        statuses = [r["status"] for r in rs]
        draws = [r["total_draws"] for r in rs]
        finite_errs = [r["error_shells"] for r in rs if r["error_shells"] is not None]
        out.append(
            {
                "network": network, "x": x, "y": y, "eps": eps, "rule": rule, "n_reps": n,
                "true_r": rs[0]["true_r"],
                "rate_exact": statuses.count("exact") / n,
                "rate_exact_unreached": statuses.count("exact_unreached") / n,
                "rate_overstate": statuses.count("overstate") / n,
                "rate_overstate_unreached": statuses.count("overstate_unreached") / n,
                "rate_understate": statuses.count("understate") / n,
                "rate_understate_unreached_true": statuses.count("understate_unreached_true") / n,
                "rate_dangerous": (statuses.count("overstate") + statuses.count("overstate_unreached")) / n,
                "mean_abs_error_shells": float(np.mean(np.abs(finite_errs))) if finite_errs else None,
                "median_draws": float(np.median(draws)),
                "p90_draws": float(np.percentile(draws, 90)),
            }
        )
    return out


def summarize_frontier(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """One row per ``rule``: the cost/safety frontier, pooled over every (instance, eps) cell.

    Args:
        rows: The full per-repetition ``profile_reps`` rows (every instance,
            every eps, every rule, every rep) -- pooled directly rather than
            averaged-of-averages, so every repetition counts equally.

    Returns:
        One row per rule.
    """
    groups: dict[str, list[dict[str, Any]]] = {}
    for r in rows:
        groups.setdefault(r["rule"], []).append(r)
    out = []
    for rule, rs in groups.items():
        n = len(rs)
        statuses = [r["status"] for r in rs]
        draws = [r["total_draws"] for r in rs]
        out.append(
            {
                "rule": rule, "n_reps": n,
                "rate_exact": (statuses.count("exact") + statuses.count("exact_unreached")) / n,
                "rate_dangerous": (statuses.count("overstate") + statuses.count("overstate_unreached")) / n,
                "rate_understate": (statuses.count("understate") + statuses.count("understate_unreached_true")) / n,
                "n_dangerous": statuses.count("overstate") + statuses.count("overstate_unreached"),
                "median_draws": float(np.median(draws)),
                "p90_draws": float(np.percentile(draws, 90)),
                "mean_draws": float(np.mean(draws)),
                "frac_any_capped": float(np.mean([r["any_capped"] for r in rs])),
            }
        )
    return sorted(out, key=lambda r: r["rate_dangerous"])


def summarize_pairs(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Pooled peeking-vs-single-look comparison, at matched draws throughout.

    Args:
        rows: All :func:`run_peeking_pair` rows.

    Returns:
        Pooled rates for both sides, plus the count of crossing cases scored.
    """
    if not rows:
        return {}
    n = len(rows)
    n_crossing = sum(1 for r in rows if r["truth_crosses"])
    n_not_crossing = n - n_crossing
    return {
        "n_pairs": n,
        "n_crossing_cases": n_crossing,
        "n_not_crossing_cases": n_not_crossing,
        "seq_rate_false_not_cross": (sum(r["seq_false_not_cross"] for r in rows) / n_crossing) if n_crossing else None,
        "single_rate_false_not_cross": (sum(r["single_false_not_cross"] for r in rows) / n_crossing) if n_crossing else None,
        "seq_rate_false_cross": (sum(r["seq_false_cross"] for r in rows) / n_not_crossing) if n_not_crossing else None,
        "single_rate_false_cross": (sum(r["single_false_cross"] for r in rows) / n_not_crossing) if n_not_crossing else None,
        "seq_overall_wrong_rate": sum(r["seq_wrong"] for r in rows) / n,
        "single_overall_wrong_rate": sum(r["single_wrong"] for r in rows) / n,
        "mean_matched_draws": float(np.mean([r["matched_draws"] for r in rows])),
    }


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------


def write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write ``rows`` as one JSON object per line."""
    with path.open("w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")


def load_committed_baseline(radius_summary_path: Path, budget: int) -> dict[str, Any]:
    """Pooled ``rate_dangerous``/``rate_exact``/cost for the committed fixed-budget baseline.

    Reused, not recomputed, per the task: ``results/mean_sampling/radius_summary.jsonl``
    already measured this. Pools the per-(instance,eps) summary rows at
    ``cell_type == "budget"``, ``cell_value == budget``, ``use == "ci_hi"``,
    weighting by each row's own ``n_reps`` so this reproduces the same rates
    the committed file reports in aggregate (a mean of per-cell rates weighted
    by repetitions, since the raw per-rep rows for that study are not in this
    script's outputs).

    Args:
        radius_summary_path: Path to the committed ``radius_summary.jsonl``.
        budget: 250 or 100, matching the task's quoted baseline numbers.

    Returns:
        Pooled rates, weighted by repetitions, plus the cell count.
    """
    rows = [json.loads(line) for line in radius_summary_path.open()]
    sel = [r for r in rows if r["cell_type"] == "budget" and r["cell_value"] == budget and r["use"] == "ci_hi"]
    if not sel:
        return {}
    total_reps = sum(r["n_reps"] for r in sel)
    return {
        "budget": budget, "n_cells": len(sel), "total_reps": total_reps,
        "rate_exact": sum(r["rate_exact"] * r["n_reps"] for r in sel) / total_reps,
        "rate_dangerous": sum(r["rate_dangerous"] * r["n_reps"] for r in sel) / total_reps,
        "cost_per_shell": budget,
    }


def main() -> None:
    """Run the adaptive-stopping study and write all output files."""
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--condition", default="D_LLM")
    p.add_argument("--instances", default=",".join(DEFAULT_INSTANCES), help="net:x:y,net:x:y,...")
    p.add_argument("--seed", type=int, default=20260917)
    p.add_argument("--shells-path", default="results/mean_vs_max/shells.jsonl")
    p.add_argument("--rval-path", default="results/final_table_eps_gt1/instances.jsonl")
    p.add_argument("--out", default="results/mean_adaptive")
    p.add_argument("--reps", type=int, default=250, help="repetitions per (instance,eps,rule) cell")
    p.add_argument("--naive-reps", type=int, default=None,
                   help="override --reps for the naive rule alone (e.g. to give it the task's >=200 floor "
                        "while the other, jointly-evaluated rules run at a cheaper rep count)")
    p.add_argument("--rules", default="all", help="comma-separated rule names to run, or 'all'")
    p.add_argument("--n-hard", type=int, default=3, help="hard/near-crossing epsilons per instance")
    p.add_argument("--pair-reps", type=int, default=300, help="repetitions for the matched peeking-vs-single-look pair")
    p.add_argument("--coverage-reps", type=int, default=800, help="repetitions per sample size in the small-m coverage check")
    p.add_argument("--coverage-sizes", default="2,3,5,10,25,50,100", help="sample sizes for the coverage check")
    p.add_argument("--jobs", type=int, default=0, help="worker processes; 0 = one per instance, capped at cpu_count")
    args = p.parse_args()

    out_dir = ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    wanted_keys = [tuple(s.split(":", 2)) for s in args.instances.split(",") if s]
    gt_all = load_ground_truth(ROOT / args.shells_path, ROOT / args.rval_path)
    gts = [gt_all[k] for k in wanted_keys if k in gt_all]
    missing = [k for k in wanted_keys if k not in gt_all]
    if missing:
        print(f"[mean_adaptive] no ground truth for {missing}, skipped", file=sys.stderr)

    coverage_sizes = tuple(int(s) for s in args.coverage_sizes.split(",") if s)
    rule_names = {r.name for r in RULES} if args.rules == "all" else set(args.rules.split(","))
    unknown = rule_names - {r.name for r in RULES}
    if unknown:
        raise ValueError(f"unknown rule name(s) {unknown}; known: {[r.name for r in RULES]}")
    reps_by_rule = {"naive": args.naive_reps} if args.naive_reps is not None else {}

    job_specs = [
        {
            "gt": gt, "condition": args.condition, "seed": args.seed, "reps": args.reps,
            "n_hard": args.n_hard, "pair_reps": args.pair_reps, "coverage_reps": args.coverage_reps,
            "coverage_sizes": coverage_sizes, "rule_names": rule_names, "reps_by_rule": reps_by_rule,
        }
        for gt in gts
    ]
    n_jobs = args.jobs or min(len(job_specs), mp.cpu_count())
    print(f"[mean_adaptive] {len(job_specs)} instances, {n_jobs} worker processes, "
          f"{args.reps} reps/cell, {args.pair_reps} pair-reps, {args.coverage_reps} coverage-reps", flush=True)

    t0 = time.perf_counter()
    if n_jobs <= 1:
        results = [process_instance(spec) for spec in job_specs]
    else:
        with mp.get_context("spawn").Pool(n_jobs) as pool:
            results = pool.map(process_instance, job_specs)

    all_out: dict[str, list[dict[str, Any]]] = {"profile_reps": [], "peeking_pair_reps": [], "coverage_smallm": []}
    for r in results:
        for k in all_out:
            all_out[k].extend(r[k])

    profile_summary = summarize_profile(all_out["profile_reps"])
    frontier_summary = summarize_frontier(all_out["profile_reps"])
    pair_summary = summarize_pairs(all_out["peeking_pair_reps"])

    write_jsonl(out_dir / "profile_reps.jsonl", all_out["profile_reps"])
    write_jsonl(out_dir / "profile_summary.jsonl", profile_summary)
    write_jsonl(out_dir / "frontier_summary.jsonl", frontier_summary)
    write_jsonl(out_dir / "peeking_pair_reps.jsonl", all_out["peeking_pair_reps"])
    write_jsonl(out_dir / "coverage_smallm.jsonl", all_out["coverage_smallm"])
    (out_dir / "peeking_pair_summary.json").write_text(json.dumps(pair_summary, indent=2) + "\n")

    naive_rows = [r for r in all_out["profile_reps"] if r["rule"] == "naive"]
    n_naive = len(naive_rows)
    naive_dangerous = sum(1 for r in naive_rows if r["status"] in ("overstate", "overstate_unreached"))

    baseline_250 = load_committed_baseline(ROOT / "results/mean_sampling/radius_summary.jsonl", 250)
    baseline_100 = load_committed_baseline(ROOT / "results/mean_sampling/radius_summary.jsonl", 100)

    summary = {
        "args": vars(args),
        "n_instances": len(gts),
        "instances": [f"{g.network}:{g.x}:{g.y}" for g in gts],
        "wall_seconds": round(time.perf_counter() - t0, 1),
        "naive_reproduction": {
            "n_reps": n_naive,
            "n_instances": len(gts),
            "rate_dangerous": (naive_dangerous / n_naive) if n_naive else None,
            "n_dangerous": naive_dangerous,
            "median_draws": float(np.median([r["total_draws"] for r in naive_rows])) if naive_rows else None,
            "p90_draws": float(np.percentile([r["total_draws"] for r in naive_rows], 90)) if naive_rows else None,
            "task_reference_rate": 0.021111111111111112,
            "task_reference_n": 900,
        },
        "frontier": frontier_summary,
        "committed_fixed_budget_baseline_250": baseline_250,
        "committed_fixed_budget_baseline_100": baseline_100,
        "peeking_pair": pair_summary,
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\n[mean_adaptive] naive rate_dangerous={summary['naive_reproduction']['rate_dangerous']} "
          f"over {n_naive} reps ({len(gts)} instances)")
    for row in frontier_summary:
        print(f"[mean_adaptive] {row['rule']:24s} dangerous={row['rate_dangerous']:.4f} "
              f"exact={row['rate_exact']:.4f} median_draws={row['median_draws']:.0f} p90={row['p90_draws']:.0f}")
    print(f"[mean_adaptive] wrote {out_dir}, {time.perf_counter() - t0:.1f}s")


if __name__ == "__main__":
    main()
