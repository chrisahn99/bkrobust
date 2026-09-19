r"""Sampled estimation of ``mu(d)``: accuracy, CI coverage, and radius fidelity.

``mu(d)`` (:mod:`bkrobust.epsilon.meanprofile`) is the subset-weighted mean bias
over retraction shell ``d``. Shell ``d`` has ``C(|K_G0|, d)`` subsets, which
explodes, so in practice ``mu(d)`` must be estimated by sampling. This script
answers the question a practitioner actually asks: *if I can only afford x% of
a shell (or N draws), what do I get, and how sure am I?*

It reuses two pieces of already-computed ground truth rather than
re-enumerating anything:

* ``results/mean_vs_max/shells.jsonl`` -- exhaustive per-shell statistics
  (``mean_subsets`` is exactly ``mu(d)``) for the 15 instances studied there.
* ``results/final_table_eps_gt1/instances.jsonl`` -- ``r_val`` per instance.

Four things are measured, matched to the task:

1. **Estimator accuracy.** :func:`shell_stat` (the module under test) is called
   many independent times per (instance, shell, sampling fraction or budget)
   cell, and the resulting ``mu(d)`` estimates are compared to the exhaustive
   truth: bias and RMSE.
2. **CI coverage.** Two 95% interval constructions are compared against the
   truth over the same repetitions: the module's own normal interval with a
   finite-population correction, and a bootstrap percentile interval built
   from the same raw draws. Coverage this script measures is the fraction of
   nominal-95% intervals that actually contain the exhaustive ``mu(d)``.
3. **Radius fidelity.** The deliverable that actually gets reported is not
   ``mu(d)`` but ``r_mean(eps)``. Full sampled profiles are built (all shells
   from ``r_val`` to ``|K_G0|``) and compared, per instance and per epsilon,
   against the radius read off the exhaustive profile -- exact-match rate,
   signed error in shells, and, critically, which direction errors run:
   overstating robustness (dangerous) or understating it (safe).
4. **A conservative reporting rule**, plus an exploratory adaptive/sequential
   sampler (stop early once the CI resolves whether ``mu(d) > eps``).

House rules this script obeys throughout:

* ``UNREACHED`` (``-1``) is a status, never averaged or treated as a large
  number; radius comparisons classify it explicitly (see
  :func:`classify_radius`).
* Degenerate instances (``r_val == 0``) are excluded (none of the 15 studied
  instances are degenerate, verified at load time).
* Every number in the summary files traces back to a repetition recorded in
  the paired ``*_reps.jsonl`` file, which in turn traces back to a
  ``shell_stat`` call or a from-ground-truth exact profile.
* Biases are SEM-conditional: one seeded linear-Gaussian SEM per instance,
  drawn exactly as in ``experiments/final_table.py`` / ``mean_vs_max_shell.py``.
* Condition ``D_LLM`` only.

Two independent verification checks run before the main study and are written
to ``<out>/verification.json``: (a) that :func:`shell_stat`'s sampling and a
hand-rolled replica using the same RNG seed agree exactly on the drawn subsets
and the resulting mean, and (b) that the finite-population-corrected standard
error the module reports matches the true sampling variance of the mean under
without-replacement sampling, computed from one fully-enumerated shell.

Usage::

    PYTHONPATH=src .venv/bin/python experiments/mean_sampling_study.py [options]

Writes, under ``<out>`` (default ``results/mean_sampling``):
    verification.json, shell_reps.jsonl, shell_summary.jsonl,
    radius_reps.jsonl, radius_summary.jsonl, adaptive_reps.jsonl,
    adaptive_summary.jsonl, summary.json.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import multiprocessing as mp
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag, random_sem  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.epsilon.bias import bias_at, knowledge_of, make_context  # noqa: E402
from bkrobust.epsilon.meanprofile import (  # noqa: E402
    ShellStat,
    _sample_subsets,  # the exact WOR sampling scheme under test; reused, not reimplemented
    r_mean,
    shell_stat,
)

import final_table as ft  # noqa: E402
from mean_vs_max_shell import DEFAULT_INSTANCES  # noqa: E402

Z95 = 1.959963984540054  # scipy.stats.norm.ppf(0.975); matches meanprofile's default z=1.96 closely

#: Sampling fractions of a shell's population, and fixed sample-size budgets.
DEFAULT_FRACTIONS = (0.01, 0.02, 0.05, 0.10, 0.25, 0.50)
DEFAULT_BUDGETS = (25, 50, 100, 250, 500, 1000)

#: Fractions of each instance's own beta_top (= mu at full retraction) used as
#: the per-instance epsilon grid for the radius study. Data-driven rather than
#: a fixed absolute grid, since mu(d)'s magnitude varies by >50x across the 15
#: instances (results/mean_vs_max/shells.jsonl).
EPS_QUANTILES = (0.10, 0.25, 0.50, 0.75, 0.90)


def seed_for(*parts: object) -> int:
    """Deterministic seed from arbitrary parts, stable across processes.

    Args:
        *parts: Anything ``str()``-able; joined with ``:``.

    Returns:
        A seed in ``[0, 2**32)``.
    """
    payload = ":".join(str(p) for p in parts).encode()
    digest = hashlib.sha256(payload).digest()
    return int.from_bytes(digest[:8], "big") % (2**32)


# --------------------------------------------------------------------------
# Instance setup and ground truth
# --------------------------------------------------------------------------


@dataclass
class GroundTruth:
    """One instance's exhaustive shell records, keyed by depth.

    Attributes:
        network: Network name.
        x: Treatment.
        y: Outcome.
        r_val: Validity radius (``results/final_table_eps_gt1``).
        n: ``|K_G0|``.
        shells: ``d -> record`` from ``results/mean_vs_max/shells.jsonl``,
            covering every ``d`` in ``0..n``.
    """

    network: str
    x: str
    y: str
    r_val: int
    n: int
    shells: dict[int, dict[str, Any]]

    def profile(self) -> list[ShellStat]:
        """The exact profile as ``ShellStat`` objects, ``d = 0..n``.

        Built directly from the exhaustive ground truth file -- no
        recomputation -- so ``se == 0`` and ``exhaustive is True`` throughout,
        matching what :func:`shell_stat` itself returns when it enumerates.
        """
        out = []
        for d in range(self.n + 1):
            r = self.shells[d]
            out.append(
                ShellStat(
                    d=d,
                    n_subsets_total=r["n_subsets"],
                    n_evaluated=r["n_subsets"],
                    n_states=r["n_states"],
                    mean=r["mean_subsets"],
                    se=0.0,
                    beta_up=r["max"],
                    min_nonzero=0.0,
                    frac_nonzero=r["frac_nonzero_subsets"],
                    exhaustive=True,
                )
            )
        return out

    @property
    def beta_top(self) -> float:
        """``mu(n) = beta_up(n)``: the mean at full retraction, single subset."""
        return self.shells[self.n]["mean_subsets"]


def load_ground_truth(shells_path: Path, rval_path: Path) -> dict[tuple[str, str, str], GroundTruth]:
    """Load the exhaustive shell study and the validity radii, joined by key.

    Args:
        shells_path: ``results/mean_vs_max/shells.jsonl``.
        rval_path: ``results/final_table_eps_gt1/instances.jsonl``.

    Returns:
        ``(network, x, y) -> GroundTruth``, excluding any instance with
        ``r_val == 0`` (degenerate; house rule) or missing an ``r_val``.
    """
    rvals: dict[tuple[str, str, str], int] = {}
    for line in rval_path.open():
        r = json.loads(line)
        rvals[(r["network"], r["x"], r["y"])] = r["r_val"]

    raw: dict[tuple[str, str, str], dict[int, dict[str, Any]]] = {}
    for line in shells_path.open():
        r = json.loads(line)
        key = (r["network"], r["x"], r["y"])
        raw.setdefault(key, {})[r["d"]] = r

    out: dict[tuple[str, str, str], GroundTruth] = {}
    for key, shells in raw.items():
        rv = rvals.get(key)
        if rv is None:
            print(f"[mean_sampling] {key}: no r_val on record, skipping", file=sys.stderr)
            continue
        if rv == 0:
            print(f"[mean_sampling] {key}: r_val == 0 (degenerate), excluded", file=sys.stderr)
            continue
        n = max(shells)
        out[key] = GroundTruth(network=key[0], x=key[1], y=key[2], r_val=rv, n=n, shells=shells)
    return out


@dataclass
class InstanceCtx:
    """Everything needed to draw and evaluate subsets for one instance."""

    gt: GroundTruth
    cpdag: Any
    g0: Any
    ctx: Any
    scale: float
    k0: list[tuple[str, str]]
    n: int


def build_instance(gt: GroundTruth, condition: str, seed: int) -> InstanceCtx | None:
    """Reconstruct the bias context for one instance, per the standard recipe.

    Args:
        gt: Ground truth record for this instance (gives network/x/y).
        condition: Knowledge condition, always ``"D_LLM"`` here.
        seed: Base seed for the per-instance SEM draw.

    Returns:
        An :class:`InstanceCtx`, or ``None`` if the instance is unavailable or
        degenerate under this rebuild (should not happen given ``gt`` already
        passed the same checks in ``mean_vs_max_shell.py``, but checked again
        rather than assumed).
    """
    know = ft.load_knowledge(condition)
    parsed = ft.load_networks({gt.network}, set(), 0)
    if gt.network not in parsed or gt.network not in know:
        return None
    P = parsed[gt.network]
    cpdag, dag = P["cpdag"], P["dag"]
    g0 = apply_orientations(cpdag, list(know[gt.network]))
    if g0 is None:
        return None
    z = optimal_adjustment_set_mpdag(g0, gt.x, gt.y)
    if z is None:
        return None
    rng = np.random.default_rng(ft.derive_seed(seed, gt.network, gt.x, gt.y))
    ctx = make_context(random_sem(dag, rng), cpdag, gt.x, gt.y, frozenset(z))
    scale = abs(ctx.theta_z)
    if scale < 1e-12:
        return None
    k0 = knowledge_of(cpdag, g0)
    if len(k0) != gt.n:
        print(
            f"[mean_sampling] {gt.network}:{gt.x}:{gt.y}: rebuilt |K_G0|={len(k0)} "
            f"!= ground truth n={gt.n}, skipping",
            file=sys.stderr,
        )
        return None
    return InstanceCtx(gt=gt, cpdag=cpdag, g0=g0, ctx=ctx, scale=scale, k0=k0, n=gt.n)


# --------------------------------------------------------------------------
# Raw sampling replica (for bootstrap, which shell_stat's summary can't give)
# --------------------------------------------------------------------------


def draw_raw(inst: InstanceCtx, d: int, m: int, rng: np.random.Generator, cache: dict[str, float]) -> np.ndarray:
    """Draw ``m`` subsets of shell ``d`` and evaluate ``B`` on each.

    Mirrors :func:`bkrobust.epsilon.meanprofile.shell_stat`'s internal loop
    exactly (same sampling call, same cache-by-edge-string, same
    ``bias_at(..., method="semilocal").worst / scale``), so that the returned
    array is what ``shell_stat`` computes internally, just not thrown away.
    Verified once against the real ``shell_stat`` in :func:`verify_replica`.

    Args:
        inst: The instance.
        d: Retraction depth.
        m: Subsets to draw (``>=` population enumerates exhaustively).
        rng: Randomness.
        cache: Edge-string to bias cache, shared across calls for this
            instance.

    Returns:
        Array of ``B`` values, one per (not-necessarily-distinct) subset.
    """
    n = inst.n
    total = math.comb(n, d)
    if m >= total:
        from bkrobust.epsilon.meanprofile import _subsets

        picks = list(_subsets(n, d))
    else:
        picks = _sample_subsets(n, d, m, rng)
    values = np.empty(len(picks), dtype=float)
    k0 = inst.k0
    for i, drop in enumerate(picks):
        dropped = set(drop)
        keep = [k0[j] for j in range(n) if j not in dropped]
        h = apply_orientations(inst.cpdag, keep)
        assert h is not None, "subset of a consistent set must be consistent"
        es = h.edge_string()
        if es not in cache:
            cache[es] = bias_at(inst.ctx, h, method="semilocal").worst / (inst.scale or 1.0)
        values[i] = cache[es]
    return values


def normal_ci(mean: float, arr: np.ndarray, total: int, z: float = Z95) -> tuple[float, float, float]:
    """The module's own normal, finite-population-corrected interval.

    Reproduces ``ShellStat.se`` / ``ShellStat.ci`` exactly (same formula),
    from a raw array rather than a ``ShellStat``, so it can sit next to the
    bootstrap interval built from the same draws.

    Args:
        mean: Sample mean.
        arr: Raw drawn values.
        total: Population size ``C(n, d)``.
        z: Normal quantile.

    Returns:
        ``(se, lo, hi)``.
    """
    m = arr.size
    if m < 2 or m >= total:
        return 0.0, mean, mean
    var = float(arr.var(ddof=1))
    se = math.sqrt(max(0.0, var / m) * max(0.0, 1.0 - m / total))
    return se, max(0.0, mean - z * se), mean + z * se


def bootstrap_ci(arr: np.ndarray, rng: np.random.Generator, b: int = 1000, z_lo: float = 2.5, z_hi: float = 97.5) -> tuple[float, float]:
    """A percentile bootstrap interval for the mean, from a WOR sample.

    Resamples ``arr`` with replacement ``b`` times (the standard bootstrap;
    ignoring the finite-population correction on the resample step is the
    usual approximation and is what is being tested against coverage here).

    Args:
        arr: Raw drawn values.
        rng: Randomness.
        b: Bootstrap resamples.
        z_lo: Lower percentile.
        z_hi: Upper percentile.

    Returns:
        ``(lo, hi)``, clipped at zero below.
    """
    m = arr.size
    if m < 2:
        mean = float(arr.mean()) if m else 0.0
        return mean, mean
    idx = rng.integers(0, m, size=(b, m))
    means = arr[idx].mean(axis=1)
    lo, hi = np.percentile(means, [z_lo, z_hi])
    return max(0.0, float(lo)), float(hi)


# --------------------------------------------------------------------------
# Verification (run once, small scale, before the main study)
# --------------------------------------------------------------------------


def verify_replica(inst: InstanceCtx, d: int, m: int, seed: int) -> dict[str, Any]:
    """Confirm the hand-rolled raw-draw loop reproduces ``shell_stat`` exactly.

    Both calls start from ``np.random.default_rng(seed)`` with the same seed;
    ``_sample_subsets`` is deterministic given the RNG state, so the two
    should draw the identical subsets and land on the identical mean and SE.

    Args:
        inst: The instance.
        d: Shell depth to check.
        m: Sample size.
        seed: Shared seed for both draws.

    Returns:
        A record with both means/SEs and whether they matched.
    """
    cache_a: dict[str, float] = {}
    official = shell_stat(inst.ctx, inst.cpdag, inst.g0, d, scale=inst.scale, sample=m, rng=np.random.default_rng(seed), cache=cache_a)
    cache_b: dict[str, float] = {}
    arr = draw_raw(inst, d, m, np.random.default_rng(seed), cache_b)
    se, lo, hi = normal_ci(float(arr.mean()), arr, math.comb(inst.n, d))
    ok = (
        official.n_evaluated == arr.size
        and math.isclose(official.mean, float(arr.mean()), rel_tol=1e-9, abs_tol=1e-12)
        and math.isclose(official.se, se, rel_tol=1e-9, abs_tol=1e-12)
    )
    return {
        "network": inst.gt.network, "x": inst.gt.x, "y": inst.gt.y, "d": d, "m": m,
        "official_mean": official.mean, "replica_mean": float(arr.mean()),
        "official_se": official.se, "replica_se": se, "match": bool(ok),
    }


def verify_fpc(inst: InstanceCtx, d: int, sample_sizes: tuple[int, ...], reps: int, seed: int) -> dict[str, Any]:
    """Check the finite-population-corrected SE against the true WOR sampling variance.

    Fully enumerates shell ``d`` (cheap for the shell this is called on -- pick
    one with a modest population), giving the exact population mean and
    variance. For each sample size, draws ``reps`` independent without-
    replacement samples and compares the *empirical* standard deviation of the
    resulting sample means to what ``shell_stat``'s formula predicts
    (``sqrt(S^2/m * (1 - m/N))`` with ``S^2`` the exact population variance,
    ``ddof=1``, matching the module).

    Args:
        inst: The instance.
        d: Shell to fully enumerate.
        sample_sizes: Sample sizes to check.
        reps: Independent WOR draws per sample size.
        seed: Base seed.

    Returns:
        A record with, per sample size, the predicted vs. empirical SE of the
        mean and their ratio.
    """
    cache: dict[str, float] = {}
    total = math.comb(inst.n, d)
    full = draw_raw(inst, d, total, np.random.default_rng(seed), cache)
    pop_mean = float(full.mean())
    pop_var = float(full.var(ddof=1))  # sample variance of the *population* of subsets
    rows = []
    for m in sample_sizes:
        if m >= total:
            continue
        predicted_se = math.sqrt(max(0.0, pop_var / m) * max(0.0, 1.0 - m / total))
        means = np.empty(reps)
        for r in range(reps):
            arr = draw_raw(inst, d, m, np.random.default_rng(seed_for(seed, d, m, r)), cache)
            means[r] = arr.mean()
        empirical_se = float(means.std(ddof=1))
        rows.append(
            {
                "m": m, "predicted_se": predicted_se, "empirical_se": empirical_se,
                "ratio": empirical_se / predicted_se if predicted_se > 0 else None,
                "empirical_bias": float(means.mean()) - pop_mean,
            }
        )
    return {
        "network": inst.gt.network, "x": inst.gt.x, "y": inst.gt.y, "d": d,
        "n_subsets_total": total, "pop_mean": pop_mean, "pop_var": pop_var, "reps": reps,
        "sample_sizes": rows,
    }


# --------------------------------------------------------------------------
# Part 1+2: per-shell estimator accuracy and CI coverage
# --------------------------------------------------------------------------


def reps_for(sample_size: int, base_reps: int, cap_draws: int, min_reps: int) -> int:
    """Scale repetitions down for expensive cells, up to ``base_reps``.

    Keeps ``reps * sample_size`` (the draws spent on one cell) under
    ``cap_draws``, so a handful of very large shells or budgets do not blow up
    the wall-clock budget; documents the actual count used in every output
    row rather than silently reducing precision.

    Args:
        sample_size: Draws per repetition.
        base_reps: Repetitions when affordable (the task's ">= 200").
        cap_draws: Total draws this cell may spend.
        min_reps: Floor, so coverage is still estimable.

    Returns:
        Repetitions to run, in ``[min_reps, base_reps]`` (or ``0`` if
        ``sample_size`` is ``0``).
    """
    if sample_size <= 0:
        return 0
    return int(max(min_reps, min(base_reps, cap_draws // sample_size)))


def select_shells_for_shell_study(inst: InstanceCtx, min_pop: int, max_per_instance: int) -> list[int]:
    """Pick up to two representative shells with enough population to sample.

    Chooses the smallest- and largest-population shell (in ``d in [r_val, n]``)
    whose population is at least ``min_pop``; a shell smaller than that would
    make every requested fraction/budget trivially exhaustive.

    Args:
        inst: The instance.
        min_pop: Minimum ``C(n, d)`` to qualify.
        max_per_instance: Cap (2 gives min- and max-population shells).

    Returns:
        Sorted, deduplicated list of depths.
    """
    candidates = [(d, math.comb(inst.n, d)) for d in range(inst.gt.r_val, inst.n + 1)]
    candidates = [(d, pop) for d, pop in candidates if pop >= min_pop]
    if not candidates:
        return []
    by_pop = sorted(candidates, key=lambda t: t[1])
    picked = {by_pop[0][0], by_pop[-1][0]}
    return sorted(picked)[:max_per_instance] if max_per_instance < 2 else sorted(picked)


def run_shell_cell(
    inst: InstanceCtx, d: int, cell_type: str, cell_value: float, sample_size: int,
    reps: int, base_seed: int, boot_b: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Run ``reps`` independent samples of shell ``d`` at a fixed sample size.

    Args:
        inst: The instance.
        d: Shell depth.
        cell_type: ``"fraction"`` or ``"budget"``.
        cell_value: The fraction or budget defining ``sample_size``.
        sample_size: Subsets drawn per repetition.
        reps: Repetitions.
        base_seed: Seeds every repetition deterministically via :func:`seed_for`.
        boot_b: Bootstrap resamples per repetition.

    Returns:
        ``(per_rep_rows, summary_row)``.
    """
    true_mu = inst.gt.shells[d]["mean_subsets"]
    total = math.comb(inst.n, d)
    cache: dict[str, float] = {}
    rows = []
    for rep in range(reps):
        seed = seed_for(base_seed, inst.gt.network, inst.gt.x, inst.gt.y, d, cell_type, cell_value, rep)
        rng = np.random.default_rng(seed)
        arr = draw_raw(inst, d, sample_size, rng, cache)
        mean_hat = float(arr.mean())
        se, lo_n, hi_n = normal_ci(mean_hat, arr, total)
        lo_b, hi_b = bootstrap_ci(arr, np.random.default_rng(seed_for(seed, "boot")), b=boot_b)
        rows.append(
            {
                "network": inst.gt.network, "x": inst.gt.x, "y": inst.gt.y, "d": d,
                "cell_type": cell_type, "cell_value": cell_value, "sample_size": sample_size,
                "n_subsets_total": total, "rep": rep, "true_mu": true_mu,
                "mean_hat": mean_hat, "se_normal": se, "ci_lo_normal": lo_n, "ci_hi_normal": hi_n,
                "ci_lo_boot": lo_b, "ci_hi_boot": hi_b,
                "covered_normal": bool(lo_n <= true_mu <= hi_n),
                "covered_boot": bool(lo_b <= true_mu <= hi_b),
                "width_normal": hi_n - lo_n, "width_boot": hi_b - lo_b,
                "error": mean_hat - true_mu,
            }
        )
    errs = np.array([r["error"] for r in rows])
    summary = {
        "network": inst.gt.network, "x": inst.gt.x, "y": inst.gt.y, "d": d,
        "cell_type": cell_type, "cell_value": cell_value, "sample_size": sample_size,
        "n_subsets_total": total, "true_mu": true_mu, "n_reps": reps,
        "bias": float(errs.mean()), "rmse": float(np.sqrt((errs**2).mean())),
        "coverage_normal": float(np.mean([r["covered_normal"] for r in rows])),
        "coverage_boot": float(np.mean([r["covered_boot"] for r in rows])),
        "mean_width_normal": float(np.mean([r["width_normal"] for r in rows])),
        "mean_width_boot": float(np.mean([r["width_boot"] for r in rows])),
    }
    return rows, summary


# --------------------------------------------------------------------------
# Part 3+4: radius fidelity
# --------------------------------------------------------------------------


def classify_radius(true_r: int, hat_r: int) -> tuple[str, int | None]:
    """Classify one radius comparison, handling ``UNREACHED`` as a status.

    ``r_mean`` is larger when the profile takes longer to cross ``eps`` --
    i.e. a larger radius claims robustness extends further. ``UNREACHED``
    (never crosses) is the extreme of "large". So:

    * both UNREACHED -> ``"exact_unreached"``.
    * true finite, hat finite, equal -> ``"exact"``.
    * true finite, hat finite, hat > true -> ``"overstate"`` (dangerous: claims
      more robustness than is actually there).
    * true finite, hat finite, hat < true -> ``"understate"`` (safe: claims
      less robustness than is actually there).
    * true finite, hat UNREACHED -> ``"overstate_unreached"`` (the most
      dangerous case: claims the bias never crosses eps when it does).
    * true UNREACHED, hat finite -> ``"understate_unreached_true"`` (safe: a
      finite radius reported where none was needed).

    Args:
        true_r: Radius from the exhaustive profile (may be ``UNREACHED``).
        hat_r: Radius from the sampled profile (may be ``UNREACHED``).

    Returns:
        ``(status, signed_error_in_shells)``. The error is ``None`` whenever
        either side is ``UNREACHED`` -- it is never given a numeric value.
    """
    if true_r == UNREACHED and hat_r == UNREACHED:
        return "exact_unreached", None
    if true_r == UNREACHED and hat_r != UNREACHED:
        return "understate_unreached_true", None
    if true_r != UNREACHED and hat_r == UNREACHED:
        return "overstate_unreached", None
    err = hat_r - true_r
    if err == 0:
        return "exact", 0
    return ("overstate" if err > 0 else "understate"), err


def eps_grid_for(inst: InstanceCtx) -> list[float]:
    """The per-instance epsilon grid: fractions of ``mu`` at full retraction.

    Args:
        inst: The instance.

    Returns:
        ``len(EPS_QUANTILES)`` epsilon values, sorted, deduplicated.
    """
    top = inst.gt.beta_top
    vals = sorted({round(q * top, 10) for q in EPS_QUANTILES})
    return [v for v in vals if v > 0]


def run_radius_cell(
    inst: InstanceCtx, cell_type: str, cell_value: float, reps: int, base_seed: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Build ``reps`` full sampled profiles at one fraction/budget and score radii.

    For each repetition, every shell ``d in [r_val, n]`` is sampled at
    ``sample_size(d)`` derived from ``cell_type``/``cell_value`` (a shell
    smaller than that size is enumerated exactly, matching ``shell_stat``'s
    own fallback). ``r_mean`` is then read off the resulting profile under
    ``use in {"mean", "ci_lo", "ci_hi"}`` and compared to the true radius for
    every epsilon in :func:`eps_grid_for`.

    Args:
        inst: The instance.
        cell_type: ``"fraction"`` or ``"budget"``.
        cell_value: Fraction or fixed budget.
        reps: Independent full-profile repetitions.
        base_seed: Seeds every repetition.

    Returns:
        ``(per_rep_per_eps_rows, per_rep_profile_rows)`` -- the second gives
        one row per repetition with the achieved sample sizes, for auditing.
    """
    true_profile = inst.gt.profile()
    eps_grid = eps_grid_for(inst)
    true_radii = {eps: r_mean(true_profile, eps) for eps in eps_grid}
    ds = list(range(inst.gt.r_val, inst.n + 1))

    cache: dict[str, float] = {}
    eps_rows: list[dict[str, Any]] = []
    profile_rows: list[dict[str, Any]] = []
    for rep in range(reps):
        sampled: list[ShellStat] = list(true_profile[: inst.gt.r_val])  # zero shells, exact, free
        total_draws = 0
        for d in ds:
            total = math.comb(inst.n, d)
            if cell_type == "fraction":
                m = max(2, round(cell_value * total))
            else:
                m = int(cell_value)
            m = min(m, total)
            seed = seed_for(base_seed, inst.gt.network, inst.gt.x, inst.gt.y, cell_type, cell_value, rep, d)
            arr = draw_raw(inst, d, m, np.random.default_rng(seed), cache)
            mean_hat = float(arr.mean())
            se, lo, hi = normal_ci(mean_hat, arr, total)
            sampled.append(
                ShellStat(
                    d=d, n_subsets_total=total, n_evaluated=arr.size, n_states=0,
                    mean=mean_hat, se=se, beta_up=float(arr.max()) if arr.size else 0.0,
                    min_nonzero=0.0, frac_nonzero=float((arr > 0).mean()) if arr.size else 0.0,
                    exhaustive=(arr.size >= total),
                )
            )
            total_draws += arr.size
        profile_rows.append(
            {
                "network": inst.gt.network, "x": inst.gt.x, "y": inst.gt.y,
                "cell_type": cell_type, "cell_value": cell_value, "rep": rep,
                "total_draws": total_draws,
            }
        )
        for eps in eps_grid:
            true_r = true_radii[eps]
            for use in ("mean", "ci_lo", "ci_hi"):
                hat_r = r_mean(sampled, eps, use=use)
                status, err = classify_radius(true_r, hat_r)
                eps_rows.append(
                    {
                        "network": inst.gt.network, "x": inst.gt.x, "y": inst.gt.y,
                        "cell_type": cell_type, "cell_value": cell_value, "rep": rep,
                        "eps": eps, "use": use, "true_r": true_r, "hat_r": hat_r,
                        "status": status, "error_shells": err,
                    }
                )
    return eps_rows, profile_rows


def summarize_radius(eps_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Aggregate per-rep radius comparisons into rates, grouped by every axis but rep.

    Args:
        eps_rows: Rows from (potentially many) :func:`run_radius_cell` calls.

    Returns:
        One summary row per ``(network,x,y,cell_type,cell_value,eps,use)``.
    """
    groups: dict[tuple, list[dict[str, Any]]] = {}
    for r in eps_rows:
        key = (r["network"], r["x"], r["y"], r["cell_type"], r["cell_value"], r["eps"], r["use"])
        groups.setdefault(key, []).append(r)
    out = []
    for key, rs in groups.items():
        network, x, y, cell_type, cell_value, eps, use = key
        n = len(rs)
        statuses = [r["status"] for r in rs]
        finite_errs = [r["error_shells"] for r in rs if r["error_shells"] is not None]
        out.append(
            {
                "network": network, "x": x, "y": y, "cell_type": cell_type, "cell_value": cell_value,
                "eps": eps, "use": use, "n_reps": n,
                "true_r": rs[0]["true_r"],
                "rate_exact": statuses.count("exact") / n,
                "rate_exact_unreached": statuses.count("exact_unreached") / n,
                "rate_overstate": statuses.count("overstate") / n,
                "rate_overstate_unreached": statuses.count("overstate_unreached") / n,
                "rate_understate": statuses.count("understate") / n,
                "rate_understate_unreached_true": statuses.count("understate_unreached_true") / n,
                "rate_dangerous": (statuses.count("overstate") + statuses.count("overstate_unreached")) / n,
                "mean_abs_error_shells": float(np.mean(np.abs(finite_errs))) if finite_errs else None,
                "n_finite_both": len(finite_errs),
            }
        )
    return out


def summarize_adaptive(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Aggregate per-rep adaptive-rule outcomes, grouped by ``(instance, eps)``.

    Args:
        rows: Rows from (potentially many) :func:`run_adaptive_instance` calls.

    Returns:
        One summary row per ``(network,x,y,eps)``, with rates in the same
        vocabulary as :func:`summarize_radius` plus the average total draws
        spent finding the radius.
    """
    groups: dict[tuple, list[dict[str, Any]]] = {}
    for r in rows:
        key = (r["network"], r["x"], r["y"], r["eps"])
        groups.setdefault(key, []).append(r)
    out = []
    for key, rs in groups.items():
        network, x, y, eps = key
        n = len(rs)
        statuses = [r["status"] for r in rs]
        finite_errs = [r["error_shells"] for r in rs if r["error_shells"] is not None]
        out.append(
            {
                "network": network, "x": x, "y": y, "eps": eps, "n_reps": n,
                "true_r": rs[0]["true_r"],
                "rate_exact": statuses.count("exact") / n,
                "rate_exact_unreached": statuses.count("exact_unreached") / n,
                "rate_dangerous": (statuses.count("overstate") + statuses.count("overstate_unreached")) / n,
                "rate_understate": statuses.count("understate") / n,
                "mean_abs_error_shells": float(np.mean(np.abs(finite_errs))) if finite_errs else None,
                "mean_total_draws": float(np.mean([r["total_draws"] for r in rs])),
            }
        )
    return out


# --------------------------------------------------------------------------
# Part 5: exploratory adaptive/sequential sampling
# --------------------------------------------------------------------------


def adaptive_shell_decision(
    inst: InstanceCtx, d: int, eps: float, rng: np.random.Generator, cache: dict[str, float],
    step: int, cap: int,
) -> dict[str, Any]:
    """Sample shell ``d`` at growing sample sizes until the CI resolves ``mu(d)`` vs ``eps``.

    Checks at cumulative sizes ``step, 2*step, 3*step, ...`` up to ``cap``
    (clipped to the population). Each checkpoint is drawn as a **fresh**
    without-replacement sample of that cumulative size -- not the previous
    checkpoint's draws plus a new increment -- because ``_sample_subsets``
    only guarantees distinctness *within* one call; concatenating independent
    increments can double-count a subset drawn twice, which would silently
    corrupt the "``arr.size >= total``, so this is exact" check the ``eps``
    decision leans on at the end of a shell's population. Redrawing the full
    checkpoint size each time costs some repeated ``apply_orientations``
    calls (bias itself is still cached by edge string) but is unambiguously a
    valid WOR sample of that size, so the exactness check is trustworthy.

    After each checkpoint, resolves if ``ci_hi < eps`` (confidently below --
    decide "does not cross") or ``ci_lo > eps`` (confidently above -- decide
    "crosses"); a checkpoint that reaches the full population is exact by
    construction and is decided from the point estimate directly. Otherwise
    continues to the next checkpoint, capped at ``cap``.

    Args:
        inst: The instance.
        d: Shell depth.
        eps: Threshold.
        rng: Randomness (one ``rng.integers`` draw seeds each checkpoint's
            fresh sample, so the sequence is still reproducible from ``rng``).
        cache: Shared bias cache.
        step: Checkpoint spacing.
        cap: Largest checkpoint size before forcing a decision from the
            point estimate.

    Returns:
        A record with the decision, draws spent at the resolving checkpoint,
        and whether it was capped rather than resolved.
    """
    total = math.comb(inst.n, d)
    ceiling = min(cap, total)
    target = 0
    while target < ceiling:
        target = min(target + step, ceiling)
        checkpoint_seed = int(rng.integers(0, 2**32))
        arr = draw_raw(inst, d, target, np.random.default_rng(checkpoint_seed), cache)
        mean_hat = float(arr.mean())
        se, lo, hi = normal_ci(mean_hat, arr, total)
        if arr.size >= total or hi < eps or lo > eps:
            crosses = mean_hat > eps if arr.size >= total else (lo > eps)
            return {
                "d": d, "eps": eps, "decision": bool(crosses), "n_drawn": int(arr.size),
                "capped": arr.size < total and not (hi < eps or lo > eps), "mean_hat": mean_hat,
            }
    # Hit the cap without resolving: fall back to the point-estimate decision
    # from the last (largest) checkpoint drawn.
    return {"d": d, "eps": eps, "decision": bool(mean_hat > eps), "n_drawn": int(arr.size), "capped": True, "mean_hat": mean_hat}


def run_adaptive_instance(inst: InstanceCtx, eps: float, reps: int, base_seed: int, step: int, cap: int) -> list[dict[str, Any]]:
    """Run the adaptive radius-finding rule ``reps`` times for one (instance, eps).

    Walks shells from ``r_val`` upward, using :func:`adaptive_shell_decision`
    at each, stopping at the first shell decided to cross (or exhausting the
    profile -> UNREACHED). Records total draws spent, for comparison against
    a fixed budget's accuracy at a matched total draw count.

    Args:
        inst: The instance.
        eps: Threshold.
        reps: Repetitions.
        base_seed: Seeds every repetition.
        step: Increment size.
        cap: Per-shell cap.

    Returns:
        One row per repetition.
    """
    true_profile = inst.gt.profile()
    true_r = r_mean(true_profile, eps)
    ds = list(range(inst.gt.r_val, inst.n + 1))
    rows = []
    for rep in range(reps):
        cache: dict[str, float] = {}
        total_draws = 0
        hat_r = UNREACHED
        for d in ds:
            seed = seed_for(base_seed, inst.gt.network, inst.gt.x, inst.gt.y, "adaptive", eps, rep, d)
            rec = adaptive_shell_decision(inst, d, eps, np.random.default_rng(seed), cache, step, cap)
            total_draws += rec["n_drawn"]
            if rec["decision"]:
                hat_r = d
                break
        status, err = classify_radius(true_r, hat_r)
        rows.append(
            {
                "network": inst.gt.network, "x": inst.gt.x, "y": inst.gt.y, "eps": eps, "rep": rep,
                "true_r": true_r, "hat_r": hat_r, "status": status, "error_shells": err,
                "total_draws": total_draws,
            }
        )
    return rows


# --------------------------------------------------------------------------
# Per-instance worker (multiprocessing unit)
# --------------------------------------------------------------------------


def process_instance(args: dict[str, Any]) -> dict[str, Any]:
    """Everything for one instance: verification (if selected), Parts 1-3, and Part 5.

    Runs in its own process; the bias cache is per-instance and never shared
    across processes, which is fine since states rarely coincide across
    instances and the cache only saves work within one.

    Args:
        args: Bundled configuration (see :func:`main`), plus the
            ``GroundTruth`` (serialized as a dict for ``multiprocessing``).

    Returns:
        Dict of row lists: ``verification``, ``shell_reps``, ``shell_summary``,
        ``radius_reps``, ``radius_summary``, ``profile_audit``,
        ``adaptive_reps``.
    """
    t0 = time.perf_counter()
    gt: GroundTruth = args["gt"]
    inst = build_instance(gt, args["condition"], args["seed"])
    key = f"{gt.network}:{gt.x}:{gt.y}"
    if inst is None:
        print(f"[mean_sampling] {key}: could not rebuild instance, skipping", file=sys.stderr)
        return {
            "verification": [], "shell_reps": [], "shell_summary": [],
            "radius_reps": [], "radius_summary": [], "profile_audit": [], "adaptive_reps": [],
            "adaptive_summary": [],
        }

    out: dict[str, list[dict[str, Any]]] = {
        "verification": [], "shell_reps": [], "shell_summary": [],
        "radius_reps": [], "radius_summary": [], "profile_audit": [], "adaptive_reps": [],
        "adaptive_summary": [],
    }

    # -- Verification (only for the instances chosen to carry it) -------
    if args["verify"]:
        d_mid = max(inst.gt.r_val, inst.n // 2)
        pop_mid = math.comb(inst.n, d_mid)
        m_check = max(2, min(pop_mid - 1, 50))
        out["verification"].append(verify_replica(inst, d_mid, m_check, seed_for(args["seed"], "verify_replica", key)))
        if pop_mid <= args["fpc_max_pop"]:
            sizes = tuple(sorted({s for s in (5, 10, 25, 50, pop_mid // 2) if 1 < s < pop_mid}))
            out["verification"].append(
                {"fpc_check": verify_fpc(inst, d_mid, sizes, args["fpc_reps"], seed_for(args["seed"], "verify_fpc", key))}
            )

    # -- Part 1+2: per-shell accuracy and coverage -----------------------
    shells = select_shells_for_shell_study(inst, args["min_pop"], args["max_shells_per_instance"])
    for d in shells:
        total = math.comb(inst.n, d)
        for frac in args["fractions"]:
            m = max(2, round(frac * total))
            if m >= total:
                continue
            reps = reps_for(m, args["reps_shell"], args["cap_draws_shell"], args["min_reps"])
            rows, summ = run_shell_cell(inst, d, "fraction", frac, m, reps, args["seed"], args["boot_b"])
            out["shell_reps"].extend(rows)
            out["shell_summary"].append(summ)
        for budget in args["budgets"]:
            if budget >= total:
                continue
            reps = reps_for(budget, args["reps_shell"], args["cap_draws_shell"], args["min_reps"])
            rows, summ = run_shell_cell(inst, d, "budget", budget, budget, reps, args["seed"], args["boot_b"])
            out["shell_reps"].extend(rows)
            out["shell_summary"].append(summ)

    # -- Part 3+4: radius fidelity, full profile -------------------------
    total_profile_pop = sum(math.comb(inst.n, d) for d in range(inst.gt.r_val, inst.n + 1))
    for cell_type, grid in (("fraction", args["fractions"]), ("budget", args["budgets"])):
        for cell_value in grid:
            if cell_type == "fraction":
                approx_draws_per_rep = round(cell_value * total_profile_pop)
            else:
                approx_draws_per_rep = sum(min(cell_value, math.comb(inst.n, d)) for d in range(inst.gt.r_val, inst.n + 1))
            reps = reps_for(max(1, approx_draws_per_rep), args["reps_profile"], args["cap_draws_profile"], args["min_reps_profile"])
            eps_rows, prof_rows = run_radius_cell(inst, cell_type, cell_value, reps, args["seed"])
            out["radius_reps"].extend(eps_rows)
            out["profile_audit"].extend(prof_rows)

    # -- Part 5: exploratory adaptive sampling ---------------------------
    if args["adaptive_eps_quantiles"]:
        eps_grid = eps_grid_for(inst)
        chosen = [eps_grid[i] for i in args["adaptive_eps_quantiles"] if i < len(eps_grid)]
        for eps in chosen:
            out["adaptive_reps"].extend(
                run_adaptive_instance(inst, eps, args["adaptive_reps"], args["seed"], args["adaptive_step"], args["adaptive_cap"])
            )

    print(
        f"[mean_sampling] {key} done in {time.perf_counter() - t0:.1f}s: "
        f"shell_reps={len(out['shell_reps'])} radius_reps={len(out['radius_reps'])} "
        f"adaptive_reps={len(out['adaptive_reps'])}",
        flush=True,
    )
    return out


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------


def write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write ``rows`` as one JSON object per line."""
    with path.open("w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")


def main() -> None:
    """Run the sampled mu(d) study and write all output files."""
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--condition", default="D_LLM")
    p.add_argument("--instances", default=",".join(DEFAULT_INSTANCES), help="net:x:y,net:x:y,...")
    p.add_argument("--seed", type=int, default=20260917)
    p.add_argument("--shells-path", default="results/mean_vs_max/shells.jsonl")
    p.add_argument("--rval-path", default="results/final_table_eps_gt1/instances.jsonl")
    p.add_argument("--out", default="results/mean_sampling")
    p.add_argument("--fractions", default=",".join(str(f) for f in DEFAULT_FRACTIONS))
    p.add_argument("--budgets", default=",".join(str(b) for b in DEFAULT_BUDGETS))
    p.add_argument("--reps-shell", type=int, default=200, help="repetitions per shell cell when affordable")
    p.add_argument("--cap-draws-shell", type=int, default=35_000, help="max draws spent on one shell cell")
    p.add_argument("--min-reps", type=int, default=30, help="floor on shell-cell repetitions")
    p.add_argument("--reps-profile", type=int, default=100, help="repetitions per radius cell when affordable")
    p.add_argument("--cap-draws-profile", type=int, default=70_000, help="max draws spent on one radius cell")
    p.add_argument("--min-reps-profile", type=int, default=20, help="floor on radius-cell repetitions")
    p.add_argument("--min-pop", type=int, default=20, help="min shell population to include in the per-shell study")
    p.add_argument("--max-shells-per-instance", type=int, default=2)
    p.add_argument("--boot-b", type=int, default=1000, help="bootstrap resamples")
    p.add_argument("--fpc-max-pop", type=int, default=2000, help="skip the FPC verification's full enumeration above this population")
    p.add_argument("--fpc-reps", type=int, default=500)
    p.add_argument(
        "--verify-keys",
        default="asia:bronc:dysp,Schipf_2010:A:TT,barley:dg25:s2225,child:CO2:DuctFlow",
        help="net:x:y,... to carry the verification checks (small + a couple with real population)",
    )
    p.add_argument("--adaptive-eps-quantiles", default="1,2,3", help="indices into each instance's eps grid to run adaptively, comma-separated, empty to skip")
    p.add_argument("--adaptive-reps", type=int, default=60)
    p.add_argument("--adaptive-step", type=int, default=25)
    p.add_argument("--adaptive-cap", type=int, default=500)
    p.add_argument(
        "--adaptive-keys",
        default="barley:dg25:s2225,magic-irri:G3212:FT,child:CO2:DuctFlow,ecoli70:cspG:hupB,paths:11:14",
        help="net:x:y,... to run the exploratory adaptive/sequential study on",
    )
    p.add_argument("--jobs", type=int, default=0, help="worker processes; 0 = one per instance, capped at cpu_count")
    args = p.parse_args()

    out_dir = ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    wanted_keys = [tuple(s.split(":", 2)) for s in args.instances.split(",") if s]
    gt_all = load_ground_truth(ROOT / args.shells_path, ROOT / args.rval_path)
    gts = [gt_all[k] for k in wanted_keys if k in gt_all]
    missing = [k for k in wanted_keys if k not in gt_all]
    if missing:
        print(f"[mean_sampling] no ground truth for {missing}, skipped", file=sys.stderr)

    fractions = tuple(float(s) for s in args.fractions.split(",") if s)
    budgets = tuple(int(s) for s in args.budgets.split(",") if s)
    adaptive_eps_idx = [int(s) for s in args.adaptive_eps_quantiles.split(",") if s]
    verify_keys = {tuple(s.split(":", 2)) for s in args.verify_keys.split(",") if s}
    adaptive_keys = {tuple(s.split(":", 2)) for s in args.adaptive_keys.split(",") if s}

    job_specs = []
    for gt in gts:
        key = (gt.network, gt.x, gt.y)
        job_specs.append(
            {
                "gt": gt, "condition": args.condition, "seed": args.seed,
                "verify": key in verify_keys, "fpc_max_pop": args.fpc_max_pop, "fpc_reps": args.fpc_reps,
                "min_pop": args.min_pop, "max_shells_per_instance": args.max_shells_per_instance,
                "fractions": fractions, "budgets": budgets,
                "reps_shell": args.reps_shell, "cap_draws_shell": args.cap_draws_shell, "min_reps": args.min_reps,
                "reps_profile": args.reps_profile, "cap_draws_profile": args.cap_draws_profile,
                "min_reps_profile": args.min_reps_profile, "boot_b": args.boot_b,
                "adaptive_eps_quantiles": adaptive_eps_idx if key in adaptive_keys else [],
                "adaptive_reps": args.adaptive_reps, "adaptive_step": args.adaptive_step, "adaptive_cap": args.adaptive_cap,
            }
        )

    n_jobs = args.jobs or min(len(job_specs), mp.cpu_count())
    print(f"[mean_sampling] {len(job_specs)} instances, {n_jobs} worker processes", flush=True)

    all_out: dict[str, list[dict[str, Any]]] = {
        "verification": [], "shell_reps": [], "shell_summary": [],
        "radius_reps": [], "radius_summary": [], "profile_audit": [], "adaptive_reps": [],
        "adaptive_summary": [],
    }
    t0 = time.perf_counter()
    if n_jobs <= 1:
        results = [process_instance(spec) for spec in job_specs]
    else:
        with mp.get_context("spawn").Pool(n_jobs) as pool:
            results = pool.map(process_instance, job_specs)
    for r in results:
        for k in all_out:
            all_out[k].extend(r[k])

    all_out["radius_summary"] = summarize_radius(all_out["radius_reps"])
    all_out["adaptive_summary"] = summarize_adaptive(all_out["adaptive_reps"])

    write_jsonl(out_dir / "shell_reps.jsonl", all_out["shell_reps"])
    write_jsonl(out_dir / "shell_summary.jsonl", all_out["shell_summary"])
    write_jsonl(out_dir / "radius_reps.jsonl", all_out["radius_reps"])
    write_jsonl(out_dir / "radius_summary.jsonl", all_out["radius_summary"])
    write_jsonl(out_dir / "profile_audit.jsonl", all_out["profile_audit"])
    write_jsonl(out_dir / "adaptive_reps.jsonl", all_out["adaptive_reps"])
    write_jsonl(out_dir / "adaptive_summary.jsonl", all_out["adaptive_summary"])
    (out_dir / "verification.json").write_text(json.dumps(all_out["verification"], indent=2) + "\n")

    # Top-level summary: coverage by cell type/value, averaged over shells and
    # instances present (each row already an (instance,shell,cell) mean, so
    # this is a mean of means -- fine as an orientation number; the report
    # must still cite the per-query rows for any real claim).
    def _grid_avg(rows: list[dict[str, Any]], key: str, field: str) -> dict[str, float]:
        by_key: dict[Any, list[float]] = {}
        for r in rows:
            by_key.setdefault(r[key], []).append(r[field])
        return {str(k): float(np.mean(v)) for k, v in sorted(by_key.items())}

    summary = {
        "args": {**vars(args), "fractions": fractions, "budgets": budgets},
        "n_instances": len(gts),
        "instances": [f"{g.network}:{g.x}:{g.y}" for g in gts],
        "wall_seconds": round(time.perf_counter() - t0, 1),
        "n_verification_checks": len(all_out["verification"]),
        "n_verification_mismatches": sum(
            1 for v in all_out["verification"] if "match" in v and not v["match"]
        ),
        "shell_study": {
            "n_cells": len(all_out["shell_summary"]),
            "coverage_normal_by_fraction": _grid_avg(
                [r for r in all_out["shell_summary"] if r["cell_type"] == "fraction"], "cell_value", "coverage_normal"
            ),
            "coverage_boot_by_fraction": _grid_avg(
                [r for r in all_out["shell_summary"] if r["cell_type"] == "fraction"], "cell_value", "coverage_boot"
            ),
            "coverage_normal_by_budget": _grid_avg(
                [r for r in all_out["shell_summary"] if r["cell_type"] == "budget"], "cell_value", "coverage_normal"
            ),
            "coverage_boot_by_budget": _grid_avg(
                [r for r in all_out["shell_summary"] if r["cell_type"] == "budget"], "cell_value", "coverage_boot"
            ),
        },
        "radius_study": {
            "n_cells": len(all_out["radius_summary"]),
            "rate_dangerous_by_fraction_use": {
                use: _grid_avg(
                    [r for r in all_out["radius_summary"] if r["cell_type"] == "fraction" and r["use"] == use],
                    "cell_value", "rate_dangerous",
                )
                for use in ("mean", "ci_lo", "ci_hi")
            },
            "rate_exact_by_fraction_use": {
                use: _grid_avg(
                    [r for r in all_out["radius_summary"] if r["cell_type"] == "fraction" and r["use"] == use],
                    "cell_value", "rate_exact",
                )
                for use in ("mean", "ci_lo", "ci_hi")
            },
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"[mean_sampling] wrote outputs to {out_dir}, total wall time {time.perf_counter() - t0:.1f}s")


if __name__ == "__main__":
    main()
