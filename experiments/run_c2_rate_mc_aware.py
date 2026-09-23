"""How often is the incumbent mean-based bias non-monotone, *after* removing Monte-Carlo noise?

The first pass at this comparison (``results/epsilon/counterexamples/c2_rate_sweep.json``)
counted an ordered pair ``G <= H`` as violating monotonicity whenever
``mean_abs_bias(G) > mean_abs_bias(H)`` by more than ``1e-6``, at 80 coefficient
draws per state. That threshold is far below the Monte-Carlo standard error of
an 80-draw mean, so the resulting rate mixes two different claims: that the
*population* mean is non-monotone, and that an 80-draw *estimate* of it happens
to be. Both are bad news for a radius built on the quantity, but they are not the
same finding, and only the first is a statement about the estimand.

This script separates them. For every comparable pair it draws the mean twice,
under independent seeds, at a large draw count, and counts a violation only when

    mean(G) - mean(H)  >  Z_MULTIPLIER * combined MC standard error

in **both** repetitions. That is a deliberately conservative rule: it will miss
real violations whose margin is small, so the rate it reports is a *lower bound*
on the population violation rate -- which is the safe direction for a claim that
the incumbent is broken.

The same pairs are scored for ``B``, which has no Monte-Carlo error at all
(it is a closed-form maximum at a fixed covariance) and is tested at ``1e-9``.
Running both on identical pairs is the point: two separate sweeps would leave
open that the instances differed.
"""

from __future__ import annotations

import json
import platform
import subprocess
import time
from pathlib import Path

import numpy as np

from bkrobust.demo.evaluate import random_sem
from bkrobust.demo.graph import MPDAG
from bkrobust.epsilon.bias import BiasContext, bias_at, make_context
from bkrobust.epsilon.counterexamples import (
    _mean_abs_bias,
    enumerate_comparable_pairs,
    mc_standard_error,
    sweep_candidates,
)

#: How many combined standard errors a deficit must exceed to count. 3 is
#: conservative for a two-sided question asked twice; the rate is reported as a
#: lower bound, so erring high here is the safe direction.
Z_MULTIPLIER = 3.0

#: Draws per state per repetition. Large enough that the standard error is small
#: relative to the margins the flagship witness exhibits (which run at 11-13
#: combined SEs), so the rule is not simply rejecting everything.
N_DRAWS = 300

N_REPETITIONS = 2
ROOT_SEED = 20260916
MAX_SPACE = 24
MAX_PAIRS_PER_INSTANCE = 120
TARGET_INSTANCES = 45

OUT = Path("results/epsilon/counterexamples")


def _git_commit() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True
        ).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def main() -> None:
    """Score comparable pairs for both quantities and write the comparison."""
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    candidates: list = []
    for n in (4, 5):
        cands, _ = sweep_candidates(n, ROOT_SEED + n, 40)
        for c in cands:
            if 2 <= len(c.space) <= MAX_SPACE:
                candidates.append(c)
        if len(candidates) >= TARGET_INSTANCES:
            break
    candidates = candidates[:TARGET_INSTANCES]

    n_pairs = 0
    n_viol_mean_naive = 0
    n_viol_mean_mc = 0
    n_viol_b = 0
    per_instance: list[dict] = []
    margins_mc: list[float] = []

    for idx, cand in enumerate(candidates):
        pairs = enumerate_comparable_pairs(cand.space)[:MAX_PAIRS_PER_INSTANCE]
        if not pairs:
            continue
        sem = random_sem(cand.dag, np.random.default_rng(ROOT_SEED + 5000 + idx))
        ctx = make_context(sem, cand.cpdag, cand.treatment, cand.outcome, cand.z)

        b_cache: dict[str, float] = {}

        def b_of(g: MPDAG, _cache: dict[str, float] = b_cache, _ctx: BiasContext = ctx) -> float:
            # The cache and context are bound as defaults rather than captured:
            # both are rebuilt per instance, and a late-binding closure here
            # would silently score one instance's pairs against another's
            # covariance.
            key = g.edge_string()
            if key not in _cache:
                _cache[key] = bias_at(_ctx, g).worst
            return _cache[key]

        inst_naive = inst_mc = inst_b = 0
        for g, h in pairs:
            n_pairs += 1

            # B: closed form, no sampling, so a plain tolerance is the right test.
            if b_of(g) > b_of(h) + 1e-9:
                inst_b += 1

            # mean_abs_bias: sampled, so the deficit is judged against its own
            # Monte-Carlo error, in two independent repetitions.
            deficits, thresholds, naive_hits = [], [], 0
            for rep in range(N_REPETITIONS):
                sg = (ROOT_SEED, idx, rep, 1, hash(g.edge_string()) % 10**6)
                sh = (ROOT_SEED, idx, rep, 2, hash(h.edge_string()) % 10**6)
                mg = _mean_abs_bias(cand.z, g, cand.treatment, cand.outcome, sg, N_DRAWS)
                mh = _mean_abs_bias(cand.z, h, cand.treatment, cand.outcome, sh, N_DRAWS)
                deficit = mg.mean_abs_bias - mh.mean_abs_bias
                if deficit > 1e-6:
                    naive_hits += 1
                se_g = mc_standard_error(mg, g, cand.z, cand.treatment, cand.outcome, sg, N_DRAWS)
                se_h = mc_standard_error(mh, h, cand.z, cand.treatment, cand.outcome, sh, N_DRAWS)
                combined = float(np.hypot(se_g, se_h))
                deficits.append(deficit)
                thresholds.append(Z_MULTIPLIER * combined)

            if naive_hits == N_REPETITIONS:
                inst_naive += 1
            if all(d > t for d, t in zip(deficits, thresholds, strict=True)):
                inst_mc += 1
                margins_mc.append(
                    float(
                        min(d / t * Z_MULTIPLIER for d, t in zip(deficits, thresholds, strict=True))
                    )
                )

        n_viol_mean_naive += inst_naive
        n_viol_mean_mc += inst_mc
        n_viol_b += inst_b
        per_instance.append(
            {
                **cand.tag(),
                "n_pairs": len(pairs),
                "viol_mean_naive": inst_naive,
                "viol_mean_mc": inst_mc,
                "viol_B": inst_b,
            }
        )

    payload = {
        "design": {
            "rule": (
                f"a pair counts as a mean_abs_bias violation only if the deficit exceeds "
                f"{Z_MULTIPLIER} combined Monte-Carlo standard errors in ALL "
                f"{N_REPETITIONS} independent repetitions, at {N_DRAWS} draws per state. "
                "This is conservative, so the rate is a LOWER BOUND on the population rate."
            ),
            "z_multiplier": Z_MULTIPLIER,
            "n_draws": N_DRAWS,
            "n_repetitions": N_REPETITIONS,
            "B_tolerance": 1e-9,
            "root_seed": ROOT_SEED,
            "n_instances": len(per_instance),
            "max_space_size": MAX_SPACE,
            "max_pairs_per_instance": MAX_PAIRS_PER_INSTANCE,
        },
        "results": {
            "n_pairs_total": n_pairs,
            "mean_abs_bias_violations_naive_1e-6": n_viol_mean_naive,
            "mean_abs_bias_violations_mc_aware": n_viol_mean_mc,
            "B_violations": n_viol_b,
            "rate_mean_naive": round(n_viol_mean_naive / max(1, n_pairs), 4),
            "rate_mean_mc_aware_LOWER_BOUND": round(n_viol_mean_mc / max(1, n_pairs), 4),
            "rate_B": round(n_viol_b / max(1, n_pairs), 6),
            "mc_margin_in_SEs": {
                "min": round(min(margins_mc), 2) if margins_mc else None,
                "median": round(float(np.median(margins_mc)), 2) if margins_mc else None,
                "max": round(max(margins_mc), 2) if margins_mc else None,
            },
        },
        "per_instance": per_instance,
        "provenance": {
            "git_commit": _git_commit(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "seconds": round(time.perf_counter() - t0, 2),
        },
    }
    (OUT / "c2_rate_sweep_mc_aware.json").write_text(json.dumps(payload, indent=2))
    print(json.dumps(payload["results"], indent=2))
    print(f"instances={len(per_instance)} seconds={payload['provenance']['seconds']}")


if __name__ == "__main__":
    main()
