"""Does the epsilon radius survive having to estimate the covariance?

``r_val`` is graphical: it reads no numbers, so it has no sampling distribution.
``r_eps`` reads a covariance, so it does, and a plug-in ``r_eps`` reported as
though it were exact would be a confidence statement with the confidence removed.
This script measures the damage and checks that the conservative construction of
:mod:`bkrobust.epsilon.finite_sample` repairs it in the right direction.

Three quantities per instance and sample size:

* **population** ``r_eps``, from the exact covariance of the true SEM. This is
  the estimand; it is not available to an analyst and is used only to score.
* **plug-in** ``r_eps``, from the sample covariance. A point estimate.
* **conservative** ``r_eps``, the 5% quantile of the bootstrap distribution of
  the plug-in. The reportable number.

The claim being tested is one-sided and is the only one that matters for a
certificate: a radius that comes back **too large** overstates robustness, and
is the error to avoid; a radius that comes back too small understates it, and is
safe. So the headline is not "how often is the plug-in right" but **how often is
it too large**, and whether the conservative version drives that rate below the
nominal 5%.

Determinism: one generator per instance, one per sample draw, one for the
bootstrap, all seeded from a fixed root; the global NumPy RNG is never touched.
"""

from __future__ import annotations

import json
import platform
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from bkrobust.core.conventions import UNREACHED
from bkrobust.core.oracle import is_valid
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag, random_sem
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations
from bkrobust.epsilon.bias import make_context
from bkrobust.epsilon.finite_sample import bootstrap_radius, sample_data
from bkrobust.epsilon.profile import epsilon_grid
from bkrobust.hybrid import breakdown_radius
from bkrobust.synth.generators import GENERATORS
from bkrobust.synth.knowledge import draw_k_true, flip
from bkrobust.synth.runner import _pick_treatment_outcome, gate

ROOT_SEED = 20260916
SAMPLE_SIZES: tuple[int, ...] = (100, 500, 2_000, 10_000)
EPS_RELATIVE: tuple[float, ...] = (0.05, 0.25)
N_BOOT = 120
ALPHA = 0.05
TARGET_INSTANCES = 120
MAX_ATTEMPTS = 8_000
MAX_UNDIRECTED = 6

OUT = Path("results/epsilon/finite_sample")


def _git_commit() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True
        ).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def _instance(seed: int) -> tuple[tuple | None, str]:
    """One screened instance, or ``None`` with the rejection reason."""
    rng = np.random.default_rng(seed)
    n = int(rng.choice([5, 6, 7]))
    dag = GENERATORS["erdos_renyi"](n, rng, edge_prob=float(rng.choice([0.3, 0.4, 0.5])))
    cpdag = dag_to_cpdag(dag)
    if len(cpdag.undirected_edges) > MAX_UNDIRECTED:
        return None, "cpdag_too_large"
    x, y = _pick_treatment_outcome("erdos_renyi", dag, rng)
    ok, reason = gate(dag, cpdag, x, y)
    if not ok:
        return None, reason
    k = flip(draw_k_true(dag, cpdag, rng, 0.8), rng, float(rng.choice([0.0, 0.2])))
    g0 = apply_orientations(cpdag, k)
    if g0 is None:
        return None, "k_inconsistent"
    z = optimal_adjustment_set_mpdag(g0, x, y)
    if z is None:
        return None, "z_not_identified"
    zf = frozenset(z)
    if not is_valid(zf, g0, x, y):
        return None, "z_invalid_at_g0"
    sem = random_sem(dag, rng)
    return (dag, cpdag, g0, x, y, zf, sem), "ok"


def main() -> None:
    """Run the coverage study and write the results."""
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    rejections: Counter = Counter()
    rows: list[dict] = []
    accepted = 0
    for attempt in range(MAX_ATTEMPTS):
        if accepted >= TARGET_INSTANCES:
            break
        built, reason = _instance(ROOT_SEED + attempt)
        if built is None:
            rejections[reason] += 1
            continue
        dag, cpdag, g0, x, y, zf, sem = built

        hres = breakdown_radius(cpdag, None, x, y, zf, g0=g0)
        if hres.radius == UNREACHED or not hres.exact:
            rejections["r_val_unreached_or_inexact"] += 1
            continue
        r_val = hres.radius

        ctx = make_context(sem, cpdag, x, y, zf)
        if abs(ctx.theta_z) < 1e-8:
            rejections["zero_estimate"] += 1
            continue

        # The epsilon grid is fixed in EFFECT units from the population estimate,
        # so every sample size and every bootstrap replicate is judged against
        # the same bar. Re-deriving a relative threshold inside each replicate
        # would move the bar with the noise and confound the two effects.
        eps_effect = [e * abs(ctx.theta_z) for e in EPS_RELATIVE]
        pop, _ = epsilon_grid(ctx, g0, eps_effect, r_val=r_val)

        row: dict = {
            "attempt": attempt,
            "n_nodes": len(dag.nodes),
            "n_k0": len(set(g0.directed_edges) - set(cpdag.directed_edges)),
            "r_val": r_val,
            "theta_z": ctx.theta_z,
            "population": {
                f"{e:g}": pop[ee] for e, ee in zip(EPS_RELATIVE, eps_effect, strict=True)
            },
            "by_sample_size": {},
        }
        for n_obs in SAMPLE_SIZES:
            data = sample_data(sem, n_obs, np.random.default_rng(ROOT_SEED + 10_000 + attempt))
            dist = bootstrap_radius(
                cpdag,
                g0,
                x,
                y,
                zf,
                data,
                sem.dag.nodes,
                eps_effect,
                r_val=r_val,
                n_boot=N_BOOT,
                rng=np.random.default_rng(ROOT_SEED + 20_000 + attempt),
            )
            row["by_sample_size"][str(n_obs)] = {
                f"{e:g}": {
                    "plugin": dist[ee].point,
                    "conservative": dist[ee].conservative(ALPHA),
                }
                for e, ee in zip(EPS_RELATIVE, eps_effect, strict=True)
            }
        rows.append(row)
        accepted += 1

    # Scoring. UNREACHED is ordered ABOVE every finite radius, because it means
    # "no crossing anywhere", which is the largest answer and hence the most
    # optimistic one -- treating it as a small number would invert the direction
    # of the very error this study is about.
    def _rank(r: int) -> float:
        return float("inf") if r == UNREACHED else float(r)

    summary: dict = {}
    for n_obs in SAMPLE_SIZES:
        for e in EPS_RELATIVE:
            key = f"n={n_obs}, eps={e:g}"
            agg = defaultdict(int)
            for row in rows:
                truth = _rank(row["population"][f"{e:g}"])
                cell = row["by_sample_size"][str(n_obs)][f"{e:g}"]
                for name in ("plugin", "conservative"):
                    got = _rank(cell[name])
                    agg[f"{name}_n"] += 1
                    if got == truth:
                        agg[f"{name}_exact"] += 1
                    elif got > truth:
                        agg[f"{name}_too_large"] += 1  # the dangerous direction
                    else:
                        agg[f"{name}_too_small"] += 1
            summary[key] = {
                name: {
                    "n": agg[f"{name}_n"],
                    "exact": agg[f"{name}_exact"],
                    "too_large_OVERSTATES_ROBUSTNESS": agg[f"{name}_too_large"],
                    "too_small_safe": agg[f"{name}_too_small"],
                    "rate_exact": round(agg[f"{name}_exact"] / max(1, agg[f"{name}_n"]), 4),
                    "rate_too_large": round(agg[f"{name}_too_large"] / max(1, agg[f"{name}_n"]), 4),
                }
                for name in ("plugin", "conservative")
            }

    payload = {
        "design": {
            "root_seed": ROOT_SEED,
            "sample_sizes": list(SAMPLE_SIZES),
            "eps_relative": list(EPS_RELATIVE),
            "n_boot": N_BOOT,
            "alpha": ALPHA,
            "target_instances": TARGET_INSTANCES,
            "n_accepted": len(rows),
            "n_attempts": min(MAX_ATTEMPTS, attempt + 1),
            "rejection_counts": dict(sorted(rejections.items())),
            "note": (
                "r_val is graphical and has no sampling distribution; only r_eps is "
                "estimated here. UNREACHED is ordered above every finite radius when "
                "scoring, since it is the most optimistic possible answer."
            ),
        },
        "summary": summary,
        "rows": rows,
        "provenance": {
            "git_commit": _git_commit(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "seconds": round(time.perf_counter() - t0, 2),
        },
    }
    (OUT / "finite_sample.json").write_text(json.dumps(payload, indent=2))

    lines = [
        "# Does `r_eps` survive an estimated covariance?",
        "",
        f"{len(rows)} instances, {N_BOOT} bootstrap replicates, alpha = {ALPHA}.",
        "",
        "`r_val` is not in this table because it has no sampling distribution: it reads",
        "the CPDAG and the asserted knowledge, and no numbers. Only `r_eps` is estimated.",
        "",
        "The column that matters is **too large**. A radius that is too large claims more",
        "certified shells than the truth supports, which overstates robustness; a radius",
        "that is too small merely wastes some. The conservative construction exists to",
        "drive the first rate below the nominal alpha, and is only useful if it does.",
        "",
        "| cell | estimator | exact | too small (safe) | too large (overstates) | rate too large |",
        "|---|---|---|---|---|---|",
    ]
    for key, cell in summary.items():
        for name in ("plugin", "conservative"):
            c = cell[name]
            lines.append(
                f"| {key} | {name} | {c['exact']}/{c['n']} | {c['too_small_safe']} "
                f"| {c['too_large_OVERSTATES_ROBUSTNESS']} | {c['rate_too_large']:.1%} |"
            )
    (OUT / "SUMMARY.md").write_text("\n".join(lines))
    print(f"wrote {OUT / 'finite_sample.json'} and {OUT / 'SUMMARY.md'}  ({len(rows)} instances)")
    for key, cell in summary.items():
        print(
            f"  {key}: plugin too-large {cell['plugin']['rate_too_large']:.1%}, "
            f"conservative too-large {cell['conservative']['rate_too_large']:.1%}"
        )


if __name__ == "__main__":
    main()
