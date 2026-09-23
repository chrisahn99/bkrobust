r"""The mean bias profile on every real network, reported alongside ``r_val``.

``r_eps`` thresholds the worst case ``beta_up(d)``, which on this corpus is flat
above ``r_val``: it jumps to ``B(Chat)`` at the validity boundary and stays, so
``r_eps`` reproduces ``r_val`` and the epsilon axis says nothing
(``docs/RESULT_EPSILON_ABOVE_ONE.md``). This script computes the companion
quantity that does vary with depth,

.. math::

    \mu(d) = \mathbb{E}\big[B(\mathrm{Meek}(\hat C, K_{G_0}\setminus S))\big],
    \quad S\sim\mathrm{Uniform}\{S\subseteq K_{G_0} : |S| = d\}

-- *if exactly ``d`` of the orientations the analyst asserted are wrong, this is
the expected bias* -- and reports the radius ``r_mu(eps) = min\{d : \mu(d) >
eps\}`` next to ``r_val`` for every query.

Why this is affordable
----------------------

``mu`` is never reported alone: it accompanies ``r_val``, and Theorem A makes
every shell strictly below ``r_val`` identically zero (a maximum of zero over
non-negative values forces every state to zero). Those shells cannot cross a
non-negative threshold, so the walk starts **at** ``r_val`` and skips everything
below it. On this corpus that turns the enumeration from intractable into
49,427 subsets in total, and every instance is exhaustive -- no sampling, no
estimate, no confidence interval needed.

Instances whose ``r_val`` is itself ``UNREACHED`` carry zero bias everywhere by
the same theorem, so their profile is identically zero and their ``r_mu`` is
``UNREACHED`` at every tolerance. They are reported that way and not enumerated.

What the number means, and does not
-----------------------------------

``mu`` is an **average-case** quantity under a **uniform prior** over which
``d`` of the analyst's claims are wrong. It does **not** bound the worst case --
``beta_up`` does, and remains what any guarantee is stated against. It is also
not the mean over the full perturbation space, which is a different quantity
that disagrees with this one in both directions
(``docs/RESULT_MEAN_SOUNDNESS_GATE.md``). Monotonicity of ``mu`` in ``d`` is
observed rather than proved, so it is checked per instance and reported.

``r_val`` is purely graph-theoretic. ``mu`` and every radius derived from it are
**SEM-conditional**: the corpus ships graph structure only, so the bias is
evaluated against a seeded linear-Gaussian SEM attached to the ground-truth DAG,
seeded per instance from ``--seed`` exactly as ``final_table.py`` does.

Usage::

    PYTHONPATH=src .venv/bin/python experiments/mean_table.py [options]

Writes ``<out>/instances.jsonl`` and ``<out>/summary.json``, and prints a
per-query table.
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
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag, random_sem  # noqa: E402
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

#: Reporting grid, in relative units -- fractions of the reported effect's own
#: magnitude, so the columns read as 2%, 10%, 25%, 40%, 60%, 90%.
#:
#: Three constraints picked it, in this order. (1) It must not sit on 1.0: half
#: the informative instances have ``B(Chat)`` equal to 1.0 to within 4e-15, and
#: a strict ``>`` test at a grid point placed there splits tied instances by
#: floating-point dust (``docs/RESULT_EPSILON_ABOVE_ONE.md`` section 3).
#: (2) It must stay below the ceiling. ``mu(|K_G0|) = B(Chat)`` exactly -- the
#: top shell is the single state ``Chat`` -- so the mean and the max share a
#: ceiling and every tolerance above an instance's own ``B(Chat)`` is
#: ``UNREACHED`` by construction, wasting the column on the 39 of 53 instances
#: whose ceiling is at most 1.0. A 1.5 column was measured and rejected: it buys
#: 0.08 distinct radii per instance and triples the ``UNREACHED`` rate from 5%
#: to 17%. ``beta_top`` is reported as its own column instead, which says where
#: each instance actually saturates rather than probing one more threshold.
#: (3) Round, practitioner-legible values only, to avoid tuning the grid to this
#: corpus. Among all round-number grids this one maximises the distinct radii
#: per instance, and it wins or ties on every corpus split tested (by ceiling,
#: by ``|K_G0|``, by ``r_val``).
DEFAULT_EPSILONS = "0.02,0.10,0.25,0.40,0.60,0.90"


def _panel(path: Path) -> list[dict[str, Any]]:
    """Load the committed ``r_val`` panel.

    ``r_val`` is reused rather than recomputed: the accelerated hybrid search is
    the expensive part of the sweep, and reusing it at the same seed keeps every
    radius here directly comparable to the committed table.

    Args:
        path: Path to an ``instances.jsonl``.

    Returns:
        The records.
    """
    return [json.loads(line) for line in path.open()]


def profile_for(
    net: str,
    x: str,
    y: str,
    r_val: int,
    cpdag: Any,
    dag: Any,
    k: list[tuple[str, str]],
    seed: int,
    *,
    shell_cap: int,
    rng_seed: int,
) -> tuple[list[ShellStat], dict[str, Any]]:
    """Enumerate ``mu(d)`` from ``r_val`` upward for one instance.

    Args:
        net: Network name.
        x: Treatment.
        y: Outcome.
        r_val: The committed validity radius (finite, and at least 1).
        cpdag: The CPDAG.
        dag: The ground-truth DAG, for the seeded SEM.
        k: The elicited knowledge.
        seed: Base seed, for the per-instance SEM.
        shell_cap: Sample a shell rather than enumerate it above this many
            subsets.
        rng_seed: Seed for any sampling.

    Returns:
        ``(profile, meta)``. ``profile`` is empty when the instance is
        degenerate or its estimate is zero.
    """
    g0 = apply_orientations(cpdag, list(k))
    if g0 is None:
        return [], {"status": "inconsistent"}
    z = optimal_adjustment_set_mpdag(g0, x, y)
    if z is None:
        return [], {"status": "degenerate"}
    ctx = make_context(
        random_sem(dag, np.random.default_rng(ft.derive_seed(seed, net, x, y))),
        cpdag, x, y, frozenset(z),
    )
    scale = abs(ctx.theta_z)
    if scale < 1e-12:
        return [], {"status": "zero_estimate"}

    n = len(knowledge_of(cpdag, g0))
    start = max(0, min(r_val, n))
    rng = np.random.default_rng(rng_seed)
    cache: dict[str, float] = {}
    prof: list[ShellStat] = []
    skipped = sum(math.comb(n, d) for d in range(0, start))
    for d in range(start, n + 1):
        total = math.comb(n, d)
        prof.append(
            shell_stat(
                ctx, cpdag, g0, d, scale=scale, cache=cache,
                sample=None if total <= shell_cap else shell_cap, rng=rng,
            )
        )
    return prof, {
        "status": "ok",
        "n_knowledge": n,
        "theta_z": ctx.theta_z,
        "subsets_skipped_below_r_val": skipped,
        "subsets_evaluated": sum(s.n_evaluated for s in prof),
        "all_exhaustive": all(s.exhaustive for s in prof),
    }


def main() -> None:
    """Compute the mean profile for the whole corpus and write the table."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--condition", default="D_LLM")
    p.add_argument("--panel", default="results/final_table_eps_gt1/instances.jsonl")
    p.add_argument("--epsilons", default=DEFAULT_EPSILONS)
    p.add_argument("--seed", type=int, default=20260917)
    p.add_argument("--rng-seed", type=int, default=424242)
    p.add_argument("--shell-cap", type=int, default=20000,
                   help="sample a shell above this many subsets")
    p.add_argument("--out", default="results/mean_table")
    args = p.parse_args()

    epsilons = [float(v) for v in args.epsilons.split(",")]
    panel = _panel(ROOT / args.panel)
    knowledge = ft.load_knowledge(args.condition)
    wanted = sorted({r["network"] for r in panel if r["status"] == "ok"})
    parsed = ft.load_networks(names=set(wanted), skip=set(), max_nodes=0)

    out_dir = ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, Any]] = []
    t0 = time.perf_counter()

    for row in panel:
        if row["status"] != "ok" or row["r_val"] in (0, None):
            continue  # not ok, or degenerate (r_val == 0): ill-posed, excluded
        net, x, y = row["network"], row["x"], row["y"]
        rec: dict[str, Any] = {
            "network": net, "x": x, "y": y, "condition": args.condition,
            "r_val": row["r_val"], "n_knowledge": row["n_knowledge"],
            "beta_top": row.get("beta_top"),
        }
        if row["r_val"] == UNREACHED:
            # Theorem A with an unreachable validity radius: the bias is zero on
            # the whole space, so mu is identically zero and no tolerance is
            # ever crossed. Recorded, not enumerated.
            rec.update({
                "status": "ok", "profile": [], "zero_everywhere": True,
                "monotone": True,
                "r_mu": {str(e): UNREACHED for e in epsilons},
                "subsets_evaluated": 0,
            })
            records.append(rec)
            continue
        if net not in parsed or net not in knowledge:
            rec.update({"status": "unavailable"})
            records.append(rec)
            continue
        prof, meta = profile_for(
            net, x, y, row["r_val"], parsed[net]["cpdag"], parsed[net]["dag"],
            knowledge[net], args.seed, shell_cap=args.shell_cap, rng_seed=args.rng_seed,
        )
        rec.update(meta)
        if prof:
            rec["zero_everywhere"] = all(s.mean == 0.0 for s in prof)
            rec["monotone"] = is_monotone(prof)
            rec["r_mu"] = {str(e): r_mean(prof, e) for e in epsilons}
            rec["profile"] = [
                {
                    "d": s.d, "mean": s.mean, "beta_up": s.beta_up,
                    "frac_nonzero": s.frac_nonzero, "min_nonzero": s.min_nonzero,
                    "n_subsets_total": s.n_subsets_total, "n_evaluated": s.n_evaluated,
                    "n_states": s.n_states, "se": s.se, "exhaustive": s.exhaustive,
                    "bound_lo": mean_bounds(s)[0], "bound_hi": mean_bounds(s)[1],
                }
                for s in prof
            ]
        records.append(rec)
        print(
            f"[mean_table] {net} {x}->{y} r_val={row['r_val']} |K|={rec.get('n_knowledge')} "
            f"subsets={rec.get('subsets_evaluated')} monotone={rec.get('monotone')} "
            f"r_mu={rec.get('r_mu')}",
            flush=True,
        )

    with (out_dir / "instances.jsonl").open("w") as f:
        for r in records:
            f.write(json.dumps(r) + "\n")

    ok = [r for r in records if r["status"] == "ok"]
    summary = {
        "args": vars(args),
        "epsilons": epsilons,
        "n_records": len(records),
        "n_ok": len(ok),
        "n_zero_everywhere": sum(1 for r in ok if r.get("zero_everywhere")),
        "n_non_monotone": sum(1 for r in ok if r.get("monotone") is False),
        "subsets_evaluated": sum(r.get("subsets_evaluated") or 0 for r in ok),
        "subsets_skipped_below_r_val": sum(
            r.get("subsets_skipped_below_r_val") or 0 for r in ok
        ),
        "all_exhaustive": all(r.get("all_exhaustive", True) for r in ok),
        "seconds": round(time.perf_counter() - t0, 2),
        "assumptions": ft.ASSUMPTIONS if hasattr(ft, "ASSUMPTIONS") else "",
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\n[mean_table] {summary['n_ok']}/{summary['n_records']} ok, "
          f"{summary['subsets_evaluated']:,} subsets evaluated, "
          f"{summary['subsets_skipped_below_r_val']:,} skipped below r_val, "
          f"non-monotone {summary['n_non_monotone']}, {summary['seconds']}s")
    print(f"[mean_table] wrote {out_dir}")


if __name__ == "__main__":
    main()
