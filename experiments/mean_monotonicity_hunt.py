r"""A search for a non-monotone ``mu(d)``, which would sink the mean radius.

``r_mu(eps) = min\{d : mu(d) > eps\}`` is only a radius if ``mu`` is
non-decreasing in ``d``. Crossing a threshold on a curve that later falls back
does not mean the analyst stays across it, and the guarantee a certificate
exists to provide would not follow. Nothing in this codebase proves ``mu`` is
monotone.

That is not an idle worry. The mean of ``B`` over the **full space** sphere at
distance ``d`` *is* non-monotone, on 34 of 79 instances measured
(``docs/RESULT_MEAN_SOUNDNESS_GATE.md``): far shells contain heavily oriented
states that pin the effect down, so they carry low bias and drag the average
back. ``mu`` restricted to the retraction up-set appears to escape this, and
this script is the attempt to falsify that.

It sweeps node counts, edge densities, every ordered treatment/outcome pair and
several SEM draws per graph, enumerates the **whole** ``mu`` profile
exhaustively for each accepted instance, and records every profile that is not
non-decreasing. Finding one counterexample would settle the question against the
statistic; finding none leaves monotonicity as an observation with a measured
sample size behind it, which is what the report is entitled to claim.

Usage::

    PYTHONPATH=src .venv/bin/python experiments/mean_monotonicity_hunt.py [options]

Writes ``<out>/hunt.json``: the profile count, the configuration swept, and any
violations in full so that each can be replayed.
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from pathlib import Path
from typing import Any, Iterator

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag, random_sem  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.epsilon.bias import knowledge_of, make_context  # noqa: E402
from bkrobust.epsilon.meanprofile import is_monotone, shell_stat  # noqa: E402


def candidates(
    n_nodes: int, seed: int, p: float, sem_draws: int, max_k: int
) -> Iterator[tuple[Any, ...]]:
    """Yield every screened ``(cpdag, g0, ctx, scale, ...)`` from one random DAG.

    Args:
        n_nodes: Node count.
        seed: Graph seed.
        p: Edge probability in the upper triangle, so acyclicity is automatic.
        sem_draws: Independent SEM draws per query.
        max_k: Skip graphs whose ``|K_G0|`` exceeds this, to keep the
            exhaustive profile affordable.

    Yields:
        One tuple per (query, SEM draw) that survives screening.
    """
    rng = np.random.default_rng(seed)
    nodes = [f"V{i}" for i in range(n_nodes)]
    edges = [
        (nodes[i], nodes[j])
        for i in range(n_nodes)
        for j in range(i + 1, n_nodes)
        if rng.random() < p
    ]
    if not edges:
        return
    dag = MPDAG(nodes, directed=edges)
    try:
        cpdag = dag_to_cpdag(dag)
    except Exception:  # noqa: BLE001 - a rejected draw is not an error here
        return
    g0 = apply_orientations(cpdag, sorted(set(dag.directed_edges) - set(cpdag.directed_edges)))
    if g0 is None:
        return
    n_k = len(knowledge_of(cpdag, g0))
    if not 2 <= n_k <= max_k:
        return
    for x, y in itertools.permutations(nodes, 2):
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None:
            continue
        for draw in range(sem_draws):
            ctx = make_context(
                random_sem(dag, np.random.default_rng(seed * 97 + draw)),
                cpdag, x, y, frozenset(z),
            )
            if abs(ctx.theta_z) < 1e-12:
                continue
            yield cpdag, g0, ctx, abs(ctx.theta_z), n_k, x, y, draw


def main() -> None:
    """Sweep the configuration grid and record every non-monotone profile."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--nodes", default="4,5,6,7")
    p.add_argument("--densities", default="0.3,0.45,0.6,0.75")
    p.add_argument("--seeds", type=int, default=400)
    p.add_argument("--sem-draws", type=int, default=3)
    p.add_argument("--max-knowledge", type=int, default=12)
    p.add_argument("--time-limit", type=float, default=1800.0)
    p.add_argument("--out", default="results/mean_monotonicity")
    args = p.parse_args()

    t0 = time.perf_counter()
    checked = 0
    violations: list[dict[str, Any]] = []
    for n_nodes in (int(v) for v in args.nodes.split(",")):
        for dens in (float(v) for v in args.densities.split(",")):
            for s in range(args.seeds):
                seed = s + n_nodes * 7919 + int(dens * 1000) * 13
                for cpdag, g0, ctx, scale, n_k, x, y, draw in candidates(
                    n_nodes, seed, dens, args.sem_draws, args.max_knowledge
                ):
                    cache: dict[str, float] = {}
                    prof = [
                        shell_stat(ctx, cpdag, g0, d, scale=scale, cache=cache)
                        for d in range(n_k + 1)
                    ]
                    checked += 1
                    if not is_monotone(prof):
                        violations.append({
                            "n_nodes": n_nodes, "density": dens, "seed": seed,
                            "sem_draw": draw, "x": x, "y": y, "n_knowledge": n_k,
                            "mu": [s_.mean for s_ in prof],
                            "cpdag": cpdag.edge_string(), "g0": g0.edge_string(),
                        })
                        print(f"[hunt] NON-MONOTONE {n_nodes} {dens} {seed} {x}->{y}", flush=True)
                    if checked % 5000 == 0:
                        print(f"[hunt] {checked} profiles, {len(violations)} violations, "
                              f"{time.perf_counter() - t0:.0f}s", flush=True)
                if time.perf_counter() - t0 > args.time_limit:
                    break
            if time.perf_counter() - t0 > args.time_limit:
                break
        if time.perf_counter() - t0 > args.time_limit:
            break

    out_dir = ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "hunt.json").write_text(json.dumps({
        "args": vars(args), "profiles_checked": checked,
        "n_violations": len(violations), "violations": violations,
        "seconds": round(time.perf_counter() - t0, 1),
    }, indent=2) + "\n")
    print(f"[hunt] DONE: {checked} profiles, {len(violations)} non-monotone, "
          f"{time.perf_counter() - t0:.0f}s -> {out_dir}/hunt.json")


if __name__ == "__main__":
    main()
