"""Exhaustive shell study: mean vs max absolute bias, per retraction shell.

``r_eps`` thresholds ``beta_up(d)``, the running maximum of ``B`` over the
retraction shell at depth ``d``. This script asks what the alternative design --
thresholding the **mean** of ``B`` over the shell instead of its maximum -- would
look like, and it answers by enumerating **every** retraction state exhaustively
rather than by searching for a radius.

For each instance it walks all ``2^|K_G0|`` subsets of the analyst's asserted
orientations, closes each one under Meek, and evaluates ``B`` at every distinct
resulting state. Per shell ``d`` it then records the full distribution of ``B``:
mean, max, median, quantiles, the contaminated fraction, and the state count.
Nothing is thresholded and no radius is searched for, so the output supports any
epsilon after the fact -- including epsilons below 0.01 and above 20, which is
what the committed grid cannot answer.

Two aggregations are reported per shell and they are not the same number:

* ``mean_states``  -- over the **distinct closures** in the shell. Meek's rules
  make many subsets collapse to the same state, and this counts each state once.
* ``mean_subsets`` -- over the ``C(|K_G0|, d)`` subsets, so a state that many
  subsets collapse to is weighted by how many. This is the quantity an analyst
  who draws a uniform random set of ``d`` mistaken claims would face.

Both are exact: no Monte-Carlo, no sampling, no early stop. The per-state ``B``
is the existing exact worst case over that state's ambiguity set, so this varies
only the shell-level aggregation -- the design decision actually under review.

Usage::

    PYTHONPATH=src .venv/bin/python experiments/mean_vs_max_shell.py [options]

Writes ``<out>/shells.jsonl`` (one record per instance/shell) and
``<out>/instances.json`` (per-instance metadata).
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import statistics
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

from bkrobust.demo.evaluate import random_sem  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.epsilon.bias import bias_at, knowledge_of, make_context  # noqa: E402
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag  # noqa: E402

import final_table as ft  # noqa: E402

#: Instances to study, as ``network:x:y``. Chosen to span |K_G0| from 3 to 14,
#: r_val from 1 to 11, and the whole observed range of B(Chat) -- including the
#: mass point at exactly 1.0 and the single instance above 20.
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


def study_instance(
    net: str, x: str, y: str, cpdag, dag, k, seed: int, max_subsets: int
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Enumerate every retraction state and summarise B per shell.

    Args:
        net: Network name.
        x: Treatment.
        y: Outcome.
        cpdag: The estimated CPDAG.
        dag: The ground-truth DAG, for the seeded SEM.
        k: The elicited knowledge.
        seed: Base seed; the per-instance seed is derived as in final_table.
        max_subsets: Refuse instances needing more than this many subsets.

    Returns:
        ``(meta, shell_records)``.
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
        return {"network": net, "x": x, "y": y, "status": "too_large", "n_knowledge": n}, []

    # B depends only on the closed state, and many subsets share a closure, so
    # evaluate once per distinct edge string and reuse.
    cache: dict[str, float] = {}
    shells: dict[int, dict[str, Any]] = {
        d: {"states": {}, "n_subsets": 0} for d in range(n + 1)
    }
    for d in range(n + 1):
        for drop in itertools.combinations(range(n), d):
            keep = [k0[i] for i in range(n) if i not in set(drop)]
            h = apply_orientations(cpdag, keep)
            if h is None:
                continue
            es = h.edge_string()
            if es not in cache:
                cache[es] = bias_at(ctx, h, method="semilocal").worst / scale
            sh = shells[d]
            sh["n_subsets"] += 1
            sh["states"].setdefault(es, {"b": cache[es], "mult": 0})["mult"] += 1

    records: list[dict[str, Any]] = []
    for d in range(n + 1):
        sh = shells[d]
        vals = [s["b"] for s in sh["states"].values()]
        if not vals:
            continue
        weighted = [s["b"] for s in sh["states"].values() for _ in range(s["mult"])]
        records.append(
            {
                "network": net,
                "x": x,
                "y": y,
                "d": d,
                "n_subsets": sh["n_subsets"],
                "n_states": len(vals),
                "max": max(vals),
                "mean_states": statistics.fmean(vals),
                "mean_subsets": statistics.fmean(weighted),
                "median_states": statistics.median(vals),
                "min": min(vals),
                "std_states": statistics.pstdev(vals) if len(vals) > 1 else 0.0,
                "frac_nonzero": sum(1 for v in vals if v > 0) / len(vals),
                "frac_nonzero_subsets": sum(1 for v in weighted if v > 0) / len(weighted),
            }
        )

    meta = {
        "network": net,
        "x": x,
        "y": y,
        "status": "ok",
        "n_knowledge": n,
        "n_subsets_total": total,
        "n_distinct_states": len(cache),
        "theta_z": ctx.theta_z,
        "beta_top": max(cache.values()) if cache else None,
        "seconds": round(time.perf_counter() - t0, 3),
    }
    return meta, records


def main() -> None:
    """Run the exhaustive shell study and write its two output files."""
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
    p.add_argument("--out", default="results/mean_vs_max")
    args = p.parse_args()

    wanted = [tuple(s.split(":", 2)) for s in args.instances.split(",") if s]
    nets_wanted = {w[0] for w in wanted}
    knowledge = ft.load_knowledge(args.condition)
    parsed = ft.load_networks(names=nets_wanted, skip=set(), max_nodes=0)

    out_dir = ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    metas, all_records = [], []
    for net, x, y in wanted:
        if net not in parsed or net not in knowledge:
            print(f"[mean_vs_max] {net}: unavailable", file=sys.stderr)
            continue
        meta, recs = study_instance(
            net, x, y, parsed[net]["cpdag"], parsed[net]["dag"],
            knowledge[net], args.seed, args.max_subsets,
        )
        metas.append(meta)
        all_records.extend(recs)
        print(
            f"[mean_vs_max] {net} {x}->{y} status={meta['status']} "
            f"|K_G0|={meta.get('n_knowledge')} states={meta.get('n_distinct_states')} "
            f"seconds={meta.get('seconds')}",
            flush=True,
        )

    with (out_dir / "shells.jsonl").open("w") as f:
        for r in all_records:
            f.write(json.dumps(r) + "\n")
    (out_dir / "instances.json").write_text(
        json.dumps({"args": vars(args), "instances": metas}, indent=2) + "\n"
    )
    print(f"[mean_vs_max] wrote {out_dir}/shells.jsonl and instances.json")


if __name__ == "__main__":
    main()
