"""The soundness gate: does a mean-based radius survive leaving the up-set?

``beta_up`` is computed over the **retraction up-set** ``^G0`` rather than over
the whole space ``G_Chat``, and Theorem C(3)-(5) of ``docs/R_EPSILON_THEORY.md``
is what licenses that: the maximum over the retraction shell at ``d`` equals the
maximum over the ball at ``d``, so the first crossing of ``beta_up`` and of the
full-space profile ``beta`` coincide. The cheap search is exact.

**No analogue of that theorem is known for a mean.** An average over a subset of
the space is not the average over the space, so a mean-based radius computed the
way the code computes things today would be an average over the retraction
up-set specifically. This script measures the gap directly: it enumerates the
**whole** space, computes ``B`` at every element, and compares, per shell,

* ``max``  over the full sphere / full ball / the retraction shell, and
* ``mean`` over the full sphere / full ball / the retraction shell,

then compares the radii each rule produces over a wide epsilon grid. The max
comparison is the harness's own control: Theorem C(5) says those first crossings
must agree, so a disagreement there would indict the harness rather than the
mean.

Two things this script is careful about, both of which would silently corrupt a
mean:

* It builds the space with :func:`bkrobust.search.space_fixed.build_space_correct`,
  **not** ``core.spacelib.build_space``. The latter filters candidates through a
  chordality test that is right for CPDAGs and wrong for MPDAGs, and it drops up
  to 10% of the elements, the loss growing with density (see that module's
  docstring). A mean over a space missing a tenth of its states is wrong in an
  uncontrolled direction; a max is much less sensitive to it.
* It evaluates ``B`` with ``method="enumerate"``, the definitional evaluator,
  rather than the semi-local shortcut whose proof is stated for this setting but
  whose differential test is reported here alongside, per element.

Usage::

    PYTHONPATH=src .venv/bin/python experiments/mean_fullspace_check.py [options]

Writes ``<out>/shells.jsonl``, ``<out>/instances.json``.
"""

from __future__ import annotations

import argparse
import itertools
import json
import signal
import statistics
import sys
import zlib
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

from bkrobust.core.spacelib import distances_from  # noqa: E402
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag, random_sem  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.epsilon.bias import bias_at, knowledge_of, make_context  # noqa: E402
from bkrobust.search.space_fixed import build_space_correct  # noqa: E402

import final_table as ft  # noqa: E402

#: Real-network instances to attempt, as ``network:x:y``. Attempted in
#: increasing order of CPDAG undirected-edge count, which is the cost knob:
#: the corrected enumerator walks ``3^k`` assignments.
DEFAULT_INSTANCES = (
    "Acid_1996:x1:x10",
    "Thoemmes_2013:e0:s2",
    "asia:bronc:dysp",
    "asia:asia:either",
    "Didelez_2010:Age:HRT",
    "mediator:I:Y",
    "Sebastiani_2005:ANXA2.5:ANXA2.11",
    "magic-niab:G1217:YLD",
    "water:CKNI_12_15:CBODN_12_45",
    "Schipf_2010:A:TT",
    "barley:dg25:s2225",
    "magic-irri:G3212:FT",
)


class Timeout(Exception):
    """Raised when one instance exceeds its budget."""


def _alarm(_signum: int, _frame: Any) -> None:
    raise Timeout()


def analyse(
    tag: dict[str, Any], cpdag, g0, ctx, *, differential: bool
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Compare max and mean over the full space against the retraction shell.

    Args:
        tag: Identifying fields copied onto every record.
        cpdag: The CPDAG.
        g0: The analyst's state.
        ctx: The bias context.
        differential: Also evaluate B semi-locally and report disagreements.

    Returns:
        ``(meta, records)``, one record per shell depth.
    """
    t0 = time.perf_counter()
    scale = abs(ctx.theta_z)
    space = build_space_correct(cpdag)
    if g0 not in space.neighbours:
        return {**tag, "status": "g0_not_in_space"}, []
    dists = distances_from(space, g0)

    # B at every element of the whole space, by the definitional evaluator.
    b_full: dict[Any, float] = {}
    n_disagree = 0
    for g in space.elements:
        val = bias_at(ctx, g, method="enumerate").worst / scale
        b_full[g] = val
        if differential:
            alt = bias_at(ctx, g, method="semilocal").worst / scale
            if abs(alt - val) > 1e-9 * max(1.0, abs(val)):
                n_disagree += 1

    # The retraction shell, built exactly as profile.retraction_shell does.
    k0 = knowledge_of(cpdag, g0)
    n = len(k0)
    shells: dict[int, dict[str, Any]] = {d: {"states": {}, "n_subsets": 0} for d in range(n + 1)}
    outside = 0
    for d in range(n + 1):
        for drop in itertools.combinations(range(n), d):
            keep = [k0[i] for i in range(n) if i not in set(drop)]
            h = apply_orientations(cpdag, keep)
            if h is None:
                continue
            if h not in b_full:
                outside += 1
                continue
            sh = shells[d]
            sh["n_subsets"] += 1
            sh["states"].setdefault(h, {"b": b_full[h], "mult": 0, "dist": dists.get(h)})[
                "mult"
            ] += 1

    reachable = [g for g in space.elements if dists.get(g) is not None]
    max_d = max(dists[g] for g in reachable)
    records: list[dict[str, Any]] = []
    for d in range(max(max_d, n) + 1):
        sphere = [b_full[g] for g in reachable if dists[g] == d]
        ball = [b_full[g] for g in reachable if dists[g] <= d]
        sh = shells.get(d, {"states": {}, "n_subsets": 0})
        rvals = [s["b"] for s in sh["states"].values()]
        rw = [s["b"] for s in sh["states"].values() for _ in range(s["mult"])]
        rec: dict[str, Any] = {
            **tag,
            "d": d,
            "n_sphere": len(sphere),
            "n_ball": len(ball),
            "sphere_max": max(sphere) if sphere else None,
            "sphere_mean": statistics.fmean(sphere) if sphere else None,
            "ball_max": max(ball) if ball else None,
            "ball_mean": statistics.fmean(ball) if ball else None,
            "n_retr_subsets": sh["n_subsets"],
            "n_retr_states": len(rvals),
            "retr_max": max(rvals) if rvals else None,
            "retr_mean_states": statistics.fmean(rvals) if rvals else None,
            "retr_mean_subsets": statistics.fmean(rw) if rw else None,
            # Corollary R1: the retraction shell at d sits between the sphere at
            # d and the ball at d, so its members may be nearer than d. Recorded
            # rather than assumed.
            "retr_dist_max": max((s["dist"] for s in sh["states"].values() if s["dist"] is not None), default=None),
            "retr_dist_min": min((s["dist"] for s in sh["states"].values() if s["dist"] is not None), default=None),
        }
        records.append(rec)

    meta = {
        **tag,
        "status": "ok",
        "n_knowledge": n,
        "space_size": len(space),
        "n_reachable": len(reachable),
        "max_distance": max_d,
        "n_undirected_cpdag": len(cpdag.undirected_edges),
        "theta_z": ctx.theta_z,
        "beta_top_full": max(b_full.values()),
        "retr_states_outside_space": outside,
        "semilocal_disagreements": n_disagree if differential else None,
        "seconds": round(time.perf_counter() - t0, 2),
    }
    return meta, records


def real_instances(names: list[str], condition: str, seed: int, budget: int, differential: bool):
    """Yield ``(meta, records)`` for each requested real-network instance."""
    wanted = [tuple(s.split(":", 2)) for s in names if s]
    knowledge = ft.load_knowledge(condition)
    parsed = ft.load_networks(names={w[0] for w in wanted}, skip=set(), max_nodes=0)
    for net, x, y in wanted:
        tag = {"source": "real", "network": net, "x": x, "y": y}
        if net not in parsed or net not in knowledge:
            yield {**tag, "status": "unavailable"}, []
            continue
        cpdag, dag = parsed[net]["cpdag"], parsed[net]["dag"]
        g0 = apply_orientations(cpdag, list(knowledge[net]))
        if g0 is None:
            yield {**tag, "status": "inconsistent"}, []
            continue
        z = optimal_adjustment_set_mpdag(g0, x, y)
        if z is None:
            yield {**tag, "status": "degenerate"}, []
            continue
        rng = np.random.default_rng(ft.derive_seed(seed, net, x, y))
        ctx = make_context(random_sem(dag, rng), cpdag, x, y, frozenset(z))
        if abs(ctx.theta_z) < 1e-12:
            yield {**tag, "status": "zero_estimate"}, []
            continue
        signal.alarm(budget)
        try:
            yield analyse(tag, cpdag, g0, ctx, differential=differential)
        except Timeout:
            yield {**tag, "status": "timeout", "budget": budget}, []
        except (MemoryError, ValueError) as exc:
            yield {**tag, "status": f"error:{type(exc).__name__}"}, []
        finally:
            signal.alarm(0)


def synthetic_instances(n_values, n_seeds, seed_base, budget, differential,
                        corruptions=(0.0,), min_knowledge=1, min_space=1):
    """Yield ``(meta, records)`` for generated instances.

    Uses the counterexample module's own accepted-instance builder, which is
    already wired to the corrected space and applies the same acceptance gate
    the published counterexample search used.
    """
    from bkrobust.epsilon.counterexamples import GEN_PARAM_GRID, try_build_candidate

    idx = 0

    for n in n_values:
        for gen, grids in sorted(GEN_PARAM_GRID.items()):
            for gp in grids:
              for corruption in corruptions:
                for s in range(n_seeds):
                    # zlib.crc32, not hash(): the builtin is salted by
                    # PYTHONHASHSEED and would make the sweep irreproducible
                    # across processes.
                    gen_salt = zlib.crc32(f"{gen}:{sorted(gp.items())}".encode()) % 997
                    seed = seed_base + 1000 * n + 7 * s + gen_salt
                    cand, why = try_build_candidate(gen, n, seed, gp, 0.6, corruption)
                    if cand is None:
                        continue
                    # Skip instances that cannot discriminate: a one-orientation
                    # knowledge set has a two-shell profile and a three-element
                    # space, so every rule trivially agrees on it.
                    if len(cand.k0) < min_knowledge or len(cand.space) < min_space:
                        continue
                    idx += 1
                    tag = {
                        "source": "synthetic",
                        # The name must separate every instance the sweep builds.
                        # (gen, n, seed) is NOT enough: distinct generator
                        # parameter grids can collide onto one seed, and the
                        # corruption rate varies independently. Both go in, plus
                        # a running index as a last resort -- two instances
                        # sharing a name silently merge into one curve and
                        # manufacture non-monotonicity in quantities that cannot
                        # be non-monotone.
                        "network": f"{gen}_n{n}_s{seed}_c{corruption}_g{gen_salt}_i{idx}",
                        "gen_params": dict(gp),
                        "x": cand.treatment,
                        "y": cand.outcome,
                        "generator": gen,
                        "n": n,
                        "seed": seed,
                        "corruption": corruption,
                    }
                    rng = np.random.default_rng(seed ^ 0x5EED)
                    ctx = make_context(
                        random_sem(cand.dag, rng), cand.cpdag,
                        cand.treatment, cand.outcome, cand.z,
                    )
                    if abs(ctx.theta_z) < 1e-12:
                        continue
                    signal.alarm(budget)
                    try:
                        yield analyse(tag, cand.cpdag, cand.g0, ctx, differential=differential)
                    except Timeout:
                        yield {**tag, "status": "timeout"}, []
                    except (MemoryError, ValueError):
                        continue
                    finally:
                        signal.alarm(0)


def main() -> None:
    """Run the full-space comparison and write its outputs."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--condition", default="D_LLM")
    p.add_argument("--instances", default=",".join(DEFAULT_INSTANCES))
    p.add_argument("--seed", type=int, default=20260917)
    p.add_argument("--budget", type=int, default=600, help="seconds per instance")
    p.add_argument("--synthetic-n", default="4,5,6")
    p.add_argument("--synthetic-seeds", type=int, default=6)
    p.add_argument("--synthetic-seed-base", type=int, default=911000)
    p.add_argument("--synthetic-corruptions", default="0.0,0.25,0.5")
    p.add_argument("--min-knowledge", type=int, default=3)
    p.add_argument("--min-space", type=int, default=20)
    p.add_argument("--no-synthetic", action="store_true")
    p.add_argument("--no-differential", action="store_true")
    p.add_argument("--out", default="results/mean_fullspace")
    args = p.parse_args()

    signal.signal(signal.SIGALRM, _alarm)
    out_dir = ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    metas: list[dict[str, Any]] = []
    records: list[dict[str, Any]] = []

    streams = [real_instances(args.instances.split(","), args.condition, args.seed,
                              args.budget, not args.no_differential)]
    if not args.no_synthetic:
        streams.append(synthetic_instances(
            tuple(int(v) for v in args.synthetic_n.split(",")),
            args.synthetic_seeds, args.synthetic_seed_base,
            args.budget, not args.no_differential,
            corruptions=tuple(float(v) for v in args.synthetic_corruptions.split(",")),
            min_knowledge=args.min_knowledge, min_space=args.min_space,
        ))

    for stream in streams:
        for meta, recs in stream:
            metas.append(meta)
            records.extend(recs)
            print(
                f"[fullspace] {meta.get('source')} {meta['network']} "
                f"{meta['x']}->{meta['y']} status={meta['status']} "
                f"|space|={meta.get('space_size')} |K_G0|={meta.get('n_knowledge')} "
                f"secs={meta.get('seconds')}",
                flush=True,
            )

    with (out_dir / "shells.jsonl").open("w") as f:
        for r in records:
            f.write(json.dumps(r) + "\n")
    (out_dir / "instances.json").write_text(
        json.dumps({"args": vars(args), "instances": metas}, indent=2) + "\n"
    )
    ok = sum(1 for m in metas if m["status"] == "ok")
    print(f"[fullspace] {ok}/{len(metas)} instances ok; wrote {out_dir}")


if __name__ == "__main__":
    main()
