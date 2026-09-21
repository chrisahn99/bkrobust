r"""The mean bias profile on the paper's own geometry, not the retraction up-set.

The perturbation space this project's geometry is defined on is

    G_Chat = { G a valid MPDAG : [G] subset-or-equal [Chat] }

ordered by model inclusion, with distance the BFS hop count on the covering
graph. That is the space ``r_val`` is defined on, and Definition 5 of
``docs/R_EPSILON_THEORY.md`` defines ``r_eps`` on it too.

The **retraction up-set** ``^G0`` is a strictly smaller object. It entered as a
computational device and is *licensed* for the worst case by a theorem chain
that is now closed (Anti-Exchange PROVED, hence Lemma R, Property S, Lemma L,
and Conjecture 2 a PROVED THEOREM): the first crossing of the retraction profile
and of the full-space profile coincide, so ``r_val`` and the max-based
``r_eps`` computed over retractions are the same numbers the full space gives.
That equivalence was also checked directly here: over 79 instances the shell
maxima agreed with the ball maxima on every shell and the first crossings agreed
on all 948 cells (``results/mean_fullspace/``).

**No such equivalence holds for a mean**, and this script therefore does not use
the up-set at all. It enumerates the whole space, takes BFS distances from
``G0``, and averages ``B`` over the sphere at each distance:

    mu(d) = mean { B(G) : G in G_Chat, d(G0, G) = d }

so that every number reported sits in the same geometry as ``r_val``. An
earlier version of this statistic averaged over the retraction shell; that is a
different quantity and the two disagree in both directions
(``docs/RESULT_MEAN_SOUNDNESS_GATE.md``). This file supersedes it for reporting.

Two consequences of staying in the real geometry, both of which are reported
rather than smoothed over:

* **Cost.** The space is ``3^k`` in the ``k`` undirected edges of ``Chat``, and
  building the covering relation is quadratic in its size. The corpus is *not*
  uniformly enumerable: ``k`` runs from 1 to 25. Instances whose space cannot be
  built within the budget are reported as such and excluded, with the
  denominator stated. Nothing is estimated or sampled.
* **Shape.** ``mu`` over the full sphere is **not** monotone in ``d``. Far
  shells contain heavily oriented states that pin the effect down and so carry
  little bias, which pulls the average back down. A first crossing is therefore
  an *onset* -- the smallest perturbation size at which the expected bias
  exceeds the tolerance -- and not a threshold beyond which it stays exceeded.
  Both the onset and whether it persists are recorded.

``r_val`` is graph-theoretic. Every bias number here is **SEM-conditional**: the
corpus ships graph structure only, so ``B`` is evaluated against a seeded
linear-Gaussian SEM attached to the ground-truth DAG, seeded per instance
exactly as ``final_table.py`` does.

Usage::

    PYTHONPATH=src .venv/bin/python experiments/mean_table_fullspace.py [options]
"""

from __future__ import annotations

import argparse
import json
import signal
import statistics
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.core.spacelib import distances_from  # noqa: E402
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag, random_sem  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.epsilon.bias import bias_at, make_context  # noqa: E402
from bkrobust.search.space_fixed import build_space_correct  # noqa: E402

import final_table as ft  # noqa: E402

#: Reporting grid, in relative units: fractions of the reported effect's own
#: magnitude. Chosen in ``docs/MEAN_REPORTING.md`` section 4 and kept here so the
#: two tables are directly comparable. It deliberately straddles rather than
#: sits on 1.0, where half the instances pile up exactly.
DEFAULT_EPSILONS = "0.02,0.10,0.25,0.40,0.60,0.90"


class Timeout(Exception):
    """Raised when one network's space exceeds its build budget."""


def _alarm(_signum: int, _frame: Any) -> None:
    raise Timeout()


def profile_for(
    cpdag: Any, g0: Any, ctx: Any, space: Any, dists: dict[Any, int], scale: float
) -> list[dict[str, Any]]:
    """The full-space bias profile by BFS distance from ``G0``.

    Args:
        cpdag: The CPDAG.
        g0: The analyst's state.
        ctx: The bias context.
        space: The enumerated space.
        dists: BFS distance per reachable element.
        scale: Divide biases by this for relative units.

    Returns:
        One record per distance, in increasing ``d``.
    """
    b = {g: bias_at(ctx, g, method="enumerate").worst / scale for g in space.elements}
    reachable = [g for g in space.elements if dists.get(g) is not None]
    out: list[dict[str, Any]] = []
    running_max = 0.0
    for d in range(max(dists[g] for g in reachable) + 1):
        sphere = [b[g] for g in reachable if dists[g] == d]
        if not sphere:
            continue
        ball = [b[g] for g in reachable if dists[g] <= d]
        running_max = max(running_max, max(sphere))
        nz = [v for v in sphere if v > 0]
        out.append({
            "d": d,
            "n_sphere": len(sphere),
            "n_ball": len(ball),
            "mu": statistics.fmean(sphere),
            "mu_ball": statistics.fmean(ball),
            "sphere_max": max(sphere),
            "beta": running_max,          # max over the ball: the r_eps profile
            "frac_nonzero": len(nz) / len(sphere),
            "median": statistics.median(sphere),
        })
    return out


def r_onset(profile: list[dict[str, Any]], eps: float, key: str = "mu") -> int:
    """First distance whose statistic exceeds ``eps``, else ``UNREACHED``."""
    for rec in profile:
        if rec[key] > eps:
            return rec["d"]
    return UNREACHED


def persists(profile: list[dict[str, Any]], eps: float, key: str = "mu") -> bool | None:
    """Whether the statistic stays above ``eps`` from its onset to the outermost shell.

    ``mu`` is not monotone on this space, so a crossing need not be permanent.
    Returns ``None`` when the tolerance is never crossed.
    """
    hit = r_onset(profile, eps, key)
    if hit == UNREACHED:
        return None
    return all(rec[key] > eps for rec in profile if rec["d"] >= hit)


def main() -> None:
    """Enumerate the full space per instance and write the profile table."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--condition", default="D_LLM")
    p.add_argument("--panel", default="results/final_table_eps_gt1/instances.jsonl")
    p.add_argument("--epsilons", default=DEFAULT_EPSILONS)
    p.add_argument("--seed", type=int, default=20260917)
    p.add_argument("--budget", type=int, default=900,
                   help="seconds allowed to BUILD one network's space")
    p.add_argument("--max-undirected", type=int, default=10,
                   help="skip a network whose CPDAG has more undirected edges than this; "
                        "3^k is a hard lower bound on the enumeration, so burning the full "
                        "budget on a hopeless network only wastes wall clock")
    p.add_argument("--networks", default="", help="comma list; default all in the panel")
    p.add_argument("--out", default="results/mean_table_fullspace")
    args = p.parse_args()

    signal.signal(signal.SIGALRM, _alarm)
    epsilons = [float(v) for v in args.epsilons.split(",")]
    panel = [json.loads(line) for line in (ROOT / args.panel).open()]
    informative = [r for r in panel
                   if r["status"] == "ok" and r["r_val"] not in (0, None)]
    want = {v for v in args.networks.split(",") if v} or {r["network"] for r in informative}
    knowledge = ft.load_knowledge(args.condition)
    parsed = ft.load_networks(names=want, skip=set(), max_nodes=0)

    out_dir = ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, Any]] = []
    t_all = time.perf_counter()

    # The space depends only on the CPDAG, so build it once per network and
    # reuse it across that network's queries -- it is by far the dominant cost.
    by_net: dict[str, list[dict[str, Any]]] = {}
    for r in informative:
        by_net.setdefault(r["network"], []).append(r)

    for net in sorted(by_net):
        rows = by_net[net]
        if net not in parsed or net not in knowledge:
            for r in rows:
                records.append({**{k: r[k] for k in ("network", "x", "y", "r_val")},
                                "status": "unavailable"})
            continue
        cpdag, dag = parsed[net]["cpdag"], parsed[net]["dag"]
        k_und = len(cpdag.undirected_edges)
        if k_und > args.max_undirected:
            print(f"[fullspace] {net}: k={k_und} (3^k={3**k_und:,}) above --max-undirected "
                  f"{args.max_undirected} -- {len(rows)} instances excluded, not attempted",
                  flush=True)
            for r in rows:
                records.append({**{q: r[q] for q in ("network", "x", "y", "r_val")},
                                "status": "space_too_large", "k_undirected": k_und})
            continue
        signal.alarm(args.budget)
        t0 = time.perf_counter()
        try:
            space = build_space_correct(cpdag)
            signal.alarm(0)
        except Timeout:
            signal.alarm(0)
            print(f"[fullspace] {net}: space build exceeded {args.budget}s "
                  f"(k={k_und}, 3^k={3**k_und:,}) -- {len(rows)} instances excluded", flush=True)
            for r in rows:
                records.append({**{q: r[q] for q in ("network", "x", "y", "r_val")},
                                "status": "space_timeout", "k_undirected": k_und,
                                "budget": args.budget})
            continue
        except MemoryError:
            signal.alarm(0)
            for r in rows:
                records.append({**{q: r[q] for q in ("network", "x", "y", "r_val")},
                                "status": "space_memory", "k_undirected": k_und})
            continue
        build_s = time.perf_counter() - t0
        print(f"[fullspace] {net}: |space|={len(space):,} k={k_und} built in {build_s:.1f}s "
              f"({len(rows)} instances)", flush=True)

        for r in rows:
            x, y = r["x"], r["y"]
            rec: dict[str, Any] = {
                "network": net, "x": x, "y": y, "r_val": r["r_val"],
                "n_knowledge": r["n_knowledge"], "k_undirected": k_und,
                "space_size": len(space), "space_build_seconds": round(build_s, 2),
            }
            g0 = apply_orientations(cpdag, list(knowledge[net]))
            if g0 is None or g0 not in space.neighbours:
                rec["status"] = "g0_not_in_space"
                records.append(rec); continue
            z = optimal_adjustment_set_mpdag(g0, x, y)
            if z is None:
                rec["status"] = "degenerate"
                records.append(rec); continue
            ctx = make_context(
                random_sem(dag, np.random.default_rng(ft.derive_seed(args.seed, net, x, y))),
                cpdag, x, y, frozenset(z))
            scale = abs(ctx.theta_z)
            if scale < 1e-12:
                rec["status"] = "zero_estimate"
                records.append(rec); continue
            t1 = time.perf_counter()
            dists = distances_from(space, g0)
            prof = profile_for(cpdag, g0, ctx, space, dists, scale)
            mus = [p_["mu"] for p_ in prof]
            rec.update({
                "status": "ok",
                "theta_z": ctx.theta_z,
                "max_distance": max(p_["d"] for p_ in prof),
                "n_reachable": sum(p_["n_sphere"] for p_ in prof),
                "beta_top": max(p_["sphere_max"] for p_ in prof),
                "mu_peak": max(mus),
                "mu_peak_at": prof[mus.index(max(mus))]["d"],
                "mu_monotone": all(mus[i] >= mus[i - 1] - 1e-12 for i in range(1, len(mus))),
                "r_mu": {str(e): r_onset(prof, e) for e in epsilons},
                "r_mu_persists": {str(e): persists(prof, e) for e in epsilons},
                "r_mu_ball": {str(e): r_onset(prof, e, "mu_ball") for e in epsilons},
                "r_beta": {str(e): r_onset(prof, e, "beta") for e in epsilons},
                "profile": prof,
                "seconds": round(time.perf_counter() - t1, 2),
            })
            records.append(rec)
            print(f"[fullspace]   {net} {x}->{y} r_val={r['r_val']} "
                  f"r_mu={[rec['r_mu'][str(e)] for e in epsilons]} "
                  f"monotone={rec['mu_monotone']}", flush=True)

    with (out_dir / "instances.jsonl").open("w") as f:
        for rec in records:
            f.write(json.dumps(rec) + "\n")
    ok = [r for r in records if r.get("status") == "ok"]
    summary = {
        "args": vars(args), "epsilons": epsilons,
        "n_informative_in_panel": len(informative),
        "n_ok": len(ok),
        "n_excluded_space_timeout": sum(1 for r in records if r.get("status") == "space_timeout"),
        "n_excluded_too_large": sum(1 for r in records if r.get("status") == "space_too_large"),
        "networks_excluded": sorted({r["network"] for r in records
                                     if r.get("status") in ("space_timeout", "space_too_large")}),
        "n_non_monotone": sum(1 for r in ok if r.get("mu_monotone") is False),
        "seconds": round(time.perf_counter() - t_all, 1),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\n[fullspace] {summary['n_ok']}/{summary['n_informative_in_panel']} informative "
          f"instances enumerated on the full space; "
          f"{summary['n_excluded_space_timeout'] + summary['n_excluded_too_large']} excluded "
          f"({', '.join(summary['networks_excluded']) or 'none'}); "
          f"non-monotone mu: {summary['n_non_monotone']}; {summary['seconds']}s")


if __name__ == "__main__":
    main()
