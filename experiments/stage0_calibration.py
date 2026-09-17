"""Stage 0: what each instrument certifies, and whether the truth honoured it.

Every instrument answers from observables alone on a corrupted knowledge set.
The truth is revealed once, after all of them have answered, and only to score.
An instrument that certifies more moves than the analyst was wrong by, on a
query whose committed adjustment set is invalid at the truth, has
over-certified; that rate is its validity. The mean number of moves it
certifies is its sharpness. Neither axis alone decides anything, which is why
both are reported.

Design and predictions: ``results/stage0/PREREGISTRATION.md``, written first.

Coverage is a prefix of one seeded permutation per network rather than the
even stride of ``benchmarks.measure.select_knowledge``, because the stride is
not nested and a coverage sweep over it is therefore not a retraction
sequence. At coverage 1.0 the two agree, which is the differential test.

Writes ``results/stage0/calibration_rows.csv`` and
``results/stage0/calibration_summary.json``.

    python experiments/stage0_calibration.py --smoke
    python experiments/stage0_calibration.py --reps 5
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import itertools
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from bkrobust.benchmarks.describe import parse_file  # noqa: E402
from bkrobust.benchmarks.measure import (  # noqa: E402
    MAX_G0_UNDIRECTED_FOR_EXTENSIONS,
    component_of,
    separation,
)
from bkrobust.demo.evaluate import (  # noqa: E402
    is_valid_adjustment_set_dag,
    optimal_adjustment_set_mpdag,
)
from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag  # noqa: E402
from bkrobust.hybrid import breakdown_radius  # noqa: E402
from bkrobust.mpdag_criterion.criterion import is_amenable  # noqa: E402
from bkrobust.synth.knowledge import flip  # noqa: E402

INSTANCES = ROOT / "results" / "axisa3" / "instances.jsonl"
MODELS = ROOT / "results" / "axisa3" / "networks" / "example_models"
OUT = ROOT / "results" / "stage0"

#: Depth of the exhaustive retraction enumeration. Above it the row carries a
#: ``>3`` sentinel rather than a guess.
MAX_DEPTH = 3
#: Instruments, in ladder order. Every one reads observables only.
RUNGS = (
    "B0_always",
    "B1_component",
    "B2_free_rule",
    "B3_min_s_kg0",
    "B4_kg0",
    "B5_k",
    "B6_hop",
    "B7_claim",
)


def load(name: str) -> tuple[MPDAG, MPDAG]:
    """Parse one network file and return its DAG and CPDAG."""
    path = next(p for p in sorted(MODELS.iterdir()) if p.name.startswith(name))
    parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
    dag = MPDAG(parsed.nodes, directed=parsed.edges)
    return dag, dag_to_cpdag(dag)


def nested_selection(k_true: list, coverage: float, seed: int) -> list:
    """A prefix of one seeded permutation, so lower coverage is a subset of higher."""
    if coverage >= 1.0:
        return list(k_true)
    keep = round(len(k_true) * coverage)
    if keep <= 0:
        return []
    order = np.random.default_rng(seed).permutation(len(k_true))
    return [k_true[int(i)] for i in sorted(order[:keep])]


def breaking_depths(
    cpdag: MPDAG, k: list, x: str, y: str, z: frozenset, budget: int
) -> tuple[int | None, str, int | None, str, int]:
    """Smallest retraction that breaks validity, and the one that breaks identifiability.

    Both are found in a single walk over retraction subsets, because the
    expensive step is the Meek re-closure and both predicates read the same
    closed graph. Failure is upward-closed for each, so the first depth at
    which any subset breaks a predicate is that predicate's answer.

    ``budget`` caps the total number of subsets closed for this row. A depth
    that would exceed it is not entered, and the row carries a bracket rather
    than a guess: enumerating depth 3 of a twenty-six-claim set means 2,600
    closures of a four-hundred-node graph.

    Returns:
        ``(r_claim, claim_status, r_id, id_status, subsets_closed)``.
    """
    r_claim = r_id = None
    claim_status = id_status = f"gt_{min(len(k), MAX_DEPTH)}"
    used = 0
    for depth in range(1, MAX_DEPTH + 1):
        if depth > len(k):
            break
        width = math.comb(len(k), depth)
        if used + width > budget:
            if r_claim is None:
                claim_status = f"censored_at_depth_{depth}"
            if r_id is None:
                id_status = f"censored_at_depth_{depth}"
            break
        for subset in itertools.combinations(range(len(k)), depth):
            drop = set(subset)
            g = apply_orientations(cpdag, [e for i, e in enumerate(k) if i not in drop])
            used += 1
            if g is None:
                continue
            if r_claim is None and not is_gac_valid_mpdag(g, x, y, z):
                r_claim, claim_status = depth, "exact"
            if r_id is None and not is_amenable(g, x, y):
                r_id, id_status = depth, "exact"
            if r_claim is not None and r_id is not None:
                break
        if r_claim is not None and r_id is not None:
            break
    return r_claim, claim_status, r_id, id_status, used


def bucket(r: int | None, w_c: int, z_valid_at_truth: bool) -> str:
    """The four-way audit verdict for one instrument's answer on one row."""
    if r is None:  # certifies everything inside the budget
        covered = True
    else:
        covered = w_c <= r - 1
    if covered:
        return "dangerous" if not z_valid_at_truth else "held"
    return "safe" if not z_valid_at_truth else "slack"


def cluster_bootstrap(by_network: dict[str, list[float]], reps: int, seed: int) -> list[float]:
    """Percentile interval for a mean, resampling networks rather than rows."""
    names = sorted(by_network)
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(reps):
        pick = rng.integers(0, len(names), size=len(names))
        vals = [v for i in pick for v in by_network[names[int(i)]]]
        if vals:
            draws.append(sum(vals) / len(vals))
    draws.sort()
    if not draws:
        return [float("nan"), float("nan")]
    return [draws[int(0.025 * len(draws))], draws[int(0.975 * len(draws)) - 1]]


def _stable(key: str) -> int:
    """A seed component that does not depend on the interpreter's hash salt."""
    return int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)


def main(argv: list[str] | None = None) -> None:
    """Run the sweep, score every instrument against the truth, write the summary."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=5, help="replicates per corruption rate above 0")
    ap.add_argument("--rates", default="0.0,0.1,0.25,0.5")
    ap.add_argument("--coverages", default="1.0,0.5")
    ap.add_argument("--seed", type=int, default=20260912)
    ap.add_argument("--skip", default="pathfinder")
    ap.add_argument("--smoke", action="store_true", help="three networks, one replicate")
    ap.add_argument("--bootstrap", type=int, default=4000)
    ap.add_argument(
        "--subset-budget",
        type=int,
        default=200,
        help="retraction subsets closed per row before the radius carries a bracket",
    )
    ap.add_argument(
        "--max-undirected",
        type=int,
        default=8,
        help=(
            "cap on undirected edges in G0 before the row is recorded as "
            "extensions_intractable. Enumerating the DAG extensions of a corrupted G0 is "
            "exponential in this count, and corruption leaves more edges undirected than a "
            "truthful set does, so the wall is reached here and not in the committed sweep, "
            f"whose guard is {MAX_G0_UNDIRECTED_FOR_EXTENSIONS}."
        ),
    )
    args = ap.parse_args(argv)

    skip = {s for s in args.skip.split(",") if s}
    rates = [float(v) for v in args.rates.split(",")]
    coverages = [float(v) for v in args.coverages.split(",")]

    rows = [json.loads(line) for line in INSTANCES.open()]
    adm = [r for r in rows if r.get("admissible") and r["network"] not in skip]
    pairs = sorted({(r["network"], r["x"], r["y"]) for r in adm})
    if args.smoke:
        keep = {"child", "asia", "Schipf_2010"}
        pairs = [p for p in pairs if p[0] in keep]
        args.reps = 1

    networks = sorted({p[0] for p in pairs})
    graphs: dict[str, tuple[MPDAG, MPDAG]] = {}
    recovering: dict[str, list] = {}
    for n in networks:
        graphs[n] = load(n)
        recovering[n] = sorted(knowledge_to_recover(*graphs[n]))
    print(f"{len(pairs)} pairs over {len(networks)} networks", flush=True)

    OUT.mkdir(parents=True, exist_ok=True)
    fh = (OUT / "calibration_rows.csv").open("w", newline="")
    writer = csv.writer(fh)
    writer.writerow(
        [
            "network",
            "x",
            "y",
            "coverage",
            "flip_rate",
            "rep",
            "status",
            "n_k",
            "k_g0",
            "n_wrong_claims",
            "separation",
            "separation_status",
            "component_size",
            "z_size",
            "g0_undirected",
            "r_hop",
            "hop_method",
            "r_claim",
            "r_claim_status",
            "r_id",
            "r_id_status",
            "phi_1",
            "subsets_closed",
            "z_valid_at_truth",
            *[f"bucket_{r}" for r in RUNGS],
            *[f"cert_{r}" for r in RUNGS],
        ]
    )

    scored: list[dict] = []
    status_count: collections.Counter = collections.Counter()
    t0 = time.perf_counter()
    done = 0
    for network, x, y in pairs:
        dag, cpdag = graphs[network]
        k_true = recovering[network]
        for coverage in coverages:
            base = nested_selection(k_true, coverage, args.seed + _stable(network) % 10_000)
            if not base:
                status_count["empty_selection"] += 1
                continue
            for rate in rates:
                n_reps = 1 if rate == 0.0 else args.reps
                for rep in range(n_reps):
                    done += 1
                    rng = np.random.default_rng(
                        [
                            args.seed,
                            _stable(f"{network}|{x}|{y}"),
                            int(coverage * 100),
                            int(rate * 100),
                            rep,
                        ]
                    )
                    k = flip(base, rng, rate) if rate > 0 else sorted(base)
                    g0 = apply_orientations(cpdag, k)
                    if g0 is None:
                        status_count["knowledge_inconsistent"] += 1
                        writer.writerow(
                            [
                                network,
                                x,
                                y,
                                coverage,
                                rate,
                                rep,
                                "knowledge_inconsistent",
                                len(k),
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                *[""] * (2 * len(RUNGS)),
                            ]
                        )
                        continue
                    n_und = len(g0.undirected_edges)
                    if n_und > args.max_undirected:
                        status_count["extensions_intractable"] += 1
                        continue
                    o = optimal_adjustment_set_mpdag(g0, x, y)
                    if o is None:
                        status_count["o_not_identified"] += 1
                        continue
                    if not o:
                        status_count["o_empty"] += 1
                        continue
                    z = frozenset(o)
                    if not is_gac_valid_mpdag(g0, x, y, z):
                        status_count["z_invalid_at_g0"] += 1
                        continue
                    status_count["ok"] += 1

                    # instruments, observables only
                    sep, sep_status = separation(cpdag, x, z)
                    comp = component_of(cpdag, x)
                    k_g0 = sum(
                        1
                        for (a, b) in g0.directed_edges
                        if tuple(sorted((a, b))) in cpdag.undirected_edges
                    )
                    hop = breakdown_radius(cpdag, None, x, y, z, g0=g0, time_limit_s=120.0)
                    r_claim, r_claim_status, r_id, r_id_status, n_subsets = breaking_depths(
                        cpdag, k, x, y, z, args.subset_budget
                    )
                    n_single_break = sum(
                        1
                        for i in range(len(k))
                        if (gi := apply_orientations(cpdag, [e for j, e in enumerate(k) if j != i]))
                        is not None
                        and not is_gac_valid_mpdag(gi, x, y, z)
                    )
                    phi_1 = n_single_break / len(k)

                    cert = {
                        "B0_always": None,
                        "B1_component": len(comp) if comp else 1,
                        "B2_free_rule": sep if sep_status == "measured" else 1,
                        "B3_min_s_kg0": min(sep, k_g0) if sep_status == "measured" else 1,
                        "B4_kg0": k_g0,
                        "B5_k": len(k),
                        "B6_hop": hop.radius,
                        "B7_claim": r_claim if r_claim is not None else MAX_DEPTH + 1,
                    }

                    # the truth, revealed only now
                    w_c = sum(1 for (a, b) in k if not dag.is_directed_edge(a, b))
                    z_ok = is_valid_adjustment_set_dag(dag, x, y, z)
                    buckets = {r: bucket(cert[r], w_c, z_ok) for r in RUNGS}

                    rec = {
                        "network": network,
                        "coverage": coverage,
                        "flip_rate": rate,
                        "w_c": w_c,
                        "z_valid": z_ok,
                        "r_hop": hop.radius,
                        "r_claim": r_claim,
                        "r_claim_status": r_claim_status,
                        "r_id": r_id,
                        "r_id_status": r_id_status,
                        "phi_1": phi_1,
                        "n_k": len(k),
                        "k_g0": k_g0,
                        "cert": cert,
                        "buckets": buckets,
                    }
                    scored.append(rec)
                    writer.writerow(
                        [
                            network,
                            x,
                            y,
                            coverage,
                            rate,
                            rep,
                            "ok",
                            len(k),
                            k_g0,
                            w_c,
                            sep,
                            sep_status,
                            len(comp) if comp else "",
                            len(z),
                            n_und,
                            hop.radius,
                            hop.method,
                            r_claim,
                            r_claim_status,
                            r_id,
                            r_id_status,
                            round(phi_1, 4),
                            n_subsets,
                            int(z_ok),
                            *[buckets[r] for r in RUNGS],
                            *[cert[r] if cert[r] is not None else "inf" for r in RUNGS],
                        ]
                    )
        if done and len(scored) % 500 < 20:
            print(
                f"  {done} draws, {len(scored)} scored, {time.perf_counter() - t0:.0f} s",
                flush=True,
            )
    fh.close()

    # ---- summaries
    def ladder(subset: list[dict]) -> dict:
        out = {}
        for rung in RUNGS:
            dang = {n: [] for n in {r["network"] for r in subset}}
            sharp = {n: [] for n in dang}
            for r in subset:
                dang[r["network"]].append(1.0 if r["buckets"][rung] == "dangerous" else 0.0)
                c = r["cert"][rung]
                sharp[r["network"]].append(float("inf") if c is None else float(max(c - 1, 0)))
            flat_d = [v for vs in dang.values() for v in vs]
            flat_s = [v for vs in sharp.values() for v in vs if math.isfinite(v)]
            b = collections.Counter(r["buckets"][rung] for r in subset)
            # Conditioned on the truth, which is where the power is: the outcome
            # the certificate exists for happens on a small fraction of rows.
            invalid = b["dangerous"] + b["safe"]
            valid = b["held"] + b["slack"]
            out[rung] = {
                "dangerous_rate": round(sum(flat_d) / len(flat_d), 4) if flat_d else None,
                "dangerous_ci": [
                    round(v, 4) for v in cluster_bootstrap(dang, args.bootstrap, args.seed)
                ],
                "mean_certified_moves": round(sum(flat_s) / len(flat_s), 3) if flat_s else None,
                "detection_rate": round(b["safe"] / invalid, 4) if invalid else None,
                "false_alarm_rate": round(b["slack"] / valid, 4) if valid else None,
                "buckets": dict(b),
            }
        return out

    corrupted = [r for r in scored if r["flip_rate"] > 0]
    summary = {
        "draws": done,
        "scored": len(scored),
        "status": dict(status_count),
        "networks_scored": len({r["network"] for r in scored}),
        "settings": {
            "max_undirected": args.max_undirected,
            "subset_budget": args.subset_budget,
            "library_guard": MAX_G0_UNDIRECTED_FOR_EXTENSIONS,
            "rates": rates,
            "coverages": coverages,
            "reps": args.reps,
            "seed": args.seed,
            "max_depth": MAX_DEPTH,
            "skipped": sorted(skip),
        },
        "r_claim_distribution": dict(
            sorted(
                collections.Counter(
                    (r["r_claim"] if r["r_claim"] is not None else r["r_claim_status"])
                    for r in scored
                ).items(),
                key=lambda kv: str(kv[0]),
            )
        ),
        "r_claim_by_rate": {
            str(rate): dict(
                sorted(
                    collections.Counter(
                        (r["r_claim"] if r["r_claim"] is not None else r["r_claim_status"])
                        for r in scored
                        if r["flip_rate"] == rate
                    ).items(),
                    key=lambda kv: str(kv[0]),
                )
            )
            for rate in rates
        },
        "r_id_distribution": dict(
            sorted(
                collections.Counter(
                    (r["r_id"] if r["r_id"] is not None else r["r_id_status"]) for r in scored
                ).items(),
                key=lambda kv: str(kv[0]),
            )
        ),
        "leverage": {},
        "invariants": {
            "r_claim_gt_r_hop": sum(
                1 for r in scored if r["r_claim"] is not None and r["r_claim"] > r["r_hop"]
            ),
            "r_id_lt_r_claim": sum(
                1
                for r in scored
                if r["r_id"] is not None and r["r_claim"] is not None and r["r_id"] < r["r_claim"]
            ),
        },
        "ladder_all_scored": ladder(scored),
        "ladder_corrupted_only": ladder(corrupted) if corrupted else {},
        "ladder_by_rate": {
            str(rate): ladder([r for r in scored if r["flip_rate"] == rate])
            for rate in rates
            if any(r["flip_rate"] == rate for r in scored)
        },
        "z_invalid_at_truth": sum(1 for r in scored if not r["z_valid"]),
        "z_invalid_at_truth_by_rate": {
            str(rate): sum(1 for r in scored if r["flip_rate"] == rate and not r["z_valid"])
            for rate in rates
        },
        "seconds": round(time.perf_counter() - t0, 1),
    }
    lev = [r["r_hop"] / r["r_claim"] for r in scored if r["r_claim"]]
    if lev:
        lev.sort()
        summary["leverage"] = {
            "min": round(lev[0], 2),
            "median": round(lev[len(lev) // 2], 2),
            "max": round(lev[-1], 2),
            "above_2": sum(1 for v in lev if v > 2),
        }
    (OUT / "calibration_summary.json").write_text(json.dumps(summary, indent=1))

    print(f"\nscored {len(scored)} of {done} draws in {summary['seconds']} s")
    print("status", json.dumps(summary["status"]))
    print("r_claim by rate", json.dumps(summary["r_claim_by_rate"]))
    print("r_id", json.dumps(summary["r_id_distribution"]))
    print("leverage", json.dumps(summary["leverage"]))
    print("invariants", json.dumps(summary["invariants"]))
    inv = sum(1 for r in scored if not r["z_valid"])
    print(f"\nrows whose committed set is invalid at the truth: {inv} of {len(scored)}")
    print("ladder on corrupted rows:")
    for rung, v in (summary["ladder_corrupted_only"] or summary["ladder_all_scored"]).items():
        print(
            f"  {rung:15s} miss {v['dangerous_rate']} {v['dangerous_ci']}  "
            f"detect {v['detection_rate']}  false-alarm {v['false_alarm_rate']}  "
            f"moves {v['mean_certified_moves']}"
        )


if __name__ == "__main__":
    main()
