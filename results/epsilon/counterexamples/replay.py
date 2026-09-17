"""Recompute every headline number in this directory's JSON files from scratch.

Run from the repository root with the project's virtualenv:

    .venv/bin/python results/epsilon/counterexamples/replay.py

This does not re-run the search: it reads the persisted witnesses (CPDAG, G0,
SEM weights/noise variances, the two witness states, RNG seeds) and recomputes
the numbers those files claim, using only the frozen evaluators
(:mod:`bkrobust.epsilon.bias`, :mod:`bkrobust.core.oracle`,
:mod:`bkrobust.core.spacelib`, :mod:`bkrobust.epsilon.profile`). A referee who
does not trust the search only needs to trust this script and the four JSON
files next to it.

Exits non-zero if any recomputed number disagrees with what was stored.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
# Locate the repo root (the directory containing "src") by walking upward.
ROOT = HERE
while not (ROOT / "src").is_dir():
    if ROOT.parent == ROOT:
        raise SystemExit("could not locate the repository root (no 'src' directory found)")
    ROOT = ROOT.parent
sys.path.insert(0, str(ROOT / "src"))

from bkrobust.epsilon import counterexamples as ce  # noqa: E402
from bkrobust.epsilon.bias import bias_at, make_context  # noqa: E402
from bkrobust.core.spacelib import distances_from, radius  # noqa: E402
from bkrobust.epsilon.profile import DEFAULT_SHELL_CAP, retraction_shell  # noqa: E402
from bkrobust.search.space_fixed import build_space_correct  # noqa: E402
from bkrobust.search.space_fixed import knowledge_of as space_knowledge_of  # noqa: E402

FAILURES: list[str] = []


def check(label: str, ok: bool, detail: str = "") -> None:
    status = "OK" if ok else "FAIL"
    print(f"[{status}] {label}{(': ' + detail) if detail and not ok else ''}")
    if not ok:
        FAILURES.append(label)


# --------------------------------------------------------------------------
# C1 -- individual-state bias is not monotone in raw distance from G0
# --------------------------------------------------------------------------


def replay_c1() -> None:
    print("\n=== C1 replay ===")
    doc = json.loads((HERE / "c1_witness.json").read_text())
    if not doc["found"]:
        print("(no C1 witness was found; nothing to replay)")
        return
    w = doc["witness"]
    inst = w["instance"]
    cpdag = ce.mpdag_from_json(inst["cpdag"])
    g0 = ce.mpdag_from_json(inst["g0"])
    x, y = inst["treatment"], inst["outcome"]
    z = frozenset(inst["z"])
    sem = ce.sem_from_json(w["sem"])
    ctx = make_context(sem, cpdag, x, y, z)

    g = ce.mpdag_from_json(w["state_G"])
    h = ce.mpdag_from_json(w["state_H"])

    b_g = bias_at(ctx, g, method="semilocal").worst
    b_h = bias_at(ctx, h, method="semilocal").worst
    b_g_enum = bias_at(ctx, g, method="enumerate").worst
    b_h_enum = bias_at(ctx, h, method="enumerate").worst

    check("B(G) matches stored value", abs(b_g - w["state_G"]["B"]) < 1e-9, f"{b_g} vs {w['state_G']['B']}")
    check("B(H) matches stored value", abs(b_h - w["state_H"]["B"]) < 1e-9, f"{b_h} vs {w['state_H']['B']}")
    check("B(G) semilocal == enumerate (Theorem E)", abs(b_g - b_g_enum) < 1e-9)
    check("B(H) semilocal == enumerate (Theorem E)", abs(b_h - b_h_enum) < 1e-9)

    margin = b_g - b_h
    check(
        "margin matches stored value",
        abs(margin - w["margin_B_G_minus_B_H"]) < 1e-9,
        f"{margin} vs {w['margin_B_G_minus_B_H']}",
    )
    check("margin exceeds MARGIN_TOL", margin > ce.MARGIN_TOL, f"margin={margin}")

    # Recompute the BFS distances from G0 independently and check the claimed
    # d(G) < d(H) inequality that makes this a distance -- not order --
    # counterexample.
    space = build_space_correct(cpdag)
    dists = distances_from(space, g0)
    d_g, d_h = dists[g], dists[h]
    check("d(G0,G) matches stored distance", d_g == w["state_G"]["distance_from_G0"])
    check("d(G0,H) matches stored distance", d_h == w["state_H"]["distance_from_G0"])
    check("d(G0,G) < d(G0,H)", d_g < d_h, f"d_g={d_g}, d_h={d_h}")

    print(
        f"C1: n={inst['n']} nodes, k={inst['n_undirected_cpdag']} undirected edges in the CPDAG, "
        f"d(G)={d_g} B(G)={b_g:.6g}, d(H)={d_h} B(H)={b_h:.6g}, margin={margin:.6g}"
    )


# --------------------------------------------------------------------------
# C2 -- mean_abs_bias is not monotone under model inclusion
# --------------------------------------------------------------------------


def replay_c2_rate_sweep() -> None:
    print("\n=== C2 rate-sweep replay (recomputes the summary counts from per_instance rows) ===")
    doc = json.loads((HERE / "c2_rate_sweep.json").read_text())
    total_pairs = sum(row["n_pairs"] for row in doc["per_instance"])
    total_viol_mean = sum(row["n_violations_mean_abs_bias"] for row in doc["per_instance"])
    total_viol_b = sum(row["n_violations_B"] for row in doc["per_instance"])
    check("n_pairs_total matches sum over instances", total_pairs == doc["n_pairs_total"])
    check(
        "n_pairs_violating_mean_abs_bias matches sum over instances",
        total_viol_mean == doc["n_pairs_violating_mean_abs_bias"],
    )
    check("n_pairs_violating_B matches sum over instances", total_viol_b == doc["n_pairs_violating_B"])
    check("B has zero violations (sanity control / Theorem B)", doc["n_pairs_violating_B"] == 0)
    print(
        f"C2 rate sweep: {doc['n_instances']} instances, {doc['n_pairs_total']} ordered "
        f"comparable pairs; mean_abs_bias violation rate "
        f"{doc['violation_rate_mean_abs_bias']:.1%}, B violation rate "
        f"{doc['violation_rate_B']:.1%}."
    )
    print(
        "(Full re-evaluation of every pair's mean_abs_bias is not repeated here -- that IS "
        "the search, not a spot check on a witness; the per-instance violation counts stored "
        "in this file are what was produced by bkrobust.epsilon.counterexamples.run_c2_rate_sweep "
        "and are re-derivable by calling that function again with the same n_values/n_seeds.)"
    )


def replay_c2_flagship() -> None:
    print("\n=== C2 flagship replay (large n_draws, independent repetitions) ===")
    doc = json.loads((HERE / "c2_flagship.json").read_text())
    if not doc.get("found"):
        print("(no C2 flagship witness was found; nothing to replay)")
        return
    result = ce.replay_c2_flagship(doc)
    check(
        "reproduction count matches stored value",
        result["n_repetitions_reproducing_sign"] == doc["n_repetitions_reproducing_sign"],
        f"{result['n_repetitions_reproducing_sign']} vs {doc['n_repetitions_reproducing_sign']}",
    )
    for i, (mg, mh, mg0, mh0) in enumerate(
        zip(
            result["mean_abs_bias_G_per_rep"],
            result["mean_abs_bias_H_per_rep"],
            doc["mean_abs_bias_G_per_rep"],
            doc["mean_abs_bias_H_per_rep"],
        )
    ):
        check(f"rep {i}: mean_abs_bias(G) matches stored", abs(mg - mg0) < 1e-9, f"{mg} vs {mg0}")
        check(f"rep {i}: mean_abs_bias(H) matches stored", abs(mh - mh0) < 1e-9, f"{mh} vs {mh0}")
    print(
        f"C2 flagship: reproduced the sign mean_abs_bias(G) > mean_abs_bias(H) in "
        f"{doc['n_repetitions_reproducing_sign']}/{doc['n_repetitions']} independent repetitions "
        f"at n_draws={doc['n_draws_per_repetition']}; smallest margin/SE across repetitions = "
        f"{doc['min_margin_over_SE']:.2f}."
    )


# --------------------------------------------------------------------------
# C2 consequence -- retraction-only search vs brute-force BFS
# --------------------------------------------------------------------------


def replay_c2_consequence() -> None:
    print("\n=== C2 consequence replay ===")
    doc = json.loads((HERE / "c2_consequence.json").read_text())
    if doc["n_mismatches"] == 0:
        print("(no mismatch was found in the searched scope; nothing to replay)")
        return
    m = doc["mismatches"][0]
    inst = m["instance"]
    cpdag = ce.mpdag_from_json(inst["cpdag"])
    g0 = ce.mpdag_from_json(inst["g0"])
    x, y, z = inst["treatment"], inst["outcome"], frozenset(inst["z"])
    eps = m["eps"]
    n_draws = m["n_draws"]
    seed_base = m["mc_seed"][:2]  # (13, cand.seed) -- the per-graph hash tag is recomputed below

    mean_cache: dict = {}

    def mean_of(g) -> float:
        hit = mean_cache.get(g)
        if hit is not None:
            return hit
        seed = (seed_base[0], seed_base[1], ce._stable_hash(g.edge_string()))
        stats = ce._mean_abs_bias(z, g, x, y, seed, n_draws)
        val = stats.mean_abs_bias if stats.n_evaluations > 0 else float("inf")
        mean_cache[g] = val
        return val

    predicate = lambda g: mean_of(g) > eps  # noqa: E731

    r_retract, w_retract = ce.retraction_only_radius(
        ce.Candidate(  # only cpdag/g0/z/treatment/outcome are used by these two calls
            generator=inst["generator"],
            n=inst["n"],
            seed=inst["seed"],
            gen_params=inst["gen_params"],
            knows_fraction=inst["knows_fraction"],
            corruption_rate=inst["corruption_rate"],
            dag=ce.mpdag_from_json(inst["true_dag"]),
            cpdag=cpdag,
            treatment=x,
            outcome=y,
            k_true=[tuple(e) for e in inst["k_true"]],
            k_assumed=[tuple(e) for e in inst["k_assumed"]],
            g0=g0,
            z=z,
            space=build_space_correct(cpdag),
            dists={},
        ),
        predicate,
    )
    space = build_space_correct(cpdag)
    dists = distances_from(space, g0)
    r_brute, witness_brute = radius(space, dists, predicate)
    w_brute = witness_brute.edge_string() if witness_brute is not None else None

    check("retraction-only radius matches stored value", r_retract == m["retraction_only_radius"])
    check("brute-force radius matches stored value", r_brute == m["brute_force_radius"])
    check("the two radii actually disagree (the point of this witness)", r_retract != r_brute)
    print(
        f"C2 consequence: eps={eps:.6g}, retraction-only radius={r_retract} "
        f"(witness={w_retract}), brute-force radius={r_brute} (witness={w_brute}). "
        f"n_undirected_cpdag={inst['n_undirected_cpdag']}, n={inst['n']}."
    )


if __name__ == "__main__":
    replay_c1()
    replay_c2_rate_sweep()
    replay_c2_flagship()
    replay_c2_consequence()

    print()
    if FAILURES:
        print(f"REPLAY FAILED: {len(FAILURES)} check(s) did not reproduce: {FAILURES}")
        raise SystemExit(1)
    print("All replayed numbers match the persisted witnesses.")
