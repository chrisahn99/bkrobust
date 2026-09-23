r"""Negative results for ``docs/R_EPSILON_THEORY.md`` Remark 4 and section 7 (P/C1, C2).

This module owns exactly two deliverables, both explicitly asked for as
*falsifications of naive readings*, not as evidence for the theorems (which are
proved elsewhere and are not what is being tested here):

* **C1** -- individual-state bias ``B`` is not monotone in raw distance from
  ``G0`` (only the ball-maximum ``beta_up`` is monotone; Theorem B is
  monotonicity **along the order**, not along distance -- Remark 4).
* **C2** -- the incumbent ``bias_stats(...).mean_abs_bias``
  (:mod:`bkrobust.core.oracle`) is not monotone **under the order** ``<=``
  (model inclusion), so ``{mean_abs_bias > eps}`` is not an up-set and the
  upward/retraction-only search that the fast machinery of
  :mod:`bkrobust.epsilon.profile` relies on (and that
  :func:`bkrobust.synth.runner._radius_eps` would need if sped up the same way)
  is not exact for it.

Everything here is a *search* over honestly generated synthetic instances,
mirroring the exact screening pipeline of :mod:`bkrobust.synth.runner`
(``GENERATORS``, ``dag_to_cpdag``, ``_pick_treatment_outcome``, ``gate``,
``draw_k_true`` + ``flip``, ``optimal_adjustment_set_mpdag``, ``is_valid``,
``random_sem``), so that a counterexample is a *real* instance the pipeline
would accept, not a hand-built pathology. The corrected, fully enumerated space
(:func:`bkrobust.search.space_fixed.build_space_correct`) is used throughout,
never the chordality-defective ``enumerate_space``.

House rules for this file, per the session brief:

* This file is the only thing under ``src/bkrobust/epsilon/`` this session
  touches. :mod:`bkrobust.epsilon.bias` and :mod:`bkrobust.epsilon.profile` are
  imported, never edited; anything that looks like a bug in them is reported in
  the accompanying ``SUMMARY.md``, not patched here.
* Every headline number is persisted with enough state (CPDAG, ``G0``, the SEM
  down to its weights and noise variances, the two witness states, RNG seeds)
  that ``results/epsilon/counterexamples/replay.py`` recomputes it from the
  JSON alone -- no re-search.
* Tolerances are never widened to manufacture a result. ``MARGIN_TOL`` below is
  the bar a "real" violation must clear, chosen to sit three orders of
  magnitude above :data:`bkrobust.epsilon.bias.ZERO_TOL` (float noise in the
  normal-equation solves), and is applied uniformly.

Run as a script: ``python -m bkrobust.epsilon.counterexamples``.
"""

from __future__ import annotations

import itertools
import json
import time
import zlib
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from bkrobust.core.oracle import BiasStats, bias_stats, is_valid
from bkrobust.core.spacelib import Space, distances_from, radius
from bkrobust.demo.evaluate import (
    LinearSEM,
    optimal_adjustment_set_mpdag,
    random_sem,
)
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.epsilon.bias import (
    ZERO_TOL,
    BiasContext,
    bias_at,
    make_context,
)
from bkrobust.epsilon.profile import DEFAULT_SHELL_CAP, retraction_shell
from bkrobust.search.space_fixed import build_space_correct
from bkrobust.search.space_fixed import knowledge_of as space_knowledge_of
from bkrobust.synth.generators import GENERATORS
from bkrobust.synth.knowledge import check_consistent, draw_k_true, flip
from bkrobust.synth.runner import _pick_treatment_outcome, gate

Node = str
Edge = tuple[str, str]

#: A bias margin below this is float noise (ZERO_TOL is 1e-9; this is 1000x
#: that), never reported as a counterexample.
MARGIN_TOL: float = 1e-6

#: Hard cap on undirected edges so the corrected space stays exhaustively
#: enumerable (3**k closures before dedup), per the task brief.
MAX_UNDIRECTED_EDGES: int = 6

#: Generator parameter grids. ``decoupled_backdoor`` needs n>=7 and is added
#: separately once n reaches that size.
GEN_PARAM_GRID: dict[str, list[dict[str, Any]]] = {
    "erdos_renyi": [{"edge_prob": p} for p in (0.25, 0.35, 0.45, 0.55, 0.65, 0.75)],
    "scale_free": [{"m_attach": m} for m in (1, 2)],
    "block": [
        {"n_blocks": 2, "p_within": pw, "p_between": pb}
        for pw, pb in ((0.7, 0.15), (0.6, 0.3), (0.8, 0.1))
    ],
}

KNOWS_FRACTIONS: tuple[float, ...] = (0.3, 0.5, 0.7, 1.0)
CORRUPTION_RATES: tuple[float, ...] = (0.0, 0.25, 0.5, 0.75, 1.0)


def _stable_hash(s: str) -> int:
    """A process-independent hash of ``s``, for use as RNG-seed material.

    Python's builtin ``hash()`` on strings is salted per-process
    (``PYTHONHASHSEED``) unless explicitly fixed, which would make a seed
    derived from ``hash(g.edge_string())`` -- and hence the exact MC numbers it
    drives -- unreproducible from one interpreter invocation to the next. This
    project's house rule (see e.g. ``random_sem``'s docstring) is that nothing
    on a randomness path may depend on hash randomisation; ``zlib.crc32`` is
    used here for exactly that reason.
    """
    return zlib.crc32(s.encode("utf-8"))


# --------------------------------------------------------------------------
# Instance generation: mirrors bkrobust.synth.runner.run_instance's screening
# --------------------------------------------------------------------------


@dataclass
class Candidate:
    """One gate-accepted instance, with its corrected space already built."""

    generator: str
    n: int
    seed: int
    gen_params: dict[str, Any]
    knows_fraction: float
    corruption_rate: float
    dag: MPDAG
    cpdag: MPDAG
    treatment: str
    outcome: str
    k_true: list[Edge]
    k_assumed: list[Edge]
    g0: MPDAG
    z: frozenset[Node]
    space: Space
    dists: dict[MPDAG, int]

    @property
    def n_undirected(self) -> int:
        """Undirected edges in the CPDAG -- the size knob for the enumerable space."""
        return len(self.cpdag.undirected_edges)

    @property
    def k0(self) -> list[Edge]:
        """The knowledge set ``K_G0``: the orientations ``G0`` adds to the CPDAG."""
        return space_knowledge_of(self.cpdag, self.g0)

    def tag(self) -> dict[str, Any]:
        """A JSON-ready identifier for this instance, for the per-instance tables."""
        return {
            "generator": self.generator,
            "n": self.n,
            "seed": self.seed,
            "gen_params": self.gen_params,
            "knows_fraction": self.knows_fraction,
            "corruption_rate": self.corruption_rate,
            "n_undirected_cpdag": self.n_undirected,
            "n_knowledge_g0": len(self.k0),
            "space_size": len(self.space),
        }


def try_build_candidate(
    generator: str,
    n: int,
    seed: int,
    gen_params: dict[str, Any],
    knows_fraction: float,
    corruption_rate: float,
) -> tuple[Candidate | None, str]:
    """Build one instance exactly as ``run_instance`` would, or report why not.

    Deliberately does not import ``run_instance`` itself: that function is
    wired to the OLD, chordality-defective ``build_space``
    (``bkrobust.core.spacelib.build_space``), and this session's search needs
    the corrected space. Every other piece is reused verbatim from the runner
    module rather than reimplemented.

    Returns:
        ``(candidate, "ok")`` on acceptance, else ``(None, reason)``.
    """
    rng = np.random.default_rng(seed)
    if generator not in GENERATORS:
        return None, "unknown_generator"
    try:
        dag = GENERATORS[generator](n, rng, **gen_params)
    except ValueError as exc:
        return None, f"generator_error:{exc}"
    cpdag = dag_to_cpdag(dag)
    if len(cpdag.undirected_edges) > MAX_UNDIRECTED_EDGES:
        return None, "too_many_undirected"

    treatment, outcome = _pick_treatment_outcome(generator, dag, rng)
    accepted, reason = gate(dag, cpdag, treatment, outcome)
    if not accepted:
        return None, reason

    k_true = draw_k_true(dag, cpdag, rng, knows_fraction)
    k_assumed = flip(k_true, rng, corruption_rate) if corruption_rate > 0 else list(k_true)
    if not check_consistent(cpdag, k_assumed):
        return None, "k_assumed_inconsistent"

    g0 = apply_orientations(cpdag, k_assumed)
    assert g0 is not None
    z = optimal_adjustment_set_mpdag(g0, treatment, outcome)
    if z is None:
        return None, "z_not_identified"
    zf = frozenset(z)
    if not is_valid(zf, g0, treatment, outcome):
        return None, "z_invalid_at_g0"

    space = build_space_correct(cpdag)
    if g0 not in space.neighbours:
        return None, "g0_not_in_corrected_space"
    dists = distances_from(space, g0)

    return (
        Candidate(
            generator=generator,
            n=n,
            seed=seed,
            gen_params=gen_params,
            knows_fraction=knows_fraction,
            corruption_rate=corruption_rate,
            dag=dag,
            cpdag=cpdag,
            treatment=treatment,
            outcome=outcome,
            k_true=k_true,
            k_assumed=k_assumed,
            g0=g0,
            z=zf,
            space=space,
            dists=dists,
        ),
        "ok",
    )


def sweep_candidates(
    n: int, seed_base: int, n_seeds: int
) -> tuple[list[Candidate], dict[str, int]]:
    """Every gate-accepted candidate for one node count ``n``, deterministically.

    Returns ``(candidates, rejection_counts)``; the latter is kept so the
    search scope (how many were tried, how many rejected and why) can be
    reported honestly in ``SUMMARY.md``.
    """
    gens = dict(GEN_PARAM_GRID)
    if n >= 7:
        gens = dict(gens, decoupled_backdoor=[{"coupling": c} for c in (0.3, 0.6, 0.9)])

    out: list[Candidate] = []
    rejects: dict[str, int] = defaultdict(int)
    idx = 0
    for generator, param_list in gens.items():
        for params in param_list:
            for s in range(n_seeds):
                seed = seed_base + idx * 100003 + s
                for kf in KNOWS_FRACTIONS:
                    for cr in CORRUPTION_RATES:
                        cand, reason = try_build_candidate(generator, n, seed, params, kf, cr)
                        if cand is None:
                            rejects[reason] += 1
                        else:
                            out.append(cand)
            idx += 1
    return out, dict(rejects)


# --------------------------------------------------------------------------
# Serialization -- everything needed to replay a witness from scratch
# --------------------------------------------------------------------------


def mpdag_to_json(g: MPDAG) -> dict[str, Any]:
    """Serialise a graph so a witness can be replayed without rerunning the search."""
    return {
        "nodes": list(g.nodes),
        "directed": [list(e) for e in sorted(g.directed_edges)],
        "undirected": [list(e) for e in sorted(g.undirected_edges)],
        "edge_string": g.edge_string(),
    }


def mpdag_from_json(d: dict[str, Any]) -> MPDAG:
    """Inverse of :func:`mpdag_to_json`."""
    return MPDAG(
        d["nodes"],
        directed=[tuple(e) for e in d["directed"]],
        undirected=[tuple(e) for e in d["undirected"]],
    )


def sem_to_json(sem: LinearSEM) -> dict[str, Any]:
    """Serialise a SEM, so a replayed witness reproduces the same numbers exactly."""
    return {
        "dag": mpdag_to_json(sem.dag),
        "weights": {f"{t}|{h}": w for (t, h), w in sorted(sem.weights.items())},
        "noise_var": dict(sorted(sem.noise_var.items())),
    }


def sem_from_json(d: dict[str, Any]) -> LinearSEM:
    """Inverse of :func:`sem_to_json`."""
    dag = mpdag_from_json(d["dag"])
    weights: dict[Edge, float] = {}
    for key, w in d["weights"].items():
        t, h = key.split("|")
        weights[(t, h)] = float(w)
    return LinearSEM(
        dag=dag, weights=weights, noise_var={k: float(v) for k, v in d["noise_var"].items()}
    )


def candidate_instance_json(cand: Candidate) -> dict[str, Any]:
    """The instance-level fields shared by every witness drawn from ``cand``."""
    return {
        "generator": cand.generator,
        "n": cand.n,
        "seed": cand.seed,
        "gen_params": cand.gen_params,
        "knows_fraction": cand.knows_fraction,
        "corruption_rate": cand.corruption_rate,
        "true_dag": mpdag_to_json(cand.dag),
        "cpdag": mpdag_to_json(cand.cpdag),
        "g0": mpdag_to_json(cand.g0),
        "treatment": cand.treatment,
        "outcome": cand.outcome,
        "z": sorted(cand.z),
        "k_true": sorted(cand.k_true),
        "k_assumed": sorted(cand.k_assumed),
        "n_undirected_cpdag": cand.n_undirected,
        "n_knowledge_g0": len(cand.k0),
        "space_size": len(cand.space),
    }


# --------------------------------------------------------------------------
# C1 -- individual-state bias is not monotone in distance from G0
# --------------------------------------------------------------------------


@dataclass
class C1Witness:
    margin: float
    b_g: float
    b_h: float
    d_g: int
    d_h: int
    state_g: MPDAG
    state_h: MPDAG
    candidate: Candidate
    sem: LinearSEM

    def sort_key(self) -> tuple[int, int, float]:
        return (self.candidate.n, self.candidate.n_undirected, -self.margin)

    def to_json(self) -> dict[str, Any]:
        return {
            "instance": candidate_instance_json(self.candidate),
            "sem": sem_to_json(self.sem),
            "state_G": {**mpdag_to_json(self.state_g), "distance_from_G0": self.d_g, "B": self.b_g},
            "state_H": {**mpdag_to_json(self.state_h), "distance_from_G0": self.d_h, "B": self.b_h},
            "margin_B_G_minus_B_H": self.margin,
            "claim": (
                "d(G0,G) < d(G0,H) but B(H) < B(G) - margin: individual-state bias "
                "is not monotone in raw distance from G0 (Remark 4). Only the "
                "ball-maximum beta_up(d) is monotone, which is why B must be "
                "maximised over the shell (Theorem C) rather than read off one "
                "perturbation."
            ),
        }


def _bias_context_for(cand: Candidate, sem_seed_salt: int) -> tuple[LinearSEM, BiasContext] | None:
    rng = np.random.default_rng((cand.seed, sem_seed_salt))
    sem = random_sem(cand.dag, rng)
    try:
        ctx = make_context(sem, cand.cpdag, cand.treatment, cand.outcome, cand.z)
    except ValueError:
        return None
    return sem, ctx


def search_c1_on_candidate(cand: Candidate) -> C1Witness | None:
    """The best (largest-margin) C1 pair found within one candidate's space."""
    built = _bias_context_for(cand, sem_seed_salt=1)
    if built is None:
        return None
    sem, ctx = built

    b_values: dict[MPDAG, float] = {}
    for g, _d in cand.dists.items():
        try:
            b_values[g] = bias_at(ctx, g).worst
        except ValueError:
            continue

    by_dist: dict[int, list[MPDAG]] = defaultdict(list)
    for g, _val in b_values.items():
        by_dist[cand.dists[g]].append(g)

    dist_levels = sorted(by_dist)
    best: tuple[float, MPDAG, MPDAG, int, int] | None = None
    for i, d1 in enumerate(dist_levels):
        # Precompute the max B at d1 to prune: we want B(G) large at the smaller
        # distance and B(H) small at the larger one.
        g_candidates = by_dist[d1]
        for d2 in dist_levels[i + 1 :]:
            for g in g_candidates:
                bg = b_values[g]
                for h in by_dist[d2]:
                    margin = bg - b_values[h]
                    if margin > MARGIN_TOL and (best is None or margin > best[0]):
                        best = (margin, g, h, d1, d2)
    if best is None:
        return None
    margin, g, h, d1, d2 = best
    return C1Witness(
        margin=margin,
        b_g=b_values[g],
        b_h=b_values[h],
        d_g=d1,
        d_h=d2,
        state_g=g,
        state_h=h,
        candidate=cand,
        sem=sem,
    )


def run_c1_search(
    n_values: tuple[int, ...] = (4, 5, 6),
    n_seeds: int = 6,
    seed_base: int = 20260916,
) -> dict[str, Any]:
    """Search increasing ``n`` for the minimal C1 witness.

    Stops at the first ``n`` that yields a hit, so "minimal" is meaningful.
    """
    scope: dict[str, Any] = {"n_values_tried": [], "candidates_tried": 0, "rejection_counts": {}}
    found_at_n: list[C1Witness] = []
    chosen_n: int | None = None

    for n in n_values:
        t0 = time.perf_counter()
        candidates, rejects = sweep_candidates(n, seed_base=seed_base + n, n_seeds=n_seeds)
        scope["n_values_tried"].append(n)
        scope["candidates_tried"] += len(candidates)
        for k, v in rejects.items():
            scope["rejection_counts"][k] = scope["rejection_counts"].get(k, 0) + v

        hits: list[C1Witness] = []
        for cand in candidates:
            w = search_c1_on_candidate(cand)
            if w is not None:
                hits.append(w)
        print(
            f"[C1] n={n}: {len(candidates)} accepted instances, {len(hits)} with a "
            f"violation, {time.perf_counter() - t0:.1f}s"
        )
        if hits:
            found_at_n = hits
            chosen_n = n
            break

    if not found_at_n:
        return {
            "found": False,
            "scope": scope,
            "note": "No C1 violation found within the searched n/generator/seed grid.",
        }

    found_at_n.sort(key=lambda w: w.sort_key())
    best = found_at_n[0]
    return {
        "found": True,
        "chosen_n": chosen_n,
        "n_witnesses_at_chosen_n": len(found_at_n),
        "scope": scope,
        "witness": best.to_json(),
    }


def replay_c1(witness_json: dict[str, Any]) -> dict[str, float]:
    """Recompute B(G), B(H) and the margin from a persisted C1 witness."""
    inst = witness_json["instance"]
    cpdag = mpdag_from_json(inst["cpdag"])
    x, y = inst["treatment"], inst["outcome"]
    z = frozenset(inst["z"])
    sem = sem_from_json(witness_json["sem"])
    ctx = make_context(sem, cpdag, x, y, z)
    g = mpdag_from_json(witness_json["state_G"])
    h = mpdag_from_json(witness_json["state_H"])
    b_g = bias_at(ctx, g).worst
    b_h = bias_at(ctx, h).worst
    return {"B_G": b_g, "B_H": b_h, "margin": b_g - b_h}


# --------------------------------------------------------------------------
# C2 -- mean_abs_bias is not monotone under model inclusion
# --------------------------------------------------------------------------


def _mean_abs_bias(
    z: frozenset[Node], g: MPDAG, x: Node, y: Node, seed: tuple[int, ...], n_draws: int
) -> BiasStats:
    rng = np.random.default_rng(seed)
    return bias_stats(z, g, x, y, rng, n_draws)


def mc_standard_error(
    stats: BiasStats,
    g: MPDAG,
    z: frozenset[Node],
    x: Node,
    y: Node,
    seed: tuple[int, ...],
    n_draws: int,
) -> float:
    """Recompute the sample std of the per-draw |bias| and return the MC standard error of the mean.

    ``bias_stats`` itself only returns the mean, so this redraws the identical
    sequence (same seed) directly with the underlying primitives to get the
    per-draw values needed for a standard error -- not a second independent
    draw, the exact same one ``bias_stats`` consumed.
    """
    from bkrobust.core.oracle import extensions
    from bkrobust.demo.evaluate import adjusted_estimand
    from bkrobust.demo.evaluate import random_sem as _random_sem

    exts = extensions(g)
    rng = np.random.default_rng(seed)
    biases = []
    for d in exts:
        for _ in range(n_draws):
            sem = _random_sem(d, rng)
            tau = sem.true_total_effect(x, y)
            est = adjusted_estimand(sem, x, y, z)
            biases.append(abs(est - tau))
    arr = np.asarray(biases)
    if arr.size < 2:
        return float("nan")
    return float(arr.std(ddof=1) / np.sqrt(arr.size))


@dataclass
class RateSweepResult:
    n_instances: int
    n_pairs_total: int
    n_pairs_violating_mean: int
    n_pairs_violating_B: int  # noqa: N815 - matches the lemma labels used in THEOREMS.md
    per_instance: list[dict[str, Any]] = field(default_factory=list)


def enumerate_comparable_pairs(space: Space) -> list[tuple[MPDAG, MPDAG]]:
    """Every ordered pair ``(G, H)`` with ``[G]`` a PROPER subset of ``[H]``.

    That is, ``G`` is strictly below ``H`` in the model-inclusion order.
    """
    elems = list(space.elements)
    pairs = []
    for g in elems:
        rg = space.reps[g]
        for h in elems:
            if g is h:
                continue
            rh = space.reps[h]
            if rg < rh:
                pairs.append((g, h))
    return pairs


def run_c2_rate_sweep(
    n_values: tuple[int, ...] = (4, 5),
    n_seeds: int = 4,
    seed_base: int = 777001,
    n_draws: int = 60,
    max_space_size: int = 130,
    max_instances: int = 25,
) -> tuple[RateSweepResult, list[Candidate], dict[MPDAG, dict[str, float]]]:
    """Monotonicity-violation rate for ``mean_abs_bias`` vs ``B``, on the SAME instances.

    Kept to modest space sizes (``max_space_size``) and a moderate ``n_draws``
    (this is the RATE sweep -- the flagship pair and the consequence test below
    re-verify with far larger ``n_draws`` and independent repetitions) so the
    O(n_states^2) pair enumeration and the O(n_states * n_draws) bias sampling
    stay affordable across many instances.
    """
    instances: list[Candidate] = []
    for n in n_values:
        candidates, _ = sweep_candidates(n, seed_base=seed_base + n, n_seeds=n_seeds)
        for c in candidates:
            if len(c.space) <= max_space_size and len(c.k0) >= 2:
                instances.append(c)
        if len(instances) >= max_instances:
            break
    instances = instances[:max_instances]

    total_pairs = 0
    viol_mean = 0
    viol_b = 0
    per_instance: list[dict[str, Any]] = []
    cache: dict[MPDAG, dict[str, float]] = {}

    for idx, cand in enumerate(instances):
        built = _bias_context_for(cand, sem_seed_salt=2)
        if built is None:
            continue
        _, ctx = built

        pairs = enumerate_comparable_pairs(cand.space)
        if not pairs:
            continue

        mean_cache: dict[MPDAG, float] = {}
        b_cache: dict[MPDAG, float] = {}

        def mean_of(g: MPDAG, _cand: Candidate = cand, _idx: int = idx) -> float:
            hit = mean_cache.get(g)  # noqa: B023 - the closure is defined and consumed inside this iteration
            if hit is not None:
                return hit
            stats = _mean_abs_bias(
                _cand.z,
                g,
                _cand.treatment,
                _cand.outcome,
                (seed_base, _idx, _stable_hash(g.edge_string())),
                n_draws,
            )
            val = stats.mean_abs_bias if stats.n_evaluations > 0 else float("inf")
            mean_cache[g] = val  # noqa: B023 - the closure is defined and consumed inside this iteration
            return val

        def b_of(g: MPDAG) -> float:
            hit = b_cache.get(g)  # noqa: B023 - the closure is defined and consumed inside this iteration
            if hit is not None:
                return hit
            try:
                val = bias_at(ctx, g).worst  # noqa: B023 - the closure is defined and consumed inside this iteration
            except ValueError:
                val = float("nan")
            b_cache[g] = val  # noqa: B023 - the closure is defined and consumed inside this iteration
            return val

        inst_viol_mean = 0
        inst_viol_b = 0
        for g, h in pairs:
            mg, mh = mean_of(g), mean_of(h)
            if mg > mh + MARGIN_TOL:
                inst_viol_mean += 1
            bg, bh = b_of(g), b_of(h)
            if bg > bh + ZERO_TOL:
                inst_viol_b += 1
        total_pairs += len(pairs)
        viol_mean += inst_viol_mean
        viol_b += inst_viol_b
        per_instance.append(
            {
                **cand.tag(),
                "n_pairs": len(pairs),
                "n_violations_mean_abs_bias": inst_viol_mean,
                "n_violations_B": inst_viol_b,
            }
        )
        for g in set(mean_cache) | set(b_cache):
            cache[g] = {
                "mean_abs_bias": mean_cache.get(g, float("nan")),
                "B": b_cache.get(g, float("nan")),
            }

    result = RateSweepResult(
        n_instances=len(per_instance),
        n_pairs_total=total_pairs,
        n_pairs_violating_mean=viol_mean,
        n_pairs_violating_B=viol_b,
        per_instance=per_instance,
    )
    return result, instances, cache


@dataclass
class C2FlagshipWitness:
    candidate: Candidate
    state_g: MPDAG
    state_h: MPDAG
    n_draws: int
    seeds: list[int]
    means_g: list[float]
    means_h: list[float]
    se_g: list[float]
    se_h: list[float]

    def to_json(self) -> dict[str, Any]:
        margins = [mg - mh for mg, mh in zip(self.means_g, self.means_h, strict=True)]
        combined_se = [
            float(np.sqrt(seg**2 + seh**2)) for seg, seh in zip(self.se_g, self.se_h, strict=True)
        ]
        ratios = [
            m / se if se > 0 else float("inf") for m, se in zip(margins, combined_se, strict=True)
        ]
        n_reproduced = sum(1 for m in margins if m > 0)
        return {
            "instance": candidate_instance_json(self.candidate),
            "state_G": mpdag_to_json(self.state_g),
            "state_H": mpdag_to_json(self.state_h),
            "claim": (
                "G <= H (model inclusion, [G] proper subset of [H]) yet "
                "mean_abs_bias(G) > mean_abs_bias(H): the incumbent bias_stats "
                "mean is not monotone under the order, so {mean > eps} is not an "
                "up-set (docs/R_EPSILON_THEORY.md Remark to Theorem B / section 0)."
            ),
            "n_draws_per_repetition": self.n_draws,
            "mc_seeds": self.seeds,
            "mean_abs_bias_G_per_rep": self.means_g,
            "mean_abs_bias_H_per_rep": self.means_h,
            "mc_standard_error_G_per_rep": self.se_g,
            "mc_standard_error_H_per_rep": self.se_h,
            "margin_per_rep": margins,
            "margin_over_combined_SE_per_rep": ratios,
            "n_repetitions": len(margins),
            "n_repetitions_reproducing_sign": n_reproduced,
            "reproduction_rate": n_reproduced / len(margins) if margins else float("nan"),
            "min_margin_over_SE": min(ratios) if ratios else float("nan"),
        }


def verify_c2_flagship(
    cand: Candidate,
    g: MPDAG,
    h: MPDAG,
    n_draws: int = 400,
    n_reps: int = 8,
    seed_base: int = 990001,
) -> C2FlagshipWitness:
    """Re-verify one (G, H) pair with a large n_draws and several independent rngs."""
    seeds = [seed_base + 17 * r for r in range(n_reps)]
    means_g, means_h, se_g, se_h = [], [], [], []
    for s in seeds:
        stats_g = _mean_abs_bias(cand.z, g, cand.treatment, cand.outcome, (s, 0), n_draws)
        stats_h = _mean_abs_bias(cand.z, h, cand.treatment, cand.outcome, (s, 1), n_draws)
        means_g.append(stats_g.mean_abs_bias)
        means_h.append(stats_h.mean_abs_bias)
        se_g.append(
            mc_standard_error(stats_g, g, cand.z, cand.treatment, cand.outcome, (s, 0), n_draws)
        )
        se_h.append(
            mc_standard_error(stats_h, h, cand.z, cand.treatment, cand.outcome, (s, 1), n_draws)
        )
    return C2FlagshipWitness(
        candidate=cand,
        state_g=g,
        state_h=h,
        n_draws=n_draws,
        seeds=seeds,
        means_g=means_g,
        means_h=means_h,
        se_g=se_g,
        se_h=se_h,
    )


def find_best_c2_pair(
    cache: dict[MPDAG, dict[str, float]], instances: list[Candidate]
) -> tuple[Candidate, MPDAG, MPDAG, float] | None:
    """The largest-margin violating pair found anywhere in the rate sweep.

    Promoted to the flagship, where it is re-verified at a far larger draw count.
    """
    best: tuple[Candidate, MPDAG, MPDAG, float] | None = None
    for cand in instances:
        pairs = enumerate_comparable_pairs(cand.space)
        for g, h in pairs:
            dg, dh = cache.get(g), cache.get(h)
            if dg is None or dh is None:
                continue
            margin = dg["mean_abs_bias"] - dh["mean_abs_bias"]
            if margin > MARGIN_TOL and (best is None or margin > best[3]):
                best = (cand, g, h, margin)
    return best


# --- the consequence: retraction-only search vs brute-force BFS ------------


def retraction_only_radius(
    cand: Candidate,
    predicate: Any,
    *,
    cap: int = DEFAULT_SHELL_CAP,
) -> tuple[int, str | None]:
    """Mirrors ``epsilon.profile``'s incremental staircase, but for an arbitrary predicate.

    This is exactly the algorithm licensed for ``B`` by Theorem U/Corollary R1
    (search only the up-set ``up(G0)``, generated by dropping subsets of
    ``K_{G0}``, running max as depth increases). It is legitimate for ``B``
    because ``{B > eps}`` is an up-set. Applied to a predicate that is NOT
    upward-closed -- ``mean_abs_bias`` -- nothing in the theory licenses it;
    it is built here purely so it can be compared against a real brute-force
    search over the whole enumerated space with the identical predicate.

    Returns:
        ``(radius, witness_edge_string)``, using
        ``bkrobust.core.conventions.UNREACHED`` (-1) if no shell up to
        ``|K_{G0}|`` satisfies the predicate.
    """
    from bkrobust.core.conventions import UNREACHED

    k0 = space_knowledge_of(cand.cpdag, cand.g0)
    for d in range(0, len(k0) + 1):
        try:
            states = retraction_shell(cand.cpdag, cand.g0, d, cap=cap)
        except MemoryError:
            return UNREACHED, None
        for s in states:
            if predicate(s):
                return d, s.edge_string()
    return UNREACHED, None


def brute_force_radius(cand: Candidate, predicate: Any) -> tuple[int, str | None]:
    """The already-correct full-space BFS radius -- what ``_radius_eps`` computes today."""
    r, witness = radius(cand.space, cand.dists, predicate)
    return r, (witness.edge_string() if witness is not None else None)


@dataclass
class ConsequenceWitness:
    candidate: Candidate
    eps: float
    n_draws: int
    seed: tuple[int, ...]
    retraction_radius: int
    retraction_witness: str | None
    brute_force_radius: int
    brute_force_witness: str | None

    def to_json(self) -> dict[str, Any]:
        return {
            "instance": candidate_instance_json(self.candidate),
            "eps": self.eps,
            "n_draws": self.n_draws,
            "mc_seed": list(self.seed),
            "retraction_only_radius": self.retraction_radius,
            "retraction_only_witness": self.retraction_witness,
            "brute_force_radius": self.brute_force_radius,
            "brute_force_witness": self.brute_force_witness,
            "agree": self.retraction_radius == self.brute_force_radius,
            "claim": (
                "The retraction-only search mirrors the fast machinery that is "
                "exact for B (Theorem U/Corollary R1) but is NOT licensed for "
                "mean_abs_bias, since {mean_abs_bias > eps} is not an up-set. "
                "Comparing it against a genuine brute-force BFS over the whole "
                "enumerated space with the same predicate is the check of "
                "whether that lack of a license is merely theoretical or "
                "actually produces a wrong number."
            ),
        }


def search_c2_consequence(
    instances: list[Candidate], cache: dict[MPDAG, dict[str, float]], n_draws: int = 200
) -> dict[str, Any]:
    """Look for an instance/eps where retraction-only and brute-force radii disagree."""
    tried = 0
    mismatches: list[ConsequenceWitness] = []
    for cand in instances:
        mean_cache: dict[MPDAG, float] = {}

        def mean_of(g: MPDAG, _cand: Candidate = cand) -> float:
            hit = mean_cache.get(g)  # noqa: B023 - the closure is defined and consumed inside this iteration
            if hit is not None:
                return hit
            seed = (13, _cand.seed, _stable_hash(g.edge_string()))
            stats = _mean_abs_bias(_cand.z, g, _cand.treatment, _cand.outcome, seed, n_draws)
            val = stats.mean_abs_bias if stats.n_evaluations > 0 else float("inf")
            mean_cache[g] = val  # noqa: B023 - the closure is defined and consumed inside this iteration
            return val

        values = sorted({mean_of(g) for g in cand.dists})
        if len(values) < 2:
            continue
        # Probe eps at the midpoints between consecutive distinct observed
        # values: these are exactly the thresholds where the predicate's
        # truth value can change, so no informative eps is skipped.
        eps_grid = [(a + b) / 2.0 for a, b in itertools.pairwise(values)]
        for eps in eps_grid:
            tried += 1
            predicate = lambda g, _e=eps: mean_of(g) > _e  # noqa: E731
            r_retract, w_retract = retraction_only_radius(cand, predicate)
            r_brute, w_brute = brute_force_radius(cand, predicate)
            if r_retract != r_brute:
                mismatches.append(
                    ConsequenceWitness(
                        candidate=cand,
                        eps=eps,
                        n_draws=n_draws,
                        seed=(13, cand.seed, -1),
                        retraction_radius=r_retract,
                        retraction_witness=w_retract,
                        brute_force_radius=r_brute,
                        brute_force_witness=w_brute,
                    )
                )
    mismatches.sort(key=lambda w: (w.candidate.n, w.candidate.n_undirected))
    return {
        "n_instances_searched": len(instances),
        "n_eps_probed": tried,
        "n_mismatches": len(mismatches),
        "mismatches": [m.to_json() for m in mismatches[:5]],
    }


def replay_c2_flagship(witness_json: dict[str, Any]) -> dict[str, Any]:
    """Recompute the flagship C2 numbers from a persisted witness."""
    inst = witness_json["instance"]
    z = frozenset(inst["z"])
    x, y = inst["treatment"], inst["outcome"]
    g = mpdag_from_json(witness_json["state_G"])
    h = mpdag_from_json(witness_json["state_H"])
    n_draws = witness_json["n_draws_per_repetition"]
    seeds = witness_json["mc_seeds"]
    means_g, means_h = [], []
    for s in seeds:
        sg = _mean_abs_bias(z, g, x, y, (s, 0), n_draws)
        sh = _mean_abs_bias(z, h, x, y, (s, 1), n_draws)
        means_g.append(sg.mean_abs_bias)
        means_h.append(sh.mean_abs_bias)
    n_reproduced = sum(1 for a, b in zip(means_g, means_h, strict=True) if a > b)
    return {
        "mean_abs_bias_G_per_rep": means_g,
        "mean_abs_bias_H_per_rep": means_h,
        "n_repetitions_reproducing_sign": n_reproduced,
        "reproduction_rate": n_reproduced / len(seeds) if seeds else float("nan"),
    }


__all__ = [
    "Candidate",
    "find_best_c2_pair",
    "mpdag_from_json",
    "mpdag_to_json",
    "replay_c1",
    "replay_c2_flagship",
    "run_c1_search",
    "run_c2_rate_sweep",
    "search_c2_consequence",
    "sem_from_json",
    "sem_to_json",
    "sweep_candidates",
    "try_build_candidate",
    "verify_c2_flagship",
]


# --------------------------------------------------------------------------
# Driver: run everything, persist JSON witnesses, write SUMMARY.md
# --------------------------------------------------------------------------


def _write_json(path: Path, obj: Any) -> None:
    path.write_text(json.dumps(obj, indent=2, sort_keys=False, default=str))


def main(out_dir: str | Path = "results/epsilon/counterexamples") -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    t_start = time.perf_counter()

    print("=== C1: individual-state bias is not monotone in distance ===")
    c1 = run_c1_search(n_values=(4, 5, 6), n_seeds=100)
    _write_json(out / "c1_witness.json", c1)
    if c1["found"]:
        wj = c1["witness"]
        print(
            f"C1: n={c1['chosen_n']}, k={wj['instance']['n_undirected_cpdag']}, "
            f"margin={wj['margin_B_G_minus_B_H']:.6g}, "
            f"d(G)={wj['state_G']['distance_from_G0']} B(G)={wj['state_G']['B']:.6g}, "
            f"d(H)={wj['state_H']['distance_from_G0']} B(H)={wj['state_H']['B']:.6g}"
        )
        replayed = replay_c1(wj)
        assert abs(replayed["margin"] - wj["margin_B_G_minus_B_H"]) < 1e-9, "C1 replay mismatch"
        print(f"C1 replay check: margin={replayed['margin']:.6g} (matches stored)")
    else:
        print("C1: NO violation found within the searched scope.", c1["scope"])

    print()
    print("=== C2: mean_abs_bias is not monotone under model inclusion ===")
    rate, instances, cache = run_c2_rate_sweep(
        n_values=(4, 5, 6), n_seeds=6, n_draws=80, max_space_size=150, max_instances=40
    )
    rate_json = {
        "n_instances": rate.n_instances,
        "n_pairs_total": rate.n_pairs_total,
        "n_pairs_violating_mean_abs_bias": rate.n_pairs_violating_mean,
        "n_pairs_violating_B": rate.n_pairs_violating_B,
        "violation_rate_mean_abs_bias": (
            rate.n_pairs_violating_mean / rate.n_pairs_total if rate.n_pairs_total else float("nan")
        ),
        "violation_rate_B": (
            rate.n_pairs_violating_B / rate.n_pairs_total if rate.n_pairs_total else float("nan")
        ),
        "n_draws_per_state": 80,
        "per_instance": rate.per_instance,
    }
    _write_json(out / "c2_rate_sweep.json", rate_json)
    print(
        f"C2 rate sweep: {rate.n_instances} instances, "
        f"{rate.n_pairs_total} ordered comparable pairs; "
        f"mean_abs_bias violates in {rate.n_pairs_violating_mean} "
        f"({rate_json['violation_rate_mean_abs_bias']:.1%}), "
        f"B violates in {rate.n_pairs_violating_B} ({rate_json['violation_rate_B']:.1%})"
    )

    best = find_best_c2_pair(cache, instances)
    if best is None:
        print("C2 flagship: no violating pair found in the rate-sweep pool (unexpected).")
        flagship_json: dict[str, Any] = {"found": False}
    else:
        cand, g, h, _margin_est = best
        flagship = verify_c2_flagship(cand, g, h, n_draws=400, n_reps=8)
        flagship_json = {"found": True, **flagship.to_json()}
        print(
            f"C2 flagship: n={cand.n}, k={cand.n_undirected}, "
            f"margin~{np.mean(flagship.means_g) - np.mean(flagship.means_h):.4g} "
            f"(mean over {len(flagship.seeds)} reps), "
            f"min margin/SE={flagship_json['min_margin_over_SE']:.2f}, "
            f"reproduced in {flagship_json['n_repetitions_reproducing_sign']}"
            f"/{flagship_json['n_repetitions']} reps"
        )
        replayed = replay_c2_flagship(flagship_json)
        assert (
            replayed["n_repetitions_reproducing_sign"]
            == flagship_json["n_repetitions_reproducing_sign"]
        ), "C2 flagship replay mismatch"
        print("C2 flagship replay check: reproduction count matches stored.")
    _write_json(out / "c2_flagship.json", flagship_json)

    print()
    print("=== C2 consequence: retraction-only search vs brute-force BFS ===")
    consequence = search_c2_consequence(instances, cache, n_draws=200)
    _write_json(out / "c2_consequence.json", consequence)
    print(
        f"Consequence search: {consequence['n_instances_searched']} instances, "
        f"{consequence['n_eps_probed']} eps thresholds probed, "
        f"{consequence['n_mismatches']} mismatches found."
    )
    if consequence["n_mismatches"] == 0:
        print(
            "No instance found where retraction-only search disagrees with brute-force "
            "BFS for {mean_abs_bias > eps} within this scope -- reported honestly, not "
            "papered over."
        )
    else:
        m0 = consequence["mismatches"][0]
        print(
            f"Example mismatch: eps={m0['eps']:.6g}, retraction_only_radius="
            f"{m0['retraction_only_radius']}, brute_force_radius={m0['brute_force_radius']}"
        )

    manifest = {
        "seconds_total": time.perf_counter() - t_start,
        "margin_tol": MARGIN_TOL,
        "zero_tol": ZERO_TOL,
        "max_undirected_edges": MAX_UNDIRECTED_EDGES,
        "files": [
            "c1_witness.json",
            "c2_rate_sweep.json",
            "c2_flagship.json",
            "c2_consequence.json",
        ],
    }
    _write_json(out / "manifest.json", manifest)
    print(f"\nDone in {manifest['seconds_total']:.1f}s. Output: {out}")


if __name__ == "__main__":
    main()
