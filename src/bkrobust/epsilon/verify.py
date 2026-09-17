"""The falsification harness for ``docs/R_EPSILON_THEORY.md`` section 7.

Seven falsifiable predictions, checked exhaustively over small, fully
enumerated instances:

* **P1** zero bias inside the certified shell (Theorem A).
* **P2** monotonicity of ``B`` under model inclusion (Theorem B).
* **P3** the staircase: sphere sufficiency, direct-shell evaluation, and
  monotonicity of the raw shell maximum itself (Theorem C.1/C.3/C.4).
* **P4** the three ``r_eps`` strategies agree with each other and with an
  independent brute-force BFS (Theorem C.4/C.5).
* **P5** the nesting ``r_val <= r_0 <= r_eps`` and monotonicity in ``eps``
  (Theorem C.6).
* **P6** the practitioner's certificate: the realised error is bounded by
  ``beta_up(k)`` and by ``B`` of the specific state that contains the truth
  (Theorem D). ``k`` counts *false orientations in the closure* ``K_{G0}``,
  not asserted claims -- Meek's rules can cascade a small number of false
  claims into a larger false closure, and the certificate's ``k`` is
  denominated in the latter.
* **P7** the semi-local possible-parent-set evaluator (Theorem E) agrees with
  brute-force DAG enumeration on every element of the enumerated space.

These are theorems, not conjectures, so a violation here is a bug -- either
in this module or in :mod:`bkrobust.epsilon` -- and is persisted with a full
replayable witness rather than smoothed into an aggregate. Tolerance is
``1e-9`` absolute throughout, matching :data:`bkrobust.epsilon.bias.ZERO_TOL`,
and is never widened to make a check pass.

**Space construction.** Every enumerated space here is built with
:func:`bkrobust.search.space_fixed.build_space_correct`, never with
:func:`bkrobust.core.spacelib.build_space`: the latter calls
:func:`bkrobust.demo.space.enumerate_space`, whose chordality-based validity
filter is known to omit legitimate knowledge states (see the module
docstring of :mod:`bkrobust.search.space_fixed`). Using the buggy enumerator
here would silently under-count the space P2/P3/P4/P7 are checked against.

**Instance generation** mirrors :func:`bkrobust.synth.runner.run_instance`
(the degeneracy gate, the treatment/outcome pick, the knowledge draw and
corruption) but adds a tighter cap on the CPDAG's undirected-edge count --
7, not the runner's 9 -- so that ``3 ** k`` stays a few thousand at most and
the corrected enumerator is cheap, and corrupts with ``flip``/``omit`` (never
just ``none``) often enough that a good fraction of analysts hold knowledge
that is false of the true DAG, which is what P6 needs to exercise.

**Determinism.** Every instance and every SEM draw gets its own
``numpy.random.default_rng``, seeded from a tuple of ``(base_seed, index,
...)`` -- never the global RNG, and never a stream shared across instances or
draws.

Run as a script to execute the sweep and write
``results/epsilon/verification/{manifest,violations}.json``::

    python -m bkrobust.epsilon.verify --target 400 --sem-draws 3
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from bkrobust.core.conventions import RADIUS_CONVENTION, UNREACHED
from bkrobust.core.oracle import is_valid
from bkrobust.core.resultsio import git_sha, is_dirty
from bkrobust.core.spacelib import Space, distances_from, radius
from bkrobust.demo.evaluate import LinearSEM, optimal_adjustment_set_mpdag, random_sem
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.epsilon.bias import (
    ZERO_TOL,
    BiasContext,
    bias_at,
    make_context,
    possible_parent_sets_enumerated,
    possible_parent_sets_semilocal,
)
from bkrobust.epsilon.profile import ShellBias, bias_profile, r_epsilon
from bkrobust.search.space_fixed import build_space_correct, knowledge_of
from bkrobust.synth.generators import GENERATORS
from bkrobust.synth.knowledge import check_consistent, compound, draw_k_true, flip, omit

# Imported despite the leading underscore: the task this module implements
# calls for reusing exactly the runner's own screening logic (gate plus
# treatment/outcome selection) rather than forking a second copy that could
# drift from it. See ``bkrobust.synth.runner.run_instance``.
from bkrobust.synth.runner import _pick_treatment_outcome, gate

Edge = tuple[str, str]
Node = str

#: Absolute tolerance for every floating-point comparison in this module.
#: Matches :data:`bkrobust.epsilon.bias.ZERO_TOL`. Never widened to pass a
#: check -- a failure at this tolerance is reported, with its magnitude, as
#: either float noise or a real violation, never silently absorbed.
TOL: float = ZERO_TOL

#: Refuse instances whose CPDAG has more undirected edges than this, so that
#: ``3 ** k`` (the corrected enumerator's cost) stays in the low thousands.
MAX_UNDIRECTED_FOR_SPACE: int = 7

#: The seven falsifiable predictions, in the order ``docs/R_EPSILON_THEORY.md``
#: section 7 states them.
PREDICATES: tuple[str, ...] = ("P1", "P2", "P3", "P4", "P5", "P6", "P7")

_GENERATOR_CYCLE: tuple[str, ...] = (
    "erdos_renyi",
    "scale_free",
    "block",
    "decoupled_backdoor",
)
_N_CYCLE: tuple[int, ...] = (5, 6, 7, 8)
_ERDOS_RENYI_PROBS: tuple[float, ...] = (0.25, 0.35, 0.45)
_SCALE_FREE_M: tuple[int, ...] = (1, 2)
_BLOCK_PARAMS: tuple[dict[str, float], ...] = (
    {"n_blocks": 2, "p_within": 0.5, "p_between": 0.1},
    {"n_blocks": 2, "p_within": 0.6, "p_between": 0.2},
)
_COUPLINGS: tuple[float, ...] = (0.0, 0.5, 1.0)

#: ``(knows_fraction, corruption_name, corruption_rate)``, cycled by instance
#: index. Weighted towards ``flip`` (which introduces orientations false of
#: the truth) over ``omit`` (which only withdraws claims, never falsifies
#: one) because P6's error budget ``k`` is only exercised by false claims.
_CORRUPTION_PROFILES: tuple[tuple[float, str, float], ...] = (
    (1.0, "none", 0.0),
    (1.0, "flip", 0.25),
    (1.0, "flip", 0.5),
    (0.8, "flip", 0.4),
    (1.0, "omit", 0.3),
    (0.9, "compound", 0.3),
    (1.0, "flip", 0.75),
)


class RejectedError(Exception):
    """Raised by :func:`_generate_graph_instance` for every rejection reason.

    Mirrors :class:`bkrobust.synth.runner.GateRejectedError`, kept separate
    rather than reused so this module's rejection vocabulary (which adds
    ``"cpdag_too_large_for_space"`` for the tighter cap here) is self
    contained.
    """

    def __init__(self, reason: str) -> None:
        """Store the machine-readable rejection code.

        Args:
            reason: Short rejection code, tallied by
                :func:`run_verification`.
        """
        super().__init__(reason)
        self.reason = reason


@dataclass(frozen=True)
class GraphInstance:
    """One accepted, purely structural instance: a graph, a query, a space.

    Everything here is independent of the SEM draw (``Sigma``), which is why
    it is built once per accepted graph and then reused across
    :data:`_CORRUPTION_PROFILES`'s several SEM draws.

    Attributes:
        instance_id: ``"g{index:06d}"``.
        generator: Key into :data:`bkrobust.synth.generators.GENERATORS`.
        n: Node count.
        seed: The attempt index this instance was drawn at (its RNG seed).
        dag: The ground-truth DAG.
        cpdag: Its CPDAG.
        x: Treatment.
        y: Outcome.
        g0: The analyst's knowledge state.
        z: The adjustment set, fixed once from ``g0``.
        k_true: The analyst's true beliefs, before corruption.
        k_assumed: What the analyst actually asserts, after corruption.
        corruption: The corruption applied (``"none"``, ``"flip"``,
            ``"omit"`` or ``"compound"``).
        corruption_rate: The rate passed to that corruption.
        space: The corrected, fully enumerated space of ``cpdag``.
        dists: BFS distance from ``g0`` to every element of ``space``.
        k_g0: ``K_{G0}`` -- the orientations ``g0`` adds to ``cpdag``.
        false_orientations: The subset of ``k_g0`` that disagrees with
            ``dag`` -- what P6's error budget ``k`` counts.
        r_val_bf: The validity radius, computed by brute force directly
            against ``space`` (never via the accelerated hybrid, to keep
            this a true differential test).
    """

    instance_id: str
    generator: str
    n: int
    seed: int
    dag: MPDAG
    cpdag: MPDAG
    x: Node
    y: Node
    g0: MPDAG
    z: frozenset[Node]
    k_true: tuple[Edge, ...]
    k_assumed: tuple[Edge, ...]
    corruption: str
    corruption_rate: float
    space: Space
    dists: dict[MPDAG, int]
    k_g0: tuple[Edge, ...]
    false_orientations: tuple[Edge, ...]
    r_val_bf: int


@dataclass(frozen=True)
class ContextBundle:
    """One SEM draw's :class:`BiasContext`, kept with enough to build a witness.

    Attributes:
        context_id: ``"{instance_id}_sem{draw}"``.
        ctx: The bias context for this draw.
        sem: The linear SEM the covariance was computed from -- kept so a
            violation's witness can carry the exact weights and noise
            variances, not just the derived ``Sigma``.
        sem_seed: The seed tuple the SEM's RNG was built from.
    """

    context_id: str
    ctx: BiasContext
    sem: LinearSEM
    sem_seed: tuple[int, ...]


@dataclass
class PredicateTally:
    """How many checks a prediction underwent, and how many failed.

    Attributes:
        checked: Number of individual assertions evaluated.
        violations: Number that failed at :data:`TOL`.
    """

    checked: int = 0
    violations: int = 0


def _generator_params(name: str, n: int, cycle_index: int) -> dict[str, Any]:
    """Deterministic keyword arguments for one of :data:`GENERATORS`.

    Args:
        name: Generator key.
        n: Node count, needed to keep ``scale_free``'s ``m_attach < n``.
        cycle_index: Rotates through each generator's own parameter grid.

    Returns:
        Keyword arguments for ``GENERATORS[name]``.

    Raises:
        ValueError: On an unknown generator name.
    """
    if name == "erdos_renyi":
        return {"edge_prob": _ERDOS_RENYI_PROBS[cycle_index % len(_ERDOS_RENYI_PROBS)]}
    if name == "scale_free":
        m = _SCALE_FREE_M[cycle_index % len(_SCALE_FREE_M)]
        return {"m_attach": max(1, min(m, n - 1))}
    if name == "block":
        return dict(_BLOCK_PARAMS[cycle_index % len(_BLOCK_PARAMS)])
    if name == "decoupled_backdoor":
        return {"coupling": _COUPLINGS[cycle_index % len(_COUPLINGS)]}
    raise ValueError(f"unknown generator {name!r}")


def _generate_graph_instance(index: int, base_seed: int) -> GraphInstance:
    """Draw and screen one instance, following ``run_instance``'s pipeline.

    Args:
        index: Attempt index; also the sole seed material for this
            instance's RNG, so re-running with the same ``base_seed`` and
            ``index`` reproduces the same instance exactly.
        base_seed: The sweep's root seed.

    Returns:
        A :class:`GraphInstance`.

    Raises:
        RejectedError: On any of the rejection reasons ``run_verification``
            tallies: ``"cpdag_too_large_for_space"``, the five reasons of
            :func:`bkrobust.synth.runner.gate`, ``"k_assumed_inconsistent"``,
            ``"z_not_identified"``, ``"z_invalid_at_g0"``, or
            ``"g0_not_in_space"``.
    """
    rng = np.random.default_rng([base_seed, index, 0xA])
    generator_name = _GENERATOR_CYCLE[index % len(_GENERATOR_CYCLE)]
    n = _N_CYCLE[index % len(_N_CYCLE)]
    if generator_name == "decoupled_backdoor" and n < 7:
        n = 7 if (index // len(_GENERATOR_CYCLE)) % 2 == 0 else 8
    params = _generator_params(generator_name, n, index // len(_GENERATOR_CYCLE))

    dag = GENERATORS[generator_name](n, rng, **params)
    cpdag = dag_to_cpdag(dag)
    treatment, outcome = _pick_treatment_outcome(generator_name, dag, rng)

    if len(cpdag.undirected_edges) > MAX_UNDIRECTED_FOR_SPACE:
        raise RejectedError("cpdag_too_large_for_space")

    accepted, reason = gate(dag, cpdag, treatment, outcome)
    if not accepted:
        raise RejectedError(reason)

    knows_fraction, corruption_name, corruption_rate = _CORRUPTION_PROFILES[
        index % len(_CORRUPTION_PROFILES)
    ]
    k_true = draw_k_true(dag, cpdag, rng, knows_fraction)
    if corruption_name == "none":
        k_assumed = list(k_true)
    elif corruption_name == "flip":
        k_assumed = flip(k_true, rng, corruption_rate)
    elif corruption_name == "omit":
        k_assumed = omit(k_true, rng, corruption_rate)
    else:
        k_assumed = compound(k_true, rng, corruption_rate)

    if not check_consistent(cpdag, k_assumed):
        raise RejectedError("k_assumed_inconsistent")

    g0 = apply_orientations(cpdag, k_assumed)
    assert g0 is not None  # guaranteed by check_consistent above

    z = optimal_adjustment_set_mpdag(g0, treatment, outcome)
    if z is None:
        raise RejectedError("z_not_identified")
    zf = frozenset(z)
    if not is_valid(zf, g0, treatment, outcome):
        raise RejectedError("z_invalid_at_g0")

    space = build_space_correct(cpdag)
    if g0 not in space.neighbours:
        raise RejectedError("g0_not_in_space")
    dists = distances_from(space, g0)

    k_g0 = tuple(knowledge_of(cpdag, g0))
    true_dir = dag.directed_edges
    false_orientations = tuple(e for e in k_g0 if e not in true_dir)

    r_val_bf, _ = radius(space, dists, lambda g: not is_valid(zf, g, treatment, outcome))

    return GraphInstance(
        instance_id=f"g{index:06d}",
        generator=generator_name,
        n=n,
        seed=index,
        dag=dag,
        cpdag=cpdag,
        x=treatment,
        y=outcome,
        g0=g0,
        z=zf,
        k_true=tuple(sorted(k_true)),
        k_assumed=tuple(sorted(k_assumed)),
        corruption=corruption_name,
        corruption_rate=float(corruption_rate),
        space=space,
        dists=dists,
        k_g0=k_g0,
        false_orientations=false_orientations,
        r_val_bf=r_val_bf,
    )


def _powerset(items: Sequence[Edge]) -> Iterator[tuple[Edge, ...]]:
    """Every subset of ``items``, smallest first.

    Args:
        items: The edges to take subsets of.

    Yields:
        Each subset as a tuple, in :func:`itertools.combinations` order.
    """
    for size in range(len(items) + 1):
        yield from itertools.combinations(items, size)


def _radius_leq(a: int, b: int) -> bool:
    """``a <= b`` under the radius convention, where ``UNREACHED`` is a status.

    ``UNREACHED`` (Theorem B2's failure set being empty) is treated as
    "no perturbation reaches this failure", i.e. as larger than every finite
    radius -- consistent with a smaller ``eps`` only ever making the failure
    set bigger (Corollary B2), never the reverse.

    Args:
        a: A radius or :data:`UNREACHED`.
        b: A radius or :data:`UNREACHED`.

    Returns:
        Whether ``a`` precedes or equals ``b`` in that ordering.
    """
    if a == UNREACHED and b == UNREACHED:
        return True
    if b == UNREACHED:
        return True
    if a == UNREACHED:
        return False
    return a <= b


def _eps_grid(values: Iterable[float]) -> list[float]:
    """A small, deterministic set of thresholds that straddle ``values``.

    Built to land exactly at, and exactly between, the distinct ``B`` values
    observed at one context -- which is what forces the crossing tested by
    P4/P5 to fall at a specific, checkable depth rather than an arbitrary one.

    Args:
        values: The ``B(G)`` values observed across the enumerated space.

    Returns:
        A sorted list of at most 10 thresholds, always including ``0.0``.
    """
    distinct = sorted({round(float(v), 12) for v in values})
    candidates: set[float] = {0.0}
    if distinct:
        candidates.add(distinct[0])
        candidates.add(distinct[-1])
        candidates.add(distinct[-1] * 2.0 + 1e-6)
        candidates.add(distinct[len(distinct) // 2])
        for a, b in itertools.pairwise(distinct):
            if b > a:
                candidates.add((a + b) / 2.0)
    grid = sorted(candidates)
    if len(grid) > 10:
        positions = np.linspace(0, len(grid) - 1, 10)
        grid = sorted({grid[round(p)] for p in positions})
    return grid


def _witness(
    predicate: str,
    graph: GraphInstance,
    bundle: ContextBundle | None,
    *,
    g: MPDAG | None = None,
    h: MPDAG | None = None,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a fully replayable witness dict for a violation.

    Every field needed to reconstruct the check independently of this run is
    included: the CPDAG, ``G0``, the query, the two states involved (when
    applicable), and -- when a SEM draw is implicated -- its exact weights
    and noise variances, not just the derived covariance.

    Args:
        predicate: Which of P1-P7 this witnesses.
        graph: The structural instance.
        bundle: The SEM context, if the violation is context-dependent
            (``None`` for the purely structural P7).
        g: The first state implicated, if any.
        h: The second state implicated, if any (e.g. P2's ``H``).
        extra: Predicate-specific fields (values, depths, thresholds).

    Returns:
        A JSON-serialisable dict.
    """
    out: dict[str, Any] = {
        "predicate": predicate,
        "instance_id": graph.instance_id,
        "generator": graph.generator,
        "n": graph.n,
        "cpdag": graph.cpdag.edge_string(),
        "g0": graph.g0.edge_string(),
        "true_dag": graph.dag.edge_string(),
        "x": graph.x,
        "y": graph.y,
        "z": sorted(graph.z),
        "k_g0": [list(e) for e in graph.k_g0],
        "false_orientations": [list(e) for e in graph.false_orientations],
        "corruption": graph.corruption,
        "corruption_rate": graph.corruption_rate,
        "r_val_bf": graph.r_val_bf,
    }
    if bundle is not None:
        out["context_id"] = bundle.context_id
        out["theta_z"] = bundle.ctx.theta_z
        out["tau_true"] = bundle.ctx.tau_true
        out["sem_seed"] = list(bundle.sem_seed)
        out["sem_weights"] = [
            [tail, head, w] for (tail, head), w in sorted(bundle.sem.weights.items())
        ]
        out["sem_noise_var"] = dict(sorted(bundle.sem.noise_var.items()))
    if g is not None:
        out["g"] = g.edge_string()
    if h is not None:
        out["h"] = h.edge_string()
    if extra:
        out.update(extra)
    return out


# --- the seven checks -------------------------------------------------------


def _check_p1(
    graph: GraphInstance,
    bundle: ContextBundle,
    b_cache: dict[MPDAG, float],
    tallies: dict[str, PredicateTally],
    violations: list[dict[str, Any]],
) -> None:
    """P1: ``B(G) = 0`` for every ``G`` with ``d(G0,G) < r_val`` (Theorem A)."""
    tally = tallies["P1"]
    for g, d in graph.dists.items():
        if d >= graph.r_val_bf:
            continue
        tally.checked += 1
        if b_cache[g] > TOL:
            tally.violations += 1
            violations.append(
                _witness(
                    "P1",
                    graph,
                    bundle,
                    g=g,
                    extra={"d": d, "B_g": b_cache[g]},
                )
            )


def _check_p2(
    graph: GraphInstance,
    bundle: ContextBundle,
    b_cache: dict[MPDAG, float],
    tallies: dict[str, PredicateTally],
    violations: list[dict[str, Any]],
) -> None:
    """P2: ``B(G) <= B(H)`` for every ordered pair ``G <= H`` (Theorem B)."""
    tally = tallies["P2"]
    elements = graph.space.elements
    reps = graph.space.reps
    for g in elements:
        rg = reps[g]
        bg = b_cache[g]
        for h in elements:
            if not (rg <= reps[h]):
                continue
            tally.checked += 1
            bh = b_cache[h]
            if bg > bh + TOL:
                tally.violations += 1
                violations.append(
                    _witness("P2", graph, bundle, g=g, h=h, extra={"B_g": bg, "B_h": bh})
                )


def _check_p3(
    graph: GraphInstance,
    bundle: ContextBundle,
    profile: Sequence[ShellBias],
    b_cache: dict[MPDAG, float],
    tallies: dict[str, PredicateTally],
    violations: list[dict[str, Any]],
) -> None:
    """P3: the staircase (Theorem C.1, C.3, C.4).

    Three things are checked at every depth ``d`` in ``profile``, each a
    distinct sub-claim rather than one aggregate:

    1. brute-force sphere-max equals brute-force ball-max (Theorem C.3,
       sphere sufficiency -- no library code under test, pure definitions);
    2. ``retraction_shell``'s ``shell_max`` equals the brute-force sphere-max
       (Theorem C.4, direct shell evaluation -- this is the differential
       test on :func:`bkrobust.epsilon.profile.shell_bias`);
    3. the raw ``shell_max`` sequence itself is non-decreasing across
       consecutive depths (the stronger form Theorem C.3 licenses: the ball
       maximum is *already* attained shell by shell, before any running max
       is taken).
    """
    tally = tallies["P3"]
    reps = graph.space.reps
    g0_rep = reps[graph.g0]
    up_g0 = [g for g in graph.space.elements if g0_rep <= reps[g]]

    for step in profile:
        d = step.d
        sphere_vals = [b_cache[g] for g in up_g0 if graph.dists.get(g) == d]
        ball_vals = [
            b_cache[g] for g in up_g0 if graph.dists.get(g) is not None and graph.dists[g] <= d
        ]
        sphere_bf = max(sphere_vals) if sphere_vals else 0.0
        ball_bf = max(ball_vals) if ball_vals else 0.0

        tally.checked += 1
        if abs(sphere_bf - ball_bf) > TOL:
            tally.violations += 1
            violations.append(
                _witness(
                    "P3",
                    graph,
                    bundle,
                    extra={
                        "check": "sphere_eq_ball_bruteforce",
                        "d": d,
                        "sphere_bf": sphere_bf,
                        "ball_bf": ball_bf,
                    },
                )
            )

        tally.checked += 1
        if abs(step.shell_max - sphere_bf) > TOL:
            tally.violations += 1
            violations.append(
                _witness(
                    "P3",
                    graph,
                    bundle,
                    extra={
                        "check": "shell_max_eq_sphere_bruteforce",
                        "d": d,
                        "shell_max": step.shell_max,
                        "sphere_bf": sphere_bf,
                        "shell_witness": step.witness,
                    },
                )
            )

        tally.checked += 1
        if abs(step.beta_up - ball_bf) > TOL:
            tally.violations += 1
            violations.append(
                _witness(
                    "P3",
                    graph,
                    bundle,
                    extra={
                        "check": "beta_up_eq_ball_bruteforce",
                        "d": d,
                        "beta_up": step.beta_up,
                        "ball_bf": ball_bf,
                    },
                )
            )

    for prev, nxt in itertools.pairwise(profile):
        tally.checked += 1
        if prev.shell_max > nxt.shell_max + TOL:
            tally.violations += 1
            violations.append(
                _witness(
                    "P3",
                    graph,
                    bundle,
                    extra={
                        "check": "shell_max_monotone",
                        "d_prev": prev.d,
                        "shell_max_prev": prev.shell_max,
                        "d_next": nxt.d,
                        "shell_max_next": nxt.shell_max,
                    },
                )
            )


def _check_p4_and_p5(
    graph: GraphInstance,
    bundle: ContextBundle,
    b_cache: dict[MPDAG, float],
    tallies: dict[str, PredicateTally],
    violations: list[dict[str, Any]],
) -> None:
    """P4 (algorithm agreement) and P5 (nesting), sharing one ``eps`` grid.

    P4 (Theorem C.4/C.5): ``r_epsilon`` under ``"incremental"``,
    ``"bisection"`` and ``"hybrid"`` must agree with each other and with an
    independent brute-force BFS over the whole enumerated space (predicate
    ``B(G) > eps``, via :func:`bkrobust.core.spacelib.radius`).

    P5 (Theorem C.6): ``r_val <= r_0 <= r_eps`` and ``eps -> r_eps`` is
    non-decreasing, both checked against the brute-force radii P4 already
    computed (so P5 is not merely re-testing the staircase code).
    """
    p4, p5 = tallies["P4"], tallies["P5"]
    ctx, space, dists = bundle.ctx, graph.space, graph.dists
    eps_grid = _eps_grid(b_cache.values())

    r_bf_by_eps: dict[float, int] = {}
    for eps in eps_grid:
        r_inc = r_epsilon(ctx, graph.g0, eps, r_val=graph.r_val_bf, strategy="incremental").radius
        r_bis = r_epsilon(ctx, graph.g0, eps, r_val=graph.r_val_bf, strategy="bisection").radius
        r_hyb = r_epsilon(ctx, graph.g0, eps, r_val=graph.r_val_bf, strategy="hybrid").radius
        r_bf, bf_witness = radius(space, dists, lambda g, eps=eps: b_cache[g] > eps)
        r_bf_by_eps[eps] = r_bf

        p4.checked += 1
        if not (r_inc == r_bis == r_hyb == r_bf):
            p4.violations += 1
            violations.append(
                _witness(
                    "P4",
                    graph,
                    bundle,
                    extra={
                        "eps": eps,
                        "r_incremental": r_inc,
                        "r_bisection": r_bis,
                        "r_hybrid": r_hyb,
                        "r_bruteforce": r_bf,
                        "bruteforce_witness": (
                            bf_witness.edge_string() if bf_witness is not None else None
                        ),
                    },
                )
            )

    p5.checked += 1
    r0 = r_bf_by_eps[0.0]
    if not _radius_leq(graph.r_val_bf, r0):
        p5.violations += 1
        violations.append(
            _witness(
                "P5",
                graph,
                bundle,
                extra={"check": "r_val_leq_r0", "r_val_bf": graph.r_val_bf, "r_0": r0},
            )
        )

    ordered = sorted(r_bf_by_eps.items())
    for (eps_a, r_a), (eps_b, r_b) in itertools.pairwise(ordered):
        p5.checked += 1
        if not _radius_leq(r_a, r_b):
            p5.violations += 1
            violations.append(
                _witness(
                    "P5",
                    graph,
                    bundle,
                    extra={
                        "check": "eps_monotone",
                        "eps_a": eps_a,
                        "r_a": r_a,
                        "eps_b": eps_b,
                        "r_b": r_b,
                    },
                )
            )


def _check_p6(
    graph: GraphInstance,
    bundle: ContextBundle,
    profile: Sequence[ShellBias],
    tallies: dict[str, PredicateTally],
    violations: list[dict[str, Any]],
) -> None:
    r"""P6: the certificate audit (Theorem D) -- the important one.

    ``k`` is the number of **false orientations in the closure**
    ``K_{G0}``, not the number of asserted claims: Meek's rules can cascade a
    handful of false assertions into a larger false closure. For every
    ``S`` a superset of the false set (the only ``S`` for which
    ``Meek(cpdag, K_{G0}\\S)`` can contain the truth at all -- see the module
    docstring), both forms of the certificate are checked in **effect
    units** throughout, never relative ones:

    * ``|theta_Z - tau_true| <= beta_up(|S|)`` (the budget form), and
    * ``|theta_Z - tau_true| <= B(H)`` for ``H = Meek(cpdag, K_{G0}\\S)``
      (the stronger, per-state form).

    ``|K_{G0}| <= 7`` throughout this sweep, so enumerating every valid
    superset of the false set costs at most ``2**7 = 128`` states -- cheap
    enough to be exhaustive rather than sampling one ``S``.
    """
    tally = tallies["P6"]
    ctx = bundle.ctx
    k_g0 = graph.k_g0
    false_set = set(graph.false_orientations)
    remaining_true = [e for e in k_g0 if e not in false_set]
    true_dir = graph.dag.directed_edges
    realised = abs(ctx.theta_z - ctx.tau_true)
    beta_up_by_d = {step.d: step.beta_up for step in profile}

    for extra_drop in _powerset(remaining_true):
        s_set = false_set | set(extra_drop)
        k = len(s_set)
        keep = [e for e in k_g0 if e not in s_set]

        h = apply_orientations(graph.cpdag, keep)
        if h is None or not (set(keep) <= true_dir):
            # Would indicate a bug in this harness (S is constructed to be a
            # superset of the false set, so `keep` is a subset of the true
            # orientations and must be both consistent and truth-preserving
            # by Proposition R) rather than a theorem violation -- but is
            # tallied and witnessed rather than silently asserted, since a
            # silent assert would hide exactly this kind of harness bug.
            tally.checked += 1
            tally.violations += 1
            violations.append(
                _witness(
                    "P6",
                    graph,
                    bundle,
                    extra={
                        "check": "harness_state_or_truth_consistency",
                        "s": [list(e) for e in sorted(s_set)],
                        "k": k,
                        "keep": [list(e) for e in sorted(keep)],
                        "h_is_none": h is None,
                    },
                )
            )
            continue

        b_h = bias_at(ctx, h).worst
        beta_up_k = beta_up_by_d.get(k)

        tally.checked += 1
        if realised > b_h + TOL:
            tally.violations += 1
            violations.append(
                _witness(
                    "P6",
                    graph,
                    bundle,
                    g=h,
                    extra={
                        "check": "realised_leq_B_H",
                        "s": [list(e) for e in sorted(s_set)],
                        "k": k,
                        "realised_error": realised,
                        "B_H": b_h,
                    },
                )
            )

        tally.checked += 1
        if beta_up_k is not None and realised > beta_up_k + TOL:
            tally.violations += 1
            violations.append(
                _witness(
                    "P6",
                    graph,
                    bundle,
                    g=h,
                    extra={
                        "check": "realised_leq_beta_up_k",
                        "s": [list(e) for e in sorted(s_set)],
                        "k": k,
                        "realised_error": realised,
                        "beta_up_k": beta_up_k,
                    },
                )
            )


def _check_p7(
    graph: GraphInstance,
    tallies: dict[str, PredicateTally],
    violations: list[dict[str, Any]],
) -> None:
    """P7: semi-local possible-parent-sets agree with enumeration (Theorem E).

    Purely structural -- no SEM involved -- so this runs once per accepted
    graph instance rather than once per SEM draw, over every element of the
    enumerated space.
    """
    tally = tallies["P7"]
    for g in graph.space.elements:
        tally.checked += 1
        enumerated = possible_parent_sets_enumerated(g, graph.x)
        semilocal = possible_parent_sets_semilocal(graph.cpdag, g, graph.x)
        if enumerated != semilocal:
            tally.violations += 1
            violations.append(
                {
                    "predicate": "P7",
                    "instance_id": graph.instance_id,
                    "cpdag": graph.cpdag.edge_string(),
                    "g": g.edge_string(),
                    "x": graph.x,
                    "enumerated": sorted(sorted(s) for s in enumerated),
                    "semilocal": sorted(sorted(s) for s in semilocal),
                }
            )


# --- orchestration -----------------------------------------------------------


def run_verification(
    *,
    target_instances: int = 400,
    n_sem_draws: int = 3,
    base_seed: int = 20260916,
    max_attempts: int | None = None,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Run the full P1-P7 sweep and return the manifest and violation list.

    Args:
        target_instances: Accepted graph instances to gather before
            stopping (each contributing ``n_sem_draws`` SEM contexts).
        n_sem_draws: SEM draws per accepted graph instance.
        base_seed: Root seed; every instance's and every SEM draw's RNG is
            derived from ``(base_seed, ...)`` alone.
        max_attempts: Attempts (accepted or rejected) before giving up.
            Defaults to ``target_instances * 80``, generous headroom for
            ``decoupled_backdoor`` draws, which the generator's own
            docstring documents as almost always failing the degeneracy
            gate's sanity check by construction.

    Returns:
        ``(manifest, violations)``, both JSON-serialisable.
    """
    if max_attempts is None:
        max_attempts = target_instances * 80

    tallies: dict[str, PredicateTally] = {p: PredicateTally() for p in PREDICATES}
    violations: list[dict[str, Any]] = []
    rejection_counts: dict[str, int] = {}
    timings: dict[str, float] = {}
    generator_counts: dict[str, int] = {}
    space_sizes: list[int] = []
    k_g0_sizes: list[int] = []
    n_false_total = 0
    n_finite_r_val = 0
    r_val_sum = 0

    def _tick(name: str, seconds: float) -> None:
        timings[name] = timings.get(name, 0.0) + seconds

    t_start = time.perf_counter()
    n_accepted = 0
    attempt = 0
    while n_accepted < target_instances and attempt < max_attempts:
        index = attempt
        attempt += 1
        t0 = time.perf_counter()
        try:
            graph = _generate_graph_instance(index, base_seed)
        except RejectedError as exc:
            rejection_counts[exc.reason] = rejection_counts.get(exc.reason, 0) + 1
            continue
        finally:
            _tick("generation", time.perf_counter() - t0)

        n_accepted += 1
        generator_counts[graph.generator] = generator_counts.get(graph.generator, 0) + 1
        space_sizes.append(len(graph.space))
        k_g0_sizes.append(len(graph.k_g0))
        n_false_total += len(graph.false_orientations)
        if graph.r_val_bf != UNREACHED:
            n_finite_r_val += 1
            r_val_sum += graph.r_val_bf

        t0 = time.perf_counter()
        _check_p7(graph, tallies, violations)
        _tick("P7", time.perf_counter() - t0)

        for draw in range(n_sem_draws):
            sem_seed = (base_seed, index, draw, 0x5EED)
            sem_rng = np.random.default_rng(list(sem_seed))
            sem = random_sem(graph.dag, sem_rng)
            ctx = make_context(sem, graph.cpdag, graph.x, graph.y, graph.z)
            bundle = ContextBundle(
                context_id=f"{graph.instance_id}_sem{draw}",
                ctx=ctx,
                sem=sem,
                sem_seed=sem_seed,
            )

            t0 = time.perf_counter()
            b_cache = {g: bias_at(ctx, g).worst for g in graph.space.elements}
            _tick("bias_evaluation", time.perf_counter() - t0)

            t0 = time.perf_counter()
            _check_p1(graph, bundle, b_cache, tallies, violations)
            _tick("P1", time.perf_counter() - t0)

            t0 = time.perf_counter()
            _check_p2(graph, bundle, b_cache, tallies, violations)
            _tick("P2", time.perf_counter() - t0)

            t0 = time.perf_counter()
            profile = bias_profile(ctx, graph.g0, start=0, max_d=len(graph.k_g0))
            _tick("bias_profile", time.perf_counter() - t0)

            t0 = time.perf_counter()
            _check_p3(graph, bundle, profile, b_cache, tallies, violations)
            _tick("P3", time.perf_counter() - t0)

            t0 = time.perf_counter()
            _check_p4_and_p5(graph, bundle, b_cache, tallies, violations)
            _tick("P4_P5", time.perf_counter() - t0)

            t0 = time.perf_counter()
            _check_p6(graph, bundle, profile, tallies, violations)
            _tick("P6", time.perf_counter() - t0)

    total_seconds = time.perf_counter() - t_start

    manifest: dict[str, Any] = {
        "tolerance": TOL,
        "base_seed": base_seed,
        "target_instances": target_instances,
        "n_sem_draws_per_instance": n_sem_draws,
        "n_graph_instances_accepted": n_accepted,
        "n_contexts_total": n_accepted * n_sem_draws,
        "n_attempts": attempt,
        "max_attempts": max_attempts,
        "rejection_counts": dict(sorted(rejection_counts.items())),
        "generator_counts": dict(sorted(generator_counts.items())),
        "predicates": {
            p: {"checked": t.checked, "violations": t.violations} for p, t in tallies.items()
        },
        "n_total_violations": sum(t.violations for t in tallies.values()),
        "diagnostics": {
            "space_size_min": min(space_sizes) if space_sizes else None,
            "space_size_max": max(space_sizes) if space_sizes else None,
            "space_size_mean": (sum(space_sizes) / len(space_sizes)) if space_sizes else None,
            "k_g0_size_min": min(k_g0_sizes) if k_g0_sizes else None,
            "k_g0_size_max": max(k_g0_sizes) if k_g0_sizes else None,
            "k_g0_size_mean": (sum(k_g0_sizes) / len(k_g0_sizes)) if k_g0_sizes else None,
            "n_false_orientations_total": n_false_total,
            "n_instances_with_finite_r_val": n_finite_r_val,
            "r_val_mean_over_finite": (r_val_sum / n_finite_r_val) if n_finite_r_val else None,
            "r_val_unreached_note": (
                "UNREACHED (-1) instances are excluded from the mean above, per house "
                "style: UNREACHED is a status, never averaged."
            ),
        },
        "space_construction": (
            "bkrobust.search.space_fixed.build_space_correct (the CORRECTED "
            "enumerator). bkrobust.core.spacelib.build_space is deliberately not "
            "used: it calls bkrobust.demo.space.enumerate_space, whose chordality "
            "validity filter omits legitimate knowledge states -- see "
            "bkrobust.search.space_fixed's module docstring."
        ),
        "radius_convention": RADIUS_CONVENTION,
        "interpreter": sys.version,
        "executable": sys.executable,
        "git_sha": git_sha(short=False),
        "git_dirty": is_dirty(),
        "timings_seconds": {k: round(v, 4) for k, v in sorted(timings.items())},
        "total_seconds": round(total_seconds, 4),
    }
    return manifest, violations


def _default_out_dir() -> Path:
    """``<repo root>/results/epsilon/verification``."""
    return Path(__file__).resolve().parents[3] / "results" / "epsilon" / "verification"


def main(argv: Sequence[str] | None = None) -> None:
    """CLI entry point: run the sweep and write the manifest and violations.

    Args:
        argv: Command-line arguments, or ``None`` to use ``sys.argv``.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else "")
    parser.add_argument("--target", type=int, default=400, help="accepted graph instances")
    parser.add_argument("--sem-draws", type=int, default=3, help="SEM draws per instance")
    parser.add_argument("--seed", type=int, default=20260916, help="root seed")
    parser.add_argument("--max-attempts", type=int, default=None)
    parser.add_argument("--out", type=str, default=None, help="output directory")
    args = parser.parse_args(argv)

    out_dir = Path(args.out) if args.out else _default_out_dir()
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest, violations = run_verification(
        target_instances=args.target,
        n_sem_draws=args.sem_draws,
        base_seed=args.seed,
        max_attempts=args.max_attempts,
    )

    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
    (out_dir / "violations.json").write_text(json.dumps(violations, indent=2, sort_keys=True))
    print(json.dumps(manifest, indent=2, sort_keys=True))
    print(f"\n{len(violations)} violation(s) written to {out_dir / 'violations.json'}")


if __name__ == "__main__":
    main()
