"""The degeneracy gate and the Axis A instance runner.

Three pieces, in the order a single instance passes through them:

1. :func:`gate` -- a purely structural screen on ``(dag, cpdag, treatment,
   outcome)``, before any background knowledge is drawn. Mirrors the five
   design gates of :mod:`bkrobust.demo.example` (which screened the one
   hand-picked worked example); here the same idea runs automatically over
   every randomly generated instance.
2. :func:`run_instance` -- builds one instance end to end: draws ``K_true``,
   corrupts it to ``K_assumed``, forms ``G0``, fixes ``Z``, enumerates the
   space, and computes the radii. Raises :class:`GateRejectedError` for every
   rejection reason (see its docstring for the full list) rather than
   returning a sentinel, so its return type is honestly "an
   :class:`Instance`, always".
3. :func:`run_grid` -- sweeps a list of parameter points, catches
   :class:`GateRejectedError` (and any other exception, recorded as
   ``"run_error"`` rather than left to kill the sweep), tallies rejections by
   reason, and writes accepted instances incrementally via
   :class:`~bkrobust.core.resultsio.ResultWriter` so a crash loses at most the
   row in flight.

Determinism: the only randomness anywhere in this module is
``numpy.random.default_rng(seed)`` at the top of :func:`run_instance`, one
per instance, and every root seed in :func:`run_grid` is a deterministic
function of ``(root_seed, grid_index)`` -- never of how much of the grid has
already run -- which is what makes resuming safe.
"""

from __future__ import annotations

import csv
import time
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np

from bkrobust.core.instance import Instance, describe_graphs
from bkrobust.core.oracle import bias_stats, is_optimal, is_valid
from bkrobust.core.resultsio import ResultWriter, write_manifest
from bkrobust.core.spacelib import Space, build_space, distances_from, radius
from bkrobust.demo.evaluate import (
    all_valid_adjustment_sets_mpdag,
    is_valid_adjustment_set_mpdag,
    optimal_adjustment_set_mpdag,
)
from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover
from bkrobust.demo.graph import MPDAG, undirected_components
from bkrobust.demo.meek import apply_orientations
from bkrobust.synth.generators import GENERATORS
from bkrobust.synth.knowledge import CORRUPTIONS, check_consistent, draw_k_true, tiered

Edge = tuple[str, str]

#: Default epsilon grid for the bias-based radius, expressed as a fraction of
#: mean |true effect| would be more principled, but the oracle reports
#: mean_abs_bias directly, so these are plain absolute thresholds.
DEFAULT_EPSILONS: tuple[float, ...] = (0.02, 0.05, 0.1, 0.2)

#: Hard cap on a CPDAG's *total* undirected-edge count before it is rejected
#: rather than handed to :func:`~bkrobust.core.spacelib.build_space`.
#:
#: ``enumerate_space`` (which ``build_space`` calls) brute-forces all ``3**k``
#: orientation assignments over *every* undirected edge in the CPDAG at once
#: -- not per chordal component -- and, for each candidate, itself enumerates
#: that candidate's own consistent DAG extensions (up to ``2**k`` more work).
#: Measured on this machine (see the runner test suite and the pilot-run
#: notes), ``k=9`` costs about a second and ``k=10`` already costs ~12
#: seconds for one instance; growth from there is steep enough that an
#: unbounded ``k`` can hang a grid sweep on a single unlucky draw. ``n<=8``
#: generators keep ``k`` small on average, but not always (a near-complete
#: skeleton with few v-structures can leave most of it undirected), so this
#: is checked explicitly rather than assumed.
MAX_UNDIRECTED_EDGES: int = 9


class GateRejectedError(Exception):
    """Raised by :func:`run_instance` when an instance is rejected, at any stage.

    Attributes:
        reason: A short machine-readable rejection code. The five structural
            codes come from :func:`gate`; four more are raised by
            :func:`run_instance` itself: ``"cpdag_too_large_for_bfs"`` (before
            knowledge is even drawn -- see :data:`MAX_UNDIRECTED_EDGES`),
            ``"k_assumed_inconsistent"`` and ``"z_not_identified"`` (once
            background knowledge has been drawn), ``"z_invalid_at_g0"`` -- the
            unconditional guard on the radius convention's ``r_val >= 1``
            (see its raise site: the optimal-set formula can silently return
            an invalid empty set when there is no causal path, which
            ``gate``'s ``"no_causal_path"`` check ordinarily prevents but this
            catches regardless), and ``"g0_not_in_space"`` -- a rare (~0.1%,
            measured) mismatch where ``apply_orientations`` returns a ``G0``
            that :func:`~bkrobust.demo.space.enumerate_space` nonetheless
            omits from the space; see the note at its raise site.
            :func:`run_grid` additionally uses ``"run_error"`` for any other,
            unexpected exception from :func:`run_instance` on one grid point,
            so a single bad instance cannot kill the whole sweep.
    """

    def __init__(self, reason: str) -> None:
        super().__init__(reason)
        self.reason = reason


# --- the degeneracy gate ------------------------------------------------


def gate(dag: MPDAG, cpdag: MPDAG, treatment: str, outcome: str) -> tuple[bool, str]:
    """Screen ``(dag, cpdag, treatment, outcome)`` for the ways an instance can be vacuous.

    Checked in order, each short-circuiting the rest:

    1. ``treatment`` is not in, or adjacent to, any undirected (chordal)
       component of ``cpdag`` -- no perturbation of background knowledge can
       ever reach it, so the whole exercise is moot for this pair.
    2. ``outcome`` is not a descendant of ``treatment`` in the true DAG -- no
       causal path exists, so the true total effect is exactly zero and the
       pair is uninformative about bias scale. This is primarily prevented at
       selection time (:func:`_pick_treatment_outcome` only proposes pairs
       with a causal path), but is checked again here as a backstop: it is
       also what makes check 3 below sound, since
       ``optimal_adjustment_set_mpdag`` computes the optimal set from
       ``cn(x, y) = descendants(x) & (ancestors(y) | {y})`` and silently
       returns the *empty* set whenever ``cn`` is empty -- which is only
       trivially valid (see check 4) when there is also no confounding; with
       a causal path ruled out here, the run_instance-level
       ``z_invalid_at_g0`` check catches the remaining, rarer case where that
       still happens to be wrong.
    3. No valid adjustment set exists at all for ``(treatment, outcome)`` in
       the **true DAG**. Checked against ``dag``, not ``cpdag``: ``cpdag`` is
       typically still ambiguous about the treatment/outcome relationship
       itself (background knowledge has not been imposed yet), so requiring a
       set to be valid in *every* CPDAG extension at this stage would reject
       almost everything for the wrong reason -- lack of identification, not
       lack of genuine confounding, which is what this check is actually
       asking about. Mirrors gate G3 of :mod:`bkrobust.demo.example`, which is
       likewise phrased "in the truth".
    4. The empty set is itself a valid adjustment set in the true DAG --
       there is no confounding to speak of, so validity can never fail
       regardless of how knowledge is perturbed.
    5. Sanity: no *single* atomic perturbation of the true knowledge
       (``knowledge_to_recover(dag, cpdag)`` with one claim dropped) changes
       the validity of any candidate valid adjustment set. This mirrors gate
       G5 of :mod:`bkrobust.demo.example`, generalised to run automatically.

    Args:
        dag: The ground-truth DAG.
        cpdag: Its CPDAG.
        treatment: Treatment node.
        outcome: Outcome node.

    Returns:
        ``(True, "ok")`` if every check passes, else ``(False, reason)`` with
        one of ``"treatment_not_in_or_adjacent_to_component"``,
        ``"no_causal_path"``, ``"no_valid_adjustment_set"``,
        ``"empty_set_trivially_valid"``, or
        ``"no_atomic_perturbation_changes_validity"``.
    """
    comps = undirected_components(cpdag)
    touches_component = any(
        treatment in comp or any(cpdag.has_edge(treatment, node) for node in sorted(comp))
        for comp in comps
    )
    if not touches_component:
        return False, "treatment_not_in_or_adjacent_to_component"

    if outcome not in dag.descendants(treatment):
        return False, "no_causal_path"

    valid_sets = all_valid_adjustment_sets_mpdag(dag, treatment, outcome)
    if not valid_sets:
        return False, "no_valid_adjustment_set"
    if frozenset() in valid_sets:
        return False, "empty_set_trivially_valid"

    k_true = sorted(knowledge_to_recover(dag, cpdag))
    sanity = False
    for drop in k_true:
        remaining = [e for e in k_true if e != drop]
        g0 = apply_orientations(cpdag, remaining)
        if g0 is None:
            continue
        for z in valid_sets:
            if not is_valid_adjustment_set_mpdag(g0, treatment, outcome, z):
                sanity = True
                break
        if sanity:
            break
    if not sanity:
        return False, "no_atomic_perturbation_changes_validity"

    return True, "ok"


# --- treatment/outcome selection for the generic generators -----------------


def _pick_treatment_outcome(
    generator_name: str, dag: MPDAG, rng: np.random.Generator
) -> tuple[str, str]:
    """Choose ``(treatment, outcome)`` for a generated DAG.

    ``decoupled_backdoor_dag`` names its own treatment/outcome nodes
    (``"X"``, ``"Y"``) by construction. For the other generators, a random
    ordered pair with ``outcome`` a genuine descendant of ``treatment`` is
    drawn -- guaranteeing a well-defined non-trivial total effect to estimate
    -- by scanning a single random permutation of all ordered pairs
    deterministically until one qualifies.

    Args:
        generator_name: Key into :data:`bkrobust.synth.generators.GENERATORS`.
        dag: The generated DAG.
        rng: Sole source of randomness.

    Returns:
        ``(treatment, outcome)``.
    """
    if generator_name == "decoupled_backdoor":
        return "X", "Y"

    nodes = sorted(dag.nodes)
    n = len(nodes)
    pairs = [(i, j) for i in range(n) for j in range(n) if i != j]
    order = rng.permutation(len(pairs))
    for idx in order:
        i, j = pairs[int(idx)]
        x, y = nodes[i], nodes[j]
        if y in dag.descendants(x):
            return x, y
    # No ordered pair has a directed path at all (e.g. a near-empty graph):
    # fall back to a fixed, still-deterministic pair.
    return nodes[0], nodes[1 % n]


# --- one instance ---------------------------------------------------------


def run_instance(
    instance_id: str,
    generator_name: str,
    generator_params: dict[str, Any],
    n: int,
    seed: int,
    knows_fraction: float,
    corruption_name: str,
    corruption_rate: float,
    epsilons: Sequence[float] = DEFAULT_EPSILONS,
    n_bias_draws: int = 25,
    tier_params: dict[str, Any] | None = None,
) -> Instance:
    """Build one Axis A instance end to end.

    Pipeline: generate the DAG and its CPDAG, pick ``(treatment, outcome)``,
    run the degeneracy :func:`gate`, draw ``K_true``, corrupt it to
    ``K_assumed``, form ``G0 = apply_orientations(cpdag, K_assumed)``, fix
    ``Z = optimal_adjustment_set_mpdag(G0, x, y)``, enumerate the space, BFS
    distances from ``G0``, and compute ``r_val``, ``r_opt`` and ``r_eps``
    (keyed off ``mean_abs_bias``, never ``max_abs_bias`` -- the latter is a
    sampled, unstable proxy and is house style to never key a radius on).

    Args:
        instance_id: Identifier stored on the returned :class:`Instance`.
        generator_name: Key into
            :data:`bkrobust.synth.generators.GENERATORS`.
        generator_params: Keyword arguments for that generator (besides
            ``n`` and ``rng``).
        n: Node count. Kept small (5-8) by callers so BFS ground truth stays
            tractable; not itself validated here (the generators validate
            their own preconditions on ``n``).
        seed: Root seed for this instance's ``numpy.random.Generator``. Every
            draw in this function's call tree descends from it alone.
        knows_fraction: Passed to
            :func:`bkrobust.synth.knowledge.draw_k_true`.
        corruption_name: One of ``"none"`` (no corruption: ``K_assumed =
            K_true``), a key of
            :data:`bkrobust.synth.knowledge.CORRUPTIONS` other than
            ``"tiered"`` (called as ``fn(k_true, rng, corruption_rate)``), or
            ``"tiered"`` (called as ``tiered(dag, cpdag, rng, **tier_params)``,
            ignoring ``k_true`` and ``corruption_rate`` -- tiered knowledge is
            generative, not corruptive; see :mod:`bkrobust.synth.knowledge`).
        corruption_rate: Rate passed to the non-tiered corruptions. Ignored
            when ``corruption_name == "tiered"``.
        epsilons: The epsilon grid for ``r_eps``.
        n_bias_draws: SEM coefficient draws per extension, passed to
            :func:`bkrobust.core.oracle.bias_stats`.
        tier_params: Keyword arguments for :func:`~bkrobust.synth.knowledge.tiered`
            (``n_tiers`` and optionally ``corruption_rate``) when
            ``corruption_name == "tiered"``. Defaults to
            ``{"n_tiers": 3, "corruption_rate": corruption_rate}``.

    Returns:
        A populated :class:`~bkrobust.core.instance.Instance`.

    Raises:
        GateRejectedError: On any of the rejection reasons documented on
            :class:`GateRejectedError` itself.
        ValueError: If ``generator_name`` is not a known generator.
    """
    if generator_name not in GENERATORS:
        raise ValueError(f"unknown generator {generator_name!r}")

    rng = np.random.default_rng(seed)
    timings: dict[str, float] = {}

    t = time.perf_counter()
    dag = GENERATORS[generator_name](n, rng, **generator_params)
    cpdag = dag_to_cpdag(dag)
    treatment, outcome = _pick_treatment_outcome(generator_name, dag, rng)
    timings["graph"] = time.perf_counter() - t

    if len(cpdag.undirected_edges) > MAX_UNDIRECTED_EDGES:
        raise GateRejectedError("cpdag_too_large_for_bfs")

    t = time.perf_counter()
    accepted, reason = gate(dag, cpdag, treatment, outcome)
    timings["gate"] = time.perf_counter() - t
    if not accepted:
        raise GateRejectedError(reason)

    t = time.perf_counter()
    k_true = draw_k_true(dag, cpdag, rng, knows_fraction)
    if corruption_name == "none":
        k_assumed = list(k_true)
        applied_rate = 0.0
    elif corruption_name == "tiered":
        tp = dict(tier_params or {})
        tp.setdefault("n_tiers", 3)
        tp.setdefault("corruption_rate", corruption_rate)
        applied_rate = tp["corruption_rate"]
        k_assumed = tiered(dag, cpdag, rng, **tp)
    else:
        if corruption_name not in CORRUPTIONS:
            raise ValueError(f"unknown corruption {corruption_name!r}")
        k_assumed = CORRUPTIONS[corruption_name](k_true, rng, corruption_rate)
        applied_rate = corruption_rate
    timings["knowledge"] = time.perf_counter() - t

    consistent = check_consistent(cpdag, k_assumed)
    if not consistent:
        raise GateRejectedError("k_assumed_inconsistent")

    t = time.perf_counter()
    g0 = apply_orientations(cpdag, k_assumed)
    assert g0 is not None  # guaranteed by check_consistent above
    z = optimal_adjustment_set_mpdag(g0, treatment, outcome)
    timings["g0_and_z"] = time.perf_counter() - t
    if z is None:
        raise GateRejectedError("z_not_identified")
    zf = frozenset(z)

    # Unconditional guard on the radius convention (bkrobust.core.conventions):
    # Z is read off G0, so it must be valid there and r_val must be >= 1.
    # optimal_adjustment_set_mpdag computes Z structurally, from
    # cn(x, y) = descendants(x) & (ancestors(y) | {y}); when cn is empty (no
    # causal path -- gated above as "no_causal_path" using the true dag, which
    # ordinarily prevents this) it returns the empty set regardless of whether
    # that is actually backdoor-valid, which it is not whenever a backdoor
    # path is still open. Checked directly here, unconditionally (not just for
    # empty Z), rather than trusted from the formula.
    if not is_valid(zf, g0, treatment, outcome):
        raise GateRejectedError("z_invalid_at_g0")

    t = time.perf_counter()
    space = build_space(cpdag)
    if g0 not in space.neighbours:
        # A rare (~0.1%, see the note on GateRejectedError.reason) mismatch
        # between apply_orientations/meek_closure and enumerate_space's own
        # is_valid_mpdag filter: g0 can be Meek-closed, keep every compelled
        # edge, and satisfy model inclusion in cpdag, yet still fail
        # is_valid_mpdag's chordality check for a reason specific to that
        # check -- it inspects only the undirected subgraph, so an
        # undirected component can look like a chordless 4-cycle there while
        # its chord exists as a *directed* edge. enumerate_space then omits
        # g0 from the space outright. This is a property of the frozen core
        # (bkrobust.core / bkrobust.demo), not a bug in this module, and is
        # deliberately not "fixed" here -- see runner test suite for the
        # captured case. Gate rather than let distances_from's KeyError
        # propagate and kill an entire grid sweep over one instance.
        raise GateRejectedError("g0_not_in_space")
    dists = distances_from(space, g0)
    timings["space"] = time.perf_counter() - t

    t = time.perf_counter()
    r_val, _ = radius(space, dists, lambda g: not is_valid(zf, g, treatment, outcome))
    r_opt, _ = radius(space, dists, lambda g: not is_optimal(zf, g, treatment, outcome))
    timings["radius_val_opt"] = time.perf_counter() - t

    t = time.perf_counter()
    r_eps = _radius_eps(space, dists, zf, treatment, outcome, epsilons, rng, n_bias_draws)
    timings["radius_eps"] = time.perf_counter() - t

    descriptors = describe_graphs(cpdag, k_assumed, treatment)

    return Instance(
        instance_id=instance_id,
        generator=generator_name,
        seed=seed,
        params=dict(generator_params),
        true_dag=dag.edge_string(),
        cpdag=cpdag.edge_string(),
        k_true=tuple(sorted(k_true)),
        k_assumed=tuple(sorted(k_assumed)),
        corruption=corruption_name,
        corruption_rate=float(applied_rate),
        k_assumed_consistent=consistent,
        g0=g0.edge_string(),
        treatment=treatment,
        outcome=outcome,
        z=tuple(sorted(zf)),
        descriptors=descriptors,
        space_size=len(space),
        n_covers=len(space.covers),
        max_shell=max(dists.values()) if dists else 0,
        r_val=r_val,
        r_opt=r_opt,
        r_eps=r_eps,
        method="bfs_exact",
        timings=timings,
    )


def _radius_eps(
    space: Space,
    dists: dict[MPDAG, int],
    zf: frozenset[str],
    treatment: str,
    outcome: str,
    epsilons: Sequence[float],
    rng: np.random.Generator,
    n_bias_draws: int,
) -> dict[str, int]:
    """The epsilon-bias radius grid, computed with the mean bias evaluated once per graph.

    ``mean_abs_bias`` is sampled (SEM coefficients are drawn), so it is
    computed at most once per graph and cached, rather than redrawn for every
    epsilon threshold -- both cheaper and statistically cleaner (one draw per
    graph, thresholded at several epsilons, rather than a fresh draw per
    threshold). The cache is populated lazily in the traversal order
    :func:`~bkrobust.core.spacelib.radius` visits :attr:`Space.elements` in --
    a fixed tuple order, not a hash-affected one -- so RNG consumption stays
    reproducible under any ``PYTHONHASHSEED``.
    """
    bias_cache: dict[MPDAG, float] = {}

    def mean_bias_of(g: MPDAG) -> float:
        cached = bias_cache.get(g)
        if cached is not None:
            return cached
        stats = bias_stats(zf, g, treatment, outcome, rng, n_bias_draws)
        # [g] empty should not arise for a genuine space element (is_valid_mpdag
        # requires at least one extension), but treat it as a failure
        # defensively, matching is_valid/is_optimal's stance on empty [g].
        value = stats.mean_abs_bias if stats.n_evaluations > 0 else float("inf")
        bias_cache[g] = value
        return value

    r_eps: dict[str, int] = {}
    for eps in epsilons:
        r, _ = radius(space, dists, lambda g, eps=eps: mean_bias_of(g) > eps)
        r_eps[f"{eps:g}"] = r
    return r_eps


# --- the grid runner ------------------------------------------------------


def _instance_id(index: int) -> str:
    return f"inst_{index:06d}"


def run_grid(
    grid: Sequence[dict[str, Any]],
    out_dir: str | Path,
    root_seed: int,
    epsilons: Sequence[float] = DEFAULT_EPSILONS,
    n_bias_draws: int = 25,
    resume: bool = True,
) -> dict[str, Any]:
    """Sweep ``grid``, writing accepted instances incrementally and tallying rejections.

    Each element of ``grid`` is a dict of keyword arguments for
    :func:`run_instance`, minus ``instance_id``, ``seed`` and ``epsilons``/
    ``n_bias_draws`` (supplied by this function): ``generator_name``,
    ``generator_params``, ``n``, ``knows_fraction``, ``corruption_name``,
    ``corruption_rate``, and optionally ``tier_params``.

    Resume safety: grid point ``i`` always gets ``seed = root_seed + i`` --
    never a seed that depends on how much of the grid has run so far -- and
    is skipped only if ``inst_{i:06d}`` already has a row in
    ``out_dir/results.csv``. So a kill-and-restart recomputes at most the one
    row that was in flight when it died (:class:`~bkrobust.core.resultsio.ResultWriter`
    flushes after every row) plus every point that was rejected last time
    (cheap to redo, and their tally would otherwise be lost, since rejections
    are never written to the CSV).

    Args:
        grid: The parameter grid, one dict per instance to attempt.
        out_dir: Results directory. Gets ``results.csv`` and
            ``manifest.json``.
        root_seed: Root seed every instance's seed is derived from.
        epsilons: Epsilon grid forwarded to every :func:`run_instance` call.
        n_bias_draws: Forwarded to every :func:`run_instance` call.
        resume: If True (default) and ``results.csv`` already exists, append
            to it and skip grid points already present. If False, an
            existing ``results.csv`` is an error (see
            :class:`~bkrobust.core.resultsio.ResultWriter`).

    Returns:
        A summary dict: ``n_total``, ``n_accepted``, ``rejection_counts``
        (by reason, complete for the whole grid regardless of what was
        skipped via resume -- see above), and ``csv_path``.
    """
    out_dir = Path(out_dir)
    csv_path = out_dir / "results.csv"

    done_ids: set[str] = set()
    if resume and csv_path.exists():
        with csv_path.open(newline="") as fh:
            for row in csv.DictReader(fh):
                done_ids.add(row["instance_id"])

    rejection_counts: dict[str, int] = {}
    errors: list[dict[str, str]] = []
    n_accepted = 0
    n_total = len(grid)

    with ResultWriter(csv_path, resume=resume) as writer:
        for index, point in enumerate(grid):
            instance_id = _instance_id(index)
            if instance_id in done_ids:
                n_accepted += 1
                continue
            seed = root_seed + index
            try:
                instance = run_instance(
                    instance_id=instance_id,
                    generator_name=point["generator_name"],
                    generator_params=point.get("generator_params", {}),
                    n=point["n"],
                    seed=seed,
                    knows_fraction=point["knows_fraction"],
                    corruption_name=point["corruption_name"],
                    corruption_rate=point.get("corruption_rate", 0.0),
                    epsilons=epsilons,
                    n_bias_draws=n_bias_draws,
                    tier_params=point.get("tier_params"),
                )
            except GateRejectedError as exc:
                rejection_counts[exc.reason] = rejection_counts.get(exc.reason, 0) + 1
                continue
            except Exception as exc:
                # A single bad grid point must never kill the whole sweep --
                # that would defeat the entire point of writing rows
                # incrementally. Record it as a rejection (so the rate stays
                # visible) plus the exception text (so it stays diagnosable),
                # and move on. Anything that should actually stop the run
                # (KeyboardInterrupt, SystemExit) is not an Exception subclass
                # and is not caught here.
                rejection_counts["run_error"] = rejection_counts.get("run_error", 0) + 1
                errors.append(
                    {
                        "instance_id": instance_id,
                        "exception_type": type(exc).__name__,
                        "message": str(exc),
                    }
                )
                continue
            writer.write(instance.to_row())
            n_accepted += 1

    write_manifest(
        out_dir,
        seed=root_seed,
        grid={
            "n_points": n_total,
            "epsilons": list(epsilons),
            "n_bias_draws": n_bias_draws,
        },
        extra={
            "n_total": n_total,
            "n_accepted": n_accepted,
            "rejection_counts": rejection_counts,
            "errors": errors,
        },
    )
    return {
        "n_total": n_total,
        "n_accepted": n_accepted,
        "rejection_counts": dict(sorted(rejection_counts.items())),
        "errors": errors,
        "csv_path": str(csv_path),
    }
