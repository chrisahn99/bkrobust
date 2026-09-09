"""Survival curves under knowledge corruption (Axis Robustness, Phase 1).

An analyst holds a CPDAG ``Ĉ``, asserts background knowledge ``K``, Meek-closes
to ``G0``, and reads an adjustment set ``Z*`` off ``G0``. This module corrupts
``K`` two different ways (:func:`bkrobust.synth.knowledge.flip`,
:func:`bkrobust.synth.knowledge.tiered`) and asks whether ``Z*`` -- the
analyst's *committed* set, held fixed -- is still valid at the corrupted
state. See ``results/axis_robustness/PREREGISTRATION.md`` for the full design;
nothing here may deviate from it without a dated appendix there.

Gate discipline
----------------
The only gate used to decide admissibility of a corrupted committed graph is
:func:`bkrobust.benchmarks.measure.fast_gate` (imported, never
reimplemented); every row this module writes carries the literal string
``"fast_gate"`` in a ``gate`` column. :func:`bkrobust.synth.component_generator
.generate_instance` is used, per the pre-approved API contract, to draw the
underlying ``(dag, cpdag)`` pair and a screened truthful ``K``; that helper's
*own internal* screening step (documented on its docstring, not called out by
name here on purpose -- see ``SURVIVAL_RUN_NOTES.md`` for the full identifier
and the argument for why it is safe at these sizes) is bounded by
construction: component sizes here never exceed 12, so the enumeration it
performs is over at most ~12 candidate nodes on a *fully oriented* DAG (no
undirected edges, hence a single trivial "extension"), which is the cheap
regime that helper's module documents, not the unbounded one this project's
hard constraints forbid. This module's *own* code never imports or calls that
exponential-enumeration path, or the gate wrapping it, directly -- see the
run notes for the full argument.

RNG discipline
--------------
No global RNG is ever touched. Every stochastic draw goes through an
``np.random.Generator`` built by :func:`derived_seed` plus
``np.random.default_rng``. :func:`derived_seed` hashes its parts with
``hashlib.sha256`` -- never Python's salted built-in ``hash()`` -- so results
are identical across ``PYTHONHASHSEED`` values (checked by acceptance
criterion A3, not assumed).
"""

from __future__ import annotations

import hashlib
from typing import Any, Iterable

import numpy as np

from bkrobust.benchmarks.measure import (
    MAX_G0_UNDIRECTED_FOR_EXTENSIONS,
    O_INTRACTABLE,
    fast_gate,
)
from bkrobust.core.conventions import UNREACHED
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag
from bkrobust.demo.graph import MPDAG, canon
from bkrobust.demo.meek import apply_orientations
from bkrobust.gac import is_gac_valid_mpdag
from bkrobust.hybrid import HybridResult, breakdown_radius
from bkrobust.synth.component_generator import (
    OUTCOME,
    TREATMENT,
    ComponentSpec,
    achievable_separations,
    generate_instance,
)
from bkrobust.synth.knowledge import flip, tiered

Edge = tuple[str, str]

#: Literal gate name stamped on every output row (hard constraint #1).
GATE_NAME = "fast_gate"

#: Fixed fractional grid for the primary AUC endpoint (PREREGISTRATION.md §4).
FRAC_GRID: tuple[float, ...] = tuple(round(i / 10, 1) for i in range(1, 11))

#: Reason code used across both arms when the screened committed graph would
#: force an intractable extension enumeration downstream. Re-exported from
#: :mod:`bkrobust.benchmarks.measure` so both modules speak the same
#: vocabulary; not redefined.
REASON_OK = "ok"


# --- seeding -----------------------------------------------------------------


def derived_seed(*parts: Any) -> int:
    """A stable, ``PYTHONHASHSEED``-independent seed derived from ``parts``.

    Uses ``hashlib.sha256`` over the string form of ``parts`` -- never
    Python's built-in ``hash()``, which is salted per-process unless hash
    randomization is disabled -- so the same parts always give the same seed,
    on any interpreter, on any run.

    Args:
        *parts: Anything with a stable ``str()``. Order matters.

    Returns:
        A nonnegative integer, suitable for ``np.random.default_rng``.
    """
    payload = "|".join(str(p) for p in parts).encode("utf-8")
    digest = hashlib.sha256(payload).digest()
    return int.from_bytes(digest[:8], "big")


# --- grid helpers --------------------------------------------------------------


def separations_for(component_size: int) -> tuple[int, int, int]:
    """The (min, median, max) achievable separations for ``component_size``.

    ``achievable_separations`` returns ``1 .. component_size - 1`` ascending,
    which is always an odd-length run for the component sizes this study
    uses (6, 8, 10, 12), so the median index is exact.

    Args:
        component_size: Number of vertices in the undirected component.

    Returns:
        ``(min, median, max)``.
    """
    seps = sorted(achievable_separations(component_size))
    return (seps[0], seps[len(seps) // 2], seps[-1])


# --- graph-level predictors ----------------------------------------------------


def directed_symdiff(g1: MPDAG, g2: MPDAG) -> int:
    """``|dir(g1) Δ dir(g2)|``, the directed-edge symmetric difference."""
    return len(g1.directed_edges ^ g2.directed_edges)


def shd_cpdag_count(g0: MPDAG, cpdag: MPDAG) -> int:
    """Number of edges directed in ``g0`` but undirected in ``cpdag``."""
    undirected_keys = {canon(a, b) for a, b in cpdag.undirected_edges}
    return sum(1 for (a, b) in g0.directed_edges if canon(a, b) in undirected_keys)


def undirected_fraction_of(cpdag: MPDAG) -> float:
    """Fraction of ``cpdag``'s edges that are undirected."""
    n_dir = len(cpdag.directed_edges)
    n_und = len(cpdag.undirected_edges)
    total = n_dir + n_und
    return (n_und / total) if total else 0.0


# --- radius -----------------------------------------------------------------


def flatten_radius(result: HybridResult) -> dict[str, Any]:
    """Flatten a :class:`HybridResult` into CSV-safe columns.

    Sentinel discipline (hard constraint #5): a timed-out ladder
    (``result.exact is False``) leaves ``r_val`` **empty**, never a number,
    and stamps ``wall_until_timeout_s`` instead of any column named
    ``seconds``. ``UNREACHED`` (``-1``) is stored as the literal sentinel
    integer when the search *completed* and genuinely found no failure --
    that is a defined, meaningful outcome per
    ``bkrobust.core.conventions``, not an undefined one -- and is flagged via
    ``r_status == "unreached"`` so it is never silently averaged downstream.

    Args:
        result: The hybrid radius result.

    Returns:
        A flat dict of columns, all prefixed ``r_``.
    """
    stats = result.stats.as_dict()
    row: dict[str, Any] = {
        "r_method": result.method,
        "r_oracle": result.oracle,
        "r_exact": result.exact,
        "r_assumes": result.assumes,
        "r_witness": result.witness or "",
        "r_search_seconds": result.search_seconds,
        "r_ladder_seconds": result.ladder_seconds,
        "r_total_seconds": result.total_seconds,
        "r_stat_elements_visited": stats["elements_visited"],
        "r_stat_closures": stats["closures"],
        "r_stat_validity_checks": stats["validity_checks"],
        "r_stat_extensions_enumerated": stats["extensions_enumerated"],
    }
    if not result.exact:
        row["r_val"] = ""
        row["r_status"] = "timeout"
        row["wall_until_timeout_s"] = result.total_seconds
    elif result.radius == UNREACHED:
        row["r_val"] = UNREACHED
        row["r_status"] = "unreached"
        row["wall_until_timeout_s"] = ""
    else:
        row["r_val"] = result.radius
        row["r_status"] = "ok"
        row["wall_until_timeout_s"] = ""
    return row


def compute_r_val(
    cpdag: MPDAG, g0: MPDAG, x: str, y: str, z_star: Iterable[str], *, time_limit_s: float = 60.0
) -> dict[str, Any]:
    """Compute and flatten the breakdown radius for one instance's committed state."""
    result = breakdown_radius(cpdag, None, x, y, frozenset(z_star), g0=g0, time_limit_s=time_limit_s)
    return flatten_radius(result)


# --- instance admission (shared by both arms) ---------------------------------


def score_committed_graph(
    dag: MPDAG, cpdag: MPDAG, x: str, y: str, k: tuple[Edge, ...], g0: MPDAG | None
) -> tuple[dict[str, Any] | None, str]:
    """Screen one committed graph ``g0`` (Meek-closure of ``cpdag`` + ``k``).

    Shared admission logic for both the flip-arm (base-wrongness-corrupted
    truthful ``K``) and tiered-arm (reference tiering) instances.

    Args:
        dag: Ground-truth DAG.
        cpdag: Its CPDAG.
        x: Treatment.
        y: Outcome.
        k: The orientation claims that produced ``g0``.
        g0: ``apply_orientations(cpdag, k)``, or ``None`` if contradictory.

    Returns:
        ``(predictors, "ok")`` on admission, else ``(None, reason)``.
        ``predictors`` holds ``z_star`` (frozenset), ``n_k``, ``shd_truth``,
        ``shd_cpdag``, ``undirected_fraction``, ``z_size``.
    """
    if g0 is None:
        return None, "k_contradictory"
    if len(k) == 0:
        return None, "n_k_zero"
    if len(g0.undirected_edges) > MAX_G0_UNDIRECTED_FOR_EXTENSIONS:
        return None, O_INTRACTABLE
    z_star = optimal_adjustment_set_mpdag(g0, x, y)
    if z_star is None:
        return None, "optimal_set_undefined"
    z_star = frozenset(z_star)
    if not is_gac_valid_mpdag(g0, x, y, z_star):
        # Should not happen: the optimal set of g0 is valid at g0 by
        # construction. If this ever fires it is a bug or a mis-gated
        # degenerate instance -- surfaced as its own reason, not silently
        # dropped into "optimal_set_undefined".
        return None, "z_invalid_at_g0"
    predictors = {
        "z_star": z_star,
        "n_k": len(k),
        "shd_truth": directed_symdiff(g0, dag),
        "shd_cpdag": shd_cpdag_count(g0, cpdag),
        "undirected_fraction": undirected_fraction_of(cpdag),
        "z_size": len(z_star),
    }
    return predictors, REASON_OK


# --- flip arm ------------------------------------------------------------------


def build_flip_instance(
    component_size: int,
    separation: int,
    coverage: float,
    base_wrongness: float,
    seed: int,
) -> tuple[dict[str, Any] | None, str]:
    """Draw and screen one flip-arm instance.

    Pipeline: ``generate_instance`` (truthful ``K_g0``) -> require
    ``.accepted`` -> ``K_b = flip(K_g0, rng_base, base_wrongness)`` ->
    ``G0 = apply_orientations(cpdag, K_b)`` -> require
    ``fast_gate(dag, cpdag, x, y) == (True, "ok")`` -> :func:`score_committed_graph`.
    ``Z*`` is computed from this ``G0`` and held fixed for the whole sweep.

    Args:
        component_size: Target undirected-component size.
        separation: Target separation (must be achievable).
        coverage: Fraction of the component's undirected edges truthfully
            asserted before base-wrongness corruption.
        base_wrongness: Fraction of the truthful ``K`` reversed before the
            instance is accepted -- this is what makes ``shd_truth``
            non-degenerate.
        seed: Root seed passed straight to ``generate_instance`` (the design's
            "seeds 0, 1, 2, ...").

    Returns:
        ``(instance, "ok")`` on admission, else ``(None, reason)``. ``reason``
        is namespaced (``"generate_instance:<code>"``, ``"fast_gate:<code>"``,
        or one of :func:`score_committed_graph`'s codes) so exclusions can be
        counted by cause.
    """
    spec = ComponentSpec(component_size=component_size, separation=separation, coverage=coverage)
    inst = generate_instance(spec, seed)
    if not inst.accepted:
        return None, f"generate_instance:{inst.reject_reason}"

    dag, cpdag = inst.dag_obj, inst.cpdag_obj
    assert dag is not None and cpdag is not None
    x, y = inst.treatment, inst.outcome

    gate_ok, gate_reason = fast_gate(dag, cpdag, x, y)
    if not gate_ok:
        return None, f"fast_gate:{gate_reason}"

    rng_base = np.random.default_rng(derived_seed(inst.instance_id, "base", base_wrongness))
    k_b = tuple(flip(list(inst.k_g0), rng_base, base_wrongness))
    g0 = apply_orientations(cpdag, k_b)

    predictors, reason = score_committed_graph(dag, cpdag, x, y, k_b, g0)
    if predictors is None:
        return None, reason

    # `inst.instance_id` (from generate_instance) encodes component_size,
    # separation, coverage and seed, but NOT base_wrongness -- three flip
    # instances at the same (component_size, separation, coverage, seed) but
    # different base_wrongness are siblings drawn from the identical
    # underlying (dag, cpdag, k_g0) and would otherwise share one ID. Since
    # every downstream RNG draw for this instance (sample_flip_state) is
    # seeded from this ID, an un-qualified ID would also correlate the
    # corruption-position draws across those three siblings, not just create
    # an ambiguous join key. Base-wrongness is folded in here so the ID is
    # unique across the whole flip-arm population and every derived seed is
    # independent across base_wrongness. `base_instance_id` keeps the
    # unqualified id for anyone who wants to find the sibling instances.
    bw_tag = f"bw{round(base_wrongness * 100):03d}"
    instance_id = f"{inst.instance_id}_{bw_tag}"

    instance = {
        "instance_id": instance_id,
        "base_instance_id": inst.instance_id,
        "arm": "flip",
        "seed": seed,
        "component_size": component_size,
        "separation": separation,
        "coverage": coverage,
        "base_wrongness": base_wrongness,
        "n_tiers": "",
        "x": x,
        "y": y,
        "dag": dag,
        "cpdag": cpdag,
        "g0": g0,
        "k": k_b,
        **predictors,
    }
    return instance, REASON_OK


def sample_flip_state(instance: dict[str, Any], d: int, rep: int) -> dict[str, Any]:
    """Draw one flip-arm sample at targeted depth ``d``.

    Args:
        instance: A flip-arm instance from :func:`build_flip_instance`.
        d: Targeted number of reversed claims, ``1 <= d <= n_k``.
        rep: Repetition index (independent draw at the same ``d``).

    Returns:
        A row dict matching the unified samples schema (see
        :mod:`bkrobust.robustness.run_survival`).
    """
    k = instance["k"]
    n = instance["n_k"]
    assert 1 <= d <= n, f"d={d} out of range for n_k={n}"
    rate = d / n
    seed = derived_seed(instance["instance_id"], "flip", d, rep)
    rng = np.random.default_rng(seed)
    k_cor = flip(list(k), rng, rate)
    n_flipped = round(rate * n)
    assert n_flipped == d, f"flip rate/depth mismatch: rate={rate}, n={n}, got {n_flipped}, want {d}"

    g = apply_orientations(instance["cpdag"], k_cor)
    row: dict[str, Any] = {
        "instance_id": instance["instance_id"],
        "arm": "flip",
        "gate": GATE_NAME,
        "coverage": instance["coverage"],
        "base_wrongness": instance["base_wrongness"],
        "n_tiers": "",
        "d": d,
        "rep": rep,
        "flip_rate": rate,
        "corruption_rate": "",
        "seed": seed,
    }
    if g is None:
        row["status"] = "corrupted_k_contradictory"
        row["survived"] = ""
        row["symdiff_proxy_not_distance"] = ""
    else:
        survived = is_gac_valid_mpdag(g, instance["x"], instance["y"], instance["z_star"])
        row["status"] = "ok"
        row["survived"] = survived
        row["symdiff_proxy_not_distance"] = directed_symdiff(instance["g0"], g)
    return row


# --- tiered arm ------------------------------------------------------------------


def build_tiered_instance(
    component_size: int,
    separation: int,
    seed: int,
    n_tiers: int,
) -> tuple[dict[str, Any] | None, str]:
    """Draw and screen one tiered-arm instance.

    This is a *separate instance population* from the flip arm
    (PREREGISTRATION.md §1.1, §5.1): :func:`bkrobust.synth.knowledge.tiered`
    is generative, not corruptive, and ignores any ``K``/``coverage`` entirely.
    The underlying ``(dag, cpdag)`` pair is still drawn via
    ``generate_instance`` (coverage is irrelevant to ``dag``/``cpdag`` under
    the default ``coverage_order="spine_first"``, which consumes no RNG
    draws, so ``coverage=1.0`` is used here purely as an arbitrary constant,
    not a swept parameter) and screened by ``fast_gate`` exactly as in the
    flip arm. ``K_ref = tiered(dag, cpdag, rng0, n_tiers, 0.0)`` is this arm's
    truthful reference knowledge; ``G0_t`` and ``Z*_t`` are read off it and
    held fixed for the sweep.

    Args:
        component_size: Target undirected-component size.
        separation: Target separation (must be achievable).
        seed: Root seed passed to ``generate_instance``.
        n_tiers: Number of temporal tiers.

    Returns:
        ``(instance, "ok")`` on admission, else ``(None, reason)``.
    """
    spec = ComponentSpec(component_size=component_size, separation=separation, coverage=1.0)
    inst = generate_instance(spec, seed)
    if not inst.accepted:
        return None, f"generate_instance:{inst.reject_reason}"

    dag, cpdag = inst.dag_obj, inst.cpdag_obj
    assert dag is not None and cpdag is not None
    x, y = inst.treatment, inst.outcome

    gate_ok, gate_reason = fast_gate(dag, cpdag, x, y)
    if not gate_ok:
        return None, f"fast_gate:{gate_reason}"

    instance_id = f"{inst.instance_id}_tiered_nt{n_tiers}"
    rng0 = np.random.default_rng(derived_seed(instance_id, "tiered_ref"))
    k_ref = tuple(tiered(dag, cpdag, rng0, n_tiers, 0.0))
    g0 = apply_orientations(cpdag, k_ref)

    predictors, reason = score_committed_graph(dag, cpdag, x, y, k_ref, g0)
    if predictors is None:
        return None, reason

    instance = {
        "instance_id": instance_id,
        "base_instance_id": inst.instance_id,
        "arm": "tiered",
        "seed": seed,
        "component_size": component_size,
        "separation": separation,
        "coverage": "",
        "base_wrongness": "",
        "n_tiers": n_tiers,
        "x": x,
        "y": y,
        "dag": dag,
        "cpdag": cpdag,
        "g0": g0,
        "k": k_ref,
        **predictors,
    }
    return instance, REASON_OK


def sample_tiered_state(instance: dict[str, Any], corruption_rate: float, rep: int) -> dict[str, Any]:
    """Draw one tiered-arm sample at a given node-relocation ``corruption_rate``.

    ``d`` (claim-depth) is *measured* post hoc, not targeted: the number of
    orientation claims in ``K_cor`` absent from (including reversed relative
    to) the instance's reference ``K_ref``.

    Args:
        instance: A tiered-arm instance from :func:`build_tiered_instance`.
        corruption_rate: Fraction of nodes relocated to a random other tier.
        rep: Repetition index.

    Returns:
        A row dict matching the unified samples schema.
    """
    seed = derived_seed(instance["instance_id"], "tiered", corruption_rate, rep)
    rng = np.random.default_rng(seed)
    k_cor = tiered(instance["dag"], instance["cpdag"], rng, instance["n_tiers"], corruption_rate)
    ref_set = set(instance["k"])
    cor_set = set(k_cor)
    d_claims = len(cor_set - ref_set)

    g = apply_orientations(instance["cpdag"], k_cor)
    row: dict[str, Any] = {
        "instance_id": instance["instance_id"],
        "arm": "tiered",
        "gate": GATE_NAME,
        "coverage": instance["coverage"],
        "base_wrongness": instance["base_wrongness"],
        "n_tiers": instance["n_tiers"],
        "d": d_claims,
        "rep": rep,
        "flip_rate": "",
        "corruption_rate": corruption_rate,
        "seed": seed,
    }
    if g is None:
        row["status"] = "corrupted_k_contradictory"
        row["survived"] = ""
        row["symdiff_proxy_not_distance"] = ""
    else:
        survived = is_gac_valid_mpdag(g, instance["x"], instance["y"], instance["z_star"])
        row["status"] = "ok"
        row["survived"] = survived
        row["symdiff_proxy_not_distance"] = directed_symdiff(instance["g0"], g)
    return row


# --- curves and AUC ------------------------------------------------------------


def build_curve(samples: list[dict[str, Any]]) -> dict[int, dict[str, Any]]:
    """Bucket sample rows by depth ``d`` into per-depth survival statistics.

    Args:
        samples: Row dicts as produced by :func:`sample_flip_state` /
            :func:`sample_tiered_state`, for one ``(instance, arm)``.

    Returns:
        ``{d: {"n_eval", "n_contradictory", "S", "S_contra_as_fail"}}``. ``S``
        is ``""`` (undefined) when every sample at that depth was
        contradictory -- never a number, never silently dropped.
    """
    by_d: dict[int, list[dict[str, Any]]] = {}
    for row in samples:
        by_d.setdefault(row["d"], []).append(row)

    curve: dict[int, dict[str, Any]] = {}
    for d in sorted(by_d):
        rows = by_d[d]
        n_total = len(rows)
        contra = [r for r in rows if r["status"] == "corrupted_k_contradictory"]
        evaluable = [r for r in rows if r["status"] == "ok"]
        n_contra = len(contra)
        n_eval = len(evaluable)
        s_val: Any = (sum(1 for r in evaluable if r["survived"]) / n_eval) if n_eval else ""
        # Sensitivity: contradictions scored as failures, over ALL samples.
        n_survived_all = sum(1 for r in evaluable if r["survived"])
        s_contra_as_fail = (n_survived_all / n_total) if n_total else ""
        curve[d] = {
            "n_eval": n_eval,
            "n_contradictory": n_contra,
            "S": s_val,
            "S_contra_as_fail": s_contra_as_fail,
        }
    return curve


def auc_frac(curve: dict[int, dict[str, Any]], n_k: int) -> float | str:
    """Primary AUC endpoint: mean ``S`` on the fixed fractional grid.

    For each ``frac in FRAC_GRID``, the nearest depth *actually present* in
    ``curve`` with a defined ``S`` is used (ties broken toward the smaller
    depth) -- this is exact for the flip arm, whose depths are dense
    ``1 .. n_k`` by construction, and is the closest well-defined analogue
    for the tiered arm, whose depths are measured rather than targeted and
    may be sparse.

    Args:
        curve: Output of :func:`build_curve`.
        n_k: The instance's claim count, for the fractional-to-integer map.

    Returns:
        The mean, or ``""`` if no grid point has a defined ``S``.
    """
    defined_depths = sorted(d for d, row in curve.items() if row["S"] != "")
    if not defined_depths or n_k <= 0:
        return ""
    vals = []
    for frac in FRAC_GRID:
        target = frac * n_k
        nearest = min(defined_depths, key=lambda d: (abs(d - target), d))
        vals.append(curve[nearest]["S"])
    return sum(vals) / len(vals)


def auc_abs(curve: dict[int, dict[str, Any]]) -> float | str:
    """Secondary AUC endpoint: raw sum of defined ``S(d)`` over the swept depths."""
    vals = [row["S"] for row in curve.values() if row["S"] != ""]
    if not vals:
        return ""
    return sum(vals)


def fraction_non_monotonic(curves_by_instance: dict[Any, dict[int, dict[str, Any]]]) -> float:
    """Fraction of instances whose ``S(d)`` is not non-increasing in ``d``.

    Monotonicity is measured, never assumed or enforced (acceptance criterion
    A5): this reports what happened, it does not sort or smooth anything.

    Args:
        curves_by_instance: ``{instance_key: curve}`` as built by
            :func:`build_curve`.

    Returns:
        The fraction of instances (with at least 2 defined depths) where some
        ``S(d) > S(d')`` for ``d > d'``.
    """
    if not curves_by_instance:
        return 0.0
    n_checked = 0
    n_violating = 0
    for key in sorted(curves_by_instance, key=str):
        curve = curves_by_instance[key]
        depths = sorted(d for d, row in curve.items() if row["S"] != "")
        if len(depths) < 2:
            continue
        n_checked += 1
        s_vals = [curve[d]["S"] for d in depths]
        if any(s_vals[i] > s_vals[i - 1] for i in range(1, len(s_vals))):
            n_violating += 1
    return (n_violating / n_checked) if n_checked else 0.0
