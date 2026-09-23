"""Matched cross-arm corruption experiment (Axis Robustness, cross-arm design).

Session 7 could not compare correlated (``tiered``) against uniform (``flip``)
knowledge corruption: the two arms were gridded on different units (fraction
of nodes relocated vs. fraction of claims reversed) and, worse, ``tiered`` is
*generative* -- it builds its own ``K``, hence its own ``G0`` and ``Z*`` -- so
the two arms ran on different instance populations entirely. Session 7
explicitly retracted the cross-arm detectability claim over this
(``results/axis_robustness/PREREGISTRATION.md`` Appendix E.3).

This module fixes both defects with one design:

1. **Matched instances.** Every instance is built *once*, the way the tiered
   arm builds it: ``K_ref = tiered(dag, cpdag, rng, n_tiers, 0.0)``,
   ``G0 = apply_orientations(cpdag, K_ref)``, ``Z* =
   optimal_adjustment_set_mpdag(G0, x, y)``. That *same* ``(cpdag, K_ref, G0,
   Z*)`` is then corrupted **both ways** -- ``tiered(dag, cpdag, rng,
   n_tiers, rate)`` for the tiered arm, ``flip(K_ref, rng, rate)`` for the
   flip arm -- so any difference between arms is process, not population.

2. **Unified intensity.** The brief's suggested metric ("percentage of the
   G0 skeleton altered") is identically zero: corruption never changes
   ``cpdag``'s adjacency, only orientations move. Instead, every sampled
   state ``G = apply_orientations(cpdag, K_cor)`` is scored on

   .. code-block:: text

       intensity = |dir(G0) Delta dir(G)| / |dir(G0)|

   the normalized symmetric difference of *directed* edge sets -- defined
   identically for both processes, counting reversals, orientations lost
   (became undirected) and orientations gained.

Known asymmetry (recorded, not hidden): flip's raw symdiff is never below 2
(reversing one claim changes two orientations at minimum -- the claim itself
and, generically, at least one Meek-propagated consequence -- in the typical
case exactly 2), while tiered's is often 0 (a relocated node's crossing
edges can happen to re-derive the same orientations, or the relocation can
land entirely outside the component actually touching X/Y). The two
processes therefore do not cover the intensity axis evenly; this is measured
and reported (J5 / achieved-intensity histograms), never assumed away.

Gate discipline
----------------
The only gate used to decide admissibility is
:func:`bkrobust.benchmarks.measure.fast_gate`; every row carries the literal
string ``"fast_gate"`` in a ``gate`` column.
:func:`bkrobust.synth.component_generator.generate_instance` is used exactly
as the tiered arm of :mod:`bkrobust.robustness.survival` uses it.

RNG discipline
--------------
No global RNG is ever touched. Every stochastic draw goes through an
``np.random.Generator`` built by :func:`derived_seed` (``hashlib.sha256``,
never Python's salted built-in ``hash()``) plus ``np.random.default_rng``.
"""

from __future__ import annotations

import hashlib
import math
from typing import Any, Iterable

import numpy as np

from bkrobust.benchmarks.measure import (
    MAX_G0_UNDIRECTED_FOR_EXTENSIONS,
    O_INTRACTABLE,
    fast_gate,
)
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag
from bkrobust.demo.graph import MPDAG, canon
from bkrobust.demo.meek import apply_orientations
from bkrobust.gac import is_gac_valid_mpdag
from bkrobust.hybrid import HybridResult, breakdown_radius
from bkrobust.synth.component_generator import ComponentSpec, achievable_separations, generate_instance
from bkrobust.synth.knowledge import flip, tiered

Edge = tuple[str, str]

#: Literal gate name stamped on every output row.
GATE_NAME = "fast_gate"

REASON_OK = "ok"

#: The tiered arm's corruption-rate grid (fraction of nodes relocated).
TIERED_RATE_GRID: tuple[float, ...] = tuple(round(0.05 * i, 2) for i in range(1, 11))

#: Intensity bin width for the shared axis, over [0, 1].
BIN_WIDTH = 0.05

#: Minimum non-contradictory samples an (instance, arm, bin) cell needs to
#: count toward AUC_intensity_usable or the paired test.
MIN_USABLE_N = 30

#: Reps per (instance, arm, grid point).
N_REPS_DEFAULT = 200


# --- seeding -----------------------------------------------------------------


def derived_seed(*parts: Any) -> int:
    """A stable, ``PYTHONHASHSEED``-independent seed derived from ``parts``.

    Args:
        *parts: Anything with a stable ``str()``. Order matters.

    Returns:
        A nonnegative integer, suitable for ``np.random.default_rng``.
    """
    payload = "|".join(str(p) for p in parts).encode("utf-8")
    digest = hashlib.sha256(payload).digest()
    return int.from_bytes(digest[:8], "big")


# --- graph-level predictors ----------------------------------------------------


def directed_symdiff(g1: MPDAG, g2: MPDAG) -> int:
    """``|dir(g1) Delta dir(g2)|``, the directed-edge symmetric difference."""
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


def intensity_bin(intensity: float) -> float:
    """Left edge of the ``BIN_WIDTH``-wide bin ``intensity`` falls into.

    No ceiling: intensity can exceed 1.0 (a heavily corrupted state can have
    *more* directed edges disagreeing with ``G0`` than ``G0`` had directed to
    begin with), and that is reported honestly rather than clipped.
    """
    idx = math.floor(round(intensity / BIN_WIDTH, 9))
    return round(idx * BIN_WIDTH, 2)


# --- radius -----------------------------------------------------------------


def flatten_radius(result: HybridResult) -> dict[str, Any]:
    """Flatten a :class:`HybridResult` into CSV-safe columns (see survival.py)."""
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
    elif result.radius == -1:
        row["r_val"] = -1
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


# --- matched instance construction -------------------------------------------


def build_matched_instance(
    component_size: int,
    separation: int,
    seed: int,
    n_tiers: int,
) -> tuple[dict[str, Any] | None, str]:
    """Draw and screen one matched instance, shared by both arms.

    Pipeline: ``generate_instance`` -> require ``.accepted`` -> ``fast_gate``
    -> ``K_ref = tiered(dag, cpdag, rng, n_tiers, 0.0)`` -> ``G0 =
    apply_orientations(cpdag, K_ref)`` -> ``Z* =
    optimal_adjustment_set_mpdag(G0, x, y)`` (reject if undefined) ->
    ``is_gac_valid_mpdag(G0, x, y, Z*)`` (must hold by construction; checked,
    not assumed). ``Z*`` is then held fixed for the whole sweep, on both arms.

    Args:
        component_size: Target undirected-component size.
        separation: Target separation (must be achievable).
        seed: Root seed passed to ``generate_instance``.
        n_tiers: Number of temporal tiers for the reference tiering.

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

    instance_id = f"{inst.instance_id}_xarm_nt{n_tiers}"
    rng0 = np.random.default_rng(derived_seed(instance_id, "k_ref"))
    k_ref = tuple(tiered(dag, cpdag, rng0, n_tiers, 0.0))

    if len(k_ref) == 0:
        return None, "n_k_zero"

    g0 = apply_orientations(cpdag, k_ref)
    if g0 is None:
        # Should not happen: tiered() at corruption_rate=0.0 asserts only
        # true-DAG orientations, which are consistent with cpdag by
        # construction. Reported as its own reason if it ever fires.
        return None, "k_ref_contradictory"

    if len(g0.undirected_edges) > MAX_G0_UNDIRECTED_FOR_EXTENSIONS:
        return None, O_INTRACTABLE

    z_star = optimal_adjustment_set_mpdag(g0, x, y)
    if z_star is None:
        return None, "optimal_set_undefined"
    z_star = frozenset(z_star)

    if not is_gac_valid_mpdag(g0, x, y, z_star):
        return None, "z_invalid_at_g0"

    n_dir_g0 = len(g0.directed_edges)
    if n_dir_g0 == 0:
        # Cannot happen given fast_gate's no_causal_path / structural
        # requirements (X has a directed ancestor path to Y in every
        # accepted instance), but guarded rather than assumed -- intensity's
        # denominator must never be zero.
        return None, "dir_g0_empty"

    instance: dict[str, Any] = {
        "instance_id": instance_id,
        "base_instance_id": inst.instance_id,
        "seed": seed,
        "component_size": component_size,
        "separation": separation,
        "n_tiers": n_tiers,
        "x": x,
        "y": y,
        "dag": dag,
        "cpdag": cpdag,
        "g0": g0,
        "k_ref": k_ref,
        "z_star": z_star,
        "n_k": len(k_ref),
        "shd_truth": directed_symdiff(g0, dag),
        "shd_cpdag": shd_cpdag_count(g0, cpdag),
        "undirected_fraction": undirected_fraction_of(cpdag),
        "z_size": len(z_star),
        "n_dir_g0": n_dir_g0,
    }
    return instance, REASON_OK


# --- sampling (both arms funnel through one scorer) --------------------------


def _score_sample(
    instance: dict[str, Any],
    arm: str,
    k_cor: list[Edge],
    *,
    grid_point: float,
    rep: int,
    seed: int,
    d: Any = "",
) -> dict[str, Any]:
    """Apply ``k_cor`` to ``instance``'s cpdag and score the resulting state."""
    g = apply_orientations(instance["cpdag"], k_cor)
    row: dict[str, Any] = {
        "instance_id": instance["instance_id"],
        "arm": arm,
        "gate": GATE_NAME,
        "component_size": instance["component_size"],
        "separation": instance["separation"],
        "n_tiers": instance["n_tiers"],
        "grid_point": grid_point,
        "flip_rate": grid_point if arm == "flip" else "",
        "corruption_rate": grid_point if arm == "tiered" else "",
        "d": d,
        "rep": rep,
        "seed": seed,
    }
    if g is None:
        row["status"] = "corrupted_k_contradictory"
        row["survived"] = ""
        row["symdiff"] = ""
        row["n_dir_g0"] = instance["n_dir_g0"]
        row["intensity"] = ""
        row["intensity_bin"] = ""
    else:
        symdiff = directed_symdiff(instance["g0"], g)
        intensity = symdiff / instance["n_dir_g0"]
        survived = is_gac_valid_mpdag(g, instance["x"], instance["y"], instance["z_star"])
        row["status"] = "ok"
        row["survived"] = survived
        row["symdiff"] = symdiff
        row["n_dir_g0"] = instance["n_dir_g0"]
        row["intensity"] = intensity
        row["intensity_bin"] = intensity_bin(intensity)
    return row


def sample_tiered_state(instance: dict[str, Any], corruption_rate: float, rep: int) -> dict[str, Any]:
    """Draw one tiered-arm sample: relocate ``corruption_rate`` of nodes, re-derive K."""
    seed = derived_seed(instance["instance_id"], "tiered", corruption_rate, rep)
    rng = np.random.default_rng(seed)
    k_cor = tiered(instance["dag"], instance["cpdag"], rng, instance["n_tiers"], corruption_rate)
    return _score_sample(instance, "tiered", k_cor, grid_point=corruption_rate, rep=rep, seed=seed)


def sample_flip_state(instance: dict[str, Any], d: int, rep: int) -> dict[str, Any]:
    """Draw one flip-arm sample: reverse ``d`` of ``K_ref``'s claims (rate = d/len(K_ref))."""
    n = instance["n_k"]
    assert 1 <= d <= n, f"d={d} out of range for n_k={n}"
    rate = d / n
    seed = derived_seed(instance["instance_id"], "flip", d, rep)
    rng = np.random.default_rng(seed)
    k_cor = flip(list(instance["k_ref"]), rng, rate)
    return _score_sample(instance, "flip", k_cor, grid_point=rate, rep=rep, seed=seed, d=d)


# --- bins ----------------------------------------------------------------------


def build_intensity_bins(samples: list[dict[str, Any]]) -> dict[Any, dict[str, Any]]:
    """Bucket one (instance, arm)'s sample rows into shared intensity bins.

    Contradictions (``apply_orientations`` returned ``None``) carry no
    measured intensity of their own -- ``g`` never existed, so there is
    nothing to diff against ``G0``. Dropping them from the binning would
    silently exclude the very samples that show a process is self-revealing,
    and the design brief asks for ``contradiction_rate`` *per bin*, so they
    cannot simply live in a separate table either.

    Resolution: samples are grouped first by ``grid_point`` (the process
    parameter actually drawn -- ``corruption_rate`` for tiered,
    ``flip_rate`` for flip), each group's **representative intensity** is
    the median measured ``intensity`` of that group's non-contradictory
    draws, and the *entire* group -- contradictions included -- is filed
    under that one shared bin. This is sound because intensity is monotone
    in the process parameter within a fixed instance (more relocation/more
    reversal cannot on average produce a less-oriented disagreement), so a
    grid point's non-contradictory draws are a representative sample of
    "what this grid point does to G0" and stand in for its contradictory
    siblings, which differ only in having crossed into logical
    inconsistency rather than landing at a slightly different orientation
    state.

    A grid point with **zero** non-contradictory draws (all 200 reps
    contradictory) has no representative intensity and is filed under the
    sentinel key ``"unresolved_all_contradictory"`` instead of guessed at --
    reported, never dropped (see ``n_grid_points_unresolved`` at the
    instance level).

    Returns:
        ``{bin_key: {"n_samples", "n_contradictory", "n_noncontra",
        "contradiction_rate", "S", "S_contra_as_fail"}}``. ``bin_key`` is
        either a numeric bin-left (``float``) or the sentinel string above.
        ``S`` is ``""`` (undefined) only when a bin's non-contradictory
        count is zero.
    """
    by_gp: dict[Any, list[dict[str, Any]]] = {}
    for row in samples:
        by_gp.setdefault(row["grid_point"], []).append(row)

    grouped: dict[Any, list[dict[str, Any]]] = {}
    for gp, rows in by_gp.items():
        noncontra = [r for r in rows if r["status"] == "ok"]
        if not noncontra:
            key: Any = "unresolved_all_contradictory"
        else:
            intens = sorted(r["intensity"] for r in noncontra)
            mid = len(intens) // 2
            median = intens[mid] if len(intens) % 2 == 1 else (intens[mid - 1] + intens[mid]) / 2
            key = intensity_bin(median)
        grouped.setdefault(key, []).extend(rows)

    bins: dict[Any, dict[str, Any]] = {}
    for key, rows_here in grouped.items():
        n_total = len(rows_here)
        noncontra = [r for r in rows_here if r["status"] == "ok"]
        n_ok = len(noncontra)
        n_contra = n_total - n_ok
        n_survived = sum(1 for r in noncontra if r["survived"])
        bins[key] = {
            "n_samples": n_total,
            "n_contradictory": n_contra,
            "n_noncontra": n_ok,
            "contradiction_rate": (n_contra / n_total) if n_total else "",
            "S": (n_survived / n_ok) if n_ok else "",
            "S_contra_as_fail": (n_survived / n_total) if n_total else "",
        }
    return bins


def contradiction_rate_overall(samples: list[dict[str, Any]]) -> float:
    """Fraction of all samples (any intensity) that were contradictory."""
    if not samples:
        return 0.0
    n_contra = sum(1 for r in samples if r["status"] == "corrupted_k_contradictory")
    return n_contra / len(samples)


def auc_intensity_usable(bins: dict[float, dict[str, Any]], *, min_n: int = MIN_USABLE_N) -> float | str:
    """Mean ``S`` over intensity bins in ``[0, 1)`` retaining ``>= min_n`` non-contradictory samples."""
    vals = [
        row["S"]
        for b, row in bins.items()
        if isinstance(b, (int, float)) and 0.0 <= b < 1.0 and row["n_noncontra"] >= min_n and row["S"] != ""
    ]
    if not vals:
        return ""
    return sum(vals) / len(vals)


__all__ = [
    "GATE_NAME",
    "REASON_OK",
    "TIERED_RATE_GRID",
    "BIN_WIDTH",
    "MIN_USABLE_N",
    "N_REPS_DEFAULT",
    "derived_seed",
    "directed_symdiff",
    "shd_cpdag_count",
    "undirected_fraction_of",
    "intensity_bin",
    "flatten_radius",
    "compute_r_val",
    "build_matched_instance",
    "sample_tiered_state",
    "sample_flip_state",
    "build_intensity_bins",
    "contradiction_rate_overall",
    "auc_intensity_usable",
    "achievable_separations",
]
