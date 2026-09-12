"""Axis Robustness, dual-axis check (P5): does `r_val` (BFS hops) agree with
the survival sweep's native depth axis (`d_claims`)?

`r_val` (``bkrobust.hybrid.breakdown_radius(...).radius``) is defined as
**minimum BFS hop distance** on the covering-relation neighbour graph
(``bkrobust.core.spacelib``). The knowledge-corruption survival sweep in
:mod:`bkrobust.robustness.survival` corrupts *knowledge claims*, and reports
depth in claim units, ``d_claims``. These are different units and the map
between them is not the identity -- see
``results/axis_robustness/PREREGISTRATION.md`` §2 and prediction P5 (§6).

No cheap single-pair distance oracle exists for this design:
``search/exact_fast.py::radius_local_up_fast`` searches **upward** only, and a
flip-corrupted state is generally not a refinement of ``G0`` (neither is above
the other in the model-inclusion order), so it cannot answer "how many hops
from G0 to this specific corrupted state." The only exact route to `d_hops` is
full space enumeration: ``bkrobust.core.spacelib.build_space`` (exponential in
the CPDAG's undirected-edge count) followed by a single BFS
(``distances_from``) from ``G0``. This module does that, on a small-instance
subsample sized by a cost census (see :func:`census_one_size`), and reports
whether the two axes' rank correlations with `r_val` agree in sign.

Construction reuse
-------------------
This module is a read-only consumer of
:mod:`bkrobust.robustness.survival`'s already-vetted flip-arm construction
(``build_flip_instance``, ``score_committed_graph``, ``compute_r_val``,
``build_curve``, ``auc_frac``, ``directed_symdiff``, ``derived_seed``,
``separations_for``, ``GATE_NAME``) rather than re-implementing any of it.
Nothing here edits ``survival.py``.

Gate discipline
----------------
The only gate is :func:`bkrobust.benchmarks.measure.fast_gate`, exclusively
via ``survival.build_flip_instance``, which calls it internally. Every sample
row and instance row carries ``gate == "fast_gate"``
(:data:`bkrobust.robustness.survival.GATE_NAME`).
The exponential all-subsets adjustment-set enumerator, and the runner-level
gate wrapping it (see the hard constraints in the task brief this module
answers), are never called, directly or transitively, from this module.

RNG discipline
--------------
No global RNG is ever touched. Every stochastic draw goes through an
``np.random.Generator`` built from ``survival.derived_seed`` (SHA-256 based,
``PYTHONHASHSEED``-independent) plus ``np.random.default_rng``.
"""

from __future__ import annotations

import signal
import time
from typing import Any, Iterable

import numpy as np

from bkrobust.core.spacelib import Space, build_space, distances_from
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.gac import is_gac_valid_mpdag
from bkrobust.robustness import survival as sv
from bkrobust.synth.component_generator import (
    ComponentSpec,
    achievable_separations,
    generate_instance,
)
from bkrobust.synth.knowledge import flip

Edge = tuple[str, str]

#: Reused verbatim -- the only gate this module ever stamps.
GATE_NAME = sv.GATE_NAME

#: Fixed fractional grid, same shape as PREREGISTRATION.md §4's AUC_frac,
#: reused unmodified via ``survival.auc_frac`` for both axes.
FRAC_GRID = sv.FRAC_GRID


# --- Step 1: cost census -----------------------------------------------------


class _CensusTimeout(Exception):
    """Raised by the SIGALRM handler when ``build_space`` exceeds its budget."""


def _alarm_handler(signum: int, frame: Any) -> None:  # pragma: no cover - trivial
    raise _CensusTimeout()


def census_one_size(
    component_size: int,
    *,
    seed_budget: int = 200,
    triangle_prob: float = 0.3,
    timeout_s: float = 30.0,
) -> dict[str, Any]:
    """Time ``build_space`` on one accepted instance at ``component_size``.

    Finds the first accepted instance at a representative (median achievable)
    separation, then times the exponential cost centre exactly once. Guarded
    by a ``SIGALRM``-based wall-clock timeout so a runaway enumeration cannot
    hang the whole census; POSIX-only, fine on darwin/linux, single-threaded
    (no interaction with any RNG).

    Args:
        component_size: Target undirected-component size.
        seed_budget: Seeds tried (via ``generate_instance``) to find one
            accepted instance.
        triangle_prob: Passed to :class:`ComponentSpec`; default matches the
            survival sweep and every other user of this generator.
        timeout_s: Wall-clock cap on the ``build_space`` call itself.

    Returns:
        A row: ``component_size``, ``separation``, ``coverage``,
        ``triangle_prob``, ``seeds_tried``, ``n_undirected_cpdag``,
        ``n_space_elements``, ``build_seconds``, ``status`` (``"ok"``,
        ``"no_accepted_instance_found"``, ``f"timeout_{timeout_s}s"``, or
        ``f"error:{type}"``).
    """
    seps = achievable_separations(component_size)
    separation = seps[len(seps) // 2]
    spec = ComponentSpec(
        component_size=component_size, separation=separation, coverage=1.0, triangle_prob=triangle_prob
    )
    inst = None
    tried = 0
    for seed in range(seed_budget):
        tried += 1
        cand = generate_instance(spec, seed)
        if cand.accepted:
            inst = cand
            break

    row: dict[str, Any] = {
        "component_size": component_size,
        "separation": separation,
        "coverage": 1.0,
        "triangle_prob": triangle_prob,
        "seeds_tried": tried,
        "n_undirected_cpdag": "",
        "n_space_elements": "",
        "build_seconds": "",
        "status": "",
    }
    if inst is None:
        row["status"] = "no_accepted_instance_found"
        return row

    assert inst.cpdag_obj is not None
    row["n_undirected_cpdag"] = len(inst.cpdag_obj.undirected_edges)

    old_handler = signal.signal(signal.SIGALRM, _alarm_handler)
    signal.setitimer(signal.ITIMER_REAL, timeout_s)
    t0 = time.perf_counter()
    try:
        space = build_space(inst.cpdag_obj)
        row["build_seconds"] = time.perf_counter() - t0
        row["n_space_elements"] = len(space)
        row["status"] = "ok"
    except _CensusTimeout:
        row["build_seconds"] = timeout_s
        row["status"] = f"timeout_{timeout_s}s"
    except Exception as exc:  # noqa: BLE001 - census must never crash the run
        row["build_seconds"] = time.perf_counter() - t0
        row["status"] = f"error:{type(exc).__name__}"
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old_handler)
    return row


# --- Step 3: dual-axis sampling ----------------------------------------------


def build_space_for_instance(cpdag: MPDAG) -> tuple[Space, float]:
    """``build_space`` timed, once per instance. Not a distance oracle shortcut."""
    t0 = time.perf_counter()
    space = build_space(cpdag)
    return space, time.perf_counter() - t0


def sample_dual_axis_state(
    instance: dict[str, Any],
    d_claims: int,
    rep: int,
    *,
    space_elements: frozenset[MPDAG],
    dist: dict[MPDAG, int],
) -> dict[str, Any]:
    """Draw one flip-arm sample at targeted claim-depth ``d_claims`` and score
    it on *both* axes: claim-depth (targeted) and hop-depth (measured from
    ``dist``, the BFS table computed once per instance).

    Mirrors ``survival.sample_flip_state``'s corruption draw exactly (same
    ``flip`` call, same rate), but additionally looks the resulting state up
    in the pre-built ``Space`` and its precomputed distance table, and
    verifies space membership per hard-constraint discipline: an undefined
    quantity is an empty string plus a status, never a number.

    Args:
        instance: A flip-arm instance from ``survival.build_flip_instance``.
        d_claims: Targeted number of reversed claims, ``1 <= d_claims <= n_k``.
        rep: Repetition index.
        space_elements: ``frozenset(space.elements)`` for O(1) membership.
        dist: ``distances_from(space, instance["g0"])``.

    Returns:
        A row dict. ``status`` is one of ``"ok"`` (corrupted state exists, is
        an element of the space, and was reached by the BFS -- ``d_hops``
        defined), ``"unreachable_in_space"`` (an element of the space but not
        reached by the BFS from ``G0`` -- ``d_hops`` empty),
        ``"not_in_space_BUG"`` (the corrupted MPDAG is not even an element of
        the fully enumerated space -- should be impossible per
        ``spacelib``'s docstring; flagged loudly, not silently dropped), or
        ``"corrupted_k_contradictory"`` (``apply_orientations`` returned
        ``None``). ``gac_status`` is ``"ok"``/``"corrupted_k_contradictory"``
        only -- the coarser split that ``survival.build_curve`` expects for
        the claims-axis curve, since GAC survival is well-defined regardless
        of whether the state was BFS-reachable.
    """
    k = instance["k"]
    n = instance["n_k"]
    assert 1 <= d_claims <= n, f"d_claims={d_claims} out of range for n_k={n}"
    rate = d_claims / n
    seed = sv.derived_seed(instance["instance_id"], "hops_flip", d_claims, rep)
    rng = np.random.default_rng(seed)
    k_cor = flip(list(k), rng, rate)
    n_flipped = round(rate * n)
    assert n_flipped == d_claims, (
        f"flip rate/depth mismatch: rate={rate}, n={n}, got {n_flipped}, want {d_claims}"
    )

    g = apply_orientations(instance["cpdag"], k_cor)
    row: dict[str, Any] = {
        "instance_id": instance["instance_id"],
        "component_size": instance["component_size"],
        "separation": instance["separation"],
        "coverage": instance["coverage"],
        "base_wrongness": instance["base_wrongness"],
        "gate": GATE_NAME,
        "d_claims": d_claims,
        "rep": rep,
        "seed": seed,
    }
    if g is None:
        row["status"] = "corrupted_k_contradictory"
        row["gac_status"] = "corrupted_k_contradictory"
        row["in_space"] = ""
        row["d_hops"] = ""
        row["symdiff_proxy_not_distance"] = ""
        row["survived"] = ""
        row["g_equals_g0"] = ""
        return row

    g0 = instance["g0"]
    survived = is_gac_valid_mpdag(g, instance["x"], instance["y"], instance["z_star"])
    symdiff = sv.directed_symdiff(g0, g)
    g_equals_g0 = g == g0
    in_space = g in space_elements

    row["gac_status"] = "ok"
    row["survived"] = survived
    row["symdiff_proxy_not_distance"] = symdiff
    row["g_equals_g0"] = g_equals_g0
    row["in_space"] = in_space

    if not in_space:
        # spacelib's docstring guarantees every state apply_orientations can
        # produce from cpdag is a model-inclusion refinement, hence an
        # element of build_space(cpdag)'s enumeration. If this ever fires it
        # is a bug in enumerate_space / covering_pairs, not in this sampler.
        row["status"] = "not_in_space_BUG"
        row["d_hops"] = ""
        return row

    dh = dist.get(g)
    if dh is None:
        row["status"] = "unreachable_in_space"
        row["d_hops"] = ""
    else:
        row["status"] = "ok"
        row["d_hops"] = dh
    return row


# --- Step 4: curves and the two AUCs ------------------------------------------


def claims_curve(rows: list[dict[str, Any]]) -> dict[int, dict[str, Any]]:
    """Survival curve on the ``d_claims`` axis, via ``survival.build_curve`` verbatim.

    Every non-contradictory row counts as evaluable (``gac_status == "ok"``),
    regardless of whether it was also BFS-reachable -- GAC survival is
    well-defined at any non-contradictory corrupted state.
    """
    translated = [
        {"d": r["d_claims"], "status": r["gac_status"], "survived": r["survived"]} for r in rows
    ]
    return sv.build_curve(translated)


def hops_curve(rows: list[dict[str, Any]]) -> dict[int, dict[str, Any]]:
    """Survival curve on the ``d_hops`` axis, via ``survival.build_curve`` verbatim.

    Only rows with ``status == "ok"`` (space element, BFS-reachable, ``d_hops``
    defined) contribute; contradictions, unreachable states, and the
    not-in-space bug sentinel are excluded here (their exclusion counts are
    reported separately, never silently folded into this curve).
    """
    translated = [
        {"d": r["d_hops"], "status": "ok", "survived": r["survived"]}
        for r in rows
        if r["status"] == "ok"
    ]
    return sv.build_curve(translated)


def auc_hops(curve: dict[int, dict[str, Any]]) -> tuple[float | str, int | str]:
    """AUC on the ``d_hops`` axis, normalized by the max finite ``d_hops`` observed.

    Same fractional-grid-with-nearest-defined-depth construction as
    ``survival.auc_frac`` (reused verbatim), but the normalizer is the
    instance's own observed ``h_max`` rather than ``n_k`` -- this is a
    distance-unit analogue of ``AUC_frac``, not the same quantity.

    Args:
        curve: Output of :func:`hops_curve`.

    Returns:
        ``(AUC_hops, h_max)``; both ``""`` if no depth in ``curve`` has a
        defined ``S`` (e.g. every sample for this instance was contradictory
        or unreachable).
    """
    defined_depths = [d for d, row in curve.items() if row["S"] != ""]
    if not defined_depths:
        return "", ""
    h_max = max(defined_depths)
    return sv.auc_frac(curve, h_max), h_max


# --- symdiff-vs-true-distance gap ---------------------------------------------


def symdiff_gap_stats(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """How often the symmetric-difference proxy equals / exceeds the true ``d_hops``.

    THEOREMS.md §10-11 give only ``r <= r_E2 <= r_E1``: the directed-edge
    symmetric difference is a documented **upper-bound surrogate** for
    distance, not distance itself (Lemma R's converse is verified, not
    proved). This is the first direct per-sample measurement of that gap
    project-wide.

    Args:
        rows: Sample rows with ``status == "ok"`` (``d_hops`` defined).

    Returns:
        Counts and fractions: ``n_symdiff_compared``, ``n_symdiff_eq_dhops``,
        ``n_symdiff_gt_dhops``, ``n_symdiff_lt_dhops`` (should be exactly 0;
        a nonzero count here would falsify the documented upper bound and
        must be reported loudly, not silently dropped), plus the ``frac_*``
        versions (``""`` if ``n_symdiff_compared == 0``).
    """
    ok_rows = [r for r in rows if r["status"] == "ok"]
    n = len(ok_rows)
    n_eq = sum(1 for r in ok_rows if r["symdiff_proxy_not_distance"] == r["d_hops"])
    n_gt = sum(1 for r in ok_rows if r["symdiff_proxy_not_distance"] > r["d_hops"])
    n_lt = sum(1 for r in ok_rows if r["symdiff_proxy_not_distance"] < r["d_hops"])
    return {
        "n_symdiff_compared": n,
        "n_symdiff_eq_dhops": n_eq,
        "n_symdiff_gt_dhops": n_gt,
        "n_symdiff_lt_dhops": n_lt,
        "frac_symdiff_eq_dhops": (n_eq / n) if n else "",
        "frac_symdiff_gt_dhops": (n_gt / n) if n else "",
        "frac_symdiff_lt_dhops": (n_lt / n) if n else "",
    }


def dclaims_dhops_relationship(rows: Iterable[dict[str, Any]]) -> dict[int, dict[str, Any]]:
    """Mean/sd of ``d_hops`` per ``d_claims``, over ``status == "ok"`` rows only.

    Args:
        rows: Sample rows (any instance/cell mix; caller decides pooling).

    Returns:
        ``{d_claims: {"n", "mean_d_hops", "sd_d_hops"}}``, sorted by
        ``d_claims``. Empty buckets are omitted, never zero-filled.
    """
    by_d: dict[int, list[int]] = {}
    for r in rows:
        if r["status"] != "ok":
            continue
        by_d.setdefault(r["d_claims"], []).append(r["d_hops"])
    out: dict[int, dict[str, Any]] = {}
    for d in sorted(by_d):
        vals = by_d[d]
        n = len(vals)
        mean = sum(vals) / n
        var = sum((v - mean) ** 2 for v in vals) / n if n > 1 else 0.0
        out[d] = {"n": n, "mean_d_hops": mean, "sd_d_hops": var**0.5}
    return out
