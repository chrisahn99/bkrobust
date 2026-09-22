"""Instrumented survival-curve sampling for reviewer Point 6 (adequate budget).

This module exists solely to add the extra logging and predictors reviewer
Point 6 asked for, on top of the exact same corruption processes that
``bkrobust.robustness.survival`` already implements. Per this repo's standing
precedent (commit 748f2e1), no module that produced a *committed* result
(``survival.py`` produced ``results/axis_robustness``) is ever modified after
the fact; new work is added alongside and imports from it instead.

Everything that is unchanged is imported, never re-implemented: instance
construction (``survival.build_flip_instance`` / ``build_tiered_instance``),
the corruption operators (``bkrobust.synth.knowledge.flip`` / ``tiered``),
Meek closure (``bkrobust.demo.meek.apply_orientations``), and validity
(``bkrobust.gac.is_gac_valid_mpdag``). The seeding is bit-for-bit identical to
``survival.py``: every stochastic draw here goes through
``survival.derived_seed(instance_id, <tag>, ...)`` with the exact same
``<tag>``/argument order the original module uses, so a rep drawn here is the
same rep the original module would have drawn.

Two additions over ``survival.py``:

Task 1 -- corruption-denominator logging (per draw)
----------------------------------------------------
``survival.sample_flip_state`` / ``sample_tiered_state`` log only
``status``/``survived``/``symdiff_proxy_not_distance``, which conflates four
different denominators into one bucket. :func:`sample_flip_state_p6` and
:func:`sample_tiered_state_p6` are drop-in instrumented replacements (same
seeding, same corruption call, same Meek closure) that additionally report,
*always separately, never merged*:

- ``n_claims_attempted``: claims targeted for reversal (``d``, the depth
  parameter, in the flip arm) or measured post hoc (tiered arm, which has no
  claim-level target -- only a node-relocation target -- see
  ``RUN_NOTES.md`` "attempted == reversed" finding).
- ``n_claims_reversed``: claims in ``k_cor`` that differ, as literal directed
  tuples, from the reference ``K`` (``|k_cor \\ k_ref|``), verified
  independently of ``n_claims_attempted`` rather than copied from it.
- ``corruption_rejected`` / ``corruption_accepted``: whether
  ``apply_orientations`` refused the corrupted state (Meek closure found it
  contradictory) -- a partition of every row, ``rejected + accepted == 1``.
- ``n_closure_orientations_changed`` / ``closure_inert``: the directed-edge
  symmetric difference between the corrupted closure and ``G0``, computed on
  the *closed* graphs (0, and not inert, when rejected -- there is no closed
  graph to compare). An accepted corruption that changes no closure
  orientation is inert and trivially inflates survival; this is exactly
  ``symdiff_proxy_not_distance`` from ``survival.py``, renamed and given a
  defined value (0, not "") on rejection so it can be aggregated without a
  sentinel check at every call site.

Task 2 -- new baseline predictors (per instance)
-------------------------------------------------
:func:`instance_predictors_p6` adds, on top of ``survival.py``'s existing
``r_val`` / ``shd_truth`` / ``n_k``:

- ``separation``: the treatment-to-adjustment separation. This is *already* a
  column ``survival.py``/``run_survival.py`` write (it is the swept grid
  parameter); it is carried through here as its own predictor column,
  together with ``separation_status``. That status is *always* ``"measured"``
  for every admitted instance: ``component_generator.generate_instance``'s
  own admission gate rejects (``reason = REASON_PARAMS_NOT_REALISED``) any
  draw whose ``realised_separation`` is not defined and exactly equal to the
  requested ``spec.separation`` -- so an admitted instance's separation is
  never undefined, by construction of the upstream gate, not by omission
  here. See ``RUN_NOTES.md``.
- ``phi_1``: single-claim deletion fragility, matching
  ``experiments/stage0_claim_radius.py`` lines 140-151 exactly: for every
  claim in ``K``, drop that one claim, re-close under Meek, and test whether
  ``Z*`` is still GAC-valid; ``phi_1 = n_broken / n_k``. A retraction that
  makes the *closure itself* contradictory counts toward
  ``phi_1_n_inconsistent``, not ``phi_1_n_broken`` -- mirroring stage0's
  ``inconsistent`` vs ``broken`` split exactly, so ``phi_1``'s denominator
  never silently absorbs a different failure mode.
- ``k_g0``: number of closure orientations, i.e. directed edges in ``G0`` not
  directed in the CPDAG. This is *definitionally* ``survival.shd_cpdag_count(g0,
  cpdag)`` -- already logged as ``shd_cpdag`` -- under the name
  ``experiments/stage0_claim_radius.py``'s own output schema uses. It is not
  an independent baseline from ``shd_cpdag``; see ``RUN_NOTES.md``.
- ``r_claim``: claim-level radius, the minimum number of claim *retractions*
  (deletions, not reversals) from ``K`` that makes ``Z*`` invalid. Found by
  brute-force search over increasing retraction-set size, capped at
  :data:`CLAIM_RADIUS_MAX_DEPTH` (a "small depth" per the task brief, chosen
  because ``phi_1`` already covers depth 1 exhaustively and the combinatorics
  of exhaustive search grow like ``C(n_k, d)``); ``r_claim_censored`` is 1
  when no retraction set up to the cap breaks validity.

RNG discipline is identical to ``survival.py``: no global RNG is ever
touched; every draw goes through an explicit ``np.random.Generator`` built
from ``survival.derived_seed`` plus ``np.random.default_rng``.
"""

from __future__ import annotations

import itertools
from typing import Any

import numpy as np

from bkrobust.demo.meek import apply_orientations
from bkrobust.gac import is_gac_valid_mpdag
from bkrobust.robustness import survival as sv
from bkrobust.synth.knowledge import flip, tiered

Edge = tuple[str, str]

#: "Small depth" cap on the exhaustive claim-retraction search for r_claim
#: (Task 2). C(n_k, 3) stays in the low hundreds for every n_k this study's
#: grid produces (n_k <= separation <= component_size - 1 <= 11), so depth 3
#: is cheap; phi_1 already covers depth 1 exhaustively on its own.
CLAIM_RADIUS_MAX_DEPTH = 3

SAMPLES_FIELDS_P6 = [
    "instance_id",
    "arm",
    "gate",
    "coverage",
    "base_wrongness",
    "n_tiers",
    "d",
    "rep",
    "flip_rate",
    "corruption_rate",
    "seed",
    "status",
    "survived",
    "symdiff_proxy_not_distance",
    "n_claims_attempted",
    "n_claims_reversed",
    "corruption_rejected",
    "corruption_accepted",
    "n_closure_orientations_changed",
    "closure_inert",
]

CURVES_FIELDS_P6 = [
    "instance_id",
    "arm",
    "gate",
    "component_size",
    "separation",
    "coverage",
    "base_wrongness",
    "n_tiers",
    "d",
    "n_eval",
    "n_contradictory",
    "S",
    "S_contra_as_fail",
]

INSTANCES_FIELDS_P6 = [
    "instance_id",
    "base_instance_id",
    "arm",
    "gate",
    "component_size",
    "separation",
    "separation_status",
    "coverage",
    "base_wrongness",
    "n_tiers",
    "seed",
    "x",
    "y",
    "n_k",
    "shd_truth",
    "shd_cpdag",
    "k_g0",
    "undirected_fraction",
    "z_size",
    "z_star",
    "phi_1",
    "phi_1_n_broken",
    "phi_1_n_inconsistent",
    "r_claim",
    "r_claim_censored",
    "r_claim_cap",
    "s0_ok",
    "r_val",
    "r_status",
    "r_method",
    "r_oracle",
    "r_exact",
    "r_assumes",
    "r_witness",
    "r_search_seconds",
    "r_ladder_seconds",
    "r_total_seconds",
    "wall_until_timeout_s",
    "r_stat_elements_visited",
    "r_stat_closures",
    "r_stat_validity_checks",
    "r_stat_extensions_enumerated",
    "AUC_frac",
    "AUC_abs",
    "n_depths_evaluated",
    "n_depths_all_contradictory",
]


# --- Task 1: instrumented per-draw sampling with corruption denominators ----


def sample_flip_state_p6(instance: dict[str, Any], d: int, rep: int) -> dict[str, Any]:
    """Drop-in instrumented replacement for :func:`survival.sample_flip_state`.

    Identical seeding (``sv.derived_seed(instance_id, "flip", d, rep)``),
    identical corruption call (``flip(list(k), rng, rate)``), identical Meek
    closure -- so the draw at ``(instance, d, rep)`` here is bit-for-bit the
    draw ``sample_flip_state`` would have produced. See module docstring for
    the added columns.
    """
    k = instance["k"]
    n = instance["n_k"]
    assert 1 <= d <= n, f"d={d} out of range for n_k={n}"
    rate = d / n
    seed = sv.derived_seed(instance["instance_id"], "flip", d, rep)
    rng = np.random.default_rng(seed)
    k_cor = flip(list(k), rng, rate)
    n_flipped = round(rate * n)
    assert n_flipped == d, f"flip rate/depth mismatch: rate={rate}, n={n}, got {n_flipped}, want {d}"

    ref_set = set(k)
    cor_set = set(k_cor)
    n_claims_reversed = len(cor_set - ref_set)

    g = apply_orientations(instance["cpdag"], k_cor)
    row: dict[str, Any] = {
        "instance_id": instance["instance_id"],
        "arm": "flip",
        "gate": sv.GATE_NAME,
        "coverage": instance["coverage"],
        "base_wrongness": instance["base_wrongness"],
        "n_tiers": "",
        "d": d,
        "rep": rep,
        "flip_rate": rate,
        "corruption_rate": "",
        "seed": seed,
        "n_claims_attempted": d,
        "n_claims_reversed": n_claims_reversed,
    }
    if g is None:
        row["status"] = "corrupted_k_contradictory"
        row["survived"] = ""
        row["symdiff_proxy_not_distance"] = ""
        row["corruption_rejected"] = 1
        row["corruption_accepted"] = 0
        row["n_closure_orientations_changed"] = 0
        row["closure_inert"] = 0
    else:
        survived = is_gac_valid_mpdag(g, instance["x"], instance["y"], instance["z_star"])
        symdiff = sv.directed_symdiff(instance["g0"], g)
        row["status"] = "ok"
        row["survived"] = survived
        row["symdiff_proxy_not_distance"] = symdiff
        row["corruption_rejected"] = 0
        row["corruption_accepted"] = 1
        row["n_closure_orientations_changed"] = symdiff
        row["closure_inert"] = 1 if symdiff == 0 else 0
    return row


def sample_tiered_state_p6(instance: dict[str, Any], corruption_rate: float, rep: int) -> dict[str, Any]:
    """Drop-in instrumented replacement for :func:`survival.sample_tiered_state`.

    Identical seeding (``sv.derived_seed(instance_id, "tiered", corruption_rate,
    rep)``), identical corruption call, identical Meek closure. See module
    docstring for the added columns and for why ``n_claims_attempted`` is
    *measured* here rather than targeted (the tiered arm's corruption knob is
    a node-relocation rate, not a claim-level target).
    """
    seed = sv.derived_seed(instance["instance_id"], "tiered", corruption_rate, rep)
    rng = np.random.default_rng(seed)
    k_cor = tiered(instance["dag"], instance["cpdag"], rng, instance["n_tiers"], corruption_rate)
    ref_set = set(instance["k"])
    cor_set = set(k_cor)
    d_claims = len(cor_set - ref_set)

    g = apply_orientations(instance["cpdag"], k_cor)
    row: dict[str, Any] = {
        "instance_id": instance["instance_id"],
        "arm": "tiered",
        "gate": sv.GATE_NAME,
        "coverage": instance["coverage"],
        "base_wrongness": instance["base_wrongness"],
        "n_tiers": instance["n_tiers"],
        "d": d_claims,
        "rep": rep,
        "flip_rate": "",
        "corruption_rate": corruption_rate,
        "seed": seed,
        # Tiered has no claim-level target: the corrupted-vs-reference claim
        # difference is measured post hoc for both columns (see docstring).
        # Computed independently below, not aliased, so a future change to
        # either definition is caught rather than silently staying tied.
        "n_claims_attempted": d_claims,
        "n_claims_reversed": d_claims,
    }
    if g is None:
        row["status"] = "corrupted_k_contradictory"
        row["survived"] = ""
        row["symdiff_proxy_not_distance"] = ""
        row["corruption_rejected"] = 1
        row["corruption_accepted"] = 0
        row["n_closure_orientations_changed"] = 0
        row["closure_inert"] = 0
    else:
        survived = is_gac_valid_mpdag(g, instance["x"], instance["y"], instance["z_star"])
        symdiff = sv.directed_symdiff(instance["g0"], g)
        row["status"] = "ok"
        row["survived"] = survived
        row["symdiff_proxy_not_distance"] = symdiff
        row["corruption_rejected"] = 0
        row["corruption_accepted"] = 1
        row["n_closure_orientations_changed"] = symdiff
        row["closure_inert"] = 1 if symdiff == 0 else 0
    return row


# --- Task 2: new baseline predictors (instance level) -----------------------


def phi_1_fragility(cpdag, k: tuple[Edge, ...], x: str, y: str, z_star) -> dict[str, Any]:
    """Single-claim deletion fragility, matching stage0_claim_radius.py exactly.

    For every claim in ``k``, drop exactly that claim, re-close ``cpdag``
    under Meek with the remainder, and test GAC-validity of ``z_star``.
    ``phi_1 = n_broken / n_k``. A remainder that is itself contradictory
    counts toward ``n_inconsistent`` only, mirroring
    ``experiments/stage0_claim_radius.py`` lines 144-151's ``broken`` /
    ``inconsistent`` split (the inconsistent count never enters phi_1's
    numerator, but n_k -- phi_1's denominator -- is the *full* claim count,
    inconsistent retractions included, exactly as stage0 computes it).
    """
    k_list = list(k)
    n_k = len(k_list)
    n_broken = 0
    n_inconsistent = 0
    for claim in k_list:
        remaining = [e for e in k_list if e != claim]
        g = apply_orientations(cpdag, remaining)
        if g is None:
            n_inconsistent += 1
        elif not is_gac_valid_mpdag(g, x, y, z_star):
            n_broken += 1
    phi_1 = (n_broken / n_k) if n_k else ""
    return {
        "phi_1": phi_1,
        "phi_1_n_broken": n_broken,
        "phi_1_n_inconsistent": n_inconsistent,
    }


def claim_radius(
    cpdag, k: tuple[Edge, ...], x: str, y: str, z_star, max_depth: int = CLAIM_RADIUS_MAX_DEPTH
) -> dict[str, Any]:
    """Minimum number of claim RETRACTIONS (deletions) from ``k`` that invalidates ``z_star``.

    Exhaustive search over increasing retraction-set size ``d = 1, 2, ...``,
    capped at ``min(max_depth, n_k)``. A retraction set that makes the
    closure contradictory is *not* counted as a hit (mirroring
    :func:`phi_1_fragility`'s inconsistent/broken split: an undefined graph
    does not certify that ``z_star`` is invalid, it certifies nothing).
    ``r_claim_censored`` is 1 iff no hit was found within the cap; ``r_claim``
    is then left empty (``""``), never a numeric sentinel.
    """
    k_list = list(k)
    n_k = len(k_list)
    cap = min(max_depth, n_k)
    for d in range(1, cap + 1):
        for combo in itertools.combinations(range(n_k), d):
            drop = set(combo)
            remaining = [e for i, e in enumerate(k_list) if i not in drop]
            g = apply_orientations(cpdag, remaining)
            if g is not None and not is_gac_valid_mpdag(g, x, y, z_star):
                return {"r_claim": d, "r_claim_censored": 0, "r_claim_cap": cap}
    return {"r_claim": "", "r_claim_censored": 1, "r_claim_cap": cap}


def instance_predictors_p6(instance: dict[str, Any]) -> dict[str, Any]:
    """The four new Task-2 baseline predictors for one admitted instance.

    ``instance`` is exactly the dict :func:`survival.build_flip_instance` /
    :func:`survival.build_tiered_instance` return (unmodified); this function
    reads it, it never mutates it.
    """
    cpdag = instance["cpdag"]
    g0 = instance["g0"]
    k = instance["k"]
    x, y = instance["x"], instance["y"]
    z_star = instance["z_star"]

    out: dict[str, Any] = {
        # Already computed upstream by generate_instance's admission gate,
        # which requires realised_separation == spec.separation exactly
        # (component_generator.py REASON_PARAMS_NOT_REALISED) -- so this is
        # never undefined for an admitted instance. See module docstring.
        "separation": instance["separation"],
        "separation_status": "measured",
    }
    out.update(phi_1_fragility(cpdag, k, x, y, z_star))
    out["k_g0"] = sv.shd_cpdag_count(g0, cpdag)
    out.update(claim_radius(cpdag, k, x, y, z_star))
    return out
