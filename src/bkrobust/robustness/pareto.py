r"""Phase 2: the Pareto frontier of `r_val` (breakdown-radius robustness) against
statistical utility, and against empirical bias under 10% domain-error injection.

See ``results/axis_robustness/PREREGISTRATION.md`` Section 9 (predictions P6,
P7) for the design this module implements. Do not deviate from it without
updating that document.

Live-tree only: every import below resolves to a *live* module
(``bkrobust.demo.*``, ``bkrobust.synth.*``, ``bkrobust.hybrid``,
``bkrobust.benchmarks.*``, ``bkrobust.core.*``). Nothing here imports the dead
top-level scaffold (``bkrobust.graphs``, ``bkrobust.knowledge``,
``bkrobust.metrics``, ``bkrobust.theory``, ``bkrobust.estimation``,
``bkrobust.data``, ``bkrobust.cfm``, ``bkrobust.representations``,
``bkrobust.utils``).

Gate discipline: the only gate ever called directly by this module is
:func:`bkrobust.benchmarks.measure.fast_gate`. The exponential all-extensions
enumerator over in ``bkrobust.demo.evaluate`` and the five-check gate living
in the synth-data ``runner`` submodule are never invoked from here. (The
generator this module calls to build candidate CPDAGs uses the latter
internally as part of its own construction pipeline -- that is the library's
business, not a decision this module makes; every acceptance decision *this
module* makes is re-derived independently from ``fast_gate``.)

Randomness: every draw is seeded from an explicit ``numpy.random.Generator``,
itself derived deterministically (via SHA-256, which is insensitive to the
interpreter's string-hash seed) from small integer/string keys. The global
NumPy or stdlib random modules' seeding entry points are never touched, and
the stdlib random module is never imported.
"""

from __future__ import annotations

import hashlib
import itertools
import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from bkrobust.benchmarks.measure import fast_gate
from bkrobust.demo.evaluate import (
    LinearSEM,
    asymptotic_variance,
    bias,
    optimal_adjustment_set_dag,
    optimal_adjustment_set_mpdag,
    random_sem,
)
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.hybrid import breakdown_radius
from bkrobust.synth.component_generator import (
    MAX_G0_UNDIRECTED_FOR_EXTENSIONS,
    ComponentSpec,
    achievable_separations,
    generate_instance,
)
from bkrobust.synth.knowledge import flip

Edge = tuple[str, str]
Node = str

GATE_NAME = "fast_gate"

#: Cap on proposals sampled per base CPDAG (design doc, Section 9).
PROPOSAL_CAP = 60
#: Independent SEM draws per base CPDAG, shared by every proposal in that
#: stratum, for both the B2 sanity check and the utility measurement.
N_SEM_DRAWS = 20
#: Bias-injection replications per proposal.
N_BIAS_REPS = 50
#: Domain-error injection rate for the bias measurement.
BIAS_FLIP_RATE = 0.10
#: Per-proposal radius budget (design doc, Section 9).
RADIUS_TIME_LIMIT_S = 60.0
#: Same tractability cap the generator itself uses for extension enumeration;
#: reused here for every ``optimal_adjustment_set_mpdag`` call this module
#: makes on a graph it built from a proposal (not from the generator), so a
#: pathological proposal with many surviving undirected edges cannot hang the
#: run. Proposals over this cap are marked ``optimal_set_untractable`` rather
#: than measured.
MAX_UNDIRECTED_FOR_EXTENSIONS = MAX_G0_UNDIRECTED_FOR_EXTENSIONS


def derive_seed(*parts: Any) -> int:
    """A stable, ``PYTHONHASHSEED``-independent seed derived from ``parts``.

    Uses SHA-256 (not Python's built-in ``hash()``, which is randomised per
    process unless ``PYTHONHASHSEED`` is fixed) over a canonical string of the
    parts, so the same ``parts`` always yield the same seed across processes
    and interpreter-hash-seed settings.
    """
    s = "|".join(str(p) for p in parts)
    digest = hashlib.sha256(s.encode("utf-8")).hexdigest()
    return int(digest[:8], 16) % (2**31 - 1)


# ---------------------------------------------------------------------------
# Base CPDAG selection
# ---------------------------------------------------------------------------


@dataclass
class BaseCpdag:
    """One admissible base CPDAG (stratum)."""

    base_id: str
    component_size: int
    separation: int
    seed: int
    treatment: str
    outcome: str
    dag_obj: MPDAG
    cpdag_obj: MPDAG


def pick_separations(component_size: int, k: int = 4) -> list[int]:
    """``k`` roughly evenly spaced achievable separations for ``component_size``.

    Pure arithmetic, no randomness: deterministic regardless of seed or
    ``PYTHONHASHSEED``.
    """
    achievable = list(achievable_separations(component_size))
    if len(achievable) <= k:
        return achievable
    idxs = sorted({round(i * (len(achievable) - 1) / (k - 1)) for i in range(k)})
    return [achievable[i] for i in idxs]


def select_base_cpdags(
    component_sizes: tuple[int, ...] = (6, 8, 10),
    seps_per_size: int = 4,
    accepted_per_cell: int = 2,
    max_attempts_per_cell: int = 10,
    root_seed: int = 0,
    limit: int | None = None,
) -> list[BaseCpdag]:
    """Build the fixed stratum set: base CPDAGs accepted by the generator AND
    ``fast_gate``.

    For each ``(component_size, separation)`` cell, up to ``accepted_per_cell``
    admissible instances are collected by trying successive deterministically
    derived seeds, stopping early once the target is reached. A cell that
    exhausts ``max_attempts_per_cell`` without reaching the target simply
    contributes fewer (possibly zero) base CPDAGs; this is recorded by the
    caller, not silently hidden.
    """
    bases: list[BaseCpdag] = []
    for c in component_sizes:
        for s in pick_separations(c, seps_per_size):
            found = 0
            for attempt in range(max_attempts_per_cell):
                if found >= accepted_per_cell:
                    break
                seed = derive_seed(root_seed, "base", c, s, attempt)
                spec = ComponentSpec(component_size=c, separation=s)
                inst = generate_instance(spec, seed=seed)
                if not inst.accepted:
                    continue
                ok, reason = fast_gate(inst.dag_obj, inst.cpdag_obj, inst.treatment, inst.outcome)
                if not ok:
                    continue
                bases.append(
                    BaseCpdag(
                        base_id=inst.instance_id,
                        component_size=c,
                        separation=s,
                        seed=seed,
                        treatment=inst.treatment,
                        outcome=inst.outcome,
                        dag_obj=inst.dag_obj,
                        cpdag_obj=inst.cpdag_obj,
                    )
                )
                found += 1
            if limit is not None and len(bases) >= limit:
                return bases[:limit]
    return bases


# ---------------------------------------------------------------------------
# Proposal enumeration
# ---------------------------------------------------------------------------


def truthful_orientation_menu(base: BaseCpdag) -> list[Edge]:
    """The undirected edges of ``base.cpdag_obj``, each oriented as in the true DAG.

    Sorted by the edges' canonical (undirected) order, so index ``i`` is stable
    across calls given the same base.
    """
    menu: list[Edge] = []
    for a, b in sorted(base.cpdag_obj.undirected_edges):
        if base.dag_obj.is_directed_edge(a, b):
            menu.append((a, b))
        elif base.dag_obj.is_directed_edge(b, a):
            menu.append((b, a))
        else:
            raise AssertionError(
                f"true DAG does not orient CPDAG-undirected edge {(a, b)} in {base.base_id}"
            )
    return menu


def classify_edge(base: BaseCpdag, edge: Edge, desc_x: set[Node], anc_y: set[Node]) -> bool:
    """Whether ``edge`` (a directed ``tail -> head`` pair) lies on some directed
    path from treatment to outcome in the true DAG.

    ``desc_x``/``anc_y`` are precomputed once per base (``descendants(x) | {x}``,
    ``ancestors(y) | {y}``) and passed in to avoid recomputation per edge.
    """
    tail, head = edge
    return tail in desc_x and head in anc_y


@dataclass
class Proposal:
    proposal_id: str
    base_id: str
    orientations: tuple[Edge, ...]
    size: int
    label: str  # "aggressive" | "conservative"
    n_on_path: int


def sample_proposals(
    base: BaseCpdag,
    menu: list[Edge],
    desc_x: set[Node],
    anc_y: set[Node],
    root_seed: int = 0,
    cap: int = PROPOSAL_CAP,
) -> list[Proposal]:
    """Diverse, deterministic proposals over ``menu``, sizes 1..len(menu).

    Sampling method (see PARETO_RUN_NOTES.md for the same text): the ``cap``
    proposal budget is split evenly across the ``m = len(menu)`` possible
    proposal sizes as ``per_size_target = max(1, cap // m)``. For a size ``k``
    whose ``C(m, k)`` is at or below that target, every subset of that size is
    taken (enumerated via ``itertools.combinations`` over sorted edge indices).
    Otherwise, ``per_size_target`` distinct size-``k`` index subsets are drawn
    by rejection sampling: ``rng.choice(m, size=k, replace=False)``, sorted,
    deduplicated in a set, seeded once per base via
    ``derive_seed(root_seed, "proposals", base.base_id)``. Proposals are then
    returned in ascending ``(size, indices)`` order; if the union across sizes
    still exceeds ``cap`` (rounding), the list is truncated at that point --
    which favours smaller proposal sizes slightly. This is a documented design
    choice, not an attempt to manufacture a trade-off (the *set* of sizes
    sampled, not the direction of any particular pair, is what determines
    whether P6 can show a trade-off).
    """
    m = len(menu)
    if m == 0:
        return []
    rng = np.random.default_rng(derive_seed(root_seed, "proposals", base.base_id))
    per_size_target = max(1, cap // m)
    chosen: dict[tuple[int, ...], list[Edge]] = {}
    for k in range(1, m + 1):
        total_combos = math.comb(m, k)
        if total_combos <= per_size_target:
            combos = list(itertools.combinations(range(m), k))
        else:
            combos_set: set[tuple[int, ...]] = set()
            attempts = 0
            max_attempts = per_size_target * 20 + 100
            while len(combos_set) < per_size_target and attempts < max_attempts:
                attempts += 1
                idx = tuple(sorted(int(i) for i in rng.choice(m, size=k, replace=False)))
                combos_set.add(idx)
            combos = sorted(combos_set)
        for idx in combos:
            chosen[idx] = [menu[i] for i in idx]

    ordered_keys = sorted(chosen.keys(), key=lambda t: (len(t), t))
    if len(ordered_keys) > cap:
        ordered_keys = ordered_keys[:cap]

    proposals: list[Proposal] = []
    for n, idx in enumerate(ordered_keys):
        orientations = tuple(chosen[idx])
        on_path_flags = [classify_edge(base, e, desc_x, anc_y) for e in orientations]
        n_on_path = sum(on_path_flags)
        label = "aggressive" if n_on_path > 0 else "conservative"
        proposals.append(
            Proposal(
                proposal_id=f"{base.base_id}__p{n:04d}_k{len(idx):02d}",
                base_id=base.base_id,
                orientations=orientations,
                size=len(idx),
                label=label,
                n_on_path=n_on_path,
            )
        )
    return proposals


def encode_orientations(orientations: tuple[Edge, ...]) -> str:
    return ";".join(f"{a}->{b}" for a, b in sorted(orientations))


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------


@dataclass
class SemDraw:
    k: int
    sem: LinearSEM


def draw_sems(base: BaseCpdag, root_seed: int, n: int = N_SEM_DRAWS) -> list[SemDraw]:
    """``n`` independent SEMs on the TRUE dag, shared by every proposal in this
    stratum (utility measurement and the B2 sanity check both use these same
    draws)."""
    out = []
    for k in range(n):
        rng = np.random.default_rng(derive_seed(root_seed, "sem", base.base_id, k))
        out.append(SemDraw(k=k, sem=random_sem(base.dag_obj, rng)))
    return out


def b2_sanity_check(base: BaseCpdag, sems: list[SemDraw]) -> list[float]:
    """``bias`` of the true DAG's own optimal adjustment set, for every SEM draw.

    Returns the list of bias values (should all be ~0). Does not raise; the
    caller (``run_pareto.py``) decides what "STOP" means operationally.
    """
    z_true = optimal_adjustment_set_dag(base.dag_obj, base.treatment, base.outcome)
    return [bias(d.sem, base.treatment, base.outcome, z_true) for d in sems]


@dataclass
class ProposalResult:
    row: dict[str, Any]
    sem_rows: list[dict[str, Any]]


def _safe_avar(sem: LinearSEM, x: str, y: str, z: set[str]) -> float | None:
    try:
        v = asymptotic_variance(sem, x, y, z)
    except (FloatingPointError, ZeroDivisionError, np.linalg.LinAlgError):
        return None
    if not math.isfinite(v) or v <= 0:
        return None
    return v


def measure_proposal(
    base: BaseCpdag,
    proposal: Proposal,
    sems: list[SemDraw],
    root_seed: int,
    radius_time_limit_s: float = RADIUS_TIME_LIMIT_S,
    n_bias_reps: int = N_BIAS_REPS,
    bias_flip_rate: float = BIAS_FLIP_RATE,
) -> ProposalResult:
    x, y = base.treatment, base.outcome
    sem_rows: list[dict[str, Any]] = []

    base_row: dict[str, Any] = {
        "base_id": base.base_id,
        "component_size": base.component_size,
        "separation": base.separation,
        "base_seed": base.seed,
        "treatment": x,
        "outcome": y,
        "proposal_id": proposal.proposal_id,
        "proposal_size": proposal.size,
        "label": proposal.label,
        "n_on_path": proposal.n_on_path,
        "orientations": encode_orientations(proposal.orientations),
        "gate": GATE_NAME,
        "status": None,
        "assumes": "",
        "method": "",
        "oracle": "",
        "exact": "",
        "radius_timed_out": "",
        "r_val": "",
        "z_p": "",
        "utility_median": "",
        "utility_q25": "",
        "utility_q75": "",
        "n_sem_valid": 0,
        "n_sem_total": len(sems),
        "mean_abs_bias": "",
        "median_abs_bias": "",
        "n_bias_valid": 0,
        "n_bias_contradictory": 0,
        "n_bias_optimal_undefined": 0,
        "n_bias_untractable": 0,
        "n_bias_total": n_bias_reps,
        "front_utility": False,
        "front_bias": False,
    }

    g0p = apply_orientations(base.cpdag_obj, proposal.orientations)
    if g0p is None:
        base_row["status"] = "proposal_contradictory"
        return ProposalResult(row=base_row, sem_rows=sem_rows)

    if len(g0p.undirected_edges) > MAX_UNDIRECTED_FOR_EXTENSIONS:
        base_row["status"] = "optimal_set_untractable"
        return ProposalResult(row=base_row, sem_rows=sem_rows)

    z_p = optimal_adjustment_set_mpdag(g0p, x, y)
    if z_p is None:
        base_row["status"] = "optimal_set_undefined"
        return ProposalResult(row=base_row, sem_rows=sem_rows)

    base_row["status"] = "ok"
    base_row["z_p"] = ";".join(sorted(z_p))

    # --- r_val ---
    result = breakdown_radius(
        base.cpdag_obj,
        list(proposal.orientations),
        x,
        y,
        frozenset(z_p),
        g0=g0p,
        time_limit_s=radius_time_limit_s,
    )
    base_row["assumes"] = result.assumes
    base_row["method"] = result.method
    base_row["oracle"] = result.oracle
    base_row["exact"] = bool(result.exact)
    base_row["radius_timed_out"] = not result.exact
    base_row["r_val"] = int(result.radius)

    # --- utility: 20 independent SEM draws, shared across proposals ---
    utilities: list[float] = []
    for d in sems:
        avar = _safe_avar(d.sem, x, y, z_p)
        status = "ok" if avar is not None else "avar_undefined"
        utility = (1.0 / avar) if avar is not None else None
        sem_rows.append(
            {
                "base_id": base.base_id,
                "proposal_id": proposal.proposal_id,
                "sem_k": d.k,
                "status": status,
                "avar": avar if avar is not None else "",
                "utility": utility if utility is not None else "",
            }
        )
        if utility is not None:
            utilities.append(utility)
    base_row["n_sem_valid"] = len(utilities)
    if utilities:
        arr = np.array(utilities, dtype=float)
        base_row["utility_median"] = float(np.median(arr))
        base_row["utility_q25"] = float(np.percentile(arr, 25))
        base_row["utility_q75"] = float(np.percentile(arr, 75))

    # --- bias under 10% domain-error injection ---
    # Design decision (documented in PARETO_RUN_NOTES.md): all reps use the
    # SAME sem (sems[0]) so the 50 reps isolate variation from which Z the
    # corrupted knowledge yields, not from re-drawn SEM weights. The 20-draw
    # requirement in the design doc is specific to the utility measurement.
    bias_sem = sems[0].sem
    abs_biases: list[float] = []
    n_contra = 0
    n_optundef = 0
    n_untract = 0
    for rep in range(n_bias_reps):
        rng = np.random.default_rng(
            derive_seed(root_seed, "bias", proposal.proposal_id, rep)
        )
        k_cor = flip(list(proposal.orientations), rng, bias_flip_rate)
        g0e = apply_orientations(base.cpdag_obj, k_cor)
        if g0e is None:
            n_contra += 1
            continue
        if len(g0e.undirected_edges) > MAX_UNDIRECTED_FOR_EXTENSIONS:
            n_untract += 1
            continue
        z_e = optimal_adjustment_set_mpdag(g0e, x, y)
        if z_e is None:
            n_optundef += 1
            continue
        b = bias(bias_sem, x, y, z_e)
        abs_biases.append(abs(b))
    base_row["n_bias_valid"] = len(abs_biases)
    base_row["n_bias_contradictory"] = n_contra
    base_row["n_bias_optimal_undefined"] = n_optundef
    base_row["n_bias_untractable"] = n_untract
    if abs_biases:
        arr = np.array(abs_biases, dtype=float)
        base_row["mean_abs_bias"] = float(np.mean(arr))
        base_row["median_abs_bias"] = float(np.median(arr))

    return ProposalResult(row=base_row, sem_rows=sem_rows)


# ---------------------------------------------------------------------------
# Pareto fronts
# ---------------------------------------------------------------------------


def _r_rank(r_val: int) -> float:
    """Sortable rank for a radius under the UNREACHED convention: UNREACHED is
    the best possible outcome (no failure found in the enumerated space), so it
    ranks above every finite radius. Never used as the literal numeric value in
    any output column or average."""
    return math.inf if r_val == -1 else float(r_val)


def _non_dominated(points: list[tuple[str, float, float]]) -> set[str]:
    """Non-dominated ids under "maximize both coordinates"."""
    ids = [p[0] for p in points]
    front: set[str] = set()
    for i, (pid, ra, ua) in enumerate(points):
        dominated = False
        for j, (pid2, rb, ub) in enumerate(points):
            if i == j:
                continue
            if rb >= ra and ub >= ua and (rb > ra or ub > ua):
                dominated = True
                break
        if not dominated:
            front.add(pid)
    return front


def compute_fronts(rows: list[dict[str, Any]]) -> None:
    """Mutates ``rows`` in place, setting ``front_utility``/``front_bias`` for
    the "ok"-status rows of a single base/stratum. ``rows`` must all share one
    ``base_id``."""
    util_points = []
    bias_points = []
    for r in rows:
        if r["status"] != "ok":
            continue
        rr = _r_rank(r["r_val"])
        if r["utility_median"] != "":
            util_points.append((r["proposal_id"], rr, float(r["utility_median"])))
        if r["mean_abs_bias"] != "":
            bias_points.append((r["proposal_id"], rr, -float(r["mean_abs_bias"])))

    util_front = _non_dominated(util_points)
    bias_front = _non_dominated(bias_points)
    for r in rows:
        if r["proposal_id"] in util_front:
            r["front_utility"] = True
        if r["proposal_id"] in bias_front:
            r["front_bias"] = True
