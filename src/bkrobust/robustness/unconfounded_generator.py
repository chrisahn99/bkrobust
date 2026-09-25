r"""Instance/query generation that does NOT target separation (v2 unconfounded arm).

Why this module exists
-----------------------
``bkrobust.synth.component_generator.generate_instance`` takes ``separation``
as an *input* (``ComponentSpec.separation``): the builder lays a spine of
exactly ``separation`` edges from the anchor to the treatment ``X`` before
growing the rest of the component, and the admission gate
(``REASON_PARAMS_NOT_REALISED``) *rejects* any draw whose realised separation
does not exactly equal the requested one. Every committed survival sweep
(``results/axis_robustness``, ``results/axis_robustness_p6``) draws instances
at ``separation in {min, median, max achievable}`` -- separation is the
generator's own **design variable**, targeted by construction, not merely
measured. A reviewer can therefore say the headline "r_val beats separation"
result is comparing the breakdown radius against the very quantity the
sampler was built around, on a population selected to realise chosen values
of that quantity. This module produces a population where that is not true:
separation is measured, never targeted, and the treatment is chosen only
*after* the graph exists.

The construction (kept as close to ``component_generator`` as the deviation
demands, and no closer)
------------------------------------------------------------------------------
1. **Grow the component exactly as :func:`component_generator.build_component_dag`
   grows its non-spine vertices** -- the very same clique-attachment loop
   (vertex ``j`` attaches to a single existing vertex, or with probability
   ``triangle_prob`` to an existing *edge*, drawn from ``0..j-1``), but
   starting at ``j = 1`` instead of ``j = separation + 1``: there is no
   preallocated spine here at all, so nothing during graph construction ever
   reads a ``separation`` value. Vertex ``0`` is not special during
   construction; it only later plays the role ``component_generator`` calls
   "anchor" (see step 2). Edges are still oriented earlier -> later in
   construction order, for the identical reason ``component_generator``
   orients them that way: this keeps every vertex's parent set a clique (no
   v-structure inside the component, so the CPDAG comes back fully
   undirected) and, by a short induction on construction order (vertex ``j``
   only ever attaches to vertices in ``0..j-1``), makes vertex ``0`` an
   ancestor of *every* other component vertex -- exactly the "anchor is an
   ancestor of X" property the spine construction existed to guarantee, but
   obtained here for free, for *every* vertex, not just one chosen at
   construction time.
2. Attach the fixed outcome ``Y`` and spectator ``W`` exactly as
   ``component_generator`` does: ``anchor -> Y``, ``candidate -> Y`` for
   whichever vertex ends up playing treatment, ``W -> Y``, ``W`` adjacent to
   nothing else. ``anchor`` is fixed to be construction-vertex ``0`` (the
   role ``component_generator`` also assigns to construction-vertex ``0``);
   this is a structural role, not part of "the query" -- see point 3.
3. **Admissible-query enumeration.** For every other component vertex
   ``candidate`` (``c - 1`` of them), build the one-candidate-at-a-time DAG
   with ``candidate`` playing treatment, take its CPDAG, and test the *same*
   admissibility rule ``survival.build_tiered_instance`` itself gates the
   final instance on -- ``bkrobust.benchmarks.measure.fast_gate`` (treatment
   touches an undirected component; outcome is a genuine descendant; a valid
   adjustment set exists in the truth; the empty set is not trivially valid;
   some atomic knowledge-perturbation changes validity -- same verdict
   vocabulary as ``synth.runner.gate``, without its exhaustive enumeration
   of every valid adjustment set, which does not scale to being run once per
   candidate per seed; see :func:`admissible_queries`'s docstring), then the
   *truthful, full-coverage* committed graph (every component edge oriented
   to its true direction) is required to be tractable
   (``<= MAX_G0_UNDIRECTED_FOR_EXTENSIONS`` undirected edges) and to have an
   identified, valid optimal adjustment set -- the same tractability/
   identification check ``component_generator.generate_instance`` makes. The
   admissible set is collected in full (never short-circuited at the first
   hit) and the treatment is then drawn **uniformly at random** from it via
   the instance's own RNG: the gate decides *which* vertices are legal
   queries, but never *which one* is asked; that choice is uniform, not
   ranked by, or biased toward, any resulting separation.
4. Downstream of the chosen treatment, this module is a thin wrapper: it
   calls ``bkrobust.benchmarks.measure.fast_gate`` and
   ``bkrobust.robustness.survival.score_committed_graph`` -- imported, never
   reimplemented -- exactly as
   ``bkrobust.robustness.survival.build_tiered_instance`` does, using
   ``bkrobust.synth.knowledge.tiered`` for the analyst's committed knowledge.
   Nothing about corruption sampling, curve-building, ``r_val``, or the other
   baseline predictors changes; :mod:`bkrobust.robustness.survival_p6`'s
   ``sample_tiered_state_p6`` / ``instance_predictors_p6`` are reused
   verbatim by the driver (``run_unconfounded_v2.py``) against the instances
   this module returns.

Separation is measured, and *kept when undefined*
---------------------------------------------------
``realised_separation`` here is the BFS distance, inside the **CPDAG's**
undirected component containing the treatment, from the treatment to the
nearest member of the *tiered arm's own* ``Z* = O(G0_tiered)`` -- not the
truthful-coverage optimal set used only for step 3's admissibility check.
Because ``Z*`` is read off the *tiered* committed graph (which need not equal
the fully-truthful one), it can fail to contain ``anchor``, or contain no
in-component member at all, in which case separation is genuinely undefined
(``separation_status != "measured"``). Per this study's brief, such rows are
**kept**, with ``separation`` empty and ``separation_status`` recording why --
never dropped, and never backfilled with a sentinel number. These are exactly
the rows where ``r_val`` is the only defined structural predictor, which is
the population this module exists to expose.

Determinism: the sole source of randomness for one instance is
``np.random.default_rng(seed)``, threaded through graph construction *and*
the treatment draw (step 3) *and* on into ``tiered()`` via
``survival.derived_seed`` exactly as ``build_tiered_instance`` does -- no
global RNG is ever touched.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any

import numpy as np

from bkrobust.benchmarks.measure import fast_gate
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG, canon
from bkrobust.demo.meek import apply_orientations
from bkrobust.robustness import survival as sv
from bkrobust.synth.component_generator import (
    MAX_G0_UNDIRECTED_FOR_EXTENSIONS,
    OUTCOME,
    SPECTATOR,
    component_containing,
    distances_within_component,
)
from bkrobust.synth.knowledge import tiered

Edge = tuple[str, str]

#: Rejection reasons, namespaced so exclusions can be counted by cause,
#: mirroring ``component_generator``'s / ``survival.build_tiered_instance``'s
#: own vocabulary wherever the same check is being made.
REASON_NO_ADMISSIBLE_QUERY = "no_admissible_query"
REASON_OK = "ok"


def _component_labels(component_size: int) -> tuple[str, ...]:
    """Generic component labels, construction order. No vertex is privileged
    as "X" here -- unlike ``component_generator._component_labels``, which
    bakes the treatment's identity into the label at build time -- because
    which vertex plays treatment is decided *after* the graph exists (step 3
    of the module docstring)."""
    width = max(2, len(str(component_size - 1)))
    return tuple(f"C{j:0{width}d}" for j in range(component_size))


@dataclass(frozen=True)
class UnconfoundedLayout:
    labels: tuple[str, ...]
    anchor: str
    component_edges: tuple[tuple[str, str], ...]


def build_unconfounded_component(
    component_size: int, triangle_prob: float, rng: np.random.Generator
) -> UnconfoundedLayout:
    """Grow a connected chordal component with no preallocated spine.

    Identical clique-attachment loop to
    ``component_generator.build_component_dag``'s non-spine growth step, run
    over the *whole* component (``j = 1 .. component_size - 1``) instead of
    only its tail. See module docstring point 1 for why vertex 0 (returned as
    ``anchor``) ends up an ancestor of every other vertex despite that.
    """
    labels = _component_labels(component_size)
    adjacency: dict[int, set[int]] = {j: set() for j in range(component_size)}
    index_edges: list[tuple[int, int]] = []

    def connect(earlier: int, later: int) -> None:
        index_edges.append((earlier, later))
        adjacency[earlier].add(later)
        adjacency[later].add(earlier)

    for j in range(1, component_size):
        anchor_idx = int(rng.integers(0, j))
        attach_to = [anchor_idx]
        neighbours = sorted(adjacency[anchor_idx])
        if neighbours and rng.uniform() < triangle_prob:
            attach_to.append(neighbours[int(rng.integers(0, len(neighbours)))])
        for base in sorted(attach_to):
            connect(base, j)

    component_edges = tuple((labels[a], labels[b]) for a, b in index_edges)
    return UnconfoundedLayout(labels=labels, anchor=labels[0], component_edges=component_edges)


def _build_dag_for_candidate(layout: UnconfoundedLayout, candidate: str) -> MPDAG:
    """The full DAG (component + Y + W) with ``candidate`` playing treatment."""
    nodes = [*layout.labels, OUTCOME, SPECTATOR]
    directed = [
        *layout.component_edges,
        (candidate, OUTCOME),
        (layout.anchor, OUTCOME),
        (SPECTATOR, OUTCOME),
    ]
    return MPDAG(nodes, directed=directed)


def _truthful_g0_identified_and_valid(dag: MPDAG, cpdag: MPDAG, candidate: str) -> bool:
    """Coverage-1.0 admissibility check: same predicate
    ``component_generator.generate_instance`` applies (tractability +
    identification + validity of the fully-truthful committed graph), minus
    the separation-realisation clause.
    """
    undirected_keys = {canon(a, b) for a, b in cpdag.undirected_edges}
    k_full = tuple(sorted((a, b) for (a, b) in dag.directed_edges if canon(a, b) in undirected_keys))
    g0 = apply_orientations(cpdag, k_full)
    if g0 is None:
        return False
    if len(g0.undirected_edges) > MAX_G0_UNDIRECTED_FOR_EXTENSIONS:
        return False
    from bkrobust.core.oracle import is_valid
    from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag

    optimal_set = optimal_adjustment_set_mpdag(g0, candidate, OUTCOME)
    if optimal_set is None:
        return False
    return is_valid(frozenset(optimal_set), g0, candidate, OUTCOME)


def admissible_queries(layout: UnconfoundedLayout) -> list[str]:
    """Every component vertex (other than ``anchor``) that is a legal query.

    Enumerates the *whole* admissible set before any random choice is made
    (module docstring point 3) -- required for "uniform among admissible",
    which a reject-and-redraw scheme with early stopping would not give.

    Screening uses ``bkrobust.benchmarks.measure.fast_gate`` -- the *same*
    gate ``survival.build_tiered_instance`` itself admits the final instance
    on (not ``synth.runner.gate``, the exhaustive-enumeration one
    ``component_generator.generate_instance`` uses internally, which does not
    scale to being called once per candidate per seed: enumerating *all*
    valid adjustment sets, as that gate does, is a materially different cost
    profile from reading off the optimal one, as ``fast_gate`` does, and
    calling the exhaustive gate a component_size-fold times per seed made the
    pilot draw intractable within seconds). ``fast_gate`` carries the same
    verdict vocabulary per its own docstring, so this is a cost substitution
    on the *screening* step, not a relaxation of the *rule*: identification
    (via ``optimal_adjustment_set_dag``), non-triviality of the empty set,
    and the atomic-perturbation sanity check are all still enforced. The
    tractability/identification check on the truthful-coverage committed
    graph below is the same check ``generate_instance`` makes.
    """
    out: list[str] = []
    for candidate in layout.labels:
        if candidate == layout.anchor:
            continue
        dag = _build_dag_for_candidate(layout, candidate)
        cpdag = dag_to_cpdag(dag)
        ok, _reason = fast_gate(dag, cpdag, candidate, OUTCOME)
        if not ok:
            continue
        if not _truthful_g0_identified_and_valid(dag, cpdag, candidate):
            continue
        out.append(candidate)
    return out


def _measure_separation(cpdag: MPDAG, x: str, z_star: frozenset) -> tuple[Any, str]:
    """BFS separation from ``x`` to the nearest ``z_star`` member inside the
    CPDAG's undirected component containing ``x`` -- mirrors
    ``component_generator.measure_realised``'s separation logic exactly, but
    against the *tiered arm's own* ``z_star``, not a truthful-coverage one.
    """
    component = component_containing(cpdag, x)
    if component is None:
        return None, "treatment_not_in_component"
    members = sorted(v for v in z_star if v in component)
    if not members:
        return None, "no_optimal_member_in_component"
    distances = distances_within_component(cpdag, component, x)
    reachable = [distances[v] for v in members if v in distances]
    if not reachable:
        return None, "no_optimal_member_in_component"
    return int(min(reachable)), "measured"


def build_unconfounded_tiered_instance(
    component_size: int,
    seed: int,
    n_tiers: int,
    *,
    triangle_prob: float = 0.3,
) -> tuple[dict[str, Any] | None, str]:
    """Draw and screen one unconfounded tiered-arm instance.

    Mirrors ``bkrobust.robustness.survival.build_tiered_instance`` downstream
    of graph construction (same ``fast_gate`` call, same
    ``score_committed_graph``, same ``tiered()`` reference knowledge, same
    instance dict shape plus ``separation``/``separation_status``), but the
    graph and the query come from this module's own construction (module
    docstring), not from ``component_generator.generate_instance``.

    Returns:
        ``(instance, "ok")`` on admission, else ``(None, reason)``.
        ``instance`` carries the same keys as ``build_tiered_instance``'s
        output plus ``separation`` (``int`` or ``None``) and
        ``separation_status``.
    """
    started = time.perf_counter()
    rng = np.random.default_rng(seed)
    layout = build_unconfounded_component(component_size, triangle_prob, rng)

    candidates = admissible_queries(layout)
    if not candidates:
        return None, REASON_NO_ADMISSIBLE_QUERY
    x_idx = int(rng.integers(0, len(candidates)))
    x = candidates[x_idx]
    y = OUTCOME

    dag = _build_dag_for_candidate(layout, x)
    cpdag = dag_to_cpdag(dag)

    gate_ok, gate_reason = fast_gate(dag, cpdag, x, y)
    if not gate_ok:
        return None, f"fast_gate:{gate_reason}"

    base_instance_id = (
        f"u2_c{component_size:02d}_seed{seed:08d}_x{x}"
    )
    instance_id = f"{base_instance_id}_tiered_nt{n_tiers}"
    rng0 = np.random.default_rng(sv.derived_seed(instance_id, "tiered_ref"))
    k_ref = tuple(tiered(dag, cpdag, rng0, n_tiers, 0.0))
    g0 = apply_orientations(cpdag, k_ref)

    predictors, reason = sv.score_committed_graph(dag, cpdag, x, y, k_ref, g0)
    if predictors is None:
        return None, reason

    separation, separation_status = _measure_separation(cpdag, x, predictors["z_star"])

    instance = {
        "instance_id": instance_id,
        "base_instance_id": base_instance_id,
        "arm": "tiered_unconfounded",
        "seed": seed,
        "component_size": component_size,
        "separation": separation,
        "separation_status": separation_status,
        "coverage": "",
        "base_wrongness": "",
        "n_tiers": n_tiers,
        "x": x,
        "y": y,
        "anchor": layout.anchor,
        "n_admissible_queries": len(candidates),
        "dag": dag,
        "cpdag": cpdag,
        "g0": g0,
        "k": k_ref,
        **predictors,
        "elapsed_seconds_generation": time.perf_counter() - started,
    }
    return instance, REASON_OK
