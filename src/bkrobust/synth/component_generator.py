r"""A designed generator for CPDAGs with a large, perturbable undirected component.

Why this module exists
----------------------
The Axis A generators in :mod:`bkrobust.synth.generators` sample DAGs and take
whatever chordal components fall out. That is the wrong instrument for the
question this study asks. In the session-4 sweep of 256 instances the *largest*
undirected component at ``n = 24`` was six vertices and 159 of 256 instances had
components of size 2-4; denser Erdos-Renyi DAGs make *more compelled* edges, not
bigger chain components, so raising ``n`` does not help. The quantity that
plausibly governs the breakdown radius is not ``n`` at all: it is the
**separation** ``s``, the graph distance *inside the undirected component* from
the treatment ``X`` to the nearest member of the optimal adjustment set
``O(G0)``. The project's worked example ``Smoke - BMI - Chol - Statin`` has
``s = 3`` and a breakdown radius of 3, the witness being exactly those three
retractions. A three-vertex component cannot exhibit ``s = 3`` at any ``n``. So
the component has to be *built*, not sampled.

The construction
----------------
1. Draw a connected **chordal** graph on ``c`` vertices by repeated
   clique-attachment: vertex ``j`` is attached to a clique among ``0..j-1``.
   Chordal-and-v-structure-free is exactly what survives into a CPDAG as an
   undirected chain component. Attaching to a clique also *cannot* create a
   shortcut between two existing vertices (their common neighbour's endpoints
   are already adjacent), so a spine laid down first keeps its length.
2. The first ``s + 1`` vertices are that spine: ``anchor = p_0, p_1, ..., p_s =
   X``. Edges are oriented **earlier -> later** in construction order, which
   makes every vertex's parent set a clique (no v-structure anywhere inside the
   component) and makes ``anchor -> ... -> X`` a directed path, so ``X`` has a
   genuine ancestor that is also a parent of ``Y``: the back-door path is open
   and the empty set is *not* trivially valid.
3. ``anchor -> Y``, ``X -> Y``, and ``W -> Y`` for a spectator ``W`` adjacent to
   nothing else.

The trap this is designed against
---------------------------------
Session 1's designed family failed twice, in opposite directions, and the two
failures are the two horns of one dilemma: **identification wants compelled
edges near the target, perturbability wants undirected ones.**

* Horn 1 -- leave ``X - Y`` undirected and the effect is not identified at all:
  no adjustment set is valid both in an "X causes Y" world and a "Y causes X"
  world, so there are *zero* valid sets and the key experiment is silently
  untestable.
* Horn 2 -- force ``X -> Y`` with a v-structure whose other parent also sits on
  the confounding structure, and the compulsion propagates: the confounder edges
  freeze too, ``O`` becomes structurally invariant across the whole perturbation
  space, no atomic perturbation changes anything, and 100% of instances are
  gated out.

The escape used here is the placement of ``W``. ``W`` supplies the v-structure
at ``Y`` (``X -> Y <- W`` and ``anchor -> Y <- W`` are both unshielded, so all
three edges into ``Y`` are compelled) while being adjacent to **nothing else**,
so Meek propagation cannot reach the component: R1 needs a directed edge whose
head still carries an undirected edge, and the only head is ``Y``, which carries
none; R2-R4 need directed paths of length two, and ``Y`` is a sink. The
component therefore comes back **fully undirected** at the intended size, and
``O(G0) = pa(Y) \ {X} = {anchor, W}`` sits at distance exactly ``s`` from ``X``.

Realised parameters are measured, never assumed
-----------------------------------------------
Every instance is returned with the realised component size, separation,
``|K_{G0}|``, whether ``X`` really is in a component, whether ``Y`` really is a
descendant of ``X``, and the gate outcome -- all read off the *constructed*
CPDAG via :func:`measure_realised`, not copied from the request. When no member
of ``O(G0)`` lies in the component, ``realised_separation`` is ``None`` and
``realised_separation_status`` carries the reason; a missing separation is never
smuggled in as a number such as ``-1`` alongside real distances.

Determinism: the only randomness is ``numpy.random.default_rng(seed)``, created
once per instance; no global RNG is touched and no hashed container is iterated
in an order that reaches an output.

Written to run on Python 3.9 as well as the repository's 3.11 target: no
``zip(strict=)``, no ``itertools.pairwise``, and PEP 604 unions appear only in
annotations, which ``from __future__ import annotations`` defers.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from bkrobust.core.oracle import is_valid
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG, canon, undirected_components
from bkrobust.demo.meek import apply_orientations
from bkrobust.synth.runner import gate

#: Name of the treatment vertex. Always inside the built component.
TREATMENT: str = "X"

#: Name of the outcome vertex. Always a sink, always a child of ``TREATMENT``.
OUTCOME: str = "Y"

#: Name of the spectator vertex that supplies the v-structure at ``OUTCOME``.
#: Adjacent to ``OUTCOME`` and to nothing else -- see the module docstring.
SPECTATOR: str = "W"

#: Hard cap on the number of undirected edges left in ``G0`` before
#: :func:`bkrobust.demo.evaluate.optimal_adjustment_set_mpdag` is called on it.
#:
#: That function goes through :func:`~bkrobust.demo.meek.enumerate_dag_extensions`,
#: which brute-forces ``2**k`` orientations of the ``k`` undirected edges. At
#: ``coverage = 1.0`` (the default) ``G0`` is fully oriented and ``k = 0``, but
#: a low-coverage analyst leaves most of a size-12 component undirected, and an
#: unbounded ``k`` would hang a sweep on a single draw. Instances above the cap
#: are rejected with ``"g0_too_ambiguous"`` rather than attempted.
MAX_G0_UNDIRECTED_FOR_EXTENSIONS: int = 12

#: Rejection code emitted when ``G0`` exceeds :data:`MAX_G0_UNDIRECTED_FOR_EXTENSIONS`.
REASON_G0_TOO_AMBIGUOUS: str = "g0_too_ambiguous"

#: Rejection code emitted when ``apply_orientations`` FAILs on ``K_{G0}``.
REASON_G0_INCONSISTENT: str = "g0_inconsistent"

#: Rejection code emitted when the optimal set is not identified by ``G0`` alone.
REASON_O_NOT_IDENTIFIED: str = "o_g0_not_identified"

#: Rejection code emitted when ``Z = O(G0)`` is not actually valid at ``G0``.
#: Mirrors the identically named guard in :func:`bkrobust.synth.runner.run_instance`.
REASON_Z_INVALID_AT_G0: str = "z_invalid_at_g0"

#: Rejection code emitted when the realised ``(c, s)`` miss the requested ones.
REASON_PARAMS_NOT_REALISED: str = "params_not_realised"

#: Accepting outcome, matching :func:`bkrobust.synth.runner.gate`'s vocabulary.
REASON_OK: str = "ok"


def achievable_separations(component_size: int) -> tuple[int, ...]:
    """The separations a component of ``component_size`` vertices can exhibit.

    A separation of ``s`` needs a spine of ``s + 1`` distinct vertices inside
    the component, so ``s`` ranges over ``1 .. component_size - 1``. ``s = 0``
    is excluded because the optimal adjustment set never contains the treatment
    itself, so distance zero is not a reachable value.

    Args:
        component_size: Requested number of vertices in the undirected
            component, at least 2.

    Returns:
        The achievable separations, ascending.

    Raises:
        ValueError: If ``component_size`` is below 2.
    """
    if component_size < 2:
        raise ValueError(f"component_size must be >= 2, got {component_size}")
    return tuple(range(1, component_size))


@dataclass(frozen=True)
class ComponentSpec:
    """The *intended* parameters of one instance.

    Attributes:
        component_size: Target number of vertices ``c`` in the undirected
            component of the CPDAG (the treatment is one of them).
        separation: Target distance ``s``, within that component, from the
            treatment to the nearest member of ``O(G0)``.
        coverage: Fraction of the component's undirected edges the analyst
            orients (truthfully) to form ``K_{G0}``. ``1.0`` orients all of
            them, which makes ``G0`` the true DAG and ``O(G0)`` identified by
            construction.
        triangle_prob: Probability that a non-spine vertex is attached to an
            *edge* (making a triangle) rather than to a single vertex. Higher
            values give denser, more clique-y chordal components; ``0.0`` gives
            a tree.
        coverage_order: ``"spine_first"`` fills ``K_{G0}`` starting from the
            component edge incident to the treatment and working outward along
            the spine before covering the rest; ``"random"`` draws a uniformly
            random subset. Irrelevant when ``coverage == 1.0``.
    """

    component_size: int
    separation: int
    coverage: float = 1.0
    triangle_prob: float = 0.3
    coverage_order: str = "spine_first"

    def __post_init__(self) -> None:
        """Validate the requested parameters.

        Raises:
            ValueError: On an out-of-range size, separation, coverage,
                triangle probability, or coverage order.
        """
        if self.component_size < 2:
            raise ValueError(f"component_size must be >= 2, got {self.component_size}")
        if self.separation not in achievable_separations(self.component_size):
            raise ValueError(
                f"separation {self.separation} is not achievable at "
                f"component_size {self.component_size}; achievable: "
                f"{achievable_separations(self.component_size)}"
            )
        if not 0.0 <= self.coverage <= 1.0:
            raise ValueError(f"coverage must be in [0, 1], got {self.coverage}")
        if not 0.0 <= self.triangle_prob <= 1.0:
            raise ValueError(f"triangle_prob must be in [0, 1], got {self.triangle_prob}")
        if self.coverage_order not in ("spine_first", "random"):
            raise ValueError(
                f"coverage_order must be 'spine_first' or 'random', got {self.coverage_order!r}"
            )

    def as_dict(self) -> dict[str, Any]:
        """The spec as a plain JSON-serialisable dict."""
        return {
            "component_size": self.component_size,
            "separation": self.separation,
            "coverage": self.coverage,
            "triangle_prob": self.triangle_prob,
            "coverage_order": self.coverage_order,
        }


@dataclass(frozen=True)
class ComponentLayout:
    """What the builder decided, kept so the measurement step never has to guess.

    Attributes:
        labels: Component vertex labels in construction order. ``labels[0]`` is
            the anchor, ``labels[separation]`` is the treatment.
        anchor: The intended ``O(G0)`` member inside the component.
        spine: The intended shortest path, treatment first, anchor last.
        component_edges: Component edges as ``(tail, head)`` in their *true*
            orientation (earlier -> later in construction order).
    """

    labels: tuple[str, ...]
    anchor: str
    spine: tuple[str, ...]
    component_edges: tuple[tuple[str, str], ...]


@dataclass(frozen=True)
class GeneratedInstance:
    """One generated instance, with intended and realised parameters side by side.

    Attributes:
        instance_id: Stable identifier.
        seed: The instance's root seed; the sole source of randomness.
        intended: The requested parameters (see :meth:`ComponentSpec.as_dict`).
        realised: Parameters *measured on the constructed CPDAG*; see
            :func:`measure_realised` for the exact keys.
        accepted: Whether the instance passed every gate and realised its
            intended ``(c, s)``.
        reject_reason: :data:`REASON_OK` when accepted, else the first failing
            check's code.
        treatment: The treatment vertex.
        outcome: The outcome vertex.
        anchor: The intended in-component member of ``O(G0)``.
        true_dag: ``MPDAG.edge_string()`` of the ground-truth DAG.
        cpdag: ``MPDAG.edge_string()`` of its CPDAG.
        g0: ``MPDAG.edge_string()`` of the analyst's graph, or ``""`` if
            ``K_{G0}`` was inconsistent.
        k_g0: The analyst's asserted orientations, sorted.
        optimal_set: ``O(G0)`` sorted, or ``None`` if not identified.
        elapsed_seconds: Wall-clock cost of this draw.
        dag_obj: The ground-truth DAG object (not serialised).
        cpdag_obj: The CPDAG object (not serialised).
        g0_obj: The analyst's graph object, or ``None`` (not serialised).
    """

    instance_id: str
    seed: int
    intended: dict[str, Any]
    realised: dict[str, Any]
    accepted: bool
    reject_reason: str
    treatment: str
    outcome: str
    anchor: str
    true_dag: str
    cpdag: str
    g0: str
    k_g0: tuple[tuple[str, str], ...]
    optimal_set: tuple[str, ...] | None
    elapsed_seconds: float = 0.0
    dag_obj: MPDAG | None = field(default=None, repr=False, compare=False)
    cpdag_obj: MPDAG | None = field(default=None, repr=False, compare=False)
    g0_obj: MPDAG | None = field(default=None, repr=False, compare=False)

    def fingerprint(self) -> str:
        """A stable digest of everything except wall-clock timing.

        Two instances drawn from the same ``(spec, seed)`` must have equal
        fingerprints, under any ``PYTHONHASHSEED``. Timing is excluded because
        it is the one field that legitimately varies between runs.

        Returns:
            A hex SHA-256 digest.
        """
        payload = json.dumps(self.to_row(), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    def to_row(self) -> dict[str, Any]:
        """A flat, JSON-serialisable record of the instance, timing excluded."""
        return {
            "instance_id": self.instance_id,
            "seed": self.seed,
            "intended": dict(self.intended),
            "realised": dict(self.realised),
            "accepted": self.accepted,
            "reject_reason": self.reject_reason,
            "treatment": self.treatment,
            "outcome": self.outcome,
            "anchor": self.anchor,
            "true_dag": self.true_dag,
            "cpdag": self.cpdag,
            "g0": self.g0,
            "k_g0": [list(e) for e in self.k_g0],
            "optimal_set": None if self.optimal_set is None else list(self.optimal_set),
        }


# --- construction ---------------------------------------------------------


def _component_labels(component_size: int, separation: int) -> tuple[str, ...]:
    """Labels for the component's vertices, in construction order.

    Position ``separation`` is the treatment; everything else gets a
    zero-padded ``C`` label whose string sort matches its construction index,
    since :class:`~bkrobust.demo.graph.MPDAG` sorts node names as strings.

    Args:
        component_size: Number of component vertices.
        separation: Construction index of the treatment.

    Returns:
        The labels, in construction order.
    """
    width = max(2, len(str(component_size - 1)))
    return tuple(TREATMENT if j == separation else f"C{j:0{width}d}" for j in range(component_size))


def build_component_dag(
    spec: ComponentSpec, rng: np.random.Generator
) -> tuple[MPDAG, ComponentLayout]:
    """Build the ground-truth DAG for one instance.

    The component is grown by clique-attachment (see the module docstring): the
    spine ``anchor = p_0 - ... - p_s = X`` is laid down first, then the
    remaining ``c - s - 1`` vertices each attach to a clique of what already
    exists -- a single vertex, or an existing edge with probability
    ``spec.triangle_prob``. Attaching to a clique preserves chordality and
    cannot shorten any existing distance, so ``d(X, anchor)`` stays exactly
    ``s``. Every edge is oriented earlier -> later in construction order, which
    makes each vertex's parents a clique (no v-structure inside the component)
    and makes ``anchor`` an ancestor of ``X``.

    Args:
        spec: The intended parameters.
        rng: Sole source of randomness.

    Returns:
        The ground-truth DAG and the layout the builder committed to.
    """
    c = spec.component_size
    s = spec.separation
    labels = _component_labels(c, s)

    adjacency: dict[int, set[int]] = {j: set() for j in range(c)}
    index_edges: list[tuple[int, int]] = []

    def connect(earlier: int, later: int) -> None:
        index_edges.append((earlier, later))
        adjacency[earlier].add(later)
        adjacency[later].add(earlier)

    for j in range(1, s + 1):
        connect(j - 1, j)

    for j in range(s + 1, c):
        anchor_idx = int(rng.integers(0, j))
        attach_to = [anchor_idx]
        neighbours = sorted(adjacency[anchor_idx])
        if neighbours and rng.uniform() < spec.triangle_prob:
            attach_to.append(neighbours[int(rng.integers(0, len(neighbours)))])
        for base in sorted(attach_to):
            connect(base, j)

    component_edges = tuple((labels[a], labels[b]) for a, b in index_edges)
    anchor = labels[0]
    spine = tuple(labels[j] for j in range(s, -1, -1))

    nodes = [*labels, OUTCOME, SPECTATOR]
    directed = [
        *component_edges,
        (TREATMENT, OUTCOME),
        (anchor, OUTCOME),
        (SPECTATOR, OUTCOME),
    ]
    dag = MPDAG(nodes, directed=directed)
    layout = ComponentLayout(
        labels=labels,
        anchor=anchor,
        spine=spine,
        component_edges=component_edges,
    )
    return dag, layout


def select_knowledge(
    dag: MPDAG,
    cpdag: MPDAG,
    layout: ComponentLayout,
    spec: ComponentSpec,
    rng: np.random.Generator,
) -> tuple[tuple[str, str], ...]:
    """Choose the analyst's background knowledge ``K_{G0}``.

    ``K_{G0}`` is always *truthful*: each claim is the orientation the edge
    actually has in ``dag``. Only how *many* claims, and which, is under the
    caller's control -- corrupting them is the job of
    :mod:`bkrobust.synth.knowledge`, downstream of this generator.

    Args:
        dag: The ground-truth DAG, read for the true orientations.
        cpdag: Its CPDAG, read for which edges are still undirected.
        layout: The builder's layout, for the spine order.
        spec: The intended parameters, for ``coverage`` and ``coverage_order``.
        rng: Sole source of randomness (used only for ``"random"`` ordering).

    Returns:
        The asserted orientations as sorted ``(tail, head)`` pairs.
    """
    undirected = {canon(a, b) for a, b in cpdag.undirected_edges}
    candidates = [e for e in layout.component_edges if canon(e[0], e[1]) in undirected]
    if not candidates:
        return ()

    if spec.coverage_order == "spine_first":
        spine_keys = [
            canon(layout.spine[i], layout.spine[i + 1]) for i in range(len(layout.spine) - 1)
        ]
        rank = {key: pos for pos, key in enumerate(spine_keys)}
        ordered = sorted(
            candidates,
            key=lambda e: (rank.get(canon(e[0], e[1]), len(rank)), e),
        )
    else:
        permutation = rng.permutation(len(candidates))
        ordered = [candidates[int(i)] for i in permutation]

    n_keep = round(spec.coverage * len(candidates))
    n_keep = max(0, min(len(candidates), n_keep))
    chosen = ordered[:n_keep]
    return tuple(sorted((tail, head) for tail, head in chosen if dag.is_directed_edge(tail, head)))


# --- measurement ----------------------------------------------------------


def component_containing(graph: MPDAG, node: str) -> frozenset[str] | None:
    """The undirected component of ``graph`` that contains ``node``, if any.

    Args:
        graph: Any MPDAG.
        node: The vertex to locate.

    Returns:
        The component as a frozen set, or ``None`` if ``node`` has no
        undirected edge at all.
    """
    for component in undirected_components(graph):
        if node in component:
            return component
    return None


def distances_within_component(
    graph: MPDAG, component: frozenset[str], source: str
) -> dict[str, int]:
    """BFS distances from ``source`` over the undirected edges inside ``component``.

    Args:
        graph: The graph whose undirected edges define the adjacency.
        component: The vertex set to stay inside.
        source: The vertex to measure from.

    Returns:
        A mapping from reachable vertex to hop count; ``source`` maps to 0.
    """
    adjacency: dict[str, list[str]] = {v: [] for v in sorted(component)}
    for a, b in sorted(graph.undirected_edges):
        if a in component and b in component:
            adjacency[a].append(b)
            adjacency[b].append(a)
    for key in adjacency:
        adjacency[key].sort()

    distances = {source: 0}
    frontier = [source]
    while frontier:
        nxt: list[str] = []
        for current in frontier:
            for neighbour in adjacency[current]:
                if neighbour not in distances:
                    distances[neighbour] = distances[current] + 1
                    nxt.append(neighbour)
        frontier = sorted(nxt)
    return distances


def measure_realised(
    dag: MPDAG,
    cpdag: MPDAG,
    g0: MPDAG | None,
    optimal_set: set[str] | None,
    k_g0: Sequence[tuple[str, str]],
    treatment: str,
    outcome: str,
) -> dict[str, Any]:
    """Measure the realised parameters on the *constructed* graphs.

    Nothing here reads the request: every value is read off ``cpdag``, ``dag``
    and ``g0``. The separation is measured inside the CPDAG's undirected
    component containing the treatment -- that component, not ``G0``'s residual
    one, is the arena the perturbation space ranges over, so it is the arena the
    hypothesis "radius tracks separation" is about. ``G0``'s own component is
    reported separately as ``realised_g0_component_size``; at ``coverage = 1.0``
    it is empty by construction, because the analyst has oriented everything.

    Args:
        dag: The ground-truth DAG.
        cpdag: Its CPDAG.
        g0: The analyst's graph, or ``None`` if ``K_{G0}`` was inconsistent.
        optimal_set: ``O(G0)``, or ``None`` if not identified.
        k_g0: The analyst's asserted orientations.
        treatment: The treatment vertex.
        outcome: The outcome vertex.

    Returns:
        A dict with keys ``realised_component_size``,
        ``realised_component_edges``, ``realised_separation`` (an ``int`` or
        ``None`` -- never a numeric sentinel),
        ``realised_separation_status``, ``realised_k_g0``,
        ``realised_g0_undirected_edges``, ``realised_g0_component_size``,
        ``realised_n_nodes``, ``realised_cpdag_undirected_edges``,
        ``x_in_component``, ``y_descendant_of_x``,
        ``o_members_in_component`` and ``optimal_set_identified``.
    """
    component = component_containing(cpdag, treatment)
    x_in_component = component is not None
    component_size = 0 if component is None else len(component)
    component_edges = (
        0
        if component is None
        else sum(1 for a, b in cpdag.undirected_edges if a in component and b in component)
    )

    separation: int | None = None
    status: str
    members_in_component: list[str] = []
    if optimal_set is None:
        status = "optimal_set_not_identified"
    elif component is None:
        status = "treatment_not_in_component"
    else:
        members_in_component = sorted(v for v in optimal_set if v in component)
        if not members_in_component:
            status = "no_optimal_member_in_component"
        else:
            distances = distances_within_component(cpdag, component, treatment)
            reachable = [distances[v] for v in members_in_component if v in distances]
            if not reachable:
                status = "no_optimal_member_in_component"
            else:
                separation = int(min(reachable))
                status = "measured"

    g0_component = None if g0 is None else component_containing(g0, treatment)

    return {
        "realised_component_size": component_size,
        "realised_component_edges": component_edges,
        "realised_separation": separation,
        "realised_separation_status": status,
        "realised_k_g0": len(k_g0),
        "realised_g0_undirected_edges": 0 if g0 is None else len(g0.undirected_edges),
        "realised_g0_component_size": 0 if g0_component is None else len(g0_component),
        "realised_n_nodes": len(cpdag.nodes),
        "realised_cpdag_undirected_edges": len(cpdag.undirected_edges),
        "x_in_component": x_in_component,
        "y_descendant_of_x": outcome in dag.descendants(treatment),
        "o_members_in_component": members_in_component,
        "optimal_set_identified": optimal_set is not None,
    }


# --- one instance ---------------------------------------------------------


def generate_instance(
    spec: ComponentSpec, seed: int, instance_id: str | None = None
) -> GeneratedInstance:
    """Draw one instance and screen it, returning it either way.

    Unlike :func:`bkrobust.synth.runner.run_instance`, this never raises on a
    rejection: a rejected instance is still returned, with its realised
    parameters measured, so a pilot can report *why* a cell fails rather than
    just that it did.

    Screening order, each short-circuiting the rest:

    1. :func:`bkrobust.synth.runner.gate` on ``(dag, cpdag, X, Y)`` -- the five
       structural degeneracy checks, reused verbatim rather than reimplemented.
    2. ``K_{G0}`` is consistent (``apply_orientations`` did not FAIL).
    3. ``G0`` is tractable for extension enumeration
       (:data:`MAX_G0_UNDIRECTED_FOR_EXTENSIONS`).
    4. ``O(G0)`` is identified, and ``Z = O(G0)`` is genuinely valid at ``G0``
       (mirroring ``run_instance``'s unconditional ``z_invalid_at_g0`` guard).
    5. The realised ``(c, s)`` equal the intended ones.

    Args:
        spec: The intended parameters.
        seed: Root seed; the sole source of randomness.
        instance_id: Identifier to stamp on the result. Defaults to a
            deterministic function of the spec and seed.

    Returns:
        The instance, accepted or not.
    """
    started = time.perf_counter()
    rng = np.random.default_rng(seed)

    dag, layout = build_component_dag(spec, rng)
    cpdag = dag_to_cpdag(dag)
    k_g0 = select_knowledge(dag, cpdag, layout, spec, rng)

    g0 = apply_orientations(cpdag, k_g0)
    optimal_set: set[str] | None = None
    reason = REASON_OK

    accepted, gate_reason = gate(dag, cpdag, TREATMENT, OUTCOME)
    if not accepted:
        reason = gate_reason

    if reason == REASON_OK and g0 is None:
        reason = REASON_G0_INCONSISTENT
    if reason == REASON_OK:
        assert g0 is not None  # narrowed by the check above
        if len(g0.undirected_edges) > MAX_G0_UNDIRECTED_FOR_EXTENSIONS:
            reason = REASON_G0_TOO_AMBIGUOUS
        else:
            optimal_set = optimal_adjustment_set_mpdag(g0, TREATMENT, OUTCOME)
            if optimal_set is None:
                reason = REASON_O_NOT_IDENTIFIED
            elif not is_valid(frozenset(optimal_set), g0, TREATMENT, OUTCOME):
                reason = REASON_Z_INVALID_AT_G0

    realised = measure_realised(dag, cpdag, g0, optimal_set, k_g0, TREATMENT, OUTCOME)

    if reason == REASON_OK and (
        realised["realised_component_size"] != spec.component_size
        or realised["realised_separation"] != spec.separation
    ):
        reason = REASON_PARAMS_NOT_REALISED

    if instance_id is None:
        instance_id = (
            f"c{spec.component_size:02d}_s{spec.separation:02d}"
            f"_cov{round(spec.coverage * 100):03d}_seed{seed:08d}"
        )

    return GeneratedInstance(
        instance_id=instance_id,
        seed=seed,
        intended=spec.as_dict(),
        realised=realised,
        accepted=reason == REASON_OK,
        reject_reason=reason,
        treatment=TREATMENT,
        outcome=OUTCOME,
        anchor=layout.anchor,
        true_dag=dag.edge_string(),
        cpdag=cpdag.edge_string(),
        g0="" if g0 is None else g0.edge_string(),
        k_g0=k_g0,
        optimal_set=None if optimal_set is None else tuple(sorted(optimal_set)),
        elapsed_seconds=time.perf_counter() - started,
        dag_obj=dag,
        cpdag_obj=cpdag,
        g0_obj=g0,
    )


# --- the pilot ------------------------------------------------------------


def _cell_report(spec: ComponentSpec, instances: Sequence[GeneratedInstance]) -> dict[str, Any]:
    """Summarise one ``(c, s)`` cell's draws.

    Args:
        spec: The cell's intended parameters.
        instances: Every draw attempted in the cell.

    Returns:
        A JSON-serialisable per-cell record.
    """
    attempts = len(instances)
    accepted = [i for i in instances if i.accepted]
    reasons: dict[str, int] = {}
    for inst in instances:
        if inst.accepted:
            continue
        reasons[inst.reject_reason] = reasons.get(inst.reject_reason, 0) + 1

    c_match = sum(
        1 for i in instances if i.realised["realised_component_size"] == spec.component_size
    )
    s_match = sum(1 for i in instances if i.realised["realised_separation"] == spec.separation)
    separations = [
        i.realised["realised_separation"]
        for i in instances
        if i.realised["realised_separation"] is not None
    ]
    total_seconds = sum(i.elapsed_seconds for i in instances)
    undirected = [i.realised["realised_cpdag_undirected_edges"] for i in instances]

    return {
        "component_size": spec.component_size,
        "separation": spec.separation,
        "attempts": attempts,
        "accepted": len(accepted),
        "acceptance_rate": (len(accepted) / attempts) if attempts else 0.0,
        "rejection_reasons": dict(sorted(reasons.items())),
        "intended_c_equals_realised_c": c_match,
        "intended_s_equals_realised_s": s_match,
        "max_realised_separation": max(separations) if separations else None,
        "n_separation_unmeasurable": attempts - len(separations),
        "mean_cpdag_undirected_edges": (sum(undirected) / attempts) if attempts else 0.0,
        "max_cpdag_undirected_edges": max(undirected) if undirected else 0,
        "seconds_total": total_seconds,
        "seconds_per_attempt": (total_seconds / attempts) if attempts else 0.0,
        "seconds_per_accepted": (total_seconds / len(accepted)) if accepted else None,
        "example_fingerprint": accepted[0].fingerprint() if accepted else None,
    }


def run_pilot(
    component_sizes: Sequence[int] = tuple(range(2, 13)),
    budget: int = 24,
    root_seed: int = 20260906,
    triangle_prob: float = 0.3,
    coverage: float = 1.0,
    coverage_sweep: Sequence[float] = (0.0, 0.25, 0.5, 0.75, 1.0),
    coverage_sweep_size: int = 8,
    coverage_sweep_separation: int = 4,
    coverage_sweep_budget: int = 8,
    radius_spotcheck_max_size: int = 7,
) -> dict[str, Any]:
    """Sweep every ``(c, s)`` cell with a fixed budget and report what happened.

    Every cell of the product ``c in component_sizes`` by
    ``s in achievable_separations(c)`` is attempted ``budget`` times. Nothing is
    resampled around a failing cell: a cell that cannot be filled shows up as an
    acceptance rate of 0 with its rejection reasons intact.

    Args:
        component_sizes: The component sizes to sweep.
        budget: Draws attempted per cell.
        root_seed: Root seed; every draw's seed is a fixed function of it and
            the draw's position in the sweep, never of how much has already run.
        triangle_prob: Passed to every :class:`ComponentSpec`.
        coverage: Passed to every ``(c, s)`` cell's :class:`ComponentSpec`.
        coverage_sweep: Coverages for the secondary knowledge-coverage sweep.
        coverage_sweep_size: Component size for that secondary sweep.
        coverage_sweep_separation: Separation for that secondary sweep.
        coverage_sweep_budget: Draws per (coverage, order) point.
        radius_spotcheck_max_size: Largest component size handed to
            :func:`_run_radius_spotcheck`.

    Returns:
        A JSON-serialisable pilot report; see :func:`main` for where it lands.
    """
    started = time.perf_counter()
    cells: list[dict[str, Any]] = []
    all_instances: list[GeneratedInstance] = []
    cell_index = 0

    for size in component_sizes:
        for sep in achievable_separations(size):
            spec = ComponentSpec(
                component_size=size,
                separation=sep,
                coverage=coverage,
                triangle_prob=triangle_prob,
            )
            draws: list[GeneratedInstance] = []
            for draw in range(budget):
                seed = root_seed + 1000 * cell_index + draw
                draws.append(generate_instance(spec, seed))
            all_instances.extend(draws)
            cells.append(_cell_report(spec, draws))
            cell_index += 1

    overall_reasons: dict[str, int] = {}
    for inst in all_instances:
        if inst.accepted:
            continue
        overall_reasons[inst.reject_reason] = overall_reasons.get(inst.reject_reason, 0) + 1

    accepted = [i for i in all_instances if i.accepted]
    achieved: dict[int, int] = {}
    for inst in accepted:
        sep = inst.realised["realised_separation"]
        size = inst.realised["realised_component_size"]
        if sep is None:
            continue
        if sep > achieved.get(size, -1):
            achieved[size] = sep

    timings = {}
    for size in (6, 9, 12):
        subset = [i for i in accepted if i.intended["component_size"] == size]
        pool = [i for i in all_instances if i.intended["component_size"] == size]
        timings[str(size)] = {
            "n_accepted": len(subset),
            "n_attempts": len(pool),
            "seconds_per_accepted_instance": (
                sum(i.elapsed_seconds for i in pool) / len(subset) if subset else None
            ),
        }

    coverage_report = _run_coverage_sweep(
        size=coverage_sweep_size,
        separation=coverage_sweep_separation,
        coverages=coverage_sweep,
        budget=coverage_sweep_budget,
        triangle_prob=triangle_prob,
        root_seed=root_seed + 900_000,
    )

    spotcheck_report = _run_radius_spotcheck(
        max_size=radius_spotcheck_max_size, root_seed=root_seed + 800_000
    )

    digest_payload = json.dumps(
        [i.fingerprint() for i in all_instances], sort_keys=True, separators=(",", ":")
    )

    return {
        "config": {
            "component_sizes": list(component_sizes),
            "budget_per_cell": budget,
            "root_seed": root_seed,
            "triangle_prob": triangle_prob,
            "coverage": coverage,
            "max_g0_undirected_for_extensions": MAX_G0_UNDIRECTED_FOR_EXTENSIONS,
        },
        "cells": cells,
        "totals": {
            "attempts": len(all_instances),
            "accepted": len(accepted),
            "acceptance_rate": len(accepted) / len(all_instances) if all_instances else 0.0,
            "rejection_reasons": dict(sorted(overall_reasons.items())),
        },
        "realisation_fidelity": {
            "n_attempts": len(all_instances),
            "intended_c_equals_realised_c": sum(
                1
                for i in all_instances
                if i.realised["realised_component_size"] == i.intended["component_size"]
            ),
            "intended_s_equals_realised_s": sum(
                1
                for i in all_instances
                if i.realised["realised_separation"] == i.intended["separation"]
            ),
            "n_separation_unmeasurable": sum(
                1 for i in all_instances if i.realised["realised_separation"] is None
            ),
            "separation_status_counts": _count(
                i.realised["realised_separation_status"] for i in all_instances
            ),
        },
        "separation_range": {
            "max_realised_separation_overall": max(
                (
                    i.realised["realised_separation"]
                    for i in accepted
                    if i.realised["realised_separation"] is not None
                ),
                default=None,
            ),
            "max_realised_separation_by_component_size": {
                str(k): v for k, v in sorted(achieved.items())
            },
            "unfilled_cells": [
                {
                    "component_size": cell["component_size"],
                    "separation": cell["separation"],
                    "rejection_reasons": cell["rejection_reasons"],
                }
                for cell in cells
                if cell["accepted"] == 0
            ],
        },
        "timings": {
            "wall_seconds_total": time.perf_counter() - started,
            "by_component_size": timings,
        },
        "coverage_sweep": coverage_report,
        "radius_spotcheck": spotcheck_report,
        "determinism": {
            "digest": hashlib.sha256(digest_payload.encode("utf-8")).hexdigest(),
            "note": (
                "SHA-256 over every attempted instance's fingerprint, in sweep order. "
                "Timing is excluded, so this digest must be identical under any "
                "PYTHONHASHSEED."
            ),
        },
    }


def _count(values: Any) -> dict[str, int]:
    """Tally an iterable of hashable labels into a sorted dict."""
    out: dict[str, int] = {}
    for value in values:
        key = str(value)
        out[key] = out.get(key, 0) + 1
    return dict(sorted(out.items()))


def _run_coverage_sweep(
    size: int,
    separation: int,
    coverages: Sequence[float],
    budget: int,
    triangle_prob: float,
    root_seed: int,
) -> dict[str, Any]:
    """Sweep the analyst's knowledge coverage at one fixed ``(c, s)``.

    Coverage is the third control knob (alongside ``c`` and ``s``) and drives
    ``|K_{G0}|``. It is swept separately from the ``(c, s)`` grid because a
    low-coverage ``G0`` leaves the component largely undirected, where
    ``O(G0)`` is typically *not* identified -- an honest finding about the knob,
    not a defect to sample around.

    Args:
        size: Component size to hold fixed.
        separation: Separation to hold fixed.
        coverages: The coverage values to sweep.
        budget: Draws per (coverage, order) point.
        triangle_prob: Passed through to the spec.
        root_seed: Root seed for the sweep.

    Returns:
        A JSON-serialisable record, one row per (coverage, coverage_order).
    """
    rows: list[dict[str, Any]] = []
    point = 0
    for order in ("spine_first", "random"):
        for cov in coverages:
            spec = ComponentSpec(
                component_size=size,
                separation=separation,
                coverage=cov,
                triangle_prob=triangle_prob,
                coverage_order=order,
            )
            draws = [generate_instance(spec, root_seed + 100 * point + d) for d in range(budget)]
            point += 1
            accepted = [i for i in draws if i.accepted]
            reasons: dict[str, int] = {}
            for inst in draws:
                if not inst.accepted:
                    reasons[inst.reject_reason] = reasons.get(inst.reject_reason, 0) + 1
            rows.append(
                {
                    "coverage": cov,
                    "coverage_order": order,
                    "attempts": len(draws),
                    "accepted": len(accepted),
                    "acceptance_rate": len(accepted) / len(draws) if draws else 0.0,
                    "mean_realised_k_g0": (
                        sum(i.realised["realised_k_g0"] for i in draws) / len(draws)
                        if draws
                        else 0.0
                    ),
                    "mean_g0_undirected_edges": (
                        sum(i.realised["realised_g0_undirected_edges"] for i in draws) / len(draws)
                        if draws
                        else 0.0
                    ),
                    "rejection_reasons": dict(sorted(reasons.items())),
                }
            )
    return {
        "component_size": size,
        "separation": separation,
        "budget_per_point": budget,
        "rows": rows,
    }


def _run_radius_spotcheck(
    max_size: int,
    root_seed: int,
    triangle_prob: float = 0.0,
    max_undirected: int = 7,
) -> dict[str, Any]:
    """Check, on small cells only, that the generated instances are actually usable.

    The point of the generator is not merely to realise ``(c, s)`` on paper but
    to hand the census instances whose breakdown radius is *measurable and
    non-trivial*. This computes the exact BFS radius ``r_val`` on the smallest
    cells, where :func:`~bkrobust.core.spacelib.build_space`'s ``3**k`` brute
    force is affordable, and records it next to the realised separation. It is
    strictly bounded: tree components only (``triangle_prob = 0``), one draw per
    cell, and any CPDAG with more than ``max_undirected`` undirected edges is
    skipped rather than attempted.

    Args:
        max_size: Largest component size to spot-check.
        root_seed: Root seed for the spot-check draws.
        triangle_prob: Passed to the spec; ``0`` keeps ``k`` at ``c - 1``.
        max_undirected: Skip any CPDAG above this many undirected edges.

    Returns:
        A JSON-serialisable record, one row per spot-checked cell.
    """
    from bkrobust.core.spacelib import build_space, distances_from, radius

    rows: list[dict[str, Any]] = []
    point = 0
    for size in range(3, max_size + 1):
        for sep in achievable_separations(size):
            spec = ComponentSpec(component_size=size, separation=sep, triangle_prob=triangle_prob)
            inst = generate_instance(spec, root_seed + point)
            point += 1
            row: dict[str, Any] = {
                "component_size": size,
                "separation": sep,
                "accepted": inst.accepted,
                "realised_separation": inst.realised["realised_separation"],
                "cpdag_undirected_edges": inst.realised["realised_cpdag_undirected_edges"],
                "r_val": None,
                "space_size": None,
                "skipped": None,
            }
            if not inst.accepted:
                row["skipped"] = f"rejected:{inst.reject_reason}"
            elif inst.realised["realised_cpdag_undirected_edges"] > max_undirected:
                row["skipped"] = "too_many_undirected_edges"
            else:
                assert inst.cpdag_obj is not None and inst.g0_obj is not None
                space = build_space(inst.cpdag_obj)
                if inst.g0_obj not in space.neighbours:
                    row["skipped"] = "g0_not_in_space"
                else:
                    distances = distances_from(space, inst.g0_obj)
                    z = frozenset(inst.optimal_set or ())
                    r_val, _ = radius(
                        space,
                        distances,
                        lambda g, z=z: not is_valid(z, g, TREATMENT, OUTCOME),
                    )
                    row["r_val"] = r_val
                    row["space_size"] = len(space)
            rows.append(row)

    checked = [r for r in rows if r["r_val"] is not None]
    return {
        "note": (
            "Exact BFS breakdown radius on the small cells, to confirm the "
            "generated instances are non-vacuous downstream, not merely "
            "correctly parameterised."
        ),
        "max_component_size": max_size,
        "n_rows": len(rows),
        "n_checked": len(checked),
        "n_r_val_equals_separation": sum(
            1 for r in checked if r["r_val"] == r["realised_separation"]
        ),
        "rows": rows,
    }


def main(argv: Sequence[str] | None = None) -> int:
    """Run the pilot and write it to ``results/axisa2/generator_pilot.json``.

    Args:
        argv: Command-line arguments; ``None`` means ``sys.argv[1:]``.

    Returns:
        Process exit code (0).
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--out",
        default="results/axisa2/generator_pilot.json",
        help="Where to write the pilot report.",
    )
    parser.add_argument("--budget", type=int, default=24, help="Draws attempted per (c, s) cell.")
    parser.add_argument("--root-seed", type=int, default=20260906, help="Root seed.")
    parser.add_argument("--max-size", type=int, default=12, help="Largest component size swept.")
    parser.add_argument(
        "--digest-only",
        action="store_true",
        help="Print only the determinism digest; do not write the report.",
    )
    args = parser.parse_args(argv)

    report = run_pilot(
        component_sizes=tuple(range(2, args.max_size + 1)),
        budget=args.budget,
        root_seed=args.root_seed,
    )
    if args.digest_only:
        print(report["determinism"]["digest"])
        return 0

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"wrote {out_path}")
    print(
        "accepted {accepted}/{attempts} ({rate:.1%})".format(
            accepted=report["totals"]["accepted"],
            attempts=report["totals"]["attempts"],
            rate=report["totals"]["acceptance_rate"],
        )
    )
    print("digest", report["determinism"]["digest"])
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
