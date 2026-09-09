r"""Phase 2, continued: P6' (achievable efficiency over a bounded menu) and
P7' (aggressive background knowledge on real benchmark networks).

See ``results/axis_robustness/PREREGISTRATION.md`` **Appendix A** for why this
module exists. Round 1 (:mod:`bkrobust.robustness.pareto`) measured the
efficiency of the *optimal* adjustment set of every truthful proposal and
found it constant within every stratum (0/24 strata varied): the optimal set
is an invariant of the MPDAG whenever it is identified at all, so no
efficiency trade-off could exist to detect. P6' replaces "the efficiency of
the optimal set" with "the efficiency of the best set the analyst can *justify
from G0*", measured over a small, fixed, auditable menu of candidate sets
built once per base CPDAG from the true DAG -- **not** the exhaustive
enumerator over every subset of V that stays off every decision path in this
project. P7' repeats the aggressive/conservative contrast on real benchmark
CPDAGs, because the synthetic generator used in round 1 structurally cannot
produce an undirected edge lying on a treatment-to-outcome path.

This module deliberately imports round 1's machinery
(:mod:`bkrobust.robustness.pareto`) rather than re-deriving it: base-CPDAG
selection, the truthful-proposal enumerator, and SEM handling are read-only
reused from there. Nothing in :mod:`bkrobust.robustness.pareto` is edited by
this module.

Live-tree only: every import below resolves to a *live* module
(``bkrobust.demo.*``, ``bkrobust.synth.*`` [only ``knowledge``, transitively
via ``pareto``], ``bkrobust.hybrid``, ``bkrobust.benchmarks.*``,
``bkrobust.gac``, ``bkrobust.core.*``). Nothing here imports the dead
top-level scaffold (``bkrobust.graphs``, ``bkrobust.knowledge``,
``bkrobust.metrics``, ``bkrobust.theory``, ``bkrobust.estimation``,
``bkrobust.data``, ``bkrobust.cfm``, ``bkrobust.representations``,
``bkrobust.utils``).

Gate discipline: the only gate ever called directly by this module (or by
:mod:`bkrobust.robustness.pareto`, which it imports) is
:func:`bkrobust.benchmarks.measure.fast_gate`. The exponential all-DAG-subset
adjustment-set enumerator in ``bkrobust.demo.evaluate`` and the five-check
gate living in the synth-data generator's own ``runner`` submodule are never
invoked from here. Validity of a candidate menu member against a partially
oriented graph is decided by the polynomial graphical criterion
(``bkrobust.gac.is_gac_valid_mpdag``), which does not enumerate DAG
extensions, so the menu-validity filter below carries none of the scaling
risk that motivated dropping ``optimal_adjustment_set_mpdag`` as the primary
endpoint.

Randomness: every draw is seeded from an explicit ``numpy.random.Generator``,
itself derived deterministically (via SHA-256, insensitive to the
interpreter's string-hash seed) from small integer/string keys, reusing
:func:`bkrobust.robustness.pareto.derive_seed`. Neither of the two global
random-number entry points in the standard library or NumPy is touched
anywhere in this module, and the standard library's random-number module is
never imported.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from bkrobust.benchmarks.describe import parse_file, to_mpdag
from bkrobust.benchmarks.measure import fast_gate
from bkrobust.demo.evaluate import LinearSEM, bias, optimal_adjustment_set_dag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.gac import is_gac_valid_mpdag
from bkrobust.hybrid import breakdown_radius
from bkrobust.robustness import pareto as P

Node = str
Edge = tuple[str, str]

GATE_NAME = "fast_gate"

#: Menu cap, per the design doc (Appendix A.3): "a fixed candidate menu of at
#: most 12 sets constructed once per base CPDAG".
MENU_CAP = 12
#: Up to this many O_true +/- v perturbation sets each, per the design doc.
MENU_PERTURB_CAP = 4

NETWORKS_DIR = Path("results/axisa3/networks/example_models")
#: Excluded from the corpus per the task brief: an ADMG (bidirected edges),
#: not a DAG.
EXCLUDED_NETWORK_FILES = {"M-bias.txt"}
#: Ordered-pair census budget per network (Task B census; structural only,
#: no gate call).
CENSUS_PAIR_CAP = 200
#: Real-network follow-up: how many of the smallest on-path-bearing networks
#: to gate and measure, and how many (x, y) pairs per chosen network.
FOLLOWUP_MAX_NETWORKS = 6
FOLLOWUP_MIN_NETWORKS = 3
FOLLOWUP_PAIRS_PER_NETWORK = 3


# ---------------------------------------------------------------------------
# Task A -- the bounded candidate menu
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class MenuEntry:
    source: str
    members: frozenset[str]


def build_menu(dag: MPDAG, x: Node, y: Node) -> list[MenuEntry]:
    """The fixed, auditable candidate menu for ``(x, y)``, built once from the
    TRUE dag (never from a proposal's ``G0``), per Appendix A.3:

    1. ``O_true`` -- the true DAG's optimal adjustment set.
    2. The empty set.
    3. ``pa_true(X)``.
    4. ``an_true(Y) \\ ({X} u de_true(X))``.
    5. Up to 4 sets ``O_true u {v}``, ``v`` ranging over the non-descendants
       of ``X`` in the true DAG, sorted, restricted to ``v`` not already in
       ``O_true u {x, y}`` (so up to 4 *new* sets are actually offered, not
       up to 4 candidates before filtering -- documented design choice, see
       ``PARETO2_RUN_NOTES.md``).
    6. Up to 4 sets ``O_true \\ {v}``, ``v`` ranging over the first 4 members
       of ``sorted(O_true)``.

    Deduplicated by member set (a set produced by two constructions keeps
    both source labels, semicolon-joined), sorted canonically by
    ``tuple(sorted(members))``, and capped at :data:`MENU_CAP`.
    """
    o_true = frozenset(optimal_adjustment_set_dag(dag, x, y))
    pa_true = frozenset(dag.parents(x))
    an_y = dag.ancestors(y) | {y}
    de_x = dag.descendants(x) | {x}
    an_minus = frozenset(an_y - de_x)

    all_nodes = set(dag.nodes)
    non_desc_x = sorted(all_nodes - dag.descendants(x) - {x})
    plus_candidates = [v for v in non_desc_x if v not in o_true and v not in (x, y)]
    plus_sets = [frozenset(o_true | {v}) for v in plus_candidates[:MENU_PERTURB_CAP]]

    o_sorted = sorted(o_true)
    minus_sets = [frozenset(o_true - {v}) for v in o_sorted[:MENU_PERTURB_CAP]]

    raw: list[tuple[str, frozenset]] = [
        ("O_true", o_true),
        ("empty", frozenset()),
        ("pa_true", pa_true),
        ("an_minus_de", an_minus),
    ]
    for i, s in enumerate(plus_sets):
        raw.append((f"O_plus_v{i}", s))
    for i, s in enumerate(minus_sets):
        raw.append((f"O_minus_v{i}", s))

    grouped: dict[frozenset, list[str]] = {}
    for src, s in raw:
        grouped.setdefault(s, []).append(src)

    canon_order = sorted(grouped.keys(), key=lambda fs: tuple(sorted(fs)))
    canon_order = canon_order[:MENU_CAP]
    return [MenuEntry(source=";".join(sorted(set(grouped[s]))), members=s) for s in canon_order]


def menu_rows(base_id: str, origin: str, network: str, x: Node, y: Node, menu: list[MenuEntry]) -> list[dict[str, Any]]:
    rows = []
    for i, entry in enumerate(menu):
        rows.append(
            {
                "origin": origin,
                "base_id": base_id,
                "network": network,
                "treatment": x,
                "outcome": y,
                "menu_index": i,
                "source": entry.source,
                "members": ";".join(sorted(entry.members)),
                "size": len(entry.members),
            }
        )
    return rows


# ---------------------------------------------------------------------------
# Task A -- per-proposal measurement over the menu
# ---------------------------------------------------------------------------


def _safe_bias(sem: LinearSEM, x: Node, y: Node, z: frozenset[Node]) -> float | None:
    try:
        b = bias(sem, x, y, z)
    except (FloatingPointError, ZeroDivisionError, np.linalg.LinAlgError):
        return None
    if not math.isfinite(b):
        return None
    return b


def measure_proposal_menu(
    base: P.BaseCpdag,
    menu: list[MenuEntry],
    proposal: P.Proposal,
    sems: list[P.SemDraw],
    root_seed: int,
    origin: str,
    network: str,
    radius_time_limit_s: float = P.RADIUS_TIME_LIMIT_S,
) -> dict[str, Any]:
    """Measure one proposal's achievable efficiency over the fixed menu.

    Returns a single row. As a side effect (recorded in the row itself as
    ``max_abs_bias_check``), every menu member found GAC-valid on this
    proposal's ``G0`` is bias-checked against every SEM draw -- this is the
    per-row contribution to acceptance criterion C2 (menu-validity/bias
    wiring sanity): every proposal here asserts only truthful orientations
    (round 1's proposal enumerator, reused unmodified), so the true DAG is
    always a member of ``G0``'s extension set, and GAC-validity on ``G0``
    therefore implies unbiasedness against the true DAG's SEM.
    """
    x, y = base.treatment, base.outcome
    row: dict[str, Any] = {
        "origin": origin,
        "network": network,
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
        "orientations": P.encode_orientations(proposal.orientations),
        "gate": GATE_NAME,
        "status": None,
        "n_menu_members": len(menu),
        "n_valid_menu_members": 0,
        "n_sem_valid": 0,
        "n_sem_total": len(sems),
        "achievable_utility_median": "",
        "achievable_utility_q25": "",
        "achievable_utility_q75": "",
        "winning_member_index": "",
        "winning_member_source": "",
        "winning_members_encoded": "",
        "representative_sem_k": "",
        "assumes": "",
        "method": "",
        "oracle": "",
        "exact": "",
        "radius_timed_out": "",
        "r_val": "",
        "max_abs_bias_check": "",
        "front_utility": False,
    }

    g0p = apply_orientations(base.cpdag_obj, proposal.orientations)
    if g0p is None:
        row["status"] = "proposal_contradictory"
        return row

    valid_idx = [i for i, m in enumerate(menu) if is_gac_valid_mpdag(g0p, x, y, m.members)]
    row["n_valid_menu_members"] = len(valid_idx)

    max_abs_bias = 0.0
    any_bias_checked = False
    for i in valid_idx:
        members = menu[i].members
        for d in sems:
            b = _safe_bias(d.sem, x, y, members)
            if b is not None:
                any_bias_checked = True
                max_abs_bias = max(max_abs_bias, abs(b))
    row["max_abs_bias_check"] = max_abs_bias if any_bias_checked else ""

    if not valid_idx:
        row["status"] = "no_menu_member_valid"
        return row

    draw_results: list[tuple[int, float, int]] = []  # (sem_k, utility, winning_menu_index)
    for d in sems:
        avars: list[tuple[float, int]] = []
        for i in valid_idx:
            av = P._safe_avar(d.sem, x, y, menu[i].members)
            if av is not None:
                avars.append((av, i))
        if not avars:
            continue
        avars.sort(key=lambda t: (t[0], t[1]))
        best_av, best_i = avars[0]
        draw_results.append((d.k, 1.0 / best_av, best_i))

    row["n_sem_valid"] = len(draw_results)
    if not draw_results:
        row["status"] = "avar_undefined"
        return row

    utilities = np.array([u for _, u, _ in draw_results], dtype=float)
    median_u = float(np.median(utilities))
    # Representative draw for winning_member_index / r_val: the draw whose
    # utility is closest to the reported median, ties broken by smallest
    # sem_k. Documented in PARETO2_RUN_NOTES.md: menu-validity is structural
    # (SEM-independent), only *which* valid member wins a given draw depends
    # on that draw's random weights, so a single representative committed set
    # is chosen this way for the radius measurement, which needs one fixed Z.
    rep_k, rep_u, rep_i = sorted(draw_results, key=lambda t: (abs(t[1] - median_u), t[0]))[0]

    row["status"] = "ok"
    row["achievable_utility_median"] = median_u
    row["achievable_utility_q25"] = float(np.percentile(utilities, 25))
    row["achievable_utility_q75"] = float(np.percentile(utilities, 75))
    row["winning_member_index"] = rep_i
    row["winning_member_source"] = menu[rep_i].source
    row["winning_members_encoded"] = ";".join(sorted(menu[rep_i].members))
    row["representative_sem_k"] = rep_k

    z_committed = frozenset(menu[rep_i].members)
    result = breakdown_radius(
        base.cpdag_obj,
        list(proposal.orientations),
        x,
        y,
        z_committed,
        g0=g0p,
        time_limit_s=radius_time_limit_s,
    )
    row["assumes"] = result.assumes
    row["method"] = result.method
    row["oracle"] = result.oracle
    row["exact"] = bool(result.exact)
    row["radius_timed_out"] = not result.exact
    row["r_val"] = int(result.radius)

    return row


def compute_front(rows: list[dict[str, Any]]) -> None:
    """Mutates ``rows`` in place, setting ``front_utility`` on the
    non-dominated (maximize ``r_val``, maximize ``achievable_utility_median``)
    subset of the "ok"-status rows of a single stratum. ``rows`` must all
    share one ``base_id``. Mirrors
    :func:`bkrobust.robustness.pareto.compute_fronts`'s utility front exactly
    (same rank convention for ``UNREACHED``, same non-domination test),
    reused via the private helpers rather than duplicated, since both are
    pure functions with no state.
    """
    points = []
    for r in rows:
        if r["status"] != "ok":
            continue
        rr = P._r_rank(r["r_val"])
        if r["achievable_utility_median"] != "":
            points.append((r["proposal_id"], rr, float(r["achievable_utility_median"])))
    front = P._non_dominated(points)
    for r in rows:
        if r["proposal_id"] in front:
            r["front_utility"] = True


# ---------------------------------------------------------------------------
# Task A -- one stratum end to end
# ---------------------------------------------------------------------------


def run_one_base_task_a(base: P.BaseCpdag, root_seed: int) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """Reuses round 1's SEM draws, truthful-orientation menu, and proposal
    sampler unmodified; only the per-proposal measurement (menu-based
    achievable efficiency, not the MPDAG-invariant optimal set) is new.

    Returns (proposal_rows, menu_rows, info).
    """
    import time

    t0 = time.perf_counter()
    x, y = base.treatment, base.outcome
    sems = P.draw_sems(base, root_seed, n=P.N_SEM_DRAWS)

    menu = build_menu(base.dag_obj, x, y)
    mrows = menu_rows(base.base_id, "synthetic", "", x, y, menu)

    truth_menu = P.truthful_orientation_menu(base)
    desc_x = base.dag_obj.descendants(x) | {x}
    anc_y = base.dag_obj.ancestors(y) | {y}
    proposals = P.sample_proposals(base, truth_menu, desc_x, anc_y, root_seed=root_seed, cap=P.PROPOSAL_CAP)

    prows: list[dict[str, Any]] = []
    for prop in proposals:
        row = measure_proposal_menu(base, menu, prop, sems, root_seed, origin="synthetic", network="")
        prows.append(row)

    ok_rows = [r for r in prows if r["status"] == "ok"]
    distinct_utility = {r["achievable_utility_median"] for r in ok_rows}
    distinct_winner = {r["winning_member_index"] for r in ok_rows}
    utility_varies = len(distinct_utility) > 1

    if utility_varies:
        compute_front(prows)

    elapsed = time.perf_counter() - t0
    info = {
        "base_id": base.base_id,
        "n_menu_members": len(menu),
        "n_undirected_cpdag": len(truth_menu),
        "n_proposals": len(proposals),
        "n_ok": len(ok_rows),
        "n_no_menu_member_valid": sum(1 for r in prows if r["status"] == "no_menu_member_valid"),
        "n_avar_undefined": sum(1 for r in prows if r["status"] == "avar_undefined"),
        "n_proposal_contradictory": sum(1 for r in prows if r["status"] == "proposal_contradictory"),
        "n_distinct_achievable_utility": len(distinct_utility),
        "n_distinct_winning_member": len(distinct_winner),
        "utility_varies": utility_varies,
        "max_abs_bias_check": max(
            (r["max_abs_bias_check"] for r in prows if r["max_abs_bias_check"] != ""), default=0.0
        ),
        "elapsed_s": elapsed,
    }
    return prows, mrows, info


# ---------------------------------------------------------------------------
# Task B -- real benchmark networks: loading
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class NetworkData:
    name: str
    source_file: str
    sha256: str
    source_format: str
    dag: MPDAG
    cpdag: MPDAG


def sha256_of_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def discover_network_files() -> list[Path]:
    exts = (".bif.gz", ".bif", ".txt", ".json")
    files = [p for p in NETWORKS_DIR.iterdir() if p.name.endswith(exts)]
    return sorted(files, key=lambda p: p.name)


def load_networks() -> tuple[list[NetworkData], list[dict[str, Any]]]:
    """Parses every benchmark file, excluding ``M-bias`` (an ADMG, not a DAG,
    per the task brief -- 39 networks, not 40) and counting every other
    exclusion out loud (C5): a parse failure, a graph with bidirected edges
    or a directed cycle (neither is a DAG), or a ``dag_to_cpdag`` failure.

    Returns (loaded networks, exclusion rows).
    """
    loaded: list[NetworkData] = []
    exclusions: list[dict[str, Any]] = []
    for path in discover_network_files():
        if path.name in EXCLUDED_NETWORK_FILES:
            exclusions.append({"file": path.name, "reason": "excluded_admg_not_dag_by_design"})
            continue
        sha = sha256_of_file(path)
        try:
            network = parse_file(path, sha)
        except ValueError as exc:
            exclusions.append({"file": path.name, "reason": f"parse_failed: {exc}"})
            continue
        if network.bidirected:
            exclusions.append({"file": path.name, "reason": "bidirected_edges_not_a_dag"})
            continue
        dag = to_mpdag(network)
        if not dag.is_dag():
            exclusions.append({"file": path.name, "reason": "not_a_dag_has_cycle"})
            continue
        try:
            cpdag = dag_to_cpdag(dag)
        except ValueError as exc:
            exclusions.append({"file": path.name, "reason": f"cpdag_failed: {exc}"})
            continue
        loaded.append(
            NetworkData(
                name=network.name,
                source_file=network.source_file,
                sha256=network.sha256,
                source_format=network.source_format,
                dag=dag,
                cpdag=cpdag,
            )
        )
    return loaded, exclusions


# ---------------------------------------------------------------------------
# Task B -- census (structural only, never calls fast_gate)
# ---------------------------------------------------------------------------


def _sample_ordered_pairs(nodes: list[Node], root_seed: int, network_name: str, cap: int = CENSUS_PAIR_CAP) -> list[Edge]:
    n = len(nodes)
    total_possible = n * (n - 1)
    if total_possible <= cap:
        return sorted((a, b) for a in nodes for b in nodes if a != b)
    rng = np.random.default_rng(P.derive_seed(root_seed, "census_pairs", network_name))
    seen: set[Edge] = set()
    max_attempts = cap * 50 + 2000
    attempts = 0
    while len(seen) < cap and attempts < max_attempts:
        attempts += 1
        i = int(rng.integers(0, n))
        j = int(rng.integers(0, n))
        if i == j:
            continue
        seen.add((nodes[i], nodes[j]))
    return sorted(seen)


def census_one_network(net: NetworkData, root_seed: int) -> tuple[dict[str, Any], list[Edge]]:
    """Structural-only census: no ``fast_gate`` call (that would cost ~29s/pair
    at 109 nodes and dominate). Counts, among up to :data:`CENSUS_PAIR_CAP`
    ordered pairs with a causal path, how many have some CPDAG-undirected
    edge lying on a directed treatment-to-outcome path.

    Returns (row, on_path_pairs) -- the second element feeds the follow-up
    stage and is not written to the census CSV (which is one row per
    network, per the task brief).
    """
    dag, cpdag = net.dag, net.cpdag
    nodes = list(dag.nodes)
    candidates = _sample_ordered_pairs(nodes, root_seed, net.name)

    oriented_undirected: list[Edge] = []
    for a, b in sorted(cpdag.undirected_edges):
        if dag.is_directed_edge(a, b):
            oriented_undirected.append((a, b))
        elif dag.is_directed_edge(b, a):
            oriented_undirected.append((b, a))
        else:
            raise AssertionError(f"{net.name}: true DAG does not orient CPDAG edge {(a, b)}")

    kept: list[Edge] = []
    on_path_pairs: list[Edge] = []
    for x, y in candidates:
        if y not in dag.descendants(x):
            continue
        kept.append((x, y))
        desc_x = dag.descendants(x) | {x}
        anc_y = dag.ancestors(y) | {y}
        if any(tail in desc_x and head in anc_y for tail, head in oriented_undirected):
            on_path_pairs.append((x, y))

    row = {
        "network": net.name,
        "source_file": net.source_file,
        "sha256": net.sha256,
        "source_format": net.source_format,
        "n_nodes": len(nodes),
        "n_undirected_edges": len(cpdag.undirected_edges),
        "n_candidate_pairs_drawn": len(candidates),
        "n_pairs_sampled": len(kept),
        "n_pairs_with_on_path_undirected_edge": len(on_path_pairs),
    }
    return row, sorted(on_path_pairs)


def run_census(networks: list[NetworkData], root_seed: int) -> tuple[list[dict[str, Any]], dict[str, list[Edge]]]:
    rows: list[dict[str, Any]] = []
    on_path_by_network: dict[str, list[Edge]] = {}
    for net in sorted(networks, key=lambda n: n.name):
        row, on_path_pairs = census_one_network(net, root_seed)
        rows.append(row)
        on_path_by_network[net.name] = on_path_pairs
    return rows, on_path_by_network


# ---------------------------------------------------------------------------
# Task B -- follow-up: gate + measure on the smallest on-path-bearing networks
# ---------------------------------------------------------------------------


def select_followup_targets(
    census_rows: list[dict[str, Any]],
    on_path_by_network: dict[str, list[Edge]],
    max_networks: int = FOLLOWUP_MAX_NETWORKS,
) -> list[tuple[str, list[Edge]]]:
    """The smallest (by ``n_nodes``) networks with at least one on-path pair,
    per the task brief ("pick 3-6 of the SMALLEST such networks")."""
    positive = [r for r in census_rows if r["n_pairs_with_on_path_undirected_edge"] > 0]
    positive_sorted = sorted(positive, key=lambda r: (r["n_nodes"], r["network"]))
    chosen = positive_sorted[:max_networks]
    return [(r["network"], on_path_by_network[r["network"]]) for r in chosen]


def run_followup(
    networks_by_name: dict[str, NetworkData],
    targets: list[tuple[str, list[Edge]]],
    root_seed: int,
    pairs_per_network: int = FOLLOWUP_PAIRS_PER_NETWORK,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """Gates the selected (network, x, y) pairs with ``fast_gate`` and, for
    every admissible one, runs the Task A menu-based measurement.

    Returns (proposal_rows, menu_rows, gate_exclusion_rows, info).
    """
    import time

    t0 = time.perf_counter()
    prows: list[dict[str, Any]] = []
    mrows: list[dict[str, Any]] = []
    gate_exclusions: list[dict[str, Any]] = []
    n_bases_measured = 0

    for net_name, on_path_pairs in targets:
        net = networks_by_name[net_name]
        for x, y in on_path_pairs[:pairs_per_network]:
            ok, reason = fast_gate(net.dag, net.cpdag, x, y)
            if not ok:
                gate_exclusions.append({"network": net_name, "treatment": x, "outcome": y, "reason": reason})
                continue
            base_id = f"real__{net_name}__{x}__{y}"
            base = P.BaseCpdag(
                base_id=base_id,
                component_size=-1,
                separation=-1,
                seed=P.derive_seed(root_seed, "real_base", net_name, x, y),
                treatment=x,
                outcome=y,
                dag_obj=net.dag,
                cpdag_obj=net.cpdag,
            )
            sems = P.draw_sems(base, root_seed, n=P.N_SEM_DRAWS)
            menu = build_menu(base.dag_obj, x, y)
            mrows.extend(menu_rows(base_id, "real_network", net_name, x, y, menu))

            truth_menu = P.truthful_orientation_menu(base)
            desc_x = base.dag_obj.descendants(x) | {x}
            anc_y = base.dag_obj.ancestors(y) | {y}
            proposals = P.sample_proposals(base, truth_menu, desc_x, anc_y, root_seed=root_seed, cap=P.PROPOSAL_CAP)

            base_rows = []
            for prop in proposals:
                row = measure_proposal_menu(
                    base, menu, prop, sems, root_seed, origin="real_network", network=net_name
                )
                base_rows.append(row)

            ok_rows = [r for r in base_rows if r["status"] == "ok"]
            if len({r["achievable_utility_median"] for r in ok_rows}) > 1:
                compute_front(base_rows)

            prows.extend(base_rows)
            n_bases_measured += 1

    elapsed = time.perf_counter() - t0
    info = {
        "n_targets_considered": sum(len(pairs[:pairs_per_network]) for _, pairs in targets),
        "n_gate_excluded": len(gate_exclusions),
        "n_bases_measured": n_bases_measured,
        "n_proposal_rows": len(prows),
        "elapsed_s": elapsed,
    }
    return prows, mrows, gate_exclusions, info
