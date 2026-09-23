"""Tests for the designed component generator.

Three things are asserted here, in order of how badly a regression in each would
hurt the study:

1. **Determinism** -- the same ``(spec, seed)`` gives a bit-identical instance,
   including from a fresh interpreter under a different ``PYTHONHASHSEED``.
2. **Realised parameters are measured, not assumed** -- the realised block is
   recomputed from the constructed CPDAG by an independent BFS written here, and
   ``measure_realised`` is fed deliberately mismatched inputs to confirm it
   reports the graph in front of it rather than the request that produced it.
3. **A generated instance passes the degeneracy gate** -- the shared
   :func:`bkrobust.synth.runner.gate`, not a local reimplementation of it.

Everything is bounded: no test sweeps more than a couple of dozen small draws.
"""

from __future__ import annotations

import os
import subprocess
import sys
from collections import deque
from pathlib import Path

import numpy as np
import pytest

from bkrobust.demo.graph import MPDAG, undirected_components
from bkrobust.synth.component_generator import (
    OUTCOME,
    SPECTATOR,
    TREATMENT,
    ComponentSpec,
    achievable_separations,
    build_component_dag,
    component_containing,
    generate_instance,
    measure_realised,
    run_pilot,
)
from bkrobust.synth.runner import gate

REPO_ROOT = Path(__file__).resolve().parents[2]


# --- independent measurement helpers (deliberately not the module's own) ----


def independent_component_of(cpdag: MPDAG, node: str) -> set[str]:
    """Flood-fill the undirected neighbourhood of ``node``, written from scratch."""
    seen = {node}
    queue = deque([node])
    while queue:
        current = queue.popleft()
        for a, b in sorted(cpdag.undirected_edges):
            other = None
            if a == current:
                other = b
            elif b == current:
                other = a
            if other is not None and other not in seen:
                seen.add(other)
                queue.append(other)
    return seen


def independent_distance(cpdag: MPDAG, source: str, target: str) -> int | None:
    """Shortest undirected-edge path length from ``source`` to ``target``, or None."""
    dist = {source: 0}
    queue = deque([source])
    while queue:
        current = queue.popleft()
        if current == target:
            return dist[current]
        for a, b in sorted(cpdag.undirected_edges):
            for here, other in ((a, b), (b, a)):
                if here == current and other not in dist:
                    dist[other] = dist[current] + 1
                    queue.append(other)
    return dist.get(target)


SMALL_CELLS = [(2, 1), (4, 2), (6, 3), (6, 5), (8, 4), (9, 8)]


# --- determinism -----------------------------------------------------------


def test_same_seed_gives_identical_instance():
    spec = ComponentSpec(component_size=8, separation=4)
    first = generate_instance(spec, 12345)
    second = generate_instance(spec, 12345)
    assert first.fingerprint() == second.fingerprint()
    assert first.to_row() == second.to_row()


def test_different_seeds_explore_different_components():
    spec = ComponentSpec(component_size=9, separation=3)
    prints = {generate_instance(spec, seed).fingerprint() for seed in range(12)}
    assert len(prints) > 1, "the generator is not actually using its seed"


@pytest.mark.parametrize("hash_seed", ["0", "1", "12345"])
def test_fingerprint_is_stable_across_pythonhashseed(hash_seed):
    spec = ComponentSpec(component_size=7, separation=3)
    expected = generate_instance(spec, 2024).fingerprint()

    env = os.environ.copy()
    env["PYTHONHASHSEED"] = hash_seed
    existing = env.get("PYTHONPATH", "")
    src = str(REPO_ROOT / "src")
    env["PYTHONPATH"] = src if not existing else f"{src}{os.pathsep}{existing}"
    script = (
        "from bkrobust.synth.component_generator import ComponentSpec, generate_instance;"
        "print(generate_instance("
        "ComponentSpec(component_size=7, separation=3), 2024).fingerprint())"
    )
    completed = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        env=env,
        cwd=str(REPO_ROOT),
        check=True,
    )
    assert completed.stdout.strip() == expected


# --- realised parameters are measured, not assumed -------------------------


@pytest.mark.parametrize(("size", "separation"), SMALL_CELLS)
def test_realised_matches_an_independent_measurement(size, separation):
    spec = ComponentSpec(component_size=size, separation=separation)
    inst = generate_instance(spec, 77)
    cpdag = inst.cpdag_obj
    assert cpdag is not None

    component = independent_component_of(cpdag, TREATMENT)
    assert inst.realised["realised_component_size"] == len(component)

    optimal = inst.optimal_set or ()
    in_component = [v for v in optimal if v in component]
    assert sorted(in_component) == inst.realised["o_members_in_component"]
    expected = min(
        d
        for d in (independent_distance(cpdag, TREATMENT, v) for v in in_component)
        if d is not None
    )
    assert inst.realised["realised_separation"] == expected


def test_measure_realised_reads_the_graph_not_the_request():
    """Feed one instance's CPDAG the *other* instance's parameters."""
    short = generate_instance(ComponentSpec(component_size=7, separation=2), 5)
    long = generate_instance(ComponentSpec(component_size=7, separation=6), 5)
    assert short.cpdag_obj is not None and long.cpdag_obj is not None

    # Measuring the long instance's graph must report 6, even though everything
    # about the call is otherwise borrowed from the short instance's request.
    measured = measure_realised(
        long.dag_obj,
        long.cpdag_obj,
        long.g0_obj,
        set(long.optimal_set or ()),
        short.k_g0,
        TREATMENT,
        OUTCOME,
    )
    assert measured["realised_separation"] == 6
    assert measured["realised_separation_status"] == "measured"
    assert measured["realised_k_g0"] == len(short.k_g0)


def test_missing_optimal_member_uses_a_sentinel_not_a_number():
    inst = generate_instance(ComponentSpec(component_size=6, separation=3), 11)
    assert inst.cpdag_obj is not None
    measured = measure_realised(
        inst.dag_obj,
        inst.cpdag_obj,
        inst.g0_obj,
        {SPECTATOR},  # W is never inside the component
        inst.k_g0,
        TREATMENT,
        OUTCOME,
    )
    assert measured["realised_separation"] is None
    assert measured["realised_separation_status"] == "no_optimal_member_in_component"
    assert measured["o_members_in_component"] == []

    unidentified = measure_realised(
        inst.dag_obj, inst.cpdag_obj, inst.g0_obj, None, inst.k_g0, TREATMENT, OUTCOME
    )
    assert unidentified["realised_separation"] is None
    assert unidentified["realised_separation_status"] == "optimal_set_not_identified"
    assert unidentified["optimal_set_identified"] is False


# --- the gate --------------------------------------------------------------


@pytest.mark.parametrize(("size", "separation"), SMALL_CELLS)
def test_generated_instance_passes_the_shared_gate(size, separation):
    spec = ComponentSpec(component_size=size, separation=separation)
    inst = generate_instance(spec, 3)
    assert inst.accepted, inst.reject_reason
    assert inst.dag_obj is not None and inst.cpdag_obj is not None
    passed, reason = gate(inst.dag_obj, inst.cpdag_obj, TREATMENT, OUTCOME)
    assert passed, reason
    assert inst.realised["x_in_component"] is True
    assert inst.realised["y_descendant_of_x"] is True
    assert inst.realised["realised_component_size"] == size
    assert inst.realised["realised_separation"] == separation


# --- structural invariants the construction promises -----------------------


@pytest.mark.parametrize(("size", "separation"), SMALL_CELLS)
def test_component_survives_into_the_cpdag_fully_undirected(size, separation):
    spec = ComponentSpec(component_size=size, separation=separation)
    dag, layout = build_component_dag(spec, np.random.default_rng(4))
    inst = generate_instance(spec, 4)
    cpdag = inst.cpdag_obj
    assert cpdag is not None
    assert dag.edge_string() == inst.true_dag

    component = component_containing(cpdag, TREATMENT)
    assert component is not None
    assert set(component) == set(layout.labels)
    # every component edge came back undirected
    for tail, head in layout.component_edges:
        assert cpdag.is_undirected_edge(tail, head)
    # exactly one undirected component: nothing else leaked in
    assert len(undirected_components(cpdag)) == 1
    # the edges into the outcome are compelled, and W touches only the outcome
    assert cpdag.is_directed_edge(TREATMENT, OUTCOME)
    assert cpdag.is_directed_edge(layout.anchor, OUTCOME)
    assert cpdag.is_directed_edge(SPECTATOR, OUTCOME)
    assert cpdag.adjacent(SPECTATOR) == {OUTCOME}


def test_spec_rejects_unachievable_parameters():
    assert achievable_separations(2) == (1,)
    assert achievable_separations(5) == (1, 2, 3, 4)
    with pytest.raises(ValueError):
        achievable_separations(1)
    with pytest.raises(ValueError):
        ComponentSpec(component_size=4, separation=4)
    with pytest.raises(ValueError):
        ComponentSpec(component_size=4, separation=1, coverage=1.5)
    with pytest.raises(ValueError):
        ComponentSpec(component_size=4, separation=1, coverage_order="whatever")


# --- the pilot report ------------------------------------------------------


def test_run_pilot_is_bounded_and_reports_per_cell():
    report = run_pilot(
        component_sizes=(3, 4),
        budget=2,
        root_seed=1,
        coverage_sweep=(1.0,),
        coverage_sweep_size=4,
        coverage_sweep_separation=2,
        coverage_sweep_budget=1,
        radius_spotcheck_max_size=4,
    )
    cells = {(c["component_size"], c["separation"]) for c in report["cells"]}
    assert cells == {(3, 1), (3, 2), (4, 1), (4, 2), (4, 3)}
    assert report["totals"]["attempts"] == 10
    for cell in report["cells"]:
        assert cell["attempts"] == 2
        assert 0.0 <= cell["acceptance_rate"] <= 1.0
        assert "rejection_reasons" in cell
    assert report["separation_range"]["max_realised_separation_overall"] == 3
    assert report["determinism"]["digest"]
