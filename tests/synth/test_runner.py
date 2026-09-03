"""Tests for bkrobust.synth.runner.

What each test discriminates:

* ``gate`` tests use small hand-built graphs (verified by hand in the
  docstrings/comments here) to pin down each of the three most reachable
  rejection reasons, plus one accepted case -- rather than trusting the gate
  only from its behaviour on random generator output.
* ``run_instance`` tests check a non-degenerate accepted instance really has
  ``r_val >= 1`` (per the radius convention: Z is read off G0, so it is valid
  there by construction) and that a genuinely degenerate hand-built instance
  raises ``GateRejectedError`` rather than silently returning something.
* ``run_grid`` tests check resume-safety directly: run a grid, simulate a
  kill by truncating the results file after some rows, resume, and check no
  row is duplicated and the run completes.
* The mandatory PYTHONHASHSEED-invariance test covers ``run_instance``,
  since it is the function whose RNG-consumption order is most at risk (it
  threads one generator through graph generation, knowledge simulation and
  bias sampling in sequence).
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from bkrobust.core.conventions import UNREACHED
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.synth.runner import (
    MAX_UNDIRECTED_EDGES,
    GateRejectedError,
    gate,
    run_grid,
    run_instance,
)

# --- gate: hand-built cases -----------------------------------------------


def test_gate_rejects_treatment_not_adjacent_to_component():
    """Two confounders both point into X and Y with no other ambiguity at all.

    Z1, Z2 -> X, Y each triggers a v-structure at X and at Y (Z1, Z2 are not
    adjacent), so the whole CPDAG ends up fully directed: no undirected edge
    anywhere for X to be in or next to.
    """
    dag = MPDAG(
        ["X", "Y", "Z1", "Z2"],
        directed=[("Z1", "X"), ("Z1", "Y"), ("Z2", "X"), ("Z2", "Y")],
    )
    cpdag = dag_to_cpdag(dag)
    assert cpdag.undirected_edges == frozenset()
    assert gate(dag, cpdag, "X", "Y") == (False, "treatment_not_in_or_adjacent_to_component")


def test_gate_rejects_empty_set_trivially_valid():
    """X -> Y with no other nodes at all: there is no confounding to speak of."""
    dag = MPDAG(["X", "Y"], directed=[("X", "Y")])
    cpdag = dag_to_cpdag(dag)
    assert gate(dag, cpdag, "X", "Y") == (False, "empty_set_trivially_valid")


def test_gate_rejects_no_atomic_perturbation_case():
    """A -> B is the only undirected CPDAG edge, and B is irrelevant to X, Y.

    True DAG: A->B, A->X, A->Y, C->X, X->Y. A and C are non-adjacent parents
    of X, so A->X and C->X are compelled by a v-structure; A and X are
    adjacent, so A->Y is not compelled by one there. The only edge left
    undirected in the CPDAG is A-B (verified below) -- and since B does not
    touch X, Y, or any valid adjustment set, dropping that one known claim
    changes nothing about validity of any candidate Z: the sanity gate must
    reject it.
    """
    dag = MPDAG(
        ["X", "Y", "A", "B", "C"],
        directed=[("A", "B"), ("A", "X"), ("A", "Y"), ("C", "X"), ("X", "Y")],
    )
    cpdag = dag_to_cpdag(dag)
    assert cpdag.undirected_edges == frozenset({("A", "B")})
    assert gate(dag, cpdag, "X", "Y") == (False, "no_atomic_perturbation_changes_validity")


def test_gate_accepts_a_genuine_non_degenerate_case():
    """Z -> X, Z -> Y with Z's own ambiguous parent: the textbook non-degenerate shape."""
    dag = MPDAG(
        ["X", "Y", "Z", "W"],
        directed=[("W", "Z"), ("Z", "X"), ("Z", "Y")],
    )
    cpdag = dag_to_cpdag(dag)
    assert gate(dag, cpdag, "X", "Y") == (True, "ok")


# --- run_instance -----------------------------------------------------


def test_run_instance_r_val_at_least_one_on_nondegenerate_case():
    """The radius convention: Z is valid at G0 by construction, so r_val >= 1."""
    accepted = False
    for seed in range(50):
        try:
            inst = run_instance(
                instance_id=f"probe_{seed}",
                generator_name="decoupled_backdoor",
                generator_params={"coupling": 0.0},
                n=6,
                seed=seed,
                knows_fraction=0.8,
                corruption_name="none",
                corruption_rate=0.0,
            )
        except GateRejectedError:
            continue
        accepted = True
        assert inst.r_val >= 1
        assert inst.r_val != UNREACHED or inst.r_val >= 1  # UNREACHED == -1, excluded either way
        assert inst.method == "bfs_exact"
        assert inst.z  # decoupled_backdoor's optimal set is never empty
        assert set(inst.timings) >= {"graph", "gate", "knowledge", "g0_and_z", "space"}
        break
    assert accepted, "expected at least one seed to clear the gate"


def test_run_instance_r_opt_and_r_eps_are_populated():
    for seed in range(50):
        try:
            inst = run_instance(
                instance_id=f"probe2_{seed}",
                generator_name="decoupled_backdoor",
                generator_params={"coupling": 0.0},
                n=6,
                seed=seed,
                knows_fraction=0.8,
                corruption_name="none",
                corruption_rate=0.0,
                epsilons=(0.05, 0.2),
            )
        except GateRejectedError:
            continue
        assert inst.r_opt >= 1
        assert set(inst.r_eps) == {"0.05", "0.2"}
        for r in inst.r_eps.values():
            assert r == UNREACHED or r >= 1
        return
    pytest.fail("expected at least one seed to clear the gate")


def test_run_instance_raises_gate_rejected_on_degenerate_generator_params():
    """decoupled_backdoor at coupling=1.0, n=8 always exceeds MAX_UNDIRECTED_EDGES."""
    with pytest.raises(GateRejectedError) as excinfo:
        run_instance(
            instance_id="degenerate",
            generator_name="decoupled_backdoor",
            generator_params={"coupling": 1.0},
            n=8,
            seed=0,
            knows_fraction=0.8,
            corruption_name="none",
            corruption_rate=0.0,
        )
    assert excinfo.value.reason == "cpdag_too_large_for_bfs"


def test_run_instance_rejects_unknown_generator():
    with pytest.raises(ValueError):
        run_instance(
            instance_id="bad",
            generator_name="not_a_real_generator",
            generator_params={},
            n=6,
            seed=0,
            knows_fraction=0.5,
            corruption_name="none",
            corruption_rate=0.0,
        )


def test_run_instance_tiered_corruption_path():
    """corruption_name='tiered' takes the generative path (no k_true corruption)."""
    for seed in range(50):
        try:
            inst = run_instance(
                instance_id=f"tiered_{seed}",
                generator_name="decoupled_backdoor",
                generator_params={"coupling": 0.0},
                n=6,
                seed=seed,
                knows_fraction=0.8,
                corruption_name="tiered",
                corruption_rate=0.0,
                tier_params={"n_tiers": 5, "corruption_rate": 0.0},
            )
        except GateRejectedError:
            continue
        assert inst.corruption == "tiered"
        return
    pytest.fail("expected at least one seed to clear the gate with tiered corruption")


def _non_timing_row(instance) -> dict:
    return {k: v for k, v in instance.to_row().items() if not k.startswith("time_")}


def test_run_instance_reproducible_from_seed():
    for seed in range(50):
        kwargs = dict(
            instance_id="repro",
            generator_name="decoupled_backdoor",
            generator_params={"coupling": 0.0},
            n=6,
            seed=seed,
            knows_fraction=0.8,
            corruption_name="flip",
            corruption_rate=0.3,
        )
        try:
            a = run_instance(**kwargs)
            b = run_instance(**kwargs)
        except GateRejectedError:
            continue
        # Wall-clock timings legitimately differ between the two calls; every
        # other field must match exactly.
        assert _non_timing_row(a) == _non_timing_row(b)
        return
    pytest.fail("expected at least one seed to clear the gate")


MAX_UNDIRECTED_EDGES_SENTINEL = MAX_UNDIRECTED_EDGES  # touch the constant so linters see it used


def test_max_undirected_edges_is_a_positive_int():
    assert isinstance(MAX_UNDIRECTED_EDGES_SENTINEL, int)
    assert MAX_UNDIRECTED_EDGES_SENTINEL > 0


# --- run_grid: resume safety -------------------------------------------


def _small_grid() -> list[dict]:
    """A small mixed grid: some points accept, some are gate-rejected -- both matter for resume."""
    grid = []
    for seed_offset in range(6):
        grid.append(
            {
                "generator_name": "decoupled_backdoor",
                "generator_params": {"coupling": 0.0},
                "n": 6,
                "knows_fraction": 0.6 + 0.05 * seed_offset,
                "corruption_name": "flip",
                "corruption_rate": 0.2,
            }
        )
    # A couple of guaranteed-degenerate points, to exercise rejection bookkeeping.
    grid.append(
        {
            "generator_name": "decoupled_backdoor",
            "generator_params": {"coupling": 1.0},
            "n": 8,
            "knows_fraction": 0.5,
            "corruption_name": "none",
            "corruption_rate": 0.0,
        }
    )
    grid.append(
        {
            "generator_name": "erdos_renyi",
            "generator_params": {"edge_prob": 0.0},
            "n": 6,
            "knows_fraction": 0.5,
            "corruption_name": "none",
            "corruption_rate": 0.0,
        }
    )
    return grid


def test_run_grid_writes_rows_and_manifest(tmp_path: Path):
    grid = _small_grid()
    summary = run_grid(grid, tmp_path, root_seed=100, epsilons=(0.1,), n_bias_draws=5)

    assert summary["n_total"] == len(grid)
    assert summary["n_accepted"] + sum(summary["rejection_counts"].values()) == len(grid)

    csv_path = tmp_path / "results.csv"
    assert csv_path.exists()
    with csv_path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == summary["n_accepted"]

    manifest_path = tmp_path / "manifest.json"
    assert manifest_path.exists()
    import json

    manifest = json.loads(manifest_path.read_text())
    assert manifest["seed"] == 100
    assert "radius_convention" in manifest
    assert manifest["rejection_counts"] == summary["rejection_counts"]


def test_run_grid_resume_survives_a_simulated_kill(tmp_path: Path):
    """Kill-and-restart must not duplicate rows and must reach the same final state."""
    grid = _small_grid()

    # Full run, uninterrupted, as the ground truth to compare against.
    ground_truth_dir = tmp_path / "full"
    full_summary = run_grid(grid, ground_truth_dir, root_seed=200, epsilons=(0.1,), n_bias_draws=5)

    # Simulated kill: run once, then truncate results.csv to simulate the
    # process dying after some rows were flushed but before the run finished
    # bookkeeping (a manifest write, say). Then resume.
    killed_dir = tmp_path / "killed"
    run_grid(grid, killed_dir, root_seed=200, epsilons=(0.1,), n_bias_draws=5)
    csv_path = killed_dir / "results.csv"
    with csv_path.open(newline="") as fh:
        all_rows = list(csv.DictReader(fh))
    if len(all_rows) >= 2:
        keep = all_rows[:-1]  # drop the "in-flight" row
        with csv_path.open() as fh:
            fieldnames = fh.readline().strip().split(",")
        with csv_path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            for row in keep:
                writer.writerow(row)

    resumed_summary = run_grid(
        grid, killed_dir, root_seed=200, epsilons=(0.1,), n_bias_draws=5, resume=True
    )

    with csv_path.open(newline="") as fh:
        resumed_rows = list(csv.DictReader(fh))

    # No duplicate instance_ids.
    ids = [row["instance_id"] for row in resumed_rows]
    assert len(ids) == len(set(ids))

    # Same final accept/reject tally as the uninterrupted run.
    assert resumed_summary["n_accepted"] == full_summary["n_accepted"]
    assert resumed_summary["rejection_counts"] == full_summary["rejection_counts"]
    assert len(resumed_rows) == full_summary["n_accepted"]


def test_run_grid_without_resume_refuses_to_clobber(tmp_path: Path):
    grid = _small_grid()[:2]
    run_grid(grid, tmp_path, root_seed=300, epsilons=(0.1,), n_bias_draws=5)
    with pytest.raises(FileExistsError):
        run_grid(grid, tmp_path, root_seed=300, epsilons=(0.1,), n_bias_draws=5, resume=False)


# --- hash-seed determinism (mandatory invariant test) ------------------------


def _run_under_hashseeds(code: str) -> None:
    import subprocess
    import sys

    outs = set()
    for hs in ("0", "1", "12345"):
        env = {"PYTHONPATH": "src", "PYTHONHASHSEED": hs, "PATH": "/usr/bin:/bin"}
        res = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, env=env, check=True
        )
        outs.add(res.stdout.strip())
    assert len(outs) == 1, f"non-deterministic across PYTHONHASHSEED: {outs}"


def test_run_instance_is_hash_invariant():
    code = (
        "from bkrobust.synth.runner import run_instance, GateRejectedError;"
        "out = None;"
        "\nfor seed in range(50):\n"
        "    try:\n"
        "        inst = run_instance(\n"
        "            instance_id='h', generator_name='decoupled_backdoor',\n"
        "            generator_params={'coupling': 0.0}, n=6, seed=seed,\n"
        "            knows_fraction=0.8, corruption_name='flip', corruption_rate=0.3,\n"
        "            epsilons=(0.05, 0.2), n_bias_draws=5)\n"
        "    except GateRejectedError:\n"
        "        continue\n"
        "    row = inst.to_row()\n"
        "    out = {k: v for k, v in row.items() if not k.startswith('time_')}\n"
        "    break\n"
        "print(repr(sorted(out.items())))"
    )
    _run_under_hashseeds(code)


def test_gate_is_hash_invariant():
    code = (
        "import numpy as np;"
        "from bkrobust.synth.generators import erdos_renyi_dag;"
        "from bkrobust.demo.example import dag_to_cpdag;"
        "from bkrobust.synth.runner import gate;"
        "out = [];"
        "\nfor seed in range(10):\n"
        "    rng = np.random.default_rng(seed)\n"
        "    dag = erdos_renyi_dag(7, rng, 0.4)\n"
        "    cpdag = dag_to_cpdag(dag)\n"
        "    out.append(gate(dag, cpdag, sorted(dag.nodes)[0], sorted(dag.nodes)[-1]))\n"
        "print(repr(out))"
    )
    _run_under_hashseeds(code)
