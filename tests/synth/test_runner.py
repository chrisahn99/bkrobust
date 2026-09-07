"""Tests for bkrobust.synth.runner.

What each test discriminates:

* ``gate`` tests use small hand-built graphs (verified by hand in the
  docstrings/comments here) to pin down each reachable rejection reason,
  plus one accepted case -- rather than trusting the gate only from its
  behaviour on random generator output.
* ``run_instance`` tests check a non-degenerate accepted instance really has
  ``r_val >= 1`` (per the radius convention: Z is read off G0, so it is
  valid there by construction) -- both on a single case and, after a review
  found this silently violated for 9% of a real run, over a broad sample of
  many seeds and generators. A genuinely degenerate hand-built instance must
  raise ``GateRejectedError`` rather than silently returning something.
* Two regression tests reproduce real bugs found by review of an actual
  pilot run rather than invented cases: ``z_invalid_at_g0`` (a captured seed
  where the optimal-set formula returned an invalid empty set) and
  ``g0_not_in_space`` (a captured structural mismatch between
  ``apply_orientations`` and ``enumerate_space``, reconstructed as a minimal
  4-node example since the exact case as pasted was internally
  contradictory -- see that test's docstring).
* ``run_grid`` tests check resume-safety directly (kill-and-restart must not
  duplicate rows) and resilience: an unexpected exception on one grid point
  must be recorded, not crash the sweep.
* The mandatory PYTHONHASHSEED-invariance test covers ``run_instance``,
  since it is the function whose RNG-consumption order is most at risk (it
  threads one generator through graph generation, knowledge simulation and
  bias sampling in sequence).

``decoupled_backdoor`` is deliberately absent from the "accepted instance"
tests below: after the fix documented in ``generators.py`` (a dedicated node
forcing X->Y to be compelled, which was needed for CPDAG-level
identification), every instance from that generator is gate-rejected with
``"no_atomic_perturbation_changes_validity"`` -- see
``tests/synth/test_generators.py`` for that finding in full and
``test_decoupled_backdoor_is_always_gate_rejected`` below for the
runner-level confirmation. ``erdos_renyi`` is used instead wherever an
accepted, non-degenerate instance is needed.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

import bkrobust.synth.runner as runner
from bkrobust.core.conventions import UNREACHED
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions, is_valid_mpdag
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


def test_gate_rejects_no_causal_path():
    """Y -> X: X cannot possibly cause Y, so the true total effect is zero.

    This is the structural check added after a review found run_instance
    silently accepting exactly this shape (via optimal_adjustment_set_mpdag
    returning an empty, invalid set) -- see test_run_instance_z_invalid_at_g0
    below for the run_instance-level regression this generalises from.
    """
    dag = MPDAG(["X", "Y", "Z"], directed=[("Y", "X"), ("Z", "Y"), ("Z", "X")])
    cpdag = dag_to_cpdag(dag)
    assert gate(dag, cpdag, "X", "Y") == (False, "no_causal_path")


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
    """X -> Y with confounder Z (itself ambiguous via its own parent W): the textbook shape."""
    dag = MPDAG(
        ["X", "Y", "Z", "W"],
        directed=[("W", "Z"), ("Z", "X"), ("Z", "Y"), ("X", "Y")],
    )
    cpdag = dag_to_cpdag(dag)
    assert gate(dag, cpdag, "X", "Y") == (True, "ok")


# --- run_instance: healthy generator (erdos_renyi) ---------------------


def test_run_instance_r_val_at_least_one_on_nondegenerate_case():
    """The radius convention: Z is valid at G0 by construction, so r_val >= 1."""
    inst = run_instance(
        instance_id="probe",
        generator_name="erdos_renyi",
        generator_params={"edge_prob": 0.5},
        n=7,
        seed=3,
        knows_fraction=0.8,
        corruption_name="none",
        corruption_rate=0.0,
    )
    assert inst.r_val >= 1
    assert inst.r_val != UNREACHED
    assert inst.method == "bfs_exact"
    assert set(inst.timings) >= {"graph", "gate", "knowledge", "g0_and_z", "space"}


def test_run_instance_r_opt_and_r_eps_are_populated():
    inst = run_instance(
        instance_id="probe2",
        generator_name="erdos_renyi",
        generator_params={"edge_prob": 0.5},
        n=7,
        seed=3,
        knows_fraction=0.8,
        corruption_name="none",
        corruption_rate=0.0,
        epsilons=(0.05, 0.2),
    )
    assert inst.r_opt >= 1
    assert set(inst.r_eps) == {"0.05", "0.2"}
    for r in inst.r_eps.values():
        assert r == UNREACHED or r >= 1


def test_run_instance_raises_gate_rejected_on_degenerate_generator_params():
    """erdos_renyi(n=8, edge_prob=0.5, seed=22) reliably exceeds MAX_UNDIRECTED_EDGES."""
    with pytest.raises(GateRejectedError) as excinfo:
        run_instance(
            instance_id="degenerate",
            generator_name="erdos_renyi",
            generator_params={"edge_prob": 0.5},
            n=8,
            seed=22,
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
    inst = run_instance(
        instance_id="tiered_probe",
        generator_name="erdos_renyi",
        generator_params={"edge_prob": 0.5},
        n=7,
        seed=5,
        knows_fraction=0.8,
        corruption_name="tiered",
        corruption_rate=0.0,
        tier_params={"n_tiers": 3, "corruption_rate": 0.0},
    )
    assert inst.corruption == "tiered"
    assert inst.r_val >= 1


def _non_timing_row(instance) -> dict:
    return {k: v for k, v in instance.to_row().items() if not k.startswith("time_")}


def test_run_instance_reproducible_from_seed():
    kwargs = dict(
        instance_id="repro",
        generator_name="erdos_renyi",
        generator_params={"edge_prob": 0.5},
        n=7,
        seed=3,
        knows_fraction=0.8,
        corruption_name="flip",
        corruption_rate=0.3,
    )
    a = run_instance(**kwargs)
    b = run_instance(**kwargs)
    # Wall-clock timings legitimately differ between the two calls; every
    # other field must match exactly.
    assert _non_timing_row(a) == _non_timing_row(b)


def test_max_undirected_edges_is_a_positive_int():
    assert isinstance(MAX_UNDIRECTED_EDGES, int)
    assert MAX_UNDIRECTED_EDGES > 0


# --- run_instance: decoupled_backdoor is always gate-rejected --------------


def test_decoupled_backdoor_is_always_gate_rejected():
    """Runner-level confirmation of the finding documented in test_generators.py.

    Forcing X->Y (and, as a consequence, C1_last/C2_last's edges into X and
    Y) to be compelled -- required for CPDAG-level identification, see
    generators.py's "Revision history" -- makes the optimal adjustment set
    structurally invariant across the whole perturbation space. Every
    instance is therefore rejected by the sanity gate. Swept over several
    seeds and both coupling endpoints to make sure this is not a one-seed
    fluke.
    """
    reasons = set()
    for seed in range(10):
        for coupling in (0.0, 1.0):
            with pytest.raises(GateRejectedError) as excinfo:
                run_instance(
                    instance_id=f"decoupled_{seed}_{coupling}",
                    generator_name="decoupled_backdoor",
                    generator_params={"coupling": coupling},
                    n=7,
                    seed=seed,
                    knows_fraction=0.8,
                    corruption_name="flip",
                    corruption_rate=0.3,
                )
            reasons.add(excinfo.value.reason)
    assert reasons == {"no_atomic_perturbation_changes_validity"}


# --- regression: z_invalid_at_g0 (real bug, captured from a pilot run) ------


def test_run_instance_z_invalid_at_g0_regression():
    """Captured case: erdos_renyi(n=7, edge_prob=0.5, seed=143) used to return r_val=0.

    A review of an actual pilot run found 9% of accepted rows with
    ``r_val == 0``, which conventions.py declares impossible (Z is read off
    G0 and must be valid there). The cause: optimal_adjustment_set_mpdag
    computes O from ``cn(x, y) = descendants(x) & (ancestors(y) | {y})`` and
    returns the *empty* set whenever ``cn`` is empty, regardless of whether
    that is actually a valid (backdoor-blocking) choice. This seed is one
    concrete case where it was not. It must now be gate-rejected, with
    ``z_invalid_at_g0`` as the backstop reason (``no_causal_path`` is the
    primary defence, at selection time, but this specific case has a causal
    path and would only be caught by the direct validity check).
    """
    with pytest.raises(GateRejectedError) as excinfo:
        run_instance(
            instance_id="z_invalid_regression",
            generator_name="erdos_renyi",
            generator_params={"edge_prob": 0.5},
            n=7,
            seed=143,
            knows_fraction=0.7,
            corruption_name="flip",
            corruption_rate=0.3,
        )
    assert excinfo.value.reason == "z_invalid_at_g0"


def test_run_instance_r_val_is_always_at_least_one_over_many_seeds():
    """The invariant that was silently violated: check it broadly, not just once.

    Every accepted instance, across many seeds and two different generators,
    must have r_val >= 1 (never 0, and UNREACHED is a distinct, allowed
    sentinel that this loop also tolerates as neither >= 1 nor a violation).
    """
    n_accepted = 0
    for generator_name, params, n in (
        ("erdos_renyi", {"edge_prob": 0.5}, 7),
        ("scale_free", {"m_attach": 2}, 7),
        ("block", {"n_blocks": 2, "p_within": 0.6, "p_between": 0.1}, 7),
    ):
        for seed in range(80):
            try:
                inst = run_instance(
                    instance_id=f"{generator_name}_{seed}",
                    generator_name=generator_name,
                    generator_params=params,
                    n=n,
                    seed=seed,
                    knows_fraction=0.7,
                    corruption_name="flip",
                    corruption_rate=0.3,
                )
            except GateRejectedError:
                continue
            n_accepted += 1
            assert inst.r_val == UNREACHED or inst.r_val >= 1, (
                f"convention violated: {generator_name} seed={seed} r_val={inst.r_val}"
            )
    assert n_accepted >= 20, "expected a healthy number of accepted instances to check"


# --- regression: g0_not_in_space (real bug, reconstructed minimally) --------


def _k4_minus_one_edge_case() -> tuple[MPDAG, MPDAG]:
    """A minimal, independently-verified reconstruction of the captured g0_not_in_space bug.

    The bug report's exact 7-node example, as pasted, contained both
    "V1->V6" and "V6->V1" -- an internal contradiction (almost certainly a
    transcription slip) that MPDAG's constructor rejects outright. Rather
    than guess which one was meant, this reconstructs the same underlying
    mechanism from scratch and verifies it independently: a CPDAG whose
    undirected component is K4 minus one edge (a 4-cycle A-B-C-D-A plus the
    chord A-C) is chordal, hence a valid CPDAG. Asserting *only* the chord's
    orientation (A->C) leaves the 4-cycle A-B-C-D-A undirected -- Meek's
    rules do not fire on it (no adjacent pair becomes non-adjacent) -- so the
    resulting g0 has an undirected component that is a chordless 4-cycle
    when only its undirected edges are inspected, even though the chord
    still exists as a directed edge. is_chordal_components looks only at
    undirected edges, so is_valid_mpdag(g0) is False and enumerate_space
    omits it -- despite g0 having 5 consistent DAG extensions of its own
    (verified below), i.e. it is a perfectly meaningful graph, just one the
    space's own filter excludes.
    """
    nodes = ["A", "B", "C", "D"]
    undirected = [("A", "B"), ("B", "C"), ("C", "D"), ("A", "D"), ("A", "C")]
    cpdag = MPDAG(nodes, directed=[], undirected=undirected)
    assert is_valid_mpdag(cpdag)
    g0 = apply_orientations(cpdag, [("A", "C")])
    assert g0 is not None
    assert not is_valid_mpdag(g0), "the captured mechanism needs g0 itself to fail is_valid_mpdag"
    return cpdag, g0


def test_g0_not_in_space_mechanism_is_reproduced():
    """The bare mechanism: g0 has consistent extensions but is excluded from the space anyway."""
    from bkrobust.core.spacelib import build_space

    cpdag, g0 = _k4_minus_one_edge_case()
    assert len(enumerate_dag_extensions(g0)) > 0, "g0 must be meaningful, not an empty shell"
    space = build_space(cpdag)
    assert g0 not in space.neighbours


def test_run_instance_gates_on_g0_not_in_space(monkeypatch: pytest.MonkeyPatch):
    """run_instance must GateRejectedError("g0_not_in_space"), not raise the underlying KeyError.

    gate(), _pick_treatment_outcome, draw_k_true, optimal_adjustment_set_mpdag
    and is_valid are all monkeypatched to bypass the ordinary pipeline for
    this one call -- the 4-node K4-minus-edge graph above has no natural
    treatment/outcome/adjustment-set story of its own (in particular, A's
    optimal set is not identified across g0's own 5 extensions, since B and D
    play different structural roles in different ones -- a fact about this
    toy graph, not about the mechanism under test). The point of this test is
    narrowly the crash-vs-gate behaviour once G0 has been computed and found
    to lack extensions in the space, which is a separate concern from
    identification.
    """
    cpdag, _g0 = _k4_minus_one_edge_case()
    dag = enumerate_dag_extensions(cpdag)[0]

    monkeypatch.setitem(runner.GENERATORS, "erdos_renyi", lambda n, rng, **kw: dag)
    monkeypatch.setattr(runner, "gate", lambda *a, **kw: (True, "ok"))
    monkeypatch.setattr(runner, "_pick_treatment_outcome", lambda *a, **kw: ("A", "C"))
    monkeypatch.setattr(runner, "draw_k_true", lambda *a, **kw: [("A", "C")])
    monkeypatch.setattr(runner, "optimal_adjustment_set_mpdag", lambda *a, **kw: set())
    monkeypatch.setattr(runner, "is_valid", lambda *a, **kw: True)

    with pytest.raises(GateRejectedError) as excinfo:
        run_instance(
            instance_id="g0_not_in_space_case",
            generator_name="erdos_renyi",
            generator_params={},
            n=4,
            seed=0,
            knows_fraction=1.0,
            corruption_name="none",
            corruption_rate=0.0,
        )
    assert excinfo.value.reason == "g0_not_in_space"


# --- run_grid: resume safety and error resilience ---------------------


def _small_grid() -> list[dict]:
    """A small mixed grid: some points accept, some are gate-rejected -- both matter for resume.

    15 erdos_renyi points at a fixed, verified-non-degenerate parameterisation
    give a handful of accepted rows regardless of which root_seed the tests
    below use (checked directly for root_seed in {100, 200, 300}: 2-4 of the
    15 accept in each case) -- acceptance depends on the drawn graph, not on
    knows_fraction/corruption_rate, so those are held fixed here rather than
    swept, to keep the accept/reject mix a property of the seed alone.
    """
    grid = []
    for _seed_offset in range(15):
        grid.append(
            {
                "generator_name": "erdos_renyi",
                "generator_params": {"edge_prob": 0.5},
                "n": 7,
                "knows_fraction": 0.7,
                "corruption_name": "flip",
                "corruption_rate": 0.2,
            }
        )
    # A guaranteed-degenerate point (decoupled_backdoor is always rejected --
    # see test_decoupled_backdoor_is_always_gate_rejected), to exercise
    # rejection bookkeeping.
    grid.append(
        {
            "generator_name": "decoupled_backdoor",
            "generator_params": {"coupling": 1.0},
            "n": 7,
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
    assert summary["errors"] == []

    csv_path = tmp_path / "results.csv"
    assert csv_path.exists()
    with csv_path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == summary["n_accepted"]
    for row in rows:
        assert int(row["r_val"]) == UNREACHED or int(row["r_val"]) >= 1

    manifest_path = tmp_path / "manifest.json"
    assert manifest_path.exists()
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


def test_run_grid_records_run_error_and_continues(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """An unexpected exception on one grid point must not kill the sweep.

    run_instance is monkeypatched to raise a plain RuntimeError (something
    that is deliberately *not* a GateRejectedError) for exactly one grid
    point; run_grid must record it under rejection_counts["run_error"] with
    the exception text in errors, and still process every other point
    normally.
    """
    grid = _small_grid()
    real_run_instance = runner.run_instance
    boom_index = 2

    def flaky(*args: object, **kwargs: object) -> object:
        if kwargs.get("instance_id") == runner._instance_id(boom_index):
            raise RuntimeError("synthetic failure for test coverage")
        return real_run_instance(*args, **kwargs)

    monkeypatch.setattr(runner, "run_instance", flaky)
    summary = run_grid(grid, tmp_path, root_seed=100, epsilons=(0.1,), n_bias_draws=5)

    assert summary["rejection_counts"].get("run_error") == 1
    assert len(summary["errors"]) == 1
    assert summary["errors"][0]["instance_id"] == runner._instance_id(boom_index)
    assert summary["errors"][0]["exception_type"] == "RuntimeError"
    assert "synthetic failure" in summary["errors"][0]["message"]
    # Every other point was still attempted: total accounts for all of them.
    assert summary["n_accepted"] + sum(summary["rejection_counts"].values()) == len(grid)


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
        "from bkrobust.synth.runner import run_instance;"
        "inst = run_instance(\n"
        "    instance_id='h', generator_name='erdos_renyi',\n"
        "    generator_params={'edge_prob': 0.5}, n=7, seed=3,\n"
        "    knows_fraction=0.8, corruption_name='flip', corruption_rate=0.3,\n"
        "    epsilons=(0.05, 0.2), n_bias_draws=5)\n"
        "row = inst.to_row()\n"
        "out = {k: v for k, v in row.items() if not k.startswith('time_')}\n"
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
