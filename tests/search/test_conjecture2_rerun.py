"""Tests for the Conjecture 2 re-run on the CORRECTED space.

These pin down the things that would silently make the re-run untrustworthy:
that it reduces to the old sweep exactly where the two spaces coincide (n=3),
that it actually exercises the extra states the correction adds (n=4), that a
saved combination is reproducible from its recorded fields alone, that output
is byte-identical across PYTHONHASHSEED, and that a crash mid-sweep loses at
most the one CPDAG in flight.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path

from bkrobust.core.oracle import is_valid
from bkrobust.core.spacelib import distances_from
from bkrobust.demo.evaluate import all_valid_adjustment_sets_mpdag
from bkrobust.demo.space import enumerate_space
from bkrobust.search.conjecture2_rerun import (
    CorrectedStudyTotals,
    build_old_vs_new,
    identify_density_gap,
    reproduce_conjecture2,
    run_dense_closure,
    run_study_corrected,
    write_corrected_summary,
)
from bkrobust.search.conjecture_study import all_cpdags
from bkrobust.search.conjectures import check_conjecture2
from bkrobust.search.space_fixed import build_space_correct, enumerate_space_correct

# --- deliverable 1: the main sweep reduces to the old one where spaces agree


def test_n3_reproduces_old_zero_counterexamples(tmp_path):
    """n=3 has no k>=5 CPDAG, so the defect never fires: old and new spaces.

    are identical there, and the corrected sweep must reproduce the old
    sweep's totals and its zero counterexamples exactly.
    """
    totals, ces = run_study_corrected(3, tmp_path, verbose=False)
    assert ces == []
    assert totals.c1_violations == 0
    assert totals.c2_violations == 0
    # Matches results/search/conjectures/n3.json's totals exactly.
    assert totals.n_dags == 25
    assert totals.n_cpdags == 11
    assert totals.n_spaces == 7
    assert totals.n_combos == 150
    assert totals.n_radius_comparisons == 150
    assert totals.skipped_cpdag_no_undirected == 4
    # The harness is what's under test here, not the space: with identical
    # spaces there is nothing for the correction to add.
    assert totals.n_states_added_vs_old == 0
    assert totals.n_combos_new == 0


def test_n3_every_space_matches_the_old_enumeration_exactly():
    """Independent of the sweep harness: at n=3 old and new spaces coincide."""
    for cpdag in all_cpdags(3):
        if not cpdag.undirected_edges:
            continue
        old = {g.edge_string() for g in enumerate_space(cpdag)}
        new = {g.edge_string() for g in enumerate_space_correct(cpdag)}
        assert old == new


# --- the correction is actually exercised at n=4 ----------------------------


def test_n4_corrected_sweep_examines_strictly_more_space_elements():
    """A k=5 CPDAG at n=4 -- where the Task 0 defect starts to bite -- must.

    yield a strictly bigger corrected space than the old enumerator gives it.
    """
    affected = [c for c in all_cpdags(4) if len(c.undirected_edges) == 5]
    assert affected, "n=4 must contain at least one k=5 CPDAG to test against"
    found_strict_increase = False
    for cpdag in affected:
        old_n = len(enumerate_space(cpdag))
        new_n = len(enumerate_space_correct(cpdag))
        assert new_n >= old_n, "the correction must never REMOVE a state"
        if new_n > old_n:
            found_strict_increase = True
    assert found_strict_increase, "expected the defect to bite at k=5, n=4"


def test_n4_sweep_totals_report_added_states(tmp_path):
    """The full n=4 sweep's n_states_added_vs_old is positive, and its skip.

    counts match the CPDAG structure independently known for n=4.
    """
    totals, _ces = run_study_corrected(4, tmp_path, verbose=False)
    assert totals.n_dags == 543
    assert totals.n_cpdags == 185
    assert totals.n_states_added_vs_old > 0
    # from n4.json: 59 CPDAGs have no undirected edges, none exceed k=6.
    assert totals.skipped_cpdag_no_undirected == 59
    assert totals.skipped_cpdag_too_many_undirected == 0


# --- reproducibility of a saved combination ---------------------------------


def _find_combo_with_valid_z() -> tuple | None:
    """Scan for any (cpdag, g0, x, y, z) with a valid Z (most have none)."""
    for cpdag in all_cpdags(4):
        if len(cpdag.undirected_edges) < 2:
            continue
        space = build_space_correct(cpdag)
        nodes = list(cpdag.nodes)
        for g0 in space.elements:
            for x, y in ((nodes[0], nodes[1]), (nodes[1], nodes[0])):
                zs = all_valid_adjustment_sets_mpdag(g0, x, y)
                valid_zs = [zz for zz in zs if is_valid(zz, g0, x, y)]
                if valid_zs:
                    return cpdag, space, g0, x, y, valid_zs[0]
    return None


def test_saved_combination_is_reproducible_from_its_fields_alone():
    """Whatever check_conjecture2 records for a real combination must come.

    back byte-for-byte from ONLY (cpdag, g0, x, y, z) -- the property that
    lets a saved counterexample be trusted from its JSON record alone.
    """
    found = _find_combo_with_valid_z()
    assert found is not None, "need at least one valid Z to exercise this path"
    cpdag, space, g0, x, y, z = found

    dists = distances_from(space, g0)

    def fails(g: object, _z: frozenset = z, _x: str = x, _y: str = y) -> bool:
        return not is_valid(_z, g, _x, _y)

    original = check_conjecture2(space, g0, fails, dists)
    rec = {"cpdag": cpdag, "g0": g0, "x": x, "y": y, "z": sorted(z)}
    reproduced = reproduce_conjecture2(rec)

    assert reproduced.r_full == original.r_full
    assert reproduced.r_up == original.r_up
    assert reproduced.holds == original.holds
    assert reproduced.witness_full == original.witness_full
    assert reproduced.witness_up == original.witness_up


def test_any_recorded_counterexample_reproduces(tmp_path):
    """Any counterexample the n=3 sweep records (expected: none) must.

    reproduce independently from its saved fields; a no-op guard when clean.
    """
    totals, ces = run_study_corrected(3, tmp_path, verbose=False)
    assert isinstance(totals, CorrectedStudyTotals)
    by_string = {c.edge_string(): c for c in all_cpdags(3)}
    for ce in ces:
        if ce["conjecture"] != 2:
            continue
        cpdag = by_string[ce["cpdag"]]
        space = build_space_correct(cpdag)
        by_g = {g.edge_string(): g for g in space.elements}
        rec = {
            "cpdag": cpdag,
            "g0": by_g[ce["g0"]],
            "x": ce["x"],
            "y": ce["y"],
            "z": ce["z"],
        }
        reproduced = reproduce_conjecture2(rec)
        assert reproduced.r_full == ce["r_full"]
        assert reproduced.r_up == ce["r_up"]
        assert reproduced.witness_full == ce["witness_full"]
        assert reproduced.witness_up == ce["witness_up"]


# --- PYTHONHASHSEED invariance (mandatory) ----------------------------------


def test_output_is_identical_across_pythonhashseed(tmp_path):
    """No global RNG and no set/dict iteration order leaking into output:.

    the sweep is byte-identical across PYTHONHASHSEED, checked via real
    subprocesses since a live process's own hash seed can't be changed.
    """
    script = tmp_path / "run_one.py"
    script.write_text(
        "import json, sys\n"
        "from dataclasses import asdict\n"
        "from bkrobust.search.conjecture2_rerun import run_study_corrected\n"
        "totals, ces = run_study_corrected(3, sys.argv[1], verbose=False)\n"
        "print(json.dumps({'totals': asdict(totals), 'ces': ces}, sort_keys=True))\n"
    )
    outputs = []
    for seed in (0, 1, 12345):
        out_dir = tmp_path / f"seed_{seed}"
        full_env = dict(os.environ)
        full_env["PYTHONHASHSEED"] = str(seed)
        full_env["PYTHONPATH"] = "src"
        result = subprocess.run(
            [sys.executable, str(script), str(out_dir)],
            capture_output=True,
            text=True,
            env=full_env,
            cwd=str(Path(__file__).resolve().parents[2]),
            timeout=120,
        )
        assert result.returncode == 0, result.stderr
        outputs.append(result.stdout.strip())
    assert outputs[0] == outputs[1] == outputs[2], "output diverges across PYTHONHASHSEED"


# --- resume: a crash mid-sweep loses at most one CPDAG ----------------------


def test_resume_after_simulated_crash_matches_full_run(tmp_path):
    """Truncating the checkpoint to all-but-the-last CPDAG and resuming must.

    reach exactly the full run's totals, without rewriting untouched lines.
    """
    full_dir = tmp_path / "full"
    full_totals, full_ces = run_study_corrected(3, full_dir, verbose=False)

    full_jsonl = full_dir / "n3_corrected.jsonl"
    lines = full_jsonl.read_text().splitlines()
    assert len(lines) == 11, "expected one checkpoint line per CPDAG at n=3"

    resumed_dir = tmp_path / "resumed"
    resumed_dir.mkdir()
    truncated = resumed_dir / "n3_corrected.jsonl"
    truncated.write_text("\n".join(lines[:-1]) + "\n")

    resumed_totals, resumed_ces = run_study_corrected(3, resumed_dir, verbose=False)

    assert asdict(resumed_totals) == asdict(full_totals)
    assert resumed_ces == full_ces

    resumed_lines = truncated.read_text().splitlines()
    assert len(resumed_lines) == 11
    # Every line but the last is byte-identical to the un-truncated run --
    # i.e. resume did not re-run or rewrite anything it already had.
    assert resumed_lines[:-1] == lines[:-1]

    summary_path = resumed_dir / "n3_corrected.json"
    assert summary_path.exists(), "a fully-covered checkpoint must auto-write the summary"


def test_resume_is_a_noop_on_an_already_complete_checkpoint(tmp_path):
    """A second call on an already-complete checkpoint rewrites nothing."""
    totals1, ces1 = run_study_corrected(3, tmp_path, verbose=False)
    jsonl = tmp_path / "n3_corrected.jsonl"
    before = jsonl.read_text()
    totals2, ces2 = run_study_corrected(3, tmp_path, verbose=False)
    after = jsonl.read_text()
    assert before == after
    assert asdict(totals1) == asdict(totals2)
    assert ces1 == ces2


# --- density gap machinery (deliverable 2), exercised at small scale -------


def test_identify_density_gap_and_dense_closure_on_n4(tmp_path):
    """n=4 has exactly one CPDAG (the fully-undirected K4 skeleton) whose.

    corrected space exceeds max_space=200; the dense closure must find it,
    examine it fully, and its summary must account for it.
    """
    run_study_corrected(4, tmp_path, max_space=200, verbose=False)
    gap = identify_density_gap(4, tmp_path)
    assert len(gap) == 1
    assert gap[0]["k_undirected"] == 6

    totals, _ces, infeasible = run_dense_closure(4, gap, tmp_path, verbose=False)
    assert infeasible == []
    assert totals.n_cpdags == 1
    assert totals.n_spaces == 1
    dense_summary = json.loads((tmp_path / "n4_dense.json").read_text())
    assert dense_summary["n_cpdags_completed"] == 1
    assert dense_summary["n_cpdags_infeasible"] == 0


def test_dense_closure_records_infeasible_cpdags_rather_than_dropping_them(tmp_path):
    """A CPDAG over max_space_dense is named, with its size and the reason,.

    in the return value and the summary file -- never silently skipped.
    """
    run_study_corrected(4, tmp_path, max_space=200, verbose=False)
    gap = identify_density_gap(4, tmp_path)
    assert gap

    _totals, _ces, infeasible = run_dense_closure(
        4, gap, tmp_path, max_space_dense=100, verbose=False
    )
    assert len(infeasible) == 1
    assert infeasible[0]["cpdag"] == gap[0]["cpdag"]
    assert infeasible[0]["reason"]
    summary = json.loads((tmp_path / "n4_dense.json").read_text())
    assert summary["n_cpdags_infeasible"] == 1
    assert summary["infeasible"][0]["cpdag"] == gap[0]["cpdag"]


# --- deliverable 3 machinery -------------------------------------------------


def test_build_old_vs_new_reads_old_verbatim_and_never_recomputes(tmp_path):
    """The comparison reads the OLD n3.json as-is and never recomputes it."""
    run_study_corrected(3, tmp_path, verbose=False)
    comparison = build_old_vs_new((3,), new_dir=tmp_path)
    entry = comparison["per_n"]["3"]
    assert entry["old"]["totals"]["n_combos"] == 150  # from the committed n3.json
    assert entry["new"]["totals"]["n_combos"] == 150
    assert "verdict" in comparison


def test_write_corrected_summary_is_idempotent_given_the_same_checkpoint(tmp_path):
    """Re-summarizing an unchanged checkpoint produces byte-identical output."""
    run_study_corrected(3, tmp_path, verbose=False)
    p1 = write_corrected_summary(3, tmp_path)
    text1 = p1.read_text()
    p2 = write_corrected_summary(3, tmp_path)
    text2 = p2.read_text()
    assert text1 == text2
