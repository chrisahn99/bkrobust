"""Regression tests for the declarative encodings.

Every one of these pins something that was once wrong, or that would silently
produce plausible-but-false numbers if it broke. In particular the two encoding
bugs found during development are pinned directly:

* the closure encoding must forbid Meek *firing configurations*, not merely
  imply the consequent -- dropping R3/R4's undirectedness premises made the
  encoding too strong and lost 67 of 133 CPDAGs;
* a collider on a d-connecting path needs **both** arrows pointing in, not
  either -- the weaker clause let a chain through ``Z`` pass as a collider and
  reported failures that were not failures, making radii too small.
"""

from __future__ import annotations

import itertools

import pytest

from bkrobust.core.conventions import UNREACHED
from bkrobust.core.oracle import is_valid
from bkrobust.core.spacelib import distances_from, radius
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.scenario import OUTCOME, TREATMENT, scenarios, true_dag
from bkrobust.sat.closure import solutions
from bkrobust.sat.e1 import radius_e1
from bkrobust.sat.e2 import radius_e2
from bkrobust.sat.e3 import radius_e3
from bkrobust.sat.verify import verify_walk
from bkrobust.search.conjecture_study import all_cpdags
from bkrobust.search.exact import radius_local_up
from bkrobust.search.space_fixed import build_space_correct


def _instances(n: int, limit: int) -> list[tuple[MPDAG, MPDAG, str, str, frozenset[str]]]:
    """Yield ``(cpdag, g0, x, y, z)`` with a non-degenerate optimal adjustment set."""
    out = []
    for cpdag in all_cpdags(n):
        if not cpdag.undirected_edges:
            continue
        for g0 in build_space_correct(cpdag).elements:
            for x, y in itertools.permutations(sorted(cpdag.nodes), 2):
                o = optimal_adjustment_set_mpdag(g0, x, y)
                if o is None:
                    continue
                z = frozenset(o)
                if not z or not is_valid(z, g0, x, y):
                    continue
                out.append((cpdag, g0, x, y, z))
                if len(out) >= limit:
                    return out
    return out


@pytest.mark.parametrize("n", [3, 4])
def test_closure_encoding_is_the_corrected_space(n: int) -> None:
    """Stage 1: the CP-SAT solution set is exactly the enumerated corrected space."""
    from bkrobust.search.space_fixed import enumerate_space_correct

    for cpdag in list(all_cpdags(n))[:25]:
        got = {frozenset(s) for s in solutions(cpdag)}
        want = {frozenset(g.directed_edges) for g in enumerate_space_correct(cpdag)}
        assert got == want, cpdag.edge_string()


def test_worked_example_radii() -> None:
    """All three encodings reproduce the published 3 / 3 / 2 on the worked example."""
    cpdag = dag_to_cpdag(true_dag())
    want = {"A": 3, "B": 3, "C": 2}
    for name, spec in scenarios().items():
        g0 = apply_orientations(cpdag, spec["knowledge"])
        z = frozenset(optimal_adjustment_set_mpdag(g0, TREATMENT, OUTCOME))
        assert radius_e1(cpdag, g0, TREATMENT, OUTCOME, z).radius == want[name]
        assert radius_e2(cpdag, g0, TREATMENT, OUTCOME, z).radius == want[name]
        assert radius_e3(cpdag, g0, TREATMENT, OUTCOME, z, max_k=6).radius == want[name]


def test_encodings_reproduce_bfs_and_local_up() -> None:
    """Stage 3: E1, E2 and E3 all agree with brute-force BFS and with ``local_up``."""
    for cpdag, g0, x, y, z in _instances(4, 60):
        space = build_space_correct(cpdag)
        dists = distances_from(space, g0)
        fails = lambda g, _z=z, _x=x, _y=y: not is_valid(_z, g, _x, _y)  # noqa: E731
        r_bfs, _ = radius(space, dists, fails)
        assert radius_local_up(cpdag, g0, fails).radius == r_bfs
        assert radius_e1(cpdag, g0, x, y, z).radius == r_bfs
        assert radius_e2(cpdag, g0, x, y, z).radius == r_bfs
        r3 = radius_e3(cpdag, g0, x, y, z, max_k=6)
        if r_bfs == UNREACHED:
            # E3 cannot prove unreachability; exhausting its budget is all it can do.
            assert r3.radius == UNREACHED and r3.exhausted_to == 6 and not r3.exact
        else:
            assert r3.radius == r_bfs


def test_every_e3_walk_is_a_certified_covering_walk() -> None:
    """E3's soundness claim is checked by replay, not taken from the solver."""
    for cpdag, g0, x, y, z in _instances(4, 40):
        r3 = radius_e3(cpdag, g0, x, y, z, max_k=6)
        if not r3.walk:
            continue
        cert = verify_walk(cpdag, g0, r3.walk, x, y, z)
        assert cert.ok, cert.problems
        assert cert.length == r3.radius


def test_e3_never_beats_e1_here_which_is_evidence_not_proof() -> None:
    """No Conjecture 2 counterexample in this scope.

    ``r_E3 < r_E1`` would refute Conjecture 2 unconditionally. Equality does not
    confirm it -- a shorter path through a multi-orientation cover is invisible
    to both encodings -- so this test pins the absence of a refutation, nothing
    more.
    """
    for cpdag, g0, x, y, z in _instances(4, 40):
        r1 = radius_e1(cpdag, g0, x, y, z).radius
        r3 = radius_e3(cpdag, g0, x, y, z, max_k=6).radius
        if r3 != UNREACHED:
            assert r3 >= r1


def test_e2_optimum_never_needs_a_down_leg() -> None:
    """The join-distance optimum is always a pure retraction, as Conjecture 2 predicts."""
    for cpdag, g0, x, y, z in _instances(4, 40):
        res = radius_e2(cpdag, g0, x, y, z)
        if res.radius != UNREACHED:
            assert res.down_steps == 0


def test_results_are_bit_identical_across_processes() -> None:
    """Same radii, witnesses and solver counters regardless of hash seed.

    The counters are included deliberately: matching radii alone would not catch
    a search that explored a different tree and happened to land on the same
    answer.
    """
    import json
    import os
    import subprocess
    import sys

    code = (
        "import json,itertools;"
        "from bkrobust.demo.example import dag_to_cpdag;"
        "from bkrobust.demo.scenario import true_dag,scenarios,TREATMENT,OUTCOME;"
        "from bkrobust.demo.meek import apply_orientations;"
        "from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag;"
        "from bkrobust.sat.e1 import radius_e1;"
        "from bkrobust.sat.e3 import radius_e3;"
        "c=dag_to_cpdag(true_dag());o=[]\n"
        "for n,s in scenarios().items():\n"
        "    g=apply_orientations(c,s['knowledge']);"
        "z=frozenset(optimal_adjustment_set_mpdag(g,TREATMENT,OUTCOME));"
        "a=radius_e1(c,g,TREATMENT,OUTCOME,z);b=radius_e3(c,g,TREATMENT,OUTCOME,z,max_k=5);"
        "o.append([n,a.radius,a.conflicts,a.branches,b.radius,b.conflicts,b.branches,"
        "[sorted(list(e) for e in w) for w in b.walk]])\n"
        "print(json.dumps(o))"
    )
    outs = []
    for seed in ("0", "1", "12345"):
        env = {**os.environ, "PYTHONHASHSEED": seed}
        p = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, env=env, check=True
        )
        outs.append(json.loads(p.stdout))
    assert outs[0] == outs[1] == outs[2]
