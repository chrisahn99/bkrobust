"""Tests for the missing back-door lower bound (Axis B, session 2).

Three things must be nailed down, in this order of importance:

1. Admissibility, exhaustively (an inadmissible lower bound silently produces
   wrong radii -- this is the whole point of ``conventions.py``'s warning).
2. ``UNREACHED`` never gets treated as 0 or infinity when combining modes.
3. The two hand-built discriminating cases: mode (B) binding with mode (A)
   structurally impossible, and vice versa -- each derived by hand in the
   docstrings below, not just asserted.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys

from bkrobust.core.conventions import UNREACHED
from bkrobust.demo.graph import MPDAG
from bkrobust.search.backdoor_bound import (
    backdoor_lower_bound,
    combined_lower_bound,
    enumerate_simple_paths,
    run_admissibility_study,
)
from bkrobust.search.exact import descendant_lower_bound

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# --------------------------------------------------------------------------------
# Hand-built discriminating cases
# --------------------------------------------------------------------------------


def _mbias_case() -> tuple[MPDAG, MPDAG]:
    """X - M - W - Y where M's edge to X is compelled ``M -> X``.

    ``M -> X`` compelled means X's only adjacency is an edge pointing away
    from X that can never be reversed, so **no** directed path X -> ... -> W
    can ever exist: mode (A) (``descendant_lower_bound``) must be UNREACHED.

    But that same ``M -> X`` orientation is exactly what a back-door path
    needs at its first edge (arrow into X), for free. With ``W`` the adjustment
    set member, the path is unblocked once W becomes a collider: M -> W <- Y.
    Currently M - W is oriented ``W -> M`` (wrong way, but retractable) and
    W - Y is oriented ``Y -> W`` (already the right way). So the minimum is 1
    retraction (flip W -> M to M -> W), and mode (B) is finite.
    """
    cpdag = MPDAG(
        ["X", "M", "W", "Y"],
        directed=[("M", "X")],
        undirected=[("M", "W"), ("W", "Y")],
    )
    g0 = MPDAG(
        ["X", "M", "W", "Y"],
        directed=[("M", "X"), ("W", "M"), ("Y", "W")],
        undirected=[],
    )
    return cpdag, g0


def test_mode_b_binds_mode_a_impossible():
    """The M-bias case: descendant mode is UNREACHED, back-door mode is finite."""
    cpdag, g0 = _mbias_case()
    z = frozenset({"W"})
    assert descendant_lower_bound(cpdag, g0, z, "X") == UNREACHED
    assert backdoor_lower_bound(cpdag, g0, z, "X", "Y") == 1
    assert combined_lower_bound(cpdag, g0, z, "X", "Y") == 1


def _chain_case() -> tuple[MPDAG, MPDAG]:
    """X - M - Y, M -> Y compelled, X - M oriented ``M -> X`` (retractable).

    ``M -> X`` is already the direction a back-door path needs at its first
    edge, for free -- but ``M`` is *also* the adjustment-set member, and the
    vertex immediately adjacent to ``X`` on any back-door path can never be a
    collider (its edge to ``X`` is forced to point away from it, into ``X``,
    which is a tail at ``M``, not an arrowhead). So the path can never be
    unblocked at ``M``: mode (B) is UNREACHED, structurally, not just because
    this particular graph lacks retraction budget.

    Mode (A) is finite: M is already one retraction away from being a forward
    descendant of X (flip ``M -> X`` to ``X -> M``).
    """
    cpdag = MPDAG(
        ["X", "M", "Y"],
        directed=[("M", "Y")],
        undirected=[("X", "M")],
    )
    g0 = MPDAG(
        ["X", "M", "Y"],
        directed=[("M", "X"), ("M", "Y")],
        undirected=[],
    )
    return cpdag, g0


def test_mode_a_binds_mode_b_impossible():
    """The chain case: back-door mode is UNREACHED, descendant mode is finite."""
    cpdag, g0 = _chain_case()
    z = frozenset({"M"})
    assert descendant_lower_bound(cpdag, g0, z, "X") == 1
    assert backdoor_lower_bound(cpdag, g0, z, "X", "Y") == UNREACHED
    assert combined_lower_bound(cpdag, g0, z, "X", "Y") == 1


def test_first_internal_vertex_can_never_be_a_collider():
    """Corollary of the chain case, checked directly on the path enumerator.

    No matter how the other edges are set, the vertex adjacent to X on a
    back-door path is always a non-collider (its edge to X is forced into X,
    which is a tail at that vertex). If that vertex is in Z, EVERY path
    through it must fail to unblock -- this is what makes the chain case's
    UNREACHED structural, not incidental.
    """
    cpdag, g0 = _chain_case()
    paths = enumerate_simple_paths(cpdag, "X", "Y")
    assert paths == [["X", "M", "Y"]]
    # And Z = {M} kills the only path, whatever the other edge did.
    assert backdoor_lower_bound(cpdag, g0, frozenset({"M"}), "X", "Y") == UNREACHED


# --------------------------------------------------------------------------------
# UNREACHED handling
# --------------------------------------------------------------------------------


def test_both_modes_unreached_combined_is_unreached_not_zero():
    """A fully disconnected X: neither mode can say anything, combined must not fabricate 0."""
    cpdag = MPDAG(["X", "Y", "W"], directed=[], undirected=[("Y", "W")])
    g0 = MPDAG(["X", "Y", "W"], directed=[("W", "Y")], undirected=[])
    z = frozenset({"W"})
    assert descendant_lower_bound(cpdag, g0, z, "X") == UNREACHED
    assert backdoor_lower_bound(cpdag, g0, z, "X", "Y") == UNREACHED
    result = combined_lower_bound(cpdag, g0, z, "X", "Y")
    assert result == UNREACHED
    assert result != 0
    assert result >= UNREACHED  # not something more negative either


def test_combined_never_smaller_than_the_smaller_defined_mode():
    """When exactly one mode is defined, combined must equal it exactly -- never 0, never negative."""
    cpdag, g0 = _mbias_case()
    z = frozenset({"W"})
    a = descendant_lower_bound(cpdag, g0, z, "X")
    b = backdoor_lower_bound(cpdag, g0, z, "X", "Y")
    assert a == UNREACHED and b != UNREACHED
    c = combined_lower_bound(cpdag, g0, z, "X", "Y")
    assert c == b
    assert c > 0


def test_combined_is_the_min_when_both_defined():
    """When both modes are defined, combined must be their min, not either alone.

    X - V - Y, X - V undirected (so the descendant route can pick V -> ...
    forward for free) and V - Y currently ``Y -> V`` and retractable. Z = {V}.
    Mode A: reach V from X -- already adjacent via an undirected edge, cost 0.
    Mode B: the only path is X-V-Y; V is the first internal vertex, so (as
    established above) it can never be a collider -- and V is in Z, so this
    path can never unblock: mode B is UNREACHED here. To get BOTH modes
    defined with different finite values, extend one hop further so the
    back-door path's collider vertex is not the one adjacent to X.
    """
    cpdag = MPDAG(
        ["X", "V", "W", "Y"],
        directed=[("V", "X")],
        undirected=[("V", "W"), ("W", "Y")],
    )
    g0 = MPDAG(
        ["X", "V", "W", "Y"],
        directed=[("V", "X"), ("V", "W"), ("W", "Y")],
        undirected=[],
    )
    z = frozenset({"W"})
    a = descendant_lower_bound(cpdag, g0, z, "X")
    b = backdoor_lower_bound(cpdag, g0, z, "X", "Y")
    c = combined_lower_bound(cpdag, g0, z, "X", "Y")
    assert a != UNREACHED and b != UNREACHED
    assert c == min(a, b)
    assert a != b, "both modes finite but equal would not discriminate min() from either alone"


# --------------------------------------------------------------------------------
# Exhaustive admissibility -- the critical test
# --------------------------------------------------------------------------------


def test_admissibility_exhaustive_n3(tmp_path):
    """Every instance at n=3: L_combined must never exceed the exact BFS radius.

    Full n=4 (185 CPDAGs, thousands of instances) is run and archived under
    ``results/axisb2/bounds/`` separately (too slow to repeat on every pytest
    run); this is the fast, always-on regression guard. Any violation here is
    a bug in this module, not a tuning parameter.
    """
    totals = run_admissibility_study(n=3, out_dir=tmp_path / "n3")
    assert totals.n_cases > 20, "the sweep must actually exercise something"
    assert totals.n_l_violations == 0, totals.violations


def test_admissibility_exhaustive_n4_subset(tmp_path):
    """A real (not tiny) slice of n=4: still zero violations, still non-trivial scope."""
    totals = run_admissibility_study(n=4, out_dir=tmp_path / "n4")
    assert totals.n_cases > 1000
    assert totals.n_l_violations == 0, totals.violations
    # The whole point of this module: the back-door mode must be picking up
    # cases the descendant mode alone could never certify.
    assert totals.n_l_back_only > 0


# --------------------------------------------------------------------------------
# PYTHONHASHSEED invariance
# --------------------------------------------------------------------------------


def test_pythonhashseed_invariance():
    """The n=3 study must be bit-identical across PYTHONHASHSEED values.

    Run as a subprocess (not by mutating ``os.environ`` in-process, which does
    not affect an already-running interpreter's hash randomisation) so this is
    a real test of the claim, not a no-op.
    """
    script = (
        "import json\n"
        "from bkrobust.search.backdoor_bound import run_admissibility_study, totals_as_dict\n"
        "import tempfile, pathlib\n"
        "with tempfile.TemporaryDirectory() as d:\n"
        "    totals = run_admissibility_study(n=3, out_dir=pathlib.Path(d) / 'study')\n"
        "    rows = (pathlib.Path(d) / 'study' / 'rows.csv').read_text()\n"
        "print(json.dumps({'totals': totals_as_dict(totals), 'rows': rows}))\n"
    )
    outputs = []
    for seed in ("0", "1", "1337"):
        env = dict(os.environ)
        env["PYTHONHASHSEED"] = seed
        env["PYTHONPATH"] = os.path.join(REPO_ROOT, "src")
        proc = subprocess.run(
            [sys.executable, "-c", script],
            cwd=REPO_ROOT,
            env=env,
            capture_output=True,
            text=True,
            check=True,
        )
        outputs.append(json.loads(proc.stdout))

    first = outputs[0]
    for other in outputs[1:]:
        assert other["rows"] == first["rows"]
        assert other["totals"] == first["totals"]
