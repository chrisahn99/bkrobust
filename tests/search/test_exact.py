"""The accelerated exact search must reproduce BFS radii, always.

Any disagreement is a bug. These are the differential tests the validation
protocol requires, run at pytest scale; the full exhaustive sweep lives in
`bkrobust.search.differential` and its results are committed under
`results/search/`.
"""

from __future__ import annotations

import itertools

import pytest

from bkrobust.core.oracle import is_valid
from bkrobust.core.spacelib import build_space, distances_from, radius
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.scenario import OUTCOME, TREATMENT, scenarios, true_dag
from bkrobust.search.conjecture_study import all_cpdags
from bkrobust.search.exact import (
    SearchStats,
    descendant_lower_bound,
    local_up_covers,
    radius_local_up,
    retractable_edges,
)


@pytest.mark.parametrize("label", ["A", "B", "C"])
def test_matches_bfs_on_worked_example(label):
    """The published radii 3/3/2 must come back from the space-free method too."""
    cpdag = dag_to_cpdag(true_dag())
    space = build_space(cpdag)
    g0 = apply_orientations(cpdag, scenarios()[label]["knowledge"])
    z = frozenset(optimal_adjustment_set_mpdag(g0, TREATMENT, OUTCOME) or ())

    def fails(g: MPDAG) -> bool:
        return not is_valid(z, g, TREATMENT, OUTCOME)

    r_bfs, _ = radius(space, distances_from(space, g0), fails)
    res = radius_local_up(cpdag, g0, fails)
    assert res.radius == r_bfs
    assert res.exact


def test_differential_against_bfs_on_all_three_node_cpdags():
    """Exhaustive differential test at n=3. Small, but it must be perfect."""
    n_cases = 0
    for cpdag in all_cpdags(3):
        if not cpdag.undirected_edges:
            continue
        space = build_space(cpdag)
        for g0 in space.elements:
            dists = distances_from(space, g0)
            for x, y in itertools.permutations(list(cpdag.nodes), 2):
                o = optimal_adjustment_set_mpdag(g0, x, y)
                if o is None:
                    continue
                z = frozenset(o)
                if not is_valid(z, g0, x, y):
                    continue

                def fails(g: MPDAG, _z: frozenset[str] = z, _x: str = x, _y: str = y) -> bool:
                    return not is_valid(_z, g, _x, _y)

                r_bfs, _ = radius(space, dists, fails)
                assert radius_local_up(cpdag, g0, fails).radius == r_bfs
                n_cases += 1
    assert n_cases > 50, "the sweep must actually exercise something"


def test_local_covers_match_space_derived_covers():
    """Locally generated upper covers must equal the covers the space knows about.

    This is the correctness risk in avoiding space enumeration: if local
    generation missed a cover, distances would come out too large.
    """
    cpdag = dag_to_cpdag(true_dag())
    space = build_space(cpdag)
    for g in space.elements:
        from_space = {h.edge_string() for (lo, h) in space.covers if lo == g}
        from_local = {h.edge_string() for h in local_up_covers(cpdag, g)}
        assert from_local == from_space, g.edge_string()


def test_retractable_excludes_compelled_edges():
    """Data-compelled edges are frozen and must never be offered as retractable."""
    cpdag = dag_to_cpdag(true_dag())
    g0 = apply_orientations(cpdag, scenarios()["B"]["knowledge"])
    retr = set(retractable_edges(cpdag, g0))
    for a, b in cpdag.directed_edges:
        assert (a, b) not in retr


def test_anytime_budget_returns_a_bound_not_a_wrong_answer():
    """Stopping early must yield a labelled bound, never an exact-looking number."""
    cpdag = dag_to_cpdag(true_dag())
    g0 = apply_orientations(cpdag, scenarios()["B"]["knowledge"])
    z = frozenset(optimal_adjustment_set_mpdag(g0, TREATMENT, OUTCOME) or ())

    def fails(g: MPDAG) -> bool:
        return not is_valid(z, g, TREATMENT, OUTCOME)

    res = radius_local_up(cpdag, g0, fails, max_depth=1)
    assert not res.exact
    assert res.method.endswith("budget")


def test_stats_are_counted():
    cpdag = dag_to_cpdag(true_dag())
    g0 = apply_orientations(cpdag, scenarios()["B"]["knowledge"])
    z = frozenset(optimal_adjustment_set_mpdag(g0, TREATMENT, OUTCOME) or ())
    st = SearchStats()
    radius_local_up(cpdag, g0, lambda g: not is_valid(z, g, TREATMENT, OUTCOME), stats=st)
    assert st.closures > 0 and st.validity_checks > 0


def test_descendant_lower_bound_is_zero_when_z_already_reachable():
    """Sanity: a Z member already downstream of X costs nothing to reach."""
    g = MPDAG(["X", "M", "Y"], directed=[("X", "M"), ("M", "Y")])
    cpdag = dag_to_cpdag(g)
    assert descendant_lower_bound(cpdag, g, frozenset({"M"}), "X") == 0
