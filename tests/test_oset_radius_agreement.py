"""Regression test for check 5: r_backdoor vs r_complete for O(G0).

Fast on purpose -- the full exhaustive-plus-sample sweep behind
``results/axisa2/oset_radius_agreement.json`` takes tens of seconds and is run
by ``python -m bkrobust.gac.sweep --check5``, not by the test suite. This file
checks two things cheaply:

1. The reviewer's counterexample (``Z`` valid but not optimal) genuinely makes
   the two radii differ: r_backdoor = 1, r_complete = 2.
2. On the exhaustive 3-node scope, where ``Z`` *is* the optimal set O(G0), the
   two radii always agree -- the paper's actual claim.
"""

from __future__ import annotations

from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag
from bkrobust.gac.sweep import _corrected_space, oset_radius, run_check5
from bkrobust.mpdag_criterion.criterion import is_valid_mpdag


def test_reviewer_counterexample_radii_diverge():
    """Z = {W} is valid but not optimal on the chain W-X-Y; radii must differ."""
    cpdag = MPDAG(["W", "X", "Y"], directed=[], undirected=[("W", "X"), ("X", "Y")])
    g0 = MPDAG(["W", "X", "Y"], directed=[("W", "X"), ("X", "Y")], undirected=[])
    space = _corrected_space(cpdag)
    x, y, z = "X", "Y", frozenset({"W"})

    # Z is valid, but not the optimal set -- O(G0) is empty here.
    assert is_valid_mpdag(g0, x, y, z)
    assert is_gac_valid_mpdag(g0, x, y, z)
    assert optimal_adjustment_set_mpdag(g0, x, y) == set()

    r_backdoor = oset_radius(cpdag, space, g0, x, y, z, is_valid_mpdag)
    r_complete = oset_radius(cpdag, space, g0, x, y, z, is_gac_valid_mpdag)

    assert r_backdoor == 1
    assert r_complete == 2


def test_run_check5_counterexample_record_matches():
    """run_check5's own asserted-and-recorded counterexample agrees with the above."""
    record = run_check5(sizes=(3,), n5_step=10**9, n5_max_undirected=0)
    counter = record["counterexample"]
    assert counter["r_backdoor"] == 1
    assert counter["r_complete"] == 2
    assert counter["optimal_set_at_g0"] == []


def test_oset_radius_agrees_on_exhaustive_n3_scope():
    """On O(G0), the closed-form radius agrees under both predicates -- n=3 exhaustive.

    ``n5_step``/``n5_max_undirected`` are set so the (expensive) 5-node sample
    scope is empty, keeping this test fast; the full sweep, including n=4 and
    the 5-node sample, is run separately by ``sweep.main_check5``.
    """
    record = run_check5(sizes=(3,), n5_step=10**9, n5_max_undirected=0)
    scope = record["by_scope"]["n3_n4_exhaustive"]

    assert scope["n_cpdags"] == 7
    assert scope["n_queries"] == 96
    assert scope["radius_equal"] == 96
    assert scope["radius_differ"] == 0
    assert scope["backdoor_invalid_at_g0_but_gac_valid"] == 0
    assert record["radius_disagreement_examples"] == []

    n5_scope = record["by_scope"]["n5_sample"]
    assert n5_scope["n_cpdags"] == 0
