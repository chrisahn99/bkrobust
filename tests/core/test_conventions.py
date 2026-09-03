"""Pin the normative radius convention and the frozen core against known values.

If any of these fail, every number produced this session shifts. They are the
reason the convention is written down rather than left implicit.
"""

from __future__ import annotations

import pytest

from bkrobust.core.conventions import (
    DEGENERATE,
    RADIUS_CONVENTION,
    UNREACHED,
    is_certified_clean,
    safe_moves,
)
from bkrobust.core.oracle import is_valid
from bkrobust.core.spacelib import Space, build_space, distances_from, radius
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.scenario import OUTCOME, TREATMENT, scenarios, true_dag


def _setup(label: str) -> tuple[Space, MPDAG, frozenset[str], dict[MPDAG, int]]:
    truth = true_dag()
    cpdag = dag_to_cpdag(truth)
    space = build_space(cpdag)
    g0 = apply_orientations(cpdag, scenarios()[label]["knowledge"])
    z = frozenset(optimal_adjustment_set_mpdag(g0, TREATMENT, OUTCOME) or ())
    dists = distances_from(space, g0)
    return space, g0, z, dists


# --- the convention itself --------------------------------------------------


def test_safe_moves_is_radius_minus_one():
    """The practitioner reading is r-1, not r. Getting this wrong is off-by-one."""
    assert safe_moves(3) == 2
    assert safe_moves(1) == 0


def test_safe_moves_rejects_unreached():
    """UNREACHED has no finite safe-move count and must not silently become -2."""
    with pytest.raises(ValueError):
        safe_moves(UNREACHED)


def test_certified_clean_covers_shells_below_radius():
    assert is_certified_clean(3, 0) and is_certified_clean(3, 2)
    assert not is_certified_clean(3, 3)
    assert is_certified_clean(UNREACHED, 99), "nothing fails, so every shell is clean"


def test_convention_string_is_recorded():
    """Manifests embed this, so results files stay self-describing."""
    assert "shells 0..r-1" in RADIUS_CONVENTION
    assert "safe-moves = r - 1" in RADIUS_CONVENTION


def test_sentinels_are_distinct():
    assert UNREACHED != DEGENERATE


# --- the convention applied -------------------------------------------------


@pytest.mark.parametrize(("label", "expected"), [("A", 3), ("B", 3), ("C", 2)])
def test_core_reproduces_published_radii(label, expected):
    """The frozen core must agree with the numbers already in report.md.

    Under the chosen convention these carry over unchanged; had the other
    convention been adopted every one of them would have shifted by one.
    """
    space, _g0, z, dists = _setup(label)
    r, _ = radius(space, dists, lambda g: not is_valid(z, g, TREATMENT, OUTCOME))
    assert r == expected


@pytest.mark.parametrize("label", ["A", "B", "C"])
def test_shells_below_radius_are_genuinely_clean(label):
    """THE convention, asserted directly: 0..r-1 clean, and r really does fail.

    This is what distinguishes the two candidate conventions, so it is checked
    element by element rather than trusted.
    """
    space, _g0, z, dists = _setup(label)
    r, witness = radius(space, dists, lambda g: not is_valid(z, g, TREATMENT, OUTCOME))
    assert r != UNREACHED and witness is not None

    for g, d in dists.items():
        if d < r:
            assert is_valid(z, g, TREATMENT, OUTCOME), f"failure inside certified shell {d}"
    assert not is_valid(z, witness, TREATMENT, OUTCOME)
    assert dists[witness] == r


@pytest.mark.parametrize("label", ["A", "B", "C"])
def test_radius_at_least_one(label):
    """Z is read off G0, so it is valid there and the radius cannot be 0."""
    space, g0, z, dists = _setup(label)
    r, _ = radius(space, dists, lambda g: not is_valid(z, g, TREATMENT, OUTCOME))
    assert r >= 1
    assert is_valid(z, g0, TREATMENT, OUTCOME)


def test_distance_to_self_is_zero():
    _space, g0, _z, dists = _setup("B")
    assert dists[g0] == 0


def test_unreached_when_nothing_fails():
    """A predicate that never fires yields UNREACHED, not a large number."""
    space, _g0, _z, dists = _setup("B")
    r, w = radius(space, dists, lambda g: False)
    assert r == UNREACHED and w is None
