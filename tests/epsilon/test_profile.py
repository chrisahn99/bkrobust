"""The staircase, the radius, and the guarantees each strategy is supposed to carry."""

from __future__ import annotations

import numpy as np
import pytest

from bkrobust.core.conventions import UNREACHED
from bkrobust.core.oracle import extensions
from bkrobust.core.spacelib import distances_from
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag, random_sem
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.epsilon.bias import BiasContext, bias_at, knowledge_of, make_context
from bkrobust.epsilon.profile import (
    bias_profile,
    certified_band,
    epsilon_grid,
    greedy_chain_bound,
    r_epsilon,
    r_epsilon_from_profile,
    retraction_shell,
    top_bias,
)
from bkrobust.hybrid import breakdown_radius
from bkrobust.search.space_fixed import build_space_correct


def _instance(seed: int = 5) -> tuple[MPDAG, MPDAG, frozenset[str], BiasContext]:
    """A screened instance with a finite ``r_val`` and a non-trivial space."""
    rng = np.random.default_rng(seed)
    dag = MPDAG(
        ["V0", "V1", "V2", "V3", "V4"],
        directed=[
            ("V0", "V1"),
            ("V0", "V2"),
            ("V1", "V3"),
            ("V2", "V3"),
            ("V3", "V4"),
            ("V2", "V4"),
        ],
    )
    cpdag = dag_to_cpdag(dag)
    g0 = apply_orientations(cpdag, sorted(set(dag.directed_edges) - set(cpdag.directed_edges)))
    assert g0 is not None
    z = optimal_adjustment_set_mpdag(g0, "V2", "V4")
    assert z is not None
    sem = random_sem(dag, rng)
    ctx = make_context(sem, cpdag, "V2", "V4", frozenset(z))
    return cpdag, g0, frozenset(z), ctx


def test_retraction_shell_is_deduplicated_and_bounded() -> None:
    """Distinct subsets can share a closure, so the shell is smaller than ``C(|K|, d)``."""
    cpdag, g0, _, _ = _instance()
    k0 = knowledge_of(cpdag, g0)
    assert len(k0) >= 2
    for d in range(len(k0) + 1):
        shell = retraction_shell(cpdag, g0, d)
        strings = [g.edge_string() for g in shell]
        assert len(strings) == len(set(strings)), "shell contains duplicates"
    assert retraction_shell(cpdag, g0, 0) == [g0]
    assert retraction_shell(cpdag, g0, len(k0)) == [cpdag]
    with pytest.raises(ValueError):
        retraction_shell(cpdag, g0, len(k0) + 1)


def test_shell_cap_raises_rather_than_truncating() -> None:
    """A silently partial shell would give a radius that is too large -- the unsafe direction."""
    cpdag, g0, _, _ = _instance()
    k0 = knowledge_of(cpdag, g0)
    if len(k0) < 2:
        pytest.skip("needs at least two orientations to overflow a cap of 1")
    with pytest.raises(MemoryError, match="above the cap"):
        retraction_shell(cpdag, g0, 1, cap=0)


def test_profile_is_non_decreasing_and_zero_below_r_val() -> None:
    """Theorem C.1 and Theorem A, on the computed staircase."""
    cpdag, g0, z, ctx = _instance()
    r_val = breakdown_radius(cpdag, None, "V2", "V4", z, g0=g0).radius
    profile = bias_profile(ctx, g0, start=0)
    values = [s.beta_up for s in profile]
    assert values == sorted(values), "beta_up must be non-decreasing"
    if r_val != UNREACHED:
        for step in profile:
            if step.d < r_val:
                assert step.beta_up == 0.0


def test_shell_max_is_already_the_ball_max() -> None:
    """Theorem C.3: the maximum over the ball is attained on its sphere."""
    _cpdag, g0, _, ctx = _instance()
    profile = bias_profile(ctx, g0, start=0)
    raw = [s.shell_max for s in profile if s.n_states > 0]
    assert raw == sorted(raw), (
        "the raw per-shell maximum is itself non-decreasing; if it were not, "
        "sphere sufficiency would be false"
    )


def test_strategies_agree_with_each_other_and_with_brute_force() -> None:
    """Theorem C.4/C.5: three code paths, one answer.

    Brute force here is a BFS over the *whole* corrected space with the predicate
    ``B > eps`` -- it does not know about up-sets, retraction, or the staircase,
    so agreement is a genuine cross-check rather than a tautology.
    """
    cpdag, g0, _, ctx = _instance()
    space = build_space_correct(cpdag)
    dists = distances_from(space, g0)
    cache: dict[str, float] = {}

    def brute(eps: float) -> int:
        best = UNREACHED
        for g in space.elements:
            d = dists.get(g)
            if d is None or (best != UNREACHED and d >= best):
                continue
            key = g.edge_string()
            if key not in cache:
                cache[key] = bias_at(ctx, g).worst
            if cache[key] > eps:
                best = d
        return best

    ceiling = top_bias(ctx).worst
    for frac in (0.0, 0.1, 0.5, 0.9, 1.5):
        eps = frac * ceiling
        want = brute(eps)
        for strategy in ("incremental", "bisection", "hybrid"):
            got = r_epsilon(ctx, g0, eps, strategy=strategy).radius
            assert got == want, f"{strategy} gave {got}, brute force gave {want} at eps={eps}"


def test_radius_is_non_decreasing_in_epsilon_and_at_least_r_val() -> None:
    """Theorem C.6, including the ordering against the validity radius."""
    cpdag, g0, z, ctx = _instance()
    r_val = breakdown_radius(cpdag, None, "V2", "V4", z, g0=g0).radius
    ceiling = top_bias(ctx).worst
    grid = [f * ceiling for f in (0.0, 0.2, 0.4, 0.6, 0.8)]
    radii, _ = epsilon_grid(ctx, g0, grid, r_val=max(0, r_val))
    ordered = [radii[e] for e in grid]
    ranked = [float("inf") if r == UNREACHED else r for r in ordered]
    assert ranked == sorted(ranked), "r_eps must be non-decreasing in eps"
    if r_val != UNREACHED:
        for r in ordered:
            assert r == UNREACHED or r >= r_val


def test_infinity_test_short_circuits_above_the_ceiling() -> None:
    """Bound U2: one evaluation at the top settles every eps above ``B(Chat)``."""
    _cpdag, g0, _, ctx = _instance()
    ceiling = top_bias(ctx).worst
    got = r_epsilon(ctx, g0, ceiling * 1.01, strategy="hybrid")
    assert got.radius == UNREACHED
    assert got.strategy == "top_only"
    assert got.exact
    assert got.states_evaluated == 1


def test_greedy_chain_is_an_upper_bound_and_says_so() -> None:
    """Bound U1: a witness bounds the radius above, and is never reported as exact."""
    _cpdag, g0, _, ctx = _instance()
    ceiling = top_bias(ctx).worst
    eps = 0.5 * ceiling
    exact = r_epsilon(ctx, g0, eps, strategy="incremental").radius
    chain = greedy_chain_bound(ctx, g0, eps)
    assert not chain.exact
    if chain.radius != UNREACHED and exact != UNREACHED:
        assert chain.radius >= exact, "a chain witness cannot beat the exact radius"


def test_grid_lookup_matches_a_per_epsilon_search() -> None:
    """One traversal, every eps: the lookup must equal the search it replaces."""
    _cpdag, g0, _, ctx = _instance()
    ceiling = top_bias(ctx).worst
    grid = [f * ceiling for f in (0.1, 0.3, 0.7)]
    radii, profile = epsilon_grid(ctx, g0, grid)
    for eps in grid:
        assert radii[eps] == r_epsilon_from_profile(profile, eps)
        assert radii[eps] == r_epsilon(ctx, g0, eps, strategy="incremental").radius


def test_certified_band_brackets_the_budget() -> None:
    """The band is (largest reachable eps, smallest unreachable eps]."""
    radii = {0.01: 2, 0.05: 2, 0.10: 4, 0.25: UNREACHED}
    assert certified_band(radii, 1) == (None, 0.01)
    assert certified_band(radii, 2) == (0.05, 0.10)
    assert certified_band(radii, 4) == (0.10, 0.25)
    assert certified_band(radii, 99) == (0.10, 0.25)


def test_extensions_are_never_empty_in_the_space() -> None:
    """The premise Theorem A leans on, checked rather than assumed."""
    cpdag, _, _, _ = _instance()
    for g in build_space_correct(cpdag).elements:
        assert extensions(g), g.edge_string()
