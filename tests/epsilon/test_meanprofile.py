"""The mean bias profile: what it is, what it is not, and what must not drift.

``mu`` is an average over the retraction shell, so two mistakes would be easy to
make and both would be serious. Reporting it as if it bounded the worst case
would turn an average-case statement into a guarantee it cannot support, and
reporting a sampled estimate without its uncertainty would hide that the number
is a draw rather than a fact. Both have dedicated tests here.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from bkrobust.core.conventions import UNREACHED
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag, random_sem
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.epsilon.bias import BiasContext, knowledge_of, make_context
from bkrobust.epsilon.meanprofile import (
    ShellStat,
    is_monotone,
    mean_bounds,
    r_mean,
    shell_stat,
)
from bkrobust.epsilon.profile import bias_profile


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
    ctx = make_context(random_sem(dag, rng), cpdag, "V2", "V4", frozenset(z))
    return cpdag, g0, frozenset(z), ctx


def _big_instance() -> tuple[MPDAG, MPDAG, frozenset[str], BiasContext]:
    """A denser instance with ``|K_G0| = 9``, so shells are large enough to sample.

    The five-node instance above has only two orientations, which makes every
    shell smaller than any sampling budget worth testing.
    """
    rng = np.random.default_rng(2)
    nodes = [f"V{i}" for i in range(7)]
    edges = [
        (nodes[i], nodes[j])
        for i in range(7)
        for j in range(i + 1, 7)
        if rng.random() < 0.45
    ]
    dag = MPDAG(nodes, directed=edges)
    cpdag = dag_to_cpdag(dag)
    g0 = apply_orientations(cpdag, sorted(set(dag.directed_edges) - set(cpdag.directed_edges)))
    assert g0 is not None
    z = optimal_adjustment_set_mpdag(g0, "V1", "V6")
    assert z is not None
    ctx = make_context(random_sem(dag, np.random.default_rng(11)), cpdag, "V1", "V6", frozenset(z))
    return cpdag, g0, frozenset(z), ctx


def _profile(cpdag: MPDAG, g0: MPDAG, ctx: BiasContext, **kw: object) -> list[ShellStat]:
    """Every shell of one instance, sharing a bias cache."""
    n = len(knowledge_of(cpdag, g0))
    cache: dict[str, float] = {}
    return [shell_stat(ctx, cpdag, g0, d, cache=cache, **kw) for d in range(n + 1)]  # type: ignore[arg-type]


def test_exhaustive_shell_counts_every_subset() -> None:
    """The population is ``C(|K_G0|, d)`` subsets, not the distinct closures.

    Several subsets Meek-close to the same state; each is a distinct way for the
    analyst to be wrong, so each is counted. Collapsing to distinct states would
    silently change the estimand.
    """
    cpdag, g0, _, ctx = _instance()
    n = len(knowledge_of(cpdag, g0))
    for d in range(n + 1):
        s = shell_stat(ctx, cpdag, g0, d)
        assert s.n_subsets_total == math.comb(n, d)
        assert s.n_evaluated == s.n_subsets_total
        assert s.exhaustive
        assert s.n_states <= s.n_evaluated


def test_mean_never_exceeds_the_worst_case() -> None:
    """``mu(d) <= beta_up(d)``. The mean is a companion to the certificate, not one."""
    cpdag, g0, _, ctx = _instance()
    for s in _profile(cpdag, g0, ctx):
        assert s.mean <= s.beta_up + 1e-12


def test_profile_is_zero_below_r_val() -> None:
    """Theorem A forces a maximum of zero, and a zero max forces a zero mean."""
    cpdag, g0, _, ctx = _instance()
    beta = bias_profile(ctx, g0, start=0)
    first_nonzero = next((s.d for s in beta if s.beta_up > 0), None)
    if first_nonzero is None or first_nonzero == 0:
        pytest.skip("instance has no strictly clean shell to check")
    for s in _profile(cpdag, g0, ctx)[:first_nonzero]:
        assert s.mean == 0.0
        assert s.frac_nonzero == 0.0


def test_bounds_bracket_the_mean() -> None:
    """``p*min+ <= mu <= p*beta_up``, the bracket the cheap estimates rest on."""
    cpdag, g0, _, ctx = _instance()
    for s in _profile(cpdag, g0, ctx):
        lo, hi = mean_bounds(s)
        assert lo <= s.mean + 1e-12
        assert s.mean <= hi + 1e-12


def test_exhaustive_has_no_standard_error() -> None:
    """Enumerating the population leaves nothing to be uncertain about."""
    cpdag, g0, _, ctx = _instance()
    for s in _profile(cpdag, g0, ctx):
        assert s.se == 0.0
        assert s.ci() == (s.mean, s.mean)


def test_sampling_reproduces_the_exhaustive_mean_at_full_budget() -> None:
    """A sample at least as large as the population is the population."""
    cpdag, g0, _, ctx = _big_instance()
    n = len(knowledge_of(cpdag, g0))
    rng = np.random.default_rng(0)
    for d in range(n + 1):
        exact = shell_stat(ctx, cpdag, g0, d)
        full = shell_stat(ctx, cpdag, g0, d, sample=math.comb(n, d), rng=rng)
        assert full.exhaustive
        assert full.mean == pytest.approx(exact.mean)


def test_sampled_shell_reports_uncertainty_and_stays_in_range() -> None:
    """A partial draw must carry a standard error and a usable interval."""
    cpdag, g0, _, ctx = _big_instance()
    n = len(knowledge_of(cpdag, g0))
    d = max(1, n // 2)
    total = math.comb(n, d)
    if total < 8:
        pytest.skip("shell too small to sample meaningfully")
    s = shell_stat(ctx, cpdag, g0, d, sample=max(2, total // 4), rng=np.random.default_rng(1))
    assert not s.exhaustive
    assert s.n_evaluated < s.n_subsets_total
    assert s.se >= 0.0
    lo, hi = s.ci()
    assert lo >= 0.0 and lo <= s.mean <= hi


def test_sampled_subsets_are_distinct() -> None:
    """Sampling is without replacement, which the finite population correction assumes."""
    cpdag, g0, _, ctx = _big_instance()
    n = len(knowledge_of(cpdag, g0))
    d = max(1, n // 2)
    total = math.comb(n, d)
    if total < 8:
        pytest.skip("shell too small to sample meaningfully")
    want = max(2, total // 3)
    s = shell_stat(ctx, cpdag, g0, d, sample=want, rng=np.random.default_rng(2))
    assert s.n_evaluated == want


def test_radius_is_a_first_crossing_and_reports_unreached() -> None:
    """Above the ceiling nothing crosses, and that is a status rather than a number."""
    cpdag, g0, _, ctx = _instance()
    prof = _profile(cpdag, g0, ctx)
    ceiling = max(s.mean for s in prof)
    assert r_mean(prof, ceiling * 10 + 1.0) == UNREACHED
    r = r_mean(prof, 0.0)
    if r != UNREACHED:
        assert prof[r].mean > 0.0
        assert all(s.mean <= 0.0 for s in prof[:r])


def test_radius_is_non_decreasing_in_epsilon() -> None:
    """A looser tolerance can only be crossed later, or not at all."""
    cpdag, g0, _, ctx = _instance()
    prof = _profile(cpdag, g0, ctx)
    seen = [r_mean(prof, e) for e in (0.0, 0.01, 0.1, 0.5, 1.0, 5.0)]
    finite = [(i, v) for i, v in enumerate(seen) if v != UNREACHED]
    for (i, a), (_, b) in zip(finite, finite[1:]):  # noqa: B905 - pairwise
        assert b >= a, f"radius fell as epsilon rose: {seen}"


def test_confidence_bounds_order_the_radius() -> None:
    """``ci_lo`` crosses no earlier than the mean, so it gives the larger radius."""
    stats = [
        ShellStat(d=0, n_subsets_total=1, n_evaluated=1, n_states=1, mean=0.0, se=0.0,
                  beta_up=0.0, min_nonzero=0.0, frac_nonzero=0.0, exhaustive=True),
        ShellStat(d=1, n_subsets_total=10, n_evaluated=5, n_states=3, mean=0.50, se=0.10,
                  beta_up=1.0, min_nonzero=0.5, frac_nonzero=0.5, exhaustive=False),
        ShellStat(d=2, n_subsets_total=10, n_evaluated=5, n_states=4, mean=0.90, se=0.10,
                  beta_up=1.0, min_nonzero=0.5, frac_nonzero=0.9, exhaustive=False),
    ]
    assert r_mean(stats, 0.55, use="ci_hi") <= r_mean(stats, 0.55, use="mean")
    assert r_mean(stats, 0.55, use="mean") <= r_mean(stats, 0.55, use="ci_lo")
    with pytest.raises(ValueError):
        r_mean(stats, 0.5, use="nonsense")


def test_is_monotone_detects_a_fall() -> None:
    """Monotonicity of ``mu`` is observed, not proved, so it is checked and reported."""
    up = [
        ShellStat(d=0, n_subsets_total=1, n_evaluated=1, n_states=1, mean=0.0, se=0.0,
                  beta_up=0.0, min_nonzero=0.0, frac_nonzero=0.0, exhaustive=True),
        ShellStat(d=1, n_subsets_total=3, n_evaluated=3, n_states=2, mean=0.4, se=0.0,
                  beta_up=1.0, min_nonzero=0.4, frac_nonzero=0.4, exhaustive=True),
    ]
    assert is_monotone(up)
    down = [*up, ShellStat(d=2, n_subsets_total=3, n_evaluated=3, n_states=2, mean=0.2,
                           se=0.0, beta_up=1.0, min_nonzero=0.2, frac_nonzero=0.2,
                           exhaustive=True)]
    assert not is_monotone(down)


def test_out_of_range_depth_and_missing_rng_are_refused() -> None:
    """Failing loudly beats returning a shell that is not the one asked for."""
    cpdag, g0, _, ctx = _instance()
    n = len(knowledge_of(cpdag, g0))
    with pytest.raises(ValueError):
        shell_stat(ctx, cpdag, g0, n + 1)
    with pytest.raises(ValueError):
        shell_stat(ctx, cpdag, g0, -1)
    cpdag, g0, _, ctx = _big_instance()
    n = len(knowledge_of(cpdag, g0))
    assert math.comb(n, 2) > 1, "fixture must have a shell larger than one subset"
    with pytest.raises(ValueError, match="rng"):
        shell_stat(ctx, cpdag, g0, 2, sample=1)
