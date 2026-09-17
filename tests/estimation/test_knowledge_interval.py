"""The knowledge interval's enumeration, against brute-force DAG extensions.

Every shortcut the interval takes is checked here on small random graphs where
the answer can be computed the slow way: the parent sets of the treatment over
every DAG extension of the full analyst graph, which
:func:`~bkrobust.demo.meek.enumerate_dag_extensions` gives by trying all
orientations. The shortcuts under test are the restriction to the treatment's
chain component, the semi-local enumeration over cliques of its undirected
neighbours, the invariance of the ball to claims outside the component, and the
direction restriction through possible descendants.
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from bkrobust.demo.evaluate import LinearSEM, adjusted_estimand
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.estimation.knowledge_interval import (
    chain_component,
    closure_by_components,
    component_graph,
    local_claims,
    possible_parent_sets,
    regression_beta_se,
    retraction_ball,
    y_possible_descendant,
)


def _random_dag(rng: np.random.Generator, n: int, p: float) -> MPDAG:
    order = list(rng.permutation(n))
    edges = [
        (f"v{order[i]}", f"v{order[j]}")
        for i in range(n)
        for j in range(i + 1, n)
        if rng.random() < p
    ]
    return MPDAG(nodes=[f"v{i}" for i in range(n)], directed=edges)


def _instances(seed: int, count: int) -> list[tuple[MPDAG, MPDAG, list]]:
    """Random (DAG, CPDAG, claims) with some claims reversed, Meek-consistent only."""
    rng = np.random.default_rng(seed)
    out: list[tuple[MPDAG, MPDAG, list]] = []
    while len(out) < count:
        dag = _random_dag(rng, int(rng.integers(5, 9)), float(rng.uniform(0.25, 0.5)))
        cpdag = dag_to_cpdag(dag)
        und = sorted(cpdag.undirected_edges)
        if not 1 <= len(und) <= 9:
            continue
        pick = rng.permutation(len(und))[: int(rng.integers(0, len(und) + 1))]
        claims = []
        for t in pick:
            a, b = und[int(t)]
            e = (a, b) if dag.is_directed_edge(a, b) else (b, a)
            claims.append(e if rng.random() < 0.7 else (e[1], e[0]))
        if apply_orientations(cpdag, claims) is None:
            continue
        out.append((dag, cpdag, claims))
    return out


def _brute_parent_sets(g: MPDAG, x: str, y: str) -> tuple[set, set]:
    """Parent sets of ``x`` over all extensions, and those with ``y`` a descendant."""
    every, towards = set(), set()
    for d in enumerate_dag_extensions(g):
        pa = frozenset(d.parents(x))
        every.add(pa)
        if y in d.descendants(x):
            towards.add(pa)
    return every, towards


def _semi_local(cpdag: MPDAG, claims: list, x: str, y: str) -> tuple[set, set] | None:
    comp = chain_component(cpdag, x)
    outside = frozenset(cpdag.parents(x))
    if comp is None:
        return {outside}, ({outside} if y in _possde_full(cpdag, claims, x) else set())
    base = component_graph(cpdag, comp)
    k_loc = local_claims(cpdag, claims, comp)
    local = apply_orientations(base, k_loc) if k_loc else base
    if local is None:
        return None
    ps = possible_parent_sets(local, x)
    assert not ps.censored
    every = {outside | p for p in ps.local}
    towards = {
        outside | p
        for p, closed in zip(ps.local, ps.closed, strict=True)
        if y_possible_descendant(cpdag, comp, closed, x, y)
    }
    return every, towards


def _possde_full(cpdag: MPDAG, claims: list, x: str) -> set:
    from bkrobust.mpdag_criterion.paths import possible_descendants

    g = apply_orientations(cpdag, claims)
    assert g is not None
    return possible_descendants(g, x)


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_semi_local_parent_sets_equal_brute_force(seed: int) -> None:
    checked = ambiguous = 0
    for dag, cpdag, claims in _instances(seed, 40):
        g = apply_orientations(cpdag, claims)
        assert g is not None
        for x in dag.nodes:
            y = sorted(set(dag.nodes) - {x})[0]
            every, _ = _brute_parent_sets(g, x, y)
            got = _semi_local(cpdag, claims, x, y)
            assert got is not None
            assert got[0] == every, (g, x, got[0], every)
            checked += 1
            ambiguous += len(every) > 1
    # the comparison is not vacuous: many treatments have several parent sets
    assert checked > 100
    assert ambiguous > 20


@pytest.mark.parametrize("seed", [4, 5])
def test_direction_restriction_equals_brute_force(seed: int) -> None:
    strict = 0
    for dag, cpdag, claims in _instances(seed, 30):
        g = apply_orientations(cpdag, claims)
        assert g is not None
        for x, y in itertools.permutations(dag.nodes, 2):
            every, towards = _brute_parent_sets(g, x, y)
            got = _semi_local(cpdag, claims, x, y)
            assert got is not None
            assert got[1] == towards, (g, x, y, got[1], towards)
            strict += bool(towards) and towards != every
    assert strict > 20


def test_local_closure_agrees_with_full_closure_on_consistency() -> None:
    rng = np.random.default_rng(7)
    seen = 0
    while seen < 150:
        dag = _random_dag(rng, int(rng.integers(5, 9)), 0.4)
        cpdag = dag_to_cpdag(dag)
        und = sorted(cpdag.undirected_edges)
        if not und:
            continue
        claims = [e if rng.random() < 0.5 else (e[1], e[0]) for e in und[: rng.integers(1, 6)]]
        full = apply_orientations(cpdag, claims)
        locals_ok = True
        for comp in {chain_component(cpdag, a) for a, _ in claims}:
            assert comp is not None
            k_loc = local_claims(cpdag, claims, comp)
            if apply_orientations(component_graph(cpdag, comp), k_loc) is None:
                locals_ok = False
        assert (full is not None) == locals_ok
        seen += 1


def test_closure_by_components_equals_full_closure() -> None:
    rng = np.random.default_rng(13)
    seen = consistent = 0
    while seen < 300:
        dag = _random_dag(rng, int(rng.integers(5, 10)), float(rng.uniform(0.25, 0.5)))
        cpdag = dag_to_cpdag(dag)
        und = sorted(cpdag.undirected_edges)
        if not und:
            continue
        pick = rng.permutation(len(und))[: int(rng.integers(1, len(und) + 1))]
        claims = [und[int(t)] if rng.random() < 0.5 else und[int(t)][::-1] for t in pick]
        full = apply_orientations(cpdag, claims)
        assembled = closure_by_components(cpdag, claims, cache={})
        assert full == assembled, (cpdag, claims, full, assembled)
        consistent += full is not None
        seen += 1
    assert consistent > 50


def test_ball_is_monotone_and_ignores_outside_claims() -> None:
    for dag, cpdag, claims in _instances(11, 40):
        for x in dag.nodes:
            comp = chain_component(cpdag, x)
            if comp is None:
                continue
            base = component_graph(cpdag, comp)
            k_loc = local_claims(cpdag, claims, comp)
            ball = retraction_ball(base, k_loc, max_depth=3, subset_budget=10_000)
            by_depth: dict[int, set] = {}
            for el in ball.elements:
                # retracting from a consistent claim set never makes it inconsistent
                assert el.graph is not None
                sets = set(possible_parent_sets(el.graph, x).local)
                by_depth.setdefault(len(el.retracted), set()).update(sets)
            depths = sorted(by_depth)
            for shallow, deep in itertools.pairwise(depths):
                # the sets reachable at one depth contain every shallower one
                assert by_depth[shallow] <= by_depth[deep]
            # retracting a claim in another component changes nothing here
            outside = [e for e in claims if e not in k_loc]
            for e in outside:
                kept = [c for c in claims if c != e]
                assert local_claims(cpdag, kept, comp) == k_loc


def test_population_effect_of_true_parents_is_the_path_sum() -> None:
    rng = np.random.default_rng(3)
    for _ in range(20):
        dag = _random_dag(rng, 7, 0.4)
        weights = {
            e: float(rng.choice([-1, 1]) * rng.uniform(0.5, 1.5))
            for e in sorted(dag.directed_edges)
        }
        sem = LinearSEM(dag=dag, weights=weights, noise_var={v: 1.0 for v in dag.nodes})
        sigma = sem.covariance()[None, :, :]
        idx = {v: i for i, v in enumerate(dag.nodes)}
        for x, y in itertools.permutations(dag.nodes, 2):
            pa = dag.parents(x)
            beta, se = regression_beta_se(sigma, idx, x, y, pa, None)
            truth = sem.true_total_effect(x, y) if y not in pa else 0.0
            assert beta[0] == pytest.approx(truth, abs=1e-9)
            assert se[0] == 0.0
            assert beta[0] == pytest.approx(
                adjusted_estimand(sem, x, y, pa) if y not in pa else 0.0, abs=1e-9
            )


def test_standard_error_matches_ordinary_least_squares() -> None:
    rng = np.random.default_rng(5)
    n = 400
    data = rng.normal(size=(n, 4))
    data[:, 1] += 0.8 * data[:, 0]
    data[:, 3] += 0.5 * data[:, 1] - 0.7 * data[:, 2]
    idx = {"a": 0, "x": 1, "c": 2, "y": 3}
    sigma = np.cov(data, rowvar=False)[None, :, :]
    beta, se = regression_beta_se(sigma, idx, "x", "y", ["a", "c"], n)
    design = np.column_stack([np.ones(n), data[:, 1], data[:, 0], data[:, 2]])
    coef, *_ = np.linalg.lstsq(design, data[:, 3], rcond=None)
    resid = data[:, 3] - design @ coef
    s2 = resid @ resid / (n - design.shape[1])
    cov = s2 * np.linalg.inv(design.T @ design)
    assert beta[0] == pytest.approx(coef[1], rel=1e-9)
    assert se[0] == pytest.approx(np.sqrt(cov[1, 1]), rel=1e-9)
