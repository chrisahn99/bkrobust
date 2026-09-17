"""The bias functional: hand-checkable cases, and the two invariants it must have.

Kept fast. The exhaustive falsification of P1-P7 over hundreds of random
instances lives in :mod:`bkrobust.epsilon.verify` and its committed results; what
is here is the part that must stay green on every commit.
"""

from __future__ import annotations

import numpy as np
import pytest

from bkrobust.demo.evaluate import LinearSEM, random_sem
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.epsilon.bias import (
    ZERO_TOL,
    bias_at,
    make_context,
    make_context_from_covariance,
    normalise,
    possible_parent_sets_enumerated,
    possible_parent_sets_semilocal,
    realised_error,
)
from bkrobust.search.space_fixed import enumerate_space_correct


def _chain_sem() -> tuple[MPDAG, LinearSEM]:
    """``C -> X -> Y`` with ``C -> Y``: one confounder, effects known by hand."""
    dag = MPDAG(["C", "X", "Y"], directed=[("C", "X"), ("X", "Y"), ("C", "Y")])
    sem = LinearSEM(
        dag=dag,
        weights={("C", "X"): 2.0, ("X", "Y"): 3.0, ("C", "Y"): 5.0},
        noise_var={"C": 1.0, "X": 1.0, "Y": 1.0},
    )
    return dag, sem


def _pinned_sem() -> tuple[MPDAG, LinearSEM]:
    """A DAG whose CPDAG compels every parent of the treatment.

    ``A -> X <- B`` is an unshielded collider, so both arrows are compelled, and
    R1 then compels ``X -> Y``. ``P -- Q`` is a separate two-node component that
    the data cannot orient, so the space has more than one element while the
    parents of ``X`` are the same in all of them -- which is the situation
    Theorem A describes and a single-element space would not test.
    """
    dag = MPDAG(
        ["A", "B", "P", "Q", "X", "Y"],
        directed=[("A", "X"), ("B", "X"), ("X", "Y"), ("P", "Q")],
    )
    sem = LinearSEM(
        dag=dag,
        weights={("A", "X"): 2.0, ("B", "X"): -1.0, ("X", "Y"): 3.0, ("P", "Q"): 0.7},
        noise_var=dict.fromkeys(["A", "B", "P", "Q", "X", "Y"], 1.0),
    )
    return dag, sem


def test_total_effect_matches_the_path_product() -> None:
    """``tau`` read off the covariance equals the structural path product."""
    dag, sem = _chain_sem()
    cpdag = dag_to_cpdag(dag)
    ctx = make_context_from_covariance(
        sem.covariance(), dag.nodes, cpdag, "X", "Y", frozenset({"C"})
    )
    # Adjusting for the true parent set {C} identifies the X->Y effect, = 3.0.
    assert ctx.total_effect_under(frozenset({"C"})) == pytest.approx(3.0)
    # Not adjusting leaves the confounding path C->X, C->Y open, so it does not.
    assert ctx.total_effect_under(frozenset()) != pytest.approx(3.0)


def test_outcome_as_a_parent_of_treatment_gives_exactly_zero() -> None:
    """The IDA convention: if ``Y`` is a parent of ``X`` then ``X`` does not cause ``Y``."""
    dag, sem = _chain_sem()
    cpdag = dag_to_cpdag(dag)
    ctx = make_context_from_covariance(
        sem.covariance(), dag.nodes, cpdag, "X", "Y", frozenset({"C"})
    )
    assert ctx.total_effect_under(frozenset({"Y"})) == 0.0


def test_bias_is_zero_at_every_state_when_the_parents_are_compelled() -> None:
    """Theorem A: where ``Z`` is valid throughout, the ambiguity set is a singleton.

    Checked at *every* element of the space, not just ``G0``: the claim is about
    the whole certified ball, and the undirected ``P -- Q`` component makes that
    ball non-trivial while leaving the parents of ``X`` alone.
    """
    dag, sem = _pinned_sem()
    cpdag = dag_to_cpdag(dag)
    ctx = make_context(sem, cpdag, "X", "Y", frozenset())
    assert ctx.theta_z == pytest.approx(3.0)
    assert realised_error(ctx) == pytest.approx(0.0, abs=1e-9)

    states = enumerate_space_correct(cpdag)
    assert len(states) > 1, "the space must be non-trivial for this to test anything"
    for g in states:
        got = bias_at(ctx, g)
        assert got.worst == 0.0, g.edge_string()
        assert got.tau_min == pytest.approx(got.tau_max)


def test_an_unoriented_triangle_leaves_real_ambiguity() -> None:
    """The complementary case, so the zero above is not mistaken for a trivial one.

    ``C -> X -> Y`` with ``C -> Y`` has a complete, v-structure-free skeleton, so
    its CPDAG is entirely undirected and the data pin down nothing: adjusting for
    ``C`` is right under the true orientation and wrong under others, and ``B``
    says so rather than reporting zero.
    """
    dag, sem = _chain_sem()
    cpdag = dag_to_cpdag(dag)
    assert not cpdag.directed_edges
    ctx = make_context(sem, cpdag, "X", "Y", frozenset({"C"}))
    got = bias_at(ctx, cpdag)
    assert got.worst > 1.0
    assert got.tau_min < got.tau_max


def test_semilocal_parent_sets_agree_with_enumeration() -> None:
    """Theorem E, differentially tested -- the check that licenses the fast path.

    Run over *every* element of the corrected space of several CPDAGs, not just
    ``G0``: the semi-local rule has to hold everywhere the search will visit it.
    """
    rng = np.random.default_rng(11)
    checked = 0
    for _ in range(8):
        dag = MPDAG(
            [f"V{i}" for i in range(5)],
            directed=[
                (f"V{i}", f"V{j}")
                for i in range(5)
                for j in range(i + 1, 5)
                if rng.uniform() < 0.45
            ],
        )
        cpdag = dag_to_cpdag(dag)
        if not cpdag.undirected_edges:
            continue
        for g in enumerate_space_correct(cpdag):
            for x in dag.nodes:
                assert possible_parent_sets_semilocal(cpdag, g, x) == (
                    possible_parent_sets_enumerated(g, x)
                ), f"semi-local disagreed at {g.edge_string()} for {x}"
                checked += 1
    assert checked > 100, "the differential test did not actually exercise anything"


def test_bias_is_monotone_over_every_ordered_pair() -> None:
    """Theorem B on a full small space: a superset of models cannot lower the max."""
    rng = np.random.default_rng(3)
    dag = MPDAG(
        ["V0", "V1", "V2", "V3"],
        directed=[("V0", "V1"), ("V0", "V2"), ("V1", "V3"), ("V2", "V3")],
    )
    cpdag = dag_to_cpdag(dag)
    sem = random_sem(dag, rng)
    ctx = make_context(sem, cpdag, "V0", "V3", frozenset())
    elems = enumerate_space_correct(cpdag)
    values = {g.edge_string(): bias_at(ctx, g).worst for g in elems}
    from bkrobust.core.oracle import extensions

    reps = {g.edge_string(): frozenset(extensions(g)) for g in elems}
    pairs = 0
    for g in elems:
        for h in elems:
            if g is h or not reps[g.edge_string()] <= reps[h.edge_string()]:
                continue
            pairs += 1
            assert values[g.edge_string()] <= values[h.edge_string()] + 1e-9
    assert pairs > 0


def test_relative_units_refuse_a_zero_estimate() -> None:
    """A relative bias against a zero estimate is undefined and must not become infinity."""
    dag, sem = _chain_sem()
    cpdag = dag_to_cpdag(dag)
    ctx = make_context(sem, cpdag, "X", "Y", frozenset({"C"}))
    zeroed = type(ctx)(**{**ctx.__dict__, "theta_z": 0.0})
    with pytest.raises(ValueError, match="relative units are undefined"):
        normalise(zeroed, 1.0, units="relative")
    assert normalise(ctx, 1.0, units="absolute") == 1.0
    assert normalise(ctx, 1.0, units="standardised") == pytest.approx(ctx.scale_std)


def test_context_rejects_treatment_or_outcome_inside_z() -> None:
    """A set containing the query's own endpoints is a caller error, not a silent no-op."""
    dag, sem = _chain_sem()
    cpdag = dag_to_cpdag(dag)
    with pytest.raises(ValueError, match="must not contain"):
        make_context(sem, cpdag, "X", "Y", frozenset({"X"}))


def test_zero_tolerance_is_applied_not_merely_documented() -> None:
    """Sub-tolerance bias is reported as exactly 0.0, since Theorem A predicts algebraic zero."""
    dag, sem = _pinned_sem()
    cpdag = dag_to_cpdag(dag)
    ctx = make_context(sem, cpdag, "X", "Y", frozenset())
    g0 = apply_orientations(cpdag, [])
    assert g0 is not None
    value = bias_at(ctx, g0).worst
    assert value == 0.0 and not (0.0 < value < ZERO_TOL)
