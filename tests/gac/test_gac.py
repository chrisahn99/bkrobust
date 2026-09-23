"""Differential and property tests for the generalised adjustment criterion."""

from __future__ import annotations

import itertools
import time

import pytest

from bkrobust.core.oracle import extensions
from bkrobust.demo.evaluate import is_valid_adjustment_set_dag
from bkrobust.demo.graph import MPDAG
from bkrobust.gac import (
    GAC_REASONS,
    causal_nodes_dag,
    clear_cache,
    forbidden_set_dag,
    gac_forbidden_set,
    is_gac_valid_dag,
    is_gac_valid_dag_by_paths,
    is_gac_valid_mpdag,
    why_invalid_gac,
)
from bkrobust.gac.sweep import (
    candidate_sets,
    classify_disagreement,
    dense_instance,
    enumerated_gac_semantics,
    mpdag_scope,
    run_check1,
    time_calls,
)
from bkrobust.mpdag_criterion import criterion as criterion_module
from bkrobust.mpdag_criterion import paths as paths_module
from bkrobust.mpdag_criterion.criterion import backdoor_forbidden_set, is_valid_mpdag
from bkrobust.search.conjecture_study import all_dags

# --------------------------------------------------------------------------------
# Scope helpers, kept small so the whole module runs in well under 90 s
# --------------------------------------------------------------------------------

DAGS4 = all_dags(4)


def _dag_cases(dags, max_z=3):  # noqa: ANN202 - a plain test-local generator
    for dag in dags:
        for x, y in itertools.permutations(sorted(dag.nodes), 2):
            for z in candidate_sets(dag.nodes, x, y, max_z):
                yield dag, x, y, z


def _mpdag_cases(graphs, max_z=3):  # noqa: ANN202 - a plain test-local generator
    for g in graphs:
        for x, y in itertools.permutations(sorted(g.nodes), 2):
            for z in candidate_sets(g.nodes, x, y, max_z):
                yield g, x, y, z


# --------------------------------------------------------------------------------
# The contrast case: what the GAC buys over back-door
# --------------------------------------------------------------------------------


def test_the_x_to_y_x_to_w_contrast_case():
    """``Z = {W}`` in ``X -> Y, X -> W`` is GAC-valid and back-door-invalid.

    ``W`` is a descendant of ``X``, which back-door bans outright, but it is not
    on or below the causal route ``X -> Y``, so the GAC permits it -- and it
    genuinely is a valid adjustment set, since summing ``P(y | x, w) P(w)`` over
    ``w`` telescopes to ``P(y | x) = P(y | do(x))``.
    """
    dag = MPDAG(["W", "X", "Y"], directed=[("X", "Y"), ("X", "W")])
    z = frozenset({"W"})

    assert is_gac_valid_dag(dag, "X", "Y", z)
    assert is_gac_valid_dag_by_paths(dag, "X", "Y", z)
    assert not is_valid_adjustment_set_dag(dag, "X", "Y", z)

    assert causal_nodes_dag(dag, "X", "Y") == {"Y"}
    assert forbidden_set_dag(dag, "X", "Y") == {"X", "Y"}
    assert "W" in dag.descendants("X")

    # The same contrast survives the lift: the DAG is its own only extension.
    assert is_gac_valid_mpdag(dag, "X", "Y", z)
    assert not is_valid_mpdag(dag, "X", "Y", z)
    assert why_invalid_gac(dag, "X", "Y", z) == ""


def test_a_descendant_on_the_causal_route_is_still_forbidden():
    """``X -> M -> Y`` with ``Z = {M}``: a mediator is banned by both criteria."""
    dag = MPDAG(["M", "X", "Y"], directed=[("X", "M"), ("M", "Y")])
    assert forbidden_set_dag(dag, "X", "Y") == {"M", "X", "Y"}
    assert not is_gac_valid_dag(dag, "X", "Y", frozenset({"M"}))
    assert not is_valid_adjustment_set_dag(dag, "X", "Y", frozenset({"M"}))


def test_a_descendant_of_a_mediator_is_forbidden_too():
    """``X -> M -> Y``, ``M -> D``: ``D`` hangs off the causal route, so it is banned."""
    dag = MPDAG(["D", "M", "X", "Y"], directed=[("X", "M"), ("M", "Y"), ("M", "D")])
    assert forbidden_set_dag(dag, "X", "Y") == {"D", "M", "X", "Y"}
    assert not is_gac_valid_dag(dag, "X", "Y", frozenset({"D"}))


# --------------------------------------------------------------------------------
# Layer 1: the reduction against the path enumeration, and against back-door
# --------------------------------------------------------------------------------


def test_layer1_reduction_equals_the_path_enumeration_on_every_4_node_dag():
    """The single d-separation test must equal the direct path-by-path check.

    :func:`is_gac_valid_dag` decides condition 2 by deleting the first edges of
    the proper causal paths and asking for d-separation, an equivalence that is
    *derived* in the module docstring under the hypothesis that condition 1
    holds. :func:`is_gac_valid_dag_by_paths` uses none of that. This is the test
    that turns the derivation into a measured claim.
    """
    checked = 0
    for dag, x, y, z in _dag_cases(DAGS4):
        assert is_gac_valid_dag(dag, x, y, z) == is_gac_valid_dag_by_paths(dag, x, y, z), (
            dag.edge_string(),
            x,
            y,
            sorted(z),
        )
        checked += 1
    assert checked > 20000


def test_backdoor_validity_implies_gac_validity_on_every_4_node_dag():
    """The GAC accepts a superset: back-door valid => GAC valid, never the reverse.

    Both halves are asserted here -- the implication on every case, and that the
    converse genuinely fails somewhere, so the superset is strict rather than an
    equality dressed up as one.
    """
    strict_witnesses = 0
    for dag, x, y, z in _dag_cases(DAGS4):
        gac = is_gac_valid_dag(dag, x, y, z)
        backdoor = is_valid_adjustment_set_dag(dag, x, y, z)
        assert not (backdoor and not gac), (dag.edge_string(), x, y, sorted(z))
        if gac and not backdoor:
            strict_witnesses += 1
    assert strict_witnesses > 0


def test_every_layer1_disagreement_is_the_expected_kind():
    """Disagreements must all be "``z`` holds a descendant of ``x`` off the causal route"."""
    for dag, x, y, z in _dag_cases(DAGS4):
        if is_gac_valid_dag(dag, x, y, z) == is_valid_adjustment_set_dag(dag, x, y, z):
            continue
        assert classify_disagreement(dag, x, y, z) == "z_has_descendant_of_x_off_the_causal_route"


def test_check1_on_a_bounded_slice_reports_no_unexplainable_disagreement():
    """The sweep function itself, on a slice small enough for the test suite."""
    record = run_check1(step5=97)
    assert record["n_cases"] > 70000
    assert record["n_disagree"] > 0
    assert record["n_unexplainable_disagreements"] == 0
    assert record["backdoor_valid_implies_gac_valid"]
    assert set(record["disagreement_classes"]) == {"z_has_descendant_of_x_off_the_causal_route"}


def test_the_degenerate_query_is_rejected():
    dag = MPDAG(["X", "Y"], directed=[("X", "Y")])
    assert not is_gac_valid_dag(dag, "X", "X", frozenset())
    assert why_invalid_gac(dag, "X", "X", frozenset()) == "degenerate_query"


# --------------------------------------------------------------------------------
# Layer 2: the semantic anchor
# --------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def scope_graphs():
    """Every CPDAG on 3 and 4 nodes with an undirected edge, plus their spaces."""
    return mpdag_scope((3, 4))


def test_layer2_equals_gac_validity_in_every_extension(scope_graphs):
    """The acceptance test: ``is_gac_valid_mpdag`` == ``all(is_gac_valid_dag(d, ...))``.

    Exhaustive over the same graph scope on which session 4's back-door
    criterion was verified, every ordered ``(x, y)`` and every candidate ``z``.
    """
    n_cases = 0
    for g, x, y, z in _mpdag_cases(scope_graphs):
        assert is_gac_valid_mpdag(g, x, y, z) == enumerated_gac_semantics(g, x, y, z), (
            g.edge_string(),
            x,
            y,
            sorted(z),
        )
        n_cases += 1
    assert n_cases > 70000


def test_the_hard_strata_are_actually_present(scope_graphs):
    """A zero-disagreement headline is only informative if the hard cases were in scope."""
    non_amenable = sum(
        1
        for g, x, y, z in _mpdag_cases(scope_graphs)
        if why_invalid_gac(g, x, y, z) == "not_amenable"
    )
    assert non_amenable > 1000


def test_reasons_are_from_the_documented_vocabulary_and_mean_invalid(scope_graphs):
    for g, x, y, z in _mpdag_cases(scope_graphs[::7]):
        reason = why_invalid_gac(g, x, y, z)
        assert reason == "" or reason in GAC_REASONS
        assert is_gac_valid_mpdag(g, x, y, z) == (reason == "")


def test_the_gac_forbidden_set_is_inside_the_backdoor_one(scope_graphs):
    """``forb(x, y, G) subset possde(x, G) u {x, y}``, and strictly so somewhere.

    Provable, and worth pinning: every anchor ``w`` of the GAC forbidden set is a
    child of ``x`` in ``G``, hence a child of ``x`` in *every* extension, so
    ``de(w, D) subset de(x, D)`` in each one and the inclusive possible-descendant
    sets nest. The strictness witness is what makes the two criteria differ.
    """
    strict = 0
    for g in scope_graphs:
        closed = criterion_module._closed(g)
        for x, y in itertools.permutations(sorted(g.nodes), 2):
            gac_forb = gac_forbidden_set(closed, x, y)
            bd_forb = backdoor_forbidden_set(closed, x, y)
            assert gac_forb <= bd_forb, (g.edge_string(), x, y)
            if gac_forb < bd_forb:
                strict += 1
    assert strict > 0


def test_backdoor_validity_implies_gac_validity_at_mpdag_level(scope_graphs):
    for g, x, y, z in _mpdag_cases(scope_graphs):
        if is_valid_mpdag(g, x, y, z):
            assert is_gac_valid_mpdag(g, x, y, z), (g.edge_string(), x, y, sorted(z))


def test_an_empty_extension_set_certifies_nothing():
    """The chordless undirected 4-cycle has no DAG extension; the convention is False.

    This is the documented incompleteness of the emptiness detection: the graph
    is acyclic and Meek-closed, so the criterion answers as though it had
    extensions and reports True where the enumeration reports False. It is not
    an MPDAG (a chain component of an essential graph is chordal), so it is
    outside the stated assumption -- pinned here so the limit stays visible
    rather than being rediscovered.
    """
    cycle = MPDAG(
        ["V0", "V1", "V2", "V3"],
        undirected=[("V0", "V1"), ("V1", "V2"), ("V2", "V3"), ("V0", "V3")],
    )
    assert extensions(cycle) == ()
    assert not enumerated_gac_semantics(cycle, "V0", "V1", frozenset())
    assert is_gac_valid_mpdag(cycle, "V0", "V1", frozenset())

    # A graph whose closure genuinely FAILs is caught, and matches the semantics.
    contradiction = MPDAG(
        ["V0", "V1", "V2"],
        directed=[("V0", "V1"), ("V1", "V2"), ("V2", "V0")],
    )
    assert why_invalid_gac(contradiction, "V0", "V1", frozenset()) == "no_extensions"


# --------------------------------------------------------------------------------
# Determinism and cost
# --------------------------------------------------------------------------------


def test_results_do_not_depend_on_the_caches(scope_graphs):
    graphs = scope_graphs[::11]
    first = [is_gac_valid_mpdag(g, x, y, z) for g, x, y, z in _mpdag_cases(graphs)]
    clear_cache()
    second = [is_gac_valid_mpdag(g, x, y, z) for g, x, y, z in _mpdag_cases(graphs)]
    assert first == second


def test_the_decision_never_calls_a_path_enumeration(monkeypatch):
    """A regression guard: the exponential walks must stay off the decision path."""

    def boom(*args: object, **kwargs: object) -> None:
        raise AssertionError("the decision reached a path-enumerating function")

    monkeypatch.setattr(criterion_module, "simple_paths", boom)
    monkeypatch.setattr(criterion_module, "possibly_causal_paths", boom)
    monkeypatch.setattr(paths_module, "simple_paths", boom)
    monkeypatch.setattr(paths_module, "possibly_causal_paths", boom)
    monkeypatch.setattr(paths_module, "unshielded_possibly_causal_paths", boom)
    monkeypatch.setattr(paths_module, "_unshielded_reach_by_paths", boom)
    clear_cache()
    for g, x, y, z in _mpdag_cases(mpdag_scope((3,))):
        why_invalid_gac(g, x, y, z)


def test_is_gac_valid_mpdag_meets_the_performance_bar_at_n_20():
    """Under 5 ms per cold-cache call on dense Erdos-Renyi instances at ``n = 20``.

    The bar exists because this predicate is on a hot path, and because session
    4's first criterion was exhaustively verified *and* exponential -- only
    correctness had been specified.
    """
    samples: list[float] = []
    for seed in (0, 1):
        for prob in (0.7, 0.85):
            for half in (False, True):
                g = dense_instance(20, seed, prob, half_knowledge=half)
                nodes = sorted(g.nodes)
                queries = [
                    (x, nodes[-1], frozenset(nodes[1:5])) for x in nodes[:4] if x != nodes[-1]
                ]
                samples.extend(time_calls(g, queries))
    worst = max(samples)
    assert worst < 5.0, f"worst cold-cache call took {worst:.2f} ms at n=20"


def test_the_whole_module_is_fast_enough_to_stay_in_ci():
    """A sanity bound on the fixtures themselves, not on the criterion."""
    started = time.perf_counter()
    mpdag_scope((3,))
    assert time.perf_counter() - started < 30.0
