"""Component and differential tests for the MPDAG adjustment criterion."""

from __future__ import annotations

import itertools

import pytest

from bkrobust.core.oracle import extensions, is_valid
from bkrobust.demo.graph import MPDAG
from bkrobust.mpdag_criterion import (
    REASONS,
    backdoor_forbidden_set,
    causal_nodes,
    clear_cache,
    definite_status_non_causal_paths,
    forbidden_set,
    is_amenable,
    is_blocked,
    is_collider,
    is_definite_non_collider,
    is_definite_status_node,
    is_definite_status_path,
    is_non_causal,
    is_possibly_causal,
    is_unshielded,
    is_valid_mpdag,
    possible_descendants,
    possibly_causal_paths,
    simple_paths,
    unshielded_possibly_causal_paths,
    why_invalid,
)
from bkrobust.mpdag_criterion.sweep import candidate_sets, run_sweep
from bkrobust.search.conjecture_study import all_cpdags
from bkrobust.search.space_fixed import build_space_correct

# --------------------------------------------------------------------------------
# Fixtures: small hand-checked graphs
# --------------------------------------------------------------------------------

#: The witness that "possibly causal" is not enough on an MPDAG carrying
#: background knowledge: <V1, V2, V0> is possibly causal but shielded, and no
#: extension makes V0 a descendant of V1 (that would close V0 -> V1 -> V2 -> V0).
KNOWN_TRIANGLE = MPDAG(
    ["V0", "V1", "V2"],
    directed=[("V0", "V1")],
    undirected=[("V0", "V2"), ("V1", "V2")],
)

#: A chain CPDAG: every edge undirected, no v-structure.
CHAIN = MPDAG(["V0", "V1", "V2"], undirected=[("V0", "V1"), ("V1", "V2")])

#: A v-structure DAG (also its own CPDAG).
COLLIDER = MPDAG(["V0", "V1", "V2"], directed=[("V0", "V1"), ("V2", "V1")])


def graph_scope(sizes: tuple[int, ...], step: int = 1) -> list[MPDAG]:
    """CPDAGs with an undirected edge plus their corrected-space elements."""
    seen: dict[tuple[tuple[str, ...], str], MPDAG] = {}
    for n in sizes:
        cpdags = [c for c in all_cpdags(n) if c.undirected_edges]
        for cpdag in cpdags[::step]:
            seen[(cpdag.nodes, cpdag.edge_string())] = cpdag
            for element in build_space_correct(cpdag).elements:
                seen[(element.nodes, element.edge_string())] = element
    return [seen[k] for k in sorted(seen)]


def differential(graphs: list[MPDAG], max_z: int = 3) -> tuple[int, int, list[str]]:
    """Compare the criterion to the oracle; return (cases, disagreements, witnesses)."""
    cases = 0
    bad: list[str] = []
    for g in graphs:
        for x, y in itertools.permutations(sorted(g.nodes), 2):
            for z in candidate_sets(g.nodes, x, y, max_z):
                cases += 1
                if is_valid_mpdag(g, x, y, z) != is_valid(z, g, x, y):
                    bad.append(f"{g.edge_string()} | x={x} y={y} z={sorted(z)}")
    return cases, len(bad), bad[:5]


# --------------------------------------------------------------------------------
# Possibly causal paths
# --------------------------------------------------------------------------------


def test_possibly_causal_ignores_forward_and_undirected_edges():
    assert is_possibly_causal(CHAIN, ("V0", "V1", "V2"))
    assert not is_non_causal(CHAIN, ("V0", "V1", "V2"))


def test_possibly_causal_rejects_a_backward_edge():
    assert not is_possibly_causal(COLLIDER, ("V0", "V1", "V2"))
    assert is_non_causal(COLLIDER, ("V0", "V1", "V2"))
    # ...but the same path read from the other end starts forward at V2? No:
    # V2 -> V1 <- V0 is backward at the far end too.
    assert is_non_causal(COLLIDER, ("V2", "V1", "V0"))


def test_simple_paths_are_deterministic_and_complete():
    g = MPDAG(["A", "B", "C"], undirected=[("A", "B"), ("B", "C"), ("A", "C")])
    assert simple_paths(g, "A", "C") == [("A", "B", "C"), ("A", "C")]
    assert simple_paths(g, "A", "A") == []


def test_possibly_causal_paths_are_a_subset_of_simple_paths():
    for g in graph_scope((3,)):
        for x, y in itertools.permutations(sorted(g.nodes), 2):
            allp = set(simple_paths(g, x, y))
            pcp = set(possibly_causal_paths(g, x, y))
            assert pcp <= allp
            assert pcp == {p for p in allp if is_possibly_causal(g, p)}


def test_unshielded_paths_are_a_subset_of_possibly_causal_paths():
    for g in graph_scope((3,)):
        for x, y in itertools.permutations(sorted(g.nodes), 2):
            pcp = set(possibly_causal_paths(g, x, y))
            usp = set(unshielded_possibly_causal_paths(g, x, y))
            assert usp == {p for p in pcp if is_unshielded(g, p)}


# --------------------------------------------------------------------------------
# Possible descendants
# --------------------------------------------------------------------------------


def test_possible_descendants_is_inclusive():
    assert "V0" in possible_descendants(CHAIN, "V0")
    assert possible_descendants(CHAIN, "V0") == {"V0", "V1", "V2"}


def test_possible_descendants_excludes_the_shielded_witness():
    # <V1, V2, V0> is possibly causal but shielded, and unrealisable: every
    # extension keeps V0 -> V1, so V0 is never a descendant of V1.
    assert possible_descendants(KNOWN_TRIANGLE, "V1") == {"V1", "V2"}
    assert possible_descendants(KNOWN_TRIANGLE, "V0") == {"V0", "V1", "V2"}


def test_possible_descendants_equals_the_union_over_extensions():
    for g in graph_scope((3,)) + graph_scope((4,), step=11):
        for x in g.nodes:
            union = {x}
            for d in extensions(g):
                union |= d.descendants(x)
            assert possible_descendants(g, x) == union, g.edge_string()


def test_possible_descendants_of_a_set_is_the_union_over_members():
    g = KNOWN_TRIANGLE
    assert possible_descendants(g, ["V1", "V2"]) == possible_descendants(
        g, "V1"
    ) | possible_descendants(g, "V2")
    assert possible_descendants(g, []) == set()


# --------------------------------------------------------------------------------
# Definite status
# --------------------------------------------------------------------------------


def test_collider_is_recognised():
    path = ("V0", "V1", "V2")
    assert is_collider(COLLIDER, path, 1)
    assert not is_definite_non_collider(COLLIDER, path, 1)
    assert is_definite_status_path(COLLIDER, path)


def test_endpoints_are_always_definite_status():
    path = ("V0", "V1", "V2")
    assert is_definite_status_node(COLLIDER, path, 0)
    assert is_definite_status_node(COLLIDER, path, 2)
    assert not is_collider(COLLIDER, path, 0)
    assert not is_definite_non_collider(COLLIDER, path, 0)


def test_definite_non_collider_via_an_edge_out_of_the_node():
    g = MPDAG(["V0", "V1", "V2"], directed=[("V0", "V1"), ("V1", "V2")])
    assert is_definite_non_collider(g, ("V0", "V1", "V2"), 1)
    # and read backwards: V1 -> V0 is an edge out of V1 on the path
    assert is_definite_non_collider(g, ("V2", "V1", "V0"), 1)


def test_definite_non_collider_via_the_unshielded_undirected_pattern():
    assert is_definite_non_collider(CHAIN, ("V0", "V1", "V2"), 1)


def test_shielded_undirected_pattern_is_indefinite_status():
    # V0 - V1 - V2 with V0, V2 adjacent: an extension may orient both edges into
    # V1 without creating an *unshielded* collider, so V1's status is undecided.
    g = MPDAG(["V0", "V1", "V2"], undirected=[("V0", "V1"), ("V1", "V2"), ("V0", "V2")])
    assert not is_definite_non_collider(g, ("V0", "V1", "V2"), 1)
    assert not is_collider(g, ("V0", "V1", "V2"), 1)
    assert not is_definite_status_node(g, ("V0", "V1", "V2"), 1)
    assert not is_definite_status_path(g, ("V0", "V1", "V2"))


def test_mixed_shielded_pattern_is_indefinite_status():
    # V0 -> V1 - V2 with V0, V2 adjacent. Meek's rule 1 does not fire (the
    # triple is shielded), and V1 is a collider in some extensions only.
    assert not is_definite_status_node(KNOWN_TRIANGLE, ("V0", "V1", "V2"), 1)


# --------------------------------------------------------------------------------
# Blocking
# --------------------------------------------------------------------------------


def test_collider_blocks_when_nothing_below_it_is_conditioned_on():
    assert is_blocked(COLLIDER, ("V0", "V1", "V2"), frozenset())
    assert not is_blocked(COLLIDER, ("V0", "V1", "V2"), frozenset({"V1"}))


def test_definite_non_collider_blocks_exactly_when_conditioned_on():
    assert not is_blocked(CHAIN, ("V0", "V1", "V2"), frozenset())
    assert is_blocked(CHAIN, ("V0", "V1", "V2"), frozenset({"V1"}))


def test_collider_is_opened_by_a_possible_descendant():
    # V0 -> V1 <- V2 with V1 - V3 undirected, shielded so Meek's rule 1 cannot
    # orient it. Conditioning on V3 opens the collider in the extension that
    # orients V1 -> V3.
    g = MPDAG(
        ["V0", "V1", "V2", "V3"],
        directed=[("V0", "V1"), ("V2", "V1")],
        undirected=[("V1", "V3"), ("V0", "V3"), ("V2", "V3")],
    )
    assert "V3" in possible_descendants(g, "V1")
    assert not is_blocked(g, ("V0", "V1", "V2"), frozenset({"V3"}))


# --------------------------------------------------------------------------------
# Amenability
# --------------------------------------------------------------------------------


def test_a_bare_undirected_edge_is_not_amenable():
    g = MPDAG(["V0", "V1"], undirected=[("V0", "V1")])
    assert not is_amenable(g, "V0", "V1")
    assert why_invalid(g, "V0", "V1", frozenset()) == "not_amenable"
    assert not is_valid(frozenset(), g, "V0", "V1")


def test_a_directed_edge_is_amenable():
    g = MPDAG(["V0", "V1"], directed=[("V0", "V1")])
    assert is_amenable(g, "V0", "V1")
    assert is_valid_mpdag(g, "V0", "V1", frozenset())


def test_amenability_counts_shielded_realisable_paths():
    # (V0, V1): the extension V0 -> V2 -> V1 is a causal path starting on the
    # undirected edge V0 - V2, so the graph is NOT amenable -- even though
    # <V0, V2, V1> is shielded. The unshielded shortcut would wrongly say yes.
    assert not is_amenable(KNOWN_TRIANGLE, "V0", "V1")


def test_amenability_ignores_unrealisable_shielded_paths():
    # (V1, V0): <V1, V2, V0> is possibly causal and starts undirected, but no
    # extension realises it, so the graph IS amenable. The naive
    # "any possibly causal path" rule would wrongly say no.
    assert is_amenable(KNOWN_TRIANGLE, "V1", "V0")


def test_amenability_matches_its_semantics_over_extensions():
    for g in graph_scope((3,)) + graph_scope((4,), step=11):
        for x, y in itertools.permutations(sorted(g.nodes), 2):
            semantic = True
            for d in extensions(g):
                for v in sorted(d.children(x)):
                    if g.is_undirected_edge(x, v) and (v == y or y in d.descendants(v)):
                        semantic = False
            assert is_amenable(g, x, y) == semantic, f"{g.edge_string()} {x} {y}"


# --------------------------------------------------------------------------------
# Forbidden sets
# --------------------------------------------------------------------------------


def test_causal_nodes_exclude_the_treatment_and_include_the_outcome():
    assert causal_nodes(CHAIN, "V0", "V2") == {"V1", "V2"}
    assert causal_nodes(COLLIDER, "V0", "V2") == set()


def test_backdoor_forbidden_set_is_larger_than_textbook_forb():
    # V0 -> V1 and V0 -> V2, with V2 off every causal path to V1. The textbook
    # forbidden set does not bar V2 -- adjusting for it really is valid -- but
    # Pearl's back-door criterion does, and the oracle implements back-door.
    g = MPDAG(["V0", "V1", "V2"], directed=[("V0", "V1"), ("V0", "V2")])
    assert forbidden_set(g, "V0", "V1") == {"V0", "V1"}
    assert backdoor_forbidden_set(g, "V0", "V1") == {"V0", "V1", "V2"}
    assert why_invalid(g, "V0", "V1", frozenset({"V2"})) == "z_hits_forbidden"
    assert not is_valid(frozenset({"V2"}), g, "V0", "V1")


def test_forbidden_set_matches_its_stated_definition():
    # forb = possde(cn) u {x}, with cn read off the possibly causal paths.
    for g in graph_scope((3,)) + graph_scope((4,), step=11):
        for x, y in itertools.permutations(sorted(g.nodes), 2):
            cn = {v for p in possibly_causal_paths(g, x, y) for v in p[1:]}
            expected = (possible_descendants(g, cn) if cn else set()) | {x}
            assert forbidden_set(g, x, y) == expected


def test_textbook_forb_is_not_contained_in_the_backdoor_forbidden_set():
    # Recorded deliberately: the DAG-level containment forb <= de(x) u {x} does
    # NOT lift to an MPDAG carrying background knowledge, because cn counts
    # shielded possibly causal paths that no extension realises. Here no
    # extension has any causal path from V1 to V0, yet cn(V1, V0) = {V3, V0}.
    g = MPDAG(
        ["V0", "V1", "V2", "V3"],
        directed=[("V0", "V1"), ("V0", "V2")],
        undirected=[("V0", "V3"), ("V1", "V3")],
    )
    assert causal_nodes(g, "V1", "V0") == {"V3", "V0"}
    assert forbidden_set(g, "V1", "V0") == {"V0", "V1", "V2", "V3"}
    assert backdoor_forbidden_set(g, "V1", "V0") == {"V0", "V1", "V3"}
    assert not forbidden_set(g, "V1", "V0") <= backdoor_forbidden_set(g, "V1", "V0")
    assert all("V0" not in d.descendants("V1") for d in extensions(g))


# --------------------------------------------------------------------------------
# Reasons
# --------------------------------------------------------------------------------


def test_open_noncausal_path_is_reported():
    # V0 <- V1 -> V2 : an unblocked confounding path.
    g = MPDAG(["V0", "V1", "V2"], directed=[("V1", "V0"), ("V1", "V2")])
    assert definite_status_non_causal_paths(g, "V0", "V2") == [("V0", "V1", "V2")]
    assert why_invalid(g, "V0", "V2", frozenset()) == "open_noncausal_path"
    assert why_invalid(g, "V0", "V2", frozenset({"V1"})) == ""
    assert is_valid_mpdag(g, "V0", "V2", frozenset({"V1"}))


def test_no_extensions_is_reported_for_a_directed_cycle():
    g = MPDAG(["V0", "V1", "V2"], directed=[("V0", "V1"), ("V1", "V2"), ("V2", "V0")])
    assert extensions(g) == ()
    assert why_invalid(g, "V0", "V1", frozenset()) == "no_extensions"
    assert not is_valid_mpdag(g, "V0", "V1", frozenset())


def test_no_extensions_is_reported_for_a_meek_conflict():
    g = MPDAG(
        ["V0", "V1", "V2", "V3"],
        directed=[("V0", "V1"), ("V3", "V2")],
        undirected=[("V1", "V2")],
    )
    assert extensions(g) == ()
    assert why_invalid(g, "V0", "V1", frozenset()) == "no_extensions"


def test_degenerate_query_when_treatment_equals_outcome():
    assert why_invalid(CHAIN, "V0", "V0", frozenset()) == "degenerate_query"
    assert not is_valid_mpdag(CHAIN, "V0", "V0", frozenset())
    assert not is_valid(frozenset(), CHAIN, "V0", "V0")


def test_every_reason_returned_is_declared():
    seen = set()
    for g in graph_scope((3,)):
        for x, y in itertools.permutations(sorted(g.nodes), 2):
            for z in candidate_sets(g.nodes, x, y, 3):
                reason = why_invalid(g, x, y, z)
                assert reason == "" or reason in REASONS
                seen.add(reason)
    assert {"", "not_amenable", "z_hits_forbidden", "open_noncausal_path"} <= seen


def test_reason_is_empty_exactly_when_valid():
    for g in graph_scope((3,)):
        for x, y in itertools.permutations(sorted(g.nodes), 2):
            for z in candidate_sets(g.nodes, x, y, 3):
                assert is_valid_mpdag(g, x, y, z) == (why_invalid(g, x, y, z) == "")


# --------------------------------------------------------------------------------
# Differential agreement with the enumerating oracle
# --------------------------------------------------------------------------------


def test_agrees_with_oracle_exhaustively_at_n3():
    graphs = graph_scope((3,))
    cases, disagreements, witnesses = differential(graphs)
    assert cases > 0
    assert disagreements == 0, witnesses


def test_agrees_with_oracle_on_a_sampled_n4_scope():
    graphs = graph_scope((4,), step=3)
    cases, disagreements, witnesses = differential(graphs)
    assert cases > 0
    assert disagreements == 0, witnesses


def test_non_amenable_graphs_are_present_and_agree_at_n4():
    graphs = graph_scope((4,), step=3)
    seen = agreed = 0
    for g in graphs:
        for x, y in itertools.permutations(sorted(g.nodes), 2):
            if is_amenable(g, x, y):
                continue
            for z in candidate_sets(g.nodes, x, y, 3):
                seen += 1
                # Non-amenable means no set is valid; the oracle must concur.
                assert not is_valid_mpdag(g, x, y, z)
                if not is_valid(z, g, x, y):
                    agreed += 1
    assert seen > 0, "no non-amenable cases in scope -- the stratum is untested"
    assert agreed == seen


@pytest.mark.slow
def test_full_sweep_reports_zero_disagreements():
    record = run_sweep()
    assert record["n_disagree"] == 0
    assert record["n_agree"] == record["n_cases"]
    assert record["non_amenable_agree"] == record["non_amenable_cases"]


# --------------------------------------------------------------------------------
# Determinism
# --------------------------------------------------------------------------------


def test_results_do_not_depend_on_the_closure_cache():
    graphs = graph_scope((3,))
    first = [
        is_valid_mpdag(g, x, y, z)
        for g in graphs
        for x, y in itertools.permutations(sorted(g.nodes), 2)
        for z in candidate_sets(g.nodes, x, y, 3)
    ]
    clear_cache()
    second = [
        is_valid_mpdag(g, x, y, z)
        for g in graphs
        for x, y in itertools.permutations(sorted(g.nodes), 2)
        for z in candidate_sets(g.nodes, x, y, 3)
    ]
    assert first == second


def test_path_enumeration_order_is_stable():
    g = graph_scope((4,), step=11)[0]
    nodes = sorted(g.nodes)
    assert simple_paths(g, nodes[0], nodes[-1]) == simple_paths(g, nodes[0], nodes[-1])
