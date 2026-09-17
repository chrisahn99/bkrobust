"""Guards for the polynomial committed set.

Three properties, each the reason a caller may trust the verdict it gets back:
the set is never returned invalid on the graph it was computed on; a query the
analyst's graph says is reversed is refused with its own verdict rather than
scored; and on a fully oriented graph the closed form is the DAG-level optimal
set exactly.
"""

from __future__ import annotations

from bkrobust.demo.evaluate import optimal_adjustment_set_dag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag
from bkrobust.mpdag_criterion.optimal import committed_adjustment_set


def _chain_with_confounder() -> MPDAG:
    """``U -> X -> M -> Y`` with ``U -> Y``: the optimal set for ``(X, Y)`` is ``{U}``."""
    nodes = ("U", "X", "M", "Y")
    return MPDAG(nodes, directed=(("U", "X"), ("X", "M"), ("M", "Y"), ("U", "Y")))


def test_fully_oriented_graph_matches_dag_optimal_set() -> None:
    dag = _chain_with_confounder()
    cs = committed_adjustment_set(dag, "X", "Y")
    assert cs.verdict == "closed_form"
    assert cs.z == frozenset(optimal_adjustment_set_dag(dag, "X", "Y"))


def test_returned_set_is_valid_on_the_graph_it_was_computed_on() -> None:
    dag = _chain_with_confounder()
    cpdag = dag_to_cpdag(dag)
    for x, y in (("X", "Y"), ("X", "M"), ("U", "Y"), ("M", "Y")):
        g0 = apply_orientations(cpdag, [("U", "X")])
        assert g0 is not None
        cs = committed_adjustment_set(g0, x, y)
        if cs.z is not None:
            assert is_gac_valid_mpdag(g0, x, y, cs.z), (x, y, cs)


def test_reversed_query_is_a_verdict_not_a_certificate() -> None:
    dag = _chain_with_confounder()
    cs = committed_adjustment_set(dag, "Y", "X")
    assert cs.verdict == "y_not_possible_descendant"
    assert cs.z is None


def test_no_edge_between_pair_is_not_a_crash() -> None:
    g = MPDAG(("A", "B", "C"), directed=(("A", "B"),))
    cs = committed_adjustment_set(g, "A", "C")
    assert cs.verdict == "y_not_possible_descendant"
