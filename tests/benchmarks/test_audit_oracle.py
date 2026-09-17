"""The audit oracle is the generalised criterion, and differs from the back-door one.

A descendant of the treatment that sits on no causal path to the outcome may be
adjusted for. The back-door criterion forbids it; the generalised criterion,
which is what the certificate decides on the analyst's graph, does not. The
audit must use the latter, or it scores the instrument against a different
definition of validity than the instrument used.
"""

from __future__ import annotations

from bkrobust.benchmarks.audit import is_valid_adjustment_set_gac_dag
from bkrobust.demo.evaluate import is_valid_adjustment_set_dag
from bkrobust.demo.graph import MPDAG
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag


def _fork() -> MPDAG:
    """``U -> X -> C -> Y`` with ``U -> Y`` and a dangling child ``X -> A``."""
    nodes = ("U", "X", "C", "Y", "A")
    return MPDAG(nodes, directed=(("U", "X"), ("U", "Y"), ("X", "C"), ("C", "Y"), ("X", "A")))


def test_descendant_off_the_causal_path_is_allowed_by_the_criterion() -> None:
    dag = _fork()
    assert is_valid_adjustment_set_gac_dag(dag, "X", "Y", {"U", "A"})
    assert is_gac_valid_mpdag(dag, "X", "Y", frozenset({"U", "A"}))
    assert not is_valid_adjustment_set_dag(dag, "X", "Y", {"U", "A"})


def test_the_three_oracles_agree_off_that_case() -> None:
    dag = _fork()
    for z, expected in (({"U"}, True), (set(), False), ({"C"}, False), ({"U", "C"}, False)):
        assert is_valid_adjustment_set_gac_dag(dag, "X", "Y", z) is expected, z
        assert is_gac_valid_mpdag(dag, "X", "Y", frozenset(z)) is expected, z
        assert is_valid_adjustment_set_dag(dag, "X", "Y", z) is expected, z


def test_endpoints_are_never_an_adjustment_set() -> None:
    dag = _fork()
    assert not is_valid_adjustment_set_gac_dag(dag, "X", "Y", {"X"})
    assert not is_valid_adjustment_set_gac_dag(dag, "X", "Y", {"Y"})
