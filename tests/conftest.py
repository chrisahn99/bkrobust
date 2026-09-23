"""Shared fixtures: small canonical graphs, hardcoded as edge lists.

These are the graphs every test reasons about. They are small enough to verify
by hand, which is the point -- a Meek implementation is checked against examples
whose answers are known independently of the code, not against its own output.

Each fixture returns a plain dict of edge lists rather than an
:class:`~bkrobust.graphs.mpdag.MPDAG`, so the fixtures work before ``MPDAG``
does and a constructor bug produces one failure rather than a collection error.
Tests build graphs from them via the ``mpdag`` factory fixture.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import pytest


@pytest.fixture
def chain() -> dict[str, Any]:
    """A -> B -> C.

    CPDAG: fully undirected (A - B - C). The three-node chain, the fork and the
    collider share a skeleton; only the collider has a v-structure, so the chain
    and the fork are Markov equivalent and the collider is not. Adjustment for
    the effect of A on C needs the empty set: B is a mediator and adjusting for
    it would block the very path being measured.
    """
    return {
        "nodes": ["A", "B", "C"],
        "directed": [("A", "B"), ("B", "C")],
        "undirected": [],
        "cpdag_directed": [],
        "cpdag_undirected": [("A", "B"), ("B", "C")],
        "treatment": "A",
        "outcome": "C",
        "optimal_adjustment_set": set(),
    }


@pytest.fixture
def fork() -> dict[str, Any]:
    """B -> A, B -> C: a common cause.

    CPDAG: fully undirected, identical to the chain's -- these two are Markov
    equivalent, which is why background knowledge has anything to do. The total
    effect of A on C is zero, and adjusting for B is what reveals that.
    """
    return {
        "nodes": ["A", "B", "C"],
        "directed": [("B", "A"), ("B", "C")],
        "undirected": [],
        "cpdag_directed": [],
        "cpdag_undirected": [("A", "B"), ("B", "C")],
        "treatment": "A",
        "outcome": "C",
        "optimal_adjustment_set": {"B"},
    }


@pytest.fixture
def collider() -> dict[str, Any]:
    """A -> B <- C: a v-structure.

    CPDAG: fully directed -- the v-structure is identifiable from data, so both
    edges are oriented and the CPDAG equals the DAG. Nothing here is left for
    background knowledge to orient, which makes this the fixture that pins down
    what *cannot* be perturbed: asserting ``B -> A`` is inconsistent, not
    consistent-but-false, and the perturbation samplers must never propose it.
    """
    return {
        "nodes": ["A", "B", "C"],
        "directed": [("A", "B"), ("C", "B")],
        "undirected": [],
        "cpdag_directed": [("A", "B"), ("C", "B")],
        "cpdag_undirected": [],
        "treatment": "A",
        "outcome": "B",
        "optimal_adjustment_set": {"C"},
    }


@pytest.fixture
def diamond() -> dict[str, Any]:
    """A -> B -> D, A -> C -> D: two parallel directed paths.

    The smallest graph where Meek's R2 matters: orienting ``A - D`` requires
    noticing the directed path through B or C. The effect of A on D runs through
    both mediators, so neither may be adjusted for and the optimal set is empty.
    """
    return {
        "nodes": ["A", "B", "C", "D"],
        "directed": [("A", "B"), ("A", "C"), ("B", "D"), ("C", "D")],
        "undirected": [],
        "cpdag_directed": [],
        "cpdag_undirected": [("A", "B"), ("A", "C"), ("B", "D"), ("C", "D")],
        "treatment": "A",
        "outcome": "D",
        "optimal_adjustment_set": set(),
    }


@pytest.fixture
def nontrivial_optimal() -> dict[str, Any]:
    """A graph where O* differs from the canonical back-door set.

    ``Z1 -> X``, ``Z1 -> Y``, ``Z2 -> Y``, ``X -> M -> Y``, ``X -> Y``.

    The canonical back-door set is ``{Z1}``: it blocks the confounding path and
    is valid. ``O*`` is ``{Z1, Z2}``. ``Z2`` is not a confounder -- it does not
    touch the treatment -- so adding it changes nothing about validity, but it
    explains outcome variance and so shrinks the residual variance of the
    adjusted estimator. That is the whole content of "optimal": it is a variance
    claim, not an identification claim, and this fixture is where a test can
    tell the two apart.

    ``M`` is a mediator and therefore forbidden: adjusting for it would block
    part of the effect being estimated.
    """
    return {
        "nodes": ["Z1", "Z2", "X", "M", "Y"],
        "directed": [
            ("Z1", "X"),
            ("Z1", "Y"),
            ("Z2", "Y"),
            ("X", "M"),
            ("M", "Y"),
            ("X", "Y"),
        ],
        "undirected": [],
        # The CPDAG of this graph is NOT hardcoded: working it out by hand is
        # error-prone, and a wrong hardcoded value would silently weaken every
        # test that used it. Derive it with `bkrobust.data.synthetic.dag_to_cpdag`
        # and assert against that instead.
        "cpdag_directed": None,
        "cpdag_undirected": None,
        "treatment": "X",
        "outcome": "Y",
        "optimal_adjustment_set": {"Z1", "Z2"},
        "canonical_adjustment_set": {"Z1"},
        "forbidden_set": {"X", "M", "Y"},
    }


@pytest.fixture
def canonical_graphs(
    chain: dict[str, Any],
    fork: dict[str, Any],
    collider: dict[str, Any],
    diamond: dict[str, Any],
    nontrivial_optimal: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    """All canonical graphs by name, for parameterised sweeps over them."""
    return {
        "chain": chain,
        "fork": fork,
        "collider": collider,
        "diamond": diamond,
        "nontrivial_optimal": nontrivial_optimal,
    }


@pytest.fixture
def mpdag() -> Callable[[dict[str, Any]], Any]:
    """Factory turning a fixture dict into an :class:`~bkrobust.graphs.mpdag.MPDAG`.

    Returns:
        A callable taking a fixture dict and an optional ``as_cpdag`` flag,
        which selects the ``cpdag_*`` edge lists instead of the DAG ones.
    """

    def _build(spec: dict[str, Any], as_cpdag: bool = False) -> Any:
        from bkrobust.graphs.mpdag import MPDAG

        if as_cpdag:
            if spec["cpdag_directed"] is None:
                raise ValueError(
                    "this fixture has no hardcoded CPDAG; derive it with "
                    "bkrobust.data.synthetic.dag_to_cpdag"
                )
            return MPDAG(
                nodes=spec["nodes"],
                directed=spec["cpdag_directed"],
                undirected=spec["cpdag_undirected"],
            )
        return MPDAG(
            nodes=spec["nodes"],
            directed=spec["directed"],
            undirected=spec["undirected"],
        )

    return _build


@pytest.fixture
def rng() -> Any:
    """A seeded NumPy generator. Tests must never touch global RNG state."""
    import numpy as np

    return np.random.default_rng(20260919)
