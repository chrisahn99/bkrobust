"""The generalised adjustment criterion (GAC), at the DAG and the MPDAG level.

Why this package exists alongside the back-door code it does not touch
----------------------------------------------------------------------
The repository's oracle
(:func:`bkrobust.demo.evaluate.is_valid_adjustment_set_dag`, lifted by
:func:`bkrobust.core.oracle.is_valid`) implements **Pearl's back-door
criterion**, which is sufficient but not necessary for adjustment. It bars every
descendant of the treatment; what actually has to be barred is only what lies on
or below the causal route from the treatment to the outcome. In ``X -> Y``,
``X -> W`` the set ``Z = {W}`` is a perfectly valid adjustment set that back-door
rejects.

The **generalised adjustment criterion** is the sound *and complete* condition
(Shpitser, VanderWeele & Robins 2010; Perkovic, Textor, Kalisch & Maathuis
2018). This package implements it in two layers, leaving every back-door
implementation in the repository untouched so that both remain available as
independent oracles:

* :mod:`bkrobust.gac.dag_level` -- ``is_gac_valid_dag``, the criterion on a
  fully oriented DAG, plus a second path-enumerating implementation kept as its
  differential reference.
* :mod:`bkrobust.gac.mpdag_level` -- ``is_gac_valid_mpdag``, the genuine lift to
  a partially oriented graph: amenability, the GAC forbidden set lifted to
  MPDAGs, and blocking of proper definite-status non-causal paths, all decided
  in **polynomial** time on session 4's state-reachability machinery
  (:mod:`bkrobust.mpdag_criterion`).

The semantics of the MPDAG layer, and its acceptance test, is the enumeration

.. code-block:: python

    all(is_gac_valid_dag(d, x, y, z) for d in extensions(g))

with the repository's convention that an empty ``[g]`` certifies nothing. The
agreement numbers, their scopes and the timings are produced by
:mod:`bkrobust.gac.sweep` and written to
``results/axisa2/gac_agreement.json``.

Entry points::

    from bkrobust.gac import is_gac_valid_dag, is_gac_valid_mpdag, why_invalid_gac

    is_gac_valid_dag(dag, "V0", "V2", frozenset({"V1"}))  # -> bool
    is_gac_valid_mpdag(g, "V0", "V2", frozenset({"V1"}))  # -> bool
    why_invalid_gac(g, "V0", "V2", frozenset({"V1"}))  # -> "" or a reason
"""

from bkrobust.gac.dag_level import (
    causal_nodes_dag,
    forbidden_set_dag,
    is_gac_valid_dag,
    is_gac_valid_dag_by_paths,
)
from bkrobust.gac.mpdag_level import (
    GAC_REASONS,
    clear_cache,
    gac_causal_children,
    gac_forbidden_set,
    is_gac_valid_mpdag,
    why_invalid_gac,
)

__all__ = [
    "GAC_REASONS",
    "causal_nodes_dag",
    "clear_cache",
    "forbidden_set_dag",
    "gac_causal_children",
    "gac_forbidden_set",
    "is_gac_valid_dag",
    "is_gac_valid_dag_by_paths",
    "is_gac_valid_mpdag",
    "why_invalid_gac",
]
