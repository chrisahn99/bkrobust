"""Adjustment-set validity decided on an MPDAG directly, without enumerating ``[G]``.

The rest of this repository answers "is ``Z`` a valid adjustment set for
``(X, Y)`` in the MPDAG ``G``?" by enumerating every DAG extension of ``G`` and
checking each one, which costs ``2**k`` back-door tests in the number ``k`` of
undirected edges and dominates every sweep. This package answers the same
question graphically, from ``G`` alone, via the generalised adjustment criterion
of Perkovic, Textor, Kalisch and Maathuis -- adapted to the *back-door*
semantics that :func:`bkrobust.core.oracle.is_valid` actually implements. See
:mod:`bkrobust.mpdag_criterion.criterion` for the derivation of that adaptation
and for the empty-``[G]`` policy.

Layout:

* :mod:`bkrobust.mpdag_criterion.paths` -- possibly causal paths, possible
  descendants, definite status, and path blocking. Each notion is a separate,
  separately testable function.
* :mod:`bkrobust.mpdag_criterion.criterion` -- amenability, the forbidden sets,
  and the public predicates :func:`is_valid_mpdag` and :func:`why_invalid`.

Entry points::

    from bkrobust.mpdag_criterion import is_valid_mpdag, why_invalid

    is_valid_mpdag(g, "V0", "V2", frozenset({"V1"}))  # -> bool
    why_invalid(g, "V0", "V2", frozenset({"V1"}))  # -> "" or a reason string
"""

from bkrobust.mpdag_criterion.criterion import (
    REASONS,
    backdoor_forbidden_set,
    causal_nodes,
    clear_cache,
    definite_status_non_causal_paths,
    forbidden_set,
    is_amenable,
    is_valid_mpdag,
    open_non_causal_path,
    why_invalid,
)
from bkrobust.mpdag_criterion.paths import (
    Path,
    clear_index_cache,
    has_open_definite_status_non_causal_path,
    is_blocked,
    is_collider,
    is_definite_non_collider,
    is_definite_status_node,
    is_definite_status_path,
    is_non_causal,
    is_possibly_causal,
    is_unshielded,
    open_definite_status_non_causal_path,
    possible_descendants,
    possibly_causal_paths,
    simple_paths,
    unshielded_possibly_causal_paths,
    unshielded_reachable,
)

__all__ = [
    "REASONS",
    "Path",
    "backdoor_forbidden_set",
    "causal_nodes",
    "clear_cache",
    "clear_index_cache",
    "definite_status_non_causal_paths",
    "forbidden_set",
    "has_open_definite_status_non_causal_path",
    "is_amenable",
    "is_blocked",
    "is_collider",
    "is_definite_non_collider",
    "is_definite_status_node",
    "is_definite_status_path",
    "is_non_causal",
    "is_possibly_causal",
    "is_unshielded",
    "is_valid_mpdag",
    "open_definite_status_non_causal_path",
    "open_non_causal_path",
    "possible_descendants",
    "possibly_causal_paths",
    "simple_paths",
    "unshielded_possibly_causal_paths",
    "unshielded_reachable",
    "why_invalid",
]
