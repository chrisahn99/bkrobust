"""Valid adjustment sets, b-adjustment, and the optimal adjustment set O*.

On an MPDAG the covariate sets that identify a total effect are characterised
by the *generalised (b-)adjustment criterion*: a set ``Z`` is valid for
``(X, Y)`` when it contains no forbidden node (nothing on or below a possibly
causal path from ``X`` to ``Y``) and blocks every proper non-causal
possibly-directed path from ``X`` to ``Y``.

Among the valid sets there is one -- the *optimal adjustment set* ``O*`` of
Henckel, Perkovic and Maathuis -- whose adjusted estimator has the smallest
asymptotic variance, uniformly over the SCMs compatible with the graph. It is
read off the graph as the parents of the *causal nodes*, minus the causal nodes
and the treatment itself.

``O*`` is a function of the graph. Change the graph -- which is what background
knowledge does -- and ``O*`` changes with it. Two distinct things can then go
wrong, and separating them is the structural contribution of this project:

* ``O*`` computed on the perturbed MPDAG is still valid but no longer optimal:
  the estimate stays unbiased and loses efficiency (radius ``delta_opt``);
* ``O*`` computed on the perturbed MPDAG is not valid at all: the estimate is
  biased (radius ``delta_valid``).
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from bkrobust.graphs.mpdag import MPDAG, Node


def is_valid_adjustment_set(
    graph: MPDAG,
    treatment: Node | Iterable[Node],
    outcome: Node | Iterable[Node],
    z: Iterable[Node],
) -> bool:
    """Whether ``z`` satisfies the generalised adjustment criterion on ``graph``.

    Args:
        graph: A DAG, CPDAG or MPDAG.
        treatment: Treatment node, or node set for a joint intervention.
        outcome: Outcome node or node set.
        z: The candidate adjustment set.

    Returns:
        ``True`` if adjusting for ``z`` identifies the total effect of
        ``treatment`` on ``outcome`` in every DAG represented by ``graph``.

    Raises:
        ValueError: If any named node is absent from ``graph``, or if
            ``treatment`` and ``outcome`` overlap.
    """
    raise NotImplementedError


def optimal_adjustment_set(
    graph: MPDAG,
    treatment: Node | Iterable[Node],
    outcome: Node | Iterable[Node],
) -> set[Node]:
    """Construct ``O*``, the optimal adjustment set of Henckel-Perkovic-Maathuis.

    Args:
        graph: A DAG, CPDAG or MPDAG.
        treatment: Treatment node or node set.
        outcome: Outcome node or node set.

    Returns:
        ``O*`` as a set of nodes. May be empty, which is a legitimate answer,
        not a failure.

    Raises:
        ValueError: If no valid adjustment set exists for ``(treatment,
            outcome)`` in ``graph`` -- the effect is then not identifiable by
            covariate adjustment and ``O*`` is undefined.
    """
    raise NotImplementedError


def causal_nodes(
    graph: MPDAG,
    treatment: Node | Iterable[Node],
    outcome: Node | Iterable[Node],
) -> set[Node]:
    """Return ``cn(X, Y)``: nodes on a proper possibly-causal path from X to Y.

    The building block of both the forbidden set and ``O*``.
    """
    raise NotImplementedError


def forbidden_set(
    graph: MPDAG,
    treatment: Node | Iterable[Node],
    outcome: Node | Iterable[Node],
) -> set[Node]:
    """Return the nodes no valid adjustment set may contain.

    The possible descendants of the causal nodes, together with ``treatment``.
    """
    raise NotImplementedError


def canonical_adjustment_set(
    graph: MPDAG,
    treatment: Node | Iterable[Node],
    outcome: Node | Iterable[Node],
) -> set[Node]:
    """Return the canonical (Pearl-style back-door) valid adjustment set.

    Valid whenever any valid set exists, and generally larger -- hence less
    efficient -- than ``O*``. Serves as the baseline against which the
    efficiency gap is reported.
    """
    raise NotImplementedError


def all_valid_adjustment_sets(
    graph: MPDAG,
    treatment: Node | Iterable[Node],
    outcome: Node | Iterable[Node],
    max_size: int | None = None,
) -> list[set[Node]]:
    """Enumerate every valid adjustment set, smallest first.

    Exponential in the number of nodes; intended only for the small canonical
    graphs used in tests and in the exact path of
    :mod:`bkrobust.theory.radius`.

    Args:
        graph: A DAG, CPDAG or MPDAG.
        treatment: Treatment node or node set.
        outcome: Outcome node or node set.
        max_size: Skip candidate sets larger than this.

    Returns:
        Valid adjustment sets, ordered by size then by sorted node labels so the
        output is deterministic.
    """
    raise NotImplementedError


def adjustment_set_exists(
    graph: MPDAG,
    treatment: Node | Iterable[Node],
    outcome: Node | Iterable[Node],
) -> bool:
    """Whether the total effect is identifiable by covariate adjustment.

    Cheaper than enumeration: checks the criterion against the canonical set,
    which is valid whenever any valid set is.
    """
    raise NotImplementedError
