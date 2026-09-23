"""A taxonomy of the ways background knowledge can be wrong.

Two orthogonal axes.

**Kind** -- what sort of claim is wrong:

``ORIENTATION``
    An edge is asserted in the wrong direction. The atomic case; every other
    kind ultimately expresses itself as some set of wrong orientations after
    the closure runs.
``ANCESTRAL``
    A reachability claim is wrong. Weaker than an orientation claim locally,
    but its consequences surface far from the pair named, because the closure
    must route a directed path.
``TIER``
    A variable is placed in the wrong tier. One statement, many forbidden
    edges; the highest ratio of expert confidence to blast radius.
``SPURIOUS``
    A constraint is asserted that the true DAG does not entail -- the
    over-confident expert.
``MISSING``
    A constraint the true DAG does entail is withheld -- the under-informative
    expert. Note this arm is consistent *and true*: it costs efficiency but
    never validity, which makes it the right control against which to read the
    other kinds.

**Severity** -- what the wrong claim costs, which is what the two radii
separate:

``BENIGN``
    ``O*`` is unchanged; the knowledge is wrong and nothing downstream notices.
``EFFICIENCY``
    ``O*`` changes but stays valid. Unbiased, larger variance. The regime
    between ``delta_opt`` and ``delta_valid``.
``VALIDITY``
    ``O*`` is no longer a valid adjustment set. Biased at any sample size. At
    and beyond ``delta_valid``.
``UNIDENTIFIABLE``
    No valid adjustment set survives in the perturbed MPDAG. The estimand
    cannot be recovered by adjustment at all -- rarer, but it is the one case
    an analyst could in principle detect without knowing the truth.
"""

from __future__ import annotations

from enum import StrEnum
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from bkrobust.graphs.mpdag import MPDAG, Node
    from bkrobust.knowledge.base import BackgroundKnowledge


class MisspecificationKind(StrEnum):
    """What sort of claim the knowledge gets wrong. Values match config strings."""

    ORIENTATION = "orientation"
    ANCESTRAL = "ancestral"
    TIER = "tier"
    SPURIOUS = "spurious"
    MISSING = "missing"


class Severity(StrEnum):
    """What the wrong claim costs downstream, ordered from cheapest to worst."""

    BENIGN = "benign"
    EFFICIENCY = "efficiency"
    VALIDITY = "validity"
    UNIDENTIFIABLE = "unidentifiable"


def classify(
    cpdag: MPDAG,
    true_graph: MPDAG,
    bk: BackgroundKnowledge,
    treatment: Node,
    outcome: Node,
) -> Severity:
    """Assign a severity to a knowledge set, given the truth.

    Imposes ``bk`` on ``cpdag``, reads ``O*`` off the result, and asks what that
    set does in ``true_graph``: unchanged (``BENIGN``), changed but still valid
    (``EFFICIENCY``), invalid (``VALIDITY``), or no valid set left at all
    (``UNIDENTIFIABLE``).

    Args:
        cpdag: The CPDAG the knowledge is imposed on.
        true_graph: The ground-truth DAG.
        bk: The knowledge to classify.
        treatment: Treatment node.
        outcome: Outcome node.

    Returns:
        The severity.

    Raises:
        ValueError: If ``bk`` is inconsistent with ``cpdag`` -- inconsistent
            knowledge is out of scope, it fails the check practitioners already
            run.
    """
    raise NotImplementedError


def kind_of(bk: BackgroundKnowledge, true_graph: MPDAG) -> set[MisspecificationKind]:
    """Which misspecification kinds a knowledge set exhibits.

    A knowledge set can exhibit several at once; the sweeps sample one kind at a
    time, but elicited knowledge from :mod:`bkrobust.knowledge.elicit` will not
    be so tidy.

    Args:
        bk: The knowledge to inspect.
        true_graph: The ground-truth DAG.

    Returns:
        The kinds present; empty if ``bk`` is entirely true of ``true_graph``.
    """
    raise NotImplementedError


def severity_order(severity: Severity) -> int:
    """Map a severity to an integer rank for sorting and plotting.

    ``BENIGN`` is 0 and ``UNIDENTIFIABLE`` is 3.
    """
    raise NotImplementedError
