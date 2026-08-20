"""Variance inflation in the band between delta_opt and delta_valid.

Between the two radii the analyst is in an odd position: their knowledge is
wrong, their adjustment set is wrong, and their estimate is nonetheless
consistent. The estimator converges to the right number. What they lose is
efficiency -- the wrong set has larger asymptotic variance than ``O*`` -- and
nothing in the output signals it. There is no diagnostic to run, because
nothing is broken in the sense that diagnostics test for. The estimate is
simply noisier than it needed to be, and the analyst has no way to know by how
much.

This module quantifies that loss, in the unit that makes it concrete: the
**effective sample size ratio**. A variance inflation of 1.4 means the analyst
threw away 29% of their data by believing something false, and stating it that
way is what turns an asymptotic-variance comparison into a claim a practitioner
can act on.

See ``docs/THEORY.md`` (target **T3**).
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    from bkrobust.graphs.mpdag import MPDAG, Node


@dataclass(frozen=True)
class EfficiencyGap:
    """Efficiency lost by adjusting on the wrong-but-valid set.

    Attributes:
        variance_ratio: ``asVar(O*_perturbed) / asVar(O*_true)``; ``>= 1`` when
            the perturbed set is valid, since ``O*`` is the minimiser.
        effective_sample_size_ratio: ``1 / variance_ratio`` -- the fraction of
            the sample that survives the mistake.
        delta: The radius at which this gap was measured.
        perturbed_set: The adjustment set chosen under false knowledge.
        optimal_set: ``O*`` under true knowledge.
        both_valid: Whether both sets are valid in the true graph. When
            ``False`` the ratio is not an efficiency comparison at all -- one of
            the estimators is biased -- and the number must not be reported as
            one.
    """

    variance_ratio: float
    effective_sample_size_ratio: float
    delta: int
    perturbed_set: frozenset[Node]
    optimal_set: frozenset[Node]
    both_valid: bool


def efficiency_gap(
    cpdag: MPDAG,
    true_graph: MPDAG,
    treatment: Node,
    outcome: Node,
    scm: Any,
    delta: int,
    *,
    distance: str = "orientation_flips",
    aggregate: str = "worst",
) -> EfficiencyGap:
    """Efficiency lost at radius ``delta``, over the ball of knowledge sets.

    Args:
        cpdag: The CPDAG the knowledge is imposed on.
        true_graph: The ground-truth DAG.
        treatment: Treatment node.
        outcome: Outcome node.
        scm: The SCM supplying the variances.
        delta: The radius.
        distance: Metric the radius is measured in.
        aggregate: ``"worst"``, ``"mean"`` or ``"median"`` over the ball.
            ``"worst"`` matches how the radii are defined and is the default;
            the others describe the typical case.

    Returns:
        An :class:`EfficiencyGap`.

    Raises:
        ValueError: If no knowledge set at radius ``delta`` yields a valid
            adjustment set -- ``delta`` is then past ``delta_valid`` and the
            question is about bias, not efficiency.
    """
    raise NotImplementedError


def variance_inflation_curve(
    cpdag: MPDAG,
    true_graph: MPDAG,
    treatment: Node,
    outcome: Node,
    scm: Any,
    delta_grid: Iterable[int],
    *,
    distance: str = "orientation_flips",
    aggregate: str = "worst",
) -> dict[int, EfficiencyGap]:
    """Variance inflation at each radius, over the band.

    Paired with :func:`bkrobust.theory.bias_bound.bias_at_radius_curve` this
    gives the paper's summary figure: one curve that lifts off zero at
    ``delta_opt`` (variance) and a second that lifts off at ``delta_valid``
    (bias), with the band between them shaded.

    Returns:
        Mapping from radius to gap. Radii past ``delta_valid`` are omitted
        rather than reported as infinite.
    """
    raise NotImplementedError


def band_width(
    cpdag: MPDAG,
    true_graph: MPDAG,
    treatment: Node,
    outcome: Node,
    scm: Any,
    *,
    distance: str = "orientation_flips",
    method: str = "exact",
) -> int:
    """Width of the silent-degradation band, ``delta_valid - delta_opt``.

    How much wrong knowledge a practitioner can hold while paying only in
    variance -- that is, while their diagnostics stay clean. A wide band means
    the failure mode is common and invisible; a narrow one means wrong knowledge
    tends to go straight to bias.

    Returns:
        The width, or :data:`bkrobust.theory.radius.UNREACHED` if either radius
        was not found.
    """
    raise NotImplementedError


def maximum_efficiency_loss(
    cpdag: MPDAG,
    true_graph: MPDAG,
    treatment: Node,
    outcome: Node,
    scm: Any,
    *,
    distance: str = "orientation_flips",
) -> float:
    """Largest variance ratio achievable while staying valid.

    The worst the analyst can do without becoming biased -- the ceiling on the
    price of silent degradation.
    """
    raise NotImplementedError
