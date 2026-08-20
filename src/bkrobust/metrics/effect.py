"""Effect-estimation metrics: bias, PEHE, coverage.

Reporting conventions this module enforces, because each has a way of going
quietly wrong in aggregate:

* **Signed and absolute bias, both.** Signed bias averages toward zero when the
  perturbations are symmetric, which reads as "no bias" and is not. Absolute
  bias alone hides direction. Report both.
* **Bias in effect units and in standard-error units.** The first is what the
  domain cares about; the second is what determines whether anyone would notice.
* **Coverage against a nominal level.** A procedure with 60% coverage at nominal
  95% is broken in a way that a bias number alone does not convey.
* **Never aggregate a missing interval as a non-covering one.** Most foundation
  models report no interval; counting those as failures would fabricate a
  finding.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    import numpy as np

    from bkrobust.estimation.base import EffectEstimate


def ate_bias(estimates: Sequence[float], truth: float) -> float:
    """Mean signed bias, ``mean(estimates) - truth``."""
    raise NotImplementedError


def ate_absolute_bias(estimates: Sequence[float], truth: float) -> float:
    """Mean absolute bias, ``mean(|estimate - truth|)``."""
    raise NotImplementedError


def ate_rmse(estimates: Sequence[float], truth: float) -> float:
    """Root mean squared error about the truth, combining bias and variance."""
    raise NotImplementedError


def standardized_bias(estimates: Sequence[float], truth: float, se: float) -> float:
    """Bias in units of the estimator's own standard error.

    The visibility scale: a bias below one standard error is invisible to the
    analyst, one above three is not, and the same absolute bias can fall on
    either side depending on sample size.

    Raises:
        ValueError: If ``se`` is not positive.
    """
    raise NotImplementedError


def pehe(predicted_ite: np.ndarray, true_ite: np.ndarray) -> float:
    """Root precision in estimating heterogeneous effects.

    Raises:
        ValueError: If the arrays have different lengths.
    """
    raise NotImplementedError


def coverage(estimates: Sequence[EffectEstimate], truth: float) -> float | None:
    """Fraction of intervals containing ``truth``.

    Returns:
        The coverage rate, or ``None`` if no estimate reported an interval.
        Estimates without intervals are excluded from the denominator rather
        than counted as failures.
    """
    raise NotImplementedError


def interval_width(estimates: Sequence[EffectEstimate]) -> float | None:
    """Mean interval width, over the estimates that reported one.

    Read next to :func:`coverage`: coverage bought by wide intervals is not the
    same achievement as coverage at a competitive width.
    """
    raise NotImplementedError


def effect_spread(estimates: Sequence[float]) -> dict[str, float]:
    """Spread of estimates across a ball of knowledge sets.

    The real-data substitute for bias. Ground truth is unavailable there, but
    the range of conclusions reachable by knowledge that all passes the
    consistency check is observable, and it is the number that should accompany
    any effect estimate derived from elicited knowledge.

    Returns:
        Keys ``"min"``, ``"max"``, ``"range"``, ``"iqr"``, ``"std"``.
    """
    raise NotImplementedError


def bias_summary(
    estimates: Sequence[EffectEstimate],
    truth: float,
) -> dict[str, float | None]:
    """Every effect metric at once, for one experimental cell.

    Returns:
        Keys ``"bias"``, ``"abs_bias"``, ``"rmse"``, ``"standardized_bias"``,
        ``"coverage"``, ``"interval_width"``, ``"n"``.
    """
    raise NotImplementedError
