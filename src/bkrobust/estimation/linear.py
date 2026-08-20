"""OLS adjustment for the total effect.

Regress the outcome on the treatment and the adjustment set; the coefficient on
the treatment is the estimate. Under a linear-Gaussian SCM and a valid
adjustment set this is exactly the estimator the theory in
:mod:`bkrobust.graphs.optimality` computes the asymptotic variance of, which
makes it the workhorse of the simulation studies: predicted variance and
realised variance are directly comparable, so a mismatch is a bug in one of the
two and not a modelling judgement call.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    import numpy as np
    import pandas as pd

    from bkrobust.estimation.base import EffectEstimate
    from bkrobust.graphs.mpdag import Node


class LinearAdjustmentEstimator:
    """Total effect by ordinary least squares with covariate adjustment.

    Args:
        fit_intercept: Include an intercept.
        robust_se: Heteroscedasticity-consistent standard errors --
            ``"HC0"``-``"HC3"``, or ``None`` for the classical formula.
        alpha: Confidence level for the reported interval.
    """

    name: str = "ols_adjustment"
    requires_scm: bool = False

    def __init__(
        self,
        fit_intercept: bool = True,
        robust_se: str | None = "HC3",
        alpha: float = 0.05,
    ) -> None:
        raise NotImplementedError

    def fit(
        self,
        data: pd.DataFrame,
        treatment: Node,
        outcome: Node,
        adjustment_set: Iterable[Node],
        **kwargs: Any,
    ) -> LinearAdjustmentEstimator:
        """Fit ``outcome ~ treatment + adjustment_set``.

        Raises:
            ValueError: If a column is missing, if the design matrix is rank
                deficient, or if ``adjustment_set`` overlaps treatment or outcome.
        """
        raise NotImplementedError

    def estimate(self) -> EffectEstimate:
        """Return the treatment coefficient as the effect estimate.

        Raises:
            RuntimeError: If called before :meth:`fit`.
        """
        raise NotImplementedError

    def estimate_cate(self, covariates: pd.DataFrame) -> np.ndarray:
        """Not available: the model is additive and its CATE is constant.

        Raises:
            NotImplementedError: Always. Use
                :class:`~bkrobust.estimation.doubly_robust.DoublyRobustEstimator`
                with a flexible outcome model for heterogeneous effects.
        """
        raise NotImplementedError

    @property
    def residual_variance(self) -> float:
        """Residual variance of the fitted regression.

        Exposed so the realised variance can be checked against the asymptotic
        variance predicted by :mod:`bkrobust.graphs.optimality`.

        Raises:
            RuntimeError: If called before :meth:`fit`.
        """
        raise NotImplementedError


def total_effect_ols(
    data: pd.DataFrame,
    treatment: Node,
    outcome: Node,
    adjustment_set: Iterable[Node],
    *,
    robust_se: str | None = "HC3",
) -> EffectEstimate:
    """One-shot convenience wrapper around :class:`LinearAdjustmentEstimator`."""
    raise NotImplementedError
