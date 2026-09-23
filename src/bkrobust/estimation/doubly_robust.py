"""AIPW / DML with cross-fitting.

Included to close off an objection the paper will otherwise attract. A reader
can reasonably ask whether the bias attributed to a wrong adjustment set is
really bias from a misspecified outcome model, since the simulations are
linear-Gaussian and OLS is the obvious estimator. A doubly robust estimator
answers it: it is consistent when *either* the outcome model or the propensity
model is right, and it is still biased on an invalid adjustment set, because
double robustness is robustness to functional form, never to conditioning on
the wrong variables.

So the expected result is that the doubly robust arm tracks the OLS arm past
``delta_valid``. That is not a null result -- it is the demonstration that the
failure is structural, and that reaching for a better estimator does not fix it.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    import numpy as np
    import pandas as pd

    from bkrobust.estimation.base import EffectEstimate
    from bkrobust.graphs.mpdag import Node


class DoublyRobustEstimator:
    """Cross-fitted AIPW / DML estimator of the ATE.

    Args:
        method: ``"aipw"`` or ``"dml"``.
        outcome_model: Sklearn-style regressor, or ``None`` for the config
            default.
        propensity_model: Sklearn-style classifier, or ``None``.
        n_folds: Cross-fitting folds.
        propensity_clip: ``(low, high)`` bounds on estimated propensities.
            Without clipping, near-zero propensities blow the variance up; with
            it, the estimand shifts slightly. The bounds are reported in the
            diagnostics rather than buried.
        alpha: Confidence level.
        random_state: Seed for the fold split. Must be derived from the
            experiment seed, never left to a global default.
    """

    name: str = "doubly_robust"
    requires_scm: bool = False

    def __init__(
        self,
        method: str = "aipw",
        outcome_model: Any | None = None,
        propensity_model: Any | None = None,
        n_folds: int = 5,
        propensity_clip: tuple[float, float] = (0.01, 0.99),
        alpha: float = 0.05,
        random_state: int | None = None,
    ) -> None:
        raise NotImplementedError

    def fit(
        self,
        data: pd.DataFrame,
        treatment: Node,
        outcome: Node,
        adjustment_set: Iterable[Node],
        **kwargs: Any,
    ) -> DoublyRobustEstimator:
        """Cross-fit the nuisance models and form the influence-function scores.

        Raises:
            ValueError: If the treatment is not binary -- the propensity model
                requires it -- or if a named column is missing.
        """
        raise NotImplementedError

    def estimate(self) -> EffectEstimate:
        """Return the AIPW estimate with its influence-function standard error.

        Diagnostics carry the propensity overlap summary and the per-fold
        nuisance scores, so poor overlap is visible rather than silently
        inflating the interval.

        Raises:
            RuntimeError: If called before :meth:`fit`.
        """
        raise NotImplementedError

    def estimate_cate(self, covariates: pd.DataFrame) -> np.ndarray:
        """Conditional effects from the fitted pseudo-outcome regression.

        Raises:
            RuntimeError: If called before :meth:`fit`.
        """
        raise NotImplementedError

    @property
    def influence_function(self) -> np.ndarray:
        """Per-observation influence-function values.

        Raises:
            RuntimeError: If called before :meth:`fit`.
        """
        raise NotImplementedError
