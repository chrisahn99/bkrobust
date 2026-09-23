"""Conformal coverage wrapper -- CONDITIONAL component, may be cut.

Split-conformal intervals around an effect estimate, to ask whether coverage
degrades gracefully or collapses as the adjustment set goes wrong.

The honest framing, which the docstrings here should not soften: conformal
prediction guarantees coverage of an *observable* quantity under exchangeability.
A causal effect is not observable, so no conformal method can guarantee coverage
of it without an identification assumption -- and the identification assumption
is exactly what a wrong adjustment set breaks. What this module can measure is
how badly a nominal-95% procedure actually does once its identifying assumption
is false, and how visible the failure is from inside. If the answer is that
coverage collapses silently, that is a finding worth reporting and a reason for
caution about conformal causal intervals in general. If it degrades gradually,
that is also worth reporting.

Either way this is a secondary result. Do not let it crowd out the main claim.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    import numpy as np
    import pandas as pd

    from bkrobust.estimation.base import EffectEstimate, Estimator
    from bkrobust.graphs.mpdag import Node


class ConformalWrapper:
    """Wrap any estimator in split-conformal prediction intervals.

    Args:
        base: The estimator to wrap.
        alpha: Miscoverage level; ``0.05`` targets 95% coverage.
        calibration_fraction: Fraction of the data held out for calibration.
        score: Conformity score -- ``"absolute"`` or ``"normalized"``.
        random_state: Seed for the split.
    """

    name: str = "conformal"
    requires_scm: bool = False

    def __init__(
        self,
        base: Estimator,
        alpha: float = 0.05,
        calibration_fraction: float = 0.3,
        score: str = "absolute",
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
    ) -> ConformalWrapper:
        """Fit the base estimator on the training split and calibrate on the rest."""
        raise NotImplementedError

    def estimate(self) -> EffectEstimate:
        """Return the base point estimate with a conformal interval.

        Raises:
            RuntimeError: If called before :meth:`fit`.
        """
        raise NotImplementedError

    def estimate_cate(self, covariates: pd.DataFrame) -> np.ndarray:
        """Delegate to the base estimator."""
        raise NotImplementedError

    def coverage(self, data: pd.DataFrame, truth: np.ndarray) -> float:
        """Empirical coverage on held-out data with known ground truth.

        Simulation only -- ``truth`` is unavailable in any real application,
        which is the whole difficulty this module documents.
        """
        raise NotImplementedError
