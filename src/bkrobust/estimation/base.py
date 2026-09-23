"""The Estimator protocol.

Every estimator in this package -- classical or amortized -- takes data, a
treatment, an outcome and an adjustment set, and returns a point estimate with
an uncertainty interval. Causal foundation models are wrapped to this same
interface in :mod:`bkrobust.cfm.adapters` so that the audit compares like with
like: whatever a foundation model does internally, from the outside it is a
function from (data, knowledge) to an effect estimate, and that is the function
under test.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Protocol, runtime_checkable

if TYPE_CHECKING:  # pragma: no cover
    import numpy as np
    import pandas as pd

    from bkrobust.graphs.mpdag import Node


@dataclass(frozen=True)
class EffectEstimate:
    """A point estimate of a causal effect with its uncertainty.

    Attributes:
        point: The point estimate.
        std_error: Standard error, or ``None`` if the estimator reports none.
            Most foundation-model checkpoints do not.
        ci_lower: Lower confidence bound, or ``None``.
        ci_upper: Upper confidence bound, or ``None``.
        alpha: The level the interval was computed at.
        estimand: ``"ate"``, ``"att"``, ``"cate"`` or ``"total_effect"``.
        adjustment_set: The set actually adjusted on. Recorded on the estimate
            itself so a number can never drift apart from the set that produced
            it once results are aggregated.
        n_samples: Sample size used.
        diagnostics: Estimator-specific extras -- fold scores, propensity
            overlap, convergence flags.
    """

    point: float
    std_error: float | None = None
    ci_lower: float | None = None
    ci_upper: float | None = None
    alpha: float = 0.05
    estimand: str = "ate"
    adjustment_set: frozenset[Node] = frozenset()
    n_samples: int = 0
    diagnostics: dict[str, Any] = field(default_factory=dict)

    def covers(self, truth: float) -> bool | None:
        """Whether the interval contains ``truth``.

        Returns:
            ``None`` when the estimator reported no interval, so that
            "no interval" is never silently aggregated as "did not cover".
        """
        raise NotImplementedError


@runtime_checkable
class Estimator(Protocol):
    """Anything that estimates an effect given data and an adjustment set."""

    name: str
    requires_scm: bool

    def fit(
        self,
        data: pd.DataFrame,
        treatment: Node,
        outcome: Node,
        adjustment_set: Iterable[Node],
        **kwargs: Any,
    ) -> Estimator:
        """Fit on ``data``, adjusting for ``adjustment_set``.

        Args:
            data: Observational data; columns are variable names.
            treatment: Treatment column.
            outcome: Outcome column.
            adjustment_set: Covariates to adjust for. May be empty.
            **kwargs: Estimator-specific options.

        Returns:
            ``self``, fitted.

        Raises:
            ValueError: If a named column is missing, or the adjustment set
                overlaps treatment or outcome.
        """
        ...

    def estimate(self) -> EffectEstimate:
        """Return the effect estimate.

        Raises:
            RuntimeError: If called before :meth:`fit`.
        """
        ...

    def estimate_cate(self, covariates: pd.DataFrame) -> np.ndarray:
        """Conditional effects at the given covariate values.

        Raises:
            NotImplementedError: If the estimator targets only the ATE.
        """
        ...


def build_estimator(cfg: Any) -> Estimator:
    """Instantiate the estimator named by a Hydra config node.

    Args:
        cfg: An ``estimator`` config group -- see ``configs/estimator/``.

    Returns:
        An unfitted estimator.

    Raises:
        ImportError: If the config selects a checkpoint-backed estimator and the
            ``[cfm]`` extra is not installed. The message must name the extra,
            since this is the most likely install-time failure for anyone
            reproducing the results.
        KeyError: If ``cfg.target`` names no known estimator.
    """
    raise NotImplementedError


def compare_estimators(
    estimators: Sequence[Estimator],
    data: pd.DataFrame,
    treatment: Node,
    outcome: Node,
    adjustment_set: Iterable[Node],
) -> Any:
    """Fit several estimators on identical data and adjustment set.

    Returns:
        A ``pandas.DataFrame``, one row per estimator.
    """
    raise NotImplementedError
