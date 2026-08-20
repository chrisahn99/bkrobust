"""Effect estimators evaluated under a supplied adjustment set.

The adjustment set is an *input* here, never a choice. It is chosen upstream by
:mod:`bkrobust.graphs.adjustment` from an MPDAG that background knowledge --
possibly false -- produced. Keeping that separation strict is what lets the
experiments attribute error to the knowledge rather than to the estimator: the
same estimator, the same data, only the set changes.

Module map:

``base``
    The :class:`~bkrobust.estimation.base.Estimator` protocol.
``linear``
    OLS adjustment for the total effect.
``doubly_robust``
    AIPW / DML, which is robust to model misspecification and, importantly, not
    to adjustment-set misspecification -- a point the experiments make directly.
``conformal``
    Conditional component: coverage under a wrong adjustment set.
"""

from __future__ import annotations

__all__: list[str] = []
