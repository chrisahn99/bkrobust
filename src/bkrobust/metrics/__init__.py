"""Evaluation metrics, split by what they measure.

``effect``
    How wrong the estimate is -- ATE bias, PEHE, interval coverage.
``structural``
    How wrong the graph and the adjustment set are -- validity flip rate,
    optimality loss, cascade size.

The pairing is the point. A structural change that never reaches ``O*`` costs
nothing measurable, and an effect error with no structural cause is a bug in the
estimator. Reporting both together is what lets a result be attributed.
"""

from __future__ import annotations

__all__: list[str] = []
