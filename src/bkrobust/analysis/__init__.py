"""Hypothesis analysis for the synthetic-ensemble session.

Consumes the frozen instance schema (:mod:`bkrobust.core.instance`) and computes
the pre-registered statistics for H1-H6. Kept separate from the generators so
that the analysis cannot be tuned to the data by whoever wrote the generator.
"""

from __future__ import annotations
