"""Declarative (CP-SAT) encodings of the breakdown radius.

Item B2.8 from the first session's list: the one scaling direction never
attempted. The point is tractability -- computing an exact radius without ever
building the perturbation space, which is what limits every enumerative method
here well before the search itself does.
"""

from __future__ import annotations
