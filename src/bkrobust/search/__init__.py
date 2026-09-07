"""Axis B: exact search for the breakdown radius, faster than shell-by-shell BFS.

Exactness is non-negotiable. Every method here either returns the BFS answer or
an explicitly labelled bound.
"""

from __future__ import annotations
