"""Shared plumbing: seeding, results IO, logging, plotting.

Nothing here is about causality. It is the machinery that makes a result
reproducible and a figure consistent, and it is centralised so that no
experiment quietly does its own thing with a seed or a results directory.
"""

from __future__ import annotations

__all__: list[str] = []
