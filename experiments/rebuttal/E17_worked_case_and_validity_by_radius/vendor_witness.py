"""Claim-retraction witness, vendored for E17.

Source: experiments/witness_table.py (functions ``witness`` and
``render``, constants MAX_DEPTH and SUBSET_BUDGET), copied unchanged except that ``render``
writes raw node ids (the released ``witness`` column uses raw ids).
Source repository HEAD: ba7ddc12c446f1087e58dc428720c84d33b1b622 (git -C bkrobust (main branch) rev-parse HEAD)

Every import below exists in the submitted src/ (rebuttal branch).
"""

from __future__ import annotations

import itertools

from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag
from bkrobust.mpdag_criterion.criterion import is_amenable

MAX_DEPTH = 2
SUBSET_BUDGET = 400


def witness(
    cpdag: MPDAG, k: list, x: str, y: str, z: frozenset
) -> tuple[int | None, tuple, int, str]:
    """The smallest retraction that invalidates ``z``, and how many there are at that depth.

    Returns ``(depth, the first such subset, how many subsets break it at that
    depth, what the retraction costs)``. The cost is ``validity`` when the set
    merely stops being valid and ``identifiability`` when the query also stops
    being amenable.
    """
    used = 0
    for depth in range(1, MAX_DEPTH + 1):
        if depth > len(k):
            break
        hits: list[tuple] = []
        cost = "validity"
        for subset in itertools.combinations(range(len(k)), depth):
            if used >= SUBSET_BUDGET:
                break
            drop = set(subset)
            kept = [e for i, e in enumerate(k) if i not in drop]
            g = apply_orientations(cpdag, kept)
            used += 1
            if g is None:
                continue
            if not is_gac_valid_mpdag(g, x, y, z):
                hits.append(tuple(k[i] for i in subset))
                if not hits[:-1] and not is_amenable(g, x, y):
                    cost = "identifiability"
        if hits:
            return depth, hits[0], len(hits), cost
    return None, (), 0, ""


def render(edges: tuple) -> str:
    """The retraction set as the released records write it."""
    return " and ".join(f"{a}->{b}" for a, b in edges)
