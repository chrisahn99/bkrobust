"""Shared helpers for the rebuttal experiments.

Everything here calls the submitted code (branch ``iclr_2027``, ``5421cab``):
the hybrid radius, the GAC predicate, Meek closure and the published-network
loaders. Nothing in ``src/`` is modified.
"""

from __future__ import annotations

import json
import sys
import time
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag  # noqa: E402
from bkrobust.hybrid import breakdown_radius  # noqa: E402

Edge = tuple[str, str]
OUT = ROOT / "results" / "rebuttal"
PAPER = ROOT.parent / "6ab5a4476c152c4289d9bb25"
RELEASE = PAPER / "release"


def k_of(cpdag: MPDAG, g: MPDAG) -> frozenset[Edge]:
    """Closure orientations of a state: its directed edges not compelled in the CPDAG."""
    return frozenset(set(g.directed_edges) - set(cpdag.directed_edges))


def dist(cpdag: MPDAG, g: MPDAG, h: MPDAG) -> int:
    """Corollary 1(a): d(G, H) = |K_G symmetric-difference K_H|."""
    return len(k_of(cpdag, g) ^ k_of(cpdag, h))


def is_closed(cpdag: MPDAG, s: frozenset[Edge]) -> bool:
    g = apply_orientations(cpdag, sorted(s))
    return g is not None and k_of(cpdag, g) == s


def radius_bruteforce_retractions(cpdag, g0, x, y, z, max_k=None):
    """Least number of retractions from G0 that invalidates Z, by plain subset
    enumeration over closed subsets of K_G0 (independent of the hybrid code).
    Returns (radius or None if none within max_k, list of witness sets)."""
    kg0 = sorted(k_of(cpdag, g0))
    top = len(kg0) if max_k is None else min(max_k, len(kg0))
    for k in range(0, top + 1):
        wit = []
        for drop in combinations(kg0, k):
            s = frozenset(kg0) - set(drop)
            g = apply_orientations(cpdag, sorted(s))
            if g is None or k_of(cpdag, g) != s:
                continue
            if not is_gac_valid_mpdag(g, x, y, z):
                wit.append(tuple(sorted(drop)))
        if wit:
            return k, wit
    return None, []


def claim_radius(cpdag, k_claims, x, y, z, max_k=None):
    """Fewest asserted claims whose withdrawal (then Meek closure) breaks Z."""
    k_claims = sorted(k_claims)
    top = len(k_claims) if max_k is None else min(max_k, len(k_claims))
    for k in range(1, top + 1):
        wit = []
        for drop in combinations(k_claims, k):
            g = apply_orientations(cpdag, sorted(set(k_claims) - set(drop)))
            if g is not None and not is_gac_valid_mpdag(g, x, y, z):
                wit.append(drop)
        if wit:
            return k, wit
    return None, []


def hybrid(cpdag, g0, x, y, z, **kw):
    return breakdown_radius(cpdag, None, x, y, frozenset(z), g0=g0, **kw)


def rad_str(r) -> str:
    return "inf" if r == UNREACHED else str(r)


def dump(obj, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, default=str))


class Timer:
    def __enter__(self):
        self.t = time.perf_counter()
        return self

    def __exit__(self, *a):
        self.s = time.perf_counter() - self.t
