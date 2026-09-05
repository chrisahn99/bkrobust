"""One entry point for the exact breakdown radius.

Sessions 1-4 produced several methods with different strengths, and choosing
between them required knowing their measured behaviour. This module removes that
burden from the caller: :func:`breakdown_radius` takes the problem and returns
the radius together with the method and oracle that produced it and the
assumptions that answer carries.

**Why a hybrid at all.** Session 3 measured a clean split, and session 4's
profiling explained it:

* When a failure is near, upward BFS finds it in one expansion. It is
  unbeatable there, and the SAT encoding must still build an ``O(n^4)``-clause
  model before it can say anything -- roughly 37x slower.
* When there is no failure anywhere, upward BFS must exhaust the entire up-set.
  It timed out on 26 of 49 such instances at a 60 s cap, while E1 answered all
  of them, the worst in 8.7 s.

The discriminator -- has a failure been found yet -- is free at runtime. So the
dispatch is: run the search under a small depth budget; if it finds a failure,
that is the answer. If it exhausts the budget without one, hand the instance to
the E1 ladder, whose UNSAT rungs certify the clean shells directly.

**Two accelerations, both applied by default and both differentially tested.**

* :mod:`bkrobust.search.exact_fast` decides cover minimality by set comparison
  (Lemma O) rather than by enumerating DAG extensions -- 70.5% of the frozen
  search's measured time, and 3.77x in aggregate.
* :mod:`bkrobust.mpdag_criterion` decides validity on the MPDAG directly rather
  than by enumerating ``[G]`` and checking every extension.

**What the answer assumes.** Both search methods are upward searches, so both
are exact *iff Conjecture 2 holds*, which rests on Anti-Exchange Case B --
verified, not proved (``THEOREMS.md`` sections 4c and 6). The error is
one-sided: they can only ever return a radius that is too **large**, never too
small, i.e. they can overstate robustness but never understate it. Every result
carries this in :attr:`HybridResult.assumes`; do not drop it when the number is
copied into a table.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

from bkrobust.core.conventions import UNREACHED
from bkrobust.core.oracle import is_valid
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.mpdag_criterion import is_valid_mpdag
from bkrobust.sat.e1 import radius_e1
from bkrobust.search.exact import SearchStats
from bkrobust.search.exact_fast import radius_local_up_fast

Edge = tuple[str, str]

#: Depth explored by the search before falling through to the ladder. Session 3
#: measured ~80% of finite radii at 1 and essentially all at <= 3, so a budget
#: of 3 catches the cases the search is good at while bounding the up-set it can
#: be forced to exhaust.
DEFAULT_SEARCH_BUDGET = 3


@dataclass
class HybridResult:
    """A radius, how it was obtained, and what it assumes.

    Attributes:
        radius: The exact radius, or :data:`~bkrobust.core.conventions.UNREACHED`
            when no reachable perturbation makes ``Z`` invalid. ``UNREACHED`` is
            a sentinel, not a number: never average or plot it numerically.
        method: ``"local_up_fast"`` when the bounded search found the failure,
            ``"e1_ladder"`` when the search exhausted its budget and the
            encoding answered, ``"degenerate"`` when ``Z`` already fails at
            ``G0``.
        oracle: ``"mpdag_criterion"`` or ``"enumeration"``.
        exact: False only if the ladder timed out; the radius is then not a
            completed search.
        assumes: The assumption chain the answer inherits.
        witness: A nearest failing graph as an edge string, when one was found.
        search_seconds: Time in the bounded upward search.
        ladder_seconds: Time in the E1 ladder, 0.0 if it was not needed.
        stats: Cost counters from the search leg.
    """

    radius: int
    method: str
    oracle: str
    exact: bool = True
    assumes: str = "Conjecture 2 (hence Anti-Exchange Case B, verified not proved)"
    witness: str | None = None
    search_seconds: float = 0.0
    ladder_seconds: float = 0.0
    stats: SearchStats = field(default_factory=SearchStats)

    @property
    def total_seconds(self) -> float:
        """Search plus ladder."""
        return self.search_seconds + self.ladder_seconds


def breakdown_radius(
    cpdag: MPDAG,
    knowledge: list[Edge] | None,
    x: str,
    y: str,
    z: frozenset[str],
    *,
    g0: MPDAG | None = None,
    search_budget: int = DEFAULT_SEARCH_BUDGET,
    use_criterion: bool = True,
    time_limit_s: float = 300.0,
) -> HybridResult:
    """Compute the exact breakdown radius, dispatching between the two methods.

    Args:
        cpdag: The estimated CPDAG.
        knowledge: The analyst's asserted orientations, used to build ``G0`` as
            ``Meek(cpdag, knowledge)``. Ignored when ``g0`` is given.
        x: Treatment.
        y: Outcome.
        z: The adjustment set, fixed once from ``G0`` and then held.
        g0: The analyst's graph, if already built. Supplying it skips the
            closure and is the cheaper path when it is already to hand.
        search_budget: Depth explored before falling through to the ladder.
        use_criterion: Decide validity on the MPDAG directly. Setting this False
            falls back to the enumeration oracle, which is far slower but shares
            no code with the criterion -- useful for differential testing.
        time_limit_s: Per-rung limit for the ladder leg.

    Returns:
        A :class:`HybridResult`.

    Raises:
        ValueError: If ``knowledge`` is inconsistent with ``cpdag``, so that no
            ``G0`` exists.
    """
    if g0 is None:
        built = apply_orientations(cpdag, knowledge or [])
        if built is None:
            raise ValueError("knowledge is inconsistent with the CPDAG; no G0 exists")
        g0 = built

    if use_criterion:

        def fails(g: MPDAG) -> bool:
            return not is_valid_mpdag(g, x, y, z)

        oracle = "mpdag_criterion"
    else:

        def fails(g: MPDAG) -> bool:
            return not is_valid(z, g, x, y)

        oracle = "enumeration"

    stats = SearchStats()
    t0 = time.perf_counter()
    found = radius_local_up_fast(cpdag, g0, fails, max_depth=search_budget, stats=stats)
    search_s = time.perf_counter() - t0

    if found.radius != UNREACHED:
        method = "degenerate" if found.radius == 0 else "local_up_fast"
        return HybridResult(
            radius=found.radius,
            method=method,
            oracle=oracle,
            exact=True,
            witness=found.witness,
            search_seconds=search_s,
            stats=stats,
        )

    if found.exact:
        # The budget was never hit: the search exhausted the whole up-set and
        # proved there is no failure. No need for the ladder.
        return HybridResult(
            radius=UNREACHED,
            method="local_up_fast",
            oracle=oracle,
            exact=True,
            search_seconds=search_s,
            stats=stats,
        )

    t0 = time.perf_counter()
    lad = radius_e1(cpdag, g0, x, y, z, time_limit_s=time_limit_s)
    ladder_s = time.perf_counter() - t0
    return HybridResult(
        radius=lad.radius,
        method="e1_ladder",
        oracle=oracle,
        exact=lad.exact,
        witness=(lad.witness_orientations and str(lad.witness_orientations)) or None,
        search_seconds=search_s,
        ladder_seconds=ladder_s,
        stats=stats,
    )


def describe(result: HybridResult) -> str:
    """A one-line, practitioner-readable summary of a result.

    Args:
        result: What :func:`breakdown_radius` returned.

    Returns:
        A sentence stating the radius, what it means operationally, and what it
        assumes.
    """
    if not result.exact:
        return (
            "The search did not complete within its budget, so no exact radius is "
            "available; the partial answer must not be read as a robustness claim."
        )
    if result.radius == UNREACHED:
        return (
            "No perturbation of the analyst's background knowledge makes this "
            "adjustment set invalid: it is robust throughout the reachable space "
            f"(established by {result.method} in {result.total_seconds:.3f}s, "
            f"assuming {result.assumes})."
        )
    if result.radius == 0:
        return "The adjustment set is already invalid at G0; the question is degenerate."
    return (
        f"The adjustment set survives any {result.radius - 1} of the analyst's "
        f"orientation claims being wrong, and fails at {result.radius} "
        f"(found by {result.method} in {result.total_seconds:.3f}s, "
        f"assuming {result.assumes})."
    )
