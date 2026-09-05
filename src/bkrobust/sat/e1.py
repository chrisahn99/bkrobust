r"""E1 -- the retraction-only encoding, solved as an incremental ladder.

By Conjecture 2 (proved in session 2 modulo Anti-Exchange Case B),

    ``r_val = min |K_G0 \ S|`` over closed ``S ⊆ K_G0`` with ``Z`` failing at
    ``Meek(Ĉ, S)``.

so the search ranges over subsets of the analyst's own orientations. The model is
the closure encoding restricted to orientations ``G0`` already has, plus a
failing witness DAG.

**Assumption status: E1 is exact iff Conjecture 2 holds.** Every radius it
returns is conditional on that, and every result row is tagged accordingly. If
Conjecture 2 were false, E1 would return radii that are too **large** -- it
searches only the up-set, so it can miss a nearer incomparable failure. It cannot
return one that is too small.

Why a ladder rather than one optimisation
-----------------------------------------

Solve "is there a failing state within ``k`` retractions?" for ``k = 0, 1, 2, …``
and stop at the first SAT. Three reasons, all borne out by the previous
sessions' data:

* radii are small -- about 80% of instances have ``r_val = 1`` -- so the ladder
  almost always terminates on the first or second rung;
* every UNSAT rung is a **certificate that shells 0…k are clean**, which is the
  statement a practitioner actually wants, and it is produced without ever
  enumerating those shells;
* when a budget runs out the partial result is still meaningful: ``r ≥ k+1``.

SAT and UNSAT times are recorded separately because they behave differently, and
the UNSAT side is what proves robustness.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any

from ortools.sat.python import cp_model

from bkrobust.core.conventions import UNREACHED
from bkrobust.demo.graph import MPDAG, canon
from bkrobust.sat.closure import add_meek_closure, add_orientation_vars
from bkrobust.sat.failure import add_failure, add_witness_dag

Edge = tuple[str, str]


@dataclass
class LadderResult:
    """Outcome of an E1 ladder solve.

    Attributes:
        radius: The radius, or ``UNREACHED`` if no failing state exists within
            the retraction budget, or ``-2`` if the budget was exhausted first.
        exact: Whether the ladder ran to a definite answer.
        encoding: Always ``"E1"``; carried so result rows are self-describing.
        assumes: The assumption the number depends on.
        witness_orientations: The failing MPDAG's orientations; empty if no
            rung was satisfiable.
        rung_seconds: Wall-clock per rung, keyed by ``k``.
        rung_status: ``"SAT"`` / ``"UNSAT"`` / ``"UNKNOWN"`` per rung.
        build_seconds: Model construction time, which is often the dominant cost
            and is routinely overlooked.
        n_variables: Model size.
        n_constraints: Model size.
        conflicts: Solver conflicts summed over rungs.
        branches: Solver branches summed over rungs.
    """

    radius: int
    exact: bool
    encoding: str = "E1"
    assumes: str = "Conjecture 2"
    witness_orientations: tuple[Edge, ...] = ()
    rung_seconds: dict[int, float] = field(default_factory=dict)
    rung_status: dict[int, str] = field(default_factory=dict)
    build_seconds: float = 0.0
    n_variables: int = 0
    n_constraints: int = 0
    conflicts: int = 0
    branches: int = 0

    @property
    def sat_seconds(self) -> float:
        """Time spent on the single satisfiable rung, if any."""
        return sum(t for k, t in self.rung_seconds.items() if self.rung_status.get(k) == "SAT")

    @property
    def unsat_seconds(self) -> float:
        """Time spent proving the lower rungs clean -- the robustness certificate."""
        return sum(t for k, t in self.rung_seconds.items() if self.rung_status.get(k) == "UNSAT")

    @property
    def total_seconds(self) -> float:
        """Build plus all rungs."""
        return self.build_seconds + sum(self.rung_seconds.values())


def _build(
    cpdag: MPDAG,
    g0: MPDAG,
    x: str,
    y: str,
    z: frozenset[str],
) -> tuple[cp_model.CpModel, Any, list[Any], dict[Edge, Any]]:
    """Build the shared E1 model; returns (model, ov, retracted_literals, d)."""
    model = cp_model.CpModel()
    ov = add_orientation_vars(model, cpdag)
    add_meek_closure(model, ov)

    # retraction-only: G may only carry orientations G0 already has
    retracted: list[Any] = []
    for a, b in ov.free:
        for u, v in ((a, b), (b, a)):
            if not g0.is_directed_edge(u, v):
                model.Add(ov.d[(u, v)] == 0)
        if g0.is_directed_edge(a, b) or g0.is_directed_edge(b, a):
            u, v = (a, b) if g0.is_directed_edge(a, b) else (b, a)
            keep = ov.d[(u, v)]
            drop = model.NewBoolVar(f"drop_{u}_{v}")
            model.Add(drop == 1 - keep)
            retracted.append(drop)

    e = add_witness_dag(model, ov)
    add_failure(model, ov, e, x, y, z)
    return model, ov, retracted, ov.d


def radius_e1(
    cpdag: MPDAG,
    g0: MPDAG,
    x: str,
    y: str,
    z: frozenset[str],
    *,
    max_k: int | None = None,
    time_limit_s: float = 60.0,
    seed: int = 0,
) -> LadderResult:
    """Compute the breakdown radius by the E1 ladder.

    Args:
        cpdag: The estimated CPDAG.
        g0: The analyst's graph, ``Meek(cpdag, K)``.
        x: Treatment.
        y: Outcome.
        z: The adjustment set, fixed once from ``g0``.
        max_k: Highest rung to try; defaults to the number of retractable
            orientations, beyond which no further retraction exists.
        time_limit_s: Per-rung solver time limit.
        seed: Solver seed. Combined with single-threaded operation this makes
            the search deterministic.

    Returns:
        A :class:`LadderResult`.
    """
    t0 = time.perf_counter()
    model, _ov, retracted, d = _build(cpdag, g0, x, y, z)
    build_s = time.perf_counter() - t0

    res = LadderResult(radius=UNREACHED, exact=True, build_seconds=build_s)
    proto = model.Proto()
    res.n_variables = len(proto.variables)
    res.n_constraints = len(proto.constraints)

    top = len(retracted) if max_k is None else min(max_k, len(retracted))
    for k in range(0, top + 1):
        sub = model.Clone()
        # Rebuild the budget on the clone by index, since literals are per-model.
        idx = [r.Index() for r in retracted]
        if idx:
            sub.Add(sum(sub.GetBoolVarFromProtoIndex(i) for i in idx) <= k)
        solver = cp_model.CpSolver()
        solver.parameters.num_workers = 1
        solver.parameters.random_seed = seed
        solver.parameters.max_time_in_seconds = time_limit_s
        t = time.perf_counter()
        status = solver.Solve(sub)
        dt = time.perf_counter() - t
        res.rung_seconds[k] = dt
        res.conflicts += solver.NumConflicts()
        res.branches += solver.NumBranches()
        if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            res.rung_status[k] = "SAT"
            got = []
            for (u, v), var in sorted(d.items()):
                if solver.Value(sub.GetBoolVarFromProtoIndex(var.Index())):
                    got.append((u, v))
            res.witness_orientations = tuple(got)
            res.radius = k
            return res
        if status == cp_model.INFEASIBLE:
            res.rung_status[k] = "UNSAT"
            continue
        res.rung_status[k] = "UNKNOWN"
        res.exact = False
        res.radius = -2
        return res
    return res


def certified_clean_through(res: LadderResult) -> int:
    """The largest ``k`` for which the ladder proved every shell 0…k clean."""
    unsat = [k for k, s in res.rung_status.items() if s == "UNSAT"]
    return max(unsat) if unsat else -1


def canon_edges(edges: tuple[Edge, ...]) -> frozenset[Edge]:
    """Canonicalise a witness's orientation set for comparison."""
    return frozenset(canon(a, b) for a, b in edges)
