r"""E2 -- join distance. Assumes the up-then-down normalisation, not Conjecture 2.

E1 restricts the failing state to a *retraction* of ``G0``; E2 lifts that
restriction and lets the failing state be any closed knowledge state ``G``,
measuring distance through the join:

    d(G0, G)  =  (|K_G0| - |K_J|) + (|K_G| - |K_J|),      K_J = cl(K_G0 ∩ K_G).

**The join needs no closure.** Theorem 2 identifies the join ``G0 v G`` with
``K_G0 ∩ K_G``, and the intersection of two closed sets is closed
(``cl(A∩B) ⊆ cl(A) ∩ cl(B) = A ∩ B``). So ``K_J`` *is* ``K_G0 ∩ K_G``, and with
``K_G0`` a constant the whole objective collapses to a symmetric difference:

    d(G0, G)  =  |K_G0 \\ K_G| + |K_G \\ K_G0|  =  |K_G0 Δ K_G|,

which is one linear constraint over a single copy of the orientation variables.
The brief budgets an extra copy for ``J``; it is not needed, because ``J`` is
determined rather than searched for. The saving is real but the assumptions are
unchanged, and they are what matters:

* **Theorem 2** (proved) for the join;
* **Lemma R** / gradedness for reading a rank difference off a set difference --
  which rests on anti-exchange Case B, verified and not proved;
* the **up-then-down normalisation**: that a shortest path may be taken up to the
  join and then down.

E2 is therefore weaker in assumption than E1 (it does not assume the down-leg is
empty, only that it comes second) and wider in reach: it searches every closed
state, not just the retractions. Because it admits ``G`` outside the up-set,
``r_E2 <= r_E1`` always. Like E1 and E3 it is an upper bound on the true radius.

A strict gap ``r_E2 < r_E1`` says the nearest failure is *not* a retraction of
``G0``, which contradicts Conjecture 2 -- but only under E2's own assumptions,
so unlike an E3 gap it is not a proof. Adjudicate any such gap with E3, which
assumes nothing.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any

from ortools.sat.python import cp_model

from bkrobust.core.conventions import UNREACHED
from bkrobust.demo.graph import MPDAG
from bkrobust.sat.closure import Edge, add_meek_closure, add_orientation_vars
from bkrobust.sat.failure import add_failure, add_witness_dag


@dataclass
class JoinResult:
    r"""An E2 radius with its witness and cost.

    Attributes:
        radius: The rung at which the ladder first became satisfiable, or
            :data:`~bkrobust.core.conventions.UNREACHED`, or ``-2`` on timeout.
        exact: False if any rung returned UNKNOWN.
        encoding: Always ``"E2"``.
        assumes: The assumption chain, carried on the result so it cannot be
            dropped when the number is copied into a table.
        witness_orientations: ``K_G`` for the failing state found.
        up_steps: ``|K_G0 \\ K_G|`` -- the retraction leg.
        down_steps: ``|K_G \\ K_G0|`` -- the orientation leg. Conjecture 2
            predicts this is always 0 at the optimum.
        rung_seconds: Solve seconds per rung.
        rung_status: ``SAT`` / ``UNSAT`` / ``UNKNOWN`` per rung.
        build_seconds: Model build time.
        n_variables: Variables in the model.
        n_constraints: Constraints in the model.
        conflicts: Solver conflicts summed over rungs.
        branches: Solver branches summed over rungs.
    """

    radius: int
    exact: bool = True
    encoding: str = "E2"
    assumes: str = "Theorem 2 + Lemma R (gradedness) + up-then-down normalisation"
    witness_orientations: tuple[Edge, ...] = ()
    up_steps: int = 0
    down_steps: int = 0
    rung_seconds: dict[int, float] = field(default_factory=dict)
    rung_status: dict[int, str] = field(default_factory=dict)
    build_seconds: float = 0.0
    n_variables: int = 0
    n_constraints: int = 0
    conflicts: int = 0
    branches: int = 0

    @property
    def sat_seconds(self) -> float:
        """Time on the satisfiable rung."""
        return sum(t for k, t in self.rung_seconds.items() if self.rung_status.get(k) == "SAT")

    @property
    def unsat_seconds(self) -> float:
        """Time spent certifying the clean shells."""
        return sum(t for k, t in self.rung_seconds.items() if self.rung_status.get(k) == "UNSAT")

    @property
    def total_seconds(self) -> float:
        """Build plus all rungs."""
        return self.build_seconds + sum(self.rung_seconds.values())


def _build(
    cpdag: MPDAG, g0: MPDAG, x: str, y: str, z: frozenset[str]
) -> tuple[cp_model.CpModel, list[tuple[Any, bool]], dict[Edge, Any]]:
    """Build the shared E2 model.

    Returns:
        ``(model, sym, d)`` where ``sym`` pairs each free orientation variable
        with whether it enters the symmetric difference *negated*. The sign is
        carried separately rather than as a ``Not()`` literal because rungs are
        added to per-rung clones, and a negated literal has no proto index to
        look up on the clone.
    """
    model = cp_model.CpModel()
    ov = add_orientation_vars(model, cpdag)
    add_meek_closure(model, ov)

    # |K_G0 Δ K_G|, one literal per free orientation. K_G0 is constant, so each
    # term is either the variable or its negation -- no reification needed.
    sym: list[tuple[Any, bool]] = []
    for a, b in ov.free:
        for u, v in ((a, b), (b, a)):
            sym.append((ov.d[(u, v)], g0.is_directed_edge(u, v)))

    e = add_witness_dag(model, ov)
    add_failure(model, ov, e, x, y, z)
    return model, sym, ov.d


def radius_e2(
    cpdag: MPDAG,
    g0: MPDAG,
    x: str,
    y: str,
    z: frozenset[str],
    *,
    max_k: int | None = None,
    time_limit_s: float = 60.0,
    seed: int = 0,
) -> JoinResult:
    """Compute the radius by the E2 join-distance ladder.

    Args:
        cpdag: The estimated CPDAG.
        g0: The analyst's graph.
        x: Treatment.
        y: Outcome.
        z: The adjustment set, fixed once from ``g0``.
        max_k: Highest rung to try; defaults to the number of free orientations.
        time_limit_s: Per-rung solver time limit.
        seed: Solver seed; with a single worker this makes the search
            deterministic.

    Returns:
        A :class:`JoinResult`. ``down_steps > 0`` at the optimum would mean the
        nearest failure is not a retraction -- report it, and adjudicate with E3.
    """
    t0 = time.perf_counter()
    model, sym, d = _build(cpdag, g0, x, y, z)
    build_s = time.perf_counter() - t0

    res = JoinResult(radius=UNREACHED, build_seconds=build_s)
    proto = model.Proto()
    res.n_variables = len(proto.variables)
    res.n_constraints = len(proto.constraints)

    top = len(sym) if max_k is None else min(max_k, len(sym))
    spec = [(var.Index(), neg) for var, neg in sym]
    for k in range(0, top + 1):
        sub = model.Clone()
        if spec:
            terms = []
            for i, neg in spec:
                v = sub.GetBoolVarFromProtoIndex(i)
                terms.append(1 - v if neg else v)
            sub.Add(sum(terms) <= k)
        solver = cp_model.CpSolver()
        solver.parameters.num_workers = 1
        solver.parameters.random_seed = seed
        solver.parameters.max_time_in_seconds = time_limit_s
        t = time.perf_counter()
        status = solver.Solve(sub)
        res.rung_seconds[k] = time.perf_counter() - t
        res.conflicts += solver.NumConflicts()
        res.branches += solver.NumBranches()

        if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            res.rung_status[k] = "SAT"
            got = [
                uv
                for uv, var in sorted(d.items())
                if solver.Value(sub.GetBoolVarFromProtoIndex(var.Index()))
            ]
            res.witness_orientations = tuple(got)
            kg = set(got)
            k0 = set(g0.directed_edges)
            res.up_steps = len(k0 - kg)
            res.down_steps = len(kg - k0)
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
