"""E3 -- assumption-free path unrolling. Does **not** assume Conjecture 2.

E1 searches only *retractions* of ``G0``. That is exact only if Conjecture 2
holds, and its error is one-sided: retraction-only can return a radius that is
too **large**, never too small. E3 exists to attack that from the other side.

**The encoding.** Unroll ``k + 1`` copies of the Meek-closure encoding,
``S_0 = K_{G0}, S_1, ..., S_k``, each an independently closed knowledge state of
the same CPDAG, and require of each consecutive pair that its symmetric
difference is **exactly one orientation**. The witness DAG and the failure
predicate are attached to the last copy only. The ladder returns the first ``k``
that is satisfiable.

**Why the step relation is exactly one orientation.** Two *sets* differing by a
single element have nothing strictly between them, so two closed states differing
by a single orientation are a covering pair *unconditionally* -- no appeal to
Lemma R, to anti-exchange, or to gradedness. Every E3 walk is therefore a genuine
walk in the covering graph, and a walk of length ``k`` to a failing state is a
constructive proof that the true BFS radius is at most ``k``.

**The asymmetry, which must not be blurred.** Covering steps that remove more
than one orientation at once are *not* representable here. (Lemma R says there
are none, but Lemma R rests on anti-exchange Case B, which is verified and not
proved -- and assuming it here would defeat the purpose.) So E3 searches a subset
of the covering walks and its radius is an **upper** bound on the true radius,
exactly as E1's is:

    r_true  <=  r_E3        and        r_true  <=  r_E1

* ``r_E3 < r_E1`` **refutes Conjecture 2**, soundly and with an explicit witness
  walk: a non-monotone path reaches failure sooner than any retraction can.
* ``r_E3 == r_E1`` is **consistent with** Conjecture 2 and is not a proof of it,
  because a shorter path through a multi-orientation cover would be invisible to
  both encodings.

E3 confirming E1 is evidence. Only E3 *disagreeing* is a theorem, and it is a
theorem in the negative direction. Do not report agreement as verification.

Unlike E1, E3 permits **down**-moves (adding an orientation), which is the move
Conjecture 2 asserts is never needed. That is the whole point of the encoding:
the search is free to descend and the conjecture predicts it never pays to.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

from ortools.sat.python import cp_model

from bkrobust.core.conventions import UNREACHED
from bkrobust.demo.graph import MPDAG
from bkrobust.sat.closure import Edge, OrientationVars, add_meek_closure, add_orientation_vars
from bkrobust.sat.failure import add_failure, add_witness_dag


@dataclass
class WalkResult:
    """An E3 radius together with its witness walk and its cost.

    Attributes:
        radius: The rung at which the unrolling first became satisfiable, or
            :data:`~bkrobust.core.conventions.UNREACHED` if no rung up to
            ``max_k`` was, or ``-2`` on solver timeout.
        exact: False if any rung returned UNKNOWN, and also false when the
            ladder was exhausted without a SAT answer -- unlike E1, whose ladder
            has a natural top (every orientation retracted), E3's ``max_k`` is a
            budget, so exhausting it proves only ``radius > max_k``. Read
            ``exhausted_to`` for that bound and do not read UNREACHED as "no
            failure exists".
        exhausted_to: The highest rung proved UNSAT when the ladder ran out,
            else ``None``. The anytime statement is ``radius > exhausted_to``.
        encoding: Always ``"E3"``.
        assumes: Always ``"nothing"`` -- this is the point of E3.
        walk: The witness walk as one frozenset of orientations per state,
            ``S_0 .. S_k``, ``S_0`` being ``G0``.
        rung_seconds: Solve seconds per rung.
        rung_status: ``SAT`` / ``UNSAT`` / ``UNKNOWN`` per rung.
        build_seconds: Total time spent building models.
        n_variables: Variables in the satisfied (or largest built) model.
        n_constraints: Constraints likewise.
        conflicts: Solver conflicts, summed over rungs (machine-independent).
        branches: Solver branches, summed over rungs (machine-independent).
    """

    radius: int
    exact: bool = True
    exhausted_to: int | None = None
    encoding: str = "E3"
    assumes: str = "nothing"
    walk: tuple[frozenset[Edge], ...] = ()
    rung_seconds: dict[int, float] = field(default_factory=dict)
    rung_status: dict[int, str] = field(default_factory=dict)
    build_seconds: float = 0.0
    n_variables: int = 0
    n_constraints: int = 0
    conflicts: int = 0
    branches: int = 0

    @property
    def sat_seconds(self) -> float:
        """Time on the one satisfiable rung, if any."""
        return sum(t for k, t in self.rung_seconds.items() if self.rung_status.get(k) == "SAT")

    @property
    def unsat_seconds(self) -> float:
        """Time proving the lower rungs unsatisfiable -- usually the dominant cost."""
        return sum(t for k, t in self.rung_seconds.items() if self.rung_status.get(k) == "UNSAT")

    @property
    def total_seconds(self) -> float:
        """Build plus all rungs."""
        return self.build_seconds + sum(self.rung_seconds.values())


def _pin(model: cp_model.CpModel, ov: OrientationVars, g0: MPDAG) -> None:
    """Pin a copy to ``G0`` exactly."""
    for a, b in ov.free:
        for u, v in ((a, b), (b, a)):
            model.Add(ov.d[(u, v)] == (1 if g0.is_directed_edge(u, v) else 0))


def _one_orientation_apart(
    model: cp_model.CpModel, lo: OrientationVars, hi: OrientationVars
) -> None:
    """Symmetric difference of the two knowledge sets is exactly one orientation.

    Only the free edges can differ; the compelled ones are pinned identically in
    both copies by :func:`add_orientation_vars`.
    """
    diffs = []
    for a, b in lo.free:
        for u, v in ((a, b), (b, a)):
            p, q = lo.d[(u, v)], hi.d[(u, v)]
            t = model.NewBoolVar(f"x_{u}_{v}")
            # t <-> (p != q)
            model.AddBoolOr([p, q, t.Not()])
            model.AddBoolOr([p.Not(), q.Not(), t.Not()])
            model.AddBoolOr([p.Not(), q, t])
            model.AddBoolOr([p, q.Not(), t])
            diffs.append(t)
    model.Add(sum(diffs) == 1)


def _build(
    cpdag: MPDAG, g0: MPDAG, x: str, y: str, z: frozenset[str], k: int
) -> tuple[cp_model.CpModel, list[OrientationVars]]:
    """Build the ``k``-step unrolling. Returns the model and the per-state vars."""
    model = cp_model.CpModel()
    states = []
    for _ in range(k + 1):
        ov = add_orientation_vars(model, cpdag)
        add_meek_closure(model, ov)
        states.append(ov)
    _pin(model, states[0], g0)
    for i in range(k):
        _one_orientation_apart(model, states[i], states[i + 1])
    last = states[-1]
    e = add_witness_dag(model, last)
    add_failure(model, last, e, x, y, z)
    return model, states


def radius_e3(
    cpdag: MPDAG,
    g0: MPDAG,
    x: str,
    y: str,
    z: frozenset[str],
    *,
    max_k: int = 6,
    time_limit_s: float = 60.0,
    seed: int = 0,
) -> WalkResult:
    """Compute an assumption-free upper bound on the radius by unrolling the walk.

    Args:
        cpdag: The estimated CPDAG.
        g0: The analyst's graph.
        x: Treatment.
        y: Outcome.
        z: The adjustment set, fixed once from ``g0``.
        max_k: Highest rung to unroll. Each rung costs one more copy of the
            closure encoding, so this is the real budget knob.
        time_limit_s: Per-rung solver time limit.
        seed: Solver seed; with a single worker this makes the search
            deterministic.

    Returns:
        A :class:`WalkResult`. Compare its radius with E1's: strictly smaller
        refutes Conjecture 2, equal is consistent with it and proves nothing.
    """
    res = WalkResult(radius=UNREACHED)
    for k in range(0, max_k + 1):
        t0 = time.perf_counter()
        model, states = _build(cpdag, g0, x, y, z, k)
        res.build_seconds += time.perf_counter() - t0
        proto = model.Proto()
        res.n_variables = len(proto.variables)
        res.n_constraints = len(proto.constraints)

        solver = cp_model.CpSolver()
        solver.parameters.num_workers = 1
        solver.parameters.random_seed = seed
        solver.parameters.max_time_in_seconds = time_limit_s
        t = time.perf_counter()
        status = solver.Solve(model)
        res.rung_seconds[k] = time.perf_counter() - t
        res.conflicts += solver.NumConflicts()
        res.branches += solver.NumBranches()

        if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            res.rung_status[k] = "SAT"
            res.walk = tuple(
                frozenset(uv for uv, var in sorted(ov.d.items()) if solver.Value(var))
                for ov in states
            )
            res.radius = k
            return res
        if status == cp_model.INFEASIBLE:
            res.rung_status[k] = "UNSAT"
            continue
        res.rung_status[k] = "UNKNOWN"
        res.exact = False
        res.radius = -2
        return res
    # Every rung UNSAT: an anytime lower bound, not a proof of unreachability.
    res.exhausted_to = max_k
    res.exact = False
    return res
