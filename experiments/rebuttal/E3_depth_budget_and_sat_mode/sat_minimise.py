"""E1 solved as ONE CP-SAT minimisation instead of a ladder of feasibility rungs.

The model is exactly the submitted E1 model (``bkrobust.sat.e1._build``): same
closure, same retraction-only restriction, same failure encoding. The only change
is the solve: minimise the number of retracted orientations in one call. OPTIMAL
with value r certifies that no failing state needs fewer than r retractions,
which is the statement the ladder's UNSAT rungs make; INFEASIBLE means no failing
state exists in the retraction-only space.
"""

from __future__ import annotations

import time

from ortools.sat.python import cp_model

from bkrobust.core.conventions import UNREACHED
from bkrobust.sat.e1 import _build


def radius_minimise(cpdag, g0, x, y, z, *, time_limit_s=300.0, seed=0, criterion="gac"):
    """Return a dict shaped like the fields of LadderResult that run.py records."""
    t0 = time.perf_counter()
    model, _ov, retracted, d = _build(cpdag, g0, x, y, z, criterion=criterion)
    build_s = time.perf_counter() - t0
    if retracted:
        model.Minimize(sum(retracted))

    solver = cp_model.CpSolver()
    # Same parameters as radius_e1, so the two modes differ only in the solve.
    solver.parameters.num_workers = 1
    solver.parameters.random_seed = seed
    solver.parameters.max_time_in_seconds = time_limit_s
    t = time.perf_counter()
    status = solver.Solve(model)
    solve_s = time.perf_counter() - t

    out = {
        "radius": UNREACHED,
        "exact": True,
        "status": solver.StatusName(status),
        "build_seconds": build_s,
        "solve_seconds": solve_s,
        "total_seconds": build_s + solve_s,
        "conflicts": solver.NumConflicts(),
        "branches": solver.NumBranches(),
        "witness_orientations": (),
    }
    if status == cp_model.OPTIMAL:
        out["radius"] = int(round(solver.ObjectiveValue())) if retracted else 0
        out["witness_orientations"] = tuple(
            sorted((u, v) for (u, v), var in d.items() if solver.Value(var))
        )
    elif status == cp_model.INFEASIBLE:
        pass
    else:
        # FEASIBLE (an upper bound only) or UNKNOWN: not a radius.
        out["exact"] = False
        out["radius"] = -2
    return out
