"""Phase 3: breakdown radii on real network structure.

Applies the pre-registered selection policy mechanically
(``results/axisa3/preregistration.md`` §2). Nothing here chooses a pair, a
network or a sample on the basis of any computed radius.

Admissibility, all four conditions required:

1. ``Y`` is a descendant of ``X`` in the true DAG — a causal path exists;
2. ``X`` lies in, or is adjacent to, an undirected component of the CPDAG;
3. the instance survives the degeneracy gate;
4. ``O(G0)`` is identified and genuinely valid at ``G0``.

The gate is the cheap equivalent introduced in session 5, which replaces the
three checks that call ``all_valid_adjustment_sets_mpdag`` — an enumeration over
every subset of V, which at these sizes would not terminate — with equivalents
costing one back-door test each. It is **strictly stricter** than
``synth.runner.gate`` on the perturbation check (16 disagreements in 46,800
cases, 0.034%, all one-directional), and the instances it drops are ones where
``O`` is robust to every atomic perturbation. That biases against finding large
radii, which is the conservative direction here.

Sentinels, kept distinct because sessions 3, 4 and 5 each shipped a bug in this
family: an undefined separation is ``None`` with a status string and is never a
number; ``UNREACHED`` is a status, not a radius; a censored run carries
``wall_until_timeout_s`` and never the key a measurement uses.
"""

from __future__ import annotations

import time
from dataclasses import asdict, dataclass, field
from typing import Any

from bkrobust.core.conventions import UNREACHED
from bkrobust.demo.evaluate import (
    is_valid_adjustment_set_dag,
    is_valid_adjustment_set_mpdag,
    optimal_adjustment_set_dag,
    optimal_adjustment_set_mpdag,
)
from bkrobust.demo.example import knowledge_to_recover
from bkrobust.demo.graph import MPDAG, undirected_components
from bkrobust.demo.meek import apply_orientations
from bkrobust.gac import is_gac_valid_mpdag
from bkrobust.hybrid import breakdown_radius

Edge = tuple[str, str]


def fast_gate(dag: MPDAG, cpdag: MPDAG, x: str, y: str) -> tuple[bool, str]:
    """The degeneracy gate, without the subset enumeration that cannot scale.

    Same verdict vocabulary as :func:`bkrobust.synth.runner.gate`.

    Args:
        dag: The ground-truth DAG.
        cpdag: Its CPDAG.
        x: Treatment.
        y: Outcome.

    Returns:
        ``(True, "ok")`` or ``(False, reason)``.
    """
    comps = undirected_components(cpdag)
    if not any(x in c or any(cpdag.has_edge(x, v) for v in sorted(c)) for c in comps):
        return False, "treatment_not_in_or_adjacent_to_component"
    if y not in dag.descendants(x):
        return False, "no_causal_path"
    o = frozenset(optimal_adjustment_set_dag(dag, x, y))
    if not is_valid_adjustment_set_dag(dag, x, y, o):
        return False, "no_valid_adjustment_set"
    if is_valid_adjustment_set_dag(dag, x, y, frozenset()):
        return False, "empty_set_trivially_valid"
    k_true = sorted(knowledge_to_recover(dag, cpdag))
    for drop in k_true:
        g0 = apply_orientations(cpdag, [e for e in k_true if e != drop])
        if g0 is None:
            continue
        if not is_valid_adjustment_set_mpdag(g0, x, y, o):
            return True, "ok"
    return False, "no_atomic_perturbation_changes_validity"


def select_knowledge(dag: MPDAG, cpdag: MPDAG, coverage: float) -> list[Edge]:
    """The analyst's asserted orientations, at a given coverage.

    Deterministic: the true orientations of the CPDAG's undirected edges, sorted,
    with an evenly spaced ``coverage`` fraction retained. No RNG, so the same
    network and coverage always give the same ``K``.

    Args:
        dag: The ground-truth DAG, read for the true orientations.
        cpdag: Its CPDAG, read for which edges are undirected.
        coverage: Fraction of undirected edges the analyst asserts, in ``[0, 1]``.

    Returns:
        The asserted orientations as ``(tail, head)`` pairs.
    """
    k_true = sorted(knowledge_to_recover(dag, cpdag))
    if coverage >= 1.0:
        return k_true
    keep = round(len(k_true) * coverage)
    if keep <= 0:
        return []
    step = len(k_true) / keep
    return [k_true[min(len(k_true) - 1, int(i * step))] for i in range(keep)]


def component_of(graph: MPDAG, node: str) -> frozenset[str] | None:
    """The undirected component containing ``node``, or None if it is in none."""
    for comp in undirected_components(graph):
        if node in comp:
            return comp
    return None


def separation(cpdag: MPDAG, x: str, z: frozenset[str]) -> tuple[int | None, str]:
    """Distance inside ``X``'s component to the nearest member of ``Z``.

    Returns:
        ``(distance, "measured")``, or ``(None, reason)``. The distance is
        **never** encoded as a number when it is undefined.
    """
    comp = component_of(cpdag, x)
    if comp is None:
        return None, "x_not_in_component"
    dist = {x: 0}
    frontier = [x]
    while frontier:
        nxt = []
        for u in sorted(frontier):
            for v in sorted(cpdag.neighbors(u)):
                if v in comp and v not in dist:
                    dist[v] = dist[u] + 1
                    nxt.append(v)
        frontier = nxt
    reach = [dist[m] for m in sorted(z) if m in dist]
    if not reach:
        return None, "no_z_member_in_component"
    return min(reach), "measured"


@dataclass
class InstanceResult:
    """One admissible instance, or one rejection, with its provenance."""

    network: str
    x: str
    y: str
    coverage: float
    admissible: bool
    reject_reason: str = ""
    component_size: int = 0
    largest_component_size: int = 0
    separation: int | None = None
    separation_status: str = ""
    k_g0: int = 0
    z_size: int = 0
    radius: int = UNREACHED
    method: str = ""
    oracle: str = ""
    exact: bool = True
    seconds: float = 0.0
    stats: dict[str, int] = field(default_factory=dict)

    def as_row(self) -> dict[str, Any]:
        """Flatten for a JSONL row."""
        return asdict(self)


def evaluate(
    network: str,
    dag: MPDAG,
    cpdag: MPDAG,
    x: str,
    y: str,
    coverage: float,
    *,
    time_limit_s: float = 300.0,
) -> InstanceResult:
    """Screen one ``(X, Y)`` pair and, if admissible, compute its radius.

    Args:
        network: Network name, for the row.
        dag: The ground-truth DAG.
        cpdag: Its CPDAG.
        x: Treatment.
        y: Outcome.
        coverage: Knowledge coverage level.
        time_limit_s: Ladder time limit.

    Returns:
        An :class:`InstanceResult`, admissible or not.
    """
    comps = undirected_components(cpdag)
    largest = max((len(c) for c in comps), default=0)
    res = InstanceResult(
        network=network,
        x=x,
        y=y,
        coverage=coverage,
        admissible=False,
        largest_component_size=largest,
    )
    ok, reason = fast_gate(dag, cpdag, x, y)
    if not ok:
        res.reject_reason = reason
        return res

    k = select_knowledge(dag, cpdag, coverage)
    g0 = apply_orientations(cpdag, k)
    if g0 is None:
        res.reject_reason = "knowledge_inconsistent"
        return res
    o = optimal_adjustment_set_mpdag(g0, x, y)
    if not o:
        res.reject_reason = "o_g0_not_identified"
        return res
    z = frozenset(o)
    if not is_gac_valid_mpdag(g0, x, y, z):
        res.reject_reason = "z_invalid_at_g0"
        return res

    comp = component_of(cpdag, x)
    sep, status = separation(cpdag, x, z)
    res.admissible = True
    res.reject_reason = "ok"
    res.component_size = len(comp) if comp else 0
    res.separation, res.separation_status = sep, status
    res.k_g0 = sum(
        1 for (a, b) in g0.directed_edges if tuple(sorted((a, b))) in cpdag.undirected_edges
    )
    res.z_size = len(z)

    t = time.perf_counter()
    out = breakdown_radius(cpdag, None, x, y, z, g0=g0, time_limit_s=time_limit_s)
    res.seconds = round(time.perf_counter() - t, 5)
    res.radius, res.method, res.oracle, res.exact = (out.radius, out.method, out.oracle, out.exact)
    res.stats = out.stats.as_dict()
    return res
