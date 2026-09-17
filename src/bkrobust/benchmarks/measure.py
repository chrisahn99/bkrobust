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

import numpy as np

from bkrobust.core.conventions import UNREACHED
from bkrobust.demo.evaluate import (
    is_valid_adjustment_set_dag,
    optimal_adjustment_set_dag,
    optimal_adjustment_set_mpdag,
)
from bkrobust.demo.example import knowledge_to_recover
from bkrobust.demo.graph import MPDAG, undirected_components
from bkrobust.demo.meek import apply_orientations
from bkrobust.gac import is_gac_valid_mpdag
from bkrobust.hybrid import DEFAULT_SEARCH_BUDGET, breakdown_radius
from bkrobust.mpdag_criterion.criterion import is_amenable

Edge = tuple[str, str]

#: ``optimal_adjustment_set_mpdag`` enumerates the DAG extensions of ``G0``, which
#: is exponential in the number of undirected edges ``G0`` still carries. At low
#: knowledge coverage on a large real component that does not terminate. Matches
#: :data:`bkrobust.synth.component_generator.MAX_G0_UNDIRECTED_FOR_EXTENSIONS`.
#:
#: Exceeding it is a **measurement limit, not a structural rejection**, and gets
#: its own status so the two can never be confused in the analysis: a pair that
#: is genuinely degenerate and a pair we could not afford to evaluate are
#: different facts.
MAX_G0_UNDIRECTED_FOR_EXTENSIONS: int = 12
O_INTRACTABLE: str = "o_g0_extensions_intractable"


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
        # The polynomial GAC predicate, not the extension-enumerating one. By
        # Theorem 14 these coincide exactly on ``Z = O(...)``, which is what ``o``
        # is, so this is a cost substitution and not a change of predicate.
        # Verified rather than assumed: 2,240 pairs across child, insurance, asia,
        # sachs and water, identical verdicts, 8.2x faster
        # (results/axisa3/gate_predicate_swap.json). The enumerating version is
        # what makes the gate the bottleneck at these sizes.
        if not is_gac_valid_mpdag(g0, x, y, o):
            return True, "ok"
    return False, "no_atomic_perturbation_changes_validity"


def select_knowledge(
    dag: MPDAG,
    cpdag: MPDAG,
    coverage: float,
    *,
    mode: str = "stride",
    seed: int = 20260912,
) -> list[Edge]:
    """The analyst's asserted orientations, at a given coverage.

    Deterministic under both modes: no global RNG, so a network and a coverage
    always give the same ``K``.

    ``stride`` is the original rule and stays the default so every committed
    result reproduces: sort the recovering set and keep an evenly spaced
    ``coverage`` fraction of it.

    ``nested`` keeps a prefix of one seeded permutation instead. The stride is
    **not nested** -- the indices it keeps at 0.25 are not a subset of those it
    keeps at 0.5 -- so a coverage sweep under ``stride`` is not a retraction
    sequence, and on this corpus that fails on four of the twenty-five
    instance-yielding networks (diabetes, ecoli70, magic-irri, munin1). Use
    ``nested`` for anything that reads a coverage sweep as retraction. The two
    modes agree at coverage 1.0.

    Args:
        dag: The ground-truth DAG, read for the true orientations.
        cpdag: Its CPDAG, read for which edges are undirected.
        coverage: Fraction of undirected edges the analyst asserts, in ``[0, 1]``.
        mode: ``"stride"`` (committed behaviour) or ``"nested"``.
        seed: Permutation seed, used by ``nested`` only.

    Returns:
        The asserted orientations as ``(tail, head)`` pairs.

    Raises:
        ValueError: If ``mode`` is neither ``"stride"`` nor ``"nested"``.
    """
    if mode not in ("stride", "nested"):
        raise ValueError(f"unknown selection mode {mode!r}")
    k_true = sorted(knowledge_to_recover(dag, cpdag))
    if coverage >= 1.0:
        return k_true
    keep = round(len(k_true) * coverage)
    if keep <= 0:
        return []
    if mode == "nested":
        order = np.random.default_rng(seed).permutation(len(k_true))
        return [k_true[int(i)] for i in sorted(order[:keep])]
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
    g0_undirected_edges: int = 0
    #: Asserted claims, as distinct from ``k_g0``, which counts the Meek closure.
    n_k: int = 0
    #: ``identified_nonempty`` / ``identified_empty`` / ``not_identified`` / ``""``.
    o_verdict: str = ""
    #: Whether the query is amenable at ``G0``; ``None`` where it was not reached.
    amenable: bool | None = None
    #: Which leg of the hybrid answered, printed so the estimator and the
    #: answer are not confounded in the analysis.
    dispatch_leg: str = ""
    search_budget: int = 0
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
    selection_mode: str = "stride",
    search_budget: int | None = None,
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
        selection_mode: Passed to :func:`select_knowledge`. ``"stride"``
            reproduces every committed row; ``"nested"`` makes the coverage
            sweep a retraction sequence.
        search_budget: Depth budget for the bounded search before the ladder
            takes over. ``None`` keeps the library default and is what every
            committed row used.

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

    k = select_knowledge(dag, cpdag, coverage, mode=selection_mode)
    res.n_k = len(k)
    g0 = apply_orientations(cpdag, k)
    if g0 is None:
        res.reject_reason = "knowledge_inconsistent"
        return res
    # Assigned on every branch from here down. It used to be written only in the
    # intractable branch, so on every other row the value was the dataclass
    # default and could be read as a measurement of zero.
    res.g0_undirected_edges = len(g0.undirected_edges)
    res.k_g0 = sum(
        1 for (a, b) in g0.directed_edges if tuple(sorted((a, b))) in cpdag.undirected_edges
    )
    if len(g0.undirected_edges) > MAX_G0_UNDIRECTED_FOR_EXTENSIONS:
        res.reject_reason = O_INTRACTABLE
        return res
    o = optimal_adjustment_set_mpdag(g0, x, y)
    # Three-valued, because ``if not o`` cannot tell "the extensions disagree"
    # (None) from "the optimal set is empty" (frozenset()), and the two are
    # different facts about the query.
    if o is None:
        res.o_verdict = "not_identified"
        # Whether the effect is identifiable by adjustment at all is the
        # substantive question here, and it is what the old single label
        # ``o_g0_not_identified`` was silently reporting.
        res.amenable = is_amenable(g0, x, y)
        res.reject_reason = "not_amenable" if not res.amenable else "o_g0_not_identified"
        return res
    if not o:
        res.o_verdict = "identified_empty"
        res.reject_reason = "o_g0_empty"
        return res
    res.o_verdict = "identified_nonempty"
    res.amenable = True
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
    res.z_size = len(z)

    t = time.perf_counter()
    kwargs = {} if search_budget is None else {"search_budget": search_budget}
    out = breakdown_radius(cpdag, None, x, y, z, g0=g0, time_limit_s=time_limit_s, **kwargs)
    res.seconds = round(time.perf_counter() - t, 5)
    res.radius, res.method, res.oracle, res.exact = (out.radius, out.method, out.oracle, out.exact)
    # ``method`` is a deterministic function of the answer: the bounded search
    # answers everything inside its budget and the ladder answers the rest, so
    # the estimator and the radius are confounded unless the leg is printed.
    res.dispatch_leg = out.method
    res.search_budget = search_budget if search_budget is not None else DEFAULT_SEARCH_BUDGET
    res.stats = out.stats.as_dict()
    return res
