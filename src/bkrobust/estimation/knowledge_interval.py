"""The knowledge interval: effects left possible by every graph in a retraction ball.

An analyst's claims ``K`` orient a CPDAG ``C`` into ``G0 = Meek(C, K)``. The
retraction ball ``B_r(K)`` holds the Meek re-closure of ``K`` minus ``S`` for every
``|S| <= r``. Each graph in the ball leaves a set of total effects of ``X`` on
``Y`` possible, and those are computed IDA-style: enumerate the parent sets of
``X`` that some DAG extension of the graph realises, and regress ``Y`` on ``X``
and that parent set. The knowledge interval ``I_r`` is the hull of those effects
over the ball; with a sample covariance it is the hull of the per-effect
confidence intervals.

Three facts keep this cheap, and each is checked by
``tests/estimation/test_knowledge_interval.py`` against brute-force
enumeration of DAG extensions rather than assumed:

* **Only the treatment's chain component matters.** In a CPDAG every node of a
  chain component has the same parents outside it, so Meek's rules fired by a
  claim inside the component never reach outside it and nothing outside ever
  fires a rule inside it. The closure of ``C`` with ``K`` restricted to the
  component equals the full closure restricted to it, and the parent sets of
  ``X`` are its fixed outside parents plus a set read off the component alone.
* **Semi-local enumeration.** A set ``P`` of undirected neighbours of ``X`` is a
  possible parent set iff orienting ``P -> X`` and ``X`` into every other neighbour is
  accepted by :func:`~bkrobust.demo.meek.apply_orientations`. Only cliques are
  tried, since two non-adjacent parents would form a new v-structure.
* **Claims outside the component leave the interval unchanged**, so the ball is
  enumerated over the component's claims only.

The population coefficient of ``X`` for a parent set that contains ``Y`` is zero
by the IDA convention (``Y`` is then a parent, not a descendant, of ``X``).
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import numpy as np

from bkrobust.demo.graph import MPDAG, canon, induced_subgraph, undirected_components
from bkrobust.demo.meek import apply_orientations
from bkrobust.mpdag_criterion.paths import possible_descendants

Node = str
Edge = tuple[str, str]

#: Two-sided 95 % normal quantile.
Z975 = 1.959963984540054


def chain_component(cpdag: MPDAG, x: Node) -> frozenset[Node] | None:
    """The undirected (chain) component of ``cpdag`` containing ``x``, or ``None``."""
    for comp in undirected_components(cpdag):
        if x in comp:
            return comp
    return None


def local_claims(cpdag: MPDAG, claims: Iterable[Edge], comp: frozenset[Node] | None) -> list[Edge]:
    """The claims that orient an undirected edge of ``cpdag`` inside ``comp``, sorted."""
    if comp is None:
        return []
    und = cpdag.undirected_edges
    return sorted((a, b) for a, b in claims if a in comp and b in comp and canon(a, b) in und)


def component_graph(cpdag: MPDAG, comp: frozenset[Node]) -> MPDAG:
    """The subgraph of ``cpdag`` induced on one chain component (all edges undirected)."""
    return induced_subgraph(cpdag, comp)


def merge_component(base: MPDAG, comp: frozenset[Node], local: MPDAG) -> MPDAG:
    """``base`` with the edges inside ``comp`` replaced by those of ``local``."""
    directed = {(a, b) for a, b in base.directed_edges if not (a in comp and b in comp)}
    undirected = {e for e in base.undirected_edges if not (e[0] in comp and e[1] in comp)}
    return MPDAG(
        base.nodes,
        directed | set(local.directed_edges),
        undirected | set(local.undirected_edges),
    )


@dataclass(frozen=True)
class ParentSets:
    """The possible parent sets of a treatment in one graph.

    Attributes:
        local: For each possible parent set, the part inside the component
            (the rest is the fixed outside parents, identical for all of them).
        closed: For each, the component graph after orienting the treatment's
            incident edges that way and closing under Meek's rules.
        candidates: Clique candidates tried.
        censored: True when the candidate count exceeded the cap and nothing was
            enumerated; ``local`` is then empty.
    """

    local: tuple[frozenset[Node], ...]
    closed: tuple[MPDAG, ...]
    candidates: int
    censored: bool


def _sibling_cliques(g: MPDAG, sib: Sequence[Node], cap: int) -> list[frozenset[Node]] | None:
    """Every clique (including the empty set) of the siblings, or ``None`` above ``cap``."""
    out: list[frozenset[Node]] = []

    def extend(current: list[Node], start: int) -> bool:
        out.append(frozenset(current))
        if len(out) > cap:
            return False
        for i in range(start, len(sib)):
            v = sib[i]
            if all(g.has_edge(v, u) for u in current):
                current.append(v)
                if not extend(current, i + 1):
                    return False
                current.pop()
        return True

    return out if extend([], 0) else None


def possible_parent_sets(local: MPDAG, x: Node, cap: int = 4096) -> ParentSets:
    """Semi-local IDA on one component graph: the parent sets of ``x`` some extension realises.

    Args:
        local: A Meek-closed graph on the treatment's chain component.
        x: The treatment.
        cap: Maximum number of clique candidates; above it the call is censored.

    Returns:
        The possible local parent sets with their closed graphs.
    """
    if x not in local.nodes:
        return ParentSets((frozenset(),), (local,), 1, False)
    fixed = frozenset(local.parents(x))
    sib = sorted(local.neighbors(x))
    cliques = _sibling_cliques(local, sib, cap)
    if cliques is None:
        return ParentSets((), (), cap + 1, True)
    sets: list[frozenset[Node]] = []
    graphs: list[MPDAG] = []
    for p in cliques:
        orient = [(v, x) for v in sorted(p)] + [(x, v) for v in sib if v not in p]
        g = apply_orientations(local, orient) if orient else local
        if g is not None:
            sets.append(fixed | p)
            graphs.append(g)
    return ParentSets(tuple(sets), tuple(graphs), len(cliques), False)


@dataclass(frozen=True)
class BallElement:
    """One retraction: the retracted claim indices and the re-closed component graph."""

    retracted: tuple[int, ...]
    graph: MPDAG | None


@dataclass(frozen=True)
class Ball:
    """The retraction ball over the local claims, enumerated by depth.

    Attributes:
        elements: Every enumerated retraction, depth by depth, in
            :func:`itertools.combinations` order.
        depth_done: The deepest depth enumerated completely.
        status: ``exhaustive`` when every depth up to ``min(max_depth, |K|)``
            was enumerated, else ``censored_at_depth_d``.
    """

    elements: tuple[BallElement, ...]
    depth_done: int
    status: str


def retraction_ball(
    local_base: MPDAG, claims: Sequence[Edge], max_depth: int = 3, subset_budget: int = 300
) -> Ball:
    """Enumerate ``B_r`` for ``r <= max_depth`` on one component, under a subset budget.

    The budget counts non-empty retractions cumulatively over depths, as the
    ledger's claim radius does: a depth whose subsets would push the count past
    ``subset_budget`` is not started and the ball is censored there.

    Args:
        local_base: The component graph with nothing oriented (the class).
        claims: The local claims, in a fixed order.
        max_depth: Deepest retraction size.
        subset_budget: Cumulative cap on non-empty retractions.

    Returns:
        The enumerated ball.
    """
    k = list(claims)
    elements = [BallElement((), apply_orientations(local_base, k) if k else local_base)]
    used = 0
    depth_done = 0
    status = "exhaustive"
    for depth in range(1, min(max_depth, len(k)) + 1):
        n_sub = math.comb(len(k), depth)
        if used + n_sub > subset_budget:
            status = f"censored_at_depth_{depth}"
            break
        for subset in itertools.combinations(range(len(k)), depth):
            drop = set(subset)
            kept = [e for i, e in enumerate(k) if i not in drop]
            g = apply_orientations(local_base, kept) if kept else local_base
            elements.append(BallElement(tuple(subset), g))
        used += n_sub
        depth_done = depth
    return Ball(tuple(elements), depth_done, status)


def effect_budget_available(ball: Ball, n_claims: int, r: int) -> bool:
    """Whether the interval at budget ``r`` is exact under the ball's enumeration.

    Budgets past the number of claims equal the budget at the number of claims
    (every claim retracted is the class itself), so they are available once the
    ball reached that depth.
    """
    return min(r, n_claims) <= ball.depth_done


def regression_beta_se(
    sigmas: np.ndarray,
    idx: dict[Node, int],
    x: Node,
    y: Node,
    z: Iterable[Node],
    n: int | None,
) -> tuple[np.ndarray, np.ndarray]:
    """Coefficient of ``x`` regressing ``y`` on ``{x} + z``, and its OLS standard error.

    Args:
        sigmas: A stack of covariance matrices, shape ``(m, p, p)``: the
            population one or sample ones with divisor ``n - 1``.
        idx: Node to row index.
        x: The treatment.
        y: The outcome.
        z: The conditioning set (``x`` and ``y`` are dropped from it).
        n: Sample size, or ``None`` for the population (standard error zero).

    Returns:
        ``(beta, se)``, each of shape ``(m,)``. When ``y`` is in ``z`` the IDA
        convention applies and both are zero.
    """
    m = sigmas.shape[0]
    if y in set(z):
        return np.zeros(m), np.zeros(m)
    beta, s_y, s_x = regression_moments(sigmas, idx, x, y, z)
    if n is None:
        return beta, np.zeros(m)
    dof = n - len(set(z) - {x, y}) - 2
    with np.errstate(divide="ignore", invalid="ignore"):
        se = np.sqrt(np.maximum(s_y, 0.0) / (dof * s_x))
    return beta, se


def regression_moments(
    sigmas: np.ndarray, idx: dict[Node, int], x: Node, y: Node, z: Iterable[Node]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The coefficient of ``x``, ``Var(y | x, z)`` and ``Var(x | z)``, for a covariance stack.

    Args:
        sigmas: Covariance matrices, shape ``(m, p, p)``.
        idx: Node to row index.
        x: The treatment.
        y: The outcome, not in ``z``.
        z: The conditioning set.

    Returns:
        Three arrays of shape ``(m,)``.
    """
    zl = sorted(set(z) - {x, y})
    reg = [idx[x]] + [idx[v] for v in zl]
    a = sigmas[:, reg][:, :, reg]
    b = sigmas[:, reg, idx[y]]
    coef = np.linalg.solve(a, b[:, :, None])[:, :, 0]
    s_y = sigmas[:, idx[y], idx[y]] - np.einsum("ij,ij->i", b, coef)
    if zl:
        zi = [idx[v] for v in zl]
        az = sigmas[:, zi][:, :, zi]
        bz = sigmas[:, zi, idx[x]]
        cz = np.linalg.solve(az, bz[:, :, None])[:, :, 0]
        s_x = sigmas[:, idx[x], idx[x]] - np.einsum("ij,ij->i", bz, cz)
    else:
        s_x = sigmas[:, idx[x], idx[x]].copy()
    return coef[:, 0], s_y, s_x


def partial_r2(
    sigma: np.ndarray,
    idx: dict[Node, int],
    target: Node,
    block: Iterable[Node],
    given: Iterable[Node],
) -> float:
    """Population partial R-squared of ``target`` on ``block`` given ``given``."""

    def resid(cond: list[Node]) -> float:
        t = idx[target]
        if not cond:
            return float(sigma[t, t])
        ci = [idx[v] for v in cond]
        coef = np.linalg.solve(sigma[np.ix_(ci, ci)], sigma[ci, t])
        return float(sigma[t, t] - sigma[ci, t] @ coef)

    g = sorted(set(given) - {target})
    u = sorted(set(block) - set(g) - {target})
    if not u:
        return 0.0
    base = resid(g)
    full = resid(sorted(set(g) | set(u)))
    return 0.0 if base <= 0 else max(0.0, 1.0 - full / base)


def zero_tolerance(sigma: np.ndarray, idx: dict[Node, int], x: Node, y: Node) -> float:
    """Scale-aware tolerance for a population coefficient being zero: ``1e-6 * sd_y / sd_x``."""
    return 1e-6 * math.sqrt(float(sigma[idx[y], idx[y]]) / float(sigma[idx[x], idx[x]]))


def contains_zero(lo: float, hi: float, tol: float) -> bool:
    """Whether ``[lo, hi]`` contains zero up to ``tol`` (the rnull-real fix)."""
    return lo <= tol and hi >= -tol


def y_possible_descendant(
    full_base: MPDAG, comp: frozenset[Node], closed: MPDAG, x: Node, y: Node
) -> bool:
    """Whether ``y`` is a possible descendant of ``x`` with the component oriented as ``closed``."""
    return y in possible_descendants(merge_component(full_base, comp, closed), x)


def closure_by_components(
    cpdag: MPDAG,
    claims: Iterable[Edge],
    components: Sequence[frozenset[Node]] | None = None,
    cache: dict | None = None,
) -> MPDAG | None:
    """``Meek(cpdag, claims)`` assembled from one closure per chain component.

    Equal to :func:`~bkrobust.demo.meek.apply_orientations` on the full CPDAG when
    every claim orients an undirected edge (the component reduction; tested), and much
    cheaper on large graphs. Falls back to the full closure otherwise.

    Args:
        cpdag: The CPDAG.
        claims: Orientation claims.
        components: The CPDAG's chain components, if already computed.
        cache: Optional memo of per-component closures keyed by (component, claims).

    Returns:
        The closed graph, or ``None`` when some component's closure fails.
    """
    k = list(claims)
    und = cpdag.undirected_edges
    if any(canon(a, b) not in und for a, b in k):
        return apply_orientations(cpdag, k)
    comps = undirected_components(cpdag) if components is None else components
    directed = set(cpdag.directed_edges)
    undirected = set(und)
    for comp in comps:
        k_loc = local_claims(cpdag, k, comp)
        if not k_loc:
            continue
        key = (comp, tuple(k_loc))
        if cache is not None and key in cache:
            local = cache[key]
        else:
            local = apply_orientations(component_graph(cpdag, comp), k_loc)
            if cache is not None:
                cache[key] = local
        if local is None:
            return None
        undirected -= {e for e in und if e[0] in comp and e[1] in comp}
        directed |= set(local.directed_edges)
        undirected |= set(local.undirected_edges)
    return MPDAG(cpdag.nodes, directed, undirected)
