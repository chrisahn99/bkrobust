"""Adjustment-set validity, optimal adjustment sets, and linear-SEM bias/variance.

Everything here is defined relative to the shared :class:`~bkrobust.demo.graph.MPDAG`
container. Two layers are provided throughout:

* A *DAG-level* layer (``is_valid_adjustment_set_dag``, ``optimal_adjustment_set_dag``,
  ...) that assumes every edge is oriented, i.e. ``dag.is_dag()`` holds. This is where
  the actual causal-inference definitions live.
* An *MPDAG-level* layer (``is_valid_adjustment_set_mpdag``, ``optimal_adjustment_set_mpdag``,
  ``all_valid_adjustment_sets_mpdag``) that lifts the DAG-level notions to a partially
  oriented graph by quantifying over every consistent DAG extension -- "valid/optimal
  in every world compatible with what is currently oriented". This mirrors the standard
  treatment of adjustment for CPDAGs/MPDAGs (Perkovic et al., 2018; Henckel, Perkovic
  and Maathuis, 2022).

d-separation is implemented once, independently, via graph moralization (Lauritzen,
1996), and every other predicate in this module is expressed in terms of it or of the
plain reachability primitives already on :class:`~bkrobust.demo.graph.MPDAG`.

Written to run on Python 3.9 as well as the repository's 3.11 target: builtin
generics and ``X | Y`` unions appear only in annotations, which ``from __future__
import annotations`` defers to strings and are never evaluated at runtime.
"""

from __future__ import annotations

import itertools
from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np

from bkrobust.demo.graph import MPDAG

Node = str


def _as_node_set(value: Node | Iterable[Node]) -> set[Node]:
    """Normalise a single node or an iterable of nodes to a plain set."""
    if isinstance(value, str):
        return {value}
    return set(value)


# --------------------------------------------------------------------------------
# d-separation
# --------------------------------------------------------------------------------


def is_dseparated(
    dag: MPDAG, a: Node | Iterable[Node], b: Node | Iterable[Node], z: Iterable[Node]
) -> bool:
    """Whether node set ``a`` is d-separated from node set ``b`` given ``z`` in ``dag``.

    ``dag`` must be a fully oriented DAG (``dag.is_dag()``); d-separation is only
    defined on DAGs. ``a`` and ``b`` may each be a single node or a set of nodes.

    Implemented by the classical moralization test (Lauritzen, "Graphical Models",
    1996, Prop 3.25), which is a *second*, independent construction from the
    back-door-criterion code in this module -- it does not enumerate paths at all:

    1. Restrict to the ancestral set ``An(a, b, z)``: the nodes in ``a``, ``b``
       and ``z`` together with all of their ancestors. Every path relevant to
       separating ``a`` from ``b`` given ``z`` stays entirely inside this set.
    2. Moralize: take the skeleton of the restricted DAG and additionally connect,
       for every node, every pair of its parents ("marry the parents") -- this is
       exactly what makes v-structures' colliders behave correctly without having
       to reason about arrowheads explicitly.
    3. ``a`` and ``b`` are d-separated by ``z`` iff removing the nodes in ``z``
       disconnects ``a`` from ``b`` in this moral graph -- checked here by a BFS
       from ``a`` that treats ``z`` as blocked.

    Args:
        dag: A fully oriented DAG.
        a: A node, or set of nodes, on one side.
        b: A node, or set of nodes, on the other side.
        z: The conditioning set (may be empty).

    Returns:
        True iff ``a`` and ``b`` are d-separated given ``z``.
    """
    a_set = _as_node_set(a)
    b_set = _as_node_set(b)
    z_set = set(z)

    relevant = a_set | b_set | z_set
    ancestral: set[Node] = set(relevant)
    for n in relevant:
        ancestral |= dag.ancestors(n)

    moral_adj: dict[Node, set[Node]] = {n: set() for n in ancestral}
    for tail, head in dag.directed_edges:
        if tail in ancestral and head in ancestral:
            moral_adj[tail].add(head)
            moral_adj[head].add(tail)
    for n in ancestral:
        parents_n = [p for p in dag.parents(n) if p in ancestral]
        for p1, p2 in itertools.combinations(parents_n, 2):
            moral_adj[p1].add(p2)
            moral_adj[p2].add(p1)

    start = [n for n in a_set if n not in z_set]
    visited: set[Node] = set(start)
    stack: list[Node] = list(start)
    while stack:
        cur = stack.pop()
        if cur in b_set:
            return False
        for nxt in moral_adj.get(cur, ()):
            if nxt in z_set or nxt in visited:
                continue
            visited.add(nxt)
            stack.append(nxt)
    return True


# --------------------------------------------------------------------------------
# Back-door criterion (DAG level)
# --------------------------------------------------------------------------------


def is_valid_adjustment_set_dag(dag: MPDAG, x: Node, y: Node, z: Iterable[Node]) -> bool:
    """Whether ``z`` satisfies the back-door criterion for ``(x, y)`` in ``dag``.

    ``dag`` must be a fully oriented DAG. The back-door criterion (Pearl,
    "Causality", 2009, Def. 3.3.1) requires:

    (a) no node in ``z`` is a descendant of ``x`` (and ``x`` not in ``z``,
        ``y`` not in ``z``); and
    (b) ``z`` blocks every back-door path from ``x`` to ``y``, i.e. every path
        that leaves ``x`` via an edge *into* ``x``.

    Condition (b) is checked via ``is_dseparated`` rather than by enumerating
    paths directly, using the standard equivalence (Pearl, 2009, proof of
    Thm. 3.3.4): let ``dag_x_under`` be ``dag`` with every edge *out of* ``x``
    deleted. Deleting those edges destroys exactly the paths that leave ``x``
    forward (the non-back-door paths) while leaving every edge *into* ``x`` --
    and hence every back-door path -- unchanged. So ``z`` blocks every back-door
    path in ``dag`` iff ``x`` and ``y`` are d-separated by ``z`` in
    ``dag_x_under``. This reuses the single d-separation implementation rather
    than a second, path-enumerating one.

    Args:
        dag: A fully oriented DAG.
        x: The treatment node.
        y: The outcome node.
        z: The candidate adjustment set.

    Returns:
        True iff ``z`` is a valid back-door adjustment set for ``(x, y)``.
    """
    z_set = frozenset(z)
    if x in z_set or y in z_set:
        return False
    if z_set & dag.descendants(x):
        return False

    trimmed_directed = [(tail, head) for tail, head in dag.directed_edges if tail != x]
    dag_x_under = MPDAG(dag.nodes, trimmed_directed, dag.undirected_edges)
    return is_dseparated(dag_x_under, x, y, z_set)


# --------------------------------------------------------------------------------
# Back-door criterion (MPDAG level)
# --------------------------------------------------------------------------------


def is_valid_adjustment_set_mpdag(g: MPDAG, x: Node, y: Node, z: Iterable[Node]) -> bool:
    """Whether ``z`` is a valid back-door adjustment set in *every* DAG extension of ``g``.

    An MPDAG represents a set of DAGs (one per way of orienting its remaining
    undirected edges consistently). ``z`` is declared valid here only if it
    satisfies the back-door criterion in each and every one of them -- so validity
    is certified by ``g`` alone, without needing to know which extension is the
    true one. If ``g`` has no consistent DAG extension at all (which should not
    happen for a graph produced by ``enumerate_dag_extensions`` on a sound MPDAG),
    this conservatively returns False.

    Args:
        g: A (possibly partially oriented) MPDAG.
        x: The treatment node.
        y: The outcome node.
        z: The candidate adjustment set.

    Returns:
        True iff ``z`` is valid in every DAG extension of ``g``.
    """
    from bkrobust.demo.meek import enumerate_dag_extensions

    extensions = enumerate_dag_extensions(g)
    if not extensions:
        return False
    z_set = frozenset(z)
    return all(is_valid_adjustment_set_dag(d, x, y, z_set) for d in extensions)


def all_valid_adjustment_sets_mpdag(
    g: MPDAG, x: Node, y: Node, max_size: int | None = None
) -> list[frozenset[Node]]:
    """Every subset of ``nodes - {x, y}`` that is a valid adjustment set per ``g``.

    Enumerates candidate subsets and keeps those for which
    ``is_valid_adjustment_set_dag`` holds in every DAG extension of ``g``
    (equivalently, ``is_valid_adjustment_set_mpdag`` would return True). This is
    exponential in the number of candidate nodes, which is acceptable at the sizes
    (``n <= 10``) this demonstration targets.

    As a cheap prune, any node that is a descendant of ``x`` in *every* extension
    is dropped from consideration up front: condition (a) of the back-door
    criterion would reject any set containing it, in every extension, regardless
    of what else is in the set.

    Args:
        g: A (possibly partially oriented) MPDAG.
        x: The treatment node.
        y: The outcome node.
        max_size: If given, only subsets of size at most this are considered.

    Returns:
        Valid adjustment sets, sorted by ``(size, sorted labels)``.
    """
    from bkrobust.demo.meek import enumerate_dag_extensions

    extensions = enumerate_dag_extensions(g)
    if not extensions:
        return []

    always_descendant_of_x: set[Node] | None = None
    for d in extensions:
        desc = d.descendants(x)
        always_descendant_of_x = (
            set(desc) if always_descendant_of_x is None else always_descendant_of_x & desc
        )
    if always_descendant_of_x is None:
        always_descendant_of_x = set()

    candidates = sorted(n for n in g.nodes if n not in (x, y) and n not in always_descendant_of_x)
    limit = len(candidates) if max_size is None else min(max_size, len(candidates))

    results: list[frozenset[Node]] = []
    for size in range(limit + 1):
        for combo in itertools.combinations(candidates, size):
            z_set = frozenset(combo)
            if all(is_valid_adjustment_set_dag(d, x, y, z_set) for d in extensions):
                results.append(z_set)
    results.sort(key=lambda s: (len(s), sorted(s)))
    return results


# --------------------------------------------------------------------------------
# Optimal adjustment set (Henckel, Perkovic and Maathuis, 2022)
# --------------------------------------------------------------------------------


def optimal_adjustment_set_dag(dag: MPDAG, x: Node, y: Node) -> set[Node]:
    r"""The Henckel-Perkovic-Maathuis optimal adjustment set for ``(x, y)`` in ``dag``.

    Defined as ``O = pa(cn(x, y)) \ (cn(x, y) union {x})``, where ``cn(x, y)`` --
    the "causal nodes" -- is the set of nodes lying on some proper (simple)
    directed path from ``x`` to ``y``, excluding ``x`` itself but including ``y``
    (Henckel, Perkovic & Maathuis, "Graphical Criteria for Efficient Total Effect
    Estimation via Adjustment in Causal Linear Models", JRSS-B 2022, Def. 3 for
    ``cn`` / forb, Def. 4-5 for ``O``). Among all valid adjustment sets, ``O``
    gives the asymptotically most efficient linear estimator.

    ``cn(x, y)`` is computed as ``descendants(x) & (ancestors(y) | {y})``: a node
    ``v`` lies on some directed path from ``x`` to ``y`` iff it is reachable from
    ``x`` (a descendant of ``x``) *and* can itself reach ``y`` (an ancestor of
    ``y``, or ``y`` itself) -- in a DAG these two sub-paths concatenate into a
    single simple path with no repeated nodes, since a node cannot be both an
    ancestor and a descendant of another without a cycle.

    Args:
        dag: A fully oriented DAG.
        x: The treatment node.
        y: The outcome node.

    Returns:
        The optimal adjustment set ``O`` (possibly empty).
    """
    cn = dag.descendants(x) & (dag.ancestors(y) | {y})
    pa_cn: set[Node] = set()
    for n in cn:
        pa_cn |= dag.parents(n)
    return pa_cn - cn - {x}


def optimal_adjustment_set_mpdag(g: MPDAG, x: Node, y: Node) -> set[Node] | None:
    """The optimal adjustment set for ``(x, y)``, if identified by the MPDAG ``g`` alone.

    ``optimal_adjustment_set_dag`` is computed in every DAG extension of ``g``. If
    every extension agrees on the resulting set, that set -- being the same
    regardless of which extension is the true DAG -- is returned. If extensions
    disagree, the optimal set is not identified by ``g`` alone (different
    completions of the remaining undirected edges would demand different
    adjustment sets for optimal efficiency), and ``None`` is returned rather than
    guessing one extension's answer. This mirrors ``is_valid_adjustment_set_mpdag``'s
    "must hold in every extension" stance, but for optimality rather than mere
    validity there is no natural intersection/union fallback, so disagreement is
    reported as non-identification instead of any set.

    Args:
        g: A (possibly partially oriented) MPDAG.
        x: The treatment node.
        y: The outcome node.

    Returns:
        The optimal adjustment set if it agrees across every DAG extension of
        ``g``, else None.
    """
    from bkrobust.demo.meek import enumerate_dag_extensions

    extensions = enumerate_dag_extensions(g)
    if not extensions:
        return None

    sets = [optimal_adjustment_set_dag(d, x, y) for d in extensions]
    first = sets[0]
    if all(s == first for s in sets[1:]):
        return set(first)
    return None


# --------------------------------------------------------------------------------
# Linear-Gaussian SEM
# --------------------------------------------------------------------------------


@dataclass
class LinearSEM:
    """A linear-Gaussian structural equation model over ``dag``.

    Each node ``v`` follows ``v = sum_{p in pa(v)} weights[(p, v)] * p + noise_v``,
    with independent noises ``noise_v ~ N(0, noise_var[v])``.

    Attributes:
        dag: The (fully oriented) causal DAG.
        weights: Structural coefficients, keyed ``(tail, head)`` for edge
            ``tail -> head``, i.e. row/child convention: ``weights[(p, v)]`` is
            the coefficient of parent ``p`` in the structural equation for ``v``.
        noise_var: Exogenous noise variance for each node.
    """

    dag: MPDAG
    weights: dict[tuple[Node, Node], float]
    noise_var: dict[Node, float]

    def _weight_matrix(self) -> tuple[np.ndarray, list[Node], dict[Node, int]]:
        """The structural weight matrix ``B`` with ``B[child, parent] = weight``.

        Row index is the child (the equation being written), column index is the
        parent (the regressor). This is the convention under which
        ``Sigma = (I - B)^-1 Omega (I - B)^-T`` is the population covariance: for
        the 2-node model ``X -> Y`` with weight ``w`` and noise variances
        ``sigma_x^2``, ``sigma_y^2``, ``B = [[0, 0], [w, 0]]`` (rows/cols ordered
        X, Y since node order is sorted and X < Y), so
        ``(I - B)^-1 = [[1, 0], [w, 1]]`` and
        ``Sigma = [[sigma_x^2, w*sigma_x^2], [w*sigma_x^2, w^2*sigma_x^2 + sigma_y^2]]``,
        which is exactly ``Var(X) = sigma_x^2``, ``Cov(X, Y) = w * Var(X)``,
        ``Var(Y) = w^2 * Var(X) + sigma_y^2`` -- the textbook answer, checked by
        hand in the test suite.
        """
        nodes = list(self.dag.nodes)
        idx = {n: i for i, n in enumerate(nodes)}
        p = len(nodes)
        b = np.zeros((p, p))
        for (tail, head), w in self.weights.items():
            b[idx[head], idx[tail]] = w
        return b, nodes, idx

    def covariance(self) -> np.ndarray:
        """The exact population covariance matrix, in ``dag.nodes`` order.

        Closed form ``Sigma = (I - B)^-1 Omega (I - B)^-T``, with ``B`` the
        weight matrix (row = child) and ``Omega = diag(noise_var)``. No sampling
        is involved; this is the exact population moment.
        """
        b, nodes, _ = self._weight_matrix()
        p = len(nodes)
        omega = np.diag([self.noise_var[n] for n in nodes])
        i_minus_b_inv = np.linalg.inv(np.eye(p) - b)
        # numpy 2.0 on the macOS Accelerate BLAS backend emits spurious
        # "divide by zero"/"overflow encountered in matmul" RuntimeWarnings on
        # well-conditioned inputs. Verified spurious: over 500 draws the result
        # was finite, SPD, and equal to the Neumann-series form (I+B+...+B^n)
        # in every case, with det(I-B)=1 and condition number ~10. Suppressed
        # narrowly, with a real finiteness guard below so a genuine blow-up
        # cannot hide behind the suppression.
        with np.errstate(all="ignore"):
            sigma = i_minus_b_inv @ omega @ i_minus_b_inv.T
        if not np.all(np.isfinite(sigma)):
            raise FloatingPointError("non-finite covariance; the SEM is degenerate")
        return sigma

    def true_total_effect(self, x: Node, y: Node) -> float:
        """Sum, over every directed path ``x -> ... -> y``, of the product of edge weights.

        Computed by the standard recursive decomposition: the total effect from a
        node ``v`` to ``y`` is 1 if ``v == y``, else the sum over each child ``c``
        of ``v`` of ``weight(v, c) * total_effect(c, y)``. Memoized over nodes
        (safe since ``dag`` is acyclic, so this recursion always terminates).
        """
        memo: dict[Node, float] = {}

        def total_from(v: Node) -> float:
            if v in memo:
                return memo[v]
            if v == y:
                memo[v] = 1.0
                return 1.0
            total = 0.0
            # sorted(): children() returns a set, and floating-point addition is
            # not associative, so summing in set-iteration order makes the result
            # depend on PYTHONHASHSEED at the last bit.
            for c in sorted(self.dag.children(v)):
                total += self.weights[(v, c)] * total_from(c)
            memo[v] = total
            return total

        return total_from(x)


def random_sem(
    dag: MPDAG,
    rng: np.random.Generator,
    weight_range: tuple[float, float] = (0.3, 1.5),
    noise_range: tuple[float, float] = (0.5, 1.5),
) -> LinearSEM:
    """Draw a random :class:`LinearSEM` over ``dag``.

    Each structural coefficient has a magnitude drawn uniformly from
    ``weight_range`` (so it is bounded away from 0) and an independently random
    sign; each noise variance is drawn uniformly from ``noise_range``. Randomness
    is drawn exclusively from the passed ``rng``; the global NumPy RNG is never
    touched. Edges are visited in sorted order so the draw sequence does not
    depend on set iteration order, which varies with ``PYTHONHASHSEED`` across
    processes -- results are reproducible from the seed within and across runs.

    Args:
        dag: A fully oriented DAG.
        rng: A ``numpy.random.Generator`` to draw from.
        weight_range: ``(low, high)`` bounds on coefficient magnitude.
        noise_range: ``(low, high)`` bounds on noise variance.

    Returns:
        A randomly parameterised :class:`LinearSEM` over ``dag``.
    """
    weights: dict[tuple[Node, Node], float] = {}
    # sorted(), NOT the raw frozenset: iterating a set of strings consumes the
    # RNG in an order that depends on PYTHONHASHSEED, which differs per process.
    # Without this the same seed gives different SEMs in different runs.
    for tail, head in sorted(dag.directed_edges):
        magnitude = rng.uniform(weight_range[0], weight_range[1])
        sign = rng.choice(np.array([-1.0, 1.0]))
        weights[(tail, head)] = float(sign * magnitude)
    noise_var = {n: float(rng.uniform(noise_range[0], noise_range[1])) for n in dag.nodes}
    return LinearSEM(dag=dag, weights=weights, noise_var=noise_var)


def _regress_coeffs(sigma: np.ndarray, reg_idx: list[int], target_idx: int) -> np.ndarray:
    """Population OLS coefficients of ``target`` on the regressors ``reg_idx``.

    Solves the normal equations ``Sigma_rr @ beta = Sigma_r,target`` on the
    population covariance matrix directly -- no sampling.
    """
    sigma_rr = sigma[np.ix_(reg_idx, reg_idx)]
    sigma_r_target = sigma[reg_idx, target_idx]
    return np.linalg.solve(sigma_rr, sigma_r_target)


def adjusted_estimand(sem: LinearSEM, x: Node, y: Node, z: Iterable[Node]) -> float:
    """The population regression coefficient on ``x`` from regressing ``y`` on ``{x} + z``.

    Computed in closed form from ``sem.covariance()`` by solving the normal
    equations on the population covariance matrix -- this is the large-sample
    limit of OLS, not an estimate from sampled data.

    Args:
        sem: The linear SEM.
        x: The treatment node (its coefficient is what is returned).
        y: The outcome node.
        z: The adjustment set.

    Returns:
        The population coefficient on ``x`` in the regression of ``y`` on
        ``{x} + z``.
    """
    nodes = list(sem.dag.nodes)
    idx = {n: i for i, n in enumerate(nodes)}
    z_list = sorted(set(z) - {x, y})
    regressors = [x, *z_list]
    sigma = sem.covariance()
    reg_idx = [idx[n] for n in regressors]
    coeffs = _regress_coeffs(sigma, reg_idx, idx[y])
    return float(coeffs[0])


def bias(sem: LinearSEM, x: Node, y: Node, z: Iterable[Node]) -> float:
    """``adjusted_estimand(sem, x, y, z) - sem.true_total_effect(x, y)``."""
    return adjusted_estimand(sem, x, y, z) - sem.true_total_effect(x, y)


def asymptotic_variance(sem: LinearSEM, x: Node, y: Node, z: Iterable[Node]) -> float:
    """The n-free part of the asymptotic variance of the OLS-adjusted estimator.

    For the linear-Gaussian OLS estimator that regresses ``y`` on ``{x} + z``,
    the coefficient on ``x`` has asymptotic variance
    ``sigma^2_{y|x,z} / (n * sigma^2_{x|z})``, where ``sigma^2_{y|x,z}`` is the
    residual variance of ``y`` given ``(x, z)`` and ``sigma^2_{x|z}`` is the
    residual variance of ``x`` given ``z`` (population values; this is the
    Frisch-Waugh-Lovell partialling-out form of the OLS variance). This function
    returns the quantity with the sample size ``n`` divided out, i.e.
    ``sigma^2_{y|x,z} / sigma^2_{x|z}``: it does not depend on ``n`` because both
    residual variances above are population moments, not sample ones, so callers
    multiply by ``1 / n`` themselves for a given sample size.

    Args:
        sem: The linear SEM.
        x: The treatment node.
        y: The outcome node.
        z: The adjustment set.

    Returns:
        ``sigma^2_{y|x,z} / sigma^2_{x|z}``, the n-free asymptotic variance factor.
    """
    nodes = list(sem.dag.nodes)
    idx = {n: i for i, n in enumerate(nodes)}
    z_list = sorted(set(z) - {x, y})
    sigma = sem.covariance()

    if z_list:
        z_idx = [idx[n] for n in z_list]
        coeffs_x = _regress_coeffs(sigma, z_idx, idx[x])
        sigma_zx = sigma[z_idx, idx[x]]
        resid_var_x = float(sigma[idx[x], idx[x]] - sigma_zx @ coeffs_x)
    else:
        resid_var_x = float(sigma[idx[x], idx[x]])

    regressors = [x, *z_list]
    reg_idx = [idx[n] for n in regressors]
    coeffs_y = _regress_coeffs(sigma, reg_idx, idx[y])
    sigma_ry = sigma[reg_idx, idx[y]]
    resid_var_y = float(sigma[idx[y], idx[y]] - sigma_ry @ coeffs_y)

    return resid_var_y / resid_var_x
