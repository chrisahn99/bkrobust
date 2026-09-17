"""The worst-case bias functional ``B``, and the ambiguity set it maximises over.

The quantity this module computes is the one thing ``bias_stats`` in
:mod:`bkrobust.core.oracle` deliberately does not compute: an **exact supremum**
rather than a sampled mean.

The difference is not cosmetic. ``bias_stats`` redraws a linear SEM on every DAG
of ``[G]`` and averages the absolute bias over extensions and draws. Two things
follow, and both are fatal for a certificate:

* A mean over a *larger* set can be *smaller*, so ``G -> mean_abs_bias(G)`` is
  not monotone under model inclusion. The set ``{mean > eps}`` is then not
  upward-closed, none of ``THEOREMS.md`` applies to it, and in particular the
  upward search is not exact for it.
* The SEM is redrawn per element, so two elements are compared under different
  data-generating parameters.

What is defined here instead. The analyst holds **one** observational law, with
population covariance ``Sigma``, and reports **one** number,

    theta_Z = beta(Z; Sigma),

the coefficient on ``X`` in the population regression of ``Y`` on ``{X} + Z``.
It does not depend on which DAG is the truth. What the knowledge state ``G``
leaves open is *which* DAG is the truth, and each candidate ``D`` in ``[G]``
assigns the effect its own value

    tau_D = beta(pa_D(X); Sigma),

the classical IDA quantity (0 by convention when ``Y`` is a parent of ``X``).
The worst-case bias is then

    B(G) = max { |theta_Z - tau| : tau in T(G) },    T(G) = { tau_D : D in [G] }.

``T(G)`` grows with ``[G]``, so ``B`` is monotone under ``<=`` by construction
(``docs/R_EPSILON_THEORY.md``, Theorem B), and it is exactly zero whenever ``Z``
is valid throughout ``[G]`` (Theorem A). Crucially the *truth* lies in ``T(G)``
for any state that the analyst's surviving claims admit, so ``B`` upper-bounds
the analyst's realised error rather than a hypothetical one.

**Two evaluators, one on the hot path.** ``T(G)`` is read off the possible
parent sets of ``X``. Those are obtained either by enumerating ``[G]``
(exponential, exact, used only as a differential test) or semi-locally by
Theorem E: try each subset of the undirected neighbours of ``X`` as its extra
parents, Meek-close, and keep the subsets that do not FAIL. That is
``2 ** deg_undirected(X)`` polynomial closures and no DAG enumeration at all --
the same move Lemma O made for cover generation in
:mod:`bkrobust.search.exact_fast`.

Determinism: no randomness is used anywhere in this module. ``Sigma`` is a
population moment computed in closed form, sets are iterated in sorted order,
and float accumulation is over sorted lists.
"""

from __future__ import annotations

import itertools
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from bkrobust.demo.evaluate import LinearSEM
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions

Node = str
Edge = tuple[str, str]

#: Absolute tolerance below which a bias is reported as exactly zero. Theorem A
#: predicts algebraic zero inside the certified shell; what is observed is the
#: rounding error of solving the normal equations, measured at ~1e-16 on the
#: paper's running example (App. C, first three shells). Anything above this is
#: a real bias, anything below is float noise -- and the choice is recorded
#: here rather than inlined so that a result depending on it is auditable.
ZERO_TOL: float = 1e-9


def knowledge_of(cpdag: MPDAG, g: MPDAG) -> list[Edge]:
    """The orientations ``g`` adds to ``cpdag`` -- its knowledge set ``K_g``.

    Duplicated from :mod:`bkrobust.search.space_fixed` rather than imported, to
    keep this module's dependency surface to the frozen core.
    """
    return sorted(set(g.directed_edges) - set(cpdag.directed_edges))


def _regress_on_covariance(
    sigma: np.ndarray, regressor_idx: list[int], target_idx: int
) -> np.ndarray:
    """Population least-squares coefficients of ``target`` on ``regressor_idx``.

    Solves the normal equations on the population covariance directly. Uses
    ``lstsq`` rather than ``solve`` so that a singular regressor block -- which
    arises legitimately when two candidate parents are deterministically related
    in the drawn SEM -- returns the minimum-norm solution instead of raising.

    Args:
        sigma: Population covariance, in a fixed node order.
        regressor_idx: Column indices of the regressors; the coefficient of
            interest is the **first** of them.
        target_idx: Column index of the regressand.

    Returns:
        The coefficient vector, aligned with ``regressor_idx``.
    """
    sigma_rr = sigma[np.ix_(regressor_idx, regressor_idx)]
    sigma_rt = sigma[regressor_idx, target_idx]
    try:
        return np.linalg.solve(sigma_rr, sigma_rt)
    except np.linalg.LinAlgError:
        return np.linalg.lstsq(sigma_rr, sigma_rt, rcond=None)[0]


@dataclass(frozen=True)
class BiasContext:
    """Everything ``B`` needs that does not vary across knowledge states.

    Built once per instance by :func:`make_context` and then reused for every
    element visited, which is what makes the fixed-``Sigma`` semantics literal:
    the covariance is computed a single time, so no two states are ever compared
    under different data-generating parameters.

    Attributes:
        cpdag: The estimated CPDAG the space refines.
        x: Treatment.
        y: Outcome.
        z: The adjustment set, fixed once from ``G0``.
        sigma: Population covariance of the true SCM, in ``nodes`` order.
        nodes: Node order the covariance is indexed by.
        index: Node name to column index.
        theta_z: The analyst's estimand ``beta(Z; Sigma)``. One number.
        tau_true: The true total effect, for measuring realised error. ``nan``
            outside simulation. Never on the computation path of ``B`` -- ``B``
            is computed from quantities the analyst holds, and the truth is used
            only to check the certificate afterwards.
        scale_std: ``sd(X) / sd(Y)``, the standardising factor of Definition 3.
    """

    cpdag: MPDAG
    x: Node
    y: Node
    z: frozenset[Node]
    sigma: np.ndarray
    nodes: tuple[Node, ...]
    index: dict[Node, int]
    theta_z: float
    tau_true: float
    scale_std: float

    def regression_coefficient(self, adjust: frozenset[Node]) -> float:
        """``beta(adjust; Sigma)``: the coefficient on ``X`` regressing ``Y`` on ``{X} + adjust``.

        Args:
            adjust: The conditioning set. ``X`` and ``Y`` are removed from it,
                so passing a set that happens to contain either is harmless.

        Returns:
            The population coefficient on ``X``.
        """
        cols = [self.index[self.x]]
        cols += [self.index[n] for n in sorted(adjust - {self.x, self.y})]
        coeffs = _regress_on_covariance(self.sigma, cols, self.index[self.y])
        return float(coeffs[0])

    def total_effect_under(self, parents_of_x: frozenset[Node]) -> float:
        """``tau_D`` for any DAG ``D`` whose parent set for ``X`` is ``parents_of_x``.

        The IDA quantity. Adjusting for the parents of the treatment identifies
        the total effect in any DAG, so the value depends on ``D`` only through
        that set -- which is what makes the whole ambiguity set computable from
        possible parent sets alone.

        Args:
            parents_of_x: ``pa_D(X)``.

        Returns:
            The total effect of ``X`` on ``Y`` under ``D``; exactly ``0.0``
            when ``Y`` is a parent of ``X``, since then ``X`` does not cause
            ``Y`` at all.
        """
        if self.y in parents_of_x:
            return 0.0
        return self.regression_coefficient(frozenset(parents_of_x))


def make_context_from_covariance(
    sigma: np.ndarray,
    nodes: Sequence[Node],
    cpdag: MPDAG,
    x: Node,
    y: Node,
    z: frozenset[Node],
    *,
    tau_true: float | None = None,
) -> BiasContext:
    """Build the per-instance context from quantities an analyst actually holds.

    This is the honest entry point, and it is deliberately the one the
    certificate is built on: a covariance matrix (estimable from data), the
    estimated CPDAG, the asserted knowledge that produced ``G0``, and the set
    ``Z`` read off it. **No true DAG appears anywhere.** That is what lets the
    resulting ``r_eps`` be reported as a diagnostic computed from the analyst's
    own inputs, in the same sense the paper claims for ``r_val``.

    Args:
        sigma: Population (or estimated) covariance, indexed by ``nodes``.
        nodes: Node order the covariance is indexed by.
        cpdag: The CPDAG the analyst estimated.
        x: Treatment.
        y: Outcome.
        z: The adjustment set read off ``G0``.
        tau_true: The true effect, when it happens to be known -- in simulation.
            Recorded for auditing the certificate afterwards and never read on
            the computation path. ``None`` outside simulation.

    Returns:
        A :class:`BiasContext`.

    Raises:
        ValueError: If ``x`` or ``y`` is unknown or is in ``z``, or if ``sigma``
            does not match ``nodes``.
    """
    nodes = tuple(nodes)
    if sigma.shape != (len(nodes), len(nodes)):
        raise ValueError(f"sigma is {sigma.shape}, expected {(len(nodes), len(nodes))}")
    if x not in nodes or y not in nodes:
        raise ValueError(f"treatment/outcome not in the node set: {x!r}, {y!r}")
    if x in z or y in z:
        raise ValueError("the adjustment set must not contain the treatment or the outcome")

    index = {n: i for i, n in enumerate(nodes)}
    blank = BiasContext(
        cpdag=cpdag,
        x=x,
        y=y,
        z=frozenset(z),
        sigma=sigma,
        nodes=nodes,
        index=index,
        theta_z=0.0,
        tau_true=float("nan"),
        scale_std=1.0,
    )
    sd_x = float(np.sqrt(sigma[index[x], index[x]]))
    sd_y = float(np.sqrt(sigma[index[y], index[y]]))
    return BiasContext(
        cpdag=cpdag,
        x=x,
        y=y,
        z=frozenset(z),
        sigma=sigma,
        nodes=nodes,
        index=index,
        theta_z=blank.regression_coefficient(frozenset(z)),
        tau_true=float("nan") if tau_true is None else float(tau_true),
        scale_std=(sd_x / sd_y) if sd_y > 0 else 1.0,
    )


def make_context(
    sem: LinearSEM,
    cpdag: MPDAG,
    x: Node,
    y: Node,
    z: frozenset[Node],
) -> BiasContext:
    """Build the context from a true SEM -- the simulation entry point.

    Supplies ``Sigma`` (the law the analyst would see at infinite sample size)
    and records ``tau_true`` for auditing. Everything the certificate is
    computed from is what :func:`make_context_from_covariance` takes; the SEM
    adds only the audit quantity.

    Args:
        sem: The true linear-Gaussian SCM.
        cpdag: The CPDAG the analyst estimated.
        x: Treatment.
        y: Outcome.
        z: The adjustment set read off ``G0``.

    Returns:
        A :class:`BiasContext`.
    """
    return make_context_from_covariance(
        sem.covariance(),
        tuple(sem.dag.nodes),
        cpdag,
        x,
        y,
        z,
        tau_true=float(sem.true_total_effect(x, y)),
    )


# --- possible parent sets -------------------------------------------------


def possible_parent_sets_enumerated(g: MPDAG, x: Node) -> frozenset[frozenset[Node]]:
    """``{ pa_D(X) : D in [G] }`` by enumerating the DAG extensions.

    Exact by construction and exponential in the number of undirected edges.
    Kept as the reference implementation that
    :func:`possible_parent_sets_semilocal` is differentially tested against;
    never on the hot path.

    Args:
        g: The knowledge state.
        x: Treatment.

    Returns:
        The distinct parent sets of ``x`` across ``[g]``.
    """
    return frozenset(frozenset(d.parents(x)) for d in enumerate_dag_extensions(g))


def possible_parent_sets_semilocal(cpdag: MPDAG, g: MPDAG, x: Node) -> frozenset[frozenset[Node]]:
    """``{ pa_D(X) : D in [G] }`` by Theorem E, with no DAG enumeration.

    For each subset ``S`` of the undirected neighbours of ``x`` in ``g``, impose
    ``s -> x`` for ``s`` in ``S`` and ``x -> t`` for the rest, Meek-close against
    the CPDAG, and keep ``S`` when the closure does not FAIL. Meek's theorem
    makes non-FAIL equivalent to the existence of an extension realising exactly
    that parent set, which is the content of Theorem E in
    ``docs/R_EPSILON_THEORY.md``.

    Cost is ``2 ** |nb_g(x)|`` polynomial closures. ``nb_g(x)`` counts the
    *undirected* edges at the treatment, which the paper's corpus measurements
    put in the low single digits.

    Args:
        cpdag: The CPDAG, needed because the closure is taken against it.
        g: The knowledge state.
        x: Treatment.

    Returns:
        The distinct parent sets of ``x`` across ``[g]``.
    """
    fixed_parents = frozenset(g.parents(x))
    neighbours = sorted(g.neighbors(x))
    base = knowledge_of(cpdag, g)

    out: set[frozenset[Node]] = set()
    for size in range(len(neighbours) + 1):
        for chosen in itertools.combinations(neighbours, size):
            into = set(chosen)
            orientations = list(base)
            orientations += [(s, x) for s in sorted(into)]
            orientations += [(x, t) for t in neighbours if t not in into]
            if apply_orientations(cpdag, orientations) is not None:
                out.add(fixed_parents | frozenset(into))
    return frozenset(out)


def possible_parent_sets(
    cpdag: MPDAG, g: MPDAG, x: Node, *, method: str = "semilocal"
) -> frozenset[frozenset[Node]]:
    """Dispatch between the two evaluators.

    Args:
        cpdag: The CPDAG.
        g: The knowledge state.
        x: Treatment.
        method: ``"semilocal"`` (default, Theorem E) or ``"enumerate"``.

    Returns:
        The distinct parent sets of ``x`` across ``[g]``.

    Raises:
        ValueError: On an unknown method.
    """
    if method == "semilocal":
        return possible_parent_sets_semilocal(cpdag, g, x)
    if method == "enumerate":
        return possible_parent_sets_enumerated(g, x)
    raise ValueError(f"unknown method {method!r}")


# --- the ambiguity set and the bias functional ----------------------------


@dataclass(frozen=True)
class BiasAt:
    """``B`` at one knowledge state, with the ambiguity set that produced it.

    Attributes:
        worst: ``B(G)``, the worst-case absolute bias. Reported as exactly
            ``0.0`` when below :data:`ZERO_TOL`.
        tau_min: Smallest effect the state leaves open.
        tau_max: Largest effect the state leaves open.
        n_parent_sets: Size of ``{pa_D(X) : D in [G]}``; 1 means the state pins
            the effect down completely.
        argmax_parents: A parent set attaining the worst case -- the witness,
            kept so the number can be replayed and checked independently.
    """

    worst: float
    tau_min: float
    tau_max: float
    n_parent_sets: int
    argmax_parents: frozenset[Node]

    @property
    def interval_width(self) -> float:
        """``tau_max - tau_min``: how much the state leaves undetermined."""
        return self.tau_max - self.tau_min


def bias_at(ctx: BiasContext, g: MPDAG, *, method: str = "semilocal") -> BiasAt:
    """Evaluate ``B`` at one knowledge state.

    Args:
        ctx: The per-instance context.
        g: The knowledge state, treated as a hypothetical constraint on the
            truth -- **not** as a truth in its own right, which is the
            difference from ``bias_stats``.
        method: Passed to :func:`possible_parent_sets`.

    Returns:
        A :class:`BiasAt`.

    Raises:
        ValueError: If ``[g]`` is empty, which no element of the space is
            (``THEOREMS.md`` section 2) and which would make the maximum
            undefined rather than vacuously anything.
    """
    parent_sets = possible_parent_sets(ctx.cpdag, g, ctx.x, method=method)
    if not parent_sets:
        raise ValueError(f"no DAG extension for state {g.edge_string()!r}; B is undefined")

    # sorted(): floating-point max is order-independent, but the tie-broken
    # argmax is not, and the witness is reported.
    taus = [(ctx.total_effect_under(p), p) for p in sorted(parent_sets, key=sorted)]
    values = [t for t, _ in taus]

    worst, witness = max(
        ((abs(ctx.theta_z - t), p) for t, p in taus),
        key=lambda pair: pair[0],
    )
    return BiasAt(
        worst=0.0 if worst < ZERO_TOL else float(worst),
        tau_min=float(min(values)),
        tau_max=float(max(values)),
        n_parent_sets=len(parent_sets),
        argmax_parents=witness,
    )


def worst_case_bias(ctx: BiasContext, g: MPDAG, *, method: str = "semilocal") -> float:
    """``B(G)`` alone, for use as a predicate inside a search.

    Args:
        ctx: The per-instance context.
        g: The knowledge state.
        method: Passed to :func:`possible_parent_sets`.

    Returns:
        The worst-case absolute bias.
    """
    return bias_at(ctx, g, method=method).worst


def realised_error(ctx: BiasContext) -> float:
    """``|theta_Z - tau_true|``: the error the analyst actually carries.

    Never used to compute a certificate -- only to audit one, which is the
    role the true DAG is allowed to play in this project.
    """
    return abs(ctx.theta_z - ctx.tau_true)


def normalise(ctx: BiasContext, value: float, *, units: str = "absolute") -> float:
    """Put a bias on a reportable scale.

    Every scale here divides by a constant of the instance, never by anything
    that varies with the knowledge state. That is what preserves monotonicity
    (``docs/R_EPSILON_THEORY.md``, Remark 1); a "relative bias" normalised by
    the state's own ambiguity set would look natural and destroy the
    certificate.

    Args:
        ctx: The per-instance context.
        value: A bias in effect units.
        units: ``"absolute"``, ``"relative"`` (fraction of ``|theta_Z|``, the
            reported estimate -- this is the one that reads as a percentage),
            or ``"standardised"`` (units of ``sd(Y)/sd(X)``).

    Returns:
        The rescaled value.

    Raises:
        ValueError: On an unknown unit, or on ``"relative"`` when the reported
            estimate is zero, where the ratio is not defined and returning
            infinity would quietly propagate into a radius.
    """
    if units == "absolute":
        return value
    if units == "standardised":
        return value * ctx.scale_std
    if units == "relative":
        if abs(ctx.theta_z) < ZERO_TOL:
            raise ValueError(
                "relative units are undefined when the reported estimate is zero; "
                "use absolute or standardised units for this instance"
            )
        return value / abs(ctx.theta_z)
    raise ValueError(f"unknown units {units!r}")
