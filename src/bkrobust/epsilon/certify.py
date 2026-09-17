"""The end-to-end call: from an analyst's inputs to a banded bias certificate.

One function, :func:`certify`, takes what an analyst has -- an estimated CPDAG,
the orientations they asserted, a query, and a covariance matrix -- and returns
the statement the objective for this session asked for:

    *If at most k of your claims are wrong, your reported effect is off by at
    most eps_y, and a bias above eps_x is already reachable.*

The chain behind it, each link proved in ``docs/R_EPSILON_THEORY.md``:

1. ``r_val`` comes from the existing accelerated hybrid
   (:func:`bkrobust.hybrid.breakdown_radius`) -- unchanged, and the paper's own
   method.
2. Theorem A makes the whole ball below ``r_val`` algebraically bias-free, so no
   bias is evaluated there. The certificate's cheapest region is the one that
   matters most.
3. Theorem C builds the staircase ``beta_up`` from ``r_val`` outwards, and every
   ``r_eps`` on the grid is then a lookup on it -- one traversal, all
   thresholds.
4. Theorem D turns two adjacent radii into the band.

The true DAG is never on the computation path. It may be supplied in simulation,
in which case the realised error is recorded next to the certificate so that the
guarantee can be *audited* rather than merely asserted -- the same discipline the
paper applies to ``r_val`` ("the data-generating DAG is never on the computation
path; it is used to draw samples and, afterwards, to check whether the
certificate held").
"""

from __future__ import annotations

import time
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from bkrobust.core.conventions import UNREACHED
from bkrobust.demo.evaluate import LinearSEM, optimal_adjustment_set_mpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.epsilon.bias import (
    BiasContext,
    make_context,
    make_context_from_covariance,
    normalise,
)
from bkrobust.epsilon.profile import (
    DEFAULT_SHELL_CAP,
    ShellBias,
    certified_band,
    epsilon_grid,
    top_bias,
)
from bkrobust.hybrid import breakdown_radius

Node = str
Edge = tuple[str, str]

#: Default thresholds, as fractions of the reported estimate: 1%, 5%, 10%, 25%,
#: 50%, 100%. Chosen to straddle the range over which a practitioner's verdict
#: would plausibly change; the grid is an argument precisely because that range
#: is domain-specific and should not be a library's opinion.
DEFAULT_EPSILONS: tuple[float, ...] = (0.01, 0.05, 0.10, 0.25, 0.50, 1.00)


@dataclass
class Certificate:
    """A banded bias certificate for one query under one knowledge state.

    Attributes:
        r_val: The validity radius, from the accelerated hybrid. Shells
            ``0 .. r_val - 1`` carry exactly zero identification bias.
        r_eps: Threshold to radius. :data:`~bkrobust.core.conventions.UNREACHED`
            means no perturbation reaches that bias -- a status, not a number,
            and never to be averaged.
        units: ``"relative"``, ``"absolute"`` or ``"standardised"``.
        epsilons: The grid, in ``units``.
        theta_z: The estimand actually reported.
        profile: The staircase the radii were read off.
        beta_top: ``B(Chat)``, the largest bias any perturbation can produce --
            the ceiling on the whole certificate.
        z: The adjustment set the certificate is about.
        g0: The analyst's state, as an edge string.
        n_knowledge: ``|K_{G0}|``, the number of orientations the certificate's
            error budget is denominated in.
        realised_error: ``|theta_Z - tau_true|`` in ``units`` when a truth was
            supplied; ``nan`` otherwise. Audit only.
        eps_scale: Effect units per reported unit. ``profile`` is stored in
            **effect** units, because that is what the traversal computes;
            everything the certificate *reports* is in ``units``, and this is the
            single conversion between them. Keeping one factor rather than
            converting the staircase in place is what makes it possible to
            re-read the same staircase in different units without recomputing.
        method: How ``r_val`` was obtained.
        seconds_r_val: Time in the validity search.
        seconds_r_eps: Time in the bias traversal.
        assumes: The inherited assumption chain.
        status: ``"ok"``, or a reason the certificate is partial.
    """

    r_val: int
    r_eps: dict[float, int]
    units: str
    epsilons: tuple[float, ...]
    theta_z: float
    profile: list[ShellBias]
    beta_top: float
    z: tuple[Node, ...]
    g0: str
    n_knowledge: int
    realised_error: float = float("nan")
    method: str = ""
    seconds_r_val: float = 0.0
    seconds_r_eps: float = 0.0
    assumes: str = (
        "Lemma L of THEOREMS.md (hence Anti-Exchange Case B, verified not proved), "
        "exactly as r_val does; error is one-sided towards larger radii"
    )
    status: str = "ok"
    eps_scale: float = 1.0
    _ctx: BiasContext | None = field(default=None, repr=False)

    def band(self, k: int) -> tuple[float | None, float | None]:
        """The band ``(eps_x, eps_y]`` the bias is pinned into at error budget ``k``.

        Args:
            k: How many of the analyst's asserted orientations might be false.

        Returns:
            ``(lower, upper)`` as in
            :func:`~bkrobust.epsilon.profile.certified_band`. ``lower is None``
            means no threshold on the grid is reachable within ``k``; ``upper is
            None`` means the grid does not bound the bias from above at ``k``
            and a coarser grid is needed to say anything.
        """
        return certified_band(self.r_eps, k)

    def worst_case_at(self, k: int) -> float:
        """``beta_up(k)``: the exact worst-case bias at error budget ``k``.

        The band is the reportable summary; this is the underlying number, read
        straight off the staircase and converted into ``units``. Depths beyond
        the traversal return the largest value the traversal established, which
        is a lower bound on the truth, so callers must not read it as exact past
        ``max(s.d for s in profile)``.

        Args:
            k: The error budget, in orientations of ``K_{G0}``.

        Returns:
            ``beta_up(k)`` in the certificate's reporting ``units`` -- the same
            units as :attr:`realised_error` and :attr:`beta_top`, so the three
            are directly comparable. (They were not, once: the staircase is
            computed in effect units and reporting it unconverted made an
            audit look like a violation.)
        """
        best = 0.0
        for step in self.profile:
            if step.d <= k:
                best = max(best, step.beta_up)
        return best / self.eps_scale if self.eps_scale else best


def certify(
    cpdag: MPDAG,
    knowledge: Sequence[Edge] | None,
    x: Node,
    y: Node,
    *,
    sem: LinearSEM | None = None,
    sigma: np.ndarray | None = None,
    nodes: Sequence[Node] | None = None,
    z: frozenset[Node] | None = None,
    g0: MPDAG | None = None,
    epsilons: Sequence[float] = DEFAULT_EPSILONS,
    units: str = "relative",
    method: str = "semilocal",
    cap: int = DEFAULT_SHELL_CAP,
    search_budget: int = 3,
) -> Certificate:
    """Compute the banded bias certificate for one query.

    Args:
        cpdag: The estimated CPDAG.
        knowledge: The asserted orientations. Ignored when ``g0`` is supplied.
        x: Treatment.
        y: Outcome.
        sem: The true SCM -- the simulation path. Supplies the covariance and,
            for auditing only, the true effect.
        sigma: Covariance matrix -- the analyst path. Requires ``nodes``.
        nodes: Node order for ``sigma``.
        z: The adjustment set. Defaults to the optimal set read off ``G0``,
            which is the paper's convention.
        g0: The analyst's state, if already built.
        epsilons: Thresholds, in ``units``.
        units: ``"relative"`` (fraction of the reported estimate; reads as a
            percentage), ``"absolute"``, or ``"standardised"``.
        method: ``"semilocal"`` (Theorem E) or ``"enumerate"``.
        cap: Shell-size cap, passed through.
        search_budget: Depth budget for the ``r_val`` hybrid.

    Returns:
        A :class:`Certificate`.

    Raises:
        ValueError: If neither ``sem`` nor ``(sigma, nodes)`` is supplied, if
            the knowledge is inconsistent with the CPDAG, or if no adjustment
            set is identified at ``G0``.
    """
    if g0 is None:
        built = apply_orientations(cpdag, list(knowledge or []))
        if built is None:
            raise ValueError("knowledge is inconsistent with the CPDAG; no G0 exists")
        g0 = built

    if z is None:
        derived = optimal_adjustment_set_mpdag(g0, x, y)
        if derived is None:
            raise ValueError("no adjustment set is identified at G0; nothing to certify")
        z = frozenset(derived)

    if sem is not None:
        ctx = make_context(sem, cpdag, x, y, z)
    elif sigma is not None and nodes is not None:
        ctx = make_context_from_covariance(sigma, nodes, cpdag, x, y, z)
    else:
        raise ValueError("supply either sem= or both sigma= and nodes=")

    # r_val: the paper's own method, unchanged.
    hybrid = breakdown_radius(cpdag, None, x, y, z, g0=g0, search_budget=search_budget)
    r_val = hybrid.radius if hybrid.radius != UNREACHED else 0

    # Thresholds arrive in report units; the staircase is built in effect units,
    # so convert the grid once rather than converting every bias.
    if units == "relative":
        scale = abs(ctx.theta_z)
        if scale < 1e-12:
            raise ValueError(
                "relative units are undefined: the reported estimate is zero. "
                "Use units='absolute' or units='standardised' for this instance."
            )
    elif units == "standardised":
        scale = 1.0 / ctx.scale_std if ctx.scale_std > 0 else 1.0
    elif units == "absolute":
        scale = 1.0
    else:
        raise ValueError(f"unknown units {units!r}")

    in_effect_units = [float(e) * scale for e in epsilons]

    t0 = time.perf_counter()
    radii_effect, profile = epsilon_grid(
        ctx, g0, in_effect_units, r_val=r_val, method=method, cap=cap
    )
    seconds_eps = time.perf_counter() - t0

    r_eps = {
        float(reported): radii_effect[effect]
        for reported, effect in zip(epsilons, in_effect_units, strict=True)
    }
    top = top_bias(ctx, method=method)
    n_knowledge = len(set(g0.directed_edges) - set(cpdag.directed_edges))

    return Certificate(
        r_val=hybrid.radius,
        r_eps=r_eps,
        units=units,
        epsilons=tuple(float(e) for e in epsilons),
        theta_z=ctx.theta_z,
        profile=profile,
        beta_top=normalise(ctx, top.worst, units=units),
        z=tuple(sorted(z)),
        g0=g0.edge_string(),
        n_knowledge=n_knowledge,
        realised_error=(
            float("nan")
            if np.isnan(ctx.tau_true)
            else normalise(ctx, abs(ctx.theta_z - ctx.tau_true), units=units)
        ),
        method=hybrid.method,
        seconds_r_val=hybrid.total_seconds,
        seconds_r_eps=seconds_eps,
        status="ok" if hybrid.exact else "r_val_not_exact",
        eps_scale=scale,
        _ctx=ctx,
    )


def describe(cert: Certificate, k: int | None = None) -> str:
    """A practitioner-readable rendering of the certificate.

    Args:
        cert: What :func:`certify` returned.
        k: The error budget to report the band at. Defaults to ``r_val``, the
            first budget at which anything can go wrong at all.

    Returns:
        A short paragraph stating the zero-bias guarantee, the band, and the
        assumption the numbers carry.
    """
    if cert.status != "ok":
        return (
            "The validity radius did not complete exactly, so no bias certificate "
            "is available; the partial numbers must not be read as a guarantee."
        )
    pct = cert.units == "relative"

    def fmt(e: float) -> str:
        return f"{e:.0%}" if pct else f"{e:g}"

    if cert.r_val == UNREACHED:
        head = (
            "No perturbation of the asserted knowledge invalidates this adjustment "
            "set, so its identification bias is exactly zero throughout the space."
        )
        return f"{head} (Assumes: {cert.assumes}.)"

    budget = cert.r_val if k is None else k
    lower, upper = cert.band(budget)
    head = (
        f"The reported estimate is {cert.theta_z:.4g}. It carries exactly zero "
        f"identification bias as long as at most {cert.r_val - 1} of the "
        f"{cert.n_knowledge} asserted orientations are wrong."
    )
    if lower is None and upper is None:
        body = ""
    elif lower is None:
        body = (
            f" At a budget of {budget} wrong orientations the bias is still certified "
            f"at or below {fmt(upper)}."
        )
    elif upper is None:
        body = (
            f" At a budget of {budget} wrong orientations a bias above {fmt(lower)} is "
            "reachable, and the supplied grid does not bound it from above."
        )
    else:
        body = (
            f" At a budget of {budget} wrong orientations the bias is at most "
            f"{fmt(upper)}, and a bias above {fmt(lower)} is reachable."
        )
    ceiling = f" No perturbation whatsoever can push it past {fmt(cert.beta_top)}."
    return f"{head}{body}{ceiling} (Assumes: {cert.assumes}.)"
