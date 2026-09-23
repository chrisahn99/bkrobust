"""How large the bias can be once ``O*`` stops being valid.

Past ``delta_valid`` the adjustment set read off the perturbed MPDAG no longer
identifies the effect. Knowing that is a qualitative statement; a practitioner
needs the quantitative one -- how wrong can the answer be?

In the linear-Gaussian case the question has a clean answer with a familiar
shape. An invalid adjustment set errs in one of two ways: it omits a
confounder, or it includes a node it should not (a collider, or a descendant of
the treatment, opening a path that adjustment then fails to block). Either way
the bias factors into two pieces, exactly as in classical omitted-variable-bias
sensitivity analysis: how strongly the mishandled variable relates to the
treatment, and how strongly it relates to the outcome given the rest.

That factoring is what makes the bound usable. Both pieces are partial
:math:`R^2` quantities, so they can be reasoned about on a bounded scale even
by someone with no access to the true graph -- which is precisely the situation
of the analyst this project is written for. The certificate in
:mod:`bkrobust.theory.certificate` reports the bound in exactly these terms.

Beyond linear-Gaussian the factoring does not survive in closed form, and
:func:`worst_case_bias` falls back to numerical maximisation over the SCMs
compatible with the graph. See ``docs/THEORY.md`` (target **T2**).
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    from bkrobust.graphs.mpdag import MPDAG, Node


@dataclass(frozen=True)
class BiasBound:
    """A bias bound with the sensitivity parameters that produced it.

    Attributes:
        bound: Upper bound on ``|E[tau_hat] - tau|``.
        realised: The actual bias, when an SCM was supplied; ``None`` otherwise.
        r2_treatment: Partial :math:`R^2` of the mishandled variables with the
            treatment, given the rest of the adjustment set.
        r2_outcome: Partial :math:`R^2` of the mishandled variables with the
            outcome, given treatment and the rest of the set.
        mechanism: ``"omitted_confounder"``, ``"included_collider"``,
            ``"included_descendant"``, or ``"mixed"``.
        offending_nodes: The nodes responsible.
        tight: Whether the bound is attained by some SCM compatible with the
            graph. A loose bound is still sound but says less, and the
            distinction belongs in any table that reports one.
    """

    bound: float
    realised: float | None
    r2_treatment: float
    r2_outcome: float
    mechanism: str
    offending_nodes: frozenset[Node]
    tight: bool


def linear_gaussian_bias(
    true_graph: MPDAG,
    treatment: Node,
    outcome: Node,
    z: Iterable[Node],
    scm: Any,
) -> float:
    """Exact asymptotic bias of the ``z``-adjusted estimator under a known SCM.

    The reference quantity every bound is checked against: with the SCM in hand
    the bias is computable, not merely boundable, and any bound that a sampled
    SCM violates is wrong.

    Args:
        true_graph: The ground-truth DAG.
        treatment: Treatment node.
        outcome: Outcome node.
        z: The adjustment set actually used, valid or not.
        scm: The linear-Gaussian SCM.

    Returns:
        ``E[tau_hat] - tau``, signed. Zero iff ``z`` is valid in ``true_graph``.
    """
    raise NotImplementedError


def worst_case_bias(
    cpdag: MPDAG,
    true_graph: MPDAG,
    treatment: Node,
    outcome: Node,
    delta: int,
    *,
    scm: Any | None = None,
    r2_treatment: float | None = None,
    r2_outcome: float | None = None,
    distance: str = "orientation_flips",
) -> BiasBound:
    """Largest bias achievable by consistent-but-false knowledge at radius ``delta``.

    Maximises over the radius-``delta`` ball of knowledge sets and, when no SCM
    is supplied, over the SCMs compatible with the graph subject to the
    sensitivity parameters.

    Args:
        cpdag: The CPDAG the knowledge is imposed on.
        true_graph: The ground-truth DAG.
        treatment: Treatment node.
        outcome: Outcome node.
        delta: The perturbation radius.
        scm: When given, the bound is evaluated on this SCM alone and is exact.
        r2_treatment: Assumed bound on the treatment-side partial :math:`R^2`.
            Required when ``scm`` is ``None``.
        r2_outcome: Assumed bound on the outcome-side partial :math:`R^2`.
            Required when ``scm`` is ``None``.
        distance: Metric the radius is measured in.

    Returns:
        A :class:`BiasBound`.

    Raises:
        ValueError: If ``scm`` is ``None`` and either sensitivity parameter is
            missing -- there is no bound without an assumption, and defaulting
            one silently would produce a number that looks like a guarantee and
            is not.
    """
    raise NotImplementedError


def omitted_variable_bias(
    r2_treatment: float,
    r2_outcome: float,
    se_tau: float,
    df: int,
) -> float:
    """Classical OVB sensitivity bound, in the Cinelli-Hazlett parameterisation.

    Bias attributable to omitting variables that explain ``r2_treatment`` of the
    residual treatment variation and ``r2_outcome`` of the residual outcome
    variation, expressed in units of the estimate's own standard error so that
    it can be read without knowing the effect scale.

    Args:
        r2_treatment: Partial :math:`R^2` with the treatment, in ``[0, 1)``.
        r2_outcome: Partial :math:`R^2` with the outcome, in ``[0, 1)``.
        se_tau: Standard error of the estimated effect.
        df: Residual degrees of freedom.

    Returns:
        The bias bound on the effect scale.

    Raises:
        ValueError: If either :math:`R^2` is outside ``[0, 1)``.
    """
    raise NotImplementedError


def bias_at_radius_curve(
    cpdag: MPDAG,
    true_graph: MPDAG,
    treatment: Node,
    outcome: Node,
    scm: Any,
    delta_grid: Iterable[int],
    *,
    distance: str = "orientation_flips",
) -> dict[int, BiasBound]:
    """Worst-case bias at each radius in ``delta_grid``.

    The main figure of the estimation-level section: flat at zero below
    ``delta_valid``, then rising. Whether the rise is a step or a ramp is an
    empirical question this function is built to answer, and the shape matters
    -- a step means the failure is all-or-nothing, a ramp means partial
    misspecification carries partial cost.

    Args:
        cpdag: The CPDAG the knowledge is imposed on.
        true_graph: The ground-truth DAG.
        treatment: Treatment node.
        outcome: Outcome node.
        scm: The SCM.
        delta_grid: Radii to evaluate.
        distance: Metric the radii are measured in.

    Returns:
        Mapping from radius to bound.
    """
    raise NotImplementedError
