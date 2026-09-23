"""The epsilon-bias radius: a bias certificate on top of the validity radius.

``r_val`` says where the committed adjustment set stops being *valid*. It does
not say how wrong the answer becomes once it does, and the paper's own App. C
table shows that the two questions have different answers: past ``r_val`` the
bias is a ramp, not a step.

This package supplies the missing half. It defines a worst-case bias functional
``B`` that is exactly zero on the whole ball ``r_val`` certifies, monotone under
knowledge retraction, and an upper bound on the error the analyst actually
carries; and it computes, from one traversal, the radius ``r_eps`` at which
``B`` first exceeds any threshold. Pairing two thresholds pins the bias into a
band: *at most ``eps_y``, and demonstrably able to exceed ``eps_x``*.

The statements and proofs are in ``docs/R_EPSILON_THEORY.md``; the summary table
there says what each one assumes. The short version is that ``r_eps`` inherits
the assumption chain of ``r_val`` and adds nothing to it.

Entry points:

* :func:`~bkrobust.epsilon.bias.make_context` -- build the per-instance context.
* :func:`~bkrobust.epsilon.bias.bias_at` -- ``B`` at one state.
* :func:`~bkrobust.epsilon.profile.epsilon_grid` -- every ``r_eps`` at once.
* :func:`~bkrobust.epsilon.profile.r_epsilon` -- one radius, with dispatch.
* :func:`~bkrobust.epsilon.certify.certify` -- the end-to-end practitioner call.
"""

from bkrobust.epsilon.bias import (
    BiasAt,
    BiasContext,
    bias_at,
    make_context,
    make_context_from_covariance,
    normalise,
    possible_parent_sets,
    realised_error,
    worst_case_bias,
)
from bkrobust.epsilon.certify import Certificate, certify, describe
from bkrobust.epsilon.finite_sample import (
    RadiusDistribution,
    bootstrap_radius,
    sample_covariance,
    sample_data,
)
from bkrobust.epsilon.profile import (
    EpsilonResult,
    ShellBias,
    bias_profile,
    certified_band,
    epsilon_grid,
    greedy_chain_bound,
    r_epsilon,
    r_epsilon_from_profile,
    retraction_shell,
    top_bias,
)

__all__ = [
    "BiasAt",
    "BiasContext",
    "Certificate",
    "EpsilonResult",
    "RadiusDistribution",
    "ShellBias",
    "bias_at",
    "bias_profile",
    "bootstrap_radius",
    "certified_band",
    "certify",
    "describe",
    "epsilon_grid",
    "greedy_chain_bound",
    "make_context",
    "make_context_from_covariance",
    "normalise",
    "possible_parent_sets",
    "r_epsilon",
    "r_epsilon_from_profile",
    "realised_error",
    "retraction_shell",
    "sample_covariance",
    "sample_data",
    "top_bias",
    "worst_case_bias",
]
