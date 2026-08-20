"""Data generation and loading.

Three tiers, trading ground truth against realism:

``synthetic``
    Random DAG, sampled SCM, simulated data. The true graph is known, so both
    breakdown radii are computable exactly. Where the theory is tested.
``semisynthetic``
    Real covariates, simulated outcomes. The true outcome mechanism is known
    while the covariate distribution is not synthetic. Where the theory is
    tested against realistic marginals.
``loaders``
    Real benchmarks. No true graph, so no radii -- only the certificate, which
    is what a practitioner would have.

``leakage`` cuts across all three. Simulated data leaks its causal order through
marginal variances -- varsortability -- and a method that exploits the leak
looks good for the wrong reason. Every synthetic and semi-synthetic result must
be reported with its varsortability alongside.
"""

from __future__ import annotations

__all__: list[str] = []
