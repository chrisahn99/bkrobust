"""Breakdown radii, bias bounds, the efficiency gap, and the certificate.

The theoretical core. Given a CPDAG, a target pair ``(X, Y)`` and the true DAG,
two radii describe how much wrong background knowledge the optimal adjustment
set ``O*`` can absorb:

``delta_opt``
    The smallest perturbation at which ``O*`` computed on the perturbed MPDAG
    stops being optimal. Past it the estimate is still unbiased; it just costs
    more data.
``delta_valid``
    The smallest perturbation at which ``O*`` stops being valid. Past it the
    estimate is biased, and no sample size helps.

The conjectured ordering is ``delta_opt <= delta_valid``: optimality is the more
fragile property, so it should be the first to go. The band between them is a
regime worth naming, because in it the analyst is wrong, is paying for it, and
has no diagnostic that would tell them so.

Module map:

``radius``
    Computes both radii, exactly on small graphs and heuristically on large ones.
``bias_bound``
    How large the bias can be at and beyond ``delta_valid``.
``efficiency_gap``
    How much variance is lost inside the band.
``certificate``
    Packages all of it for someone who has knowledge, data, and no ground truth.
"""

from __future__ import annotations

__all__: list[str] = []
