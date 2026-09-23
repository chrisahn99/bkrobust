"""Varsortability and R2-sortability diagnostics.

Simulated data leaks its causal order. If each node's variance is the sum of its
parents' contributions plus its own noise, variance tends to grow along the
causal order, and a method that sorts variables by marginal variance recovers
much of the DAG without doing any causal reasoning at all. Reisach et al. named
this varsortability and showed that on standard synthetic benchmarks it is high
enough for the trivial sorting baseline to be competitive with published methods.

Why it matters here specifically. The perturbation experiments claim that wrong
background knowledge changes ``O*`` and thereby changes the estimate. If the data
is highly varsortable, discovery recovers the true order from the variances no
matter what the knowledge says, the knowledge never gets to matter, and the
measured radii are large for a reason that has nothing to do with robustness.
That would be a benchmark artefact reported as a finding.

So every synthetic and semi-synthetic result carries its varsortability. Three
responses when it is high, in decreasing order of preference: sample noise
scales to break the pattern (which is why ``noise_scale`` is a range, not a
constant); report results split by varsortability; or standardise the data,
which is the worst option -- it hides the leak from the diagnostic without
removing the information, and it changes the estimand's scale.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    import pandas as pd

    from bkrobust.graphs.mpdag import MPDAG


def varsortability(data: pd.DataFrame, dag: MPDAG) -> float:
    """Fraction of directed edges whose marginal variances increase along the edge.

    Args:
        data: Observational data; columns are variable names.
        dag: The true DAG.

    Returns:
        A value in ``[0, 1]``. ``0.5`` is chance; values near ``1`` mean the
        causal order is readable from the marginal variances alone.

    Raises:
        ValueError: If ``data`` is missing a column for some node of ``dag``.
    """
    raise NotImplementedError


def r2_sortability(data: pd.DataFrame, dag: MPDAG) -> float:
    """Varsortability's scale-invariant analogue, using explained variance.

    Survives standardisation, which plain varsortability does not. If a dataset
    is standardised and this is still high, the leak is real and not an artefact
    of measurement units.

    Returns:
        A value in ``[0, 1]``; ``0.5`` is chance.
    """
    raise NotImplementedError


def sortability_report(data: pd.DataFrame, dag: MPDAG) -> dict[str, float]:
    """Compute both diagnostics plus a trivial-baseline score.

    Returns:
        Keys ``"varsortability"``, ``"r2_sortability"``,
        ``"sorting_baseline_shd"`` -- the SHD of the DAG obtained by sorting on
        variance alone, which is the concrete version of "how much of this graph
        comes free".
    """
    raise NotImplementedError


def check_leakage(
    data: pd.DataFrame,
    dag: MPDAG,
    *,
    threshold: float = 0.8,
    report_only: bool = True,
) -> dict[str, float]:
    """Run the diagnostics and optionally refuse to proceed.

    Args:
        data: Observational data.
        dag: The true DAG.
        threshold: Varsortability above this counts as leaking.
        report_only: Warn rather than raise. Default ``True`` so that a leaking
            configuration is still measured and reported -- suppressing the run
            would remove the evidence that the benchmark family has this
            property.

    Returns:
        The report from :func:`sortability_report`.

    Raises:
        ValueError: If the threshold is exceeded and ``report_only`` is ``False``.
    """
    raise NotImplementedError
