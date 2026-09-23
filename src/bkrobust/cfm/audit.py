"""The bias-versus-delta sweep over foundation-model checkpoints.

Holds the data fixed, varies the knowledge, and watches the estimate move.

Three arms are needed for the comparison to mean anything:

1. **Unconditioned.** The checkpoint with no knowledge at all -- the baseline
   error attributable to the model rather than to the knowledge.
2. **True knowledge.** Conditioned on constraints the true DAG satisfies. If
   this arm is no better than the baseline, the model is ignoring the
   conditioning signal, and every result in the false arm is measuring noise.
   Establishing that the signal does something is a precondition for the audit,
   not an optional extra.
3. **Consistent-but-false knowledge.** The object of study, swept over delta.

Alongside these, a classical arm -- OLS on ``O*`` from the same perturbed MPDAG
-- runs on identical data, so the amortized and explicit pipelines can be
compared at each radius.

Two sanity gates run before any sweep, and both abort rather than warn:

* :func:`~bkrobust.cfm.conditioning.verify_zero_scale_identity` -- a zero-scale
  bias must not change the forward pass;
* the true-knowledge arm must beat the unconditioned arm.

Without gate two, a flat bias-versus-delta curve is uninterpretable: it looks
like robustness and is indistinguishable from the model not reading the input.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    import numpy as np
    import pandas as pd

    from bkrobust.graphs.mpdag import MPDAG, Node


@dataclass(frozen=True)
class AuditResult:
    """One checkpoint's behaviour at one radius under one conditioning mode.

    Attributes:
        checkpoint: Registry key.
        revision: Checkpoint revision, carried through so a row can always be
            traced to a specific artefact.
        conditioning_mode: How knowledge was injected.
        bias_scale: Injection magnitude.
        delta: Perturbation radius.
        arm: ``"unconditioned"``, ``"true"`` or ``"false"``.
        ate_bias: Mean signed bias of the ATE estimate.
        ate_abs_bias: Mean absolute bias.
        pehe: Root PEHE, where the checkpoint supports CATE; ``None`` otherwise.
        classical_ate_bias: Bias of the OLS-on-``O*`` arm on the same data, for
            side-by-side reading.
        agreement: Correlation between the amortized and classical estimates
            across replicates. A model that tracks the classical estimator is
            inheriting the same structural failure; one that does not is failing
            -- or succeeding -- some other way, and the distinction is the
            headline result of this section.
        n_replicates: Replicates behind these numbers.
        diagnostics: Per-checkpoint extras.
    """

    checkpoint: str
    revision: str
    conditioning_mode: str
    bias_scale: float
    delta: int
    arm: str
    ate_bias: float
    ate_abs_bias: float
    pehe: float | None
    classical_ate_bias: float | None
    agreement: float | None
    n_replicates: int
    diagnostics: dict[str, Any]


def audit_checkpoint(
    checkpoint: str,
    data: pd.DataFrame,
    cpdag: MPDAG,
    true_graph: MPDAG,
    treatment: Node,
    outcome: Node,
    true_ate: float,
    delta_grid: Iterable[int],
    rng: np.random.Generator,
    *,
    conditioning_modes: Sequence[str] = ("none", "attention_bias"),
    bias_scale_grid: Sequence[float] = (0.0, 1.0),
    n_replicates: int = 10,
    reference_estimator: str | None = "ols_adjustment",
) -> list[AuditResult]:
    """Sweep one checkpoint over radii, conditioning modes and bias scales.

    Args:
        checkpoint: Registry key.
        data: The dataset, held fixed across the whole sweep.
        cpdag: The CPDAG knowledge is sampled against.
        true_graph: The ground-truth DAG.
        treatment: Treatment column.
        outcome: Outcome column.
        true_ate: The ground-truth effect, for the bias computation.
        delta_grid: Radii to sweep.
        rng: Seeded generator for the knowledge sampler.
        conditioning_modes: Injection modes to compare.
        bias_scale_grid: Injection magnitudes to compare.
        n_replicates: Knowledge draws per cell.
        reference_estimator: Classical estimator to run alongside; ``None``
            skips the comparison.

    Returns:
        One :class:`AuditResult` per cell.

    Raises:
        ImportError: If the ``[cfm]`` extra is not installed.
        RuntimeError: If either sanity gate fails -- zero-scale identity, or the
            true-knowledge arm failing to improve on the unconditioned arm.
    """
    raise NotImplementedError


def audit_all(
    checkpoints: Sequence[str],
    *args: Any,
    **kwargs: Any,
) -> list[AuditResult]:
    """Run :func:`audit_checkpoint` over several checkpoints.

    Returns:
        The concatenated results.
    """
    raise NotImplementedError


def sanity_check_conditioning(
    checkpoint: str,
    data: pd.DataFrame,
    cpdag: MPDAG,
    true_graph: MPDAG,
    treatment: Node,
    outcome: Node,
    true_ate: float,
    rng: np.random.Generator,
) -> dict[str, Any]:
    """Run both gates and report what they found.

    Args:
        checkpoint: Registry key.
        data: The dataset.
        cpdag: The CPDAG.
        true_graph: The ground-truth DAG.
        treatment: Treatment column.
        outcome: Outcome column.
        true_ate: Ground-truth effect.
        rng: Seeded generator.

    Returns:
        Keys ``"zero_scale_identity"`` (bool), ``"true_knowledge_helps"`` (bool),
        ``"unconditioned_abs_bias"``, ``"true_knowledge_abs_bias"``, and
        ``"verdict"`` -- ``"ok"`` or a description of which gate failed.

    Note:
        Report these numbers in the paper whatever they say. A checkpoint that
        ignores its conditioning input is a finding about that checkpoint, not
        an inconvenience to be dropped from the table.
    """
    raise NotImplementedError


def results_to_frame(results: Iterable[AuditResult]) -> pd.DataFrame:
    """Flatten audit results into a DataFrame, one row per result."""
    raise NotImplementedError
