"""Probes: does the learned geometry predict delta, or separate true from false?

Two questions, and they are not equally likely to have interesting answers.

**Does the geometry predict the breakdown radius?** Plausible. The radius is a
function of the knowledge and the CPDAG, both of which the embedding sees, so
the information is present in principle and the only question is whether the
geometry organises it usefully. Bear in mind that the radius is *computable*
exactly on the relevant graph sizes, so a positive result here is a claim about
speed or about generalisation to larger graphs, not about access to something
previously unavailable. Say which one is being claimed.

**Does the geometry separate consistent-and-true from consistent-but-false
knowledge?** This one deserves suspicion. Consistent-but-false knowledge is
defined by being indistinguishable from true knowledge given the CPDAG -- that
is what consistency means. So a probe that separates them cannot be detecting
truth; it must be detecting something about how the two arms were *generated*.
The sampler draws false knowledge by perturbing away from the truth and true
knowledge by reading off the truth, and those two procedures leave different
structural fingerprints even when the results are equally consistent.

:func:`shuffle_control` exists for exactly this. It re-labels the arms at random
and re-runs the probe; a probe that still separates is reading the sampler.
Report the control alongside the probe, always, and treat a high separation
score with a clean control as a hypothesis about the sampler until proven
otherwise.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    import numpy as np


@dataclass(frozen=True)
class ProbeResult:
    """Outcome of one probe, with the baselines it must beat.

    Attributes:
        probe: Which probe -- ``"predict_delta"``, ``"separate_true_false"``,
            ``"predict_bias"``.
        metric: Metric name -- ``"r2"``, ``"auc"``, ``"rmse"``.
        score: Probe score on the test split.
        baseline_scores: Structural baselines, keyed by name. The probe is only
            interesting if it beats all of them.
        control_score: Score under :func:`shuffle_control`. Should sit at chance;
            anything else means the probe is reading the data-generating
            procedure rather than the quantity of interest.
        n_train: Training set size.
        n_test: Test set size.
        beats_baselines: Whether ``score`` exceeds every baseline by more than
            the bootstrap standard error.
    """

    probe: str
    metric: str
    score: float
    baseline_scores: dict[str, float]
    control_score: float | None
    n_train: int
    n_test: int
    beats_baselines: bool


def predict_delta(
    embeddings: np.ndarray,
    deltas: np.ndarray,
    *,
    groups: np.ndarray | None = None,
    baselines: dict[str, np.ndarray] | None = None,
    probe: str = "ridge",
    seed: int = 0,
) -> ProbeResult:
    """Regress the breakdown radius on the embedding.

    Args:
        embeddings: ``(n_samples, dim)`` embedded knowledge sets.
        deltas: ``(n_samples,)`` radii.
        groups: Graph identifiers for the split. Splitting by graph rather than
            by sample is mandatory: knowledge sets drawn from the same graph
            share its structure, so a random split leaks and the score is
            meaningless.
        baselines: Structural features to compare against, keyed by name.
        probe: ``"ridge"``, ``"mlp"`` or ``"knn"``. Keep it weak -- a strong
            probe measures its own capacity, not the representation's.
        seed: Seed for the split and the probe.

    Returns:
        A :class:`ProbeResult` with ``metric="r2"``.
    """
    raise NotImplementedError


def separate_true_false(
    embeddings: np.ndarray,
    labels: np.ndarray,
    *,
    groups: np.ndarray | None = None,
    baselines: dict[str, np.ndarray] | None = None,
    probe: str = "logistic",
    seed: int = 0,
) -> ProbeResult:
    """Classify consistent-and-true against consistent-but-false knowledge.

    Args:
        embeddings: ``(n_samples, dim)`` embedded knowledge sets.
        labels: ``(n_samples,)``; ``1`` for false, ``0`` for true.
        groups: Graph identifiers for the split.
        baselines: Structural features to compare against.
        probe: ``"logistic"``, ``"mlp"`` or ``"knn"``.
        seed: Seed.

    Returns:
        A :class:`ProbeResult` with ``metric="auc"``.

    Note:
        Always populate ``control_score`` via :func:`shuffle_control` for this
        probe. See the module docstring for why a positive result here is more
        likely to be about the sampler than about truth.
    """
    raise NotImplementedError


def predict_bias(
    embeddings: np.ndarray,
    biases: np.ndarray,
    *,
    groups: np.ndarray | None = None,
    baselines: dict[str, np.ndarray] | None = None,
    probe: str = "ridge",
    seed: int = 0,
) -> ProbeResult:
    """Regress the realised estimation bias on the embedding.

    The most practically useful of the three if it works: it would let an
    analyst rank knowledge sets by expected damage without knowing the truth.
    Also the hardest, since the bias depends on the SCM and not on the graph
    alone, so the embedding is being asked for something it does not see.

    Returns:
        A :class:`ProbeResult` with ``metric="r2"``.
    """
    raise NotImplementedError


def shuffle_control(
    probe_fn: Any,
    embeddings: np.ndarray,
    targets: np.ndarray,
    *,
    n_shuffles: int = 20,
    seed: int = 0,
    **kwargs: Any,
) -> float:
    """Re-run a probe on shuffled targets, to measure what it scores by construction.

    Args:
        probe_fn: One of the probe functions above.
        embeddings: The embeddings.
        targets: The true targets, which are permuted.
        n_shuffles: Permutations to average over.
        seed: Seed.
        **kwargs: Passed through to ``probe_fn``.

    Returns:
        Mean score across permutations. A probe scoring well above chance here
        is exploiting structure in the sampling, and its headline score should
        not be reported without this number beside it.
    """
    raise NotImplementedError


def baseline_features(
    knowledge_sets: Sequence[Any],
    cpdags: Sequence[Any],
    names: Sequence[str] = ("constraint_count", "cascade_size", "shd_to_true_cpdag"),
) -> dict[str, np.ndarray]:
    """Compute the structural baselines a learned representation must beat.

    Args:
        knowledge_sets: The sampled knowledge sets.
        cpdags: The CPDAG each was sampled against.
        names: Which features to compute.

    Returns:
        Feature arrays keyed by name.

    Raises:
        KeyError: If a name is not a known feature.
    """
    raise NotImplementedError
