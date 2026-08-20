"""CausalFM adapter.

A foundation model for causal effect estimation, pretrained across a large
family of synthetic SCMs. Where CausalPFN is purely in-context, CausalFM exposes
more of a structural interface, which makes it the checkpoint most likely to
accept an ancestral-matrix encoding directly rather than through an attention
bias -- and therefore the most informative one for separating "the model ignores
the knowledge" from "the model uses the knowledge and is hurt by it".

Requires the ``[cfm]`` extra. See ``docs/CHECKPOINTS.md``.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    import numpy as np
    import pandas as pd

    from bkrobust.estimation.base import EffectEstimate
    from bkrobust.graphs.mpdag import Node
    from bkrobust.knowledge.base import BackgroundKnowledge


class CausalFMAdapter:
    """Wrap a CausalFM checkpoint as an :class:`~bkrobust.estimation.base.Estimator`.

    Args:
        checkpoint: Registry key or local path; ``None`` uses the registry
            default.
        revision: Checkpoint revision. Required for any reported number.
        device: ``"cpu"``, ``"cuda"`` or ``"mps"``.
        dtype: Inference dtype.
        batch_size: Inference batch size.
        conditioning: The ``conditioning`` config block.

    Raises:
        ImportError: If the ``[cfm]`` extra is not installed.
    """

    name: str = "causalfm"
    requires_scm: bool = False

    def __init__(
        self,
        checkpoint: str | None = None,
        revision: str | None = None,
        device: str = "cpu",
        dtype: str = "float32",
        batch_size: int = 64,
        conditioning: dict[str, Any] | None = None,
    ) -> None:
        raise NotImplementedError

    def fit(
        self,
        data: pd.DataFrame,
        treatment: Node,
        outcome: Node,
        adjustment_set: Iterable[Node],
        *,
        knowledge: BackgroundKnowledge | None = None,
        **kwargs: Any,
    ) -> CausalFMAdapter:
        """Prepare the input batch and the conditioning tensors.

        Args:
            data: The dataset.
            treatment: Treatment column.
            outcome: Outcome column.
            adjustment_set: Recorded on the estimate for comparability with the
                classical arm; not consumed by the model.
            knowledge: Knowledge to condition on. ``None`` runs unconditioned.
            **kwargs: Reserved.

        Raises:
            ValueError: If the dataset exceeds the checkpoint's limits, or if
                the configured conditioning mode is not in the checkpoint's
                ``supports_conditioning``.
        """
        raise NotImplementedError

    def estimate(self) -> EffectEstimate:
        """Run inference and return the ATE, with revision and mode in diagnostics.

        Raises:
            RuntimeError: If called before :meth:`fit`.
        """
        raise NotImplementedError

    def estimate_cate(self, covariates: pd.DataFrame) -> np.ndarray:
        """Conditional effects at the given covariate values.

        Raises:
            NotImplementedError: If the loaded checkpoint is ATE-only.
        """
        raise NotImplementedError
