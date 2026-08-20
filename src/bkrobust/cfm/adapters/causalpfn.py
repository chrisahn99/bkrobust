"""CausalPFN adapter.

An amortized in-context estimator: the dataset is fed as context and the effect
comes out of a forward pass. There is no adjustment set anywhere in the
computation, which is exactly why the audit is interesting -- the classical
breakdown radii are defined through ``O*``, and this model has none, so any
breakdown it exhibits has to arise some other way.

Requires the ``[cfm]`` extra. See ``docs/CHECKPOINTS.md`` for the source,
revision and licence.
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


class CausalPFNAdapter:
    """Wrap a CausalPFN checkpoint as an :class:`~bkrobust.estimation.base.Estimator`.

    Args:
        checkpoint: Registry key or local path; ``None`` uses the registry
            default.
        revision: Checkpoint revision. Required for any reported number.
        device: ``"cpu"``, ``"cuda"`` or ``"mps"``.
        dtype: Inference dtype.
        batch_size: Inference batch size.
        conditioning: The ``conditioning`` config block -- mode, bias scale and
            related flags.

    Raises:
        ImportError: If the ``[cfm]`` extra is not installed.
    """

    name: str = "causalpfn"
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
    ) -> CausalPFNAdapter:
        """Load the context batch and the conditioning tensors.

        Args:
            data: The dataset, passed in context.
            treatment: Treatment column.
            outcome: Outcome column.
            adjustment_set: Ignored by the model itself, but recorded on the
                estimate so the amortized and classical arms stay comparable
                row for row.
            knowledge: Knowledge to condition on, encoded per the adapter's
                ``conditioning`` config. ``None`` runs unconditioned.
            **kwargs: Reserved.

        Raises:
            ValueError: If the dataset exceeds the checkpoint's feature or
                sample limits (see :class:`~bkrobust.cfm.registry.CheckpointSpec`).
        """
        raise NotImplementedError

    def estimate(self) -> EffectEstimate:
        """Run the forward pass and return the ATE.

        The estimate carries the checkpoint revision, the conditioning mode and
        the bias scale in its diagnostics.

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
