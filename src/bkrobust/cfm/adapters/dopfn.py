"""Do-PFN adapter.

Do-PFN answers interventional queries directly -- it is trained to predict the
post-intervention distribution rather than an average effect -- so the ATE is
obtained by averaging its conditional predictions rather than read off a single
head. That makes it the natural place to look at PEHE alongside ATE bias, and
the only checkpoint in the audit where wrong knowledge could plausibly distort
the *shape* of the effect surface while leaving its average intact.

Because the estimand is interventional by construction, this adapter is also the
sharpest test of whether attention-bias conditioning behaves like an
intervention on the model's implicit graph or merely like a hint. Do not assume
the former; the audit is designed to find out.

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


class DoPFNAdapter:
    """Wrap a Do-PFN checkpoint as an :class:`~bkrobust.estimation.base.Estimator`.

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

    name: str = "dopfn"
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
    ) -> DoPFNAdapter:
        """Prepare the in-context batch and the interventional query.

        Args:
            data: The dataset.
            treatment: The intervened-upon column.
            outcome: Outcome column.
            adjustment_set: Recorded for comparability; not consumed.
            knowledge: Knowledge to condition on. ``None`` runs unconditioned.
            **kwargs: Reserved.

        Raises:
            ValueError: If the dataset exceeds the checkpoint's limits.
        """
        raise NotImplementedError

    def estimate(self) -> EffectEstimate:
        """Return the ATE, obtained by averaging the predicted conditional effects.

        The averaging weights are recorded in the diagnostics, since an ATE
        computed under a different covariate distribution is a different
        estimand and the audit tables must not conflate them.

        Raises:
            RuntimeError: If called before :meth:`fit`.
        """
        raise NotImplementedError

    def estimate_cate(self, covariates: pd.DataFrame) -> np.ndarray:
        """Conditional effects at the given covariate values -- the native estimand.

        Raises:
            RuntimeError: If called before :meth:`fit`.
        """
        raise NotImplementedError
