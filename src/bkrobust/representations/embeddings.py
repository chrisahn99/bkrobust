"""Order, box and hyperbolic-cone embeddings of the constraint poset.

Three families, each encoding partial order as a geometric relation:

``OrderEmbedding``
    Coordinate-wise dominance in the positive orthant. Simplest, and reflexive
    and transitive by construction -- the order axioms come free.
``BoxEmbedding``
    Axis-aligned hyperrectangles ordered by containment. Boxes have volumes, so
    "how much does this constraint pin down" becomes a volume ratio, and that
    is the most direct candidate for a geometric proxy of the breakdown radius.
``HyperbolicConeEmbedding``
    Entailment cones in hyperbolic space. Best suited to deep, tree-like orders;
    whether the constraint poset is deep enough to benefit is an empirical
    question, and on shallow posets it will not beat the Euclidean options.

All three need torch, from the ``[cfm]`` extra, and import it lazily so that
``import bkrobust`` stays light. The relevant baseline is not another embedding
-- it is the scalar structural features in
``configs/experiment/knowledge_embedding.yaml``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    import numpy as np

    from bkrobust.representations.poset import ConstraintPoset


class OrderEmbedding:
    """Order embedding: entailment as coordinate-wise dominance in the positive orthant.

    Args:
        dim: Embedding dimension.
        margin: Margin in the order-violation loss.
        negative_sampling: Non-entailed pairs drawn per entailed pair.
        lr: Learning rate.
        epochs: Maximum training epochs.
        device: Torch device.
        seed: Seed for initialisation and negative sampling.

    Raises:
        ImportError: If the ``[cfm]`` extra is not installed.
    """

    def __init__(
        self,
        dim: int = 64,
        margin: float = 1.0,
        negative_sampling: int = 5,
        lr: float = 1e-3,
        epochs: int = 100,
        device: str = "cpu",
        seed: int | None = None,
    ) -> None:
        raise NotImplementedError

    def fit(
        self,
        poset: ConstraintPoset,
        val_poset: ConstraintPoset | None = None,
    ) -> OrderEmbedding:
        """Train on the cover relation, early-stopping on validation order violation.

        Args:
            poset: Training poset.
            val_poset: Validation poset; ``None`` disables early stopping.

        Returns:
            ``self``, fitted.
        """
        raise NotImplementedError

    def transform(self, poset: ConstraintPoset) -> np.ndarray:
        """Embed the elements of ``poset``.

        Returns:
            An ``(n_elements, dim)`` array.

        Raises:
            RuntimeError: If called before :meth:`fit`.
        """
        raise NotImplementedError

    def score(self, poset: ConstraintPoset) -> dict[str, float]:
        """Reconstruction quality on a held-out poset.

        Returns:
            Keys ``"order_violation"``, ``"pairwise_accuracy"``, ``"auc"``.
        """
        raise NotImplementedError


class BoxEmbedding:
    """Box embedding: entailment as hyperrectangle containment, with calibrated volumes.

    Args:
        dim: Number of box dimensions.
        volume: ``"hard"``, ``"soft"`` or ``"bessel"``. Soft volumes keep a
            gradient when boxes are disjoint, which hard volumes do not, and
            disjoint boxes are the common case early in training.
        temperature: Softplus temperature for soft volumes.
        min_side: Numerical floor on side length.
        lr: Learning rate.
        epochs: Maximum training epochs.
        device: Torch device.
        seed: Seed.

    Raises:
        ImportError: If the ``[cfm]`` extra is not installed.
    """

    def __init__(
        self,
        dim: int = 32,
        volume: str = "soft",
        temperature: float = 1.0,
        min_side: float = 1e-4,
        lr: float = 1e-3,
        epochs: int = 100,
        device: str = "cpu",
        seed: int | None = None,
    ) -> None:
        raise NotImplementedError

    def fit(self, poset: ConstraintPoset, val_poset: ConstraintPoset | None = None) -> BoxEmbedding:
        """Train by maximising the conditional probabilities of covered pairs."""
        raise NotImplementedError

    def transform(self, poset: ConstraintPoset) -> np.ndarray:
        """Embed the elements as boxes.

        Returns:
            An ``(n_elements, 2 * dim)`` array of ``(min, max)`` corners.

        Raises:
            RuntimeError: If called before :meth:`fit`.
        """
        raise NotImplementedError

    def volume(self, poset: ConstraintPoset) -> np.ndarray:
        """Per-element box volume.

        The candidate geometric proxy for informativeness: a constraint that
        pins down more should occupy less volume.
        """
        raise NotImplementedError

    def score(self, poset: ConstraintPoset) -> dict[str, float]:
        """Reconstruction quality on a held-out poset."""
        raise NotImplementedError


class HyperbolicConeEmbedding:
    """Entailment cones in the Poincare ball.

    Args:
        dim: Embedding dimension.
        curvature: Negative curvature magnitude.
        aperture: Cone aperture parameter.
        lr: Riemannian learning rate.
        epochs: Maximum training epochs.
        device: Torch device.
        seed: Seed.

    Raises:
        ImportError: If the ``[cfm]`` extra is not installed.
    """

    def __init__(
        self,
        dim: int = 32,
        curvature: float = 1.0,
        aperture: float = 0.1,
        lr: float = 1e-3,
        epochs: int = 100,
        device: str = "cpu",
        seed: int | None = None,
    ) -> None:
        raise NotImplementedError

    def fit(
        self,
        poset: ConstraintPoset,
        val_poset: ConstraintPoset | None = None,
    ) -> HyperbolicConeEmbedding:
        """Train with Riemannian SGD on the cone-violation energy."""
        raise NotImplementedError

    def transform(self, poset: ConstraintPoset) -> np.ndarray:
        """Embed the elements into the Poincare ball.

        Raises:
            RuntimeError: If called before :meth:`fit`.
        """
        raise NotImplementedError

    def score(self, poset: ConstraintPoset) -> dict[str, float]:
        """Reconstruction quality on a held-out poset."""
        raise NotImplementedError


def build_embedding(cfg: Any) -> Any:
    """Instantiate the embedding named by a ``model`` config group.

    Args:
        cfg: See ``configs/model/``.

    Returns:
        An unfitted embedding.

    Raises:
        ImportError: If the ``[cfm]`` extra is not installed.
        KeyError: If ``cfg.target`` names no known embedding.
    """
    raise NotImplementedError
