"""Encode background knowledge in the form a foundation model consumes.

Three encodings, in increasing distance from the model's native inputs:

``ancestral_matrix``
    An ``{-1, 0, 1}`` matrix over variable pairs: ``1`` for asserted ancestor,
    ``-1`` for asserted non-ancestor, ``0`` for unstated. Passed as an extra
    input where the architecture has a slot for one.
``attention_bias``
    An additive bias on attention logits between variable tokens, following the
    approach of *Use What You Know* -- positive for asserted relations, negative
    for forbidden ones. Works with any transformer-based checkpoint without
    retraining, which is why it is the default.
``prompt``
    Knowledge serialised into whatever textual or tabular prefix the model
    accepts. Most fragile, and included only for checkpoints with no better
    route.

**The scale-zero invariant.** At ``bias_scale = 0.0`` the attention-bias
encoding must reproduce ``mode="none"`` bit for bit. If it does not, the
injection is changing the forward pass through some route other than the
intended bias -- a shape bug, a masking bug, a normalisation that sees the
added term -- and every number downstream is measuring that bug instead of the
knowledge. :func:`verify_zero_scale_identity` checks it and the audit calls it
before doing anything else.

**Interpretation, and its limits.** Injecting knowledge this way is not the same
operation as imposing it on a CPDAG. Meek's closure has a semantics: the
orientations it forces are entailed. An attention bias has none -- it is a
nudge, and the model may follow it, ignore it, or overreact to it. The audit
should therefore be read as measuring what these models *do* with a conditioning
signal, not what they are guaranteed to do with knowledge, and the paper should
say so plainly.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    import numpy as np

    from bkrobust.graphs.mpdag import Node
    from bkrobust.knowledge.base import BackgroundKnowledge


def knowledge_to_ancestral_matrix(
    bk: BackgroundKnowledge,
    nodes: Sequence[Node],
) -> np.ndarray:
    """Encode knowledge as an ``{-1, 0, 1}`` ancestral matrix.

    Entry ``M[i, j]`` is ``1`` if ``bk`` asserts that ``nodes[i]`` is an
    ancestor of ``nodes[j]``, ``-1`` if it asserts the negation, and ``0`` if it
    says nothing. Required edges and tier constraints are expanded into their
    implied ancestral assertions; the diagonal is ``0``.

    Args:
        bk: The knowledge to encode.
        nodes: Variable order, which fixes the row/column order. Must match the
            column order of the data the model sees, or the encoding silently
            describes a different graph.

    Returns:
        An ``(n_nodes, n_nodes)`` array of dtype ``int8``.

    Raises:
        ValueError: If ``bk`` names a node absent from ``nodes``.
    """
    raise NotImplementedError


def knowledge_to_adjacency_bias(
    bk: BackgroundKnowledge,
    nodes: Sequence[Node],
    *,
    bias_scale: float = 1.0,
    encode_forbidden: bool = True,
    symmetrize: bool = False,
) -> np.ndarray:
    """Encode knowledge as an additive attention-bias matrix.

    Entry ``B[i, j]`` is added to the attention logit from the token for
    ``nodes[i]`` to the token for ``nodes[j]``: ``+bias_scale`` for asserted
    edges, ``-bias_scale`` for forbidden ones when ``encode_forbidden``, ``0``
    otherwise.

    Args:
        bk: The knowledge to encode.
        nodes: Variable order.
        bias_scale: Magnitude of the bias. ``0.0`` must produce an all-zero
            matrix -- see :func:`verify_zero_scale_identity`.
        encode_forbidden: Emit negative bias for forbidden relations.
        symmetrize: Apply the bias in both directions. Loses the orientation
            information, so it is only appropriate for checkpoints whose
            attention over variable tokens is symmetric by construction.

    Returns:
        An ``(n_nodes, n_nodes)`` array of dtype ``float32``.

    Raises:
        ValueError: If ``bk`` names a node absent from ``nodes``, or
            ``bias_scale`` is negative.
    """
    raise NotImplementedError


def knowledge_to_prompt(
    bk: BackgroundKnowledge,
    nodes: Sequence[Node],
    *,
    template: str | None = None,
) -> str:
    """Serialise knowledge into a textual prefix.

    Fallback for checkpoints exposing no structured conditioning slot. Whatever
    template is used must be recorded in the run manifest: phrasing changes
    results here, and an unrecorded prompt makes the number unreproducible.
    """
    raise NotImplementedError


def encode(
    bk: BackgroundKnowledge,
    nodes: Sequence[Node],
    mode: str,
    **kwargs: Any,
) -> Any:
    """Dispatch to the encoding named by ``mode``.

    Args:
        bk: The knowledge to encode.
        nodes: Variable order.
        mode: ``"none"``, ``"attention_bias"``, ``"ancestral_matrix"`` or
            ``"prompt"``.
        **kwargs: Passed through to the chosen encoder.

    Returns:
        ``None`` for ``mode="none"``, otherwise the encoded object.

    Raises:
        KeyError: If ``mode`` is unknown.
    """
    raise NotImplementedError


def verify_zero_scale_identity(
    adapter: Any,
    data: Any,
    bk: BackgroundKnowledge,
    nodes: Sequence[Node],
    *,
    tol: float = 1e-6,
) -> bool:
    """Check that ``bias_scale=0.0`` reproduces unconditioned inference exactly.

    Runs the adapter twice -- once with ``mode="none"`` and once with
    ``mode="attention_bias"`` at zero scale -- and compares the outputs.

    Args:
        adapter: The checkpoint adapter to test.
        data: A batch to run through it.
        bk: Any non-empty knowledge set; at zero scale its content must not
            matter, which is precisely what is being tested.
        nodes: Variable order.
        tol: Absolute tolerance on the comparison.

    Returns:
        ``True`` if the two agree within ``tol``.

    Note:
        Called by :mod:`bkrobust.cfm.audit` before any sweep. A ``False`` here
        invalidates every downstream number for that checkpoint, so the audit
        must abort rather than warn.
    """
    raise NotImplementedError
