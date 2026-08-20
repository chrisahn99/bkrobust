"""Adapters for background knowledge that comes from a real source.

The sampled knowledge in :mod:`bkrobust.knowledge.perturb` is the instrument:
controlled radius, controlled kind, known ground truth. It answers "what would
knowledge this wrong do". It does not answer "how wrong is the knowledge people
actually use", and that second question is where the paper's claim about
practice has to be earned.

Two sources that practitioners genuinely draw on:

* **Knowledge graphs.** Curated relations mined into required or forbidden
  edges. Precise, incomplete, and prone to conflating association with
  causation -- a co-occurrence edge read as a causal one is exactly a spurious
  constraint.
* **LLM elicitation.** Increasingly common and largely unaudited. An LLM will
  answer every orientation question asked of it, at high confidence, with no
  signal distinguishing recalled fact from plausible-sounding guess. Its output
  is consistent-but-false knowledge in the typical case, which makes it the
  most direct real-world instance of the failure mode this project studies.

Both adapters are stubs and both are optional to the main results. If they land,
they supply the empirical anchor: a measured distribution over ``delta`` for
knowledge from each source, placed against the breakdown radii computed for the
same graphs. Neither may be allowed to smuggle a network call into a test.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Protocol

if TYPE_CHECKING:  # pragma: no cover
    from bkrobust.graphs.mpdag import MPDAG, Node
    from bkrobust.knowledge.base import BackgroundKnowledge


@dataclass(frozen=True)
class ElicitedConstraint:
    """One constraint from an external source, with whatever confidence it carried.

    Attributes:
        pair: The node pair the constraint is about.
        kind: ``"required"``, ``"forbidden"``, ``"ancestor"``, ``"non_ancestor"``.
        confidence: Source-reported confidence in ``[0, 1]``, or ``None`` when
            the source reports none. Retained rather than thresholded away: a
            confidence-weighted knowledge set is a different object from a hard
            one, and whether calibration helps is an open question this scaffold
            should not foreclose.
        provenance: Free-form record of where it came from -- KG relation id,
            model name and prompt hash. Required for the run manifest.
    """

    pair: tuple[Node, Node]
    kind: str
    confidence: float | None
    provenance: dict[str, Any]


class KnowledgeElicitor(Protocol):
    """Anything that can produce constraints over a named variable set."""

    def elicit(
        self,
        nodes: list[Node],
        context: dict[str, Any] | None = None,
    ) -> list[ElicitedConstraint]:
        """Return constraints over ``nodes``.

        Args:
            nodes: The variables to ask about.
            context: Domain description, variable glossary, dataset card --
                whatever the source needs to interpret the variable names.

        Returns:
            The elicited constraints; may be empty.
        """
        ...


class KnowledgeGraphElicitor:
    """Derive constraints from a curated knowledge graph.

    Args:
        source: Identifier of the knowledge graph.
        relation_map: Mapping from KG relation types to constraint kinds. The
            interesting decisions live here: which relations are read as causal
            at all, and which are demoted to non-adjacency.
        min_confidence: Drop constraints below this confidence.
    """

    def __init__(
        self,
        source: str,
        relation_map: dict[str, str] | None = None,
        min_confidence: float = 0.0,
    ) -> None:
        raise NotImplementedError

    def elicit(
        self,
        nodes: list[Node],
        context: dict[str, Any] | None = None,
    ) -> list[ElicitedConstraint]:
        """Look up ``nodes`` in the knowledge graph and map relations to constraints."""
        raise NotImplementedError


class LLMElicitor:
    """Elicit constraints from a language model, one pairwise question at a time.

    Args:
        model: Model identifier.
        prompt_template: Template for the pairwise orientation question.
        temperature: Sampling temperature; ``0.0`` for reproducibility.
        n_votes: Ask each pair this many times and take the majority. A crude
            confidence estimate, and worth reporting as such rather than
            presenting as calibrated.
        cache_path: Where to cache responses. Required in practice -- results
            must be reproducible without re-querying a model that may have
            changed underneath the experiment.
    """

    def __init__(
        self,
        model: str,
        prompt_template: str | None = None,
        temperature: float = 0.0,
        n_votes: int = 1,
        cache_path: str | None = None,
    ) -> None:
        raise NotImplementedError

    def elicit(
        self,
        nodes: list[Node],
        context: dict[str, Any] | None = None,
    ) -> list[ElicitedConstraint]:
        """Query the model over all node pairs and parse the answers into constraints.

        Raises:
            RuntimeError: If no cache is configured and no backend is available.
                Tests must never reach a network call.
        """
        raise NotImplementedError


def to_background_knowledge(
    constraints: list[ElicitedConstraint],
    *,
    min_confidence: float = 0.0,
    drop_inconsistent: bool = True,
    cpdag: MPDAG | None = None,
) -> BackgroundKnowledge:
    """Assemble elicited constraints into a usable knowledge set.

    Args:
        constraints: The elicited constraints.
        min_confidence: Drop anything below this.
        drop_inconsistent: Drop constraints that make the set inconsistent with
            ``cpdag``, greedily in decreasing confidence order. This mirrors what
            a practitioner does when Meek's Algorithm 1 reports FAIL, and it is
            worth being explicit about what it accomplishes: it restores
            consistency, so the check passes, while doing nothing whatever about
            truth. The surviving set is the paper's object of study.
        cpdag: Required when ``drop_inconsistent`` is set.

    Returns:
        The assembled knowledge, with ``source`` set from the constraints'
        provenance and ``metadata`` recording how many were dropped and why.

    Raises:
        ValueError: If ``drop_inconsistent`` is set and ``cpdag`` is ``None``.
    """
    raise NotImplementedError
