"""The BackgroundKnowledge container.

Four kinds of constraint, in rough order of how often practitioners supply them:

``tiers``
    An ordered partition of the variables -- temporal order, measurement waves,
    an assumed causal hierarchy. Forbids every edge pointing from a later tier
    to an earlier one. Cheap to state and expensive to get wrong: one misplaced
    variable forbids a whole block of edges at once.
``forbidden_edges``
    "``a`` cannot cause ``b``." Usually from domain impossibility.
``required_edges``
    "``a`` causes ``b``." The strongest claim, and the one experts are most
    confident about and least often right about.
``ancestral``
    "``a`` is (not) an ancestor of ``b``." Weaker than an edge claim -- it
    permits mediation -- but it can still force orientations far from the pair
    named, because the closure has to make room for a directed path.

The container is immutable. Perturbations return new instances rather than
mutating, so a sampled knowledge set can be cached, hashed, and compared
against the one it was derived from.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, TypeAlias

if TYPE_CHECKING:  # pragma: no cover
    from bkrobust.graphs.mpdag import MPDAG, Edge, Node

#: An ancestral constraint: ``(a, b, holds)`` reads "a is an ancestor of b" when
#: ``holds`` is ``True``, and "a is not an ancestor of b" when ``False``.
AncestralConstraint: TypeAlias = "tuple[Node, Node, bool]"


@dataclass(frozen=True)
class BackgroundKnowledge:
    """An immutable body of background knowledge over a fixed variable set.

    Attributes:
        required_edges: Pairs ``(a, b)`` asserting the directed edge ``a -> b``.
        forbidden_edges: Pairs ``(a, b)`` asserting that ``a -> b`` is absent.
            Forbidding both directions asserts non-adjacency.
        tiers: An ordered partition; every edge from a later tier to an earlier
            one is forbidden. Nodes absent from every tier are unconstrained.
        ancestral: Ancestral constraints as ``(a, b, holds)`` triples.
        source: Provenance tag -- ``"manual"``, ``"sampled"``, ``"kg"``,
            ``"llm"``. Carried into the run manifest so a result can always be
            traced back to where its knowledge came from.
        metadata: Free-form annotations. The samplers record the perturbation
            kind, the realised radius and the rejection count here.
    """

    required_edges: frozenset[Edge] = frozenset()
    forbidden_edges: frozenset[Edge] = frozenset()
    tiers: tuple[frozenset[Node], ...] = ()
    ancestral: frozenset[AncestralConstraint] = frozenset()
    source: str = "manual"
    metadata: dict[str, object] = field(default_factory=dict, compare=False)

    # --- construction ------------------------------------------------------

    @classmethod
    def empty(cls) -> BackgroundKnowledge:
        """Return the knowledge set that asserts nothing.

        The no-knowledge control arm; imposing it leaves a CPDAG unchanged.
        """
        raise NotImplementedError

    @classmethod
    def from_dag(
        cls,
        dag: MPDAG,
        kinds: Sequence[str] = ("required_edges",),
        subset: Iterable[Edge] | None = None,
    ) -> BackgroundKnowledge:
        """Extract true knowledge from a known DAG.

        Builds the consistent-and-*true* control arm: knowledge of the same
        shape and size as a sampled perturbation, but correct. Without this arm
        an experiment cannot separate the effect of knowledge being wrong from
        the effect of knowledge being present at all.

        Args:
            dag: The ground-truth DAG.
            kinds: Which constraint kinds to extract.
            subset: Restrict extraction to these edges; ``None`` uses all.

        Returns:
            Knowledge that is true of ``dag``.
        """
        raise NotImplementedError

    @classmethod
    def from_config(cls, cfg: object) -> BackgroundKnowledge:
        """Build from a Hydra/OmegaConf node.

        Reads the ``knowledge`` block of an experiment config -- see
        ``configs/experiment/real_data.yaml`` for the expected shape.

        Raises:
            ValueError: If the config names an unknown constraint kind or
                malformed node pairs.
        """
        raise NotImplementedError

    # --- derivation --------------------------------------------------------

    def with_edge(self, tail: Node, head: Node, required: bool = True) -> BackgroundKnowledge:
        """Return a copy with one more edge constraint.

        Args:
            tail: Source endpoint.
            head: Target endpoint.
            required: Add to ``required_edges`` if ``True``, else to
                ``forbidden_edges``.

        Returns:
            A new instance; ``self`` is unchanged.
        """
        raise NotImplementedError

    def without_edge(self, tail: Node, head: Node) -> BackgroundKnowledge:
        """Return a copy with any constraint on ``(tail, head)`` removed."""
        raise NotImplementedError

    def merge(self, other: BackgroundKnowledge) -> BackgroundKnowledge:
        """Return the union of two knowledge sets.

        Raises:
            ValueError: If the two directly contradict each other -- the same
                edge required by one and forbidden by the other -- since the
                union would then be unsatisfiable by any DAG, not merely
                inconsistent with a particular CPDAG.
        """
        raise NotImplementedError

    def restricted_to(self, nodes: Iterable[Node]) -> BackgroundKnowledge:
        """Return the constraints whose endpoints all lie in ``nodes``."""
        raise NotImplementedError

    # --- expansion ---------------------------------------------------------

    def implied_forbidden_edges(self, nodes: Iterable[Node]) -> frozenset[Edge]:
        """Expand :attr:`tiers` into the explicit set of forbidden edges.

        Tiers are stored compactly but act as many constraints; the closure and
        the distance metrics need them expanded.

        Args:
            nodes: The variable set the tiers are interpreted over.

        Returns:
            Every ``(a, b)`` forbidden because ``a`` sits in a later tier than
            ``b``, together with :attr:`forbidden_edges`.
        """
        raise NotImplementedError

    def constraints(self) -> Iterator[tuple[str, tuple[Node, Node]]]:
        """Iterate every constraint as ``(kind, (a, b))``, in a deterministic order.

        Order is fixed so that Meek's Algorithm 1 -- which imposes constraints
        one at a time and can FAIL partway -- gives a reproducible verdict.

        Yields:
            ``kind`` is one of ``"required"``, ``"forbidden"``, ``"tier"``,
            ``"ancestor"``, ``"non_ancestor"``.
        """
        raise NotImplementedError

    def nodes_mentioned(self) -> frozenset[Node]:
        """Every node named by at least one constraint."""
        raise NotImplementedError

    # --- summary -----------------------------------------------------------

    def __len__(self) -> int:
        """Total number of constraints, with tiers expanded to forbidden edges."""
        raise NotImplementedError

    def __bool__(self) -> bool:
        """Whether this knowledge set asserts anything at all."""
        raise NotImplementedError

    def summary(self) -> dict[str, int]:
        """Counts per constraint kind, for logging and the run manifest."""
        raise NotImplementedError

    def __str__(self) -> str:
        """Render the constraints in a form an analyst can read and check."""
        raise NotImplementedError
