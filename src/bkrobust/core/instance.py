"""The instance schema: one serialisable record per experimental unit. Frozen.

An :class:`Instance` is everything needed to reproduce and interpret one
measurement: the graphs, the simulated knowledge, the adjustment set, the seeds,
the computed radii, the method that produced them, timings, and the graph
descriptors the analysis regresses on.

The governing cost parameter is **not** ``n``. It is the size and density of the
largest chordal component that the knowledge intersects, since that component is
what the space is enumerated over. :func:`describe_graphs` computes it and every
scaling plot uses it on the abscissa.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

from bkrobust.demo.graph import MPDAG, undirected_components

Edge = tuple[str, str]


@dataclass(frozen=True)
class GraphDescriptors:
    """Structural summary of one instance's graphs.

    Attributes:
        n_nodes: Vertices.
        n_edges: Skeleton edges.
        n_compelled: Edges the data orients (frozen; never perturbed).
        n_undirected: Edges background knowledge may orient.
        n_components: Undirected (chordal) components.
        largest_component_nodes: Vertices in the largest component.
        largest_component_edges: Edges in the largest component.
        knowledge_component_nodes: Vertices in the largest component the
            knowledge intersects.
        knowledge_component_edges: Edges in that component. **The governing cost
            parameter.**
        knowledge_component_density: Its edge density,
            ``2E / (V(V-1))``, or 0 for a component with fewer than two nodes.
        treatment_in_component: Whether the treatment lies inside an undirected
            component -- the condition that makes perturbation able to reach it.
    """

    n_nodes: int
    n_edges: int
    n_compelled: int
    n_undirected: int
    n_components: int
    largest_component_nodes: int
    largest_component_edges: int
    knowledge_component_nodes: int
    knowledge_component_edges: int
    knowledge_component_density: float
    treatment_in_component: bool


def describe_graphs(cpdag: MPDAG, knowledge: list[Edge], treatment: str) -> GraphDescriptors:
    """Compute the structural descriptors for one instance.

    Args:
        cpdag: The estimated CPDAG.
        knowledge: The asserted orientations.
        treatment: The treatment node.

    Returns:
        A :class:`GraphDescriptors`.
    """
    comps = undirected_components(cpdag)
    und = cpdag.undirected_edges

    def comp_edges(c: frozenset[str]) -> int:
        return sum(1 for a, b in und if a in c and b in c)

    largest_n = max((len(c) for c in comps), default=0)
    largest_e = max((comp_edges(c) for c in comps), default=0)

    touched = [c for c in comps if any(a in c or b in c for a, b in knowledge)]
    if touched:
        best = max(touched, key=lambda c: (comp_edges(c), len(c)))
        k_n, k_e = len(best), comp_edges(best)
    else:
        k_n, k_e = 0, 0
    density = (2.0 * k_e / (k_n * (k_n - 1))) if k_n > 1 else 0.0

    return GraphDescriptors(
        n_nodes=len(cpdag.nodes),
        n_edges=len(cpdag.skeleton()),
        n_compelled=len(cpdag.directed_edges),
        n_undirected=len(und),
        n_components=len(comps),
        largest_component_nodes=largest_n,
        largest_component_edges=largest_e,
        knowledge_component_nodes=k_n,
        knowledge_component_edges=k_e,
        knowledge_component_density=density,
        treatment_in_component=any(treatment in c for c in comps),
    )


@dataclass(frozen=True)
class Instance:
    """One experimental unit, fully serialisable.

    Graphs are stored as ``edge_string()`` renderings: compact, human-readable,
    and losslessly re-parseable given the node set.
    """

    instance_id: str
    generator: str
    seed: int
    params: dict[str, Any]

    true_dag: str
    cpdag: str
    k_true: tuple[Edge, ...]
    k_assumed: tuple[Edge, ...]
    corruption: str
    corruption_rate: float
    k_assumed_consistent: bool

    g0: str
    treatment: str
    outcome: str
    z: tuple[str, ...]

    descriptors: GraphDescriptors

    space_size: int
    n_covers: int
    max_shell: int

    r_val: int
    r_opt: int
    r_eps: dict[str, int]
    method: str

    timings: dict[str, float] = field(default_factory=dict)
    notes: str = ""

    def to_row(self) -> dict[str, Any]:
        """Flatten to a single CSV row.

        Nested descriptors are prefixed ``desc_``; the epsilon radii ``r_eps_``.
        """
        row: dict[str, Any] = {
            k: v
            for k, v in asdict(self).items()
            if k not in ("descriptors", "r_eps", "params", "timings", "k_true", "k_assumed", "z")
        }
        row["z"] = "|".join(self.z)
        row["k_true"] = "|".join(f"{a}->{b}" for a, b in self.k_true)
        row["k_assumed"] = "|".join(f"{a}->{b}" for a, b in self.k_assumed)
        for k, v in asdict(self.descriptors).items():
            row[f"desc_{k}"] = v
        for k, v in self.r_eps.items():
            row[f"r_eps_{k}"] = v
        for k, v in self.params.items():
            row[f"param_{k}"] = v
        for k, v in self.timings.items():
            row[f"time_{k}"] = v
        return row
