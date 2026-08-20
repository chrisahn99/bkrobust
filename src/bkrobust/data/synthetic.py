"""Random DAGs, SCMs, and the CPDAGs the experiments actually see.

The generation pipeline, in the order the experiments use it:

1. Sample a DAG from a random graph family.
2. Attach an SCM -- linear-Gaussian by default, additive-noise otherwise.
3. Sample observational data from it.
4. Compute the CPDAG.

Step 4 has a choice with consequences. The CPDAG can come from the *true DAG*
(the exact equivalence class, obtained by graph operations) or from *running
discovery on the sampled data* (an estimate, with its own errors). The first
isolates the effect of background knowledge from the effect of discovery error,
and is the right default for the theory experiments -- it is the only way a
measured radius is attributable to the knowledge alone. The second is what a
practitioner faces, and the real-data experiment uses it.

Both are supported; ``cpdag_from`` selects. Which one produced a number must be
recorded, because a radius computed against an estimated CPDAG confounds two
sources of error and is not comparable to one computed against the exact class.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    import numpy as np
    import pandas as pd

    from bkrobust.graphs.mpdag import MPDAG, Node


@dataclass(frozen=True)
class LinearGaussianSCM:
    """A linear-Gaussian structural causal model on a known DAG.

    Attributes:
        graph: The DAG.
        weights: Edge coefficients keyed by ``(parent, child)``.
        noise_scales: Per-node noise standard deviations. Deliberately unequal
            by default: equal variances make the causal order recoverable from
            marginal variances alone, which is a leak and not a property of real
            data. See :mod:`bkrobust.data.leakage`.
        intercepts: Per-node intercepts.
    """

    graph: MPDAG
    weights: dict[tuple[Node, Node], float]
    noise_scales: dict[Node, float]
    intercepts: dict[Node, float]

    def sample(self, n: int, rng: np.random.Generator) -> pd.DataFrame:
        """Draw ``n`` observations in topological order.

        Returns:
            A DataFrame whose columns follow ``graph.nodes``.
        """
        raise NotImplementedError

    def covariance(self) -> np.ndarray:
        """Population covariance matrix, in ``graph.nodes`` order.

        Closed form from the weight matrix and noise scales; this is what
        :mod:`bkrobust.graphs.optimality` uses, so asymptotic variances are
        exact rather than estimated.
        """
        raise NotImplementedError

    def total_effect(self, treatment: Node, outcome: Node) -> float:
        """Ground-truth total effect: the sum over directed paths of path products."""
        raise NotImplementedError

    def intervene(self, node: Node, value: float) -> LinearGaussianSCM:
        """Return the SCM under ``do(node = value)``, with the node's parents cut."""
        raise NotImplementedError


def random_dag(
    n_nodes: int,
    rng: np.random.Generator,
    *,
    generator: str = "erdos_renyi",
    edge_prob: float = 0.3,
    m_attach: int = 2,
    max_in_degree: int | None = None,
) -> MPDAG:
    """Sample a random DAG.

    Args:
        n_nodes: Number of nodes.
        rng: Seeded generator.
        generator: ``"erdos_renyi"`` or ``"scale_free"``.
        edge_prob: Edge probability, for Erdos-Renyi.
        m_attach: Attachment count, for scale-free.
        max_in_degree: Optional cap on in-degree.

    Returns:
        A fully oriented MPDAG with nodes labelled ``"X1"``..``"Xn"``.

    Raises:
        ValueError: If ``generator`` is unknown or the parameters are out of range.
    """
    raise NotImplementedError


def random_scm(
    dag: MPDAG,
    rng: np.random.Generator,
    *,
    family: str = "linear_gaussian",
    weight_range: tuple[float, float] = (0.5, 2.0),
    weight_sign: str = "random",
    noise: str = "gaussian",
    noise_scale: tuple[float, float] = (0.5, 1.5),
) -> Any:
    """Attach a randomly parameterised SCM to ``dag``.

    Args:
        dag: The DAG.
        rng: Seeded generator.
        family: ``"linear_gaussian"``, ``"anm"`` or ``"post_nonlinear"``.
        weight_range: Range for ``|coefficient|``. Bounded away from zero:
            near-zero coefficients make an edge causally inert, so the graph
            would no longer describe the distribution and the ground truth would
            be wrong in a way that quietly flatters every method.
        weight_sign: ``"random"`` or ``"positive"``.
        noise: Noise family.
        noise_scale: Range for the per-node noise standard deviations.

    Returns:
        A :class:`LinearGaussianSCM`, or the corresponding SCM object for other
        families.

    Raises:
        ValueError: If ``family`` is unknown.
    """
    raise NotImplementedError


def dag_to_cpdag(dag: MPDAG) -> MPDAG:
    """Return the exact CPDAG of ``dag``: keep the v-structures, take the Meek closure.

    The no-discovery-error path.
    """
    raise NotImplementedError


def discover_cpdag(
    data: pd.DataFrame,
    *,
    algorithm: str = "pc",
    alpha: float = 0.01,
    **kwargs: Any,
) -> MPDAG:
    """Estimate a CPDAG from data.

    Args:
        data: Observational data.
        algorithm: ``"pc"`` or ``"ges"``.
        alpha: Significance level for conditional independence tests.
        **kwargs: Algorithm-specific options.

    Returns:
        The estimated CPDAG.

    Raises:
        ImportError: If no causal-discovery backend is installed. The backend is
            an unresolved decision point -- see ``pyproject.toml``.
    """
    raise NotImplementedError


def select_target_pair(
    dag: MPDAG,
    rng: np.random.Generator,
    *,
    selection: str = "random_pair",
    require_causal_path: bool = True,
) -> tuple[Node, Node]:
    """Choose a ``(treatment, outcome)`` pair to study on ``dag``.

    Args:
        dag: The DAG.
        rng: Seeded generator.
        selection: ``"random_pair"``, ``"max_path_count"`` or ``"fixed"``.
        require_causal_path: Reject pairs with no directed path, whose true
            effect is zero and whose bias is therefore uninformative about scale.

    Returns:
        The chosen pair.

    Raises:
        ValueError: If no pair satisfies the constraints.
    """
    raise NotImplementedError


def generate_dataset(
    cfg: Any,
    rng: np.random.Generator,
) -> dict[str, Any]:
    """Run the whole pipeline from a ``graph`` config group.

    Args:
        cfg: See ``configs/graph/``.
        rng: Seeded generator.

    Returns:
        Keys ``"dag"``, ``"cpdag"``, ``"scm"``, ``"data"``, ``"treatment"``,
        ``"outcome"``, ``"true_effect"``, and ``"leakage"`` carrying the
        varsortability diagnostics.
    """
    raise NotImplementedError
