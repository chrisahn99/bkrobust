"""bkrobust -- robustness of causal inference to consistent-but-false background knowledge.

Background knowledge injected into causal discovery (required and forbidden
edges, tiers, ancestral constraints) is checked only for *consistency*: Meek's
Algorithm 1 either terminates with an MPDAG or reports FAIL. It is never
checked for *truth*. Knowledge that is consistent but false therefore passes
every check in current practice, and this package quantifies what such
knowledge costs.

The package is organised by the three levels of the claim:

``bkrobust.graphs``
    MPDAG representation, Meek rules, consistency, valid adjustment,
    the optimal adjustment set ``O*``, and distances between graphs.
``bkrobust.knowledge``
    :class:`~bkrobust.knowledge.base.BackgroundKnowledge` as a first-class
    object, a taxonomy of misspecifications, and samplers that draw knowledge
    which is consistent with an observed CPDAG yet false with respect to a
    known DAG.
``bkrobust.theory``
    The two breakdown radii ``delta_opt`` and ``delta_valid``, the bias bound
    beyond ``delta_valid``, the efficiency gap between them, and the
    practitioner-facing certificate.
``bkrobust.estimation``
    Effect estimators evaluated under a supplied adjustment set.
``bkrobust.cfm``
    Audit of causal foundation models conditioned on the same knowledge.
``bkrobust.representations``
    Conditional component: learned geometry of the constraint poset.
``bkrobust.data``, ``bkrobust.metrics``, ``bkrobust.utils``
    Generators and loaders, evaluation metrics, and shared plumbing.

Nothing in this package is implemented yet: this is a scaffold, and every
function body raises :class:`NotImplementedError`. See ``docs/THEORY.md`` for
notation and the theorem targets, and ``docs/EXPERIMENTS.md`` for what each
experiment is meant to establish.
"""

from __future__ import annotations

__version__ = "0.1.0"

__all__ = [
    "__version__",
]
