"""Background knowledge as a first-class object.

In most causal-discovery tooling background knowledge is a bag of arguments
threaded into an algorithm and then discarded. Here it is a value with its own
type, its own distances, its own samplers and its own taxonomy of ways to be
wrong -- because the object of study is the knowledge, not the algorithm.

Module map:

``base``
    :class:`~bkrobust.knowledge.base.BackgroundKnowledge`: required edges,
    forbidden edges, tiers, ancestral constraints.
``taxonomy``
    The kinds of misspecification -- orientation, ancestral, tier, and the
    missing/spurious axis -- as an enumerable, so experiments sweep over it
    rather than hardcoding one kind.
``perturb``
    Samplers that draw knowledge which is consistent with an observed CPDAG and
    false with respect to a known DAG, at a controlled radius.
``cascade``
    Meek-closure amplification: how many orientations ``k`` imposed constraints
    actually force.
``elicit``
    Adapters for knowledge that comes from somewhere real -- a knowledge graph,
    an LLM -- rather than from a sampler.
"""

from __future__ import annotations

__all__: list[str] = []
