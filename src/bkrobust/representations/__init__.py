"""Learned geometry of the constraint poset -- CONDITIONAL component, may be cut.

The hypothesis: background-knowledge constraints form a partial order under
entailment, that order has geometry, and the geometry carries information about
consequences -- so an embedding of a knowledge set might predict its breakdown
radius, or separate consistent-and-true from consistent-but-false knowledge,
without running the closure at all.

The reasons to be sceptical, stated up front so the section is not written
backwards from a hoped-for result:

* The radius is computable exactly on the graphs where it matters most, so a
  learned predictor of it is a speed optimisation, not a new capability, and
  must be justified on those terms.
* Separating true from false knowledge using only the CPDAG and the constraints
  would be a strong claim -- consistent-but-false knowledge is by construction
  indistinguishable from true knowledge *given the CPDAG*. Any separation a
  probe achieves is therefore picking up a structural regularity of how the
  perturbations were *sampled*, not truth. If this probe succeeds, the first
  hypothesis to test is that it has learned the sampler.

Both probes must beat the structural baselines in
``configs/experiment/knowledge_embedding.yaml`` -- constraint count, cascade
size, SHD -- and if they do not, this component does not go in the paper.
"""

from __future__ import annotations

__all__: list[str] = []
