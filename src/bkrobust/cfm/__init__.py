"""Audit of causal foundation models under consistent-but-false knowledge.

Amortized causal models -- CausalPFN, CausalFM, Do-PFN -- are trained on
millions of synthetic SCMs and estimate effects in a forward pass, with no
explicit graph, no explicit adjustment set, and no identification argument the
user can inspect. Increasingly they accept background knowledge as an input,
typically injected as a bias on the attention pattern.

That raises the question this module exists to answer. In the classical
pipeline, wrong knowledge does its damage through a mechanism we can name: it
changes the MPDAG, which changes ``O*``, which changes what is conditioned on.
An amortized model has no ``O*`` to change. So does wrong knowledge hurt it in
the same way, differently, or not at all? Three outcomes are possible and all
are informative -- the model may inherit the same breakdown radii, may be more
robust because it hedges across the equivalence class it was trained over, or
may be *less* robust because it over-trusts a conditioning signal that it was
trained to treat as reliable.

Module map:

``registry``
    Checkpoint sources, versions and licences.
``conditioning``
    Encode a :class:`~bkrobust.knowledge.base.BackgroundKnowledge` as the
    matrices these models consume.
``adapters``
    One wrapper per checkpoint, all conforming to
    :class:`~bkrobust.estimation.base.Estimator`.
``audit``
    The bias-versus-delta sweep.

Everything here needs the ``[cfm]`` extra. Imports are guarded so that
``import bkrobust`` works without torch, and the core theory and simulation code
must never grow a hard dependency on this subpackage.
"""

from __future__ import annotations

__all__: list[str] = []
