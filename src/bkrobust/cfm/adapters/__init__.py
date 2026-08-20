"""Per-checkpoint adapters conforming to the Estimator protocol.

Each adapter hides one checkpoint's input conventions -- feature ordering,
scaling, in-context format, where the conditioning tensor goes -- behind
``fit`` / ``estimate``. Uniformity is what makes the audit a comparison rather
than three anecdotes.

Every adapter must:

* import torch lazily, inside ``__init__``, so ``import bkrobust`` works without
  the ``[cfm]`` extra;
* refuse to load an unpinned checkpoint;
* honour ``conditioning.mode="none"`` by passing no bias at all, rather than a
  zero tensor -- the two must agree numerically, and checking that they do is
  the point of
  :func:`~bkrobust.cfm.conditioning.verify_zero_scale_identity`, which is not a
  check the adapter can perform on itself;
* record the checkpoint revision on every returned
  :class:`~bkrobust.estimation.base.EffectEstimate`.
"""

from __future__ import annotations

__all__: list[str] = []
