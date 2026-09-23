"""Audit causal foundation models under consistent-but-false knowledge.

Sweeps radius, conditioning mode and bias scale across the registered
checkpoints, with a classical OLS-on-``O*`` arm on identical data for reference.

Two gates run before the sweep and abort it if they fail. First, a zero-scale
bias must leave the forward pass unchanged; otherwise the injection is doing
something other than injecting knowledge. Second, true knowledge must beat no
knowledge; otherwise the model is not reading the conditioning input, and a flat
bias-versus-delta curve would look like robustness while meaning nothing.

Requires the ``[cfm]`` extra and downloaded checkpoints. See
``docs/CHECKPOINTS.md``.

Usage::

    python experiments/run_cfm_audit.py experiment=cfm_audit
    python experiments/run_cfm_audit.py experiment=cfm_audit compute.device=cuda
"""

from __future__ import annotations

import hydra
from omegaconf import DictConfig


@hydra.main(version_base=None, config_path="../configs", config_name="config")
def main(cfg: DictConfig) -> None:
    """Run the foundation-model audit.

    Steps once implemented:

    1. Seed and resolve the results directory.
    2. Resolve each checkpoint from the registry, refusing any without a pinned
       revision, and record the revisions in the manifest.
    3. Generate or load the audit datasets, keeping data fixed across the sweep
       so that only the knowledge varies.
    4. Run both sanity gates per checkpoint and abort on failure, recording what
       failed rather than skipping quietly.
    5. Sweep radius, conditioning mode and bias scale over all three arms.
    6. Run the classical reference arm on identical data at each radius.
    7. Write the audit table, including the gate results, whatever they say.

    Args:
        cfg: The composed Hydra config.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    raise NotImplementedError


if __name__ == "__main__":
    main()
