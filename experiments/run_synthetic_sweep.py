"""Delta sweep on synthetic MPDAGs -- the structural-level experiment.

Establishes the two breakdown radii empirically and tests the conjectured
ordering ``delta_opt <= delta_valid``. For each sampled graph and SCM it walks
the radius grid, draws consistent-but-false knowledge at each radius, and
records what happens to ``O*``, to the asymptotic variance, and to the estimate.

Runs three arms at every radius -- no knowledge, size-matched true knowledge,
and false knowledge -- because the interesting quantity is the difference
between the last two, not the level of any one.

Usage::

    python experiments/run_synthetic_sweep.py
    python experiments/run_synthetic_sweep.py graph=scale_free perturbation=tier
    python experiments/run_synthetic_sweep.py -m seed=1,2,3,4,5
"""

from __future__ import annotations

import hydra
from omegaconf import DictConfig


@hydra.main(version_base=None, config_path="../configs", config_name="config")
def main(cfg: DictConfig) -> None:
    """Run the synthetic delta sweep.

    Steps once implemented:

    1. Seed from ``cfg.seed`` and spawn one generator per replicate.
    2. Resolve the results directory and write the manifest.
    3. For each graph replicate: sample a DAG, attach an SCM, sample data,
       compute the CPDAG, and select a target pair.
    4. Record the leakage diagnostics -- a highly varsortable draw makes the
       rest of the numbers uninterpretable, so this comes before the sweep and
       not after.
    5. Compute ``delta_valid`` and ``delta_opt`` exactly where the graph is small
       enough, and by heuristic search otherwise, recording which path was used.
    6. For each radius in the grid and each arm: draw knowledge, impose it, read
       off ``O*``, estimate the effect, and record the structural and effect
       metrics together.
    7. Check the ordering and persist the full witness for any exact violation.
    8. Write results and aggregates.

    Args:
        cfg: The composed Hydra config.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    raise NotImplementedError


if __name__ == "__main__":
    main()
