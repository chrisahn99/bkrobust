"""Semi-synthetic sweep: real covariates, simulated outcomes, known truth.

Repeats the synthetic sweep on data with realistic marginals, to check that the
radii are not an artefact of the synthetic generator. The failure mode being
guarded against is specific: synthetic benchmarks are often highly varsortable,
and a robustness result that only holds on varsortable data is a result about
the benchmark.

Usage::

    python experiments/run_semisynthetic.py experiment=semisynthetic
    python experiments/run_semisynthetic.py experiment=semisynthetic experiment.generator=realcause
"""

from __future__ import annotations

import hydra
from omegaconf import DictConfig


@hydra.main(version_base=None, config_path="../configs", config_name="config")
def main(cfg: DictConfig) -> None:
    """Run the semi-synthetic sweep.

    Steps once implemented:

    1. Seed and resolve the results directory.
    2. Load the generator's datasets and their reference networks.
    3. Simulate treatment and outcome on the real covariates.
    4. Record leakage diagnostics per dataset and carry them into every row.
    5. Sweep radii and arms as in the synthetic experiment, using the heuristic
       radius path where the networks are too large for exact enumeration, and
       labelling those radii as upper bounds.
    6. Write results, aggregates, and a per-dataset comparison against the
       synthetic sweep at matched graph size.

    Args:
        cfg: The composed Hydra config.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    raise NotImplementedError


if __name__ == "__main__":
    main()
