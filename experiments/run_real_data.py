"""Real-data certificates -- what a practitioner can actually compute.

No ground-truth DAG, so no measured radii and no measured bias. What this
experiment produces is the certificate: given the knowledge an analyst supplied,
how far it can be wrong before the conclusion moves, which constraints are
load-bearing, and how much the effect estimate varies across the ball of
knowledge that passes the same consistency check.

The discovery step is bootstrapped. The CPDAG is an estimate, and a certificate
computed against a single point estimate of it would understate the uncertainty
by exactly the amount that discovery error contributes.

Usage::

    python experiments/run_real_data.py experiment=real_data
    python experiments/run_real_data.py experiment=real_data experiment.dataset=ihdp
"""

from __future__ import annotations

import hydra
from omegaconf import DictConfig


@hydra.main(version_base=None, config_path="../configs", config_name="config")
def main(cfg: DictConfig) -> None:
    """Run the real-data certificate experiment.

    Steps once implemented:

    1. Seed and resolve the results directory.
    2. Load the dataset and its variable glossary.
    3. Read the analyst's knowledge from ``cfg.experiment.knowledge``, or elicit
       it if a source is configured, and record its provenance verbatim.
    4. Bootstrap the discovery step to get a distribution over CPDAGs.
    5. Check consistency on each; a knowledge set that FAILs on some bootstrap
       replicates and not others is itself a finding worth reporting, since the
       analyst would have seen only one draw.
    6. Certify on each replicate, and report the spread of radii and verdicts.
    7. Estimate the effect over the ball of consistent knowledge and report the
       spread of conclusions it admits.
    8. Write the certificates, the aggregate table, and the rendered reports.

    Args:
        cfg: The composed Hydra config.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    raise NotImplementedError


if __name__ == "__main__":
    main()
