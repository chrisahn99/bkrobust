"""Learned constraint-poset representations -- CONDITIONAL, may be cut.

Trains an embedding of the background-knowledge constraint poset and probes
whether its geometry predicts the breakdown radius or separates
consistent-and-true from consistent-but-false knowledge.

Two guards, both non-optional. Splits are by graph, never by constraint --
knowledge sets from one graph share its structure, so a random split leaks.
And the true/false probe is always reported next to its shuffle control: that
probe is trying to detect something the CPDAG by construction does not
determine, so a strong score is evidence about the sampler until the control
says otherwise.

If the probes do not beat the structural baselines in the config, this
experiment does not appear in the paper.

Usage::

    python experiments/run_knowledge_embedding.py experiment=knowledge_embedding
    python experiments/run_knowledge_embedding.py experiment=knowledge_embedding model=box_embedding
"""

from __future__ import annotations

import hydra
from omegaconf import DictConfig


@hydra.main(version_base=None, config_path="../configs", config_name="config")
def main(cfg: DictConfig) -> None:
    """Run the knowledge-embedding experiment.

    Steps once implemented:

    1. Seed and resolve the results directory.
    2. Generate graphs, sample knowledge sets of both arms, and label each with
       its radius, truth value and realised bias.
    3. Build the constraint poset per knowledge set.
    4. Split by graph into train/val/test.
    5. Fit the embedding named by the ``model`` config group.
    6. Run the probes, plus the structural baselines and the shuffle control.
    7. Write probe results with the baseline and control scores in the same
       table, so no score can be quoted without them.

    Args:
        cfg: The composed Hydra config.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    raise NotImplementedError


if __name__ == "__main__":
    main()
