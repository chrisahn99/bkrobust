"""Meek amplification measurement -- how far one wrong assertion travels.

Measures the ratio of forced orientations to imposed constraints across graph
families, densities and constraint kinds. The mechanism behind everything else:
if a single false assertion routinely forces a dozen orientations, then small
knowledge errors have large structural consequences and the breakdown radii will
be small even on large graphs.

Separately worth knowing is *where* the forced orientations land. Amplification
that never reaches the parents of the causal nodes cannot move ``O*`` and costs
nothing, so the study records cascade reach and hit rate alongside the raw
factor.

Usage::

    python experiments/run_cascade_study.py
    python experiments/run_cascade_study.py -m graph=erdos_renyi,scale_free
"""

from __future__ import annotations

import hydra
from omegaconf import DictConfig


@hydra.main(version_base=None, config_path="../configs", config_name="config")
def main(cfg: DictConfig) -> None:
    """Run the cascade amplification study.

    Steps once implemented:

    1. Seed and resolve the results directory.
    2. For each graph family and density: sample DAGs and their CPDAGs.
    3. For ``k = 1, 2, ...`` impose ``k`` consistent constraints of each kind
       and record the cascade report -- forced count, reach, per-rule
       attribution, and whether ``O*`` moved.
    4. Aggregate the amplification factor by family, density and kind, and
       report R4's share separately: R4 only fires once knowledge is present, so
       its contribution is the part of the cascade that exists *because* the
       analyst supplied knowledge.

    Args:
        cfg: The composed Hydra config.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    raise NotImplementedError


if __name__ == "__main__":
    main()
