"""Regenerate every paper figure from the latest results.

Reads from ``results/<experiment>/`` and writes to ``paper/figures/``. Takes the
most recent run per experiment unless a run id is given, so the usual invocation
needs no arguments.

Figures are regenerated, never edited by hand. A figure that cannot be rebuilt
from a results directory and a git SHA cannot be checked by a reviewer or by
whoever picks this up in six months.

Usage::

    python experiments/make_figures.py
    python experiments/make_figures.py output.run_id=2026-09-01_12-00-00_abc1234
"""

from __future__ import annotations

import hydra
from omegaconf import DictConfig


@hydra.main(version_base=None, config_path="../configs", config_name="config")
def main(cfg: DictConfig) -> None:
    """Rebuild the paper figures.

    Figures produced once implemented:

    * ``breakdown_curve`` -- bias and variance inflation against radius, with
      both radii marked and the band shaded. The main figure.
    * ``severity_stack`` -- composition of the knowledge ball by severity class
      as the radius grows.
    * ``cascade_amplification`` -- forced orientations against imposed
      constraints, by graph family.
    * ``radius_ordering`` -- ``delta_opt`` against ``delta_valid`` across graphs,
      with the diagonal drawn. Any point above the diagonal is a counterexample
      to the conjecture and must be plotted individually rather than aggregated.
    * ``cfm_audit`` -- bias against radius per checkpoint, with the classical arm.
    * ``certificate_examples`` -- rendered certificates for the real datasets.

    Args:
        cfg: The composed Hydra config.

    Raises:
        NotImplementedError: Always -- this is a scaffold.
    """
    raise NotImplementedError


if __name__ == "__main__":
    main()
