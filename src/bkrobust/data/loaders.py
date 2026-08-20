"""Real-data loaders.

No ground-truth DAG anywhere in this module, and therefore no breakdown radii.
What these datasets support is the certificate: given the knowledge an analyst
supplies, how much of their conclusion rests on it. That is the only honest
question to ask of real data, and it happens to be the one practitioners
actually have.

Two families, with different reasons to be here:

*Structure benchmarks* -- Sachs, Perturb-seq -- have a consensus network,
sometimes with interventional data to check it against. The consensus is not
ground truth; it is expert belief that has held up, which is a different thing
and is worth stating in the paper rather than in a footnote.

*Effect benchmarks* -- IHDP, ACIC, Twins, Jobs, LaLonde -- have a known or
partly known effect, usually because the outcome was simulated or because a
randomised arm exists. They exercise the estimation path without a graph.

Every loader documents provenance and licence in ``docs/DATA.md``. None of them
downloads: they read from ``data/``, which is gitignored, and point at
``scripts/download_data.sh`` when a file is missing. Some of these datasets have
redistribution terms that a convenience auto-download would quietly violate.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    import pandas as pd


def load_sachs(
    *,
    data_dir: str | None = None,
    regime: str = "observational",
    log_transform: bool = True,
    standardize: bool = False,
) -> dict[str, Any]:
    """Load the Sachs protein-signalling data and its consensus network.

    Args:
        data_dir: Override the data directory.
        regime: ``"observational"``, ``"interventional"`` or ``"pooled"``.
        log_transform: Apply a log transform; the raw abundances are heavily
            skewed and most published analyses transform them.
        standardize: Standardise columns. Off by default -- standardising
            destroys the variance information that varsortability diagnostics
            read, and hides a leak rather than removing it.

    Returns:
        Keys ``"data"``, ``"consensus_dag"``, ``"nodes"``, ``"metadata"``.

    Raises:
        FileNotFoundError: If the data is absent.
    """
    raise NotImplementedError


def load_ihdp(
    *,
    data_dir: str | None = None,
    replication: int = 0,
) -> dict[str, Any]:
    """Load one IHDP replication -- real covariates, simulated outcomes.

    Args:
        data_dir: Override the data directory.
        replication: Which of the 1000 replications.

    Returns:
        Keys ``"data"``, ``"true_ate"``, ``"true_ite"``, ``"metadata"``.
    """
    raise NotImplementedError


def load_acic(
    *,
    data_dir: str | None = None,
    year: int = 2016,
    setting: int = 1,
) -> dict[str, Any]:
    """Load an ACIC data challenge setting.

    Returns:
        Keys ``"data"``, ``"true_ate"``, ``"metadata"``.
    """
    raise NotImplementedError


def load_twins(*, data_dir: str | None = None) -> dict[str, Any]:
    """Load the Twins dataset.

    Both potential outcomes are observed -- one twin each -- so the individual
    effect is known without simulation, which is rare and makes this the best
    PEHE benchmark available here.

    Returns:
        Keys ``"data"``, ``"true_ite"``, ``"metadata"``.
    """
    raise NotImplementedError


def load_jobs(*, data_dir: str | None = None) -> dict[str, Any]:
    """Load the Jobs dataset (randomised arm plus observational arm).

    Returns:
        Keys ``"data"``, ``"randomized_mask"``, ``"metadata"``.
    """
    raise NotImplementedError


def load_lalonde(
    *,
    data_dir: str | None = None,
    control: str = "psid",
) -> dict[str, Any]:
    """Load the LaLonde data with the requested observational control group.

    Args:
        data_dir: Override the data directory.
        control: ``"psid"``, ``"cps"`` or ``"experimental"``.

    Returns:
        Keys ``"data"``, ``"experimental_ate"``, ``"metadata"``.
    """
    raise NotImplementedError


def load_perturb_seq(
    *,
    data_dir: str | None = None,
    subset: str | None = None,
) -> dict[str, Any]:
    """Load Perturb-seq data: single-cell expression with genetic interventions.

    Large and sparse. The interventional arm gives partial structural ground
    truth -- knocking out a gene reveals its downstream targets -- which is the
    closest thing to a verified graph in this module.

    Args:
        data_dir: Override the data directory.
        subset: Gene subset identifier; ``None`` loads the full matrix, which is
            far too large for exact radius computation.

    Returns:
        Keys ``"data"``, ``"interventions"``, ``"nodes"``, ``"metadata"``.
    """
    raise NotImplementedError


def load(dataset: str, **kwargs: Any) -> dict[str, Any]:
    """Dispatch to the loader named by ``dataset``.

    Raises:
        KeyError: If unknown; the message lists the available names.
    """
    raise NotImplementedError


def variable_glossary(dataset: str) -> dict[str, str]:
    """Human-readable descriptions of a dataset's variables.

    Needed by :mod:`bkrobust.knowledge.elicit`: an LLM asked to orient ``"praf"``
    against ``"pmek"`` needs to know what those are, and the quality of the
    elicited knowledge depends on what it was told. The glossary is therefore
    part of the experimental condition and belongs in the run manifest.
    """
    raise NotImplementedError
