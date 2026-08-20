"""Semi-synthetic generators: real covariates, simulated outcomes.

The compromise tier. Fully synthetic data has known ground truth and marginals
nothing like real data; real data has realistic marginals and no ground truth.
Semi-synthetic generators take real covariates and simulate the outcome from a
specified mechanism, so the true effect is known while the covariate
distribution is not something anyone made up.

``DREAM``
    Gene-regulatory network challenges. Networks with published gold standards,
    plus simulators, at sizes where exact radius computation is still feasible
    at the small end.
``RealCause``
    Generative models fitted to real causal benchmarks, sampled to produce data
    with realistic marginals and a known effect.

The caveat that must ride along with every number from this tier: the true graph
here is the *simulator's* graph, and it is true of the simulated outcome, not of
the real process the covariates came from. A radius measured here says how
robust the method is on data with realistic marginals, not how robust it is in
the domain the covariates were collected from. See ``docs/DATA.md``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    import numpy as np


def load_dream(
    dataset: str,
    *,
    data_dir: str | None = None,
) -> dict[str, Any]:
    """Load one DREAM challenge network and its data.

    Args:
        dataset: Dataset identifier, e.g. ``"dream4_size10_1"``.
        data_dir: Override the data directory (``data/``, gitignored).

    Returns:
        Keys ``"dag"``, ``"data"``, ``"nodes"``, ``"metadata"``.

    Raises:
        FileNotFoundError: If the data is absent, with a pointer to
            ``scripts/download_data.sh``.
    """
    raise NotImplementedError


def load_realcause(
    dataset: str,
    *,
    data_dir: str | None = None,
    n_samples: int | None = None,
    seed: int = 0,
) -> dict[str, Any]:
    """Sample from a RealCause-style generator fitted to a real benchmark.

    Args:
        dataset: Which fitted generator, e.g. ``"lalonde_cps"``.
        data_dir: Override the data directory.
        n_samples: How many rows to draw; ``None`` uses the generator default.
        seed: Sampling seed.

    Returns:
        Keys ``"data"``, ``"true_ate"``, ``"true_ite"``, ``"metadata"``.

    Raises:
        FileNotFoundError: If the generator artefacts are absent.
    """
    raise NotImplementedError


def simulate_outcome(
    covariates: Any,
    dag: Any,
    rng: np.random.Generator,
    *,
    mechanism: str = "linear",
    effect_size: float = 1.0,
    noise_scale: float = 1.0,
) -> Any:
    """Simulate treatment and outcome on top of real covariates.

    Args:
        covariates: Real covariate matrix.
        dag: The DAG specifying which covariates enter which mechanism.
        rng: Seeded generator.
        mechanism: ``"linear"``, ``"nonlinear"`` or ``"interaction"``.
        effect_size: True treatment effect.
        noise_scale: Outcome noise standard deviation.

    Returns:
        A DataFrame with the covariates plus the simulated treatment and outcome.
    """
    raise NotImplementedError


def available_datasets(generator: str) -> list[str]:
    """List the dataset identifiers a generator provides.

    Raises:
        KeyError: If ``generator`` is not ``"dream"`` or ``"realcause"``.
    """
    raise NotImplementedError
