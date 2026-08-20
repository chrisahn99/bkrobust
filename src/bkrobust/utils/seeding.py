"""Seeding, and the discipline around it.

One rule: every source of randomness in an experiment descends from the single
``seed`` in the config. No module may call a global RNG -- not ``np.random.rand``,
not ``random.random``, not an unseeded sklearn default. Generators are created
here and passed down explicitly.

The reason is that this project's results are distributions over sampled
knowledge sets, and a distribution whose sampling cannot be replayed is not a
result anyone can check. The nuisance case is nested parallelism: replicates
running concurrently must not share a stream, and must not be reseeded with
``seed + i``, which correlates the streams. :func:`spawn` uses NumPy's
``SeedSequence`` spawning, which gives independent streams with a reproducible
derivation.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    import numpy as np


def set_global_seed(seed: int) -> None:
    """Seed the global RNGs of Python, NumPy and torch if present.

    A backstop for third-party code that reaches for a global RNG, not a
    substitute for passing generators explicitly. Experiment code must still
    thread a :class:`numpy.random.Generator` through.

    Args:
        seed: The seed.
    """
    raise NotImplementedError


def get_rng(seed: int) -> np.random.Generator:
    """Create a NumPy generator from a seed."""
    raise NotImplementedError


def spawn(seed: int, n: int) -> list[np.random.Generator]:
    """Derive ``n`` independent generators from one seed.

    Uses ``SeedSequence.spawn``, so the streams are statistically independent
    and the derivation is deterministic -- replicate ``i`` gets the same stream
    whether the sweep runs serially or across processes.

    Args:
        seed: The root seed.
        n: How many streams.

    Returns:
        ``n`` independent generators.

    Raises:
        ValueError: If ``n`` is not positive.
    """
    raise NotImplementedError


def seed_from_config(cfg: object) -> int:
    """Read the root seed out of a Hydra config.

    Raises:
        KeyError: If the config has no ``seed``. There is no default: an
            experiment that silently picks its own seed is unreproducible, and
            failing loudly is cheaper than discovering that later.
    """
    raise NotImplementedError


def seed_grid(seed: int, n_seeds: int) -> Sequence[int]:
    """Expand a root seed into the seeds for ``n_seeds`` full replications.

    Distinct from :func:`spawn`: this produces *seeds* for whole independent
    runs of an experiment, each of which then spawns its own generators.
    """
    raise NotImplementedError
