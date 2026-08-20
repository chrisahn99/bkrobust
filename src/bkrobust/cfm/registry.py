"""Checkpoint registry for the audited causal foundation models.

Every reported number must be traceable to a specific checkpoint at a specific
revision. These models are updated, and a result produced against "the CausalPFN
checkpoint" without a revision is not reproducible -- so :class:`CheckpointSpec`
requires a revision and the audit refuses to run without one.

Nothing here downloads anything. ``scripts/download_checkpoints.sh`` is the only
place a download may be initiated, and it is deliberately a stub: licences vary
across these checkpoints, and pulling them silently as a side effect of an
import would be wrong regardless of what the licence says. See
``docs/CHECKPOINTS.md``.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CheckpointSpec:
    """Everything needed to locate and cite one checkpoint.

    Attributes:
        key: Registry key, matching the config name -- ``"causalpfn"``,
            ``"causalfm"``, ``"dopfn"``.
        display_name: Name as it appears in the paper's tables.
        source: Where it comes from -- a model hub id or release URL.
        revision: Commit hash or release tag. Required: an unpinned checkpoint
            makes the result unreproducible.
        licence: SPDX identifier or licence name. Recorded because these
            checkpoints are not uniformly permissive and the paper must state
            terms accurately.
        estimands: Which estimands this checkpoint supports -- ``"ate"``,
            ``"cate"``.
        supports_conditioning: Which conditioning modes it accepts -- ``"none"``,
            ``"attention_bias"``, ``"ancestral_matrix"``, ``"prompt"``.
        max_features: Feature-count limit, where the architecture imposes one.
        max_samples: In-context sample limit, likewise. Both are binding for
            these models and determine which datasets are auditable at all.
        notes: Caveats -- known preprocessing requirements, expected input
            scaling, quirks worth recording next to a number.
    """

    key: str
    display_name: str
    source: str
    revision: str
    licence: str
    estimands: tuple[str, ...]
    supports_conditioning: tuple[str, ...]
    max_features: int | None
    max_samples: int | None
    notes: str = ""


def register(spec: CheckpointSpec) -> None:
    """Add a checkpoint to the registry.

    Raises:
        ValueError: If the key is already registered, or ``revision`` is empty.
    """
    raise NotImplementedError


def get(key: str) -> CheckpointSpec:
    """Look up a checkpoint by registry key.

    Raises:
        KeyError: If unknown; the message lists the registered keys.
    """
    raise NotImplementedError


def available() -> list[str]:
    """Registered checkpoint keys, sorted."""
    raise NotImplementedError


def local_path(key: str, cache_dir: str | None = None) -> str:
    """Resolve where a checkpoint should live on disk.

    Does not download. Reports where the file is expected and, if absent, points
    at ``scripts/download_checkpoints.sh``.

    Args:
        key: Registry key.
        cache_dir: Override the default cache directory (``checkpoints/``,
            gitignored).

    Returns:
        Absolute path.

    Raises:
        FileNotFoundError: If nothing is at that path, with instructions.
    """
    raise NotImplementedError


def cite(key: str) -> str:
    """Return the BibTeX key for a checkpoint's paper.

    Keeps the audit tables and ``paper/references.bib`` from drifting apart.

    Raises:
        KeyError: If the checkpoint is unknown.
    """
    raise NotImplementedError
