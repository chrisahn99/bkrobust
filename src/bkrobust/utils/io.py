"""Results directories and run manifests.

Results go to ``results/<experiment>/<timestamp>_<git-sha>/``. The git SHA is in
the path because a number without the code that produced it is not reproducible,
and a timestamp alone does not identify the code.

Every run directory holds a ``manifest.json``: the resolved config, the git SHA,
whether the working tree was dirty, package versions, platform, and wall-clock
timing. The dirty-tree flag matters more than it looks -- a result produced from
uncommitted changes cannot be reproduced from the SHA, and silently recording
the SHA anyway would be worse than recording nothing.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    import pandas as pd


def git_sha(short: bool = True) -> str:
    """Return the current commit SHA.

    Args:
        short: Return the abbreviated form.

    Returns:
        The SHA, or ``"nogit"`` if this is not a git working tree.
    """
    raise NotImplementedError


def is_dirty() -> bool:
    """Whether the working tree has uncommitted changes to tracked files."""
    raise NotImplementedError


def results_dir(
    experiment: str,
    *,
    root: str | Path = "results",
    run_id: str | None = None,
    create: bool = True,
) -> Path:
    """Resolve (and by default create) the directory for one run.

    Args:
        experiment: Experiment name, from ``cfg.experiment.name``.
        root: Results root.
        run_id: Override the generated ``<timestamp>_<sha>`` directory name.
        create: Create the directory.

    Returns:
        The path.

    Raises:
        FileExistsError: If the directory exists and holds results already.
            Overwriting a completed run silently is how results get lost.
    """
    raise NotImplementedError


def write_manifest(
    path: str | Path,
    cfg: Any,
    *,
    extra: dict[str, Any] | None = None,
) -> Path:
    """Write ``manifest.json`` into a run directory.

    Records the resolved config, git SHA and dirty flag, Python and package
    versions, platform, hostname-free environment summary, and the start time.

    Args:
        path: The run directory.
        cfg: The resolved Hydra config.
        extra: Additional fields to record.

    Returns:
        The manifest path.

    Note:
        Must not record usernames, home directories, hostnames or cluster names.
        The manifest is a tracked artefact of a double-blind submission; see
        ``docs/ANONYMITY_CHECKLIST.md``.
    """
    raise NotImplementedError


def save_results(
    path: str | Path,
    name: str,
    data: Any,
    *,
    fmt: str = "parquet",
) -> Path:
    """Write a results table into a run directory.

    Args:
        path: The run directory.
        name: Base filename without extension.
        data: A DataFrame, or anything convertible to one.
        fmt: ``"parquet"``, ``"csv"`` or ``"json"``. Parquet by default: these
            tables carry many float columns and CSV silently rounds them.

    Returns:
        The written path.
    """
    raise NotImplementedError


def load_results(path: str | Path, name: str | None = None) -> pd.DataFrame:
    """Read results back from a run directory.

    Args:
        path: The run directory.
        name: Which table; ``None`` concatenates every table in the directory,
            adding a column identifying which file each row came from.

    Returns:
        The results.

    Raises:
        FileNotFoundError: If nothing matches.
    """
    raise NotImplementedError


def latest_run(experiment: str, *, root: str | Path = "results") -> Path:
    """Return the most recent run directory for an experiment.

    Used by ``experiments/make_figures.py`` so that regenerating figures does
    not require pasting a timestamp.

    Raises:
        FileNotFoundError: If the experiment has no runs.
    """
    raise NotImplementedError
