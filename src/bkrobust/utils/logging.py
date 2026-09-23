"""Logging setup.

Hydra configures the root logger; this module gives the package a consistent
namespace under it and adds the two things the experiments actually need:
progress bars that do not fight with log lines, and a config dump at the top of
every run so a log file alone identifies what was run.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable, Iterator
from typing import TYPE_CHECKING, Any, TypeVar

if TYPE_CHECKING:  # pragma: no cover
    from pathlib import Path

T = TypeVar("T")


def get_logger(name: str) -> logging.Logger:
    """Return a logger under the ``bkrobust`` namespace.

    Args:
        name: Usually ``__name__``.

    Returns:
        The logger.
    """
    raise NotImplementedError


def setup_logging(
    level: str = "INFO",
    *,
    rich: bool = False,
    log_file: str | Path | None = None,
) -> None:
    """Configure package logging.

    Args:
        level: Log level name.
        rich: Use rich formatting when available; plain otherwise. Off by
            default -- rich output in a redirected log file is unreadable.
        log_file: Also write to this file.
    """
    raise NotImplementedError


def log_config(cfg: Any, logger: logging.Logger | None = None) -> None:
    """Log the fully resolved config at the start of a run.

    Resolved, not raw: interpolations must be expanded, or the log records what
    was written rather than what ran.
    """
    raise NotImplementedError


def progress(
    iterable: Iterable[T],
    *,
    desc: str = "",
    total: int | None = None,
    enabled: bool = True,
) -> Iterator[T]:
    """Wrap an iterable in a progress bar, respecting ``compute.progress``.

    Args:
        iterable: What to iterate.
        desc: Bar description.
        total: Length, when the iterable does not report one.
        enabled: Pass ``False`` to disable -- in tests, and in non-interactive
            runs where the bar just fills the log with control characters.

    Yields:
        The items.
    """
    raise NotImplementedError
