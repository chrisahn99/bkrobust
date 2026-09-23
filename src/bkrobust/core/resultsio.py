"""The results contract: append-only rows plus a self-describing manifest. Frozen.

Two rules, both learned the hard way:

* **Append incrementally.** A crashed run must never lose more than the row it
  was computing. :class:`ResultWriter` flushes after every row.
* **Every results directory is self-describing.** The manifest records the root
  seed, the git SHA, the environment, the parameter grid and the radius
  convention, so a CSV can always be traced to the code and settings that made
  it.
"""

from __future__ import annotations

import csv
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from bkrobust.core.conventions import RADIUS_CONVENTION


def git_sha(short: bool = True) -> str:
    """Current commit SHA, or ``"nogit"`` outside a working tree."""
    try:
        args = ["git", "rev-parse", "--short" if short else "HEAD", "HEAD"]
        if not short:
            args = ["git", "rev-parse", "HEAD"]
        out = subprocess.run(args, capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return "nogit"


def is_dirty() -> bool:
    """Whether tracked files have uncommitted changes."""
    try:
        out = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=no"],
            capture_output=True,
            text=True,
            check=True,
        )
        return bool(out.stdout.strip())
    except Exception:
        return False


def write_manifest(
    directory: str | Path,
    *,
    seed: int,
    grid: dict[str, Any],
    extra: dict[str, Any] | None = None,
) -> Path:
    """Write ``manifest.json`` describing a results directory.

    Args:
        directory: The results directory.
        seed: The root seed every draw descends from.
        grid: The parameter grid this run swept.
        extra: Anything else worth recording.

    Returns:
        The manifest path.
    """
    d = Path(directory)
    d.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        # datetime.UTC is 3.11+; this session runs on 3.9.
        "created_utc": datetime.now(timezone.utc).isoformat(),  # noqa: UP017
        "git_sha": git_sha(short=False),
        "git_dirty": is_dirty(),
        "radius_convention": RADIUS_CONVENTION,
        "seed": seed,
        "grid": grid,
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
        },
    }
    for mod in ("numpy", "networkx", "pandas", "scipy", "matplotlib"):
        try:
            payload["environment"][mod] = __import__(mod).__version__
        except Exception:
            payload["environment"][mod] = "absent"
    if extra:
        payload.update(extra)
    path = d / "manifest.json"
    path.write_text(json.dumps(payload, indent=2, default=str))
    return path


class ResultWriter:
    """Append-only CSV writer that flushes after every row.

    The header is fixed by the first row. Later rows with extra keys raise
    rather than silently dropping columns -- a schema drift mid-run would make
    the file unreadable and is always a bug.

    Args:
        path: Destination CSV.
        resume: If ``True`` and the file exists, append without rewriting the
            header. If ``False``, an existing file is an error, so a run cannot
            silently clobber earlier results.

    Raises:
        FileExistsError: If ``path`` exists and ``resume`` is ``False``.
    """

    def __init__(self, path: str | Path, *, resume: bool = False) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        if self.path.exists() and not resume:
            raise FileExistsError(f"{self.path} exists; pass resume=True to append")
        self._fh = None
        self._writer = None
        self._fields: list[str] | None = None
        if self.path.exists() and resume:
            with self.path.open() as fh:
                header = fh.readline().strip()
            if header:
                self._fields = header.split(",")

    def write(self, row: dict[str, Any]) -> None:
        """Append one row and flush.

        Raises:
            ValueError: If the row's keys do not match the established header.
        """
        if self._fields is None:
            self._fields = list(row.keys())
            self._fh = self.path.open("w", newline="")
            self._writer = csv.DictWriter(self._fh, fieldnames=self._fields)
            self._writer.writeheader()
        elif self._fh is None:
            self._fh = self.path.open("a", newline="")
            self._writer = csv.DictWriter(self._fh, fieldnames=self._fields)
        missing = set(self._fields) - set(row)
        extra = set(row) - set(self._fields)
        if missing or extra:
            raise ValueError(f"row schema drift: missing={sorted(missing)} extra={sorted(extra)}")
        self._writer.writerow(row)
        self._fh.flush()

    def close(self) -> None:
        """Close the underlying file."""
        if self._fh is not None:
            self._fh.close()
            self._fh = None
            self._writer = None

    def __enter__(self) -> ResultWriter:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()
