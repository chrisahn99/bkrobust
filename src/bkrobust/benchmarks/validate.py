"""An independent count of a BIF file's nodes and edges, for cross-checking.

:mod:`bkrobust.benchmarks.bif` builds a graph; this module counts. The two share
no code beyond the standard library, on purpose. If the BIF parser had a subtle
bug -- a regex that skipped headers spanning a newline, a parent list split on
the wrong character -- it would produce a network that looks entirely plausible
in a results table, and nothing downstream would catch it. A second, dumber
implementation that agrees on the totals is cheap insurance against exactly
that.

The counting path here is deliberately naive: iterate the decompressed text
line by line, count lines that start with ``variable``, and for lines that start
with ``probability`` count the comma-separated names after the ``|``. No shared
regexes, no shared helpers, no gzip helper in common.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class RawCounts:
    """Counts obtained by the independent path.

    Args:
        n_variable_blocks: Lines opening a ``variable`` block.
        n_probability_blocks: Lines opening a ``probability`` block.
        n_parent_references: Total parent names across all ``probability``
            headers -- i.e. the arc count.
        n_root_blocks: ``probability`` headers with no ``|``.
    """

    n_variable_blocks: int
    n_probability_blocks: int
    n_parent_references: int
    n_root_blocks: int


def count_bif(path: Path) -> RawCounts:
    """Count nodes and parent references in a BIF file, line by line.

    Args:
        path: A ``.bif`` or ``.bif.gz`` file.

    Returns:
        The counts.
    """
    path = Path(path)
    opener = gzip.open if str(path).endswith(".gz") else open
    n_var = 0
    n_prob = 0
    n_parents = 0
    n_root = 0
    with opener(path, "rt", encoding="latin-1") as handle:  # type: ignore[operator]
        for line in handle:
            text = line.strip()
            if text.startswith("variable "):
                n_var += 1
            elif text.startswith("probability"):
                n_prob += 1
                inner = text[text.index("(") + 1 : text.rindex(")")]
                if "|" not in inner:
                    n_root += 1
                    continue
                tail = inner.split("|", 1)[1]
                n_parents += len([piece for piece in tail.split(",") if piece.strip()])
    return RawCounts(
        n_variable_blocks=n_var,
        n_probability_blocks=n_prob,
        n_parent_references=n_parents,
        n_root_blocks=n_root,
    )
