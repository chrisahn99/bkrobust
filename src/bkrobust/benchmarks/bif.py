"""A structure-only parser for the BIF (Bayesian Interchange Format) files.

The 24 benchmark networks ship as gzipped ``.bif``. A full BIF parser would
read the conditional probability tables too, and every existing one we could
have used drags in a dependency tree we do not want for a structural question.
We need exactly two things:

* ``variable X { ... }``   -> ``X`` is a vertex;
* ``probability ( C | P1, P2 ) { ... }`` -> arcs ``P1 -> C``, ``P2 -> C``;
  a ``probability ( C )`` with no bar declares ``C`` a root.

So this module reads the declaration headers and ignores every table body. The
numbers are irrelevant to us and skipping them keeps the parser small enough to
audit by eye -- which matters, because a silently wrong parser would produce a
plausible-looking but fictitious network, the exact failure mode the whole
acquisition path is built to prevent.

Robustness notes, all driven by what the real files contain rather than by
guesswork: headers may carry arbitrary whitespace, ``|`` may be surrounded by
none, parent lists are comma-separated, and the corpus uses no comments (we
strip ``//`` and ``/* */`` anyway, since BIF permits them and stripping is
cheap). Every header in the corpus fits on one line, but the parser matches
across newlines so that a wrapped header would still be read correctly.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

from bkrobust.benchmarks.parsed import ParsedNetwork, build

_BLOCK_COMMENT = re.compile(r"/\*.*?\*/", re.DOTALL)
_LINE_COMMENT = re.compile(r"//[^\n]*")

#: ``variable NAME {`` -- the vertex declarations.
_VARIABLE = re.compile(r"\bvariable\s+([^\s{|,)]+)\s*\{", re.MULTILINE)

#: ``probability ( CHILD [| PARENTS] ) {`` -- the arc declarations.
_PROBABILITY = re.compile(r"\bprobability\s*\(\s*([^)]*?)\s*\)\s*\{", re.DOTALL)


def read_text(path: Path) -> str:
    """Read a ``.bif`` or ``.bif.gz`` file as text.

    Args:
        path: The file to read.

    Returns:
        The decoded contents. Latin-1 is the fallback for the few upstream
        files that are not valid UTF-8; it never fails, and node names in this
        corpus are ASCII regardless.
    """
    path = Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as handle:
            raw = handle.read()
    else:
        raw = path.read_bytes()
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")


def strip_comments(text: str) -> str:
    """Remove BIF block and line comments.

    Args:
        text: Raw BIF source.

    Returns:
        The source with comments replaced by nothing.
    """
    return _LINE_COMMENT.sub("", _BLOCK_COMMENT.sub("", text))


def _split_header(header: str) -> tuple[str, list[str]]:
    """Split a ``probability`` header's interior into child and parents.

    Args:
        header: The text between the parentheses, e.g. ``"C | P1, P2"``.

    Returns:
        ``(child, parents)``.

    Raises:
        ValueError: If the header names no child, or contains more than one
            ``|``.
    """
    parts = header.split("|")
    if len(parts) > 2:
        raise ValueError(f"probability header has multiple '|': {header!r}")
    child = parts[0].strip().strip(",").strip()
    if not child or any(ch.isspace() for ch in child):
        raise ValueError(f"unreadable child in probability header: {header!r}")
    parents: list[str] = []
    if len(parts) == 2:
        for token in parts[1].split(","):
            name = token.strip()
            if name:
                parents.append(name)
    return child, parents


def parse_bif_text(
    text: str,
    name: str,
    source_file: str,
    sha256: str,
) -> ParsedNetwork:
    """Parse BIF source into a :class:`ParsedNetwork`.

    Args:
        text: The BIF source (already decompressed).
        name: Short network name to record.
        source_file: Base name of the source file.
        sha256: Digest of the source file.

    Returns:
        The parsed network.

    Raises:
        ValueError: If no ``variable`` block is found, if a variable is
            declared twice, or if a ``probability`` header names a variable
            that was never declared.
    """
    body = strip_comments(text)

    nodes: list[str] = []
    seen: set[str] = set()
    for match in _VARIABLE.finditer(body):
        var = match.group(1)
        if var in seen:
            raise ValueError(f"variable declared twice: {var!r}")
        seen.add(var)
        nodes.append(var)
    if not nodes:
        raise ValueError("no 'variable' blocks found -- not a BIF file?")

    edges: list[tuple[str, str]] = []
    for match in _PROBABILITY.finditer(body):
        child, parents = _split_header(match.group(1))
        if child not in seen:
            raise ValueError(f"probability block for undeclared variable {child!r}")
        for parent in parents:
            if parent not in seen:
                raise ValueError(f"undeclared parent {parent!r} of {child!r}")
            edges.append((parent, child))

    return build(
        name=name,
        source_file=source_file,
        sha256=sha256,
        source_format="bif",
        raw_nodes=nodes,
        raw_edges=edges,
    )


def parse_bif_file(path: Path, sha256: str, name: str | None = None) -> ParsedNetwork:
    """Parse a ``.bif`` / ``.bif.gz`` file.

    Args:
        path: The file to parse.
        sha256: Digest of that file, from the acquisition manifest.
        name: Short network name; defaults to the file name with ``.bif``
            and ``.gz`` stripped.

    Returns:
        The parsed network.
    """
    path = Path(path)
    stem = path.name
    for suffix in (".gz", ".bif"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
    return parse_bif_text(
        read_text(path),
        name=name if name is not None else stem,
        source_file=path.name,
        sha256=sha256,
    )
