"""A parser for the four bnlearn-style JSON networks in the same corpus.

``arth150``, ``ecoli70``, ``magic-irri`` and ``magic-niab`` ship as JSON rather
than BIF because they are Gaussian networks, whose parameters do not fit BIF's
discrete table syntax. The *structure* is stored plainly:

    {"nodes": ["a", "b", ...], "arcs": [["a", "b"], ...]}

with each arc a ``[parent, child]`` pair. That is all the structural sweep
needs, so they are included: ``arth150`` (107 nodes) and ``ecoli70`` are real
gene-expression networks and widen the corpus beyond expert-elicited discrete
models.

The parser is strict about shape and refuses anything that is not this exact
schema, rather than trying to divine an alternative layout.
"""

from __future__ import annotations

import json
from pathlib import Path

from bkrobust.benchmarks.parsed import ParsedNetwork, build


class BnJsonParseError(ValueError):
    """Raised when a ``.json`` file is not the ``{nodes, arcs}`` schema."""


def parse_bnjson_text(
    text: str,
    name: str,
    source_file: str,
    sha256: str,
) -> ParsedNetwork:
    """Parse bnlearn-style JSON into a :class:`ParsedNetwork`.

    Args:
        text: The file contents.
        name: Short network name to record.
        source_file: Base name of the source file.
        sha256: Digest of the source file.

    Returns:
        The parsed network.

    Raises:
        BnJsonParseError: If the payload is not an object with a list of
            string ``nodes`` and a list of two-element string ``arcs``.
    """
    try:
        payload = json.loads(text)
    except json.JSONDecodeError as exc:
        raise BnJsonParseError(f"invalid JSON: {exc}") from exc
    if not isinstance(payload, dict) or "nodes" not in payload or "arcs" not in payload:
        raise BnJsonParseError("expected an object with 'nodes' and 'arcs'")

    raw_nodes = payload["nodes"]
    raw_arcs = payload["arcs"]
    if not isinstance(raw_nodes, list) or not all(isinstance(n, str) for n in raw_nodes):
        raise BnJsonParseError("'nodes' is not a list of strings")
    if not isinstance(raw_arcs, list):
        raise BnJsonParseError("'arcs' is not a list")

    edges: list[tuple[str, str]] = []
    for arc in raw_arcs:
        if (
            not isinstance(arc, (list, tuple))
            or len(arc) != 2
            or not all(isinstance(x, str) for x in arc)
        ):
            raise BnJsonParseError(f"arc is not a [parent, child] string pair: {arc!r}")
        edges.append((arc[0], arc[1]))

    return build(
        name=name,
        source_file=source_file,
        sha256=sha256,
        source_format="bnjson",
        raw_nodes=list(raw_nodes),
        raw_edges=edges,
    )


def parse_bnjson_file(path: Path, sha256: str, name: str | None = None) -> ParsedNetwork:
    """Parse a bnlearn-style ``.json`` file.

    Args:
        path: The file to parse.
        sha256: Digest of that file, from the acquisition manifest.
        name: Short network name; defaults to the file stem.

    Returns:
        The parsed network.
    """
    path = Path(path)
    return parse_bnjson_text(
        path.read_text(encoding="utf-8"),
        name=name if name is not None else path.stem,
        source_file=path.name,
        sha256=sha256,
    )
