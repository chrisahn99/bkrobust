"""A parser for the dagitty-style ``.txt`` DAGs from applied papers.

These twelve files are a different kind of object from the BIF benchmarks: they
are the causal diagrams drawn *by hand, by domain experts, in published applied
papers* (Shrier & Platt 2008 on sports injury, Schipf 2010 on diabetes,
Sebastiani 2005 on stroke, ...), exported from the dagitty editor. They are
small, but they are the closest thing available to "the structure an analyst
actually asserts", which is exactly the population the background-knowledge
question is about.

The format, as it appears in these files:

    dag {
    bb="-3,-0.5,2,1.2"
    Name [pos="1.0,2.0"]
    Other [exposure,pos="-2.0,1.0"]
    Name -> Other
    Name <-> Other
    }

so: an optional bounding-box line, one declaration line per node carrying
bracketed attributes (``exposure``, ``outcome``, ``latent``, ``adjusted``,
``pos``), then one edge line per arc.

``<->`` is a *bidirected* edge -- dagitty's notation for an unmeasured common
cause -- and is emphatically not a DAG arc. Exactly one file in the corpus
(``M-bias.txt``) uses it. It is parsed and recorded separately rather than
silently coerced into an arc, and the descriptive sweep then declines to treat
that file as a DAG. Coercing it would manufacture a network nobody drew.

``--`` (undirected) does not occur in these files; it is accepted and recorded
as bidirected-style non-arc information would be, i.e. it makes the graph
non-DAG, so that a file using it is never quietly mis-read as a DAG.
"""

from __future__ import annotations

import re
from pathlib import Path

from bkrobust.benchmarks.parsed import ParsedNetwork, build

#: A node declaration: a bare name, optionally followed by ``[attributes]``.
_NODE_DECL = re.compile(r'^(?P<name>(?:"[^"]+"|[^\s\[\]{}]+))\s*(?:\[(?P<attrs>[^\]]*)\])?$')

#: An edge line: ``A -> B``, ``A <- B``, ``A <-> B`` or ``A -- B``, with an
#: optional trailing ``[...]`` attribute block.
_EDGE = re.compile(
    r'^(?P<left>(?:"[^"]+"|[^\s\[\]{}]+))\s*'
    r"(?P<arrow><->|->|<-|--)\s*"
    r'(?P<right>(?:"[^"]+"|[^\s\[\]{}]+))\s*'
    r"(?:\[[^\]]*\])?$"
)

#: Graph-level attribute lines such as ``bb="-3,-0.5,2,1.2"``.
_GRAPH_ATTR = re.compile(r'^[A-Za-z_][A-Za-z0-9_]*\s*=\s*(?:"[^"]*"|\S+)$')


class DagittyParseError(ValueError):
    """Raised when a ``.txt`` file does not match the dagitty grammar above."""


def _unquote(token: str) -> str:
    """Strip surrounding double quotes from a dagitty identifier."""
    quoted = len(token) >= 2 and token.startswith('"') and token.endswith('"')
    return token[1:-1] if quoted else token


def parse_dagitty_text(
    text: str,
    name: str,
    source_file: str,
    sha256: str,
) -> ParsedNetwork:
    """Parse dagitty source into a :class:`ParsedNetwork`.

    Nodes are taken from declaration lines *and* from edge endpoints, since
    dagitty does not require a node to be declared before it is used.

    Args:
        text: The file contents.
        name: Short network name to record.
        source_file: Base name of the source file.
        sha256: Digest of the source file.

    Returns:
        The parsed network, with ``<->`` and ``--`` edges in ``bidirected``.

    Raises:
        DagittyParseError: If the file has no ``dag {`` / ``graph {`` header,
            or if any line matches none of the recognised forms. Failing loudly
            is deliberate: a skipped line is a silently missing edge.
    """
    lines = [line.strip() for line in text.splitlines()]
    lines = [line for line in lines if line and not line.startswith("#")]
    if not lines:
        raise DagittyParseError("empty file")

    header = lines[0]
    if not re.match(r"^(dag|graph|digraph|pdag|mag|pag)\b.*\{$", header):
        raise DagittyParseError(f"unrecognised header line: {header!r}")
    body = lines[1:]
    if body and body[-1] == "}":
        body = body[:-1]
    else:
        raise DagittyParseError("file does not end with '}'")

    nodes: list[str] = []
    seen: set[str] = set()
    edges: list[tuple[str, str]] = []
    bidirected: list[tuple[str, str]] = []

    def note(node: str) -> None:
        if node not in seen:
            seen.add(node)
            nodes.append(node)

    for line in body:
        stripped = line.rstrip(";").strip()
        if not stripped:
            continue
        if _GRAPH_ATTR.match(stripped):
            continue
        edge = _EDGE.match(stripped)
        if edge is not None:
            left = _unquote(edge.group("left"))
            right = _unquote(edge.group("right"))
            note(left)
            note(right)
            arrow = edge.group("arrow")
            if arrow == "->":
                edges.append((left, right))
            elif arrow == "<-":
                edges.append((right, left))
            else:
                bidirected.append((left, right))
            continue
        decl = _NODE_DECL.match(stripped)
        if decl is not None:
            note(_unquote(decl.group("name")))
            continue
        raise DagittyParseError(f"unrecognised line: {stripped!r}")

    if not nodes:
        raise DagittyParseError("no nodes found")

    return build(
        name=name,
        source_file=source_file,
        sha256=sha256,
        source_format="dagitty",
        raw_nodes=nodes,
        raw_edges=edges,
        raw_bidirected=bidirected,
    )


def parse_dagitty_file(path: Path, sha256: str, name: str | None = None) -> ParsedNetwork:
    """Parse a dagitty ``.txt`` file.

    Args:
        path: The file to parse.
        sha256: Digest of that file, from the acquisition manifest.
        name: Short network name; defaults to the file stem.

    Returns:
        The parsed network.
    """
    path = Path(path)
    return parse_dagitty_text(
        path.read_text(encoding="utf-8", errors="replace"),
        name=name if name is not None else path.stem,
        source_file=path.name,
        sha256=sha256,
    )
