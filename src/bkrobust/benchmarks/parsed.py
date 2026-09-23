"""The format-independent result of parsing one benchmark network.

Three on-disk formats (BIF, dagitty, bnlearn JSON) feed one downstream sweep,
so they agree on a single container. It carries only what the structural
question needs -- vertices and arcs -- plus the provenance (file name, sha256)
that lets a row of the results table be traced back to specific bytes.

Node-name normalisation lives here rather than in the parsers because it must
be applied identically across formats, and because it must be *reversible*: the
mapping from normalised name back to the name as it appears in the file is
recorded, so a claim about "node X of ALARM" can always be checked against the
source.
"""

from __future__ import annotations

from dataclasses import dataclass, field

#: Characters kept verbatim in a normalised node name.
_SAFE = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-")


def normalise_names(raw_names: list[str]) -> dict[str, str]:
    """Map raw node names to CSV/JSON-safe names, reversibly and deterministically.

    Any character outside ``[A-Za-z0-9_.-]`` becomes ``_``. If two raw names
    collide after substitution, later ones (in sorted order) get a ``__2``,
    ``__3``, ... suffix, so the map stays injective and therefore invertible.

    Args:
        raw_names: The names exactly as they appear in the source file.

    Returns:
        A dict from raw name to normalised name. It is the identity map when
        every name is already safe, which is the common case.

    Raises:
        ValueError: If a name is empty, or if the same raw name appears twice.
    """
    mapping: dict[str, str] = {}
    used: set[str] = set()
    for raw in sorted(set(raw_names)):
        if not raw:
            raise ValueError("empty node name")
        base = "".join(ch if ch in _SAFE else "_" for ch in raw)
        candidate = base
        suffix = 2
        while candidate in used:
            candidate = f"{base}__{suffix}"
            suffix += 1
        used.add(candidate)
        mapping[raw] = candidate
    return mapping


@dataclass(frozen=True)
class ParsedNetwork:
    """A benchmark network's structure, as read from one file.

    Args:
        name: Short network name (the file's stem, e.g. ``alarm``).
        source_file: Base name of the file it was read from.
        sha256: Digest of that file's bytes.
        source_format: ``bif``, ``dagitty`` or ``bnjson``.
        nodes: Normalised node names, sorted.
        edges: Directed arcs as ``(parent, child)``, normalised, sorted.
        bidirected: Bidirected (``<->``) pairs, canonically ordered and sorted.
            Non-empty only for dagitty files that encode latent confounding;
            such a file is parsed successfully but is not a DAG.
        name_map: Normalised name -> raw name, recorded only where the two
            differ, so normalisation is auditable and reversible.
    """

    name: str
    source_file: str
    sha256: str
    source_format: str
    nodes: tuple[str, ...]
    edges: tuple[tuple[str, str], ...]
    bidirected: tuple[tuple[str, str], ...] = ()
    name_map: dict[str, str] = field(default_factory=dict)

    @property
    def n_nodes(self) -> int:
        """Number of vertices."""
        return len(self.nodes)

    @property
    def n_edges(self) -> int:
        """Number of directed arcs."""
        return len(self.edges)


def build(
    name: str,
    source_file: str,
    sha256: str,
    source_format: str,
    raw_nodes: list[str],
    raw_edges: list[tuple[str, str]],
    raw_bidirected: list[tuple[str, str]] | None = None,
) -> ParsedNetwork:
    """Normalise names and assemble a :class:`ParsedNetwork`.

    Args:
        name: Short network name.
        source_file: Base name of the source file.
        sha256: Digest of the source file.
        source_format: Parser tag.
        raw_nodes: Vertices as named in the file.
        raw_edges: Arcs ``(parent, child)`` as named in the file.
        raw_bidirected: Bidirected pairs as named in the file.

    Returns:
        The assembled network, with sorted, de-duplicated node and edge tuples.

    Raises:
        ValueError: If an edge references a node not in ``raw_nodes``, or is a
            self-loop.
    """
    mapping = normalise_names(list(raw_nodes))
    known = set(mapping)
    edges: set[tuple[str, str]] = set()
    for parent, child in raw_edges:
        if parent not in known or child not in known:
            raise ValueError(f"edge {(parent, child)} references an undeclared node")
        if parent == child:
            raise ValueError(f"self-loop on {parent!r}")
        edges.add((mapping[parent], mapping[child]))

    bidirected: set[tuple[str, str]] = set()
    for a, b in raw_bidirected or []:
        if a not in known or b not in known:
            raise ValueError(f"bidirected edge {(a, b)} references an undeclared node")
        if a == b:
            raise ValueError(f"self-loop on {a!r}")
        x, y = mapping[a], mapping[b]
        bidirected.add((x, y) if x <= y else (y, x))

    return ParsedNetwork(
        name=name,
        source_file=source_file,
        sha256=sha256,
        source_format=source_format,
        nodes=tuple(sorted(mapping.values())),
        edges=tuple(sorted(edges)),
        bidirected=tuple(sorted(bidirected)),
        name_map={v: k for k, v in sorted(mapping.items()) if v != k},
    )
