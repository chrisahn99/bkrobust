"""The descriptive structural sweep over the benchmark corpus.

This is the phase that has to happen *before* any breakdown-radius measurement,
and be committed before it, so the corpus cannot be filtered after the fact to
whichever networks flatter a hypothesis. It computes, per network, only what is
a property of the graph itself:

* how big it is (nodes, edges) and whether it really is a DAG;
* its CPDAG -- the Markov equivalence class, i.e. what a discovery algorithm
  could have recovered without background knowledge;
* the **fraction of edges the CPDAG leaves undirected**, which is precisely the
  budget that background knowledge has to work with. No undirected edges means
  no room for false knowledge to do damage;
* the **chain components** of that undirected subgraph, with the *full* size
  distribution rather than just the maximum. The distribution is what matters:
  a network with one component of size 40 and one with twenty components of
  size 2 have the same maximum edge budget and completely different geometry.
* how many components have **three or more vertices** -- those are the only
  ones that can host a separation of 2 or more, so this count is the corpus's
  upper bound on where the interesting failures can live at all.

Cost note. ``dag_to_cpdag`` is used exactly as the repository defines it, with
no local acceleration, because the descriptive numbers must be the numbers that
module produces. Its Meek closure re-derives every rule from scratch after each
single orientation, and ``MPDAG.parents``/``adjacent`` scan the whole edge set,
so the cost grows sharply with the number of undirected edges. Each network is
therefore run under a wall-clock budget (:data:`DEFAULT_BUDGET_S`) and a network
that exceeds it is recorded as timed out, never silently dropped.
"""

from __future__ import annotations

import csv
import json
import signal
import time
from collections import Counter
from dataclasses import asdict, dataclass, field
from pathlib import Path
from types import FrameType

from bkrobust.benchmarks import bif, bnjson, dagitty
from bkrobust.benchmarks.acquire import DEFAULT_CACHE, Acquisition, acquire
from bkrobust.benchmarks.parsed import ParsedNetwork
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG, undirected_components

#: Wall-clock budget per network for the CPDAG computation, in seconds.
DEFAULT_BUDGET_S = 300.0

#: Column order of ``descriptive_structure.csv``.
CSV_COLUMNS = (
    "network",
    "source_file",
    "sha256",
    "source_format",
    "n_nodes",
    "n_edges",
    "is_dag",
    "n_undirected_edges",
    "undirected_fraction",
    "n_components",
    "max_component_size",
    "n_components_ge3",
    "n_components_ge6",
    "n_nodes_in_components",
    "seconds_to_cpdag",
)


class CpdagTimeoutError(RuntimeError):
    """Raised when a network's CPDAG exceeds its wall-clock budget."""


@dataclass(frozen=True)
class NetworkDescription:
    """The descriptive row for one network.

    Args:
        network: Short network name.
        source_file: File it was parsed from.
        sha256: Digest of that file.
        source_format: ``bif``, ``dagitty`` or ``bnjson``.
        n_nodes: Vertex count.
        n_edges: Arc count.
        is_dag: Whether the parsed graph is a fully oriented acyclic graph.
        n_undirected_edges: Undirected edges in the CPDAG.
        undirected_fraction: ``n_undirected_edges / n_edges``, or 0.0 when the
            graph has no edges.
        n_components: Number of chain components of size >= 2 (singletons are
            not components in the sense that matters here: an isolated vertex
            carries no orientation choice).
        max_component_size: Largest chain component, 0 if there are none.
        n_components_ge3: Chain components with at least three vertices -- the
            ones capable of separation >= 2.
        n_components_ge6: Chain components with at least six vertices, the
            scale at which the synthetic families topped out.
        n_nodes_in_components: Vertices lying in some chain component.
        component_sizes: Full size distribution, sorted descending.
        seconds_to_cpdag: Wall-clock seconds the CPDAG took.
        status: ``ok``, ``timeout``, or ``not_a_dag``.
        note: Human-readable reason when ``status`` is not ``ok``.
    """

    network: str
    source_file: str
    sha256: str
    source_format: str
    n_nodes: int
    n_edges: int
    is_dag: bool
    n_undirected_edges: int | None
    undirected_fraction: float | None
    n_components: int | None
    max_component_size: int | None
    n_components_ge3: int | None
    n_components_ge6: int | None
    n_nodes_in_components: int | None
    component_sizes: tuple[int, ...] = ()
    seconds_to_cpdag: float | None = None
    status: str = "ok"
    note: str = ""
    name_map: dict[str, str] = field(default_factory=dict)


def to_mpdag(network: ParsedNetwork) -> MPDAG:
    """Convert a parsed network to the repository's :class:`MPDAG`.

    Only the directed arcs are carried across. Bidirected (``<->``) edges,
    which encode latent confounding, have no representation in ``MPDAG`` and
    are deliberately *not* coerced into arcs or into undirected edges; callers
    check ``network.bidirected`` and treat such a graph as non-DAG.

    Args:
        network: The parsed network.

    Returns:
        The MPDAG with every edge directed.
    """
    return MPDAG(network.nodes, directed=network.edges)


class _Budget:
    """A SIGALRM-based wall-clock guard for a single pure-Python computation.

    ``dag_to_cpdag`` has no cancellation hook, so the budget is enforced by an
    interval timer whose handler raises inside the running code. This works
    only on the main thread of a POSIX process; elsewhere the guard degrades to
    a no-op and the caller simply waits, which is acceptable because the budget
    is a reporting device, not a correctness one.

    Args:
        seconds: The budget. Non-positive disables the guard.
        label: Name used in the raised error.
    """

    def __init__(self, seconds: float, label: str) -> None:
        self._seconds = seconds
        self._label = label
        self._armed = False

    def _fire(self, signum: int, frame: FrameType | None) -> None:
        raise CpdagTimeoutError(f"{self._label}: exceeded {self._seconds:g}s budget")

    def __enter__(self) -> _Budget:
        """Arm the timer, if this process and thread support one."""
        if self._seconds <= 0:
            return self
        try:
            self._previous = signal.signal(signal.SIGALRM, self._fire)
            signal.setitimer(signal.ITIMER_REAL, self._seconds)
            self._armed = True
        except (ValueError, AttributeError, OSError):  # pragma: no cover - platform guard
            self._armed = False
        return self

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        """Disarm the timer and restore the previous handler."""
        if self._armed:
            signal.setitimer(signal.ITIMER_REAL, 0.0)
            signal.signal(signal.SIGALRM, self._previous)
            self._armed = False


def describe_network(
    network: ParsedNetwork,
    budget_s: float = DEFAULT_BUDGET_S,
) -> NetworkDescription:
    """Compute the descriptive row for one parsed network.

    Args:
        network: The parsed network.
        budget_s: Wall-clock budget for the CPDAG computation.

    Returns:
        The description. On timeout or on a non-DAG input the CPDAG-derived
        fields are ``None`` and ``status`` says why.
    """
    graph = to_mpdag(network)
    is_dag = graph.is_dag() and not network.bidirected

    base = {
        "network": network.name,
        "source_file": network.source_file,
        "sha256": network.sha256,
        "source_format": network.source_format,
        "n_nodes": network.n_nodes,
        "n_edges": network.n_edges,
        "is_dag": is_dag,
        "name_map": dict(network.name_map),
    }

    if not is_dag:
        if network.bidirected:
            note = (
                f"{len(network.bidirected)} bidirected (<->) edge(s): latent confounding, "
                "not representable as a DAG"
            )
        else:
            note = "directed graph contains a cycle"
        return NetworkDescription(
            n_undirected_edges=None,
            undirected_fraction=None,
            n_components=None,
            max_component_size=None,
            n_components_ge3=None,
            n_components_ge6=None,
            n_nodes_in_components=None,
            seconds_to_cpdag=None,
            status="not_a_dag",
            note=note,
            **base,
        )

    start = time.perf_counter()
    try:
        with _Budget(budget_s, network.name):
            cpdag = dag_to_cpdag(graph)
    except CpdagTimeoutError:
        return NetworkDescription(
            n_undirected_edges=None,
            undirected_fraction=None,
            n_components=None,
            max_component_size=None,
            n_components_ge3=None,
            n_components_ge6=None,
            n_nodes_in_components=None,
            seconds_to_cpdag=round(time.perf_counter() - start, 3),
            status="timeout",
            note=f"dag_to_cpdag exceeded the {budget_s:g}s budget",
            **base,
        )
    elapsed = time.perf_counter() - start

    n_undirected = len(cpdag.undirected_edges)
    components = undirected_components(cpdag)
    sizes = sorted((len(c) for c in components), reverse=True)

    return NetworkDescription(
        n_undirected_edges=n_undirected,
        undirected_fraction=(n_undirected / network.n_edges) if network.n_edges else 0.0,
        n_components=len(sizes),
        max_component_size=max(sizes) if sizes else 0,
        n_components_ge3=sum(1 for s in sizes if s >= 3),
        n_components_ge6=sum(1 for s in sizes if s >= 6),
        n_nodes_in_components=sum(sizes),
        component_sizes=tuple(sizes),
        seconds_to_cpdag=round(elapsed, 3),
        status="ok",
        **base,
    )


def parse_file(path: Path, sha256: str) -> ParsedNetwork:
    """Dispatch a file to the parser for its format.

    Args:
        path: The file to parse.
        sha256: Digest of that file.

    Returns:
        The parsed network.

    Raises:
        ValueError: If the extension is not one of the three known formats.
    """
    name = path.name
    if name.endswith(".bif.gz") or name.endswith(".bif"):
        return bif.parse_bif_file(path, sha256)
    if name.endswith(".txt"):
        return dagitty.parse_dagitty_file(path, sha256)
    if name.endswith(".json"):
        return bnjson.parse_bnjson_file(path, sha256)
    raise ValueError(f"no parser for {name!r}")


def parse_corpus(
    acquisition: Acquisition,
) -> tuple[list[ParsedNetwork], list[tuple[str, str]]]:
    """Parse every network file in an acquisition.

    Args:
        acquisition: The acquisition record naming the extracted files.

    Returns:
        ``(networks, failures)`` where ``failures`` pairs a file name with the
        reason it could not be parsed. Iteration is over the acquisition's
        sorted file list, so the result is deterministic.
    """
    networks: list[ParsedNetwork] = []
    failures: list[tuple[str, str]] = []
    for record in acquisition.files:
        if not (
            record.name.endswith((".bif", ".bif.gz", ".txt", ".json"))
            and not record.name.startswith("_")
        ):
            continue
        try:
            networks.append(parse_file(record.path, record.sha256))
        except Exception as exc:  # a parse failure is data to report, not a crash
            failures.append((record.name, f"{type(exc).__name__}: {exc}"))
    return sorted(networks, key=lambda n: n.name), sorted(failures)


def sweep(
    acquisition: Acquisition,
    budget_s: float = DEFAULT_BUDGET_S,
) -> tuple[list[NetworkDescription], list[tuple[str, str]]]:
    """Parse and describe the whole corpus.

    Args:
        acquisition: The acquisition record.
        budget_s: Per-network CPDAG budget.

    Returns:
        ``(descriptions, parse_failures)``, descriptions sorted by network name.
    """
    networks, failures = parse_corpus(acquisition)
    return [describe_network(n, budget_s=budget_s) for n in networks], failures


def write_csv(descriptions: list[NetworkDescription], path: Path) -> None:
    """Write the descriptive table as CSV.

    Args:
        descriptions: The rows, in the order they should appear.
        path: Destination file; parent directories are created.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(CSV_COLUMNS))
        writer.writeheader()
        for row in descriptions:
            record = asdict(row)
            writer.writerow({k: record[k] for k in CSV_COLUMNS})


def write_json(
    descriptions: list[NetworkDescription],
    failures: list[tuple[str, str]],
    acquisition: Acquisition,
    path: Path,
) -> None:
    """Write the descriptive table as JSON, with the size distributions.

    The component-size distribution is the reason this file exists: it does not
    fit a CSV cell, and summarising it to a maximum throws away exactly the
    shape the downstream question cares about.

    Args:
        descriptions: The rows.
        failures: Files that could not be parsed, with reasons.
        acquisition: Provenance for the corpus.
        path: Destination file; parent directories are created.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for row in descriptions:
        record = asdict(row)
        record["component_sizes"] = list(row.component_sizes)
        record["component_size_counts"] = {
            str(size): count for size, count in sorted(Counter(row.component_sizes).items())
        }
        rows.append(record)
    payload = {
        "provenance": {
            "sdist_url": acquisition.sdist_url,
            "sdist_sha256": acquisition.sdist_sha256,
            "cache_dir": str(acquisition.cache_dir),
            "n_files_extracted": len(acquisition.files),
        },
        "n_networks": len(rows),
        "n_parse_failures": len(failures),
        "parse_failures": [{"file": name, "reason": reason} for name, reason in failures],
        "networks": rows,
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main(
    cache_dir: Path = DEFAULT_CACHE,
    out_dir: Path = Path("results/axisa3"),
    sdist_path: Path | None = None,
    budget_s: float = DEFAULT_BUDGET_S,
) -> list[NetworkDescription]:
    """Run the whole pipeline: acquire, parse, describe, write both tables.

    Args:
        cache_dir: Where the sdist and extracted networks are cached.
        out_dir: Where ``descriptive_structure.{csv,json}`` are written.
        sdist_path: An explicit sdist path to prefer over the cache.
        budget_s: Per-network CPDAG budget.

    Returns:
        The descriptive rows.
    """
    acquisition = acquire(cache_dir=cache_dir, sdist_path=sdist_path)
    descriptions, failures = sweep(acquisition, budget_s=budget_s)
    write_csv(descriptions, Path(out_dir) / "descriptive_structure.csv")
    write_json(descriptions, failures, acquisition, Path(out_dir) / "descriptive_structure.json")
    return descriptions


if __name__ == "__main__":  # pragma: no cover - manual entry point
    import argparse

    argument_parser = argparse.ArgumentParser(description=__doc__)
    argument_parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    argument_parser.add_argument("--out-dir", type=Path, default=Path("results/axisa3"))
    argument_parser.add_argument("--sdist", type=Path, default=None)
    argument_parser.add_argument("--budget-s", type=float, default=DEFAULT_BUDGET_S)
    args = argument_parser.parse_args()
    rows = main(
        cache_dir=args.cache_dir,
        out_dir=args.out_dir,
        sdist_path=args.sdist,
        budget_s=args.budget_s,
    )
    for row in rows:
        print(
            f"{row.network:16s} n={row.n_nodes:5d} e={row.n_edges:5d} "
            f"undir={row.n_undirected_edges} frac="
            f"{'' if row.undirected_fraction is None else format(row.undirected_fraction, '.3f')} "
            f"maxcomp={row.max_component_size} ge3={row.n_components_ge3} "
            f"t={row.seconds_to_cpdag}s [{row.status}]"
        )
