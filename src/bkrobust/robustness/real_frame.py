"""Build and freeze the real-structure survival frame ([RE-11] Phase 1).

The **primary frame** for the real-network survival campaign is every row of
``results/axisa3/instances.jsonl`` with ``admissible == true``: 831 rows over
25 networks and 543 distinct ``(network, x, y)`` pairs, at coverage
1.0 / 0.5 / 0.25 (543 / 182 / 106). This module rebuilds each network's DAG and
CPDAG once, recomputes the analyst's graph ``G0`` once per ``(network,
coverage)``, then recomputes every per-row quantity the pre-registration
(``results/axis_robustness_real/PREREGISTRATION.md`` §1) requires, including a
full reproduction of the committed breakdown radius via
:func:`bkrobust.hybrid.breakdown_radius` ("T7"). It writes the frozen frame and
a self-describing hash manifest, and touches nothing else: the committed
corpus under ``results/axisa3/`` is read-only input, never rewritten.

Sentinel discipline, carried over from :mod:`bkrobust.benchmarks.measure` and
:mod:`bkrobust.core.conventions` because this module recomputes the same
quantities and must not reintroduce a bug that family already paid for:

* :data:`bkrobust.core.conventions.UNREACHED` (``-1``) is a status, never a
  radius; a row's ``r_status_recomputed`` names it explicitly and
  ``radius_recomputed`` is never averaged as ``-1``.
* A censored T7 reproduction (``not out.exact``) carries
  ``wall_until_timeout_s`` and sets ``radius_recomputed`` to ``null`` --
  **never** a key named ``seconds``, which is the exact shape of the defect
  [RE-5b]'s neighbours warn about.
* An undefined separation is ``null`` with a status string, never ``0`` and
  never ``-1``.
* [RE-5b]: ``benchmarks/measure.py::evaluate`` never assigns
  ``g0_undirected_edges`` on the admissible path, so every committed row
  carries the dataclass default ``0`` there -- not a measurement. This module
  recomputes the real value under a distinct column name
  (``g0_undirected_edges_recomputed``) and leaves the committed file alone.

Every knowledge-dependent condition that could in principle fail (``G0``
non-existent, the optimal set not identified, a GAC check failing at ``G0``) is
recorded as a status and a count in ``frame_hash.json``'s ``checks`` object,
never used to drop a row -- the frame's denominator is always 831.

Run as::

    PYTHONPATH=src /usr/bin/python3 -m bkrobust.robustness.real_frame
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from bkrobust.benchmarks import describe as D
from bkrobust.benchmarks.acquire import FileRecord, load_cached
from bkrobust.benchmarks.measure import (
    MAX_G0_UNDIRECTED_FOR_EXTENSIONS,
    component_of,
    select_knowledge,
    separation,
)
from bkrobust.core.conventions import UNREACHED
from bkrobust.core.resultsio import git_sha, is_dirty
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag
from bkrobust.demo.graph import MPDAG, undirected_components
from bkrobust.demo.meek import apply_orientations
from bkrobust.gac import is_gac_valid_mpdag
from bkrobust.hybrid import breakdown_radius

#: Source corpus this module reads and never writes.
DEFAULT_SOURCE = Path("results/axisa3/instances.jsonl")

#: Where the frozen frame and its hash manifest are written.
DEFAULT_OUT_DIR = Path("results/axis_robustness_real")

#: Expected shape of the primary frame, asserted before anything is written.
#: See ``PREREGISTRATION.md`` §1.1.
EXPECTED_N_ROWS = 831
EXPECTED_N_NETWORKS = 25
EXPECTED_N_PAIRS = 543
EXPECTED_ROWS_BY_COVERAGE = {"1.0": 543, "0.5": 182, "0.25": 106}

#: Status used on the row and in the checks when a knowledge-dependent
#: quantity could not be computed for a reason that "should not happen" on an
#: admissible row -- see the pre-registration and the module docstring.
STATUS_G0_NONE = "g0_none"
STATUS_Z_UNDEFINED = "optimal_set_undefined"

#: Exact row schema, in order. Every key is present on every row.
ROW_KEYS: tuple[str, ...] = (
    "row_id",
    "source_row_sha256",
    "network",
    "x",
    "y",
    "coverage",
    "largest_component_size",
    "component_size_committed",
    "component_size_recomputed",
    "n_k",
    "k_g0_committed",
    "k_g0_recomputed",
    "g0_undirected_edges_committed",
    "g0_undirected_edges_recomputed",
    "shd_truth_bw0",
    "separation_committed",
    "separation_status_committed",
    "separation_recomputed",
    "separation_status_recomputed",
    "z_star",
    "z_size_committed",
    "z_size",
    "z_valid_at_g0",
    "radius_committed",
    "r_status_committed",
    "radius_recomputed",
    "r_status_recomputed",
    "wall_until_timeout_s",
    "dispatch_leg_committed",
    "dispatch_leg",
    "oracle",
    "exact",
    "gate",
    "assumes",
    "radius_matches_committed",
)


def _sha256_bytes(data: bytes) -> str:
    """Hex sha256 of a bytes object.

    Args:
        data: The bytes to hash.

    Returns:
        The hex digest.
    """
    return hashlib.sha256(data).hexdigest()


def row_id(network: str, x: str, y: str, coverage: float) -> str:
    """The frame's stable row identifier.

    Args:
        network: Network name.
        x: Treatment.
        y: Outcome.
        coverage: Knowledge coverage level.

    Returns:
        The first 16 hex characters of
        ``sha256("|".join([network, x, y, repr(coverage)]))``.
    """
    key = "|".join([network, x, y, repr(coverage)])
    return hashlib.sha256(key.encode()).hexdigest()[:16]


def committed_r_status(radius: int, exact: bool) -> str:
    """The committed row's radius status, derived the same way the recomputed one is.

    Args:
        radius: The committed radius (may be :data:`UNREACHED`).
        exact: The committed ``exact`` flag.

    Returns:
        ``"unreached"`` if ``radius == UNREACHED``, ``"timeout"`` if
        ``exact`` is False, else ``"ok"``.
    """
    if radius == UNREACHED:
        return "unreached"
    if not exact:
        return "timeout"
    return "ok"


def load_source_rows(source_path: Path) -> list[tuple[dict[str, Any], bytes]]:
    """Read the admissible rows of the source corpus, each paired with its raw line.

    Iterates the file once, in file order, so ``source_row_sha256`` is the hash
    of the exact bytes that were on disk -- never a re-serialization.

    Args:
        source_path: Path to ``results/axisa3/instances.jsonl``.

    Returns:
        ``(parsed_row, raw_line_bytes)`` pairs, one per admissible row, raw
        bytes stripped of the trailing newline.
    """
    out: list[tuple[dict[str, Any], bytes]] = []
    raw = source_path.read_bytes()
    for line in raw.split(b"\n"):
        if not line.strip():
            continue
        parsed = json.loads(line)
        if parsed.get("_meta"):
            continue
        if not parsed.get("admissible"):
            continue
        out.append((parsed, line))
    return out


def resolve_network_file(network: str, recs: dict[str, FileRecord]) -> FileRecord:
    """Find the acquired file for a network name.

    Args:
        network: Network name as it appears in the source corpus.
        recs: File records keyed by base file name.

    Returns:
        The matching :class:`FileRecord`.

    Raises:
        KeyError: If none of the three known extensions match a cached file.
    """
    for ext in (".bif.gz", ".json", ".txt"):
        fname = f"{network}{ext}"
        if fname in recs:
            return recs[fname]
    raise KeyError(f"no cached file for network {network!r} (tried .bif.gz/.json/.txt)")


def load_graphs(network: str, recs: dict[str, FileRecord]) -> tuple[MPDAG, MPDAG]:
    """Parse a network once and build its DAG and CPDAG.

    Args:
        network: Network name.
        recs: File records keyed by base file name, from an :class:`Acquisition`.

    Returns:
        ``(dag, cpdag)``.
    """
    rec = resolve_network_file(network, recs)
    parsed = D.parse_file(rec.path, rec.sha256)
    dag = D.to_mpdag(parsed)
    cpdag = dag_to_cpdag(dag)
    return dag, cpdag


class NetworkCoverageStats:
    """The network-level quantities computed once per ``(network, coverage)``.

    Attributes:
        k: The analyst's asserted orientations (the STATED claim set).
        g0: The Meek closure of ``cpdag`` under ``k``, or ``None`` on FAIL.
        n_k: ``len(k)``.
        g0_undirected_edges_recomputed: ``len(g0.undirected_edges)``, or
            ``None`` if ``g0`` is ``None``.
        k_g0_recomputed: The commitment size, or ``None`` if ``g0`` is ``None``.
        shd_truth_bw0: ``len(g0.directed_edges ^ dag.directed_edges)``, or
            ``None`` if ``g0`` is ``None``.
    """

    __slots__ = ("k", "g0", "n_k", "g0_undirected_edges_recomputed", "k_g0_recomputed", "shd_truth_bw0")

    def __init__(
        self,
        k: list[tuple[str, str]],
        g0: MPDAG | None,
        n_k: int,
        g0_undirected_edges_recomputed: int | None,
        k_g0_recomputed: int | None,
        shd_truth_bw0: int | None,
    ) -> None:
        self.k = k
        self.g0 = g0
        self.n_k = n_k
        self.g0_undirected_edges_recomputed = g0_undirected_edges_recomputed
        self.k_g0_recomputed = k_g0_recomputed
        self.shd_truth_bw0 = shd_truth_bw0


def compute_network_coverage_stats(dag: MPDAG, cpdag: MPDAG, coverage: float) -> NetworkCoverageStats:
    """Compute the once-per-``(network, coverage)`` quantities.

    Args:
        dag: The ground-truth DAG.
        cpdag: Its CPDAG.
        coverage: Knowledge coverage level.

    Returns:
        A populated :class:`NetworkCoverageStats`. If ``G0`` fails to build
        (``apply_orientations`` returns ``None``), the recomputed fields are
        ``None`` rather than raising -- this "may in principle happen" per the
        pre-registration and is recorded as a status, not a crash.
    """
    k = select_knowledge(dag, cpdag, coverage)
    g0 = apply_orientations(cpdag, k)
    if g0 is None:
        return NetworkCoverageStats(k=k, g0=None, n_k=len(k), g0_undirected_edges_recomputed=None,
                                     k_g0_recomputed=None, shd_truth_bw0=None)
    k_g0_recomputed = sum(
        1 for (a, b) in g0.directed_edges if tuple(sorted((a, b))) in cpdag.undirected_edges
    )
    shd_truth_bw0 = len(g0.directed_edges ^ dag.directed_edges)
    return NetworkCoverageStats(
        k=k,
        g0=g0,
        n_k=len(k),
        g0_undirected_edges_recomputed=len(g0.undirected_edges),
        k_g0_recomputed=k_g0_recomputed,
        shd_truth_bw0=shd_truth_bw0,
    )


def build_row(
    source: dict[str, Any],
    raw_line: bytes,
    cpdag: MPDAG,
    stats: NetworkCoverageStats,
    largest_component_size: int,
    counters: dict[str, int],
    mismatch_examples: list[dict[str, Any]],
) -> dict[str, Any]:
    """Build one frame row from a source row plus its cached network state.

    Args:
        source: The parsed source JSONL row (admissible == true).
        raw_line: The exact source line bytes, no trailing newline.
        cpdag: The network's CPDAG.
        stats: The ``(network, coverage)`` state from
            :func:`compute_network_coverage_stats`.
        largest_component_size: The network's largest undirected-component size.
        counters: Mutated in place with the running check counts.
        mismatch_examples: Mutated in place; up to 10 full mismatching rows
            are appended here.

    Returns:
        The frame row, with every key in :data:`ROW_KEYS` present.
    """
    network, x, y, coverage = source["network"], source["x"], source["y"], source["coverage"]
    row: dict[str, Any] = {
        "row_id": row_id(network, x, y, coverage),
        "source_row_sha256": _sha256_bytes(raw_line),
        "network": network,
        "x": x,
        "y": y,
        "coverage": coverage,
        "largest_component_size": largest_component_size,
        "component_size_committed": source["component_size"],
        "component_size_recomputed": len(component_of(cpdag, x) or ()),
        "n_k": stats.n_k,
        "k_g0_committed": source["k_g0"],
        "k_g0_recomputed": stats.k_g0_recomputed,
        "g0_undirected_edges_committed": source["g0_undirected_edges"],
        "g0_undirected_edges_recomputed": stats.g0_undirected_edges_recomputed,
        "shd_truth_bw0": stats.shd_truth_bw0,
        "separation_committed": source["separation"],
        "separation_status_committed": source["separation_status"],
        "z_size_committed": source["z_size"],
        "radius_committed": source["radius"],
        "r_status_committed": committed_r_status(source["radius"], source["exact"]),
        "dispatch_leg_committed": source["method"],
        "gate": "fast_gate",
    }

    if source["g0_undirected_edges"] != 0:
        counters["g0_undirected_committed_nonzero"] += 1
    if stats.g0_undirected_edges_recomputed not in (None, 0):
        counters["n_g0_undirected_recomputed_nonzero"] += 1
    if stats.k_g0_recomputed is not None and stats.k_g0_recomputed != source["k_g0"]:
        counters["n_k_g0_mismatches"] += 1

    if stats.g0 is None:
        counters["n_g0_none"] += 1
        row.update(
            separation_recomputed=None,
            separation_status_recomputed=STATUS_G0_NONE,
            z_star=None,
            z_size=None,
            z_valid_at_g0=None,
            radius_recomputed=None,
            r_status_recomputed=STATUS_G0_NONE,
            wall_until_timeout_s=None,
            dispatch_leg=None,
            oracle=None,
            exact=None,
            assumes=None,
            radius_matches_committed=False,
        )
        return {k: row[k] for k in ROW_KEYS}

    g0 = stats.g0
    assert len(g0.undirected_edges) <= MAX_G0_UNDIRECTED_FOR_EXTENSIONS, (
        f"{network}/{coverage}: G0 carries {len(g0.undirected_edges)} undirected edges, "
        f"beyond the extension-enumeration limit -- should not happen on an admissible "
        f"(network, coverage) cell"
    )

    z = optimal_adjustment_set_mpdag(g0, x, y)
    if not z:
        counters["n_optimal_set_undefined"] += 1
        row.update(
            separation_recomputed=None,
            separation_status_recomputed=STATUS_Z_UNDEFINED,
            z_star=None,
            z_size=None,
            z_valid_at_g0=None,
            radius_recomputed=None,
            r_status_recomputed=STATUS_Z_UNDEFINED,
            wall_until_timeout_s=None,
            dispatch_leg=None,
            oracle=None,
            exact=None,
            assumes=None,
            radius_matches_committed=False,
        )
        return {k: row[k] for k in ROW_KEYS}

    z_star = sorted(z)
    z_frozen = frozenset(z_star)
    z_valid = is_gac_valid_mpdag(g0, x, y, z_frozen)
    if not z_valid:
        counters["n_z_invalid_at_g0"] += 1
    sep_recomputed, sep_status_recomputed = separation(cpdag, x, z_frozen)

    out = breakdown_radius(cpdag, None, x, y, z_frozen, g0=g0, time_limit_s=300.0)
    if not out.exact:
        counters["n_censored"] += 1
        radius_recomputed: int | None = None
        r_status_recomputed = "timeout"
        wall_until_timeout_s: float | None = out.total_seconds
    else:
        radius_recomputed = out.radius
        wall_until_timeout_s = None
        if out.radius == UNREACHED:
            counters["n_unreached"] += 1
            r_status_recomputed = "unreached"
        else:
            r_status_recomputed = "ok"

    radius_matches_committed = (
        radius_recomputed == row["radius_committed"] and r_status_recomputed == row["r_status_committed"]
    )
    if radius_matches_committed:
        counters["n_radius_matches_committed"] += 1
    else:
        counters["n_radius_mismatches"] += 1

    row.update(
        separation_recomputed=sep_recomputed,
        separation_status_recomputed=sep_status_recomputed,
        z_star=z_star,
        z_size=len(z_star),
        z_valid_at_g0=z_valid,
        radius_recomputed=radius_recomputed,
        r_status_recomputed=r_status_recomputed,
        wall_until_timeout_s=wall_until_timeout_s,
        dispatch_leg=out.method,
        oracle=out.oracle,
        exact=out.exact,
        assumes=out.assumes,
        radius_matches_committed=radius_matches_committed,
    )
    full_row = {k: row[k] for k in ROW_KEYS}
    if not radius_matches_committed and len(mismatch_examples) < 10:
        mismatch_examples.append(full_row)
    return full_row


def build_frame(source_path: Path) -> tuple[list[dict[str, Any]], dict[str, Any], bytes]:
    """Build every frame row from the committed source corpus.

    Args:
        source_path: Path to ``results/axisa3/instances.jsonl``.

    Returns:
        ``(rows, checks, source_bytes)`` where ``rows`` is sorted by
        ``(network, coverage, x, y)``, ``checks`` is the partially assembled
        checks object (missing only ``mismatch_examples`` promotion, which is
        already inlined), and ``source_bytes`` is the raw source file content.

    Raises:
        AssertionError: If the admissible-row counts do not match the
            pre-registered shape (§1.1).
    """
    source_bytes = source_path.read_bytes()
    admissible = load_source_rows(source_path)

    n_rows = len(admissible)
    networks = sorted({r["network"] for r, _ in admissible})
    pairs = {(r["network"], r["x"], r["y"]) for r, _ in admissible}

    assert n_rows == EXPECTED_N_ROWS, f"expected {EXPECTED_N_ROWS} admissible rows, found {n_rows}"
    assert len(networks) == EXPECTED_N_NETWORKS, (
        f"expected {EXPECTED_N_NETWORKS} networks, found {len(networks)}"
    )
    assert len(pairs) == EXPECTED_N_PAIRS, f"expected {EXPECTED_N_PAIRS} distinct pairs, found {len(pairs)}"
    rows_by_coverage: dict[str, int] = {}
    for r, _ in admissible:
        key = str(r["coverage"])
        rows_by_coverage[key] = rows_by_coverage.get(key, 0) + 1
    for cov_key, expected_n in EXPECTED_ROWS_BY_COVERAGE.items():
        got = rows_by_coverage.get(cov_key)
        assert got == expected_n, (
            f"coverage {cov_key}: expected {expected_n} admissible rows, found {got}"
        )

    acq = load_cached()
    recs = {f.name: f for f in acq.files}

    counters = {
        "n_radius_matches_committed": 0,
        "n_radius_mismatches": 0,
        "n_k_g0_mismatches": 0,
        "n_z_invalid_at_g0": 0,
        "n_optimal_set_undefined": 0,
        "n_g0_none": 0,
        "n_censored": 0,
        "n_unreached": 0,
        "g0_undirected_committed_nonzero": 0,
        "n_g0_undirected_recomputed_nonzero": 0,
    }
    mismatch_examples: list[dict[str, Any]] = []

    # Group admissible source rows by network, then by coverage, in sorted
    # order throughout -- never raw dict/set iteration order.
    by_network: dict[str, list[tuple[dict[str, Any], bytes]]] = {}
    for pair in admissible:
        by_network.setdefault(pair[0]["network"], []).append(pair)

    rows: list[dict[str, Any]] = []
    for network in networks:
        dag, cpdag = load_graphs(network, recs)
        largest_component_size = max(
            (len(c) for c in undirected_components(cpdag)), default=0
        )
        coverages = sorted({r["coverage"] for r, _ in by_network[network]})
        stats_by_coverage = {
            cov: compute_network_coverage_stats(dag, cpdag, cov) for cov in coverages
        }
        net_rows = sorted(
            by_network[network], key=lambda pair: (pair[0]["coverage"], pair[0]["x"], pair[0]["y"])
        )
        for source, raw_line in net_rows:
            stats = stats_by_coverage[source["coverage"]]
            rows.append(
                build_row(
                    source,
                    raw_line,
                    cpdag,
                    stats,
                    largest_component_size,
                    counters,
                    mismatch_examples,
                )
            )

    rows.sort(key=lambda r: (r["network"], r["coverage"], r["x"], r["y"]))

    checks = {
        "n_radius_matches_committed": counters["n_radius_matches_committed"],
        "n_radius_mismatches": counters["n_radius_mismatches"],
        "mismatch_examples": mismatch_examples,
        "n_k_g0_mismatches": counters["n_k_g0_mismatches"],
        "n_z_invalid_at_g0": counters["n_z_invalid_at_g0"],
        "n_optimal_set_undefined": counters["n_optimal_set_undefined"],
        "n_g0_none": counters["n_g0_none"],
        "n_censored": counters["n_censored"],
        "n_unreached": counters["n_unreached"],
        "g0_undirected_committed_all_zero": counters["g0_undirected_committed_nonzero"] == 0,
        "g0_undirected_recomputed_nonzero_count": counters["n_g0_undirected_recomputed_nonzero"],
    }
    return rows, checks, source_bytes


def write_frame(rows: list[dict[str, Any]], out_dir: Path) -> tuple[Path, bytes]:
    """Write ``frame.jsonl``, one JSON object per line in schema key order.

    Args:
        rows: Frame rows, already sorted by ``(network, coverage, x, y)``.
        out_dir: Destination directory, created if missing.

    Returns:
        ``(frame_path, frame_bytes)``.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    frame_path = out_dir / "frame.jsonl"
    lines = [json.dumps(row, sort_keys=False) for row in rows]
    frame_bytes = ("\n".join(lines) + "\n").encode("utf-8") if lines else b""
    frame_path.write_bytes(frame_bytes)
    return frame_path, frame_bytes


def write_frame_hash(
    out_dir: Path,
    frame_bytes: bytes,
    source_bytes: bytes,
    rows: list[dict[str, Any]],
    checks: dict[str, Any],
) -> Path:
    """Write ``frame_hash.json``, the frame's self-describing manifest.

    Args:
        out_dir: Destination directory.
        frame_bytes: The exact bytes written to ``frame.jsonl``.
        source_bytes: The exact bytes of the source ``instances.jsonl``.
        rows: The frame rows (for the row/network/pair/coverage counts).
        checks: The checks object assembled by :func:`build_frame`.

    Returns:
        The manifest path.
    """
    networks = sorted({r["network"] for r in rows})
    pairs = {(r["network"], r["x"], r["y"]) for r in rows}
    rows_by_coverage: dict[str, int] = {}
    for r in rows:
        key = str(r["coverage"])
        rows_by_coverage[key] = rows_by_coverage.get(key, 0) + 1

    payload = {
        "frame_sha256": _sha256_bytes(frame_bytes),
        "source_instances_sha256": _sha256_bytes(source_bytes),
        "n_rows": len(rows),
        "n_networks": len(networks),
        "n_distinct_pairs": len(pairs),
        "rows_by_coverage": rows_by_coverage,
        "built_utc": datetime.now(timezone.utc).isoformat(),  # noqa: UP017 -- py3.9
        "git_sha": git_sha(short=False),
        "git_dirty": is_dirty(),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "checks": checks,
    }
    path = out_dir / "frame_hash.json"
    path.write_text(json.dumps(payload, indent=2, default=str))
    return path


def main(argv: list[str] | None = None) -> int:
    """Build, freeze and hash the real-structure survival frame.

    Args:
        argv: Command-line arguments, or ``None`` to use ``sys.argv[1:]``.

    Returns:
        ``0`` on success. Assertion failures on the pre-registered shape
        (§1.1) propagate as exceptions rather than being caught, per the
        "fail loudly" requirement.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    args = parser.parse_args(argv)

    rows, checks, source_bytes = build_frame(args.source)
    frame_path, frame_bytes = write_frame(rows, args.out_dir)
    hash_path = write_frame_hash(args.out_dir, frame_bytes, source_bytes, rows, checks)

    summary = {
        "frame_path": str(frame_path),
        "hash_path": str(hash_path),
        "n_rows": len(rows),
        "frame_sha256": _sha256_bytes(frame_bytes),
        "source_instances_sha256": _sha256_bytes(source_bytes),
        "checks": checks,
    }
    print(json.dumps(summary, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
