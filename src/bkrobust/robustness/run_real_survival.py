"""Sharded, resumable driver for the real-structure survival sweep ([RE-11]).

Design fixed by ``results/axis_robustness_real/PREREGISTRATION.md``. This module
is orchestration and I/O only; the corruption and scoring logic lives in
:mod:`bkrobust.robustness.real_survival`.

Incrementality contract, and it is the whole point of the file layout
---------------------------------------------------------------------
Work is cut into **shards** at ``(arm, network, cell, replicate)`` granularity.
A shard owns exactly **two** output files and writes them append-only, one row
at a time, flushed. When -- and only when -- a shard finishes, it writes
``_done/<shard_id>.json`` carrying its row counts and the SHA-256 of both files.

A shard with no completion marker is **not resumed**: its files are truncated
and it is redone from scratch. Partial output is therefore cleanly discardable
by construction, which is the property session 8's survival sweep lacked when a
relaunched worker overwrote its own outputs in place. Re-invoking any stage is
safe to repeat and skips completed shards.

Usage::

    PYTHONPATH=src /usr/bin/python3 -m bkrobust.robustness.run_real_survival list
    PYTHONPATH=src /usr/bin/python3 -m bkrobust.robustness.run_real_survival run --shard <id>
    PYTHONPATH=src /usr/bin/python3 -m bkrobust.robustness.run_real_survival status

No global RNG is touched anywhere. Every draw is seeded by a SHA-256 derivation
in :mod:`~bkrobust.robustness.real_survival`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import statistics
import time
from pathlib import Path
from typing import Any, Iterator

import numpy as np

from bkrobust.core.resultsio import git_sha, is_dirty
from bkrobust.robustness import real_survival as rs
from bkrobust.synth.knowledge import flip as _flip

DEFAULT_OUT = Path("results/axis_robustness_real")

#: Coverage levels that mirror a synthetic stratum. Coverage 0.25 is run and
#: reported as supplementary, never as one of the nine.
PRIMARY_COVERAGES: tuple[float, ...] = (1.0, 0.5)
SUPPLEMENTARY_COVERAGES: tuple[float, ...] = (0.25,)

#: Base wrongness levels, mirroring the synthetic design.
BASE_WRONGNESS: tuple[float, ...] = (0.0, 0.10, 0.25)

#: Analyst replicates. ``b = 0`` is the identity corruption, so a second draw
#: would be a duplicate rather than a replicate.
REPLICATES_BY_BW: dict[float, int] = {0.0: 1, 0.10: 3, 0.25: 3}

#: Absolute base-wrongness levels, in **claims reversed** rather than a rate.
#: Added by Appendix A of the pre-registration, before any survival curve
#: existed, because the rate-based mechanism reverses nothing at all in 24 of
#: the corpus's 50 ``(network, coverage)`` cells at ``b = 0.25`` and 38 of 50 at
#: ``b = 0.10``: ``|K|`` is 1 or 2 on most of this corpus and
#: ``round(b·|K|) == 0``. Without an absolute level the ``SHD(G0, truth)``
#: baseline stays identically zero on those cells, which is exactly the rigged
#: comparison the design set out to avoid. Implemented through the same
#: ``synth.knowledge.flip`` call at ``rate = n/|K|``, so no new corruption
#: operator enters the design.
BASE_WRONGNESS_ABS: tuple[int, ...] = (1,)
REPLICATES_ABS: int = 3

#: Tier counts, mirroring the synthetic design.
N_TIERS: tuple[int, ...] = (2, 3, 4)


# --- frame --------------------------------------------------------------------


def load_frame(out_dir: Path) -> list[dict[str, Any]]:
    """Read the frozen frame and verify its hash before anything runs on it.

    Args:
        out_dir: The results subtree holding ``frame.jsonl`` and
            ``frame_hash.json``.

    Returns:
        The frame rows, in file order.

    Raises:
        FileNotFoundError: If the frame or its hash file is absent.
        RuntimeError: If the frame's digest differs from the recorded one.
    """
    fpath, hpath = out_dir / "frame.jsonl", out_dir / "frame_hash.json"
    if not fpath.is_file() or not hpath.is_file():
        raise FileNotFoundError(f"frozen frame missing under {out_dir}")
    recorded = json.loads(hpath.read_text())["frame_sha256"]
    digest = hashlib.sha256(fpath.read_bytes()).hexdigest()
    if digest != recorded:
        raise RuntimeError(f"frame.jsonl digest {digest} != recorded {recorded}")
    return [json.loads(line) for line in fpath.read_text().splitlines() if line.strip()]


# --- shards -------------------------------------------------------------------


def shard_id(spec: dict[str, Any]) -> str:
    """The deterministic id, and filename stem, of a shard.

    Args:
        spec: A shard specification.

    Returns:
        The id.
    """
    kind = spec["kind"]
    if kind == "flip":
        bw = spec.get("base_wrongness")
        tag = f"bw{round(bw * 100):03d}" if bw is not None else f"bwa{spec['bw_abs']}"
        return (
            f"flip__{spec['network']}__cov{round(spec['coverage'] * 100):03d}"
            f"__{tag}__r{spec['replicate']}"
        )
    if kind == "tiered":
        return f"tiered__{spec['network']}__nt{spec['n_tiers']}"
    return f"xarm__{spec['network']}__nt{spec['n_tiers']}__{spec['xarm']}"


def enumerate_shards(frame: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Every unit of work the campaign consists of, in a deterministic order.

    Args:
        frame: The frozen frame.

    Returns:
        Shard specifications, sorted by id.
    """
    networks = sorted({r["network"] for r in frame})
    cov_by_net: dict[str, set[float]] = {}
    for r in frame:
        cov_by_net.setdefault(r["network"], set()).add(r["coverage"])

    specs: list[dict[str, Any]] = []
    for net in networks:
        for cov in sorted(cov_by_net[net]):
            for bw in BASE_WRONGNESS:
                for rep in range(REPLICATES_BY_BW[bw]):
                    specs.append(
                        {"kind": "flip", "network": net, "coverage": cov,
                         "base_wrongness": bw, "bw_abs": None, "replicate": rep}
                    )
            for n_abs in BASE_WRONGNESS_ABS:
                for rep in range(REPLICATES_ABS):
                    specs.append(
                        {"kind": "flip", "network": net, "coverage": cov,
                         "base_wrongness": None, "bw_abs": n_abs, "replicate": rep}
                    )
        for nt in N_TIERS:
            specs.append({"kind": "tiered", "network": net, "n_tiers": nt})
            for arm in ("tiered", "flip"):
                specs.append({"kind": "xarm", "network": net, "n_tiers": nt, "xarm": arm})
    return sorted(specs, key=shard_id)


def is_complete(out_dir: Path, sid: str) -> bool:
    """Whether a shard has a completion marker.

    Args:
        out_dir: The results subtree.
        sid: The shard id.

    Returns:
        ``True`` if and only if the marker exists.
    """
    return (out_dir / "_done" / f"{sid}.json").is_file()


# --- writing ------------------------------------------------------------------


class JsonlWriter:
    """Append-only JSONL writer that flushes after every row.

    A shard truncates its file on open, because a shard without a marker is
    always redone from scratch and never resumed mid-file.

    Args:
        path: Destination file.
    """

    def __init__(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        self.path = path
        self._fh = path.open("w", encoding="utf-8")
        self.n = 0

    def write(self, row: dict[str, Any]) -> None:
        """Append one row and flush.

        Args:
            row: The row. Must contain no key named ``seconds`` -- a censored
                run must never share a key with a real measurement.

        Raises:
            ValueError: If the row carries a forbidden key.
        """
        if "seconds" in row:
            raise ValueError("'seconds' is a forbidden key: see PREREGISTRATION §5.3")
        self._fh.write(json.dumps(row, default=str) + "\n")
        self._fh.flush()
        self.n += 1

    def close(self) -> None:
        """Close the file."""
        self._fh.close()

    def __enter__(self) -> JsonlWriter:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()


def sha_file(path: Path) -> str:
    """SHA-256 of a file's bytes.

    Args:
        path: The file.

    Returns:
        The hex digest.
    """
    return hashlib.sha256(path.read_bytes()).hexdigest()


# --- the three shard kinds ----------------------------------------------------


def _instance_base(spec: dict[str, Any], row: dict[str, Any]) -> dict[str, Any]:
    """The columns every instance row carries, whatever its status.

    Args:
        spec: The shard specification.
        row: A frame row.

    Returns:
        The shared columns.
    """
    return {
        "shard_id": shard_id(spec),
        "arm": spec["kind"] if spec["kind"] != "xarm" else f"xarm_{spec['xarm']}",
        "gate": rs.GATE_NAME,
        "network": row["network"],
        "x": row["x"],
        "y": row["y"],
        "coverage": row["coverage"] if spec["kind"] == "flip" else None,
        "base_wrongness": spec.get("base_wrongness"),
        "bw_abs_level": spec.get("bw_abs"),
        "analyst_replicate": spec.get("replicate"),
        "n_tiers": spec.get("n_tiers"),
        "frame_row_id": row["row_id"],
        "largest_component_size": row["largest_component_size"],
        "component_size": row["component_size_recomputed"],
        "separation": row["separation_recomputed"],
        "separation_status": row["separation_status_recomputed"],
        "radius_committed_bw0": row["radius_recomputed"],
    }


def run_flip_shard(
    spec: dict[str, Any], frame: list[dict[str, Any]], out_dir: Path, n_draws: int
) -> dict[str, Any]:
    """Run one flip-arm shard: one network, coverage, base wrongness and analyst.

    The corrupted state is drawn **once per grid point per repetition** and
    scored against every admissible pair of the network, because on this corpus
    the analyst's knowledge is a property of the network, not of the query.

    Args:
        spec: The shard specification.
        frame: The frozen frame.
        out_dir: The results subtree.
        n_draws: Draws per grid point.

    Returns:
        The completion marker's payload.
    """
    sid = shard_id(spec)
    net, cov, rep_i = spec["network"], spec["coverage"], spec["replicate"]
    bw, bw_abs = spec.get("base_wrongness"), spec.get("bw_abs")
    rows = [r for r in frame if r["network"] == net and r["coverage"] == cov]
    dag, cpdag = rs.load_network(net)

    k_true = rs.select_knowledge(dag, cpdag, cov)
    family = f"{net}|cov{cov}|bw{bw}|bwa{bw_abs}|r{rep_i}"
    if bw_abs is not None and len(k_true) > 0:
        # Absolute level: exactly `bw_abs` claims reversed, expressed as a rate
        # so the corruption operator is unchanged.
        rng = np.random.default_rng(rs.derived_seed(net, cov, "base_abs", bw_abs, rep_i))
        k_b = tuple(_flip(list(k_true), rng, min(1.0, bw_abs / len(k_true))))
    elif bw:
        rng = np.random.default_rng(rs.derived_seed(net, cov, "base", bw, rep_i))
        k_b = tuple(_flip(list(k_true), rng, bw))
    else:
        k_b = tuple(sorted(k_true))

    g0, g0_reason = rs.build_g0(cpdag, k_b)
    shared = {
        "n_k": len(k_b),
        "k_b_sha256": rs.sha_of(sorted(k_b)),
        "g0_status": g0_reason,
        "g0_sha256": rs.sha_of(sorted(g0.directed_edges)) if g0 is not None else None,
        "g0_undirected_edges": len(g0.undirected_edges) if g0 is not None else None,
        "k_g0": rs.commitment_size(g0, cpdag) if g0 is not None else None,
        "shd_truth": rs.directed_symdiff(g0, dag) if g0 is not None else None,
        "n_claims_actually_wrong": (
            sum(1 for e in k_b if e not in set(k_true)) if k_b else 0
        ),
        "bw_abs": bw_abs,
        "bw_is_inert": bool(
            (bw or bw_abs) and len(k_b) > 0 and all(e in set(k_true) for e in k_b)
        ),
        "depth_grid": "frac10",
    }

    inst_path = out_dir / "shards" / f"{sid}.instances.jsonl"
    cell_path = out_dir / "shards" / f"{sid}.cells.jsonl"
    t0 = time.perf_counter()
    scored: list[tuple[dict[str, Any], frozenset]] = []

    with JsonlWriter(inst_path) as iw:
        for row in rows:
            irow = {**_instance_base(spec, row), **shared}
            if g0 is None:
                irow.update({"status": g0_reason, "z_star": None, "z_size": None,
                             "radius": None, "r_status": None, "dispatch_leg": None,
                             "oracle": None, "exact": None, "assumes": None,
                             "wall_until_timeout_s": None, "r_search_seconds": None,
                             "r_total_seconds": None})
                iw.write(irow)
                continue
            z_star, z_reason = rs.commit_z_star(g0, row["x"], row["y"])
            if z_star is None:
                irow.update({"status": z_reason, "z_star": None, "z_size": None,
                             "radius": None, "r_status": None, "dispatch_leg": None,
                             "oracle": None, "exact": None, "assumes": None,
                             "wall_until_timeout_s": None, "r_search_seconds": None,
                             "r_total_seconds": None})
                iw.write(irow)
                continue
            irow.update({"status": "ok", "z_star": sorted(z_star), "z_size": len(z_star)})
            irow.update(rs.radius_columns(cpdag, g0, row["x"], row["y"], z_star))
            iw.write(irow)
            scored.append((irow, z_star))

    n_cells = 0
    with JsonlWriter(cell_path) as cw:
        if g0 is not None and scored:
            cache = rs.ClosureCache(cpdag)
            for d in rs.depth_grid(len(k_b)):
                n_contra = 0
                survived = [0] * len(scored)
                symdiffs: list[int] = []
                for rep in range(n_draws):
                    k_cor, seed = rs.flip_draw(k_b, d, rep, family=family)
                    key, g = cache.graph(k_cor)
                    if g is None:
                        n_contra += 1
                        continue
                    symdiffs.append(rs.directed_symdiff(g0, g))
                    for i, (irow, z_star) in enumerate(scored):
                        if cache.survived(key, g, irow["x"], irow["y"], z_star):
                            survived[i] += 1
                for i, (irow, _z) in enumerate(scored):
                    crow = {
                        "shard_id": sid, "arm": "flip", "gate": rs.GATE_NAME,
                        "network": net, "x": irow["x"], "y": irow["y"],
                        "coverage": cov, "base_wrongness": bw, "bw_abs": bw_abs,
                        "analyst_replicate": rep_i,
                        "n_tiers": None, "frame_row_id": irow["frame_row_id"],
                        "grid_kind": "d_claims", "grid_point": d,
                        "frac_of_n_k": d / len(k_b),
                        **rs.cell_statistics(n_draws, n_contra, survived[i]),
                        "symdiff_proxy_not_distance_median": (
                            statistics.median(symdiffs) if symdiffs else None
                        ),
                        "status": "ok",
                    }
                    cw.write(crow)
                    n_cells += 1
            hits = {"graph_hits": cache.graph_hits, "graph_misses": cache.graph_misses,
                    "verdict_hits": cache.verdict_hits, "verdict_misses": cache.verdict_misses}
        else:
            hits = {"graph_hits": 0, "graph_misses": 0, "verdict_hits": 0, "verdict_misses": 0}

    return {
        "shard_id": sid, "kind": "flip", "spec": spec, "n_draws": n_draws,
        "n_frame_rows": len(rows), "n_scored": len(scored), "n_cells": n_cells,
        "g0_status": g0_reason, "n_k": len(k_b),
        "n_claims_actually_wrong": shared["n_claims_actually_wrong"],
        "bw_is_inert": shared["bw_is_inert"], "cache": hits,
        "elapsed_s": round(time.perf_counter() - t0, 3),
    }


def run_tiered_shard(
    spec: dict[str, Any], frame: list[dict[str, Any]], out_dir: Path, n_draws: int
) -> dict[str, Any]:
    """Run one tiered-arm shard: one network and tier count.

    The tiered arm is a **separate instance population** from the flip arm: its
    knowledge is generated from the network's temporal order rather than drawn
    at a coverage level, so it has no coverage parameter and its denominator is
    the network's distinct ``(X, Y)`` pairs. The axis is ``corruption_rate``,
    never the tiered ``d_claims`` column, which misses removals and whose zero
    bin is contaminated (synthetic PREREGISTRATION Appendix E).

    Args:
        spec: The shard specification.
        frame: The frozen frame.
        out_dir: The results subtree.
        n_draws: Draws per grid point.

    Returns:
        The completion marker's payload.
    """
    sid = shard_id(spec)
    net, nt = spec["network"], spec["n_tiers"]
    seen: set[tuple[str, str]] = set()
    rows = []
    for r in frame:
        if r["network"] != net:
            continue
        key = (r["x"], r["y"])
        if key in seen:
            continue
        seen.add(key)
        rows.append(r)
    dag, cpdag = rs.load_network(net)

    family = f"{net}|nt{nt}"
    k_ref, _seed = rs.tiered_draw(dag, cpdag, nt, 0.0, 0, family=f"{family}|ref")
    k_ref = tuple(sorted(k_ref))
    g0, g0_reason = rs.build_g0(cpdag, k_ref)
    shared = {
        "n_k": len(k_ref),
        "k_b_sha256": rs.sha_of(sorted(k_ref)),
        "g0_status": g0_reason,
        "g0_sha256": rs.sha_of(sorted(g0.directed_edges)) if g0 is not None else None,
        "g0_undirected_edges": len(g0.undirected_edges) if g0 is not None else None,
        "k_g0": rs.commitment_size(g0, cpdag) if g0 is not None else None,
        "shd_truth": rs.directed_symdiff(g0, dag) if g0 is not None else None,
        "n_claims_actually_wrong": 0,
        "bw_abs": None, "bw_is_inert": False,
        "depth_grid": "corruption_rate11",
    }

    inst_path = out_dir / "shards" / f"{sid}.instances.jsonl"
    cell_path = out_dir / "shards" / f"{sid}.cells.jsonl"
    t0 = time.perf_counter()
    scored: list[tuple[dict[str, Any], frozenset]] = []

    with JsonlWriter(inst_path) as iw:
        for row in rows:
            irow = {**_instance_base(spec, row), **shared}
            irow["coverage"] = None
            if g0 is None:
                irow.update({"status": g0_reason, "z_star": None, "z_size": None,
                             "radius": None, "r_status": None, "dispatch_leg": None,
                             "oracle": None, "exact": None, "assumes": None,
                             "wall_until_timeout_s": None, "r_search_seconds": None,
                             "r_total_seconds": None})
                iw.write(irow)
                continue
            z_star, z_reason = rs.commit_z_star(g0, row["x"], row["y"])
            if z_star is None:
                irow.update({"status": z_reason, "z_star": None, "z_size": None,
                             "radius": None, "r_status": None, "dispatch_leg": None,
                             "oracle": None, "exact": None, "assumes": None,
                             "wall_until_timeout_s": None, "r_search_seconds": None,
                             "r_total_seconds": None})
                iw.write(irow)
                continue
            irow.update({"status": "ok", "z_star": sorted(z_star), "z_size": len(z_star)})
            irow.update(rs.radius_columns(cpdag, g0, row["x"], row["y"], z_star))
            iw.write(irow)
            scored.append((irow, z_star))

    n_cells = 0
    with JsonlWriter(cell_path) as cw:
        if g0 is not None and scored:
            cache = rs.ClosureCache(cpdag)
            for rate in rs.RATE_GRID:
                n_contra = 0
                survived = [0] * len(scored)
                symdiffs: list[int] = []
                d_claims: list[int] = []
                for rep in range(n_draws):
                    k_cor, seed = rs.tiered_draw(dag, cpdag, nt, rate, rep, family=family)
                    key, g = cache.graph(k_cor)
                    d_claims.append(len(set(key) - set(k_ref)))
                    if g is None:
                        n_contra += 1
                        continue
                    symdiffs.append(rs.directed_symdiff(g0, g))
                    for i, (irow, z_star) in enumerate(scored):
                        if cache.survived(key, g, irow["x"], irow["y"], z_star):
                            survived[i] += 1
                for i, (irow, _z) in enumerate(scored):
                    crow = {
                        "shard_id": sid, "arm": "tiered", "gate": rs.GATE_NAME,
                        "network": net, "x": irow["x"], "y": irow["y"],
                        "coverage": None, "base_wrongness": None, "bw_abs": None,
                        "analyst_replicate": None, "n_tiers": nt,
                        "frame_row_id": irow["frame_row_id"],
                        "grid_kind": "corruption_rate", "grid_point": rate,
                        "frac_of_n_k": None,
                        **rs.cell_statistics(n_draws, n_contra, survived[i]),
                        "symdiff_proxy_not_distance_median": (
                            statistics.median(symdiffs) if symdiffs else None
                        ),
                        "d_claims_defective_do_not_bin": statistics.median(d_claims),
                        "status": "ok",
                    }
                    cw.write(crow)
                    n_cells += 1
            hits = {"graph_hits": cache.graph_hits, "graph_misses": cache.graph_misses,
                    "verdict_hits": cache.verdict_hits, "verdict_misses": cache.verdict_misses}
        else:
            hits = {"graph_hits": 0, "graph_misses": 0, "verdict_hits": 0, "verdict_misses": 0}

    return {
        "shard_id": sid, "kind": "tiered", "spec": spec, "n_draws": n_draws,
        "n_frame_rows": len(rows), "n_scored": len(scored), "n_cells": n_cells,
        "g0_status": g0_reason, "n_k": len(k_ref), "cache": hits,
        "elapsed_s": round(time.perf_counter() - t0, 3),
    }


def run_xarm_shard(
    spec: dict[str, Any], frame: list[dict[str, Any]], out_dir: Path, n_draws: int
) -> dict[str, Any]:
    """Run one cross-arm shard: one network, tier count and corruption process.

    The instance is the **tiered** instance and is constructed identically in
    both of a network's two cross-arm shards. Each shard stamps
    ``instance_sha256`` over ``(K_ref, dir(G0), Z*)``; the pairing step refuses
    to pair any instance whose two shards disagree. That is a stronger guarantee
    than shared process memory and it is checked rather than assumed.

    Both processes are scored on the shared intensity axis
    ``|dir(G0) Δ dir(G)| / |dir(G0)|``.

    Args:
        spec: The shard specification.
        frame: The frozen frame.
        out_dir: The results subtree.
        n_draws: Draws per grid point.

    Returns:
        The completion marker's payload.
    """
    sid = shard_id(spec)
    net, nt, arm = spec["network"], spec["n_tiers"], spec["xarm"]
    seen: set[tuple[str, str]] = set()
    rows = []
    for r in frame:
        if r["network"] != net:
            continue
        key = (r["x"], r["y"])
        if key in seen:
            continue
        seen.add(key)
        rows.append(r)
    dag, cpdag = rs.load_network(net)

    family = f"{net}|nt{nt}"
    k_ref, _s = rs.tiered_draw(dag, cpdag, nt, 0.0, 0, family=f"{family}|ref")
    k_ref = tuple(sorted(k_ref))
    g0, g0_reason = rs.build_g0(cpdag, k_ref)
    n_dir_g0 = len(g0.directed_edges) if g0 is not None else 0
    # `dir_g0_empty` is a structural fact about a graph that EXISTS and has no
    # directed edge, so the shared-intensity denominator would be zero. It must
    # never stand in for `g0 is None`, which is a different fact -- and, when
    # the reason is `o_g0_extensions_intractable`, a *measurement limit* rather
    # than a structural one. Conflating the two is the bug fixed on 2026-09-16
    # and recorded in PREREGISTRATION.md Appendix B.
    if g0 is None:
        g0_status = g0_reason
    elif n_dir_g0 == 0:
        g0_status = "dir_g0_empty"
    else:
        g0_status = g0_reason
    shared = {
        "n_k": len(k_ref),
        "k_b_sha256": rs.sha_of(sorted(k_ref)),
        "g0_status": g0_status,
        "g0_sha256": rs.sha_of(sorted(g0.directed_edges)) if g0 is not None else None,
        "g0_undirected_edges": len(g0.undirected_edges) if g0 is not None else None,
        "k_g0": rs.commitment_size(g0, cpdag) if g0 is not None else None,
        "shd_truth": rs.directed_symdiff(g0, dag) if g0 is not None else None,
        "n_claims_actually_wrong": 0,
        "bw_abs": None, "bw_is_inert": False,
        "n_dir_g0": n_dir_g0,
        "depth_grid": "xarm",
    }
    blocked = g0 is None or n_dir_g0 == 0

    inst_path = out_dir / "shards" / f"{sid}.instances.jsonl"
    cell_path = out_dir / "shards" / f"{sid}.cells.jsonl"
    t0 = time.perf_counter()
    scored: list[tuple[dict[str, Any], frozenset]] = []

    with JsonlWriter(inst_path) as iw:
        for row in rows:
            irow = {**_instance_base(spec, row), **shared}
            irow["coverage"] = None
            if blocked:
                irow.update({"status": shared["g0_status"], "z_star": None, "z_size": None,
                             "instance_sha256": None, "radius": None, "r_status": None,
                             "dispatch_leg": None, "oracle": None, "exact": None,
                             "assumes": None, "wall_until_timeout_s": None,
                             "r_search_seconds": None, "r_total_seconds": None})
                iw.write(irow)
                continue
            z_star, z_reason = rs.commit_z_star(g0, row["x"], row["y"])
            if z_star is None:
                irow.update({"status": z_reason, "z_star": None, "z_size": None,
                             "instance_sha256": None, "radius": None, "r_status": None,
                             "dispatch_leg": None, "oracle": None, "exact": None,
                             "assumes": None, "wall_until_timeout_s": None,
                             "r_search_seconds": None, "r_total_seconds": None})
                iw.write(irow)
                continue
            irow.update({
                "status": "ok", "z_star": sorted(z_star), "z_size": len(z_star),
                "instance_sha256": rs.sha_of(
                    (sorted(k_ref), sorted(g0.directed_edges), sorted(z_star))
                ),
            })
            irow.update(rs.radius_columns(cpdag, g0, row["x"], row["y"], z_star))
            iw.write(irow)
            scored.append((irow, z_star))

    grid: tuple[float, ...] | tuple[int, ...]
    grid = rs.XARM_TIERED_RATE_GRID if arm == "tiered" else rs.depth_grid(len(k_ref)) if k_ref else ()

    n_cells = 0
    with JsonlWriter(cell_path) as cw:
        if not blocked and scored and grid:
            cache = rs.ClosureCache(cpdag)
            for gp in grid:
                n_contra = 0
                survived = [0] * len(scored)
                intens: list[float] = []
                for rep in range(n_draws):
                    if arm == "tiered":
                        k_cor, seed = rs.tiered_draw(dag, cpdag, nt, gp, rep, family=f"{family}|xarm")
                    else:
                        k_cor, seed = rs.flip_draw(k_ref, int(gp), rep, family=f"{family}|xarm")
                    key, g = cache.graph(k_cor)
                    if g is None:
                        n_contra += 1
                        continue
                    intens.append(rs.directed_symdiff(g0, g) / n_dir_g0)
                    for i, (irow, z_star) in enumerate(scored):
                        if cache.survived(key, g, irow["x"], irow["y"], z_star):
                            survived[i] += 1
                med = statistics.median(intens) if intens else None
                for i, (irow, _z) in enumerate(scored):
                    crow = {
                        "shard_id": sid, "arm": f"xarm_{arm}", "gate": rs.GATE_NAME,
                        "network": net, "x": irow["x"], "y": irow["y"],
                        "coverage": None, "base_wrongness": None, "bw_abs": None,
                        "analyst_replicate": None, "n_tiers": nt,
                        "frame_row_id": irow["frame_row_id"],
                        "instance_sha256": irow["instance_sha256"],
                        "grid_kind": "corruption_rate" if arm == "tiered" else "d_claims",
                        "grid_point": gp,
                        **rs.cell_statistics(n_draws, n_contra, survived[i]),
                        "median_intensity": med,
                        "intensity_bin": rs.intensity_bin(med) if med is not None else None,
                        "n_dir_g0": n_dir_g0,
                        "status": "ok" if med is not None else "unresolved_all_contradictory",
                    }
                    cw.write(crow)
                    n_cells += 1
            hits = {"graph_hits": cache.graph_hits, "graph_misses": cache.graph_misses,
                    "verdict_hits": cache.verdict_hits, "verdict_misses": cache.verdict_misses}
        else:
            hits = {"graph_hits": 0, "graph_misses": 0, "verdict_hits": 0, "verdict_misses": 0}

    return {
        "shard_id": sid, "kind": "xarm", "spec": spec, "n_draws": n_draws,
        "n_frame_rows": len(rows), "n_scored": len(scored), "n_cells": n_cells,
        "g0_status": shared["g0_status"], "n_k": len(k_ref), "cache": hits,
        "elapsed_s": round(time.perf_counter() - t0, 3),
    }


RUNNERS = {"flip": run_flip_shard, "tiered": run_tiered_shard, "xarm": run_xarm_shard}


#: Rough per-draw Meek-closure cost in seconds, measured in the orientation
#: pass. Used **only** to order the work queue longest-first so the tail of a
#: parallel run is not one enormous shard; it never influences a measurement.
CLOSURE_COST_S: dict[str, float] = {
    "pathfinder": 1.0, "diabetes": 0.042, "arth150": 0.038, "ecoli70": 0.006,
    "hailfinder": 0.0055, "munin1": 0.0034, "insurance": 0.0022,
    "win95pts": 0.0021, "magic-irri": 0.002, "hepar2": 0.0016,
    "barley": 0.0011, "magic-niab": 0.0011, "sachs": 0.0008,
    "Polzer_2012": 0.0007, "paths": 0.0007, "child": 0.0006,
    "Sebastiani_2005": 0.0005, "Kampen_2014": 0.0004, "water": 0.0004,
    "magic": 0.0004, "barley_": 0.0004,
}
DEFAULT_CLOSURE_COST_S = 0.001


def estimated_cost_s(spec: dict[str, Any], frame: list[dict[str, Any]], n_draws: int) -> float:
    """A rough wall-time estimate for a shard, used only to order the queue.

    Args:
        spec: The shard specification.
        frame: The frozen frame.
        n_draws: Draws per grid point.

    Returns:
        Seconds. Never used in any measurement or analysis.
    """
    net = spec["network"]
    closure = CLOSURE_COST_S.get(net, DEFAULT_CLOSURE_COST_S)
    if spec["kind"] == "flip":
        rows = [r for r in frame if r["network"] == net and r["coverage"] == spec["coverage"]]
        n_grid = min(10, max(1, rows[0]["n_k"] if rows else 1))
    else:
        rows = list({(r["x"], r["y"]): r for r in frame if r["network"] == net}.values())
        n_grid = 11 if spec["kind"] == "tiered" else 10
    return n_grid * n_draws * (closure + 0.0004 * len(rows))


def run_one(spec: dict[str, Any], frame: list[dict[str, Any]], out_dir: Path, n_draws: int) -> dict[str, Any]:
    """Run one shard and write its completion marker.

    The marker is written last, and only on success, so a killed shard leaves no
    marker and is redone from scratch next time.

    Args:
        spec: The shard specification.
        frame: The frozen frame.
        out_dir: The results subtree.
        n_draws: Draws per grid point.

    Returns:
        The marker payload.
    """
    sid = shard_id(spec)
    payload = RUNNERS[spec["kind"]](spec, frame, out_dir, n_draws)
    inst = out_dir / "shards" / f"{sid}.instances.jsonl"
    cell = out_dir / "shards" / f"{sid}.cells.jsonl"
    payload.update({
        "instances_path": str(inst), "cells_path": str(cell),
        "instances_sha256": sha_file(inst), "cells_sha256": sha_file(cell),
        "git_sha": git_sha(short=False), "git_dirty": is_dirty(),
        "completed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    })
    marker = out_dir / "_done" / f"{sid}.json"
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text(json.dumps(payload, indent=2, default=str))
    return payload


# --- CLI ----------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    """The command-line interface.

    Returns:
        The parser.
    """
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=["list", "run", "status", "pool"])
    p.add_argument("--out-dir", default=str(DEFAULT_OUT))
    p.add_argument("--shard", default=None, help="shard id, for `run`")
    p.add_argument("--kind", default=None, choices=["flip", "tiered", "xarm"])
    p.add_argument("--n-draws", type=int, default=rs.N_DRAWS_DEFAULT)
    p.add_argument("--force", action="store_true", help="rerun a completed shard")
    p.add_argument("--workers", type=int, default=8, help="parallel shards, for `pool`")
    return p


def run_pool(specs: list[dict[str, Any]], frame: list[dict[str, Any]], out_dir: Path,
             n_draws: int, workers: int) -> int:
    """Run every pending shard as its own subprocess, longest-first.

    Each shard runs in a fresh interpreter, so a shard that crashes or is killed
    takes nothing else down and leaves no completion marker -- it is simply
    redone on the next invocation. The queue is ordered by estimated cost so the
    tail of the run is not one enormous shard.

    Args:
        specs: Candidate shards.
        frame: The frozen frame.
        out_dir: The results subtree.
        n_draws: Draws per grid point.
        workers: Maximum concurrent subprocesses.

    Returns:
        The number of shards that failed.
    """
    import subprocess
    import sys

    pending = [s for s in specs if not is_complete(out_dir, shard_id(s))]
    pending.sort(key=lambda s: -estimated_cost_s(s, frame, n_draws))
    queue = [shard_id(s) for s in pending]
    print(f"pool: {len(queue)} pending, {workers} workers", flush=True)

    running: dict[str, Any] = {}
    failed: list[str] = []
    logdir = out_dir / "_logs"
    logdir.mkdir(parents=True, exist_ok=True)
    started = time.time()

    while queue or running:
        while queue and len(running) < workers:
            sid = queue.pop(0)
            log = (logdir / f"{sid}.log").open("w")
            proc = subprocess.Popen(
                [sys.executable, "-m", "bkrobust.robustness.run_real_survival", "run",
                 "--shard", sid, "--out-dir", str(out_dir), "--n-draws", str(n_draws)],
                stdout=log, stderr=subprocess.STDOUT,
            )
            running[sid] = (proc, log)
            print(f"[{time.time() - started:8.0f}s] start {sid}", flush=True)
        time.sleep(1.0)
        for sid in list(running):
            proc, log = running[sid]
            if proc.poll() is None:
                continue
            log.close()
            del running[sid]
            ok = proc.returncode == 0 and is_complete(out_dir, sid)
            if not ok:
                failed.append(sid)
            print(f"[{time.time() - started:8.0f}s] {'done ' if ok else 'FAIL '} {sid}"
                  f" (rc={proc.returncode}, {len(queue)} queued)", flush=True)

    print(json.dumps({"failed": failed, "n_failed": len(failed),
                      "elapsed_s": round(time.time() - started, 1)}), flush=True)
    return len(failed)


def main(argv: list[str] | None = None) -> int:
    """Entry point.

    Args:
        argv: Argument vector, or ``None`` for ``sys.argv``.

    Returns:
        Process exit status.
    """
    args = build_arg_parser().parse_args(argv)
    out_dir = Path(args.out_dir)
    frame = load_frame(out_dir)
    specs = enumerate_shards(frame)
    if args.kind:
        specs = [s for s in specs if s["kind"] == args.kind]

    if args.command == "list":
        for s in specs:
            sid = shard_id(s)
            if not args.force and is_complete(out_dir, sid):
                continue
            print(sid)
        return 0

    if args.command == "pool":
        return 1 if run_pool(specs, frame, out_dir, args.n_draws, args.workers) else 0

    if args.command == "status":
        done = sum(1 for s in specs if is_complete(out_dir, shard_id(s)))
        print(json.dumps({"total": len(specs), "done": done, "pending": len(specs) - done}))
        return 0

    if not args.shard:
        raise SystemExit("run requires --shard")
    by_id = {shard_id(s): s for s in specs}
    if args.shard not in by_id:
        raise SystemExit(f"unknown shard {args.shard}")
    if is_complete(out_dir, args.shard) and not args.force:
        print(json.dumps({"shard_id": args.shard, "skipped": "already complete"}))
        return 0
    payload = run_one(by_id[args.shard], frame, out_dir, args.n_draws)
    print(json.dumps({k: payload[k] for k in
                      ("shard_id", "n_frame_rows", "n_scored", "n_cells", "g0_status", "elapsed_s")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
