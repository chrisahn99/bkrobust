"""Sharded, resumable driver for the elicited-knowledge survival sweep.

Runs the survival/ranking experiment with ``K`` read from
``results/elicit/knowledge.json`` instead of built by
``benchmarks.measure.select_knowledge``. The knowledge layer is
:mod:`bkrobust.robustness.llm_survival`; the grid, seeding, closure cache,
radius and cell machinery are imported unchanged from
:mod:`bkrobust.robustness.real_survival` so the elicited arm and the committed
real arm cannot disagree about what an endpoint means.

The shard is ``(condition, network)``
--------------------------------------
Elicited ``K`` is a property of the condition and the network, exactly as the
committed real arm's ``K`` is a property of the network and the coverage. It
carries no coverage parameter, no base-wrongness dial and no analyst
replicate: the model answered once, and what it answered is the analyst's
knowledge. A shard therefore draws its corruption once per grid point per
repetition and scores it against every distinct ``(X, Y)`` pair the frozen
frame admits for that network.

Incrementality contract, inherited from ``run_real_survival``
--------------------------------------------------------------
A shard owns exactly two append-only output files, flushed per row, and writes
``_done/<shard_id>.json`` only when it finishes. A shard with no marker is
never resumed -- its files are truncated and it is redone from scratch -- so
partial output is cleanly discardable by construction. Re-invoking any stage
skips completed shards.

Every marker carries the SHA-256 of ``knowledge.json``, so a results tree can
always be tied back to the elicitation bundle it was computed from, and of the
frozen frame, which is hash-checked before anything runs on it.

Degenerate shards run to completion and are written out
--------------------------------------------------------
A condition that declined every pair (``n_k_zero``), one whose claims admit no
MPDAG (``k_contradictory``) or one whose ``G0`` is too large to enumerate
extensions for (``o_g0_extensions_intractable``) still writes a full set of
instance rows carrying that status and an empty cell file. It is a measured
outcome of the panel, not an error, and it must appear in the denominator.
``pathfinder`` is the one exception: it is absent from ``knowledge.json``
altogether, so it has no shard at all and is named in the manifest's
``networks_absent_from_bundle``.

Usage::

    PYTHONPATH=src .venv/bin/python -m bkrobust.robustness.run_llm_survival list
    PYTHONPATH=src .venv/bin/python -m bkrobust.robustness.run_llm_survival run --shard <id>
    PYTHONPATH=src .venv/bin/python -m bkrobust.robustness.run_llm_survival pool --workers 8
    PYTHONPATH=src .venv/bin/python -m bkrobust.robustness.run_llm_survival status
    PYTHONPATH=src .venv/bin/python -m bkrobust.robustness.run_llm_survival manifest

No global RNG is touched anywhere.
"""

from __future__ import annotations

import argparse
import json
import statistics
import time
from pathlib import Path
from typing import Any

from bkrobust.core.resultsio import git_sha, is_dirty
from bkrobust.robustness import llm_survival as ls
from bkrobust.robustness import real_survival as rs
from bkrobust.robustness.run_real_survival import (
    JsonlWriter,
    _censored_cell,
    load_frame,
    sha_file,
)

DEFAULT_OUT = Path("results/axis_robustness_llm")

#: The frozen frame the instance population is read from. It is the committed
#: 831-row corpus of ``results/axis_robustness_real/``; this sweep reads it and
#: never writes to that tree.
DEFAULT_FRAME_DIR = Path("results/axis_robustness_real")

#: Per-shard wall cap. A safety net, not a scientific parameter: grid points
#: are swept in ascending order, so a cap that fires truncates the
#: high-intensity end of a curve and every truncated grid point is written as a
#: ``censored_wall_cap`` cell rather than dropped.
SHARD_WALL_CAP_S: float = 18000.0


# --- shards -------------------------------------------------------------------


def shard_id(spec: dict[str, Any]) -> str:
    """The deterministic id, and filename stem, of a shard.

    Args:
        spec: A shard specification.

    Returns:
        The id, ``llm__<condition>__<network>``.
    """
    return f"llm__{spec['condition']}__{spec['network']}"


def enumerate_shards(
    frame: list[dict[str, Any]], bundle: dict[str, Any], *, naming: str | None = None
) -> list[dict[str, Any]]:
    """Every unit of work the campaign consists of, in a deterministic order.

    The cross product of the panel's conditions with the frame's networks,
    minus the ``(condition, network)`` cells the bundle never covered.

    Args:
        frame: The frozen frame.
        bundle: The loaded elicitation bundle.
        naming: Restrict to one naming arm, or ``None`` for every condition in
            the bundle (the panel and its scrambled control together).

    Returns:
        Shard specifications, sorted by id.
    """
    networks = sorted({r["network"] for r in frame})
    specs = [
        {"condition": c, "network": n}
        for c in ls.panel_conditions(bundle, naming=naming)
        for n in networks
        if ls.has_network(bundle, c, n)
    ]
    return sorted(specs, key=shard_id)


def absent_networks(frame: list[dict[str, Any]], bundle: dict[str, Any]) -> list[str]:
    """Frame networks no condition in the bundle was ever asked about.

    Args:
        frame: The frozen frame.
        bundle: The loaded elicitation bundle.

    Returns:
        The network names, sorted. Named in the manifest so a network that
        drops out of the sweep never drops out silently.
    """
    networks = sorted({r["network"] for r in frame})
    return [n for n in networks if not any(ls.has_network(bundle, c, n) for c in bundle)]


def is_complete(out_dir: Path, sid: str) -> bool:
    """Whether a shard has a completion marker.

    Args:
        out_dir: The results subtree.
        sid: The shard id.

    Returns:
        ``True`` if and only if the marker exists.
    """
    return (out_dir / "_done" / f"{sid}.json").is_file()


# --- the shard ----------------------------------------------------------------


def _instance_base(spec: dict[str, Any], row: dict[str, Any]) -> dict[str, Any]:
    """The columns every instance row carries, whatever its status.

    Args:
        spec: The shard specification.
        row: A frame row.

    Returns:
        The shared identity columns.
    """
    return {
        "shard_id": shard_id(spec),
        "arm": "llm",
        "gate": rs.GATE_NAME,
        "network": row["network"],
        "x": row["x"],
        "y": row["y"],
        # No coverage parameter exists on this arm: `K` is what the model said,
        # not a fraction of the truth. Written as None rather than omitted so
        # the schema matches the committed real arm's rows column for column.
        "coverage": None,
        "base_wrongness": None,
        "bw_abs_level": None,
        "analyst_replicate": None,
        "n_tiers": None,
        "frame_row_id": row["row_id"],
        "largest_component_size": row["largest_component_size"],
        "component_size": row["component_size_recomputed"],
        "separation": row["separation_recomputed"],
        "separation_status": row["separation_status_recomputed"],
        # The committed radius of the SAME pair under programmatic K at the
        # frame's own coverage. Kept for reference only: it is not this arm's
        # radius, and the analysis must never substitute one for the other.
        "radius_committed_programmatic_k": row["radius_recomputed"],
    }


_BLOCKED_COLUMNS: dict[str, Any] = {
    "z_star": None, "z_size": None, "radius": None, "r_status": None,
    "dispatch_leg": None, "oracle": None, "exact": None, "assumes": None,
    "wall_until_timeout_s": None, "r_search_seconds": None, "r_total_seconds": None,
}


def run_shard(
    spec: dict[str, Any],
    frame: list[dict[str, Any]],
    bundle: dict[str, Any],
    out_dir: Path,
    n_draws: int,
    wall_cap_s: float = SHARD_WALL_CAP_S,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run one ``(condition, network)`` shard end to end.

    Args:
        spec: The shard specification.
        frame: The frozen frame.
        bundle: The loaded elicitation bundle.
        out_dir: The results subtree.
        n_draws: Draws per grid point.
        wall_cap_s: Per-shard wall cap.
        verbose: Print a progress line per grid point to stdout.

    Returns:
        The completion marker's payload.
    """
    sid = shard_id(spec)
    cond, net = spec["condition"], spec["network"]
    rows = ls.distinct_pairs(frame, net)
    dag, cpdag = rs.load_network(net)

    k, g0, g0_reason, shared = ls.build_state(bundle, cond, net, dag, cpdag)
    family = ls.panel_family(cond, net)

    if verbose:
        print(
            f"[{sid}] |K|={len(k)} asked={shared['n_asked']} "
            f"acc={shared['k_accuracy']} g0={g0_reason} pairs={len(rows)} "
            f"draws={n_draws}",
            flush=True,
        )
    if shared["n_off_skeleton"]:
        # Never silently tolerated: a claim outside the CPDAG skeleton would
        # mean the answers were parsed against a different graph than the one
        # being scored here.
        print(
            f"[{sid}] WARNING off-skeleton claims: {shared['off_skeleton_claims']}",
            flush=True,
        )

    inst_path = out_dir / "shards" / f"{sid}.instances.jsonl"
    cell_path = out_dir / "shards" / f"{sid}.cells.jsonl"
    t0 = time.perf_counter()
    scored: list[tuple[dict[str, Any], frozenset]] = []
    status_counts: dict[str, int] = {}

    with JsonlWriter(inst_path) as iw:
        for row in rows:
            irow = {**_instance_base(spec, row), **shared}
            if g0 is None:
                irow.update({"status": g0_reason, **_BLOCKED_COLUMNS})
                status_counts[g0_reason] = status_counts.get(g0_reason, 0) + 1
                iw.write(irow)
                continue
            z_star, z_reason = rs.commit_z_star(g0, row["x"], row["y"])
            if z_star is None:
                irow.update({"status": z_reason, **_BLOCKED_COLUMNS})
                status_counts[z_reason] = status_counts.get(z_reason, 0) + 1
                iw.write(irow)
                continue
            irow.update({"status": "ok", "z_star": sorted(z_star), "z_size": len(z_star)})
            irow.update(rs.radius_columns(cpdag, g0, row["x"], row["y"], z_star))
            status_counts["ok"] = status_counts.get("ok", 0) + 1
            iw.write(irow)
            scored.append((irow, z_star))

    if verbose:
        print(
            f"[{sid}] instances written: {len(rows)} ({status_counts}) "
            f"scored={len(scored)} in {time.perf_counter() - t0:.1f}s",
            flush=True,
        )

    n_cells = 0
    n_censored = 0
    last_cells: list[dict[str, Any]] = []
    with JsonlWriter(cell_path) as cw:
        if g0 is not None and scored:
            cache = rs.ClosureCache(cpdag)
            grid_all = rs.depth_grid(len(k))
            for gi, d in enumerate(grid_all):
                if time.perf_counter() - t0 > wall_cap_s:
                    print(f"[{sid}] WALL CAP at depth {d}; censoring the tail", flush=True)
                    for gp in grid_all[gi:]:
                        for tmpl in last_cells:
                            cw.write(_censored_cell(tmpl, gp, time.perf_counter() - t0))
                            n_cells += 1
                            n_censored += 1
                    break
                t_d = time.perf_counter()
                n_contra = 0
                survived = [0] * len(scored)
                symdiffs: list[int] = []
                for rep in range(n_draws):
                    k_cor, _seed = rs.flip_draw(k, d, rep, family=family)
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
                        "shard_id": sid, "arm": "llm", "gate": rs.GATE_NAME,
                        "condition": cond, "model": shared["model"],
                        "family": shared["family"], "naming": shared["naming"],
                        "network": net, "x": irow["x"], "y": irow["y"],
                        "coverage": None, "base_wrongness": None, "bw_abs": None,
                        "analyst_replicate": None, "n_tiers": None,
                        "frame_row_id": irow["frame_row_id"],
                        "grid_kind": "d_claims", "grid_point": d,
                        "frac_of_n_k": d / len(k),
                        "n_k": len(k), "k_accuracy": shared["k_accuracy"],
                        **rs.cell_statistics(n_draws, n_contra, survived[i]),
                        "symdiff_proxy_not_distance_median": (
                            statistics.median(symdiffs) if symdiffs else None
                        ),
                        "status": "ok",
                    }
                    cw.write(crow)
                    n_cells += 1
                    last_cells = [
                        c for c in last_cells if c["frame_row_id"] != crow["frame_row_id"]
                    ]
                    last_cells.append(crow)
                if verbose:
                    n_eval = n_draws - n_contra
                    s_vals = [s / n_eval for s in survived] if n_eval else []
                    s_txt = (
                        f"S_median={statistics.median(s_vals):.3f}" if s_vals
                        else "S=undefined(all_contradictory)"
                    )
                    print(
                        f"[{sid}]   d={d}/{grid_all[-1]} contra={n_contra}/{n_draws} "
                        f"{s_txt} ({time.perf_counter() - t_d:.1f}s)",
                        flush=True,
                    )
            hits = {"graph_hits": cache.graph_hits, "graph_misses": cache.graph_misses,
                    "verdict_hits": cache.verdict_hits, "verdict_misses": cache.verdict_misses}
        else:
            hits = {"graph_hits": 0, "graph_misses": 0, "verdict_hits": 0, "verdict_misses": 0}

    radii = [i["radius"] for i, _ in scored if i.get("r_status") == "ok"]
    return {
        "shard_id": sid, "kind": "llm", "spec": spec, "n_draws": n_draws,
        "condition": cond, "network": net,
        "model": shared["model"], "family": shared["family"], "naming": shared["naming"],
        "n_frame_rows": len(rows), "n_scored": len(scored), "n_cells": n_cells,
        "n_cells_censored_wall_cap": n_censored, "wall_cap_s": wall_cap_s,
        "g0_status": g0_reason, "n_k": len(k),
        "n_asked": shared["n_asked"], "assert_rate": shared["assert_rate"],
        "k_accuracy": shared["k_accuracy"], "n_k_correct": shared["n_k_correct"],
        "n_k_wrong": shared["n_k_wrong"], "n_off_skeleton": shared["n_off_skeleton"],
        "k_g0": shared["k_g0"], "shd_truth": shared["shd_truth"],
        "instance_status_counts": status_counts,
        "radius_median_ok": statistics.median(radii) if radii else None,
        "n_radius_ok": len(radii),
        "cache": hits,
        "elapsed_s": round(time.perf_counter() - t0, 3),
    }


def run_one(
    spec: dict[str, Any], frame: list[dict[str, Any]], bundle: dict[str, Any],
    bundle_sha: str, frame_sha: str, out_dir: Path, n_draws: int,
    wall_cap_s: float = SHARD_WALL_CAP_S,
) -> dict[str, Any]:
    """Run one shard and write its completion marker.

    The marker is written last and only on success, so a killed shard leaves no
    marker and is redone from scratch next time.

    Args:
        spec: The shard specification.
        frame: The frozen frame.
        bundle: The loaded elicitation bundle.
        bundle_sha: SHA-256 of ``knowledge.json``.
        frame_sha: SHA-256 of ``frame.jsonl``.
        out_dir: The results subtree.
        n_draws: Draws per grid point.
        wall_cap_s: Per-shard wall cap.

    Returns:
        The marker payload.
    """
    sid = shard_id(spec)
    payload = run_shard(spec, frame, bundle, out_dir, n_draws, wall_cap_s)
    inst = out_dir / "shards" / f"{sid}.instances.jsonl"
    cell = out_dir / "shards" / f"{sid}.cells.jsonl"
    payload.update({
        "instances_path": str(inst), "cells_path": str(cell),
        "instances_sha256": sha_file(inst), "cells_sha256": sha_file(cell),
        "knowledge_json_sha256": bundle_sha, "frame_sha256": frame_sha,
        "git_sha": git_sha(short=False), "git_dirty": is_dirty(),
        "completed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    })
    marker = out_dir / "_done" / f"{sid}.json"
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text(json.dumps(payload, indent=2, default=str))
    return payload


# --- cost model ---------------------------------------------------------------


def estimated_cost_s(spec: dict[str, Any], frame: list[dict[str, Any]],
                     bundle: dict[str, Any], n_draws: int) -> float:
    """A rough wall-time estimate for a shard, used only to order the queue.

    Args:
        spec: The shard specification.
        frame: The frozen frame.
        bundle: The loaded elicitation bundle.
        n_draws: Draws per grid point.

    Returns:
        Seconds. Never used in any measurement or analysis.
    """
    from bkrobust.robustness.run_real_survival import (
        CLOSURE_COST_S,
        DEFAULT_CLOSURE_COST_S,
    )

    net = spec["network"]
    closure = CLOSURE_COST_S.get(net, DEFAULT_CLOSURE_COST_S)
    n_pairs = len(ls.distinct_pairs(frame, net))
    try:
        k, _ = ls.elicited_k(bundle, spec["condition"], net)
    except KeyError:
        return 0.0
    n_grid = min(10, max(1, len(k)))
    return n_grid * n_draws * (closure + 0.0004 * n_pairs)


# --- manifest -----------------------------------------------------------------


def build_manifest(
    frame: list[dict[str, Any]], bundle: dict[str, Any], bundle_sha: str,
    frame_sha: str, specs: list[dict[str, Any]], out_dir: Path, n_draws: int,
) -> dict[str, Any]:
    """The campaign's manifest: what was enumerated, and what was left out.

    Args:
        frame: The frozen frame.
        bundle: The loaded elicitation bundle.
        bundle_sha: SHA-256 of ``knowledge.json``.
        frame_sha: SHA-256 of ``frame.jsonl``.
        specs: The enumerated shards.
        out_dir: The results subtree.
        n_draws: Draws per grid point.

    Returns:
        The manifest payload, also written to ``manifest.json``.
    """
    networks = sorted({r["network"] for r in frame})
    conds = ls.panel_conditions(bundle, naming=None)
    covered = {(s["condition"], s["network"]) for s in specs}
    payload = {
        "campaign": "axis_robustness_llm",
        "k_source": "results/elicit/knowledge.json",
        "k_source_note": (
            "benchmarks.measure.select_knowledge is never called on this arm; "
            "K is the model's own asserted orientations."
        ),
        "knowledge_json_sha256": bundle_sha,
        "frame": str(DEFAULT_FRAME_DIR / "frame.jsonl"),
        "frame_sha256": frame_sha,
        "n_frame_rows": len(frame),
        "n_frame_networks": len(networks),
        "n_distinct_pairs": sum(len(ls.distinct_pairs(frame, n)) for n in networks),
        "conditions": [ls.condition_meta(bundle, c) for c in conds],
        "n_conditions": len(conds),
        "panel_real_naming": ls.panel_conditions(bundle, naming=ls.REAL_NAMING),
        "control_scrambled_naming": ls.panel_conditions(bundle, naming=ls.SCRAMBLED_NAMING),
        "networks_absent_from_bundle": absent_networks(frame, bundle),
        "cells_not_elicited": sorted(
            [c, n] for c in conds for n in networks if (c, n) not in covered
        ),
        "n_shards": len(specs),
        "n_draws": n_draws,
        "wall_cap_s": SHARD_WALL_CAP_S,
        "grid": "rs.depth_grid(|K|) -- frac10, identical to the committed real arm",
        "git_sha": git_sha(short=False),
        "git_dirty": is_dirty(),
        "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "manifest.json").write_text(json.dumps(payload, indent=2, default=str))
    return payload


# --- CLI ----------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    """The command-line interface.

    Returns:
        The parser.
    """
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=["list", "run", "pool", "status", "manifest"])
    p.add_argument("--out-dir", default=str(DEFAULT_OUT))
    p.add_argument("--frame-dir", default=str(DEFAULT_FRAME_DIR))
    p.add_argument("--elicit", default=str(ls.DEFAULT_ELICIT_PATH))
    p.add_argument("--shard", default=None, help="shard id, for `run`")
    p.add_argument("--condition", default=None, help="restrict to one condition")
    p.add_argument("--network", default=None, help="restrict to one network")
    p.add_argument("--naming", default=None, choices=["real", "scrambled"])
    p.add_argument("--n-draws", type=int, default=rs.N_DRAWS_DEFAULT)
    p.add_argument("--force", action="store_true", help="rerun a completed shard")
    p.add_argument("--workers", type=int, default=8, help="parallel shards, for `pool`")
    p.add_argument("--wall-cap-s", type=float, default=SHARD_WALL_CAP_S)
    return p


def run_pool(specs: list[dict[str, Any]], frame: list[dict[str, Any]],
             bundle: dict[str, Any], out_dir: Path, elicit: Path,
             n_draws: int, workers: int, wall_cap_s: float) -> int:
    """Run every pending shard as its own subprocess, longest-first.

    Each shard runs in a fresh interpreter, so one that crashes or is killed
    takes nothing else down and leaves no completion marker -- it is simply
    redone on the next invocation. A failure is logged and the sweep continues:
    the panel is never halted by one degenerate cell.

    Args:
        specs: Candidate shards.
        frame: The frozen frame.
        bundle: The loaded elicitation bundle.
        out_dir: The results subtree.
        elicit: Path to the elicitation bundle, passed to each child.
        n_draws: Draws per grid point.
        workers: Maximum concurrent subprocesses.
        wall_cap_s: Per-shard wall cap.

    Returns:
        The number of shards that failed.
    """
    import subprocess
    import sys

    pending = [s for s in specs if not is_complete(out_dir, shard_id(s))]
    pending.sort(key=lambda s: -estimated_cost_s(s, frame, bundle, n_draws))
    queue = [shard_id(s) for s in pending]
    print(f"pool: {len(queue)} pending of {len(specs)}, {workers} workers", flush=True)

    running: dict[str, Any] = {}
    failed: list[str] = []
    logdir = out_dir / "_logs"
    logdir.mkdir(parents=True, exist_ok=True)
    started = time.time()
    n_done = 0

    while queue or running:
        while queue and len(running) < workers:
            sid = queue.pop(0)
            log = (logdir / f"{sid}.log").open("w")
            proc = subprocess.Popen(
                [sys.executable, "-m", "bkrobust.robustness.run_llm_survival", "run",
                 "--shard", sid, "--out-dir", str(out_dir), "--elicit", str(elicit),
                 "--n-draws", str(n_draws), "--wall-cap-s", str(wall_cap_s)],
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
            n_done += 1
            if not ok:
                failed.append(sid)
            print(f"[{time.time() - started:8.0f}s] {'done ' if ok else 'FAIL '} {sid}"
                  f" (rc={proc.returncode}, {n_done}/{len(pending)} finished, "
                  f"{len(queue)} queued)", flush=True)

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
    frame_dir = Path(args.frame_dir)
    elicit = Path(args.elicit)

    frame = load_frame(frame_dir)
    frame_sha = json.loads((frame_dir / "frame_hash.json").read_text())["frame_sha256"]
    bundle, bundle_sha = ls.load_bundle(elicit)
    specs = enumerate_shards(frame, bundle, naming=args.naming)
    if args.condition:
        specs = [s for s in specs if s["condition"] == args.condition]
    if args.network:
        specs = [s for s in specs if s["network"] == args.network]

    if args.command == "manifest":
        payload = build_manifest(
            frame, bundle, bundle_sha, frame_sha,
            enumerate_shards(frame, bundle, naming=None), out_dir, args.n_draws,
        )
        print(json.dumps({k: payload[k] for k in (
            "n_shards", "n_conditions", "n_frame_networks", "n_distinct_pairs",
            "networks_absent_from_bundle", "cells_not_elicited",
            "knowledge_json_sha256")}, indent=2))
        return 0

    if args.command == "list":
        for s in specs:
            sid = shard_id(s)
            if not args.force and is_complete(out_dir, sid):
                continue
            print(sid)
        return 0

    if args.command == "status":
        done = sum(1 for s in specs if is_complete(out_dir, shard_id(s)))
        print(json.dumps({"total": len(specs), "done": done, "pending": len(specs) - done}))
        return 0

    if args.command == "pool":
        build_manifest(frame, bundle, bundle_sha, frame_sha,
                       enumerate_shards(frame, bundle, naming=None), out_dir, args.n_draws)
        n_failed = run_pool(specs, frame, bundle, out_dir, elicit,
                            args.n_draws, args.workers, args.wall_cap_s)
        return 1 if n_failed else 0

    if not args.shard:
        raise SystemExit("run requires --shard")
    by_id = {shard_id(s): s for s in specs}
    if args.shard not in by_id:
        raise SystemExit(f"unknown shard {args.shard}")
    if is_complete(out_dir, args.shard) and not args.force:
        print(json.dumps({"shard_id": args.shard, "skipped": "already complete"}))
        return 0
    payload = run_one(by_id[args.shard], frame, bundle, bundle_sha, frame_sha,
                      out_dir, args.n_draws, args.wall_cap_s)
    print(json.dumps({k: payload[k] for k in (
        "shard_id", "n_k", "k_accuracy", "n_frame_rows", "n_scored", "n_cells",
        "g0_status", "elapsed_s")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
