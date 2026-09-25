"""Recompute every reported breakdown radius under GAC and compare to back-door.

Every radius the paper reports was produced by ``bkrobust.hybrid.breakdown_radius``
using Pearl's back-door predicate lifted to MPDAGs
(``bkrobust.mpdag_criterion.is_valid_mpdag``). This script independently
recomputes the same radius two ways -- via the same back-door predicate (a
sanity check on the reconstruction) and via the generalised adjustment
criterion (``bkrobust.gac.mpdag_level.is_gac_valid_mpdag``) -- using a plain
breadth-first retraction search from ``G0`` to full, unlimited depth. It does
NOT call ``bkrobust.hybrid.breakdown_radius`` or anything in
``bkrobust.sat``/``bkrobust.hybrid``, which another session is concurrently
rewriting.

Method (paper's Theorem: the nearest failure is reachable by retractions):
  1. Rebuild (cpdag, G0, X, Y, committed Z) exactly as the original experiment
     did, reusing bkrobust.robustness.real_survival / llm_survival's own
     builder functions.
  2. Check the predicate at the CPDAG itself (top of the order). Failure is
     upward-closed and the CPDAG is the maximum element, so radius = inf iff Z
     is valid there -- matches bkrobust.hybrid's own top-state shortcut, but
     implemented locally against each predicate directly.
  3. Otherwise, breadth-first search upward from G0 via
     bkrobust.search.exact_fast.radius_local_up_fast, with max_depth=None
     (unlimited), and the predicate as the ``fails`` oracle.

Usage::

    PYTHONPATH=src .venv/bin/python experiments/recompute_radius_gac.py <dataset> [--limit N] [--time-budget S]

Datasets: llm, llm_v2, p6, unconfounded, mean_table, e2e, demo, all
"""

from __future__ import annotations

import argparse
import csv
import json
import multiprocessing as mp
import signal
import sys
import time
from pathlib import Path

sys.path.insert(0, "src")


class _PerUnitTimeout(Exception):
    pass


def _alarm_handler(signum, frame):  # noqa: ARG001
    raise _PerUnitTimeout()


class per_unit_timeout:
    """Context manager: raise _PerUnitTimeout if the block runs past `seconds`."""

    def __init__(self, seconds: float):
        self.seconds = seconds

    def __enter__(self):
        if self.seconds and self.seconds > 0:
            signal.signal(signal.SIGALRM, _alarm_handler)
            signal.setitimer(signal.ITIMER_REAL, self.seconds)
        return self

    def __exit__(self, exc_type, exc, tb):
        signal.setitimer(signal.ITIMER_REAL, 0)
        return False

from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag  # noqa: E402
from bkrobust.mpdag_criterion import is_valid_mpdag  # noqa: E402
from bkrobust.robustness import llm_survival as ls  # noqa: E402
from bkrobust.robustness import real_survival as rs  # noqa: E402
from bkrobust.search.exact_fast import radius_local_up_fast  # noqa: E402

OUT_DIR = Path("results/radius_gac_recompute")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def radius_under(predicate_fails, cpdag: MPDAG, g0: MPDAG, *, max_depth=None):
    """radius under an arbitrary validity-failure predicate, exact BFS.

    ``max_depth=None`` means unlimited (full up-set), matching the paper's
    theorem (nearest failure reachable by retractions). A finite max_depth is
    only ever used as a defensive budget against pathological blow-up; when it
    is hit without a definitive answer, ``exact`` is False.

    Returns (radius_or_UNREACHED, exact_bool, method_str).
    """
    if not predicate_fails(cpdag):
        return UNREACHED, True, "top_state"
    res = radius_local_up_fast(cpdag, g0, predicate_fails, max_depth=max_depth)
    return res.radius, res.exact, "bfs"


def gac_radius(cpdag, g0, x, y, z, *, max_depth=None):
    def fails(g):
        return not is_gac_valid_mpdag(g, x, y, z)

    return radius_under(fails, cpdag, g0, max_depth=max_depth)


def bd_radius(cpdag, g0, x, y, z, *, max_depth=None):
    def fails(g):
        return not is_valid_mpdag(g, x, y, z)

    return radius_under(fails, cpdag, g0, max_depth=max_depth)


#: Wall-clock cap on a single predicate's BFS, as a defence against
#: pathological blow-up (a huge up-set the search must exhaust to prove
#: UNREACHED). Units that hit it are recorded as not-recomputed with a
#: reason, never silently reported as equal or different.
PER_UNIT_TIMEOUT_S = 8.0

#: Networks above this node count are skipped for the GAC recompute (not
#: attempted at all, rather than left to a per-unit timeout): on this
#: hardware a single MPDAG-validity BFS over e.g. the 413-node ``diabetes``
#: or the 186-node ``munin1`` network measured 20-30+ seconds *per query*,
#: which Python's SIGALRM cannot reliably interrupt mid-call (a single
#: large set/graph operation is one uninterruptible C-level step), and which
#: would consume the whole time budget on a handful of networks. Recorded as
#: not_recomputed / "network_too_large_skipped", never silently dropped.
NETWORK_NODE_LIMIT = 120


def gac_valid_at_g0(g0, x, y, z) -> bool:
    return is_gac_valid_mpdag(g0, x, y, z)


def process_unit(cpdag: MPDAG, g0: MPDAG, x: str, y: str, z_star, stored_radius_val,
                  *, do_sanity: bool, timeout_s: float = PER_UNIT_TIMEOUT_S):
    """Recompute GAC radius (and, if requested, a back-door sanity recheck)
    for one already-built unit. Returns a result dict; never raises.
    """
    out = {
        "z_star": json.dumps(sorted(z_star)),
        "stored_radius": stored_radius_val,
        "gac_radius": None, "gac_exact": False, "equal": None,
        "z_gac_valid_at_g0": None, "backdoor_recheck": None,
        "sanity_status": "",
    }
    try:
        with per_unit_timeout(timeout_s):
            gac_rad, gac_exact, _ = gac_radius(cpdag, g0, x, y, z_star)
            g0_gac_valid = gac_valid_at_g0(g0, x, y, z_star)
    except _PerUnitTimeout:
        out["sanity_status"] = "not_recomputed_timeout"
        return out

    out["gac_radius"] = gac_rad
    out["gac_exact"] = gac_exact
    out["z_gac_valid_at_g0"] = g0_gac_valid
    out["equal"] = (stored_radius_val == gac_rad) if gac_exact else None

    if do_sanity:
        try:
            with per_unit_timeout(timeout_s):
                bd_rad, bd_exact, _ = bd_radius(cpdag, g0, x, y, z_star)
            out["backdoor_recheck"] = bd_rad
            if not bd_exact:
                out["sanity_status"] = "not_recomputed_timeout"
            elif stored_radius_val is not None and bd_rad == stored_radius_val:
                out["sanity_status"] = "match"
            else:
                out["sanity_status"] = "MISMATCH"
        except _PerUnitTimeout:
            out["sanity_status"] = "not_recomputed_timeout"
    return out


def run_generic_dataset(units, out_csv: Path, *, time_budget_s: float, sanity_check_n: int,
                         id_fields: list[str], total_hint: int | None = None,
                         build_timeout_s: float = 15.0):
    """units: iterable of dicts with id_fields + cpdag,g0,x,y,z_star,stored_radius.

    ``units`` is consumed lazily (it may itself do expensive per-row work,
    e.g. rebuilding a synthetic instance from its seed) so the time budget and
    per-item timeout apply to construction as well as to the GAC/back-door
    recompute -- a single pathological instance cannot silently consume the
    whole run before the loop body ever executes.

    Handles timing, aggregation and CSV writing uniformly for the datasets
    that hand over already-built (cpdag, g0, z_star) tuples (p6, unconfounded,
    mean_table, e2e, demo), as opposed to the LLM panel which must rebuild
    G0/Z from raw elicitation data per row.
    """
    t_start = time.time()
    rows_out = []
    n_units = 0
    n_equal = 0
    n_differ = 0
    n_not_recomputed = 0
    not_recomputed_reasons: dict[str, int] = {}
    differences = []
    n_sanity_ok = 0
    n_sanity_fail = 0
    sanity_examples = []

    total_str = str(total_hint) if total_hint is not None else "?"
    print(f"[{out_csv.name}] up to {total_str} units to recompute (lazily built)", flush=True)

    it = iter(units)
    i = -1
    while True:
        if time.time() - t_start > time_budget_s:
            print(f"[{out_csv.name}] time budget ({time_budget_s}s) exceeded after "
                  f"{i+1}/{total_str}; stopping.", flush=True)
            break
        i += 1
        try:
            with per_unit_timeout(build_timeout_s):
                u = next(it)
        except StopIteration:
            break
        except _PerUnitTimeout:
            n_not_recomputed += 1
            not_recomputed_reasons["instance_build_timeout"] = (
                not_recomputed_reasons.get("instance_build_timeout", 0) + 1
            )
            print(f"  [{out_csv.name}] instance build #{i} exceeded {build_timeout_s}s; "
                  f"skipping and abandoning this generator (state may be inconsistent).",
                  flush=True)
            break

        n_units += 1
        do_sanity = i < sanity_check_n
        res = process_unit(u["cpdag"], u["g0"], u["x"], u["y"], u["z_star"],
                            u["stored_radius"], do_sanity=do_sanity)

        if res["sanity_status"] == "match":
            n_sanity_ok += 1
        elif res["sanity_status"] == "MISMATCH":
            n_sanity_fail += 1
            sanity_examples.append({**{k: u.get(k) for k in id_fields},
                                     "stored": u["stored_radius"],
                                     "recomputed_backdoor": res["backdoor_recheck"]})

        if res["sanity_status"] == "not_recomputed_timeout" and res["gac_radius"] is None:
            n_not_recomputed += 1
            not_recomputed_reasons["gac_per_unit_timeout"] = (
                not_recomputed_reasons.get("gac_per_unit_timeout", 0) + 1
            )
        elif not res["gac_exact"]:
            n_not_recomputed += 1
            not_recomputed_reasons["gac_bfs_inexact"] = (
                not_recomputed_reasons.get("gac_bfs_inexact", 0) + 1
            )
        elif res["equal"]:
            n_equal += 1
        else:
            n_differ += 1
            differences.append({**{k: u.get(k) for k in id_fields},
                                 "z": json.loads(res["z_star"]),
                                 "stored_radius": res["stored_radius"],
                                 "gac_radius": res["gac_radius"]})

        rows_out.append({**{k: u.get(k) for k in id_fields}, **res})

        if (i + 1) % 10 == 0:
            elapsed = time.time() - t_start
            print(f"[{out_csv.name}] {i+1}/{total_str} done, n_equal={n_equal} "
                  f"n_differ={n_differ} elapsed={elapsed:.0f}s", flush=True)

    fieldnames = id_fields + [
        "z_star", "stored_radius", "gac_radius", "gac_exact", "equal",
        "z_gac_valid_at_g0", "backdoor_recheck", "sanity_status",
    ]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows_out:
            w.writerow(r)

    seconds = time.time() - t_start
    return {
        "n_units": n_units, "n_equal": n_equal, "n_differ": n_differ,
        "n_not_recomputed": n_not_recomputed,
        "not_recomputed_reasons": not_recomputed_reasons,
        "n_sanity_checked": n_sanity_ok + n_sanity_fail,
        "n_sanity_ok": n_sanity_ok, "n_sanity_fail": n_sanity_fail,
        "sanity_examples": sanity_examples,
        "differences": differences[:50], "n_differences_total": len(differences),
        "seconds": seconds,
    }


# --------------------------------------------------------------------------
# Dataset 1: LLM survival panel (axis_robustness_llm, axis_robustness_llm_v2)
# --------------------------------------------------------------------------


def run_llm_dataset(
    csv_path: Path,
    z_mode: str,
    out_csv: Path,
    *,
    limit: int | None,
    time_budget_s: float,
    sanity_check_n: int,
):
    """z_mode: 'commit' (z_star = rs.commit_z_star(g0,x,y)) or 'empty' (z_star = frozenset())."""
    t_start = time.time()
    bundle, _bundle_sha = ls.load_bundle()

    rows_out = []
    n_units = 0
    n_equal = 0
    n_differ = 0
    n_not_recomputed = 0
    not_recomputed_reasons: dict[str, int] = {}
    differences = []
    n_sanity_ok = 0
    n_sanity_fail = 0
    sanity_examples = []

    net_cache: dict[str, tuple] = {}  # network -> (dag, cpdag)
    g0_cache: dict[tuple, tuple] = {}  # (condition, network) -> (k, g0, g0_reason, shared)

    with open(csv_path) as f:
        reader = csv.DictReader(f)
        all_rows = list(reader)

    # Only rows the original experiment actually scored a radius for.
    scored_rows = [r for r in all_rows if r.get("status") == "ok"]
    if limit is not None:
        scored_rows = scored_rows[:limit]

    print(f"[{csv_path.name}] {len(scored_rows)} status==ok rows to recompute "
          f"(of {len(all_rows)} total)", flush=True)

    for i, row in enumerate(scored_rows):
        if time.time() - t_start > time_budget_s:
            n_not_recomputed += len(scored_rows) - i
            not_recomputed_reasons["time_budget_exceeded"] = (
                not_recomputed_reasons.get("time_budget_exceeded", 0) + len(scored_rows) - i
            )
            print(f"[{csv_path.name}] time budget ({time_budget_s}s) exceeded at row {i}/"
                  f"{len(scored_rows)}; stopping.", flush=True)
            break

        cond = row["condition"]
        net = row["network"]
        x, y = row["x"], row["y"]
        stored_radius_raw = row["radius"]

        try:
            if net not in net_cache:
                net_cache[net] = rs.load_network(net)
            dag, cpdag = net_cache[net]

            if len(cpdag.nodes) > NETWORK_NODE_LIMIT:
                n_not_recomputed += 1
                not_recomputed_reasons["network_too_large_skipped"] = (
                    not_recomputed_reasons.get("network_too_large_skipped", 0) + 1
                )
                continue

            key = (cond, net)
            if key not in g0_cache:
                g0_cache[key] = ls.build_state(bundle, cond, net, dag, cpdag)
            k, g0, g0_reason, shared = g0_cache[key]

            if g0 is None or g0_reason != "ok":
                n_not_recomputed += 1
                not_recomputed_reasons[f"g0_{g0_reason}"] = (
                    not_recomputed_reasons.get(f"g0_{g0_reason}", 0) + 1
                )
                continue

            if z_mode == "commit":
                z_star, z_reason = rs.commit_z_star(g0, x, y)
                if z_star is None:
                    n_not_recomputed += 1
                    not_recomputed_reasons[f"z_{z_reason}"] = (
                        not_recomputed_reasons.get(f"z_{z_reason}", 0) + 1
                    )
                    continue
            else:
                z_star = frozenset()

            n_units += 1

            timed_out = False
            try:
                with per_unit_timeout(PER_UNIT_TIMEOUT_S):
                    gac_rad, gac_exact, gac_method = gac_radius(cpdag, g0, x, y, z_star)
                    g0_gac_valid = gac_valid_at_g0(g0, x, y, z_star)
            except _PerUnitTimeout:
                timed_out = True
                gac_rad, gac_exact = None, False
                g0_gac_valid = None

            if timed_out:
                n_not_recomputed += 1
                not_recomputed_reasons["gac_per_unit_timeout"] = (
                    not_recomputed_reasons.get("gac_per_unit_timeout", 0) + 1
                )
                rows_out.append({
                    "network": net, "x": x, "y": y, "condition": cond,
                    "frame_row_id": row.get("frame_row_id", ""),
                    "shard_id": row.get("shard_id", ""),
                    "z_star": json.dumps(sorted(z_star)),
                    "stored_radius": (
                        int(stored_radius_raw) if stored_radius_raw not in ("", None) else None
                    ),
                    "gac_radius": None, "gac_exact": False, "equal": None,
                    "z_gac_valid_at_g0": None, "backdoor_recheck": None,
                    "sanity_status": "not_recomputed_timeout",
                })
                continue

            do_sanity = i < sanity_check_n
            bd_recomputed = None
            sanity_status = ""
            if do_sanity:
                try:
                    with per_unit_timeout(PER_UNIT_TIMEOUT_S):
                        bd_rad, bd_exact, bd_method = bd_radius(cpdag, g0, x, y, z_star)
                    bd_recomputed = bd_rad
                except _PerUnitTimeout:
                    bd_rad, bd_exact = None, False
                    bd_recomputed = None
                stored_radius_val = (
                    int(stored_radius_raw) if stored_radius_raw not in ("", None) else None
                )
                if not bd_exact:
                    sanity_status = "not_recomputed_timeout"
                elif stored_radius_val is not None and bd_rad == stored_radius_val:
                    sanity_status = "match"
                    n_sanity_ok += 1
                else:
                    sanity_status = "MISMATCH"
                    n_sanity_fail += 1
                    sanity_examples.append({
                        "network": net, "x": x, "y": y, "condition": cond,
                        "stored": stored_radius_val, "recomputed_backdoor": bd_rad,
                    })

            stored_radius_val = (
                int(stored_radius_raw) if stored_radius_raw not in ("", None) else None
            )
            equal = (stored_radius_val == gac_rad) if gac_exact else None
            if gac_exact:
                if equal:
                    n_equal += 1
                else:
                    n_differ += 1
                    differences.append({
                        "network": net, "x": x, "y": y, "condition": cond,
                        "z": sorted(z_star), "stored_radius": stored_radius_val,
                        "gac_radius": gac_rad,
                    })
            else:
                n_not_recomputed += 1
                not_recomputed_reasons["gac_bfs_inexact"] = (
                    not_recomputed_reasons.get("gac_bfs_inexact", 0) + 1
                )

            rows_out.append({
                "network": net, "x": x, "y": y, "condition": cond,
                "frame_row_id": row.get("frame_row_id", ""),
                "shard_id": row.get("shard_id", ""),
                "z_star": " ".join(f"{a}->{b}" for a, b in sorted(z_star)) if False else json.dumps(sorted(z_star)),
                "stored_radius": stored_radius_val,
                "gac_radius": gac_rad,
                "gac_exact": gac_exact,
                "equal": equal,
                "z_gac_valid_at_g0": g0_gac_valid,
                "backdoor_recheck": bd_recomputed,
                "sanity_status": sanity_status,
            })

            if (i + 1) % 10 == 0:
                elapsed = time.time() - t_start
                print(f"[{csv_path.name}] {i+1}/{len(scored_rows)} done, "
                      f"n_equal={n_equal} n_differ={n_differ} elapsed={elapsed:.0f}s", flush=True)

        except Exception as exc:  # noqa: BLE001
            n_not_recomputed += 1
            not_recomputed_reasons[f"exception_{type(exc).__name__}"] = (
                not_recomputed_reasons.get(f"exception_{type(exc).__name__}", 0) + 1
            )
            print(f"  EXCEPTION on {net}/{x}/{y}/{cond}: {exc!r}", flush=True)
            continue

    with open(out_csv, "w", newline="") as f:
        fieldnames = [
            "network", "x", "y", "condition", "frame_row_id", "shard_id",
            "z_star", "stored_radius", "gac_radius", "gac_exact", "equal",
            "z_gac_valid_at_g0", "backdoor_recheck", "sanity_status",
        ]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows_out:
            w.writerow(r)

    seconds = time.time() - t_start
    return {
        "n_units": n_units,
        "n_equal": n_equal,
        "n_differ": n_differ,
        "n_not_recomputed": n_not_recomputed,
        "not_recomputed_reasons": not_recomputed_reasons,
        "n_sanity_checked": n_sanity_ok + n_sanity_fail,
        "n_sanity_ok": n_sanity_ok,
        "n_sanity_fail": n_sanity_fail,
        "sanity_examples": sanity_examples,
        "differences": differences[:50],
        "n_differences_total": len(differences),
        "seconds": seconds,
    }


# --------------------------------------------------------------------------
# Dataset 1b: parallel recompute of a network-filtered subset of the LLM
# panel (e.g. the large networks skipped by NETWORK_NODE_LIMIT above). One
# unit per multiprocessing task; each worker process rebuilds (cpdag, g0,
# z_star) itself via the *same* rs.load_network / ls.build_state /
# rs.commit_z_star calls the sequential path uses, and scores it with the
# same gac_radius / gac_valid_at_g0 functions, so results are directly
# comparable to axis_robustness_llm.csv.
# --------------------------------------------------------------------------

_worker_bundle = None
_worker_net_cache: dict = {}
_worker_g0_cache: dict = {}


def _worker_init():
    """multiprocessing.Pool initializer: load the elicitation bundle once per
    worker process; reused (along with per-network caches below) across every
    task that worker subsequently handles."""
    global _worker_bundle
    _worker_bundle, _ = ls.load_bundle()


def _worker_recompute_unit(task):
    """Run in a worker process. task = (row_dict, per_unit_timeout_s).

    Returns a dict with the same fields run_llm_dataset produces per row,
    plus "_not_recomputed_reason" (None on success) that the caller strips
    before writing the CSV.
    """
    row, timeout_s = task
    cond, net, x, y = row["condition"], row["network"], row["x"], row["y"]
    stored_radius_raw = row.get("radius", "")
    stored_radius_val = (
        int(stored_radius_raw) if stored_radius_raw not in ("", None) else None
    )
    out = {
        "network": net, "x": x, "y": y, "condition": cond,
        "frame_row_id": row.get("frame_row_id", ""),
        "shard_id": row.get("shard_id", ""),
        "z_star": "", "stored_radius": stored_radius_val,
        "gac_radius": None, "gac_exact": False, "equal": None,
        "z_gac_valid_at_g0": None, "backdoor_recheck": None,
        "sanity_status": "",
        "_not_recomputed_reason": None,
    }
    try:
        if net not in _worker_net_cache:
            _worker_net_cache[net] = rs.load_network(net)
        dag, cpdag = _worker_net_cache[net]

        key = (cond, net)
        if key not in _worker_g0_cache:
            _worker_g0_cache[key] = ls.build_state(_worker_bundle, cond, net, dag, cpdag)
        _k, g0, g0_reason, _shared = _worker_g0_cache[key]
        if g0 is None or g0_reason != "ok":
            out["_not_recomputed_reason"] = f"g0_{g0_reason}"
            return out

        z_star, z_reason = rs.commit_z_star(g0, x, y)
        if z_star is None:
            out["_not_recomputed_reason"] = f"z_{z_reason}"
            return out
        out["z_star"] = json.dumps(sorted(z_star))

        try:
            with per_unit_timeout(timeout_s):
                gac_rad, gac_exact, _method = gac_radius(cpdag, g0, x, y, z_star)
                g0_gac_valid = gac_valid_at_g0(g0, x, y, z_star)
        except _PerUnitTimeout:
            out["sanity_status"] = "not_recomputed_timeout"
            out["_not_recomputed_reason"] = "gac_per_unit_timeout"
            return out

        out["gac_radius"] = gac_rad
        out["gac_exact"] = gac_exact
        out["z_gac_valid_at_g0"] = g0_gac_valid
        if gac_exact:
            out["equal"] = (stored_radius_val == gac_rad)
        else:
            out["_not_recomputed_reason"] = "gac_bfs_inexact"
        return out
    except Exception as exc:  # noqa: BLE001
        out["_not_recomputed_reason"] = f"exception_{type(exc).__name__}"
        return out


def run_llm_dataset_parallel(
    csv_path: Path,
    out_csv: Path,
    *,
    only_networks: set[str] | None,
    max_nodes: int,
    workers: int,
    per_unit_timeout: float,
    time_budget_s: float,
    limit: int | None = None,
):
    """Parallel counterpart to run_llm_dataset, for a network-filtered subset
    (e.g. the networks that NETWORK_NODE_LIMIT skips in the sequential run).

    max_nodes <= 0 disables the node-count skip entirely. Reuses the exact
    same rebuild (rs.load_network / ls.build_state / rs.commit_z_star) and
    scoring (gac_radius / gac_valid_at_g0) functions as run_llm_dataset, via
    _worker_recompute_unit, so results are directly comparable.
    """
    t_start = time.time()
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        all_rows = list(reader)

    scored_rows = [r for r in all_rows if r.get("status") == "ok"]
    if only_networks:
        scored_rows = [r for r in scored_rows if r["network"] in only_networks]
    if max_nodes and max_nodes > 0:
        node_counts: dict[str, int] = {}
        kept = []
        for r in scored_rows:
            net = r["network"]
            if net not in node_counts:
                _dag, cpdag = rs.load_network(net)
                node_counts[net] = len(cpdag.nodes)
            if node_counts[net] <= max_nodes:
                kept.append(r)
        scored_rows = kept
    if limit is not None:
        scored_rows = scored_rows[:limit]

    tasks = [(r, per_unit_timeout) for r in scored_rows]
    n_total = len(tasks)
    print(f"[{out_csv.name}] {n_total} units to recompute in parallel "
          f"(workers={workers}, per_unit_timeout={per_unit_timeout}s, "
          f"time_budget={time_budget_s}s)", flush=True)

    rows_out = []
    n_equal = 0
    n_differ = 0
    n_not_recomputed = 0
    not_recomputed_reasons: dict[str, int] = {}
    differences = []
    n_done = 0
    hit_time_budget = False

    pool = mp.Pool(processes=workers, initializer=_worker_init)
    try:
        for res in pool.imap_unordered(_worker_recompute_unit, tasks):
            n_done += 1
            reason = res.pop("_not_recomputed_reason", None)
            if reason is not None:
                n_not_recomputed += 1
                not_recomputed_reasons[reason] = not_recomputed_reasons.get(reason, 0) + 1
            elif res["equal"]:
                n_equal += 1
            else:
                n_differ += 1
                differences.append({
                    "network": res["network"], "x": res["x"], "y": res["y"],
                    "condition": res["condition"],
                    "z": json.loads(res["z_star"]) if res["z_star"] else None,
                    "stored_radius": res["stored_radius"], "gac_radius": res["gac_radius"],
                })
            rows_out.append(res)

            if n_done % 5 == 0 or n_done == n_total:
                elapsed = time.time() - t_start
                print(f"[{out_csv.name}] {n_done}/{n_total} done, n_equal={n_equal} "
                      f"n_differ={n_differ} n_not_recomputed={n_not_recomputed} "
                      f"elapsed={elapsed:.0f}s", flush=True)

            if time.time() - t_start > time_budget_s:
                print(f"[{out_csv.name}] overall time budget ({time_budget_s}s) exceeded "
                      f"after {n_done}/{n_total}; terminating remaining workers.", flush=True)
                hit_time_budget = True
                break
    finally:
        pool.terminate()
        pool.join()

    if hit_time_budget:
        remaining = n_total - n_done
        n_not_recomputed += remaining
        not_recomputed_reasons["time_budget_exceeded"] = (
            not_recomputed_reasons.get("time_budget_exceeded", 0) + remaining
        )

    fieldnames = [
        "network", "x", "y", "condition", "frame_row_id", "shard_id",
        "z_star", "stored_radius", "gac_radius", "gac_exact", "equal",
        "z_gac_valid_at_g0", "backdoor_recheck", "sanity_status",
    ]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows_out:
            w.writerow({k: r.get(k) for k in fieldnames})

    seconds = time.time() - t_start
    return {
        "n_units": n_total,
        "n_done": n_done,
        "n_remaining": n_total - n_done,
        "n_equal": n_equal,
        "n_differ": n_differ,
        "n_not_recomputed": n_not_recomputed,
        "not_recomputed_reasons": not_recomputed_reasons,
        "differences": differences[:50],
        "n_differences_total": len(differences),
        "seconds": seconds,
    }


# --------------------------------------------------------------------------
# Dataset 2: synthetic survival sweep (axis_robustness_p6): flip + tiered
# --------------------------------------------------------------------------


def run_p6_dataset(csv_path: Path, out_csv: Path, *, limit, time_budget_s, sanity_check_n):
    from bkrobust.demo.graph import MPDAG as _MPDAG  # noqa: F401
    from bkrobust.robustness import survival as sv

    with open(csv_path) as f:
        reader = csv.DictReader(f)
        all_rows = list(reader)
    if limit is not None:
        all_rows = all_rows[:limit]

    def gen():
        for row in all_rows:
            arm = row["arm"]
            seed = int(row["seed"])
            component_size = int(row["component_size"])
            separation = int(row["separation"])
            if arm == "flip":
                coverage = float(row["coverage"])
                base_wrongness = float(row["base_wrongness"])
                inst, reason = sv.build_flip_instance(
                    component_size, separation, coverage, base_wrongness, seed
                )
            elif arm == "tiered":
                n_tiers = int(row["n_tiers"])
                inst, reason = sv.build_tiered_instance(component_size, separation, seed, n_tiers)
            else:
                continue
            if inst is None:
                print(f"  SKIP {row['instance_id']}: rebuild failed ({reason})", flush=True)
                continue
            z_csv = frozenset(s for s in row["z_star"].split(",") if s)
            if z_csv != inst["z_star"]:
                print(f"  WARNING {row['instance_id']}: rebuilt z_star {sorted(inst['z_star'])} "
                      f"!= csv z_star {sorted(z_csv)}", flush=True)
            stored_raw = row["r_val"]
            stored_val = int(stored_raw) if stored_raw not in ("", None) else None
            yield {
                "instance_id": row["instance_id"], "arm": arm,
                "cpdag": inst["cpdag"], "g0": inst["g0"],
                "x": inst["x"], "y": inst["y"], "z_star": inst["z_star"],
                "stored_radius": stored_val,
            }

    return run_generic_dataset(
        gen(), out_csv, time_budget_s=time_budget_s, sanity_check_n=sanity_check_n,
        id_fields=["instance_id", "arm"], total_hint=len(all_rows),
    )


# --------------------------------------------------------------------------
# Dataset 3: uniform-query (unconfounded) sweep, 480 tiered instances
# --------------------------------------------------------------------------


def run_unconfounded_dataset(csv_path: Path, out_csv: Path, *, limit, time_budget_s, sanity_check_n):
    from bkrobust.robustness import unconfounded_generator as ug

    with open(csv_path) as f:
        reader = csv.DictReader(f)
        all_rows = list(reader)
    if limit is not None:
        all_rows = all_rows[:limit]

    def gen():
        for row in all_rows:
            component_size = int(row["component_size"])
            seed = int(row["seed"])
            n_tiers = int(row["n_tiers"])
            inst, reason = ug.build_unconfounded_tiered_instance(component_size, seed, n_tiers)
            if inst is None:
                print(f"  SKIP {row['instance_id']}: rebuild failed ({reason})", flush=True)
                continue
            z_csv = frozenset(s for s in row["z_star"].split(",") if s)
            if z_csv != inst["z_star"]:
                print(f"  WARNING {row['instance_id']}: rebuilt z_star {sorted(inst['z_star'])} "
                      f"!= csv z_star {sorted(z_csv)}", flush=True)
            stored_raw = row["r_val"]
            stored_val = int(stored_raw) if stored_raw not in ("", None) else None
            yield {
                "instance_id": row["instance_id"], "arm": row["arm"],
                "cpdag": inst["cpdag"], "g0": inst["g0"],
                "x": inst["x"], "y": inst["y"], "z_star": inst["z_star"],
                "stored_radius": stored_val,
            }

    return run_generic_dataset(
        gen(), out_csv, time_budget_s=time_budget_s, sanity_check_n=sanity_check_n,
        id_fields=["instance_id", "arm"], total_hint=len(all_rows),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset", choices=["llm", "llm_v2", "p6", "unconfounded", "all"])
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--time-budget", type=float, default=600.0)
    ap.add_argument("--sanity-n", type=int, default=25)
    ap.add_argument(
        "--only-networks", type=str, default=None,
        help="Comma-separated network names restricting the llm dataset "
             "(e.g. diabetes,munin1). Triggers the parallel path and writes "
             "to axis_robustness_llm_large.csv / summary key "
             "axis_robustness_llm_large instead of the default output.",
    )
    ap.add_argument(
        "--max-nodes", type=int, default=None,
        help="Override NETWORK_NODE_LIMIT for the llm dataset's parallel "
             "path; 0 disables the node-count skip entirely.",
    )
    ap.add_argument(
        "--workers", type=int, default=1,
        help="If >1, recompute the llm dataset's units in parallel across "
             "this many worker processes (one unit per task).",
    )
    ap.add_argument(
        "--per-unit-timeout", type=float, default=None,
        help="Per-unit GAC BFS timeout in seconds. Default: 8s on the "
             "sequential path, 1200s (20 min) on the parallel path.",
    )
    args = ap.parse_args()

    summary = {}

    if args.dataset in ("llm", "all"):
        only_networks = (
            {n.strip() for n in args.only_networks.split(",") if n.strip()}
            if args.only_networks else None
        )
        if only_networks or args.workers > 1:
            per_unit_timeout_val = (
                args.per_unit_timeout if args.per_unit_timeout is not None else 1200.0
            )
            max_nodes_val = args.max_nodes if args.max_nodes is not None else NETWORK_NODE_LIMIT
            out_csv = OUT_DIR / (
                "axis_robustness_llm_large.csv" if only_networks else "axis_robustness_llm.csv"
            )
            summary_key = "axis_robustness_llm_large" if only_networks else "axis_robustness_llm"
            res = run_llm_dataset_parallel(
                Path("results/axis_robustness_llm/analysis_units.csv"),
                out_csv,
                only_networks=only_networks,
                max_nodes=max_nodes_val,
                workers=max(args.workers, 1),
                per_unit_timeout=per_unit_timeout_val,
                time_budget_s=args.time_budget,
                limit=args.limit,
            )
            summary[summary_key] = res
            print(json.dumps(res, indent=2, default=str))
        else:
            res = run_llm_dataset(
                Path("results/axis_robustness_llm/analysis_units.csv"),
                "commit",
                OUT_DIR / "axis_robustness_llm.csv",
                limit=args.limit,
                time_budget_s=args.time_budget,
                sanity_check_n=args.sanity_n,
            )
            summary["axis_robustness_llm"] = res
            print(json.dumps(res, indent=2, default=str))

    if args.dataset in ("llm_v2", "all"):
        res = run_llm_dataset(
            Path("results/axis_robustness_llm_v2/empty_valid_units.csv"),
            "empty",
            OUT_DIR / "axis_robustness_llm_v2_empty_valid.csv",
            limit=args.limit,
            time_budget_s=args.time_budget,
            sanity_check_n=args.sanity_n,
        )
        summary["axis_robustness_llm_v2_empty_valid"] = res
        print(json.dumps(res, indent=2, default=str))

    if args.dataset in ("p6", "all"):
        res = run_p6_dataset(
            Path("results/axis_robustness_p6/survival_instances.csv"),
            OUT_DIR / "axis_robustness_p6.csv",
            limit=args.limit,
            time_budget_s=args.time_budget,
            sanity_check_n=args.sanity_n,
        )
        summary["axis_robustness_p6"] = res
        print(json.dumps(res, indent=2, default=str))

    if args.dataset in ("unconfounded", "all"):
        res = run_unconfounded_dataset(
            Path("results/axis_robustness_unconfounded_v2/survival_instances.csv"),
            OUT_DIR / "axis_robustness_unconfounded_v2.csv",
            limit=args.limit,
            time_budget_s=args.time_budget,
            sanity_check_n=args.sanity_n,
        )
        summary["axis_robustness_unconfounded_v2"] = res
        print(json.dumps(res, indent=2, default=str))

    summary_path = OUT_DIR / "summary.partial.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()
