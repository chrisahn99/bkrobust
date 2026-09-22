"""Claim-level baselines for the real corpus ([Point 6]: "the SHD baseline is
insufficient").

This module computes two predictors that do not exist anywhere else in the
repository, on the **frozen 831-row frame** at
``results/axis_robustness_real/frame.jsonl`` (never modified, never re-hashed
to something else):

``phi_1``
    Single-claim deletion fragility: the fraction of the analyst's asserted
    claims ``K`` whose individual retraction, followed by re-closing under
    Meek, breaks the GAC-validity of the frame's own held-fixed ``Z*`` at the
    frame's own ``G0``. Canonical definition: ``experiments/stage0_claim_radius.py``
    lines ~140-155 (the per-claim loop inside its coverage-1.0 check). That
    script only ever ran this at coverage 1.0; here it is applied verbatim to
    every one of the 831 admissible rows, at whatever coverage the row itself
    carries, using **that row's own** ``K = select_knowledge(dag, cpdag,
    coverage)`` and **that row's own** ``z_star`` -- both reconstructed
    deterministically (``select_knowledge`` takes no RNG) and checked against
    the frame's own ``n_k`` and ``z_star`` columns before anything is trusted.

``r_claim``
    Claim-level radius: the minimum number of claim **retractions**
    (deletions, never reversals) whose *joint* removal breaks ``Z*``'s
    validity. Found by exhaustive depth-first-by-depth search over subsets of
    ``K``, depth 1, 2, 3, ... A search is **exact** only up to the depth it
    actually completes; deeper depths are skipped, and the row is marked
    ``r_claim_status = "censored_depth_cap"``, whenever the combinatorial
    cost of the next depth is not affordable on this machine (see
    :data:`DEPTH_BUDGET_SECONDS` and the cost estimate in
    :func:`compute_baselines`). A row that is exhausted all the way to
    ``|K|`` without ever breaking validity is a **different, informative**
    status, ``"never_invalidated_up_to_full_retraction"`` -- not a censoring
    artifact, a fact about that query: no subset of the analyst's claims,
    however large, changes whether ``Z*`` is valid.

Known degeneracy, stated up front rather than discovered late (task brief):
at coverage 1.0, ``K`` for a given row is **exactly** the ``k_true`` that
``benchmarks.measure.fast_gate`` already tested claim-by-claim as its
admission criterion, and the frame's own ``z_star`` at that coverage is the
DAG's true optimal adjustment set. So the admission gate **is** the depth-1
``phi_1``/``r_claim`` test at coverage 1.0, and every one of the 543
admissible coverage-1.0 rows is expected to resolve at depth 1 with
``r_claim == 1`` and ``phi_1 > 0``, **by construction**, not by any property
of the radius axis. This is verified empirically below (never assumed) and
reported loudly if it fails anywhere.

Two baselines this module does **not** compute, because they already exist
and the task brief is explicit that they must not be reimplemented:

- ``separation``: already a column on every analysis unit
  (:mod:`bkrobust.robustness.real_analyse`), carried through from the shard
  instance rows. This module does not touch it; ``__main__`` below simply
  adds it to the set of predictors that get *ranked*, which the committed
  analysis never did (it only used it to stratify).
- ``k_g0``: already one of :data:`bkrobust.robustness.real_analyse.PREDICTORS`
  and already ranked in ``analysis_tau.csv``. Carried straight through into
  the extended table this module writes.

Reads
-----
``results/axis_robustness_real/frame.jsonl`` (hash-verified via
:func:`bkrobust.robustness.real_analyse.verify_frame`, imported, never
duplicated) and the network files under
``results/axisa3/networks/example_models/``. Also reads every completed
shard via :func:`bkrobust.robustness.real_analyse.load_all_shards` /
:func:`bkrobust.robustness.real_analyse.build_analysis_units` to build the
extended tau table -- read-only; nothing under ``results/axis_robustness_real/``
is ever written by this module.

Writes, all under ``--out-dir`` (default ``results/axis_robustness_real_p6``)
--------------------------------------------------------------------------
``baselines_per_instance.csv`` (831 rows, one per frame row), ``analysis_tau_extended.csv``
(the committed four predictors plus ``separation``, ``phi_1``, ``r_claim``,
ranked by the same :func:`bkrobust.robustness.real_analyse.tau_for_stratum`
the committed analysis uses -- imported, not reimplemented), and
``baselines_manifest.json`` (depth caps, timings, the empirical check of the
degeneracy above).

    PYTHONPATH=src /usr/bin/python3 src/bkrobust/robustness/real_baselines_p6.py
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

from bkrobust.benchmarks.describe import parse_file
from bkrobust.benchmarks.measure import select_knowledge
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.gac import is_gac_valid_mpdag
from bkrobust.robustness.real_analyse import (
    DEFAULT_N_BOOT,
    DEFAULT_SEED,
    PREDICTORS as EXISTING_PREDICTORS,
    build_analysis_units,
    group_units_by_stratum,
    load_all_shards,
    tau_for_stratum,
    write_csv,
    write_json,
)
from bkrobust.robustness.real_analyse import verify_frame as _verify_frame

DEFAULT_NETWORKS_DIR = "results/axisa3/networks/example_models"
DEFAULT_FRAME_DIR = "results/axis_robustness_real"
DEFAULT_OUT_DIR = "results/axis_robustness_real_p6"

#: Networks with |K| at or below this size get an **exhaustive** depth
#: search, up to depth |K| itself. Every (network, coverage) cell on this
#: corpus with |K| > 12 is coverage 1.0 (diabetes 26, arth150 29, pathfinder
#: 79 -- see ``results/axis_robustness_real_p6/NOTES.md``), where depth 1 is
#: expected to resolve every row by construction (see module docstring), so
#: this threshold costs nothing in practice and only bounds the worst case.
SMALL_K_EXHAUSTIVE_MAX = 12

#: Wall-clock budget, per (network, coverage) cell, per depth, before that
#: depth is skipped and any row still unresolved is censored. Estimated from
#: the measured per-call cost of ``apply_orientations`` on that cell's own
#: depth-1 pass times ``math.comb(|K|, depth)`` -- an approximation (deeper
#: subsets are smaller, hence typically *cheaper* per call, so this is a
#: conservative -- i.e. more-likely-to-skip -- estimate).
DEPTH_BUDGET_SECONDS = 180.0

#: The new predictors this module adds to the ranking.
NEW_PREDICTORS: tuple[str, ...] = ("phi_1", "r_claim")

#: Every predictor the extended tau table ranks: the four the committed
#: analysis already ranks, plus ``separation`` (already a unit column, never
#: previously ranked), plus the two new ones.
EXTENDED_PREDICTORS: tuple[str, ...] = (*EXISTING_PREDICTORS, "separation", *NEW_PREDICTORS)

BASELINE_COLUMNS: list[str] = [
    "row_id", "network", "x", "y", "coverage", "n_k",
    "phi_1", "n_broken_single", "n_inconsistent_single",
    "r_claim", "r_claim_status", "r_claim_censored", "r_claim_depth_cap_used",
]

EXTENDED_TAU_COLUMNS: list[str] = [
    "stratum", "arm", "predictor", "endpoint", "status",
    "n", "n_networks", "stratum_n_units", "n_excluded_unreached", "n_excluded_nan",
    "tau_b", "p_value", "ci_lo_2p5", "ci_hi_97p5", "n_boot_nan",
    "loo_tau_min", "loo_tau_max", "loo_network_at_min", "loo_network_at_max", "loo_verdict_flips",
]


# ---------------------------------------------------------------------------
# Network loading (mirrors experiments/stage0_claim_radius.py::load exactly)
# ---------------------------------------------------------------------------


def load_network(name: str, networks_dir: Path) -> tuple[MPDAG, MPDAG]:
    """Parses one network file and returns its ground-truth DAG and CPDAG.

    Args:
        name: Network name as it appears on frame rows (e.g. ``"pathfinder"``).
        networks_dir: Directory of network source files.

    Returns:
        ``(dag, cpdag)``.
    """
    path = next(p for p in sorted(networks_dir.iterdir()) if p.name.startswith(name))
    parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
    dag = MPDAG(parsed.nodes, directed=parsed.edges)
    return dag, dag_to_cpdag(dag)


# ---------------------------------------------------------------------------
# phi_1 / r_claim
# ---------------------------------------------------------------------------


def compute_baselines(
    frame_rows: list[dict[str, Any]],
    networks_dir: Path,
    *,
    small_k_exhaustive_max: int = SMALL_K_EXHAUSTIVE_MAX,
    depth_budget_s: float = DEPTH_BUDGET_SECONDS,
    log=print,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Computes ``phi_1`` and ``r_claim`` for every frame row.

    Rows are grouped into ``(network, coverage)`` cells, since ``K`` and
    every claim-dropped ``G0`` variant depend only on the cell, not the
    individual ``(x, y)`` pair -- so each ``apply_orientations`` call is made
    once per cell and reused across every row that shares it. This is the
    difference between minutes and hours on ``pathfinder`` (``|K| = 79``,
    ``apply_orientations`` at ~2s/call there).

    Args:
        frame_rows: Rows of ``frame.jsonl``.
        networks_dir: Directory of network source files.
        small_k_exhaustive_max: See :data:`SMALL_K_EXHAUSTIVE_MAX`.
        depth_budget_s: See :data:`DEPTH_BUDGET_SECONDS`.
        log: Callable for verbose progress logging (default ``print``).

    Returns:
        ``(baseline_rows, cell_log)``. ``baseline_rows`` has one dict per
        frame row with :data:`BASELINE_COLUMNS` keys. ``cell_log`` has one
        dict per ``(network, coverage)`` cell with timing and resolution
        counts, for the manifest.
    """
    cells: dict[tuple[str, float], list[dict[str, Any]]] = defaultdict(list)
    for r in frame_rows:
        cells[(r["network"], r["coverage"])].append(r)

    dag_cache: dict[str, tuple[MPDAG, MPDAG]] = {}
    out_rows: list[dict[str, Any]] = []
    cell_log: list[dict[str, Any]] = []

    for (network, coverage), rows in sorted(cells.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        t_cell = time.perf_counter()
        if network not in dag_cache:
            t0 = time.perf_counter()
            dag_cache[network] = load_network(network, networks_dir)
            log(f"[load] {network}: {time.perf_counter() - t0:.2f}s")
        dag, cpdag = dag_cache[network]

        t0 = time.perf_counter()
        K = select_knowledge(dag, cpdag, coverage)
        n_k_expected = rows[0]["n_k"]
        if len(K) != n_k_expected:
            raise AssertionError(
                f"select_knowledge mismatch for {network}@{coverage}: "
                f"recomputed |K|={len(K)} vs frame n_k={n_k_expected} -- "
                "refusing to trust a K that does not reproduce the committed frame."
            )
        log(f"[select_knowledge] {network}@{coverage}: |K|={len(K)} "
            f"({time.perf_counter() - t0:.2f}s, {len(rows)} rows in cell)")

        # --- depth 1: one apply_orientations per claim, cached for the cell ---
        t0 = time.perf_counter()
        drop_g0: dict[int, MPDAG | None] = {}
        n_inconsistent = 0
        for i, claim in enumerate(K):
            g = apply_orientations(cpdag, [e for e in K if e != claim])
            drop_g0[i] = g
            if g is None:
                n_inconsistent += 1
        depth1_time = time.perf_counter() - t0
        avg_call_s = depth1_time / max(1, len(K))
        log(f"[depth1] {network}@{coverage}: |K|={len(K)} inconsistent={n_inconsistent} "
            f"time={depth1_time:.2f}s avg_call={avg_call_s * 1000:.1f}ms")

        unresolved: list[tuple[dict[str, Any], dict[str, Any]]] = []
        for r in rows:
            z = frozenset(r["z_star"])
            broken = 0
            for i in range(len(K)):
                g = drop_g0[i]
                if g is None:
                    continue
                if not is_gac_valid_mpdag(g, r["x"], r["y"], z):
                    broken += 1
            phi_1 = (broken / len(K)) if K else None
            row_out = {
                "row_id": r["row_id"], "network": network, "x": r["x"], "y": r["y"],
                "coverage": coverage, "n_k": len(K),
                "phi_1": phi_1, "n_broken_single": broken, "n_inconsistent_single": n_inconsistent,
                "r_claim": None, "r_claim_status": None,
                "r_claim_censored": None, "r_claim_depth_cap_used": 1,
            }
            if broken > 0:
                row_out["r_claim"] = 1
                row_out["r_claim_status"] = "resolved"
                row_out["r_claim_censored"] = False
            else:
                unresolved.append((r, row_out))
            out_rows.append(row_out)

        if coverage == 1.0 and unresolved:
            log(f"  *** ANOMALY *** {network}@coverage=1.0: {len(unresolved)} row(s) NOT resolved "
                "at depth 1, contradicting the construction argument (module docstring) that the "
                "admission gate already performed this exact test. Reported, not silently deepened.")

        # --- deeper depths, only for rows still unresolved after depth 1 ---
        depth_cap_for_cell = 1
        n_resolved_deeper = 0
        n_depths_attempted = 1
        if unresolved:
            max_feasible_depth = len(K) if len(K) <= small_k_exhaustive_max else 1
            d = 2
            remaining = unresolved
            while d <= max_feasible_depth and remaining:
                n_combos = math.comb(len(K), d)
                est_s = n_combos * avg_call_s
                if est_s > depth_budget_s:
                    log(f"[depth{d}] {network}@{coverage}: SKIPPED -- estimated {est_s:.0f}s > "
                        f"budget {depth_budget_s:.0f}s ({n_combos} combos of size {d}); "
                        f"{len(remaining)} row(s) remain censored at cap={d - 1}")
                    break
                t0 = time.perf_counter()
                before_this_depth = len(remaining)
                for combo in itertools.combinations(range(len(K)), d):
                    if not remaining:
                        break
                    kept = [K[j] for j in range(len(K)) if j not in combo]
                    g = apply_orientations(cpdag, kept)
                    if g is None:
                        continue
                    next_remaining = []
                    for (r, row_out) in remaining:
                        z = frozenset(r["z_star"])
                        if not is_gac_valid_mpdag(g, r["x"], r["y"], z):
                            row_out["r_claim"] = d
                            row_out["r_claim_status"] = "resolved"
                            row_out["r_claim_censored"] = False
                            row_out["r_claim_depth_cap_used"] = d
                            n_resolved_deeper += 1
                        else:
                            next_remaining.append((r, row_out))
                    remaining = next_remaining
                depth_time = time.perf_counter() - t0
                log(f"[depth{d}] {network}@{coverage}: {n_combos} combos, {depth_time:.2f}s, "
                    f"{before_this_depth - len(remaining)} newly resolved"
                    f" -- {len(remaining)} row(s) remain unresolved")
                depth_cap_for_cell = d
                n_depths_attempted = d
                unresolved = remaining
                d += 1

            for (r, row_out) in unresolved:
                row_out["r_claim_depth_cap_used"] = depth_cap_for_cell
                if depth_cap_for_cell >= len(K):
                    row_out["r_claim_status"] = "never_invalidated_up_to_full_retraction"
                    row_out["r_claim_censored"] = False
                else:
                    row_out["r_claim_status"] = "censored_depth_cap"
                    row_out["r_claim_censored"] = True

        cell_log.append({
            "network": network, "coverage": coverage, "n_k": len(K), "n_rows": len(rows),
            "n_inconsistent_single_drop": n_inconsistent,
            "n_resolved_depth1": len(rows) - len(unresolved) - n_resolved_deeper,
            "n_resolved_deeper": n_resolved_deeper,
            "n_unresolved_final": len(unresolved),
            "depth_cap_used": depth_cap_for_cell,
            "n_depths_attempted": n_depths_attempted,
            "small_k_exhaustive": len(K) <= small_k_exhaustive_max,
            "cell_seconds": round(time.perf_counter() - t_cell, 2),
        })
        log(f"[cell done] {network}@{coverage}: {round(time.perf_counter() - t_cell, 2)}s total")

    return out_rows, cell_log


# ---------------------------------------------------------------------------
# Extended tau table (reuses real_analyse.tau_for_stratum verbatim)
# ---------------------------------------------------------------------------


def build_extended_units(
    loaded_shards: list[dict[str, Any]], baseline_by_row_id: dict[str, dict[str, Any]]
) -> list[dict[str, Any]]:
    """Builds :func:`real_analyse.build_analysis_units` output, augmented with
    ``phi_1`` and ``r_claim``.

    ``phi_1``/``r_claim`` are attached **only** to flip-arm units at
    ``base_wrongness == 0.0``: that is the only case where the unit's
    operative ``K`` is exactly the frame's own ``K`` that
    :func:`compute_baselines` used. Every other unit (flip at higher base
    wrongness, every tiered unit) is built from a *different*, corrupted or
    independently-generated knowledge state that the frozen frame does not
    capture, so attaching the frame-derived value there would silently
    compare a stale predictor to an endpoint it was never computed against.
    Those units keep ``phi_1 = r_claim = None`` -- excluded from tau, not
    imputed -- and the exclusion is total for their strata, stated in
    ``NOTES.md``, not smoothed over here.

    Args:
        loaded_shards: Output of ``real_analyse.load_all_shards``.
        baseline_by_row_id: ``{frame_row_id: baseline_row}`` from
            :func:`compute_baselines`.

    Returns:
        Every unit row, each carrying ``phi_1`` and ``r_claim`` in addition
        to the columns ``real_analyse.build_analysis_units`` already sets.
    """
    units = build_analysis_units(loaded_shards)
    for u in units:
        u["phi_1"] = None
        u["r_claim"] = None
        if u["arm"] == "flip" and u.get("base_wrongness") == 0.0:
            b = baseline_by_row_id.get(u["frame_row_id"])
            if b is not None:
                u["phi_1"] = b["phi_1"]
                u["r_claim"] = b["r_claim"]
    return units


def build_extended_tau_table(
    stratum_groups: dict[str, list[dict[str, Any]]],
    *, n_boot: int = DEFAULT_N_BOOT, seed: int = DEFAULT_SEED,
) -> list[dict[str, Any]]:
    """One row per (stratum, predictor) over :data:`EXTENDED_PREDICTORS`.

    Endpoint is chosen per stratum's arm exactly as
    ``real_analyse.build_tau_table`` does. Calls
    ``real_analyse.tau_for_stratum`` directly -- the statistic itself is not
    reimplemented, only the predictor loop is widened.

    Args:
        stratum_groups: Output of ``real_analyse.group_units_by_stratum`` on
            the extended units.
        n_boot: Bootstrap resamples per cell.
        seed: Seed for every bootstrap.

    Returns:
        Rows in stratum, then predictor order.
    """
    from bkrobust.robustness.real_analyse import ENDPOINTS_BY_ARM

    rows = []
    for stratum in sorted(stratum_groups):
        units = stratum_groups[stratum]
        arm = units[0]["arm"]
        endpoints = ENDPOINTS_BY_ARM[arm]
        for predictor in EXTENDED_PREDICTORS:
            for endpoint in endpoints:
                r = tau_for_stratum(units, predictor, endpoint, n_boot=n_boot, seed=seed)
                r["stratum"], r["arm"] = stratum, arm
                rows.append(r)
    return rows


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    """Builds the CLI parser."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--frame-dir", default=DEFAULT_FRAME_DIR)
    p.add_argument("--networks-dir", default=DEFAULT_NETWORKS_DIR)
    p.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    p.add_argument("--n-boot", type=int, default=DEFAULT_N_BOOT)
    return p


def main(argv: list[str] | None = None) -> int:
    """Entry point: computes ``phi_1``/``r_claim`` on the frame, then the
    extended tau table, writing everything under ``--out-dir``.

    Returns:
        Process exit status. Nonzero if the frame hash does not verify.
    """
    args = build_arg_parser().parse_args(argv)
    frame_dir = Path(args.frame_dir)
    networks_dir = Path(args.networks_dir)
    out_dir = Path(args.out_dir)
    n_boot = args.n_boot
    seed = DEFAULT_SEED

    t_start = time.perf_counter()
    try:
        frame_rows, frame_info = _verify_frame(frame_dir)
    except (FileNotFoundError, RuntimeError) as exc:
        print(f"REFUSING TO RUN: {exc}")
        return 1
    print(f"frame verified: {len(frame_rows)} rows, sha256={frame_info['frame_sha256'][:16]}...")

    baseline_rows, cell_log = compute_baselines(frame_rows, networks_dir, log=print)
    write_csv(out_dir / "baselines_per_instance.csv", BASELINE_COLUMNS, baseline_rows)
    print(f"wrote baselines_per_instance.csv: {len(baseline_rows)} rows")

    baseline_by_row_id = {b["row_id"]: b for b in baseline_rows}

    loaded_shards, skipped_shards = load_all_shards(frame_dir)
    print(f"loaded {len(loaded_shards)} shards ({len(skipped_shards)} skipped)")
    units = build_extended_units(loaded_shards, baseline_by_row_id)
    stratum_groups = group_units_by_stratum(units)
    tau_rows = build_extended_tau_table(stratum_groups, n_boot=n_boot, seed=seed)
    write_csv(out_dir / "analysis_tau_extended.csv", EXTENDED_TAU_COLUMNS, tau_rows)
    print(f"wrote analysis_tau_extended.csv: {len(tau_rows)} rows "
          f"({len(stratum_groups)} strata x {len(EXTENDED_PREDICTORS)} predictors)")

    # Empirical check of the depth-1-at-coverage-1.0 degeneracy claim (module docstring).
    cov1_rows = [b for b in baseline_rows if b["coverage"] == 1.0]
    n_cov1_depth1 = sum(1 for b in cov1_rows if b["r_claim"] == 1)
    degeneracy_check = {
        "n_coverage_1_rows": len(cov1_rows),
        "n_resolved_at_depth_1": n_cov1_depth1,
        "all_resolved_at_depth_1": n_cov1_depth1 == len(cov1_rows),
        "phi_1_positive_count": sum(1 for b in cov1_rows if (b["phi_1"] or 0) > 0),
    }
    print("degeneracy check (coverage=1.0):", json.dumps(degeneracy_check))

    manifest = {
        "out_dir": str(out_dir),
        "frame_dir": str(frame_dir),
        "networks_dir": str(networks_dir),
        "n_boot": n_boot,
        "seed": seed,
        "small_k_exhaustive_max": SMALL_K_EXHAUSTIVE_MAX,
        "depth_budget_seconds": DEPTH_BUDGET_SECONDS,
        "frame_sha256": frame_info["frame_sha256"],
        "n_frame_rows": len(frame_rows),
        "n_baseline_rows": len(baseline_rows),
        "cell_log": cell_log,
        "degeneracy_check_coverage_1": degeneracy_check,
        "extended_predictors": list(EXTENDED_PREDICTORS),
        "n_strata": len(stratum_groups),
        "total_seconds": round(time.perf_counter() - t_start, 1),
    }
    write_json(out_dir / "baselines_manifest.json", manifest)
    print(f"wrote baselines_manifest.json. total time {manifest['total_seconds']}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
