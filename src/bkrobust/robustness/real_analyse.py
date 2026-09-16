"""Analysis of the real-structure survival sweep ([RE-11]).

Reads
-----
``results/axis_robustness_real/frame.jsonl`` (hash-checked against
``frame_hash.json``; refuses to run on a mismatch) and every
``results/axis_robustness_real/shards/*.instances.jsonl`` /
``*.cells.jsonl`` pair whose shard has a completion marker in
``results/axis_robustness_real/_done/``. A shard without a marker is ignored
entirely (the incrementality contract of the sweep driver). Each read file's
SHA-256 is verified against the marker's own ``instances_sha256`` /
``cells_sha256`` before it is trusted; a mismatch skips that one shard and is
recorded in ``analysis_summary.json``'s ``skipped_shards`` list rather than
aborting the run, because the sweep may still be writing.

Writes, all under the ``--out-dir`` (default ``results/axis_robustness_real``)
--------------------------------------------------------------------------
``analysis_units.csv``, ``analysis_tau.csv``, ``analysis_per_network.csv``,
``analysis_effect_sizes.csv``, ``analysis_contradiction.csv``,
``analysis_separation_strata.csv``, ``analysis_summary.json``,
``xarm_bins.csv``, ``xarm_paired_test.csv``, ``xarm_intensity_histogram.csv``.

Design is fixed by ``results/axis_robustness_real/PREREGISTRATION.md``,
sections 5 and 6 and Appendices A and B. This module imports its grid and
endpoint helpers (``auc_over_grid``, ``auc_usable``, ``depth_grid``,
``FRAC_GRID``, ``RATE_GRID``, ``MIN_USABLE_N``, ``intensity_bin``) from
:mod:`bkrobust.robustness.real_survival` rather than reimplementing them, so
the sweep's and the analysis's notion of an endpoint cannot drift apart.

Why the bootstrap is over networks, not instances
--------------------------------------------------
On this corpus a network's rows are not independent: the analyst's corrupted
knowledge state is shared across every pair of a network (the flip arm) or
generated once per network (the tiered arm), so two rows of one network at
one grid point are perfectly dependent draws of the same corruption. A
per-row (per-instance) bootstrap, as the synthetic counterpart
(:mod:`bkrobust.robustness.analyse`) uses, would treat those dependent rows
as independent evidence and understate every interval. Every bootstrap here
therefore resamples **networks** with replacement and takes *all* of a drawn
network's rows, exactly as ``PREREGISTRATION.md`` Section 6.2 specifies.

``all_valid_adjustment_sets_mpdag`` is never imported or called anywhere in
this module. No global RNG is ever touched; every resampling routine builds
its own ``numpy.random.Generator`` from an explicit seed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import binomtest, kendalltau, wilcoxon

from bkrobust.core.conventions import UNREACHED
from bkrobust.robustness.real_survival import (
    FRAC_GRID,
    MIN_USABLE_N,
    RATE_GRID,
    auc_over_grid,
    auc_usable,
)

DEFAULT_OUT_DIR = "results/axis_robustness_real"
DEFAULT_N_BOOT = 10_000
DEFAULT_SEED = 0

CENSORED_STATUS = "censored_wall_cap"

#: The four predictors under test (PREREGISTRATION.md Sec. 6.1).
PREDICTORS: tuple[str, ...] = ("radius", "shd_truth", "n_k", "k_g0")

#: Which pair of (raw, usable) endpoint columns a given arm's units carry.
ENDPOINTS_BY_ARM: dict[str, tuple[str, str]] = {
    "flip": ("AUC_frac", "AUC_frac_usable"),
    "tiered": ("AUC_rate", "AUC_rate_usable"),
}

#: The coverage / base-wrongness combinations that form the nine primary
#: strata's flip half (Sec. 3). The tiered half is simply n_tiers in {2,3,4}.
PRIMARY_FLIP_COVERAGES: tuple[float, ...] = (1.0, 0.5)
PRIMARY_FLIP_BASE_WRONGNESS: tuple[float, ...] = (0.0, 0.10, 0.25)


# ---------------------------------------------------------------------------
# Hashing / JSONL IO
# ---------------------------------------------------------------------------


def sha256_file(path: Path) -> str:
    """The SHA-256 hex digest of a file's raw bytes.

    Args:
        path: File to hash.

    Returns:
        The hex digest.
    """
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    """Reads a JSONL file into a list of dicts.

    Args:
        path: The file to read.

    Returns:
        One dict per non-blank line, in file order.
    """
    rows = []
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def verify_frame(out_dir: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Loads ``frame.jsonl`` and verifies its SHA-256 against ``frame_hash.json``.

    Args:
        out_dir: The results subtree holding both files.

    Returns:
        ``(frame_rows, frame_hash_info)``.

    Raises:
        FileNotFoundError: If either file is missing.
        RuntimeError: If the frame's digest does not match the recorded one.
            The caller must refuse to run on this, per the task brief.
    """
    fpath, hpath = out_dir / "frame.jsonl", out_dir / "frame_hash.json"
    if not fpath.is_file() or not hpath.is_file():
        raise FileNotFoundError(f"frame.jsonl / frame_hash.json missing under {out_dir}")
    info = json.loads(hpath.read_text())
    digest = sha256_file(fpath)
    if digest != info["frame_sha256"]:
        raise RuntimeError(
            f"frame.jsonl sha256 mismatch: computed {digest}, "
            f"recorded {info['frame_sha256']} -- refusing to run"
        )
    return read_jsonl(fpath), info


# ---------------------------------------------------------------------------
# Shard discovery and loading
# ---------------------------------------------------------------------------


def discover_shard_markers(out_dir: Path) -> list[dict[str, Any]]:
    """Lists every completed shard's marker payload, sorted by shard id.

    A shard with no marker in ``_done/`` does not exist for this analysis --
    that is the incrementality contract the sweep driver guarantees.

    Args:
        out_dir: The results subtree.

    Returns:
        Parsed marker payloads (each carries ``shard_id``, ``kind``, ``spec``,
        ``instances_path``, ``cells_path``, ``instances_sha256``,
        ``cells_sha256``), sorted by ``shard_id``.
    """
    done_dir = out_dir / "_done"
    if not done_dir.is_dir():
        return []
    markers = []
    for p in sorted(done_dir.glob("*.json")):
        markers.append(json.loads(p.read_text()))
    markers.sort(key=lambda m: m["shard_id"])
    return markers


def load_shard(
    out_dir: Path, marker: dict[str, Any]
) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
    """Loads one shard's instance and cell rows, hash-verified against its marker.

    Args:
        out_dir: The results subtree (used only if the marker's own recorded
            paths are unusable).
        marker: The shard's completion-marker payload.

    Returns:
        ``(loaded, skip_info)``. Exactly one is ``None``. ``loaded`` is
        ``{"kind", "marker", "instances", "cells"}``. ``skip_info`` is
        ``{"shard_id", "reason"}`` when the shard is skipped -- a missing
        file or a hash mismatch -- which must not abort the whole run since
        the sweep may still be writing.
    """
    sid = marker["shard_id"]
    inst_path = Path(marker.get("instances_path") or (out_dir / "shards" / f"{sid}.instances.jsonl"))
    cell_path = Path(marker.get("cells_path") or (out_dir / "shards" / f"{sid}.cells.jsonl"))
    if not inst_path.is_file() or not cell_path.is_file():
        return None, {"shard_id": sid, "reason": "missing_shard_files"}

    inst_digest = sha256_file(inst_path)
    cell_digest = sha256_file(cell_path)
    reasons = []
    if inst_digest != marker.get("instances_sha256"):
        reasons.append("instances_sha256_mismatch")
    if cell_digest != marker.get("cells_sha256"):
        reasons.append("cells_sha256_mismatch")
    if reasons:
        return None, {"shard_id": sid, "reason": ",".join(reasons)}

    instances = read_jsonl(inst_path)
    cells = read_jsonl(cell_path)
    return {"kind": marker["kind"], "marker": marker, "instances": instances, "cells": cells}, None


def load_all_shards(out_dir: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Discovers and loads every completed, hash-verified shard.

    Args:
        out_dir: The results subtree.

    Returns:
        ``(loaded_shards, skipped_shards)``.
    """
    loaded, skipped = [], []
    for marker in discover_shard_markers(out_dir):
        shard, skip_info = load_shard(out_dir, marker)
        if skip_info is not None:
            skipped.append(skip_info)
        else:
            loaded.append(shard)
    return loaded, skipped


# ---------------------------------------------------------------------------
# Stratum labelling
# ---------------------------------------------------------------------------


def flip_stratum_tag(coverage: float, base_wrongness: float | None, bw_abs: int | None) -> str:
    """The flip arm's stratum id, in the same ``covNNN``/``bwNNN``/``bwaN``
    tag style :func:`bkrobust.robustness.run_real_survival.shard_id` uses,
    so a stratum id can be grepped straight back to its shard ids.

    Args:
        coverage: The coverage level.
        base_wrongness: The rate-based base wrongness, or ``None`` if this is
            a ``bw_abs`` unit.
        bw_abs: The absolute base-wrongness level, or ``None`` if rate-based.

    Returns:
        The stratum id.
    """
    cov_tag = f"cov{round(coverage * 100):03d}"
    bw_tag = f"bwa{bw_abs}" if bw_abs is not None else f"bw{round(base_wrongness * 100):03d}"
    return f"flip_{cov_tag}_{bw_tag}"


def tiered_stratum_tag(n_tiers: int) -> str:
    """The tiered arm's stratum id.

    Args:
        n_tiers: Tier count.

    Returns:
        The stratum id.
    """
    return f"tiered_nt{n_tiers}"


def classify_stratum(
    arm: str, coverage: float | None, base_wrongness: float | None, bw_abs: int | None
) -> str:
    """Whether a unit's stratum is one of the nine primary strata (Sec. 3) or
    supplementary.

    Args:
        arm: ``"flip"`` or ``"tiered"``.
        coverage: The unit's coverage (flip only).
        base_wrongness: The unit's rate-based base wrongness (flip only).
        bw_abs: The unit's absolute base-wrongness level (flip only).

    Returns:
        ``"primary"`` or ``"supplementary"``.
    """
    if arm == "tiered":
        return "primary"
    if bw_abs is not None:
        return "supplementary"  # Appendix A's bw_abs=1 level is never primary.
    if coverage == 0.25:
        return "supplementary"  # Sec. 3: coverage 0.25 has no synthetic counterpart.
    if coverage in PRIMARY_FLIP_COVERAGES and base_wrongness in PRIMARY_FLIP_BASE_WRONGNESS:
        return "primary"
    return "supplementary"


# ---------------------------------------------------------------------------
# Analysis units (flip + tiered arms)
# ---------------------------------------------------------------------------

UNIT_COLUMNS: list[str] = [
    "arm", "stratum", "stratum_class", "status",
    "network", "x", "y",
    "coverage", "base_wrongness", "bw_abs", "analyst_replicate", "n_tiers",
    "radius", "shd_truth", "n_k", "k_g0",
    "AUC_frac", "AUC_frac_usable", "AUC_rate", "AUC_rate_usable", "endpoint_gap",
    "n_grid_points", "n_grid_points_usable", "n_grid_points_censored", "endpoint_censored",
    "min_n_eval", "median_n_eval",
    "separation", "separation_status", "largest_component_size", "component_size",
    "n_claims_actually_wrong", "bw_is_inert", "dispatch_leg", "gate", "assumes",
    "shard_id", "frame_row_id",
]


def _build_unit(
    irow: dict[str, Any], cells_for_unit: list[dict[str, Any]], endpoint_kind: str
) -> dict[str, Any]:
    """Builds one analysis-unit row from one instance row and its joined cells.

    Args:
        irow: The instance row (one frame pair within one shard).
        cells_for_unit: The cell rows sharing this unit's ``frame_row_id``,
            within the same shard.
        endpoint_kind: ``"AUC_frac"`` (flip) or ``"AUC_rate"`` (tiered).

    Returns:
        The unit row, with every :data:`UNIT_COLUMNS` key present.
    """
    censored = [c for c in cells_for_unit if c.get("status") == CENSORED_STATUS]
    noncensored = [c for c in cells_for_unit if c.get("status") != CENSORED_STATUS]
    cells_map = {c["grid_point"]: c for c in noncensored}

    n_eval_vals = [c["n_eval"] for c in noncensored if c.get("n_eval") is not None]
    min_n_eval = min(n_eval_vals) if n_eval_vals else None
    median_n_eval = statistics.median(n_eval_vals) if n_eval_vals else None
    n_grid_points_usable = sum(
        1 for c in noncensored
        if c.get("n_eval") is not None and c["n_eval"] >= MIN_USABLE_N and c.get("S") is not None
    )

    auc_frac = auc_frac_usable = auc_rate = auc_rate_usable = None
    if irow["status"] == "ok":
        if endpoint_kind == "AUC_frac":
            n_k = irow["n_k"]
            targets = [f * n_k for f in FRAC_GRID]
            auc_frac = auc_over_grid(cells_map, targets)
            auc_frac_usable = auc_usable(cells_map, targets)
        else:
            targets = list(RATE_GRID)
            auc_rate = auc_over_grid(cells_map, targets)
            auc_rate_usable = auc_usable(cells_map)

    raw = auc_frac if endpoint_kind == "AUC_frac" else auc_rate
    usable = auc_frac_usable if endpoint_kind == "AUC_frac" else auc_rate_usable
    endpoint_gap = abs(raw - usable) if (raw is not None and usable is not None) else None

    return {
        "arm": irow["arm"],
        "network": irow["network"], "x": irow["x"], "y": irow["y"],
        "coverage": irow.get("coverage"),
        "base_wrongness": irow.get("base_wrongness"),
        "bw_abs": irow.get("bw_abs_level"),
        "analyst_replicate": irow.get("analyst_replicate"),
        "n_tiers": irow.get("n_tiers"),
        "status": irow["status"],
        "radius": irow.get("radius"),
        "shd_truth": irow.get("shd_truth"),
        "n_k": irow.get("n_k"),
        "k_g0": irow.get("k_g0"),
        "AUC_frac": auc_frac, "AUC_frac_usable": auc_frac_usable,
        "AUC_rate": auc_rate, "AUC_rate_usable": auc_rate_usable,
        "endpoint_gap": endpoint_gap,
        "n_grid_points": len(cells_for_unit),
        "n_grid_points_usable": n_grid_points_usable,
        "n_grid_points_censored": len(censored),
        "endpoint_censored": len(censored) > 0,
        "min_n_eval": min_n_eval, "median_n_eval": median_n_eval,
        "separation": irow.get("separation"), "separation_status": irow.get("separation_status"),
        "largest_component_size": irow.get("largest_component_size"),
        "component_size": irow.get("component_size"),
        "n_claims_actually_wrong": irow.get("n_claims_actually_wrong"),
        "bw_is_inert": irow.get("bw_is_inert"),
        "dispatch_leg": irow.get("dispatch_leg"),
        "gate": irow.get("gate"),
        "assumes": irow.get("assumes"),
        "shard_id": irow.get("shard_id"),
        "frame_row_id": irow.get("frame_row_id"),
    }


def build_analysis_units(loaded_shards: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Builds every analysis unit from the loaded ``flip`` and ``tiered`` shards.

    An analysis unit is one frame pair within one shard: for the flip arm,
    ``(network, x, y, coverage, base_wrongness, bw_abs, analyst_replicate)``;
    for the tiered arm, ``(network, x, y, n_tiers)``. Units whose instance
    status is not ``"ok"`` are kept with every endpoint ``None`` -- they are
    never dropped, per the task brief. ``xarm`` shards are a separate
    instance population (Sec. 7) and never contribute units here.

    Args:
        loaded_shards: Output of :func:`load_all_shards` (loaded half only).

    Returns:
        Every unit row, each carrying ``stratum`` and ``stratum_class``.
    """
    units: list[dict[str, Any]] = []
    for shard in loaded_shards:
        kind = shard["kind"]
        if kind not in ("flip", "tiered"):
            continue
        endpoint_kind = "AUC_frac" if kind == "flip" else "AUC_rate"
        cells_by_frid: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for c in shard["cells"]:
            cells_by_frid[c["frame_row_id"]].append(c)
        for irow in shard["instances"]:
            unit = _build_unit(irow, cells_by_frid.get(irow["frame_row_id"], []), endpoint_kind)
            if kind == "flip":
                unit["stratum"] = flip_stratum_tag(unit["coverage"], unit["base_wrongness"], unit["bw_abs"])
            else:
                unit["stratum"] = tiered_stratum_tag(unit["n_tiers"])
            unit["stratum_class"] = classify_stratum(
                unit["arm"], unit["coverage"], unit["base_wrongness"], unit["bw_abs"]
            )
            units.append(unit)
    return units


def group_units_by_stratum(units: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    """Groups unit rows by their ``stratum`` label.

    Args:
        units: Output of :func:`build_analysis_units`.

    Returns:
        ``{stratum: [unit, ...]}``, insertion-ordered by first appearance.
    """
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for u in units:
        groups[u["stratum"]].append(u)
    return groups


# ---------------------------------------------------------------------------
# Statistics: tau_for_stratum, cluster bootstrap over networks, LOO
# ---------------------------------------------------------------------------


def _cluster_bootstrap_tau(
    x: np.ndarray, y: np.ndarray, net_of_row: list[str], *, n_boot: int, seed: int
) -> tuple[float | None, float | None, int]:
    """Cluster bootstrap over networks for a Kendall tau-b confidence interval.

    Each resample draws ``len(distinct networks)`` networks with replacement
    and concatenates *all* rows of each drawn network (a network drawn twice
    contributes its rows twice), then recomputes tau-b.

    Args:
        x: Predictor values, one per row.
        y: Endpoint values, one per row.
        net_of_row: The network each row belongs to, same length as ``x``.
        n_boot: Number of resamples.
        seed: Seed for ``np.random.default_rng``.

    Returns:
        ``(ci_lo_2p5, ci_hi_97p5, n_boot_nan)``. The CI bounds are ``None``
        if fewer than 100 resamples were valid.
    """
    networks = sorted(set(net_of_row))
    net_to_idx: dict[str, list[int]] = defaultdict(list)
    for i, net in enumerate(net_of_row):
        net_to_idx[net].append(i)
    net_arr = np.array(networks, dtype=object)
    n_nets = len(networks)

    rng = np.random.default_rng(seed)
    taus = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        drawn = net_arr[rng.integers(0, n_nets, size=n_nets)]
        idx: list[int] = []
        for net in drawn:
            idx.extend(net_to_idx[net])
        xb, yb = x[idx], y[idx]
        if np.unique(xb).size <= 1 or np.unique(yb).size <= 1:
            taus[b] = np.nan
            continue
        t, _ = kendalltau(xb, yb, variant="b", nan_policy="propagate")
        taus[b] = t

    n_nan = int(np.isnan(taus).sum())
    valid = taus[~np.isnan(taus)]
    if len(valid) < 100:
        return None, None, n_nan
    lo, hi = np.percentile(valid, [2.5, 97.5])
    return float(lo), float(hi), n_nan


def _leave_one_network_out(
    x: np.ndarray, y: np.ndarray, net_of_row: list[str]
) -> dict[str, float]:
    """Point-estimate tau-b with each network's rows removed in turn.

    Args:
        x: Predictor values.
        y: Endpoint values.
        net_of_row: The network each row belongs to.

    Returns:
        ``{network: tau_b}`` for networks whose removal still leaves a
        computable tau (>=2 rows, non-constant predictor and endpoint).
    """
    networks = sorted(set(net_of_row))
    out: dict[str, float] = {}
    net_arr = np.array(net_of_row)
    for net in networks:
        keep = net_arr != net
        xo, yo = x[keep], y[keep]
        if len(xo) < 2 or np.unique(xo).size <= 1 or np.unique(yo).size <= 1:
            continue
        t, _ = kendalltau(xo, yo, variant="b", nan_policy="propagate")
        if not (isinstance(t, float) and math.isnan(t)):
            out[net] = float(t)
    return out


def tau_for_stratum(
    rows: list[dict[str, Any]], predictor: str, endpoint: str,
    *, n_boot: int = DEFAULT_N_BOOT, seed: int = DEFAULT_SEED,
) -> dict[str, Any]:
    """One (stratum, predictor, endpoint) row of ``analysis_tau.csv``.

    Follows ``PREREGISTRATION.md`` Sec. 6.2 exactly: Kendall tau-b, a cluster
    bootstrap over networks, and a leave-one-network-out sensitivity check.

    Args:
        rows: Every unit row of one stratum (any status; endpoint/predictor
            ``None`` rows are dropped here, not upstream).
        predictor: One of :data:`PREDICTORS`.
        endpoint: The endpoint column name for this stratum's arm.
        n_boot: Bootstrap resamples.
        seed: Seed for every ``np.random.default_rng`` this call constructs.

    Returns:
        A dict with ``predictor``, ``endpoint``, ``status``, ``n``,
        ``n_networks``, ``stratum_n_units``, ``n_excluded_unreached``,
        ``n_excluded_nan``, ``tau_b``, ``p_value``, ``ci_lo_2p5``,
        ``ci_hi_97p5``, ``n_boot_nan``, the LOO fields, and ``r_assumes``.
    """
    stratum_n_units = len(rows)
    cleaned: list[dict[str, Any]] = []
    n_excluded_nan = 0
    n_excluded_unreached = 0
    for r in rows:
        pv, ev = r.get(predictor), r.get(endpoint)
        if pv is None or ev is None:
            n_excluded_nan += 1
            continue
        if predictor == "radius" and pv == UNREACHED:
            n_excluded_unreached += 1
            continue
        cleaned.append(r)

    out: dict[str, Any] = {
        "predictor": predictor, "endpoint": endpoint,
        "n": len(cleaned), "n_networks": len(set(r["network"] for r in cleaned)),
        "stratum_n_units": stratum_n_units,
        "n_excluded_unreached": n_excluded_unreached, "n_excluded_nan": n_excluded_nan,
        "tau_b": None, "p_value": None, "ci_lo_2p5": None, "ci_hi_97p5": None, "n_boot_nan": None,
        "loo_tau_min": None, "loo_tau_max": None,
        "loo_network_at_min": None, "loo_network_at_max": None, "loo_verdict_flips": None,
        "r_assumes": "", "status": "ok",
    }

    if predictor == "radius" and cleaned:
        assumes_vals = sorted({r["assumes"] for r in cleaned if r.get("assumes")})
        if len(assumes_vals) > 1:
            raise AssertionError(
                f"r_assumes is not constant within stratum for predictor=radius: {assumes_vals}"
            )
        out["r_assumes"] = assumes_vals[0] if assumes_vals else ""

    if out["n"] < 2:
        out["status"] = "insufficient_n"
        return out

    x = np.array([float(r[predictor]) for r in cleaned])
    y = np.array([float(r[endpoint]) for r in cleaned])
    net_of_row = [r["network"] for r in cleaned]

    if np.unique(x).size <= 1:
        out["status"] = "undefined (predictor constant)"
        return out
    if np.unique(y).size <= 1:
        out["status"] = "undefined (endpoint constant)"
        return out

    tau, pval = kendalltau(x, y, variant="b", nan_policy="propagate")
    out["tau_b"] = float(tau)
    out["p_value"] = float(pval)

    lo, hi, n_boot_nan = _cluster_bootstrap_tau(x, y, net_of_row, n_boot=n_boot, seed=seed)
    out["n_boot_nan"] = n_boot_nan
    if lo is None:
        out["status"] = "bootstrap_degenerate"
    else:
        out["ci_lo_2p5"], out["ci_hi_97p5"] = lo, hi

    loo = _leave_one_network_out(x, y, net_of_row)
    if loo:
        net_min = min(loo, key=lambda k: loo[k])
        net_max = max(loo, key=lambda k: loo[k])
        out["loo_tau_min"], out["loo_network_at_min"] = loo[net_min], net_min
        out["loo_tau_max"], out["loo_network_at_max"] = loo[net_max], net_max
        full_sign = 0.0 if out["tau_b"] == 0 else math.copysign(1.0, out["tau_b"])
        opposite_sign = any(
            v != 0 and full_sign != 0 and math.copysign(1.0, v) != full_sign for v in loo.values()
        )
        ci_excludes_zero = out["ci_lo_2p5"] is not None and (
            out["ci_lo_2p5"] > 0 or out["ci_hi_97p5"] < 0
        )
        flips_by_ci = ci_excludes_zero and (loo[net_min] <= 0 <= loo[net_max])
        out["loo_verdict_flips"] = bool(opposite_sign or flips_by_ci)

    return out


def stratum_weighted_medians(
    units: list[dict[str, Any]], endpoint: str
) -> tuple[float | None, float | None]:
    """Pair-weighted and network-weighted median of one endpoint in one stratum.

    Pair-weighted is the plain median over every unit with a defined value.
    Network-weighted is each network's own median, then the mean over
    networks (never the median of medians).

    Args:
        units: The stratum's unit rows.
        endpoint: The endpoint column to summarise.

    Returns:
        ``(pair_weighted_median, network_weighted_median)``, either ``None``
        if no unit has a defined value.
    """
    vals = [u[endpoint] for u in units if u[endpoint] is not None]
    pair_weighted = statistics.median(vals) if vals else None

    by_net: dict[str, list[float]] = defaultdict(list)
    for u in units:
        if u[endpoint] is not None:
            by_net[u["network"]].append(u[endpoint])
    net_medians = [statistics.median(v) for v in by_net.values()]
    network_weighted = (sum(net_medians) / len(net_medians)) if net_medians else None
    return pair_weighted, network_weighted


TAU_COLUMNS: list[str] = [
    "stratum", "stratum_class", "arm", "predictor", "endpoint", "status",
    "n", "n_networks", "stratum_n_units", "n_excluded_unreached", "n_excluded_nan",
    "tau_b", "p_value", "ci_lo_2p5", "ci_hi_97p5", "n_boot_nan",
    "loo_tau_min", "loo_tau_max", "loo_network_at_min", "loo_network_at_max", "loo_verdict_flips",
    "r_assumes", "median_endpoint_pair_weighted", "median_endpoint_network_weighted",
]


def build_tau_table(
    stratum_groups: dict[str, list[dict[str, Any]]],
    *, n_boot: int = DEFAULT_N_BOOT, seed: int = DEFAULT_SEED,
) -> list[dict[str, Any]]:
    """Builds every row of ``analysis_tau.csv``: one per (stratum, predictor,
    endpoint), where a stratum's endpoint pair is determined by its arm.

    Args:
        stratum_groups: Output of :func:`group_units_by_stratum`.
        n_boot: Bootstrap resamples per cell.
        seed: Seed for every bootstrap.

    Returns:
        Rows in stratum, then predictor, then endpoint order.
    """
    rows = []
    for stratum in sorted(stratum_groups):
        units = stratum_groups[stratum]
        arm = units[0]["arm"]
        stratum_class = units[0]["stratum_class"]
        endpoints = ENDPOINTS_BY_ARM[arm]
        medians = {ep: stratum_weighted_medians(units, ep) for ep in endpoints}
        for predictor in PREDICTORS:
            for endpoint in endpoints:
                r = tau_for_stratum(units, predictor, endpoint, n_boot=n_boot, seed=seed)
                r["stratum"], r["stratum_class"], r["arm"] = stratum, stratum_class, arm
                pw, nw = medians[endpoint]
                r["median_endpoint_pair_weighted"] = pw
                r["median_endpoint_network_weighted"] = nw
                rows.append(r)
    return rows


# ---------------------------------------------------------------------------
# Per-network, effect-size and separation tables
# ---------------------------------------------------------------------------


def _median_iqr(values: list[float]) -> tuple[float | None, float | None, float | None]:
    """Median, Q1 and Q3 of a list of floats.

    Args:
        values: The values. May be empty.

    Returns:
        ``(median, q1, q3)``, all ``None`` if ``values`` is empty.
    """
    if not values:
        return None, None, None
    arr = np.asarray(sorted(values), dtype=float)
    med, q1, q3 = (float(v) for v in np.percentile(arr, [50, 25, 75]))
    return med, q1, q3


def _share_radius_1(radii: list[int]) -> float | None:
    return (sum(1 for r in radii if r == 1) / len(radii)) if radii else None


PER_NETWORK_COLUMNS: list[str] = [
    "stratum", "stratum_class", "arm", "network", "n_units",
    "endpoint_raw_name", "median_endpoint_raw", "q1_endpoint_raw", "q3_endpoint_raw",
    "endpoint_usable_name", "median_endpoint_usable", "q1_endpoint_usable", "q3_endpoint_usable",
    "n_distinct_radius", "share_radius_1", "median_n_k", "median_shd_truth",
]


def build_per_network_table(
    stratum_groups: dict[str, list[dict[str, Any]]]
) -> list[dict[str, Any]]:
    """Builds ``analysis_per_network.csv``: one row per (stratum, network).

    Args:
        stratum_groups: Output of :func:`group_units_by_stratum`.

    Returns:
        Rows in stratum, then network order.
    """
    rows = []
    for stratum in sorted(stratum_groups):
        units = stratum_groups[stratum]
        arm, stratum_class = units[0]["arm"], units[0]["stratum_class"]
        ep_raw, ep_usable = ENDPOINTS_BY_ARM[arm]
        by_net: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for u in units:
            by_net[u["network"]].append(u)
        for net in sorted(by_net):
            us = by_net[net]
            radii = [u["radius"] for u in us if u["radius"] is not None]
            med_r, q1_r, q3_r = _median_iqr([u[ep_raw] for u in us if u[ep_raw] is not None])
            med_u, q1_u, q3_u = _median_iqr([u[ep_usable] for u in us if u[ep_usable] is not None])
            n_k_vals = [u["n_k"] for u in us if u["n_k"] is not None]
            shd_vals = [u["shd_truth"] for u in us if u["shd_truth"] is not None]
            rows.append({
                "stratum": stratum, "stratum_class": stratum_class, "arm": arm, "network": net,
                "n_units": len(us),
                "endpoint_raw_name": ep_raw, "median_endpoint_raw": med_r,
                "q1_endpoint_raw": q1_r, "q3_endpoint_raw": q3_r,
                "endpoint_usable_name": ep_usable, "median_endpoint_usable": med_u,
                "q1_endpoint_usable": q1_u, "q3_endpoint_usable": q3_u,
                "n_distinct_radius": len(set(radii)),
                "share_radius_1": _share_radius_1(radii),
                "median_n_k": statistics.median(n_k_vals) if n_k_vals else None,
                "median_shd_truth": statistics.median(shd_vals) if shd_vals else None,
            })
    return rows


RADIUS_BUCKETS_ORDER = {"1": 1, "2": 2, "3": 3, "4": 4, "5+": 5, "UNREACHED": 6}


def radius_bucket(r: int | None) -> str | None:
    """The effect-size radius bucket a radius value falls into.

    Args:
        r: A radius value, ``UNREACHED`` (-1), or ``None``.

    Returns:
        One of ``"1"``, ``"2"``, ``"3"``, ``"4"``, ``"5+"``, ``"UNREACHED"``,
        or ``None`` if ``r`` is ``None`` (radius never computed -- excluded
        from the effect-size table, distinct from ``UNREACHED``).
    """
    if r is None:
        return None
    if r == UNREACHED:
        return "UNREACHED"
    return f"{r}" if r < 5 else "5+"


EFFECT_SIZE_COLUMNS: list[str] = [
    "stratum", "stratum_class", "arm", "radius_bucket", "n_units", "n_networks",
    "endpoint_raw_name", "median_endpoint_raw", "q1_endpoint_raw", "q3_endpoint_raw",
    "endpoint_usable_name", "median_endpoint_usable", "q1_endpoint_usable", "q3_endpoint_usable",
]


def build_effect_size_table(
    stratum_groups: dict[str, list[dict[str, Any]]]
) -> list[dict[str, Any]]:
    """Builds ``analysis_effect_sizes.csv``: one row per (stratum, radius bucket).

    ``UNREACHED`` is its own bucket, per Sec. 6.2, and is never placed on the
    numeric radius axis; a unit whose radius was never computed (``None``) is
    excluded from this table entirely, since it is a different fact from
    ``UNREACHED``.

    Args:
        stratum_groups: Output of :func:`group_units_by_stratum`.

    Returns:
        Rows in stratum, then bucket order (1, 2, 3, 4, 5+, UNREACHED).
    """
    rows = []
    for stratum in sorted(stratum_groups):
        units = stratum_groups[stratum]
        arm, stratum_class = units[0]["arm"], units[0]["stratum_class"]
        ep_raw, ep_usable = ENDPOINTS_BY_ARM[arm]
        by_bucket: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for u in units:
            b = radius_bucket(u["radius"])
            if b is not None:
                by_bucket[b].append(u)
        for bucket in sorted(by_bucket, key=lambda b: RADIUS_BUCKETS_ORDER[b]):
            us = by_bucket[bucket]
            med_r, q1_r, q3_r = _median_iqr([u[ep_raw] for u in us if u[ep_raw] is not None])
            med_u, q1_u, q3_u = _median_iqr([u[ep_usable] for u in us if u[ep_usable] is not None])
            rows.append({
                "stratum": stratum, "stratum_class": stratum_class, "arm": arm,
                "radius_bucket": bucket, "n_units": len(us),
                "n_networks": len({u["network"] for u in us}),
                "endpoint_raw_name": ep_raw, "median_endpoint_raw": med_r,
                "q1_endpoint_raw": q1_r, "q3_endpoint_raw": q3_r,
                "endpoint_usable_name": ep_usable, "median_endpoint_usable": med_u,
                "q1_endpoint_usable": q1_u, "q3_endpoint_usable": q3_u,
            })
    return rows


SEPARATION_COLUMNS: list[str] = [
    "stratum", "stratum_class", "arm", "separation_status", "n_units", "n_networks",
    "endpoint_raw_name", "median_endpoint_raw", "endpoint_usable_name", "median_endpoint_usable",
    "share_radius_1",
]


def build_separation_table(
    stratum_groups: dict[str, list[dict[str, Any]]]
) -> list[dict[str, Any]]:
    """Builds ``analysis_separation_strata.csv``: one row per (stratum,
    separation_status).

    Args:
        stratum_groups: Output of :func:`group_units_by_stratum`.

    Returns:
        Rows in stratum, then separation_status order.
    """
    rows = []
    for stratum in sorted(stratum_groups):
        units = stratum_groups[stratum]
        arm, stratum_class = units[0]["arm"], units[0]["stratum_class"]
        ep_raw, ep_usable = ENDPOINTS_BY_ARM[arm]
        by_sep: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for u in units:
            by_sep[u["separation_status"] or "NA"].append(u)
        for sep_status in sorted(by_sep):
            us = by_sep[sep_status]
            radii = [u["radius"] for u in us if u["radius"] is not None]
            med_r, *_ = _median_iqr([u[ep_raw] for u in us if u[ep_raw] is not None])
            med_u, *_ = _median_iqr([u[ep_usable] for u in us if u[ep_usable] is not None])
            rows.append({
                "stratum": stratum, "stratum_class": stratum_class, "arm": arm,
                "separation_status": sep_status, "n_units": len(us),
                "n_networks": len({u["network"] for u in us}),
                "endpoint_raw_name": ep_raw, "median_endpoint_raw": med_r,
                "endpoint_usable_name": ep_usable, "median_endpoint_usable": med_u,
                "share_radius_1": _share_radius_1(radii),
            })
    return rows


# ---------------------------------------------------------------------------
# Contradiction table (flip + tiered arms, pooled over units, not averaged)
# ---------------------------------------------------------------------------

CONTRADICTION_COLUMNS: list[str] = [
    "arm", "stratum", "grid_point", "n_cells", "n_draws", "n_contradictory", "contradiction_rate",
]


def build_contradiction_table(loaded_shards: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Builds ``analysis_contradiction.csv``.

    Pooled per (arm, stratum, grid_point): the contradiction rate is the sum
    of contradictory draws over the sum of all draws across every cell in
    that group -- never the mean of the per-cell rates.

    Args:
        loaded_shards: Loaded ``flip`` and ``tiered`` shards (``xarm`` shards
            are ignored; they are a separate instance population).

    Returns:
        Rows in (arm, stratum, grid_point) order.
    """
    groups: dict[tuple[str, str, Any], dict[str, int]] = defaultdict(
        lambda: {"n_contra": 0, "n_draws": 0, "n_cells": 0}
    )
    for shard in loaded_shards:
        if shard["kind"] not in ("flip", "tiered"):
            continue
        arm = shard["kind"]
        for c in shard["cells"]:
            if c.get("status") == CENSORED_STATUS:
                continue
            if arm == "flip":
                stratum = flip_stratum_tag(c["coverage"], c.get("base_wrongness"), c.get("bw_abs"))
            else:
                stratum = tiered_stratum_tag(c["n_tiers"])
            key = (arm, stratum, c["grid_point"])
            g = groups[key]
            g["n_contra"] += c.get("n_contradictory") or 0
            g["n_draws"] += c.get("n_draws") or 0
            g["n_cells"] += 1

    rows = []
    for (arm, stratum, gp), g in groups.items():
        rate = (g["n_contra"] / g["n_draws"]) if g["n_draws"] else None
        rows.append({
            "arm": arm, "stratum": stratum, "grid_point": gp,
            "n_cells": g["n_cells"], "n_draws": g["n_draws"],
            "n_contradictory": g["n_contra"], "contradiction_rate": rate,
        })
    return sorted(rows, key=lambda r: (r["arm"], r["stratum"], r["grid_point"]))


# ---------------------------------------------------------------------------
# Cross-arm: xarm_bins, xarm_paired_test, xarm_intensity_histogram
# ---------------------------------------------------------------------------

XARM_BINS_COLUMNS: list[str] = [
    "network", "x", "y", "n_tiers", "arm", "intensity_bin",
    "n_draws", "n_eval", "n_contradictory", "n_survived", "S",
    "contradiction_rate", "silent_failure_rate", "instance_sha256",
]


def build_xarm_bins(loaded_shards: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Builds ``xarm_bins.csv``.

    Pools every cell of one instance that lands in the same intensity bin --
    summing draw counts and recomputing rates from the sums -- rather than
    averaging per-cell rates, per the task brief.

    Args:
        loaded_shards: All loaded shards (only ``xarm`` ones contribute).

    Returns:
        Rows in (network, x, y, n_tiers, arm, intensity_bin) order.
    """
    groups: dict[tuple, dict[str, Any]] = {}
    for shard in loaded_shards:
        if shard["kind"] != "xarm":
            continue
        for c in shard["cells"]:
            if c.get("status") != "ok" or c.get("intensity_bin") is None:
                continue
            key = (c["network"], c["x"], c["y"], c["n_tiers"], c["arm"], c["intensity_bin"])
            g = groups.setdefault(
                key, {"n_draws": 0, "n_eval": 0, "n_contradictory": 0, "n_survived": 0,
                      "instance_sha256": c.get("instance_sha256")}
            )
            g["n_draws"] += c["n_draws"]
            g["n_eval"] += c["n_eval"]
            g["n_contradictory"] += c["n_contradictory"]
            g["n_survived"] += c["n_survived"]

    rows = []
    for (net, x, y, nt, arm, ibin), g in groups.items():
        n_draws, n_eval, n_contra, n_surv = g["n_draws"], g["n_eval"], g["n_contradictory"], g["n_survived"]
        S = (n_surv / n_eval) if n_eval else None
        contradiction_rate = (n_contra / n_draws) if n_draws else None
        silent_failure_rate = (
            (1.0 - contradiction_rate) * (1.0 - S) if (S is not None and contradiction_rate is not None) else None
        )
        rows.append({
            "network": net, "x": x, "y": y, "n_tiers": nt, "arm": arm, "intensity_bin": ibin,
            "n_draws": n_draws, "n_eval": n_eval, "n_contradictory": n_contra, "n_survived": n_surv,
            "S": S, "contradiction_rate": contradiction_rate, "silent_failure_rate": silent_failure_rate,
            "instance_sha256": g["instance_sha256"],
        })
    return sorted(rows, key=lambda r: (r["network"], r["x"], r["y"], r["n_tiers"], r["arm"], r["intensity_bin"]))


def _instance_sha_map(loaded_shards: list[dict[str, Any]], xarm_kind: str) -> dict[tuple, str | None]:
    """Maps ``(network, x, y, n_tiers)`` -> ``instance_sha256`` for one cross-arm process.

    Args:
        loaded_shards: All loaded shards.
        xarm_kind: ``"flip"`` or ``"tiered"`` -- the corruption process.

    Returns:
        The map, built from cell rows of the matching ``xarm`` shards.
    """
    m: dict[tuple, str | None] = {}
    for shard in loaded_shards:
        if shard["kind"] != "xarm":
            continue
        spec = shard["marker"].get("spec", {})
        if spec.get("xarm") != xarm_kind:
            continue
        for c in shard["cells"]:
            key = (c["network"], c["x"], c["y"], c["n_tiers"])
            m.setdefault(key, c.get("instance_sha256"))
    return m


def _cluster_bootstrap_mean(
    values: np.ndarray, networks: list[str], *, n_boot: int, seed: int
) -> tuple[float | None, float | None, float | None, int]:
    """Cluster bootstrap over networks for the mean of a paired difference.

    Args:
        values: One paired-difference value per matched instance.
        networks: The network each matched instance belongs to.
        n_boot: Resamples.
        seed: Seed for ``np.random.default_rng``.

    Returns:
        ``(boot_mean, ci_lo_2p5, ci_hi_97p5, n_boot_nan)``. The mean and CI
        are ``None`` if fewer than 100 resamples were valid.
    """
    if len(values) == 0:
        return None, None, None, 0
    uniq_nets = sorted(set(networks))
    net_to_idx: dict[str, list[int]] = defaultdict(list)
    for i, net in enumerate(networks):
        net_to_idx[net].append(i)
    net_arr = np.array(uniq_nets, dtype=object)
    n_nets = len(uniq_nets)

    rng = np.random.default_rng(seed)
    means = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        drawn = net_arr[rng.integers(0, n_nets, size=n_nets)]
        idx: list[int] = []
        for net in drawn:
            idx.extend(net_to_idx[net])
        means[b] = values[idx].mean() if idx else np.nan

    n_nan = int(np.isnan(means).sum())
    valid = means[~np.isnan(means)]
    if len(valid) < 100:
        return None, None, None, n_nan
    lo, hi = np.percentile(valid, [2.5, 97.5])
    return float(valid.mean()), float(lo), float(hi), n_nan


def _paired_stats(
    diffs: list[float], networks: list[str], *, n_boot: int, seed: int
) -> dict[str, Any]:
    """Sign test, Wilcoxon and cluster bootstrap for one set of paired differences.

    Args:
        diffs: ``tiered - flip`` differences, one per matched instance.
        networks: The network of each matched instance, same length as ``diffs``.
        n_boot: Bootstrap resamples.
        seed: Seed for the bootstrap's ``np.random.default_rng``.

    Returns:
        A flat dict of the paired-test columns for one metric.
    """
    arr = np.asarray(diffs, dtype=float)
    n_pos = int((arr > 0).sum())
    n_neg = int((arr < 0).sum())
    n_zero = int((arr == 0).sum())
    nonzero = arr[arr != 0]

    sign_p = None
    if len(nonzero) > 0:
        k = int((nonzero > 0).sum())
        sign_p = float(binomtest(k, len(nonzero), 0.5, alternative="two-sided").pvalue)

    wil_p = None
    if len(nonzero) >= 6:
        try:
            wil_p = float(wilcoxon(nonzero).pvalue)
        except ValueError:
            wil_p = None

    boot_mean, ci_lo, ci_hi, boot_nan = _cluster_bootstrap_mean(arr, networks, n_boot=n_boot, seed=seed)

    return {
        "mean_diff_tiered_minus_flip": float(arr.mean()) if len(arr) else None,
        "median_diff": float(np.median(arr)) if len(arr) else None,
        "n_pos": n_pos, "n_neg": n_neg, "n_zero": n_zero,
        "sign_test_p": sign_p, "wilcoxon_p": wil_p,
        "boot_mean": boot_mean, "boot_ci_lo_2p5": ci_lo, "boot_ci_hi_97p5": ci_hi, "boot_n_nan": boot_nan,
    }


_PAIRED_METRICS = ("S", "contradiction_rate", "silent_failure_rate")

XARM_PAIRED_COLUMNS: list[str] = ["intensity_bin", "n_instances_matched"] + [
    f"{metric}_{field}"
    for metric in _PAIRED_METRICS
    for field in (
        "mean_diff_tiered_minus_flip", "median_diff", "n_pos", "n_neg", "n_zero",
        "sign_test_p", "wilcoxon_p", "boot_mean", "boot_ci_lo_2p5", "boot_ci_hi_97p5", "boot_n_nan",
    )
]


def paired_cross_arm(
    loaded_shards: list[dict[str, Any]], xarm_bins: list[dict[str, Any]],
    *, n_boot: int = DEFAULT_N_BOOT, seed: int = DEFAULT_SEED,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Builds ``xarm_paired_test.csv``, per Sec. 7 of the pre-registration.

    Before pairing, asserts ``instance_sha256`` is identical between the two
    arms for each ``(network, x, y, n_tiers)`` instance; an instance whose
    two shards disagree is refused and counted, never paired.

    Args:
        loaded_shards: All loaded shards.
        xarm_bins: Output of :func:`build_xarm_bins`.
        n_boot: Bootstrap resamples per bin per metric.
        seed: Seed for every bootstrap this call performs.

    Returns:
        ``(rows, extra_summary)`` where ``extra_summary`` carries the
        instance_sha256 pairing-refusal count for ``analysis_summary.json``.
    """
    flip_sha = _instance_sha_map(loaded_shards, "flip")
    tiered_sha = _instance_sha_map(loaded_shards, "tiered")
    common = set(flip_sha) & set(tiered_sha)
    mismatched = {k for k in common if flip_sha[k] != tiered_sha[k]}
    pairable = common - mismatched

    by_key_bin: dict[tuple, dict[str, dict[str, Any]]] = defaultdict(dict)
    for row in xarm_bins:
        inst_key = (row["network"], row["x"], row["y"], row["n_tiers"])
        sub_arm = "flip" if row["arm"] == "xarm_flip" else "tiered"
        by_key_bin[(inst_key, row["intensity_bin"])][sub_arm] = row

    by_bin: dict[float, list[tuple]] = defaultdict(list)
    for (inst_key, ibin), arms in by_key_bin.items():
        if inst_key not in pairable:
            continue
        if "flip" not in arms or "tiered" not in arms:
            continue
        f, t = arms["flip"], arms["tiered"]
        if f["n_eval"] < MIN_USABLE_N or t["n_eval"] < MIN_USABLE_N:
            continue
        by_bin[ibin].append((inst_key, f, t))

    rows = []
    for ibin in sorted(by_bin):
        matched = by_bin[ibin]
        networks = [k[0] for k, _, _ in matched]
        row: dict[str, Any] = {"intensity_bin": ibin, "n_instances_matched": len(matched)}
        for metric in _PAIRED_METRICS:
            diffs = [t[metric] - f[metric] for _, f, t in matched]
            stats = _paired_stats(diffs, networks, n_boot=n_boot, seed=seed)
            for field, val in stats.items():
                row[f"{metric}_{field}"] = val
        rows.append(row)

    extra = {
        "n_instances_common_to_both_xarm_processes": len(common),
        "n_instances_instance_sha256_mismatch": len(mismatched),
        "n_instances_pairable": len(pairable),
        "n_intensity_bins_with_matched_data": len(by_bin),
    }
    return rows, extra


XARM_HIST_COLUMNS: list[str] = ["arm", "intensity_bin", "n_cells", "n_draws"]


def build_xarm_intensity_histogram(loaded_shards: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Builds ``xarm_intensity_histogram.csv``: achieved-intensity coverage per arm.

    Args:
        loaded_shards: All loaded shards (only ``xarm`` ones contribute).

    Returns:
        Rows in (arm, intensity_bin) order.
    """
    groups: dict[tuple[str, float], dict[str, int]] = defaultdict(lambda: {"n_cells": 0, "n_draws": 0})
    for shard in loaded_shards:
        if shard["kind"] != "xarm":
            continue
        for c in shard["cells"]:
            if c.get("status") != "ok" or c.get("intensity_bin") is None:
                continue
            g = groups[(c["arm"], c["intensity_bin"])]
            g["n_cells"] += 1
            g["n_draws"] += c["n_draws"]
    rows = [
        {"arm": a, "intensity_bin": b, "n_cells": g["n_cells"], "n_draws": g["n_draws"]}
        for (a, b), g in groups.items()
    ]
    return sorted(rows, key=lambda r: (r["arm"], r["intensity_bin"]))


# ---------------------------------------------------------------------------
# Summary checks
# ---------------------------------------------------------------------------


def check_decomposition_identity(loaded_shards: list[dict[str, Any]]) -> dict[str, Any]:
    """Verifies ``S_contra_as_fail = (1 - contradiction_rate) * S`` on every
    non-censored cell of every loaded shard (flip, tiered and xarm alike).

    Args:
        loaded_shards: All loaded shards.

    Returns:
        ``{"n_checked", "max_abs_deviation"}``.
    """
    max_dev = 0.0
    n_checked = 0
    for shard in loaded_shards:
        for c in shard["cells"]:
            if c.get("status") == CENSORED_STATUS:
                continue
            S, cr, scaf = c.get("S"), c.get("contradiction_rate"), c.get("S_contra_as_fail")
            if S is None or cr is None or scaf is None:
                continue
            dev = abs((1.0 - cr) * S - scaf)
            max_dev = max(max_dev, dev)
            n_checked += 1
    return {"n_checked": n_checked, "max_abs_deviation": max_dev}


def check_n_draws(loaded_shards: list[dict[str, Any]]) -> dict[str, Any]:
    """Checks every non-censored cell's ``n_draws`` equals 1000.

    Args:
        loaded_shards: All loaded shards.

    Returns:
        ``{"n_checked", "all_equal_1000", "n_exceptions", "examples"}``.
    """
    exceptions = []
    n_checked = 0
    for shard in loaded_shards:
        sid = shard["marker"]["shard_id"]
        for c in shard["cells"]:
            if c.get("status") == CENSORED_STATUS:
                continue
            n_checked += 1
            if c.get("n_draws") != 1000:
                exceptions.append({"shard_id": sid, "grid_point": c.get("grid_point"), "n_draws": c.get("n_draws")})
    return {
        "n_checked": n_checked, "all_equal_1000": len(exceptions) == 0,
        "n_exceptions": len(exceptions), "examples": exceptions[:50],
    }


def check_censored_cells(loaded_shards: list[dict[str, Any]]) -> dict[str, Any]:
    """Counts censored cells and the shards they appear in.

    Args:
        loaded_shards: All loaded shards.

    Returns:
        ``{"n_censored_cells", "shards"}``.
    """
    n = 0
    shards: set[str] = set()
    for shard in loaded_shards:
        sid = shard["marker"]["shard_id"]
        for c in shard["cells"]:
            if c.get("status") == CENSORED_STATUS:
                n += 1
                shards.add(sid)
    return {"n_censored_cells": n, "shards": sorted(shards)}


def check_tiered_rate0(loaded_shards: list[dict[str, Any]]) -> dict[str, Any]:
    """Falsification trigger T4: ``S == 1.000`` exactly at rate 0, for every
    tiered-arm unit individually.

    Args:
        loaded_shards: All loaded shards (only plain ``tiered`` shards checked
            -- ``xarm`` tiered-process shards never sweep rate 0).

    Returns:
        ``{"n_checked", "n_violations", "examples"}``.
    """
    violations = []
    n_checked = 0
    for shard in loaded_shards:
        if shard["kind"] != "tiered":
            continue
        sid = shard["marker"]["shard_id"]
        for c in shard["cells"]:
            if c.get("status") == CENSORED_STATUS or c.get("grid_point") != 0.0:
                continue
            n_checked += 1
            if c.get("S") != 1.0:
                violations.append({"shard_id": sid, "frame_row_id": c.get("frame_row_id"), "S": c.get("S")})
    return {"n_checked": n_checked, "n_violations": len(violations), "examples": violations[:10]}


# ---------------------------------------------------------------------------
# CSV / JSON writers
# ---------------------------------------------------------------------------


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    """Writes rows to a CSV with a stable header via ``csv.DictWriter``.

    Args:
        path: Destination file.
        fieldnames: The stable column order. Every row carries every column
            (missing keys become blank, NaN floats become blank).
        rows: The rows to write.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for r in rows:
            clean = {}
            for k in fieldnames:
                v = r.get(k)
                if isinstance(v, float) and math.isnan(v):
                    v = None
                clean[k] = v
            writer.writerow(clean)


def write_json(path: Path, obj: Any) -> None:
    """Writes ``obj`` as indented JSON.

    Args:
        path: Destination file.
        obj: A JSON-serialisable object.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, default=str))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    """The command-line interface.

    Returns:
        The parser.
    """
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    p.add_argument("--n-boot", type=int, default=DEFAULT_N_BOOT)
    return p


def main(argv: list[str] | None = None) -> int:
    """Entry point: reads the frame and every completed shard, writes every
    output file under ``--out-dir``.

    Args:
        argv: Argument vector, or ``None`` for ``sys.argv``.

    Returns:
        Process exit status. Nonzero if the frame hash does not verify.
    """
    args = build_arg_parser().parse_args(argv)
    out_dir = Path(args.out_dir)
    n_boot = args.n_boot
    seed = DEFAULT_SEED

    try:
        frame_rows, frame_info = verify_frame(out_dir)
    except (FileNotFoundError, RuntimeError) as exc:
        print(f"REFUSING TO RUN: {exc}")
        return 1

    loaded_shards, skipped_shards = load_all_shards(out_dir)

    units = build_analysis_units(loaded_shards)
    write_csv(out_dir / "analysis_units.csv", UNIT_COLUMNS, units)

    stratum_groups = group_units_by_stratum(units)

    tau_rows = build_tau_table(stratum_groups, n_boot=n_boot, seed=seed)
    write_csv(out_dir / "analysis_tau.csv", TAU_COLUMNS, tau_rows)

    per_network_rows = build_per_network_table(stratum_groups)
    write_csv(out_dir / "analysis_per_network.csv", PER_NETWORK_COLUMNS, per_network_rows)

    effect_size_rows = build_effect_size_table(stratum_groups)
    write_csv(out_dir / "analysis_effect_sizes.csv", EFFECT_SIZE_COLUMNS, effect_size_rows)

    contradiction_rows = build_contradiction_table(loaded_shards)
    write_csv(out_dir / "analysis_contradiction.csv", CONTRADICTION_COLUMNS, contradiction_rows)

    separation_rows = build_separation_table(stratum_groups)
    write_csv(out_dir / "analysis_separation_strata.csv", SEPARATION_COLUMNS, separation_rows)

    xarm_bins_rows = build_xarm_bins(loaded_shards)
    write_csv(out_dir / "xarm_bins.csv", XARM_BINS_COLUMNS, xarm_bins_rows)

    paired_rows, xarm_pairing_summary = paired_cross_arm(loaded_shards, xarm_bins_rows, n_boot=n_boot, seed=seed)
    write_csv(out_dir / "xarm_paired_test.csv", XARM_PAIRED_COLUMNS, paired_rows)

    hist_rows = build_xarm_intensity_histogram(loaded_shards)
    write_csv(out_dir / "xarm_intensity_histogram.csv", XARM_HIST_COLUMNS, hist_rows)

    status_counts: dict[str, int] = defaultdict(int)
    for u in units:
        status_counts[f'{u["arm"]}:{u["status"]}'] += 1
    n_unreached_units = sum(1 for u in units if u["radius"] == UNREACHED)

    checks = {
        "decomposition_identity": check_decomposition_identity(loaded_shards),
        "n_draws_1000_non_censored": check_n_draws(loaded_shards),
        "n_unreached_units": n_unreached_units,
        "censored_cells": check_censored_cells(loaded_shards),
        "tiered_rate0_S_exactly_1": check_tiered_rate0(loaded_shards),
    }

    summary = {
        "out_dir": str(out_dir),
        "n_boot": n_boot,
        "seed": seed,
        "frame": {
            "n_rows": len(frame_rows),
            "frame_sha256": frame_info["frame_sha256"],
            "source_instances_sha256": frame_info.get("source_instances_sha256"),
        },
        "n_shards_read": len(loaded_shards),
        "n_shards_skipped": len(skipped_shards),
        "skipped_shards": skipped_shards,
        "shards_read": [
            {
                "shard_id": s["marker"]["shard_id"], "kind": s["kind"],
                "instances_sha256": s["marker"].get("instances_sha256"),
                "cells_sha256": s["marker"].get("cells_sha256"),
            }
            for s in loaded_shards
        ],
        "n_units": len(units),
        "unit_status_counts": dict(sorted(status_counts.items())),
        "n_strata_with_data": len(stratum_groups),
        "strata": sorted(stratum_groups),
        "n_strata_primary": sum(
            1 for us in stratum_groups.values() if us[0]["stratum_class"] == "primary"
        ),
        "n_strata_supplementary": sum(
            1 for us in stratum_groups.values() if us[0]["stratum_class"] == "supplementary"
        ),
        "xarm_cross_arm_pairing": xarm_pairing_summary,
        "checks": checks,
    }
    write_json(out_dir / "analysis_summary.json", summary)

    print(json.dumps({
        "n_shards_read": summary["n_shards_read"],
        "n_shards_skipped": summary["n_shards_skipped"],
        "n_units": summary["n_units"],
        "n_strata_with_data": summary["n_strata_with_data"],
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
