"""Survival AUC for the "O* agrees to empty and is GAC-valid" units that
``commit_z_star``'s ``if not z:`` bug excluded from ``results/axis_robustness_llm``.

Context (see scratchpad BRIEF.md and ``experiments/llm_status_split_v2.py``,
whose output ``results/axis_robustness_llm_v2/status_split.json`` this script
reads): :func:`bkrobust.robustness.real_survival.commit_z_star` tests
``if not z:``, which is true both when the optimal adjustment set is
genuinely unidentified (``None``) and when every DAG extension of G0 agrees
it is the *empty set*. The second case never reached its own GAC validity
check and was silently folded into the ``optimal_set_undefined`` status.
``llm_status_split_v2.py`` already separated the two cases and, for the
units where the agreed-empty O* is GAC-valid at G0 ("empty_valid"), already
computed ``r_val`` with :func:`bkrobust.hybrid.breakdown_radius`.  That
script did **not** compute the survival endpoint (``AUC_frac``): the
flip-corruption sweep over the claim-depth grid that the committed
``analysis_units.csv`` reports for every other unit.  This script fills that
gap, for those units alone, forcing ``Z* = frozenset()`` (instead of the
buggy skip) and otherwise running the *exact* committed survival pipeline:

* ``K`` and ``G0``: :func:`bkrobust.robustness.llm_survival.build_state`,
  unchanged -- the same elicited-knowledge bundle, the same Meek closure.
* Corruption draws: :func:`bkrobust.robustness.real_survival.flip_draw` with
  the same seed family (:func:`bkrobust.robustness.llm_survival.panel_family`)
  and the same ``n_draws`` (1000, ``scripts/run_llm_survival_panel.sh``'s
  ``N_DRAWS`` default) as ``run_llm_survival.run_shard`` uses.
* Grid: :func:`bkrobust.robustness.real_survival.depth_grid`, unchanged.
* Survival check: :class:`bkrobust.robustness.real_survival.ClosureCache`,
  unchanged -- ``is_gac_valid_mpdag`` against the SAME corrupted graph every
  other unit in the group is scored against, only ``z_star`` is fixed to
  ``frozenset()`` instead of whatever ``commit_z_star`` would have committed.
* Endpoint: :func:`bkrobust.robustness.real_analyse._build_unit` with
  ``endpoint_kind="AUC_frac"`` and
  :func:`bkrobust.robustness.llm_analyse._contradiction_columns`, unchanged --
  the identical machinery ``llm_analyse.build_units`` calls for every
  committed unit.
* ``radius`` (r_val): :func:`bkrobust.robustness.real_survival.radius_columns`
  (which calls :func:`bkrobust.hybrid.breakdown_radius`), unchanged -- this
  recomputes r_val independently of ``status_split.json``'s own r_val, as a
  cross-check (see the ``r_val_status_split`` / ``r_val_matches_status_split``
  columns below).

No committed code path is modified. Units are grouped by ``(condition,
network)`` and, within a group, corruption draws are shared across units
exactly as ``run_llm_survival.run_shard`` shares them -- one draw per
``(depth, repetition)``, checked against every unit's fixed ``Z*`` at that
draw -- so this reproduces bit-for-bit what a shard would have written for
these rows had ``commit_z_star`` not short-circuited on ``z == set()``.

Scope: all 122 ``empty_valid`` units from ``status_split.json`` (94 on the
eight real-naming conditions of the paper's pooled panel, 28 on the
scrambled-naming control) are scored -- cheap here (31 groups, <=9 units
each, dense grids of at most 10 points, 1000 draws/point: same order of cost
as a single ordinary shard).

Output: ``results/axis_robustness_llm_v2/empty_valid_units.csv``, the same
column schema as ``results/axis_robustness_llm/analysis_units.csv``
(``bkrobust.robustness.llm_analyse.UNIT_COLUMNS``), plus a few traceability
columns appended at the end (never inserted into the middle of the schema):
``classification`` (always ``"empty_valid"``), ``r_val_status_split`` (the
r_val ``status_split.json`` already computed for this unit, for comparison),
``r_val_matches_status_split`` (whether this script's independently
recomputed ``radius`` agrees), ``source_script``.

Usage::

    PYTHONPATH=src .venv/bin/python experiments/llm_empty_valid_survival_v2.py
"""

from __future__ import annotations

import csv
import json
import statistics
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from bkrobust.robustness import llm_survival as ls  # noqa: E402
from bkrobust.robustness import real_survival as rs  # noqa: E402
from bkrobust.robustness.llm_analyse import (  # noqa: E402
    ENDPOINT,
    LLM_UNIT_COLUMNS,
    UNIT_COLUMNS,
    _contradiction_columns,
)
from bkrobust.robustness.real_analyse import _build_unit, write_csv  # noqa: E402
from bkrobust.robustness.run_llm_survival import DEFAULT_FRAME_DIR  # noqa: E402
from bkrobust.robustness.run_real_survival import load_frame  # noqa: E402

STATUS_SPLIT_JSON = REPO_ROOT / "results/axis_robustness_llm_v2/status_split.json"
SRC_ANALYSIS_UNITS_CSV = REPO_ROOT / "results/axis_robustness_llm/analysis_units.csv"
OUT_DIR = REPO_ROOT / "results/axis_robustness_llm_v2"
OUT_CSV = OUT_DIR / "empty_valid_units.csv"
OUT_JSON = OUT_DIR / "empty_valid_units.manifest.json"

#: Matches ``scripts/run_llm_survival_panel.sh``'s ``N_DRAWS`` default and
#: therefore ``run_llm_survival.run_shard``'s draws/grid-point, exactly.
N_DRAWS = 1000

EMPTY_Z = frozenset()

_START = time.time()


def _log(msg: str) -> None:
    print(f"[{time.time() - _START:7.1f}s] {msg}", flush=True)


def _balanced_triples_8of8(path: Path = SRC_ANALYSIS_UNITS_CSV) -> set[tuple[str, str, str]]:
    """``(network, x, y)`` triples every real-naming condition scored ``ok`` on
    the COMMITTED panel -- read-only, mirrors
    ``llm_analyse.balanced_triples``/``llm_paired_v2.ok_conditions_by_triple``
    at threshold 8, used only to fill this script's own ``in_balanced_panel``
    column consistently with the committed one.
    """
    ok_by_triple: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("naming") == "real" and r.get("status") == "ok":
                ok_by_triple[(r["network"], r["x"], r["y"])].add(r["condition"])
    return {t for t, cs in ok_by_triple.items() if len(cs) >= 8}


def score_group(
    condition: str, network: str, group_units: list[dict[str, Any]],
    bundle: dict[str, Any], frame: list[dict[str, Any]], balanced: set[tuple[str, str, str]],
) -> list[dict[str, Any]]:
    """Score every unit of one ``(condition, network)`` group, Z* forced empty.

    Mirrors ``run_llm_survival.run_shard`` exactly (same K/G0, same corruption
    draws shared across the group's units, same cache, same grid, same
    endpoint machinery), restricted to this group's units and with
    ``z_star = frozenset()`` fixed for every one of them instead of whatever
    ``commit_z_star`` would have committed.
    """
    dag, cpdag = rs.load_network(network)
    k, g0, g0_reason, shared = ls.build_state(bundle, condition, network, dag, cpdag)
    if g0 is None:
        raise RuntimeError(
            f"({condition}, {network}): g0 is None ({g0_reason}) but status_split.json "
            f"classified units here as empty_valid -- g0 must have existed when that "
            f"script ran. Bundle or frame may have drifted."
        )
    family = ls.panel_family(condition, network)
    pair_rows = {(r["x"], r["y"]): r for r in ls.distinct_pairs(frame, network)}

    scored: list[dict[str, Any]] = []
    for u in group_units:
        x, y = u["x"], u["y"]
        frame_row = pair_rows.get((x, y))
        if frame_row is None:
            raise RuntimeError(f"({condition}, {network}, {x}, {y}): no frame row found")
        rad_cols = rs.radius_columns(cpdag, g0, x, y, EMPTY_Z)
        irow: dict[str, Any] = {
            "shard_id": f"llm_empty_valid__{condition}__{network}",
            "arm": "llm",
            "gate": rs.GATE_NAME,
            "network": network, "x": x, "y": y,
            "coverage": None, "base_wrongness": None, "bw_abs_level": None,
            "analyst_replicate": None, "n_tiers": None,
            "frame_row_id": frame_row["row_id"],
            "largest_component_size": frame_row["largest_component_size"],
            "component_size": frame_row["component_size_recomputed"],
            "separation": frame_row["separation_recomputed"],
            "separation_status": frame_row["separation_status_recomputed"],
            "radius_committed_programmatic_k": frame_row["radius_recomputed"],
            **shared,
            "status": "ok", "z_star": [], "z_size": 0,
            **rad_cols,
        }
        scored.append(irow)

    cache = rs.ClosureCache(cpdag)
    cells_by_frid: dict[str, list[dict[str, Any]]] = defaultdict(list)
    grid = rs.depth_grid(len(k))
    for d in grid:
        n_contra = 0
        survived = [0] * len(scored)
        symdiffs: list[int] = []
        for rep in range(N_DRAWS):
            k_cor, _seed = rs.flip_draw(k, d, rep, family=family)
            key, g = cache.graph(k_cor)
            if g is None:
                n_contra += 1
                continue
            symdiffs.append(rs.directed_symdiff(g0, g))
            for i, irow in enumerate(scored):
                if cache.survived(key, g, irow["x"], irow["y"], EMPTY_Z):
                    survived[i] += 1
        for i, irow in enumerate(scored):
            crow = {
                "shard_id": irow["shard_id"], "arm": "llm", "gate": rs.GATE_NAME,
                "condition": condition, "model": shared["model"],
                "family": shared["family"], "naming": shared["naming"],
                "network": network, "x": irow["x"], "y": irow["y"],
                "coverage": None, "base_wrongness": None, "bw_abs": None,
                "analyst_replicate": None, "n_tiers": None,
                "frame_row_id": irow["frame_row_id"],
                "grid_kind": "d_claims", "grid_point": d,
                "frac_of_n_k": d / len(k),
                "n_k": len(k), "k_accuracy": shared["k_accuracy"],
                **rs.cell_statistics(N_DRAWS, n_contra, survived[i]),
                "symdiff_proxy_not_distance_median": (
                    statistics.median(symdiffs) if symdiffs else None
                ),
                "status": "ok",
            }
            cells_by_frid[irow["frame_row_id"]].append(crow)

    out_rows: list[dict[str, Any]] = []
    for u, irow in zip(group_units, scored):
        cells = cells_by_frid[irow["frame_row_id"]]
        unit = _build_unit(irow, cells, ENDPOINT)
        unit["r_status"] = irow.get("r_status")
        for col in LLM_UNIT_COLUMNS:
            unit[col] = irow.get(col)
        unit.update(_contradiction_columns(irow, cells))
        raw, cf = unit.get(ENDPOINT), unit.get("AUC_frac_contra_as_fail")
        unit["endpoint_contra_gap"] = abs(raw - cf) if (raw is not None and cf is not None) else None
        unit["in_balanced_panel"] = (network, irow["x"], irow["y"]) in balanced

        r_val_split = u.get("r_val")
        recomputed = unit.get("radius")
        matches = (
            r_val_split is not None and recomputed is not None
            and float(r_val_split) == float(recomputed)
        )
        unit["classification"] = "empty_valid"
        unit["r_val_status_split"] = r_val_split
        unit["r_val_matches_status_split"] = matches
        unit["source_script"] = "experiments/llm_empty_valid_survival_v2.py"
        out_rows.append(unit)
    return out_rows


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    split = json.loads(STATUS_SPLIT_JSON.read_text())
    real_conditions = set(split["meta"]["real_conditions"])
    empty_units = [u for u in split["units"] if u["classification"] == "empty_valid"]
    _log(f"{len(empty_units)} empty_valid units from status_split.json "
         f"({sum(1 for u in empty_units if u['condition'] in real_conditions)} real-condition)")

    groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for u in empty_units:
        groups[(u["condition"], u["network"])].append(u)
    _log(f"{len(groups)} (condition, network) groups")

    bundle, bundle_sha = ls.load_bundle()
    frame = load_frame(DEFAULT_FRAME_DIR)
    balanced = _balanced_triples_8of8()

    all_rows: list[dict[str, Any]] = []
    n_mismatch = 0
    for gi, ((condition, network), group_units) in enumerate(sorted(groups.items())):
        t0 = time.time()
        rows = score_group(condition, network, group_units, bundle, frame, balanced)
        all_rows.extend(rows)
        n_bad = sum(1 for r in rows if r["r_val_status_split"] is not None and not r["r_val_matches_status_split"])
        n_mismatch += n_bad
        _log(f"group {gi + 1}/{len(groups)} ({condition}, {network}): {len(rows)} units, "
             f"AUC_frac={[round(r['AUC_frac'], 3) if r['AUC_frac'] is not None else None for r in rows]}, "
             f"r_val mismatches={n_bad} ({time.time() - t0:.1f}s)")

    extra_cols = ["classification", "r_val_status_split", "r_val_matches_status_split", "source_script"]
    write_csv(OUT_CSV, [*UNIT_COLUMNS, *extra_cols], all_rows)
    _log(f"wrote {len(all_rows)} rows to {OUT_CSV}")
    if n_mismatch:
        _log(f"WARNING: {n_mismatch} units' recomputed radius disagrees with status_split.json's r_val")

    manifest = {
        "script": "experiments/llm_empty_valid_survival_v2.py",
        "status_split_json": str(STATUS_SPLIT_JSON.relative_to(REPO_ROOT)),
        "out_csv": str(OUT_CSV.relative_to(REPO_ROOT)),
        "n_draws": N_DRAWS,
        "n_units": len(all_rows),
        "n_units_real_naming": sum(1 for r in all_rows if r.get("naming") == "real"),
        "n_groups": len(groups),
        "n_r_val_mismatches": n_mismatch,
        "bundle_sha256": bundle_sha,
        "elapsed_s": round(time.time() - _START, 1),
    }
    OUT_JSON.write_text(json.dumps(manifest, indent=2, default=str))
    _log("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
