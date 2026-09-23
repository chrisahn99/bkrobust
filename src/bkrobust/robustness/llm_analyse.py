"""Analysis of the elicited-knowledge survival sweep.

Reads every completed shard under ``--out-dir`` (default
``results/axis_robustness_llm``) and writes the ranking tables. Endpoint and
tau machinery is **imported** from :mod:`bkrobust.robustness.real_analyse`
rather than reimplemented -- ``_build_unit``, ``tau_for_stratum``,
``load_all_shards``, ``stratum_weighted_medians`` -- so this arm and the
committed real arm cannot drift apart on what ``AUC_frac`` or a tau means.
What this module adds is the stratification: the stratum is the **model
condition**, not a coverage/base-wrongness cell.

Four stratifications, and they answer different questions
----------------------------------------------------------
``within_condition``
    One stratum per condition. Ranks the queries of a single elicited
    knowledge state against each other, across networks. This is the closest
    analogue of the committed analysis and inherits its weakness: inside one
    condition, ``|K|`` is still fixed per network.

``panel_pooled``
    One stratum over all real-naming conditions at once. This is the axis the
    sweep exists for. A unit is a ``(condition, network, X, Y)`` triple, so
    ``|K|``, ``k_g0`` and ``shd_truth`` all vary *within a network* -- one
    network contributes up to ten different knowledge states -- and the
    baselines finally have something to rank.

``panel_balanced``
    The pooled stratum restricted to the ``(network, X, Y)`` triples that
    every real-naming condition scored. See :func:`balanced_triples`: the
    pooled stratum is selected on ``|K|``, and this removes that channel.

The scrambled-naming conditions are a control arm. They get their own
``control_scrambled`` stratification and are **never** pooled into
``panel_pooled``: their claims are answers about relabelled variables.

Clustering
----------
Every bootstrap resamples **networks**, inherited unchanged. On the pooled
stratification a network's rows are dependent twice over -- one corrupted
state is shared by all of a network's pairs within a condition, and the same
network's structure recurs across conditions -- so the network cluster is the
correct and the conservative unit. It is not the instance.

``UNREACHED`` is excluded from every radius correlation by
``tau_for_stratum`` itself, and the count is reported per row.

Writes, all under ``--out-dir``
-------------------------------
``analysis_units.csv``, ``analysis_tau.csv``, ``analysis_panel.csv``,
``analysis_per_condition.csv``, ``analysis_summary.json``.

    PYTHONPATH=src .venv/bin/python -m bkrobust.robustness.llm_analyse
"""

from __future__ import annotations

import argparse
import json
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any

from bkrobust.core.conventions import UNREACHED
from bkrobust.robustness.real_analyse import (
    CENSORED_STATUS,
    DEFAULT_N_BOOT,
    DEFAULT_SEED,
    _build_unit,
    load_all_shards,
    stratum_weighted_medians,
    tau_for_stratum,
    write_csv,
    write_json,
)
from bkrobust.robustness.real_survival import FRAC_GRID, auc_over_grid

DEFAULT_OUT_DIR = "results/axis_robustness_llm"

#: The endpoint this arm reports. The flip grid is the committed one.
ENDPOINT = "AUC_frac"

#: The second endpoint, and on this arm it is not optional.
#:
#: Elicited ``K`` is dense on an almost-determined CPDAG and is Meek-consistent
#: by construction, so reversing a claim very often produces a claim set that
#: admits **no** MPDAG at all. The pilot saw ``child`` at depth 2 draw 923
#: contradictions in 1000, and every depth from 3 up draw 1000. ``S`` is
#: undefined at those grid points -- correctly, it is a conditional on
#: non-contradictory draws -- and :func:`auc_over_grid` then falls back to the
#: nearest grid point that *is* defined, so ``AUC_frac`` on this arm is carried
#: by the low-depth end of the curve.
#:
#: ``AUC_frac_contra_as_fail`` reads the same cells the other way: a corrupted
#: ``K`` that admits no MPDAG is scored as a failure rather than dropped, which
#: is the reading an analyst who would notice the contradiction and stop would
#: want. It is computed from the ``S_contra_as_fail`` column every cell already
#: carries, through the **same** :func:`auc_over_grid`, and is reported beside
#: the raw endpoint, never instead of it. The gap between them is the
#: contradiction-saturation diagnostic.
ENDPOINT_CONTRA_FAIL = "AUC_frac_contra_as_fail"

#: Predictors under test. The committed four, plus the two the panel makes
#: meaningful for the first time on real structure:
#:
#: ``n_k``
#:     on the committed corpus this is a deterministic function of
#:     ``(network, coverage)`` and cannot rank queries within a stratum. Here
#:     it is how many claims the model chose to assert.
#: ``shd_truth``
#:     identically zero at ``base_wrongness = 0`` on the committed corpus.
#:     Here it is the measured error of a real knowledge state.
#:
#: ``k_accuracy`` and ``assert_rate`` are elicitation-quality predictors that
#: exist only on this arm.
PREDICTORS: tuple[str, ...] = (
    "radius", "shd_truth", "n_k", "k_g0", "k_accuracy", "assert_rate",
)

#: The extra columns an elicited unit carries beyond the committed schema.
LLM_UNIT_COLUMNS: list[str] = [
    "condition", "model", "family", "naming", "elicit_arm",
    "n_asked", "n_asserted", "assert_rate", "k_accuracy",
    "n_k_correct", "n_k_wrong", "n_declined", "n_parse_fail", "n_not_reached",
    "n_off_skeleton", "k_source", "radius_committed_programmatic_k",
]

#: A unit belongs to more than one stratification at once -- every panel unit
#: is in both its own condition's stratum and the pooled one -- so the unit
#: table carries no single ``stratum`` column that would have to pick one.
#: ``condition`` and ``naming`` recover the grouping exactly.
UNIT_COLUMNS: list[str] = [
    "arm", "status", "in_balanced_panel",
    "network", "x", "y",
    "radius", "r_status", "shd_truth", "n_k", "k_g0",
    "AUC_frac", "AUC_frac_usable", "endpoint_gap",
    "AUC_frac_contra_as_fail", "endpoint_contra_gap",
    "n_grid_points_S_defined", "contradiction_rate_at_d1", "median_contradiction_rate",
    "n_grid_points", "n_grid_points_usable", "n_grid_points_censored", "endpoint_censored",
    "min_n_eval", "median_n_eval",
    "separation", "separation_status", "largest_component_size", "component_size",
    "dispatch_leg", "gate", "assumes", "shard_id", "frame_row_id",
    *LLM_UNIT_COLUMNS,
]


def build_units(loaded_shards: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Build one analysis unit per ``(condition, network, X, Y)``.

    The endpoint is computed by :func:`real_analyse._build_unit`, unchanged.
    Units whose instance status is not ``"ok"`` are kept with every endpoint
    ``None`` -- a condition that declined every pair is a measured outcome of
    the panel and must stay in the denominator.

    Args:
        loaded_shards: Output of :func:`real_analyse.load_all_shards`.

    Returns:
        Every unit row.
    """
    units: list[dict[str, Any]] = []
    for shard in loaded_shards:
        cells_by_frid: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for c in shard["cells"]:
            cells_by_frid[c["frame_row_id"]].append(c)
        for irow in shard["instances"]:
            cells = cells_by_frid.get(irow["frame_row_id"], [])
            unit = _build_unit(irow, cells, ENDPOINT)
            unit["r_status"] = irow.get("r_status")
            for col in LLM_UNIT_COLUMNS:
                unit[col] = irow.get(col)
            unit.update(_contradiction_columns(irow, cells))
            raw, cf = unit.get(ENDPOINT), unit.get(ENDPOINT_CONTRA_FAIL)
            unit["endpoint_contra_gap"] = (
                abs(raw - cf) if (raw is not None and cf is not None) else None
            )
            units.append(unit)
    return units


def _contradiction_columns(
    irow: dict[str, Any], cells: list[dict[str, Any]]
) -> dict[str, Any]:
    """The contra-as-fail endpoint and the saturation diagnostics for one unit.

    Args:
        irow: The instance row.
        cells: Its cell rows, censored ones included.

    Returns:
        ``AUC_frac_contra_as_fail``, ``n_grid_points_S_defined``,
        ``contradiction_rate_at_d1`` and ``median_contradiction_rate``. The
        gap between the two endpoints is filled in by the caller, which has
        both in hand.
    """
    live = [c for c in cells if c.get("status") != CENSORED_STATUS]
    rates = [
        c["contradiction_rate"] for c in live if c.get("contradiction_rate") is not None
    ]
    at_d1 = next(
        (c["contradiction_rate"] for c in sorted(live, key=lambda c: c["grid_point"])),
        None,
    )
    auc_cf = None
    if irow["status"] == "ok" and irow.get("n_k"):
        # Same helper, same nearest-grid-point rule; only the column read for
        # `S` differs, and `n_eval` becomes the full draw count because a
        # contradictory draw is now scored rather than excluded.
        remapped = {
            c["grid_point"]: {**c, "S": c.get("S_contra_as_fail"), "n_eval": c.get("n_draws")}
            for c in live
        }
        auc_cf = auc_over_grid(remapped, [f * irow["n_k"] for f in FRAC_GRID])
    return {
        ENDPOINT_CONTRA_FAIL: auc_cf,
        "n_grid_points_S_defined": sum(1 for c in live if c.get("S") is not None),
        "contradiction_rate_at_d1": at_d1,
        "median_contradiction_rate": statistics.median(rates) if rates else None,
    }


def stratifications(units: list[dict[str, Any]]) -> dict[str, dict[str, list[dict[str, Any]]]]:
    """Group the units the three ways the analysis reports.

    Args:
        units: Output of :func:`build_units`.

    Returns:
        ``{stratification: {stratum: rows}}`` for ``within_condition``,
        ``panel_pooled``, ``panel_balanced`` and ``control_scrambled``.
    """
    out: dict[str, dict[str, list[dict[str, Any]]]] = {
        "within_condition": defaultdict(list),
        "panel_pooled": defaultdict(list),
        "panel_balanced": defaultdict(list),
        "control_scrambled": defaultdict(list),
    }
    balanced = balanced_triples(units)
    for u in units:
        out["within_condition"][u["condition"]].append(u)
        if u.get("naming") == "real":
            out["panel_pooled"]["panel_all_real_naming"].append(u)
            if (u["network"], u["x"], u["y"]) in balanced:
                out["panel_balanced"]["panel_balanced_real_naming"].append(u)
        elif u.get("naming") == "scrambled":
            out["control_scrambled"]["control_all_scrambled"].append(u)
    return {k: dict(v) for k, v in out.items()}


def balanced_triples(units: list[dict[str, Any]]) -> set[tuple[str, str, str]]:
    """The ``(network, X, Y)`` triples every real-naming condition scored.

    Why this control exists, and it is not optional
    ------------------------------------------------
    A query enters the pooled stratum only when its condition's ``K`` produced
    a usable ``G0`` *and* determined an optimal adjustment set. Partial
    elicited ``K`` frequently does neither -- ``optimal_set_undefined`` and
    ``o_g0_extensions_intractable`` between them account for a large share of
    instance rows -- and a condition that asserts more claims clears both bars
    more often. So the pooled stratum is **selected on** the very predictor it
    is ranking: a positive ``n_k`` correlation there could be nothing but the
    selection.

    Restricting to the triples that are ``ok`` in *every* real-naming
    condition removes that channel entirely. Each condition then contributes
    the same queries, so ``|K|`` still varies within a network but the
    composition of the sample no longer does. The two taus are reported side
    by side; a ``n_k`` effect that survives only in the unbalanced stratum is
    a selection artefact and must be read as one.

    Args:
        units: Output of :func:`build_units`.

    Returns:
        The balanced triples. Empty if no triple is scorable everywhere, which
        is itself reportable rather than an error.
    """
    conds = {u["condition"] for u in units if u.get("naming") == "real"}
    if not conds:
        return set()
    ok_by_triple: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    for u in units:
        if u.get("naming") == "real" and u["status"] == "ok":
            ok_by_triple[(u["network"], u["x"], u["y"])].add(u["condition"])
    return {t for t, cs in ok_by_triple.items() if cs == conds}


def tau_rows(
    groups: dict[str, dict[str, list[dict[str, Any]]]], *, n_boot: int, seed: int
) -> list[dict[str, Any]]:
    """One tau row per ``(stratification, stratum, predictor)``.

    Args:
        groups: Output of :func:`stratifications`.
        n_boot: Bootstrap resamples.
        seed: Seed for every generator the call constructs.

    Returns:
        The tau rows, in a deterministic order. A stratum whose ``radius``
        rows carry more than one non-empty ``assumes`` string is recorded with
        status ``mixed_r_assumes`` and no tau, rather than aborting the run:
        the sweep is long, and a mixed-assumption stratum is a finding to
        report, not a crash.
    """
    rows: list[dict[str, Any]] = []
    for strat_name in (
        "within_condition", "panel_pooled", "panel_balanced", "control_scrambled",
    ):
        for stratum in sorted(groups.get(strat_name, {})):
            srows = groups[strat_name][stratum]
            for endpoint in (ENDPOINT, ENDPOINT_CONTRA_FAIL):
                for pred in PREDICTORS:
                    try:
                        r = tau_for_stratum(srows, pred, endpoint, n_boot=n_boot, seed=seed)
                    except AssertionError as exc:
                        print(f"  !! {strat_name}/{stratum}/{endpoint}/{pred}: {exc}",
                              flush=True)
                        r = {
                            "predictor": pred, "endpoint": endpoint,
                            "status": "mixed_r_assumes", "n": None, "n_networks": None,
                            "stratum_n_units": len(srows), "n_excluded_unreached": None,
                            "n_excluded_nan": None, "tau_b": None, "p_value": None,
                            "ci_lo_2p5": None, "ci_hi_97p5": None, "n_boot_nan": None,
                            "loo_tau_min": None, "loo_tau_max": None,
                            "loo_network_at_min": None, "loo_network_at_max": None,
                            "loo_verdict_flips": None, "r_assumes": str(exc),
                        }
                    rows.append({"stratification": strat_name, "stratum": stratum, **r})
    return rows


def per_condition_table(units: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """One descriptive row per condition: the two-dimensional panel variation.

    Args:
        units: Output of :func:`build_units`.

    Returns:
        A row per condition carrying ``|K|`` and accuracy summaries beside the
        survival endpoint and the radius distribution, so the ranking result
        can be read against the knowledge state that produced it.
    """
    by_cond: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for u in units:
        by_cond[u["condition"]].append(u)

    rows: list[dict[str, Any]] = []
    for cond in sorted(by_cond):
        us = by_cond[cond]
        ok = [u for u in us if u["status"] == "ok"]
        nets = sorted({u["network"] for u in us})
        # |K| and accuracy are per (condition, network), never per query: take
        # one value per network so a 95-pair network cannot outvote a 2-pair one.
        k_by_net = {u["network"]: u["n_k"] for u in us if u["n_k"] is not None}
        acc_by_net = {
            u["network"]: u["k_accuracy"] for u in us if u.get("k_accuracy") is not None
        }
        aucs = [u[ENDPOINT] for u in ok if u.get(ENDPOINT) is not None]
        aucs_cf = [
            u[ENDPOINT_CONTRA_FAIL] for u in ok if u.get(ENDPOINT_CONTRA_FAIL) is not None
        ]
        contra = [
            u["median_contradiction_rate"] for u in ok
            if u.get("median_contradiction_rate") is not None
        ]
        sdef = [
            u["n_grid_points_S_defined"] for u in ok
            if u.get("n_grid_points_S_defined") is not None
        ]
        radii = [u["radius"] for u in ok if u.get("r_status") == "ok"]
        pair_med, net_med = stratum_weighted_medians(us, ENDPOINT)
        statuses: dict[str, int] = defaultdict(int)
        for u in us:
            statuses[u["status"]] += 1
        rows.append({
            "condition": cond,
            "model": us[0].get("model"),
            "family": us[0].get("family"),
            "naming": us[0].get("naming"),
            "n_units": len(us),
            "n_units_ok": len(ok),
            "n_networks": len(nets),
            "sum_n_k": sum(k_by_net.values()),
            "median_n_k_per_network": statistics.median(k_by_net.values()) if k_by_net else None,
            "min_n_k": min(k_by_net.values()) if k_by_net else None,
            "max_n_k": max(k_by_net.values()) if k_by_net else None,
            "mean_accuracy_over_networks": (
                statistics.mean(acc_by_net.values()) if acc_by_net else None
            ),
            "median_accuracy_over_networks": (
                statistics.median(acc_by_net.values()) if acc_by_net else None
            ),
            "n_networks_g0_blocked": sum(
                1 for n in nets
                if all(u["status"] != "ok" for u in us if u["network"] == n)
            ),
            "median_AUC_frac_pair_weighted": pair_med,
            "median_AUC_frac_network_weighted": net_med,
            "n_AUC_defined": len(aucs),
            "median_AUC_frac_contra_as_fail": (
                statistics.median(aucs_cf) if aucs_cf else None
            ),
            "median_contradiction_rate": (
                statistics.median(contra) if contra else None
            ),
            "median_n_grid_points_S_defined": (
                statistics.median(sdef) if sdef else None
            ),
            "median_radius": statistics.median(radii) if radii else None,
            "share_radius_1": (
                sum(1 for r in radii if r == 1) / len(radii) if radii else None
            ),
            "n_radius_ok": len(radii),
            "n_radius_unreached": sum(
                1 for u in ok if u.get("radius") == UNREACHED
            ),
            "status_counts": json.dumps(dict(sorted(statuses.items()))),
        })
    return rows


def panel_variation_table(units: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """One row per ``(network, condition)``: the raw two-dimensional variation.

    This is the table that makes the design's premise checkable rather than
    asserted: within one network, does ``|K|`` actually move across conditions,
    and does the endpoint move with it?

    Args:
        units: Output of :func:`build_units`.

    Returns:
        A row per ``(network, condition)`` cell.
    """
    by_cell: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for u in units:
        by_cell[(u["network"], u["condition"])].append(u)

    rows: list[dict[str, Any]] = []
    for (net, cond) in sorted(by_cell):
        us = by_cell[(net, cond)]
        ok = [u for u in us if u["status"] == "ok"]
        aucs = [u[ENDPOINT] for u in ok if u.get(ENDPOINT) is not None]
        aucs_cf = [
            u[ENDPOINT_CONTRA_FAIL] for u in ok if u.get(ENDPOINT_CONTRA_FAIL) is not None
        ]
        radii = [u["radius"] for u in ok if u.get("r_status") == "ok"]
        rows.append({
            "network": net, "condition": cond,
            "model": us[0].get("model"), "naming": us[0].get("naming"),
            "n_pairs": len(us), "n_pairs_ok": len(ok),
            "n_k": us[0].get("n_k"), "n_asked": us[0].get("n_asked"),
            "assert_rate": us[0].get("assert_rate"),
            "k_accuracy": us[0].get("k_accuracy"),
            "n_k_wrong": us[0].get("n_k_wrong"),
            "k_g0": us[0].get("k_g0"), "shd_truth": us[0].get("shd_truth"),
            "g0_status": us[0].get("status") if not ok else "ok",
            "median_AUC_frac": statistics.median(aucs) if aucs else None,
            "median_AUC_frac_contra_as_fail": (
                statistics.median(aucs_cf) if aucs_cf else None
            ),
            "median_radius": statistics.median(radii) if radii else None,
            "n_radius_ok": len(radii),
        })
    return rows


def summarise(
    units: list[dict[str, Any]], taus: list[dict[str, Any]],
    loaded: list[dict[str, Any]], skipped: list[dict[str, Any]],
) -> dict[str, Any]:
    """The run-level summary written to ``analysis_summary.json``.

    Args:
        units: Output of :func:`build_units`.
        taus: Output of :func:`tau_rows`.
        loaded: Shards that loaded and verified.
        skipped: Shards whose digest did not match.

    Returns:
        The summary payload.
    """
    statuses: dict[str, int] = defaultdict(int)
    for u in units:
        statuses[u["status"]] += 1
    nets = sorted({u["network"] for u in units})
    conds = sorted({u["condition"] for u in units})
    k_spread = {}
    for n in nets:
        ks = sorted({u["n_k"] for u in units if u["network"] == n and u["n_k"] is not None})
        k_spread[n] = {"distinct_n_k": len(ks), "min": min(ks) if ks else None,
                       "max": max(ks) if ks else None}
    pooled = [
        t for t in taus
        if t["stratification"] == "panel_pooled" and t["endpoint"] == ENDPOINT
    ]
    balanced = [
        t for t in taus
        if t["stratification"] == "panel_balanced" and t["endpoint"] == ENDPOINT
    ]
    pooled_cf = [
        t for t in taus
        if t["stratification"] == "panel_pooled" and t["endpoint"] == ENDPOINT_CONTRA_FAIL
    ]
    return {
        "k_source": "results/elicit/knowledge.json",
        "select_knowledge_called": False,
        "endpoint": ENDPOINT,
        "n_shards_loaded": len(loaded),
        "n_shards_skipped_digest_mismatch": len(skipped),
        "skipped_shards": [s.get("shard_id") for s in skipped],
        "n_units": len(units),
        "n_units_ok": statuses.get("ok", 0),
        "unit_status_counts": dict(sorted(statuses.items())),
        "n_conditions": len(conds),
        "conditions": conds,
        "n_networks": len(nets),
        "within_network_n_k_spread": k_spread,
        "n_networks_with_n_k_spread": sum(
            1 for v in k_spread.values() if (v["distinct_n_k"] or 0) > 1
        ),
        "panel_pooled_tau": {
            t["predictor"]: {
                "tau_b": t["tau_b"], "ci_lo_2p5": t["ci_lo_2p5"],
                "ci_hi_97p5": t["ci_hi_97p5"], "n": t["n"],
                "n_networks": t["n_networks"], "status": t["status"],
                "loo_verdict_flips": t["loo_verdict_flips"],
            }
            for t in pooled
        },
        "panel_balanced_tau": {
            t["predictor"]: {
                "tau_b": t["tau_b"], "ci_lo_2p5": t["ci_lo_2p5"],
                "ci_hi_97p5": t["ci_hi_97p5"], "n": t["n"],
                "n_networks": t["n_networks"], "status": t["status"],
                "loo_verdict_flips": t["loo_verdict_flips"],
            }
            for t in balanced
        },
        "n_balanced_triples": len(balanced_triples(units)),
        "panel_pooled_tau_contra_as_fail": {
            t["predictor"]: {
                "tau_b": t["tau_b"], "ci_lo_2p5": t["ci_lo_2p5"],
                "ci_hi_97p5": t["ci_hi_97p5"], "n": t["n"],
                "n_networks": t["n_networks"], "status": t["status"],
                "loo_verdict_flips": t["loo_verdict_flips"],
            }
            for t in pooled_cf
        },
    }


def build_arg_parser() -> argparse.ArgumentParser:
    """The command-line interface.

    Returns:
        The parser.
    """
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    p.add_argument("--n-boot", type=int, default=DEFAULT_N_BOOT)
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    return p


def main(argv: list[str] | None = None) -> int:
    """Entry point.

    Args:
        argv: Argument vector, or ``None`` for ``sys.argv``.

    Returns:
        Process exit status.
    """
    args = build_arg_parser().parse_args(argv)
    out_dir = Path(args.out_dir)

    loaded, skipped = load_all_shards(out_dir)
    print(f"loaded {len(loaded)} shards, skipped {len(skipped)}", flush=True)
    units = build_units(loaded)
    print(f"built {len(units)} units", flush=True)

    groups = stratifications(units)
    for name, g in groups.items():
        print(f"  {name}: {len(g)} strata, "
              f"{sum(len(v) for v in g.values())} units", flush=True)

    taus = tau_rows(groups, n_boot=args.n_boot, seed=args.seed)
    print(f"computed {len(taus)} tau rows", flush=True)

    # Recoverable from the CSV alone: which units the selection control keeps.
    balanced = balanced_triples(units)
    for u in units:
        u["in_balanced_panel"] = (
            u.get("naming") == "real" and (u["network"], u["x"], u["y"]) in balanced
        )

    write_csv(out_dir / "analysis_units.csv", UNIT_COLUMNS, units)
    write_csv(out_dir / "analysis_tau.csv", list(taus[0]) if taus else [], taus)
    per_cond = per_condition_table(units)
    write_csv(out_dir / "analysis_per_condition.csv",
              list(per_cond[0]) if per_cond else [], per_cond)
    panel = panel_variation_table(units)
    write_csv(out_dir / "analysis_panel.csv", list(panel[0]) if panel else [], panel)
    write_json(out_dir / "analysis_summary.json",
               summarise(units, taus, loaded, skipped))
    print(f"wrote analysis tables under {out_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
