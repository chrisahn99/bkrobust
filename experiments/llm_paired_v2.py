"""Reproduction + extension of the real-network LLM-elicited-knowledge table
(Table tab:real-pooled in sections/06_experiments.tex, Table tab:balanced-full
and the paired r_val-SHD figure in sections/appendix/g_experiments.tex).

Reads the COMMITTED, read-only ``results/axis_robustness_llm/analysis_units.csv``
(never modified here) and writes everything new to
``results/axis_robustness_llm_v2/`` (never touches the committed directory).

What this does
---------------
1. Reproduces the pooled and balanced-8/8 tau_b of r_val (radius), k_g0,
   |K| (n_k) and SHD (shd_truth) against AUC_frac, using the same predictor/
   endpoint cleaning rule as ``real_analyse.tau_for_stratum`` (radius
   UNREACHED and missing predictor/endpoint dropped), and the same paired
   pooled r_val-SHD network-cluster bootstrap delta reported in the appendix
   (+0.265 [0.041, 0.456]), via ``bkrobust.robustness.paired_resample
   .paired_cluster_bootstrap`` (identical matched-cleaning + cluster-bootstrap
   machinery already used elsewhere in this repo for exactly this kind of
   paired claim).
2. Extends the paired comparison to r_val - k_g0 and r_val - |K|, and to two
   subset levels (>=7/8, >=6/8 real-naming conditions scorable) not present
   in the committed analysis_tau.csv, following the subset definitions of
   Table tab:balanced-full.
3. If a per-unit ``separation`` column already exists (it does --
   ``analysis_units.csv`` carries it, "measured" for 3060/5400 rows), computes
   tau_b and paired deltas for it on the pooled panel. ``phi_1`` does not
   exist as a column on this corpus (it is only computed for the synthetic
   flip-arm corpus by ``real_baselines_p6.build_extended_units``, which
   operates on a different unit schema) and is not cheaply derivable from the
   stored LLM-arm states, so it is reported as skipped, not recomputed.
4. A within-state stratified tau_b (pooling concordant/discordant pairs
   *within* each (network, condition) knowledge state, then bootstrapping by
   resampling networks) for r_val, to check against the paper's
   tau_b=0.42 [0.05, 0.68] (92 states, 1,973 units, 18 networks), and the same
   statistic for separation, which -- unlike |K| or k_g0 -- can vary from
   query to query inside one knowledge state and so is the only baseline that
   is a fair within-state comparator.

    PYTHONPATH=src .venv/bin/python experiments/llm_paired_v2.py
"""

from __future__ import annotations

import csv
import json
import math
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.robustness.paired_resample import (  # noqa: E402
    DEFAULT_N_BOOT,
    DEFAULT_SEED,
    paired_comparison,
)
from bkrobust.robustness.real_analyse import tau_for_stratum  # noqa: E402

SRC_UNITS_CSV = REPO_ROOT / "results" / "axis_robustness_llm" / "analysis_units.csv"
SRC_PANEL_CSV = REPO_ROOT / "results" / "axis_robustness_llm" / "analysis_panel.csv"
SRC_TAU_CSV = REPO_ROOT / "results" / "axis_robustness_llm" / "analysis_tau.csv"
OUT_DIR = REPO_ROOT / "results" / "axis_robustness_llm_v2"

ENDPOINT = "AUC_frac"
N_BOOT = DEFAULT_N_BOOT  # 10_000
SEED = DEFAULT_SEED  # 0

PAPER_POOLED = {
    "radius": (0.319, 0.122, 0.485),
    "k_g0": (0.214, 0.079, 0.339),
    "n_k": (0.185, 0.071, 0.327),
    "shd_truth": (0.054, -0.077, 0.200),
}
PAPER_BALANCED_8OF8 = {
    "radius": (0.499, 0.240, 0.624),
    "k_g0": (0.380, -0.074, 0.562),
    "n_k": (0.333, -0.125, 0.484),
    "shd_truth": (0.285, -0.167, 0.456),
}
PAPER_PAIRED_RVAL_MINUS_SHD_POOLED = (0.265, 0.041, 0.456)


def _to_float(v: Any) -> float | None:
    if v is None:
        return None
    s = str(v).strip()
    if s == "":
        return None
    try:
        f = float(s)
    except ValueError:
        return None
    if math.isnan(f):
        return None
    return f


def load_units(path: Path = SRC_UNITS_CSV) -> list[dict[str, Any]]:
    """Loads the committed per-unit table, coercing the numeric columns this
    script uses. Every other column is kept as the raw string from the CSV.
    Never writes to ``path``.
    """
    numeric_cols = [
        "radius", "shd_truth", "n_k", "k_g0", "AUC_frac", "AUC_frac_usable",
        "separation", "k_accuracy",
    ]
    rows = []
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            row = dict(r)
            for c in numeric_cols:
                row[c] = _to_float(row.get(c))
            rows.append(row)
    return rows


def real_naming(units: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [u for u in units if u.get("naming") == "real"]


def ok_conditions_by_triple(units: list[dict[str, Any]]) -> dict[tuple, set]:
    """``{(network, x, y): {conditions where status == "ok"}}`` over the
    real-naming units. Mirrors ``llm_analyse.balanced_triples`` but keeps the
    per-triple condition *set*, not just the fully-balanced (8/8) triples, so
    thresholds other than 8/8 can be derived from the same table.
    """
    d: dict[tuple, set] = defaultdict(set)
    for u in units:
        if u.get("naming") == "real" and u.get("status") == "ok":
            d[(u["network"], u["x"], u["y"])].add(u["condition"])
    return d


def subset_at_least(
    units: list[dict[str, Any]], ok_by_triple: dict[tuple, set], threshold: int
) -> list[dict[str, Any]]:
    """Real-naming units whose ``(network, x, y)`` triple is ``status=="ok"``
    in at least ``threshold`` of the real-naming conditions (Table
    tab:balanced-full's "8/8", ">=7/8", ">=6/8" levels). ``threshold=8`` is
    exactly ``llm_analyse.balanced_triples``.
    """
    qualifying = {t for t, cs in ok_by_triple.items() if len(cs) >= threshold}
    return [
        u for u in units
        if u.get("naming") == "real" and (u["network"], u["x"], u["y"]) in qualifying
    ]


# ---------------------------------------------------------------------------
# Task 1: reproduce pooled / balanced tau_b, and the paired pooled delta
# ---------------------------------------------------------------------------


def task1_reproduce(units: list[dict[str, Any]]) -> dict[str, Any]:
    pooled = real_naming(units)
    ok_by_triple = ok_conditions_by_triple(units)
    balanced_8of8 = subset_at_least(units, ok_by_triple, 8)

    out: dict[str, Any] = {"pooled": {}, "balanced_8of8": {}, "discrepancies": []}
    for label, rows, paper in (
        ("pooled", pooled, PAPER_POOLED),
        ("balanced_8of8", balanced_8of8, PAPER_BALANCED_8OF8),
    ):
        for pred in ("radius", "k_g0", "n_k", "shd_truth"):
            r = tau_for_stratum(rows, pred, ENDPOINT, n_boot=N_BOOT, seed=SEED)
            out[label][pred] = r
            paper_tau, paper_lo, paper_hi = paper[pred]
            match = (
                r["tau_b"] is not None
                and abs(r["tau_b"] - paper_tau) < 0.005
                and abs((r["ci_lo_2p5"] or float("nan")) - paper_lo) < 0.01
                and abs((r["ci_hi_97p5"] or float("nan")) - paper_hi) < 0.01
            )
            if not match:
                out["discrepancies"].append({
                    "stratum": label, "predictor": pred,
                    "paper": {"tau_b": paper_tau, "ci_lo": paper_lo, "ci_hi": paper_hi},
                    "reproduced": {
                        "tau_b": r["tau_b"], "ci_lo": r["ci_lo_2p5"], "ci_hi": r["ci_hi_97p5"],
                    },
                })

    paired = paired_comparison(pooled, "radius", "shd_truth", ENDPOINT, "network",
                                n_boot=N_BOOT, seed=SEED)
    out["paired_rval_minus_shd_pooled"] = paired
    paper_d, paper_lo, paper_hi = PAPER_PAIRED_RVAL_MINUS_SHD_POOLED
    match = (
        paired["delta_point"] is not None
        and abs(paired["delta_point"] - paper_d) < 0.01
        and abs((paired["ci_lo_2p5"] or float("nan")) - paper_lo) < 0.02
        and abs((paired["ci_hi_97p5"] or float("nan")) - paper_hi) < 0.02
    )
    if not match:
        out["discrepancies"].append({
            "stratum": "pooled", "predictor": "paired_rval_minus_shd",
            "paper": {"delta": paper_d, "ci_lo": paper_lo, "ci_hi": paper_hi},
            "reproduced": {
                "delta": paired["delta_point"], "ci_lo": paired["ci_lo_2p5"],
                "ci_hi": paired["ci_hi_97p5"],
            },
        })
    return out


# ---------------------------------------------------------------------------
# Task 2: paired deltas r_val - k_g0 and r_val - |K|, four subset levels
# ---------------------------------------------------------------------------


def task2_paired_subsets(units: list[dict[str, Any]]) -> list[dict[str, Any]]:
    ok_by_triple = ok_conditions_by_triple(units)
    subsets = {
        "pooled": real_naming(units),
        "balanced_8of8": subset_at_least(units, ok_by_triple, 8),
        "at_least_7of8": subset_at_least(units, ok_by_triple, 7),
        "at_least_6of8": subset_at_least(units, ok_by_triple, 6),
    }
    rows = []
    for level, rset in subsets.items():
        n_triples = len({(u["network"], u["x"], u["y"]) for u in rset})
        n_networks = len({u["network"] for u in rset})
        for baseline in ("k_g0", "n_k", "shd_truth"):
            r = paired_comparison(rset, "radius", baseline, ENDPOINT, "network",
                                   n_boot=N_BOOT, seed=SEED)
            r["level"] = level
            r["level_n_triples"] = n_triples
            r["level_n_networks_total"] = n_networks
            rows.append(r)
    return rows


# ---------------------------------------------------------------------------
# Task 3: separation (and phi_1, if present) on the pooled panel
# ---------------------------------------------------------------------------


def task3_separation(units: list[dict[str, Any]]) -> dict[str, Any]:
    pooled = real_naming(units)
    has_phi1 = any("phi_1" in u for u in units)
    out: dict[str, Any] = {
        "phi_1_available": has_phi1,
        "phi_1_skip_reason": (
            None if has_phi1 else
            "no phi_1 column on results/axis_robustness_llm/analysis_units.csv; "
            "phi_1 is only produced by real_baselines_p6.build_extended_units for the "
            "synthetic flip-arm corpus (flip arm, base_wrongness==0.0) via a many-minute "
            "combinatorial per-instance search over consistent extensions, and that "
            "function's unit schema does not apply to the LLM-elicited corpus's units "
            "(no base_wrongness / flip-arm structure here). Not cheaply computable from "
            "the stored states within the 30-minute budget; skipped, not recomputed."
        ),
    }
    n_sep_defined = sum(1 for u in pooled if u.get("separation") is not None)
    out["n_units_separation_defined"] = n_sep_defined
    out["n_units_pooled"] = len(pooled)
    out["tau_separation_pooled"] = tau_for_stratum(
        pooled, "separation", ENDPOINT, n_boot=N_BOOT, seed=SEED
    )
    out["paired_rval_minus_separation_pooled"] = paired_comparison(
        pooled, "radius", "separation", ENDPOINT, "network", n_boot=N_BOOT, seed=SEED
    )
    return out


# ---------------------------------------------------------------------------
# Task 4: within-state stratified tau_b (r_val, separation)
# ---------------------------------------------------------------------------


def _stratified_tau_b(
    rows_by_state: dict[Any, list[dict[str, Any]]], x_col: str, y_col: str
) -> dict[str, Any]:
    """Pooled concordant/discordant-pair tau_b across strata: only pairs of
    rows drawn from the *same* stratum (knowledge state) are counted, exactly
    mirroring the standard Kendall tau_b formula but restricted to
    within-stratum pairs. This is the natural pair-counting generalisation of
    "combine the states by stratification" the paper's within-state result
    (Sec 6 / App G) describes, and it reduces to ordinary tau_b when there is
    a single stratum.
    """
    C = D = Tx = Ty = 0
    n_states_used = 0
    for _state, rows in rows_by_state.items():
        n = len(rows)
        if n < 2:
            continue
        xs = [r[x_col] for r in rows]
        ys = [r[y_col] for r in rows]
        if len(set(xs)) < 2:
            # x constant in this state contributes no informative pairs, but
            # still counts x-ties (Ty pairs) below; harmless to include.
            pass
        used_this_state = False
        for i in range(n):
            for j in range(i + 1, n):
                dx = xs[i] - xs[j]
                dy = ys[i] - ys[j]
                if dx == 0 and dy == 0:
                    continue
                used_this_state = True
                if dx == 0:
                    Tx += 1
                elif dy == 0:
                    Ty += 1
                elif dx * dy > 0:
                    C += 1
                else:
                    D += 1
        if used_this_state:
            n_states_used += 1
    denom = math.sqrt((C + D + Tx) * (C + D + Ty))
    tau_b = (C - D) / denom if denom > 0 else None
    return {
        "tau_b": tau_b, "C": C, "D": D, "Tx": Tx, "Ty": Ty,
        "n_states_used": n_states_used, "n_states_total": len(rows_by_state),
    }


def _within_state_bootstrap(
    rows_by_state: dict[tuple, list[dict[str, Any]]], x_col: str, y_col: str,
    *, n_boot: int = N_BOOT, seed: int = SEED,
) -> tuple[float | None, float | None, int]:
    """Network-cluster bootstrap CI for :func:`_stratified_tau_b`: each
    resample draws networks with replacement and takes every state belonging
    to a drawn network (a network's states are its ``(network, condition)``
    keys; a network drawn twice contributes its states' rows twice, exactly
    as ``paired_resample.paired_cluster_bootstrap`` treats a drawn cluster).
    """
    state_keys = list(rows_by_state)
    network_of_state = {k: k[0] for k in state_keys}
    networks = sorted(set(network_of_state.values()))
    n_networks = len(networks)
    states_by_network: dict[str, list[tuple]] = defaultdict(list)
    for k, net in network_of_state.items():
        states_by_network[net].append(k)

    rng = np.random.default_rng(seed)
    net_arr = np.array(networks, dtype=object)
    vals = np.full(n_boot, np.nan, dtype=float)
    for b in range(n_boot):
        drawn = net_arr[rng.integers(0, n_networks, size=n_networks)]
        resampled: dict[Any, list[dict[str, Any]]] = {}
        for rep, net in enumerate(drawn):
            for state_key in states_by_network[net]:
                resampled[(rep, state_key)] = rows_by_state[state_key]
        res = _stratified_tau_b(resampled, x_col, y_col)
        if res["tau_b"] is not None:
            vals[b] = res["tau_b"]
    valid = vals[~np.isnan(vals)]
    if len(valid) < 100:
        return None, None, int(n_boot - len(valid))
    lo, hi = np.percentile(valid, [2.5, 97.5])
    return float(lo), float(hi), int(n_boot - len(valid))


def task4_within_state(units: list[dict[str, Any]]) -> dict[str, Any]:
    pooled = real_naming(units)
    # A knowledge state is one (network, condition) pair. Only rows with a
    # defined, non-UNREACHED radius and a defined endpoint enter the state
    # (matching tau_for_stratum's own cleaning rule).
    by_state: dict[tuple, list[dict[str, Any]]] = defaultdict(list)
    for u in pooled:
        rv, ev = u.get("radius"), u.get(ENDPOINT)
        if rv is None or ev is None or rv == UNREACHED:
            continue
        by_state[(u["network"], u["condition"])].append(u)

    # States with variation in radius (the paper's selection rule).
    varying_states = {
        k: rows for k, rows in by_state.items()
        if len({r["radius"] for r in rows}) > 1
    }
    n_k_const = all(
        len({r["n_k"] for r in rows if r.get("n_k") is not None}) <= 1
        for rows in varying_states.values()
    )
    k_g0_const = all(
        len({r["k_g0"] for r in rows if r.get("k_g0") is not None}) <= 1
        for rows in varying_states.values()
    )

    res_radius = _stratified_tau_b(varying_states, "radius", ENDPOINT)
    lo, hi, n_boot_nan = _within_state_bootstrap(varying_states, "radius", ENDPOINT)
    res_radius["ci_lo_2p5"], res_radius["ci_hi_97p5"], res_radius["n_boot_nan"] = lo, hi, n_boot_nan
    res_radius["n_units"] = sum(len(r) for r in varying_states.values())
    res_radius["n_networks"] = len({k[0] for k in varying_states})

    # Separation: cheap (already a stored column). Use the SAME states
    # (network,condition) as radius's within-state set -- these are the
    # states the paper's headline number is computed over -- so this is a
    # like-for-like comparison. separation is defined ("measured") only on a
    # subset of rows; undefined rows are dropped from this predictor only.
    sep_states: dict[tuple, list[dict[str, Any]]] = {}
    for k, rows in varying_states.items():
        sub = [r for r in rows if r.get("separation") is not None]
        if len(sub) >= 2:
            sep_states[k] = sub
    res_sep = _stratified_tau_b(sep_states, "separation", ENDPOINT)
    lo_s, hi_s, n_boot_nan_s = _within_state_bootstrap(sep_states, "separation", ENDPOINT)
    res_sep["ci_lo_2p5"], res_sep["ci_hi_97p5"], res_sep["n_boot_nan"] = lo_s, hi_s, n_boot_nan_s
    res_sep["n_units"] = sum(len(r) for r in sep_states.values())
    res_sep["n_networks"] = len({k[0] for k in sep_states})
    res_sep["n_states_with_separation_variation"] = sum(
        1 for rows in sep_states.values() if len({r["separation"] for r in rows}) > 1
    )

    return {
        "n_states_total": len(by_state),
        "n_states_with_radius_variation": len(varying_states),
        "n_k_constant_within_every_varying_state": n_k_const,
        "k_g0_constant_within_every_varying_state": k_g0_const,
        "radius": res_radius,
        "separation_same_states": res_sep,
        "paper_reference": {"tau_b": 0.42, "ci_lo": 0.05, "ci_hi": 0.68,
                             "n_states": 92, "n_units": 1973, "n_networks": 18},
    }


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def _json_default(o: Any) -> Any:
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, set):
        return sorted(o)
    return str(o)


def main() -> int:
    t0 = time.perf_counter()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    units = load_units()
    print(f"loaded {len(units)} units from {SRC_UNITS_CSV}")

    print("== task 1: reproduce pooled / balanced-8of8 tau_b + paired r_val-SHD ==")
    t1 = task1_reproduce(units)
    print(f"  discrepancies: {len(t1['discrepancies'])}")
    for d in t1["discrepancies"]:
        print(f"    MISMATCH {d['stratum']}/{d['predictor']}: paper={d['paper']} "
              f"reproduced={d['reproduced']}")

    print("== task 2: paired deltas r_val - {k_g0, n_k, shd_truth}, 4 subset levels ==")
    t2 = task2_paired_subsets(units)
    for r in t2:
        print(f"  {r['level']:16s} {r['predictor_a']}-{r['predictor_b']:10s} "
              f"n={r['n']:5d} delta={r['delta_point']}"
              f" CI=[{r['ci_lo_2p5']},{r['ci_hi_97p5']}] P(>0)={r['frac_delta_gt_0']}")

    print("== task 3: separation (+ phi_1 availability check) on pooled panel ==")
    t3 = task3_separation(units)
    print(f"  phi_1_available={t3['phi_1_available']}")
    print(f"  tau_b(separation)={t3['tau_separation_pooled']['tau_b']} "
          f"CI=[{t3['tau_separation_pooled']['ci_lo_2p5']},{t3['tau_separation_pooled']['ci_hi_97p5']}]")

    print("== task 4: within-state stratified tau_b (radius, separation) ==")
    t4 = task4_within_state(units)
    print(f"  radius: tau_b={t4['radius']['tau_b']} "
          f"CI=[{t4['radius']['ci_lo_2p5']},{t4['radius']['ci_hi_97p5']}] "
          f"n_states={t4['n_states_with_radius_variation']} n_units={t4['radius']['n_units']} "
          f"n_networks={t4['radius']['n_networks']}")
    print(f"  separation (same states): tau_b={t4['separation_same_states']['tau_b']} "
          f"CI=[{t4['separation_same_states']['ci_lo_2p5']},{t4['separation_same_states']['ci_hi_97p5']}] "
          f"n_units={t4['separation_same_states']['n_units']} "
          f"n_networks={t4['separation_same_states']['n_networks']}")

    elapsed = round(time.perf_counter() - t0, 1)
    print(f"total elapsed: {elapsed}s")

    with (OUT_DIR / "task1_reproduce.json").open("w", encoding="utf-8") as fh:
        json.dump(t1, fh, indent=2, default=_json_default)
    with (OUT_DIR / "task2_paired_subsets.json").open("w", encoding="utf-8") as fh:
        json.dump(t2, fh, indent=2, default=_json_default)
    with (OUT_DIR / "task3_separation.json").open("w", encoding="utf-8") as fh:
        json.dump(t3, fh, indent=2, default=_json_default)
    with (OUT_DIR / "task4_within_state.json").open("w", encoding="utf-8") as fh:
        json.dump(t4, fh, indent=2, default=_json_default)

    manifest = {
        "script": "experiments/llm_paired_v2.py",
        "source_units_csv": str(SRC_UNITS_CSV.relative_to(REPO_ROOT)),
        "source_units_csv_untouched": True,
        "out_dir": str(OUT_DIR.relative_to(REPO_ROOT)),
        "n_boot": N_BOOT, "seed": SEED,
        "n_units_loaded": len(units),
        "elapsed_seconds": elapsed,
        "outputs": [
            "task1_reproduce.json", "task2_paired_subsets.json",
            "task3_separation.json", "task4_within_state.json",
            "FINDINGS.md", "manifest.json",
        ],
    }
    with (OUT_DIR / "manifest.json").open("w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
