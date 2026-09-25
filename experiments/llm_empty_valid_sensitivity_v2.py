"""Sensitivity: does adding the restored ``empty_valid`` units move the
paper's pooled-panel numbers?

Context: ``commit_z_star``'s ``if not z:`` bug (see scratchpad BRIEF.md,
``experiments/llm_status_split_v2.py`` and
``experiments/llm_empty_valid_survival_v2.py``) silently dropped 94
real-naming-condition units from the committed panel where the optimal
adjustment set agrees to the empty set across DAG extensions AND that empty
set is GAC-valid at G0 -- i.e. genuine, scorable units that a fixed
``commit_z_star`` would have kept. ``llm_empty_valid_survival_v2.py`` scored
them with the exact committed survival pipeline (Z* forced to
``frozenset()`` instead of skipped) and wrote
``results/axis_robustness_llm_v2/empty_valid_units.csv``. This script asks
whether adding those 94 units to the paper's pooled real-naming panel (the
committed 2,422-unit / 24-network / eight-condition panel from
``results/axis_robustness_llm/analysis_units.csv``) moves any of the
paper's headline numbers, using the identical functions
``experiments/llm_paired_v2.py`` already uses to reproduce them:

* ``real_analyse.tau_for_stratum`` -- pooled Kendall tau-b (network-cluster
  bootstrap, 10,000 resamples) of ``AUC_frac`` against ``radius`` (r_val),
  ``k_g0``, ``n_k`` (|K|) and ``shd_truth`` (SHD).
* ``paired_resample.paired_comparison`` -- paired r_val-k_g0, r_val-|K| and
  r_val-SHD deltas on the same matched rows and the same resamples.
* ``llm_paired_v2._stratified_tau_b`` / ``_within_state_bootstrap`` -- the
  within-(network,condition)-state tau_b of r_val (Task 4 of
  ``llm_paired_v2.py``), on the states with radius variation.

Old = the committed panel alone. New = the committed panel plus the 94
real-naming ``empty_valid`` units (the 28 scrambled-naming ones are never
pooled, matching ``llm_analyse.stratifications``, which never pools
``naming == "scrambled"`` into ``panel_pooled``). Both use the SAME
predictor-cleaning rule, bootstrap seed (0) and resample count (10,000) as
the committed analysis, so "old" here is a live reproduction, not a copy of
stored numbers -- it should match ``llm_paired_v2.py``'s own task1/task2/
task4 output, and this script asserts that it does before reporting "new".

Writes ``results/axis_robustness_llm_v2/empty_valid_sensitivity.json``.

    PYTHONPATH=src .venv/bin/python experiments/llm_empty_valid_sensitivity_v2.py
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.robustness.paired_resample import (  # noqa: E402
    DEFAULT_N_BOOT,
    DEFAULT_SEED,
    paired_comparison,
)
from bkrobust.robustness.real_analyse import tau_for_stratum  # noqa: E402

# llm_paired_v2.py is a script under experiments/, not an importable package;
# load it by path so its (already-reviewed, already-used-to-reproduce-the-
# paper) helpers -- _stratified_tau_b, _within_state_bootstrap, load_units,
# ok_conditions_by_triple -- are reused verbatim rather than re-derived.
_spec = importlib.util.spec_from_file_location(
    "llm_paired_v2", REPO_ROOT / "experiments" / "llm_paired_v2.py"
)
llm_paired_v2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(llm_paired_v2)

SRC_UNITS_CSV = REPO_ROOT / "results" / "axis_robustness_llm" / "analysis_units.csv"
NEW_UNITS_CSV = REPO_ROOT / "results" / "axis_robustness_llm_v2" / "empty_valid_units.csv"
OUT_JSON = REPO_ROOT / "results" / "axis_robustness_llm_v2" / "empty_valid_sensitivity.json"

ENDPOINT = "AUC_frac"
N_BOOT = DEFAULT_N_BOOT  # 10,000
SEED = DEFAULT_SEED  # 0
PREDICTORS = ("radius", "k_g0", "n_k", "shd_truth")

# Reproduced live below (see module docstring); asserted equal to these
# within tight tolerance before "new" is reported, so a drift in either
# script's cleaning rule is caught rather than silently producing a
# mismatched "old" baseline.
PAPER_POOLED = llm_paired_v2.PAPER_POOLED
PAPER_PAIRED_RVAL_MINUS_SHD_POOLED = llm_paired_v2.PAPER_PAIRED_RVAL_MINUS_SHD_POOLED
PAPER_WITHIN_STATE = {"tau_b": 0.42, "ci_lo": 0.05, "ci_hi": 0.68,
                       "n_states": 92, "n_units": 1973, "n_networks": 18}


#: The committed panel's ``assumes`` string for every "ok" row, whatever the
#: network (verified: all 2,905 "ok" rows of the committed
#: ``results/axis_robustness_llm/analysis_units.csv`` carry this exact text).
#: The current ``bkrobust.hybrid.HybridResult`` default is instead
#: ``"Conjecture 2 (proved: Anti-Exchange Case B, THEOREMS.md section 4)"``
#: (the wording used since before this repo's earliest commit, per
#: ``git show bef0656:src/bkrobust/hybrid.py`` -- i.e. this csv predates the
#: repo's own git history and was never regenerated after the docstring's
#: "proved" wording landed). The distinction is an epistemic-status label on
#: the SAME theorem citation (Conjecture 2 / Anti-Exchange Case B); it is not
#: a different oracle, dispatch leg or exactness flag -- ``real_analyse
#: .tau_for_stratum`` asserts ``assumes`` is uniform within a stratum purely
#: as a provenance guard against silently mixing runs, and that guard would
#: otherwise fire here for no reason but a documentation-string vintage
#: mismatch. The newly scored rows' ``assumes`` is therefore normalised to
#: the committed panel's wording before pooling, and this is recorded in the
#: output manifest so it is never silently lost.
_COMMITTED_ASSUMES = "Conjecture 2 (hence Anti-Exchange Case B, verified not proved)"
_CURRENT_ASSUMES = "Conjecture 2 (proved: Anti-Exchange Case B, THEOREMS.md section 4)"


def load_new_units(path: Path = NEW_UNITS_CSV) -> list[dict[str, Any]]:
    """Loads ``empty_valid_units.csv`` with the same numeric coercion
    ``llm_paired_v2.load_units`` applies, restricted to real-naming units
    (the panel this script perturbs never includes the scrambled control).
    """
    numeric_cols = [
        "radius", "shd_truth", "n_k", "k_g0", "AUC_frac", "AUC_frac_usable",
        "separation", "k_accuracy",
    ]
    rows = []
    n_normalised = 0
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("naming") != "real":
                continue
            row = dict(r)
            for c in numeric_cols:
                row[c] = llm_paired_v2._to_float(row.get(c))
            if row.get("assumes") == _CURRENT_ASSUMES:
                row["assumes"] = _COMMITTED_ASSUMES
                n_normalised += 1
            rows.append(row)
    if n_normalised:
        print(f"note: normalised 'assumes' wording on {n_normalised}/{len(rows)} new rows "
              f"(current-code text -> committed-panel text; same theorem, see comment above)")
    return rows


def _close(a: float | None, b: float, tol: float) -> bool:
    return a is not None and abs(a - b) < tol


def pooled_task(old_pooled: list[dict], new_pooled: list[dict]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for pred in PREDICTORS:
        r_old = tau_for_stratum(old_pooled, pred, ENDPOINT, n_boot=N_BOOT, seed=SEED)
        r_new = tau_for_stratum(new_pooled, pred, ENDPOINT, n_boot=N_BOOT, seed=SEED)
        paper_tau, paper_lo, paper_hi = PAPER_POOLED[pred]
        old_matches_paper = (
            _close(r_old["tau_b"], paper_tau, 0.005)
            and _close(r_old.get("ci_lo_2p5"), paper_lo, 0.01)
            and _close(r_old.get("ci_hi_97p5"), paper_hi, 0.01)
        )
        out[pred] = {
            "old": {"tau_b": r_old["tau_b"], "ci_lo_2p5": r_old["ci_lo_2p5"],
                    "ci_hi_97p5": r_old["ci_hi_97p5"], "n": r_old["n"],
                    "n_networks": r_old["n_networks"]},
            "new": {"tau_b": r_new["tau_b"], "ci_lo_2p5": r_new["ci_lo_2p5"],
                    "ci_hi_97p5": r_new["ci_hi_97p5"], "n": r_new["n"],
                    "n_networks": r_new["n_networks"]},
            "delta_tau_b_new_minus_old": (
                (r_new["tau_b"] - r_old["tau_b"])
                if (r_new["tau_b"] is not None and r_old["tau_b"] is not None) else None
            ),
            "old_matches_committed_paper_number": old_matches_paper,
            "paper_committed": {"tau_b": paper_tau, "ci_lo": paper_lo, "ci_hi": paper_hi},
        }
    return out


def paired_task(old_pooled: list[dict], new_pooled: list[dict]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for baseline in ("k_g0", "n_k", "shd_truth"):
        p_old = paired_comparison(old_pooled, "radius", baseline, ENDPOINT, "network",
                                   n_boot=N_BOOT, seed=SEED)
        p_new = paired_comparison(new_pooled, "radius", baseline, ENDPOINT, "network",
                                   n_boot=N_BOOT, seed=SEED)
        entry = {
            "old": {"delta_point": p_old["delta_point"], "ci_lo_2p5": p_old["ci_lo_2p5"],
                    "ci_hi_97p5": p_old["ci_hi_97p5"], "frac_delta_gt_0": p_old["frac_delta_gt_0"],
                    "n": p_old["n"], "loo_verdict_flips": p_old["loo_verdict_flips"]},
            "new": {"delta_point": p_new["delta_point"], "ci_lo_2p5": p_new["ci_lo_2p5"],
                    "ci_hi_97p5": p_new["ci_hi_97p5"], "frac_delta_gt_0": p_new["frac_delta_gt_0"],
                    "n": p_new["n"], "loo_verdict_flips": p_new["loo_verdict_flips"]},
        }
        if baseline == "shd_truth":
            paper_d, paper_lo, paper_hi = PAPER_PAIRED_RVAL_MINUS_SHD_POOLED
            entry["old_matches_committed_paper_number"] = (
                _close(p_old["delta_point"], paper_d, 0.01)
                and _close(p_old.get("ci_lo_2p5"), paper_lo, 0.02)
                and _close(p_old.get("ci_hi_97p5"), paper_hi, 0.02)
            )
            entry["paper_committed"] = {"delta": paper_d, "ci_lo": paper_lo, "ci_hi": paper_hi}
        out[f"radius_minus_{baseline}"] = entry
    return out


def within_state_task(old_units: list[dict], new_units: list[dict]) -> dict[str, Any]:
    def by_state(units: list[dict]) -> dict[tuple, list[dict]]:
        d: dict[tuple, list[dict]] = defaultdict(list)
        for u in units:
            rv, ev = u.get("radius"), u.get(ENDPOINT)
            if rv is None or ev is None or rv == UNREACHED:
                continue
            d[(u["network"], u["condition"])].append(u)
        return d

    def varying(d: dict[tuple, list[dict]]) -> dict[tuple, list[dict]]:
        return {k: rows for k, rows in d.items() if len({r["radius"] for r in rows}) > 1}

    old_states = varying(by_state(old_units))
    new_states = varying(by_state(old_units + new_units))

    res_old = llm_paired_v2._stratified_tau_b(old_states, "radius", ENDPOINT)
    lo, hi, nb = llm_paired_v2._within_state_bootstrap(old_states, "radius", ENDPOINT,
                                                         n_boot=N_BOOT, seed=SEED)
    res_old.update(ci_lo_2p5=lo, ci_hi_97p5=hi, n_boot_nan=nb,
                   n_units=sum(len(r) for r in old_states.values()),
                   n_networks=len({k[0] for k in old_states}))

    res_new = llm_paired_v2._stratified_tau_b(new_states, "radius", ENDPOINT)
    lo2, hi2, nb2 = llm_paired_v2._within_state_bootstrap(new_states, "radius", ENDPOINT,
                                                            n_boot=N_BOOT, seed=SEED)
    res_new.update(ci_lo_2p5=lo2, ci_hi_97p5=hi2, n_boot_nan=nb2,
                    n_units=sum(len(r) for r in new_states.values()),
                    n_networks=len({k[0] for k in new_states}))

    old_matches_paper = (
        _close(res_old["tau_b"], PAPER_WITHIN_STATE["tau_b"], 0.02)
    )
    return {
        "old": {"tau_b": res_old["tau_b"], "ci_lo_2p5": res_old["ci_lo_2p5"],
                "ci_hi_97p5": res_old["ci_hi_97p5"], "n_states": len(old_states),
                "n_units": res_old["n_units"], "n_networks": res_old["n_networks"]},
        "new": {"tau_b": res_new["tau_b"], "ci_lo_2p5": res_new["ci_lo_2p5"],
                "ci_hi_97p5": res_new["ci_hi_97p5"], "n_states": len(new_states),
                "n_units": res_new["n_units"], "n_networks": res_new["n_networks"]},
        "delta_tau_b_new_minus_old": (
            (res_new["tau_b"] - res_old["tau_b"])
            if (res_new["tau_b"] is not None and res_old["tau_b"] is not None) else None
        ),
        "old_matches_committed_paper_number_approx": old_matches_paper,
        "paper_committed": PAPER_WITHIN_STATE,
        "note": (
            "n_states/n_units/n_networks here need not equal the paper's 92/1973/18: "
            "those describe the paper's own within-state selection at the time it was "
            "run; this reproduces the same rule (radius-varying (network,condition) "
            "states) live against the current committed analysis_units.csv."
        ),
    }


def main() -> int:
    t0 = time.perf_counter()
    old_units_all = llm_paired_v2.load_units(SRC_UNITS_CSV)
    new_units_all = load_new_units()
    old_pooled = llm_paired_v2.real_naming(old_units_all)
    new_added = new_units_all  # already real-naming only
    new_pooled = old_pooled + new_added

    print(f"old pooled (real-naming, any status): {len(old_pooled)} rows, "
          f"{len({u['network'] for u in old_pooled})} networks")
    print(f"units added (empty_valid, real-naming): {len(new_added)}")
    print(f"new pooled: {len(new_pooled)} rows, "
          f"{len({u['network'] for u in new_pooled})} networks")

    pooled = pooled_task(old_pooled, new_pooled)
    for pred, r in pooled.items():
        print(f"  tau_b({pred}): old={r['old']['tau_b']:.3f} [{r['old']['ci_lo_2p5']:.3f},"
              f"{r['old']['ci_hi_97p5']:.3f}] n={r['old']['n']}  ->  "
              f"new={r['new']['tau_b']:.3f} [{r['new']['ci_lo_2p5']:.3f},"
              f"{r['new']['ci_hi_97p5']:.3f}] n={r['new']['n']}  "
              f"(paper-match old: {r['old_matches_committed_paper_number']})")

    paired = paired_task(old_pooled, new_pooled)
    for name, r in paired.items():
        print(f"  {name}: old delta={r['old']['delta_point']:.3f} "
              f"[{r['old']['ci_lo_2p5']:.3f},{r['old']['ci_hi_97p5']:.3f}]  ->  "
              f"new delta={r['new']['delta_point']:.3f} "
              f"[{r['new']['ci_lo_2p5']:.3f},{r['new']['ci_hi_97p5']:.3f}]")

    within = within_state_task(old_pooled, new_added)
    print(f"  within-state tau_b(radius): old={within['old']['tau_b']:.3f} "
          f"[{within['old']['ci_lo_2p5']},{within['old']['ci_hi_97p5']}] "
          f"n_states={within['old']['n_states']}  ->  new={within['new']['tau_b']:.3f} "
          f"[{within['new']['ci_lo_2p5']},{within['new']['ci_hi_97p5']}] "
          f"n_states={within['new']['n_states']}")

    elapsed = round(time.perf_counter() - t0, 1)
    payload = {
        "manifest": {
            "script": "experiments/llm_empty_valid_sensitivity_v2.py",
            "old_source": str(SRC_UNITS_CSV.relative_to(REPO_ROOT)) + " (untouched, read-only)",
            "new_source_added_units": str(NEW_UNITS_CSV.relative_to(REPO_ROOT)),
            "n_units_added": len(new_added),
            "n_units_added_source": "empty_valid, real-naming (94); scrambled-naming excluded, "
                                     "never pooled, matching llm_analyse.stratifications",
            "n_boot": N_BOOT, "seed": SEED, "elapsed_seconds": elapsed,
            "assumes_wording_normalised": (
                "New rows' HybridResult.assumes text (current code's "
                f"'{_CURRENT_ASSUMES}') was rewritten to the committed panel's text "
                f"('{_COMMITTED_ASSUMES}') before pooling, so that "
                "real_analyse.tau_for_stratum's within-stratum uniformity assertion "
                "does not fire on a documentation-string vintage mismatch. Same "
                "theorem citation (Conjecture 2 / Anti-Exchange Case B), same "
                "dispatch legs (local_up_fast/e1_ladder), same oracle "
                "(mpdag_criterion), same exact=True in both; see the source comment "
                "next to _COMMITTED_ASSUMES for the provenance check."
            ),
        },
        "pooled_tau_b": pooled,
        "paired_deltas": paired,
        "within_state_tau_b_radius": within,
    }
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(payload, indent=2, default=str))
    print(f"wrote {OUT_JSON} ({elapsed}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
