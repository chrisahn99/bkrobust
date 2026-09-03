"""Write scenario results to ``results/breakdown_radius_demo/`` as CSV and JSON.

One subdirectory per scenario. Everything the report and the figures read comes
from these files, so a number in ``report.md`` can always be traced to a row on
disk.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from bkrobust.demo.pipeline import robustness_frontier, run_scenario
from bkrobust.demo.space import atomic_moves

RESULTS_ROOT = Path("results/breakdown_radius_demo")


def elements_frame(result: dict[str, Any]) -> pd.DataFrame:
    """One row per element of the space."""
    g0, shells = result["g0"], result["shells"]
    rows = []
    for g, row in zip(result["space"], result["rows"]):  # noqa: B905 - equal by construction
        rec = dict(row)
        rec["shell"] = shells.get(g, -1)
        rec["moves_from_g0"] = " ; ".join(atomic_moves(g0, g))
        rec["n_moves"] = len(atomic_moves(g0, g))
        rows.append(rec)
    return pd.DataFrame(rows).sort_values(["shell", "edge_string"]).reset_index(drop=True)


def shells_frame(elements: pd.DataFrame) -> pd.DataFrame:
    """Shell index against counts and bias summaries."""
    grp = elements[elements["shell"] >= 0].groupby("shell")
    out = grp.agg(
        n_graphs=("edge_string", "count"),
        n_z_invalid=("z_valid", lambda s: int((~s).sum())),
        n_z_suboptimal=("z_is_optimal", lambda s: int((~s).sum())),
        mean_abs_bias=("mean_abs_bias", "mean"),
        max_abs_bias=("max_abs_bias", "max"),
    ).reset_index()
    return out


def write_scenario(label: str, n_draws: int = 200) -> dict[str, Any]:
    """Run one scenario and write every artefact for it.

    Args:
        label: ``"A"``, ``"B"`` or ``"C"``.
        n_draws: Coefficient draws per DAG.

    Returns:
        The full result dict, so callers can reuse it without recomputing.
    """
    result = run_scenario(label, n_draws=n_draws)
    out = RESULTS_ROOT / f"scenario_{label}"
    out.mkdir(parents=True, exist_ok=True)

    elements = elements_frame(result)
    elements.to_csv(out / "elements.csv", index=False)
    shells_frame(elements).to_csv(out / "shells.csv", index=False)

    frontier = pd.DataFrame(robustness_frontier(result))
    frontier.to_csv(out / "adjustment_sets.csv", index=False)

    summary = {
        "label": label,
        "relation_to_k_true": result["spec"]["relation"],
        "description": result["spec"]["description"],
        "g0": result["g0"].edge_string(),
        "cpdag": result["cpdag"].edge_string(),
        "truth": result["truth"].edge_string(),
        "adjustment_set_Z": sorted(result["z"]),
        "space_size": len(result["space"]),
        "n_covering_pairs": len(result["covers"]),
        "max_shell": max(result["shells"].values()),
        "radii": result["radii"],
        "witnesses": result["witnesses"],
        "eps_grid": result["eps_grid"],
        "standard_error_n1000": result["standard_error_n1000"],
        "mean_abs_effect": result["mean_abs_effect"],
        "metric_check": {k: v for k, v in result["metric"].items() if k != "violations"},
        "metric_violations": result["metric"].get("violations", []),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    return result


def write_all(n_draws: int = 200) -> dict[str, dict[str, Any]]:
    """Run and write all three scenarios."""
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    return {lab: write_scenario(lab, n_draws) for lab in ("A", "B", "C")}
