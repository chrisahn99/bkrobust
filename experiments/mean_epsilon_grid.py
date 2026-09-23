r"""Choose a defensible epsilon reporting grid for the mean bias profile ``mu(d)``.

``r_eps`` thresholds ``beta_up(d)``, the shell **maximum**. This script is
about a companion quantity, ``mu(d)`` (``src/bkrobust/epsilon/meanprofile.py``):
the subset-weighted **mean** of ``B`` over the retraction shell at depth ``d``.
Because ``mu`` averages instead of maximising, it is systematically smaller
than ``beta_up`` shell by shell, so the epsilon grid tuned for the max
(``{0.01, 0.05, 0.10, 0.25, 0.50, 1.00}``) is not automatically right for the
mean, and this script asks what grid is.

It does three things, none of which recomputes a single Meek closure --
everything here reads exhaustive per-shell data already committed under
``results/``:

1. **Characterises the bias distribution** -- individual per-state ``B``
   (via the shell-level min/median/max/frac-nonzero fields that already
   summarise it), the shell means ``mu(d)`` pooled over every instance and
   depth, and the per-instance ceiling ``B(Chat)`` (the bias at the top of
   the retraction poset, i.e. with every asserted orientation retracted).
   It looks explicitly for **mass points** -- values that recur exactly
   across independent instances, which a threshold must straddle rather
   than sit on, since a strict ``>`` test splits a tied mass point by
   floating-point dust alone.

2. **Scores candidate grids** by how well they discriminate between
   instances: how many distinct finite radii a grid produces per instance
   (more is better), how many grid cells land ``UNREACHED`` (threshold too
   high, wasted column) or pinned at ``r_val`` (threshold too low, adds no
   information beyond ``r_val`` itself), and whether any grid point sits
   within ``1e-6`` of a mass point (disqualifying).

3. **Guards against overfitting** by scoring every grid on several
   independent splits of the corpus -- real vs. synthetic networks, and
   small vs. large ``|K_G0|`` -- and flagging any grid (in particular a
   quantile-derived one) that scores well on one half and poorly on the
   other.

Every bias number here is conditional on the seeded linear-Gaussian SEM
attached to each instance's ground-truth DAG: the corpus ships no
observational data, so ``theta_z`` and every ``B`` value are SEM-conditional,
not estimated from data. Condition ``D_LLM`` only; every claim below is
scoped to the instances actually loaded.

Data reused (not recomputed):
    ``results/mean_vs_max/shells.jsonl``     -- 15 real instances, exhaustive
        per-shell enumeration (``mean_subsets`` is ``mu(d)``; also carries
        ``min``/``median_states``/``max``/``std_states``/``frac_nonzero``,
        the richest per-state summary available).
    ``results/mean_fullspace/shells.jsonl``  -- 79 instances (10 real + 69
        synthetic), per-shell retraction statistics (``retr_mean_subsets``
        is ``mu(d)``, ``retr_max`` is ``beta_up(d)``). This is the primary
        corpus for grid scoring and the real/synthetic overfitting check.
    ``results/final_table_eps_gt1/instances.jsonl`` -- 160 real-network
        query records with ``r_val`` and ``beta_top`` (the ceiling
        ``B(Chat)``), the largest available sample of the ceiling
        distribution.

Usage::

    PYTHONPATH=src .venv/bin/python experiments/mean_epsilon_grid.py [options]

Writes ``<out>/mass_points.json``, ``<out>/grid_scores.csv`` and
``<out>/summary.json``, plus ``figures/mean_eps_distribution.png`` and
``figures/mean_eps_grid_scores.png``.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

UNREACHED = -1

#: The legacy grid, tuned for beta_up (the shell max), not mu.
LEGACY_GRID = [0.01, 0.05, 0.10, 0.25, 0.50, 1.00]

#: Round, practitioner-recognisable fractions of the reported effect,
#: deliberately straddling the mass point found at exactly 1.0 (0.75 below,
#: 1.25 above) rather than landing on it.
ROUND_NUMBER_GRID = [0.02, 0.10, 0.25, 0.50, 0.75, 1.25, 2.00]

#: Candidate values to test for being a mass point in the pooled data.
MASS_POINT_CANDIDATES = [0.0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 2.00]


def log_spaced_grid(lo: float = 0.02, hi: float = 2.00, n: int = 6) -> list[float]:
    """A log-spaced candidate grid between ``lo`` and ``hi``.

    Args:
        lo: Smallest value.
        hi: Largest value.
        n: Number of points.

    Returns:
        Sorted list of ``n`` values, rounded to 4 significant places.
    """
    return sorted(float(f"{v:.4g}") for v in np.geomspace(lo, hi, n))


def quantile_derived_grid(pooled: np.ndarray, qs: tuple[float, ...] = (10, 25, 50, 75, 90, 97)) -> list[float]:
    """A grid read directly off quantiles of the pooled ``mu(d)`` sample.

    This is the grid the overfitting check (step 3) is designed to catch:
    it is tuned to *this* corpus's empirical distribution by construction.

    Args:
        pooled: Pooled non-degenerate ``mu(d)`` values.
        qs: Percentiles to read off.

    Returns:
        Sorted, deduplicated list of quantile values.
    """
    vals = sorted({round(float(v), 4) for v in np.percentile(pooled, qs)})
    return vals


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------


def _key(rec: dict[str, Any]) -> tuple[str, str, str]:
    return (rec["network"], rec["x"], rec["y"])


def load_mean_vs_max(root: Path) -> dict[str, Any]:
    """Load the 15-instance exhaustive real-network shell study.

    Args:
        root: Repository root.

    Returns:
        ``{"profiles": {key: {d: mean_subsets}}, "shells": {key: [rows]},
        "meta": {key: instance_meta}}``.
    """
    path = root / "results" / "mean_vs_max" / "shells.jsonl"
    rows = [json.loads(line) for line in path.open()]
    shells: dict[tuple, list[dict]] = {}
    for r in rows:
        shells.setdefault(_key(r), []).append(r)
    for k in shells:
        shells[k].sort(key=lambda r: r["d"])
    profiles = {k: {r["d"]: r["mean_subsets"] for r in v} for k, v in shells.items()}
    meta_all = json.loads((root / "results" / "mean_vs_max" / "instances.json").read_text())
    meta = {_key(m): m for m in meta_all["instances"] if m.get("status") == "ok"}
    return {"profiles": profiles, "shells": shells, "meta": meta}


def load_mean_fullspace(root: Path) -> dict[str, Any]:
    """Load the 79-instance (10 real + 69 synthetic) retraction-shell study.

    Args:
        root: Repository root.

    Returns:
        ``{"profiles": {key: {d: retr_mean_subsets}}, "shells": {key: [rows]},
        "meta": {key: instance_meta}}``. ``meta[key]["source"]`` is
        ``"real"`` or ``"synthetic"``.
    """
    path = root / "results" / "mean_fullspace" / "shells.jsonl"
    rows = [json.loads(line) for line in path.open()]
    shells: dict[tuple, list[dict]] = {}
    for r in rows:
        shells.setdefault(_key(r), []).append(r)
    for k in shells:
        shells[k].sort(key=lambda r: r["d"])
    profiles = {
        k: {r["d"]: r["retr_mean_subsets"] for r in v if r["retr_mean_subsets"] is not None}
        for k, v in shells.items()
    }
    meta_all = json.loads((root / "results" / "mean_fullspace" / "instances.json").read_text())
    meta = {_key(m): m for m in meta_all["instances"] if m.get("status") == "ok"}
    return {"profiles": profiles, "shells": shells, "meta": meta}


def load_final_table(root: Path) -> list[dict[str, Any]]:
    """Load the 160-record real-network ``r_val``/``beta_top`` table.

    Args:
        root: Repository root.

    Returns:
        Records with ``status == "ok"``.
    """
    path = root / "results" / "final_table_eps_gt1" / "instances.jsonl"
    recs = [json.loads(line) for line in path.open()]
    return [r for r in recs if r["status"] == "ok"]


# --------------------------------------------------------------------------
# r_val, r_mean
# --------------------------------------------------------------------------


def derive_r_val(rows: list[dict[str, Any]], max_field: str) -> int:
    """First depth whose shell maximum is non-zero, else UNREACHED.

    Mirrors ``r_val``: Theorem A gives ``beta(d) = 0`` for ``d < r_val``, so
    the first depth carrying a non-zero shell maximum *is* ``r_val``.

    Args:
        rows: Shell rows for one instance, sorted by ``d``.
        max_field: Field name holding the shell maximum (``"max"`` or
            ``"retr_max"``).

    Returns:
        The depth, or :data:`UNREACHED` if no shell ever carries non-zero
        bias (a fully robust instance).
    """
    for r in rows:
        v = r.get(max_field)
        if v is not None and v > 1e-9:
            return r["d"]
    return UNREACHED


def r_mean_from_profile(profile: dict[int, float], r_val: int, eps: float) -> int:
    """First depth ``d >= r_val`` with ``mu(d) > eps``, else UNREACHED.

    Applies the same first-crossing rule as
    ``bkrobust.epsilon.meanprofile.r_mean`` to an already-computed ``mu(d)``
    table, rather than recomputing it.

    Args:
        profile: ``{d: mu(d)}``, exhaustive per Theorem A below ``r_val``.
        r_val: The instance's validity radius. If :data:`UNREACHED`, the
            profile is identically zero and this returns UNREACHED
            immediately (never averaged, per house rule).
        eps: Threshold.

    Returns:
        The crossing depth, or :data:`UNREACHED`.
    """
    if r_val == UNREACHED:
        return UNREACHED
    for d in sorted(profile):
        if d < r_val:
            continue
        v = profile[d]
        if v is not None and v > eps:
            return d
    return UNREACHED


# --------------------------------------------------------------------------
# Mass points
# --------------------------------------------------------------------------


def mass_point_report(values: np.ndarray, tol: float = 1e-6) -> dict[str, dict[str, float]]:
    """Count exact hits near each candidate mass point.

    Args:
        values: Pooled sample.
        tol: Matching tolerance.

    Returns:
        ``{str(candidate): {"count": int, "frac": float}}``.
    """
    out = {}
    n = len(values)
    for c in MASS_POINT_CANDIDATES:
        hit = np.sum(np.abs(values - c) < tol)
        out[str(c)] = {"count": int(hit), "frac": float(hit) / n if n else 0.0}
    return out


def disqualified_grid_points(grid: list[float], mass_points: list[float], tol: float = 1e-6) -> list[float]:
    """Grid values sitting within ``tol`` of a confirmed mass point.

    Args:
        grid: Candidate grid.
        mass_points: Confirmed mass-point values (see :func:`mass_point_report`).
        tol: Matching tolerance.

    Returns:
        The offending subset of ``grid`` (empty if none).
    """
    return [g for g in grid for m in mass_points if abs(g - m) < tol]


# --------------------------------------------------------------------------
# Grid scoring
# --------------------------------------------------------------------------


def score_grid(
    instances: list[dict[str, Any]], grid: list[float]
) -> dict[str, float]:
    """Score one grid against one set of instances.

    An instance contributes ``{"profile": {d: mu}, "r_val": int}``.
    Degenerate instances (``r_val == 0``) must already be excluded by the
    caller, per house rule.

    Args:
        instances: List of ``{"profile": ..., "r_val": ...}``.
        grid: Candidate epsilon values.

    Returns:
        Aggregate scores: ``n_instances``, ``avg_distinct_finite``,
        ``avg_n_unreached``, ``avg_n_pinned``, ``frac_all_unreached``,
        ``frac_all_pinned``, ``grid_size``.
    """
    n = len(instances)
    distinct_counts, unreached_counts, pinned_counts = [], [], []
    all_unreached, all_pinned = 0, 0
    for inst in instances:
        radii = [r_mean_from_profile(inst["profile"], inst["r_val"], eps) for eps in grid]
        finite = [r for r in radii if r != UNREACHED]
        n_unreached = sum(1 for r in radii if r == UNREACHED)
        n_pinned = sum(1 for r in radii if r == inst["r_val"])
        distinct_counts.append(len(set(finite)))
        unreached_counts.append(n_unreached)
        pinned_counts.append(n_pinned)
        if n_unreached == len(grid):
            all_unreached += 1
        if n_pinned == len(grid):
            all_pinned += 1
    return {
        "n_instances": n,
        "grid_size": len(grid),
        "avg_distinct_finite": float(np.mean(distinct_counts)) if n else 0.0,
        "avg_n_unreached": float(np.mean(unreached_counts)) if n else 0.0,
        "avg_n_pinned": float(np.mean(pinned_counts)) if n else 0.0,
        "frac_all_unreached": all_unreached / n if n else 0.0,
        "frac_all_pinned": all_pinned / n if n else 0.0,
    }


def build_instance_list(
    profiles: dict[tuple, dict[int, float]],
    shells: dict[tuple, list[dict]],
    max_field: str,
    keys: list[tuple],
) -> list[dict[str, Any]]:
    """Assemble ``{"profile", "r_val", "key"}`` records, excluding degenerate.

    Args:
        profiles: ``{key: {d: mu(d)}}``.
        shells: ``{key: [shell rows]}``, used to derive ``r_val``.
        max_field: Field name for the shell maximum (``derive_r_val``).
        keys: Which instances to include.

    Returns:
        Non-degenerate instance records (``r_val != 0``).
    """
    out = []
    for k in keys:
        rv = derive_r_val(shells[k], max_field)
        if rv == 0:
            continue  # degenerate: excluded from every summary, per house rule
        out.append({"key": k, "profile": profiles[k], "r_val": rv})
    return out


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------


def main() -> None:
    """Characterise the bias distribution, score candidate grids, and write outputs."""
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", default="results/mean_epsilon")
    p.add_argument("--figdir", default="figures")
    p.add_argument("--tol", type=float, default=1e-6, help="mass-point matching tolerance")
    args = p.parse_args()

    out_dir = ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = ROOT / args.figdir
    fig_dir.mkdir(parents=True, exist_ok=True)

    mvm = load_mean_vs_max(ROOT)
    mfs = load_mean_fullspace(ROOT)
    ft = load_final_table(ROOT)

    # ---- 1. Characterise the distribution -------------------------------

    # Ceiling B(Chat): the largest sample, from final_table_eps_gt1.
    # "Informative" = r_val != 0 (matches the corpus finding that half of
    # these sit exactly at B(Chat) == 1.0); degenerate (r_val == 0) excluded.
    informative = [r for r in ft if r["r_val"] not in (0, None)]
    ceiling_all = np.array([r["beta_top"] for r in informative])
    finite_only = [r for r in informative if r["r_val"] != UNREACHED]
    ceiling_finite = np.array([r["beta_top"] for r in finite_only])

    ceiling_mass = mass_point_report(ceiling_all, args.tol)
    ceiling_mass_finite = mass_point_report(ceiling_finite, args.tol)

    # Shell means mu(d), pooled over every instance and every non-degenerate
    # depth (d >= r_val), from the two exhaustive/near-exhaustive corpora.
    mvm_keys = list(mvm["meta"].keys())
    mvm_insts = build_instance_list(mvm["profiles"], mvm["shells"], "max", mvm_keys)
    mfs_keys = list(mfs["meta"].keys())
    mfs_insts = build_instance_list(mfs["profiles"], mfs["shells"], "retr_max", mfs_keys)

    def pooled_mu(insts: list[dict[str, Any]]) -> np.ndarray:
        vals = []
        for inst in insts:
            for d, v in inst["profile"].items():
                if d >= inst["r_val"] and v is not None:
                    vals.append(v)
        return np.array(vals)

    pooled_mvm = pooled_mu(mvm_insts)
    pooled_mfs = pooled_mu(mfs_insts)
    pooled_mu_mass = mass_point_report(pooled_mfs, args.tol)

    # Individual per-state B: no raw per-state list is stored, but the
    # exhaustive mean_vs_max shells carry min/median/max/std per shell,
    # which is the richest already-computed summary of the per-state
    # distribution available without recomputation.
    state_summary_rows = [r for k in mvm["shells"] for r in mvm["shells"][k] if r["n_states"] > 0]
    state_medians = np.array([r["median_states"] for r in state_summary_rows])
    state_maxes = np.array([r["max"] for r in state_summary_rows])
    state_frac_nonzero = np.array([r["frac_nonzero"] for r in state_summary_rows])

    distribution = {
        "ceiling_B_Chat": {
            "source": "results/final_table_eps_gt1/instances.jsonl, D_LLM, status=ok, r_val != 0",
            "n_informative": len(ceiling_all),
            "n_finite_r_val": len(ceiling_finite),
            "n_unreached_r_val": len(ceiling_all) - len(ceiling_finite),
            "quantiles_informative": {
                str(q): float(np.percentile(ceiling_all, q)) for q in (0, 10, 25, 50, 75, 90, 100)
            },
            "mass_points_informative": ceiling_mass,
            "mass_points_finite_only": ceiling_mass_finite,
        },
        "shell_means_mu_d": {
            "source_exhaustive_15": "results/mean_vs_max/shells.jsonl",
            "source_79": "results/mean_fullspace/shells.jsonl",
            "n_pooled_exhaustive_15": len(pooled_mvm),
            "n_pooled_79": len(pooled_mfs),
            "quantiles_pooled_79": {
                str(q): float(np.percentile(pooled_mfs, q)) for q in (0, 10, 25, 50, 75, 90, 95, 99, 100)
            },
            "quantiles_pooled_exhaustive_15": {
                str(q): float(np.percentile(pooled_mvm, q)) for q in (0, 10, 25, 50, 75, 90, 100)
            },
            "mass_points_pooled_79": pooled_mu_mass,
        },
        "individual_state_B_proxy": {
            "note": (
                "No raw per-state B list is stored anywhere on disk; these are "
                "shell-level min/median/max/frac_nonzero fields from the 15 "
                "exhaustive real instances (results/mean_vs_max/shells.jsonl), "
                "pooled over every shell with at least one state."
            ),
            "n_shells_pooled": len(state_summary_rows),
            "median_states_quantiles": {
                str(q): float(np.percentile(state_medians, q)) for q in (0, 25, 50, 75, 100)
            },
            "max_quantiles": {str(q): float(np.percentile(state_maxes, q)) for q in (0, 25, 50, 75, 100)},
            "frac_nonzero_quantiles": {
                str(q): float(np.percentile(state_frac_nonzero, q)) for q in (0, 25, 50, 75, 100)
            },
        },
    }
    (out_dir / "mass_points.json").write_text(json.dumps(distribution, indent=2) + "\n")

    # ---- 2 & 3. Score candidate grids on several corpora -----------------

    quantile_grid = quantile_derived_grid(pooled_mfs)
    log_grid = log_spaced_grid()
    grids = {
        "legacy_max_grid": LEGACY_GRID,
        "log_spaced": log_grid,
        "quantile_derived_79": quantile_grid,
        "round_number_proposal": ROUND_NUMBER_GRID,
    }

    confirmed_mass_points = [1.00]  # see mass_points.json: the only value

    # median split on |K_G0|
    nk = {k: mfs["meta"][k]["n_knowledge"] for k in mfs["meta"]}
    median_nk = float(np.median(list(nk.values())))

    def subset(keys_source, pred):
        keys = [inst["key"] for inst in mfs_insts if pred(inst["key"])]
        return build_instance_list(mfs["profiles"], mfs["shells"], "retr_max", keys)

    corpora = {
        "mean_fullspace_all_79": mfs_insts,
        "mean_fullspace_real_10": subset(mfs_insts, lambda k: mfs["meta"][k]["source"] == "real"),
        "mean_fullspace_synthetic_69": subset(mfs_insts, lambda k: mfs["meta"][k]["source"] == "synthetic"),
        "mean_fullspace_small_K": subset(mfs_insts, lambda k: nk[k] <= median_nk),
        "mean_fullspace_large_K": subset(mfs_insts, lambda k: nk[k] > median_nk),
        "mean_vs_max_exhaustive_real_15": mvm_insts,
    }

    rows_out = []
    for corpus_name, insts in corpora.items():
        for grid_name, grid in grids.items():
            score = score_grid(insts, grid)
            dq = disqualified_grid_points(grid, confirmed_mass_points, args.tol)
            rows_out.append(
                {
                    "corpus": corpus_name,
                    "grid": grid_name,
                    "grid_values": ",".join(f"{v:g}" for v in grid),
                    **score,
                    "disqualified": bool(dq),
                    "disqualified_values": ",".join(f"{v:g}" for v in dq),
                }
            )

    with (out_dir / "grid_scores.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()))
        writer.writeheader()
        writer.writerows(rows_out)

    summary = {
        "grids": {k: v for k, v in grids.items()},
        "confirmed_mass_points": confirmed_mass_points,
        "median_n_knowledge_mean_fullspace": median_nk,
        "corpus_sizes": {k: len(v) for k, v in corpora.items()},
        "recommendation": {
            "grid_name": "round_number_proposal",
            "values": ROUND_NUMBER_GRID,
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    for r in rows_out:
        print(
            f"[mean_epsilon_grid] {r['corpus']:32s} {r['grid']:24s} "
            f"avg_distinct={r['avg_distinct_finite']:.2f} "
            f"avg_unreached={r['avg_n_unreached']:.2f} "
            f"avg_pinned={r['avg_n_pinned']:.2f} "
            f"disqualified={r['disqualified']}",
            flush=True,
        )

    # ---- Figures -----------------------------------------------------

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    def _plot_capped(ax, data, cap, color, title, xlabel):
        n_over = int(np.sum(data > cap))
        clipped = np.clip(data, 0, cap)
        ax.hist(clipped, bins=np.linspace(0, cap, 41), color=color, alpha=0.85)
        for v in ROUND_NUMBER_GRID:
            if v <= cap:
                ax.axvline(v, color="#c44e52", linestyle="--", linewidth=1)
        suffix = f"\n({n_over} values > {cap:g} not shown, max={float(data.max()):.1f})" if n_over else ""
        ax.set_title(title + suffix, fontsize=9)
        ax.set_xlabel(xlabel)
        ax.set_xlim(0, cap)

    _plot_capped(
        axes[0], ceiling_all, 3.0, "#4c72b0",
        f"Ceiling B(Chat), D_LLM, r_val != 0 (n={len(ceiling_all)})", "B(Chat) [relative units]",
    )
    axes[0].set_ylabel("count")
    _plot_capped(
        axes[1], pooled_mfs, 3.0, "#55a868",
        f"Pooled shell means mu(d), d >= r_val (n={len(pooled_mfs)})", "mu(d) [relative units]",
    )

    fig.suptitle("Proposed grid (dashed red) against the empirical bias distribution")
    fig.tight_layout()
    fig.savefig(fig_dir / "mean_eps_distribution.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    corpus_names = list(corpora.keys())
    grid_names = list(grids.keys())
    x = np.arange(len(corpus_names))
    width = 0.8 / len(grid_names)
    for i, gname in enumerate(grid_names):
        vals = [
            next(r["avg_distinct_finite"] for r in rows_out if r["corpus"] == c and r["grid"] == gname)
            for c in corpus_names
        ]
        ax.bar(x + i * width, vals, width=width, label=gname)
    ax.set_xticks(x + width * (len(grid_names) - 1) / 2)
    ax.set_xticklabels(corpus_names, rotation=30, ha="right")
    ax.set_ylabel("avg distinct finite radii per instance")
    ax.set_title("Grid discrimination by corpus split")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(fig_dir / "mean_eps_grid_scores.png", dpi=150)
    plt.close(fig)

    print(f"[mean_epsilon_grid] wrote {out_dir}/mass_points.json, grid_scores.csv, summary.json")
    print(f"[mean_epsilon_grid] wrote {fig_dir}/mean_eps_distribution.png, mean_eps_grid_scores.png")


if __name__ == "__main__":
    main()
