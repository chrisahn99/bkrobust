"""Publication anchor table: tau_b(predictor, AUC_*_usable) for r_val,
shd_truth (SHD to ground truth) and n_k (|K|), across the nine pre-registered
strata plus the null-cell high-precision re-measurement.

Read-only against every committed input under ``results/axis_robustness/``.
Never regenerates data, never runs a sweep, never writes to an existing file.
Its only outputs are ``table_tau_comparisons.md`` (repo root) and
``results/axis_robustness/table_tau_source.csv``, both created fresh by the
caller (``scripts`` below), not by this module itself -- this module only
computes and returns in-memory frames plus a markdown string, so it can be
imported and its output diffed twice for the determinism check (L3) without
touching the filesystem more than the two intended writes.

Design:

- flip arm, 6 strata (coverage in {0.5, 1.0}) x (base_wrongness in
  {0.0, 0.10, 0.25}): endpoint is ``AUC_frac_usable``, recomputed here from
  ``survival_curves.csv`` + ``survival_instances.csv`` via the exact same
  nearest-defined-depth / tie-break-toward-smaller-depth algorithm as
  session 7's sensitivity analysis
  (:func:`bkrobust.robustness.survival.auc_frac`, restricted to depths with
  ``n_eval >= 30``) -- see PREREGISTRATION.md Appendix F.3. This is a
  read-only recomputation from source, never a copy of
  ``analysis_tau_primary.csv``'s numbers (those are used only for the L1
  cross-check in the caller).

- tiered arm, 3 strata (n_tiers in {2, 3, 4}): endpoint is ``AUC_rate_usable``,
  computed here for the first time (the gap flagged by the task) as the mean
  of ``S_rate_*`` in ``analysis_tiered_rebinned.csv`` restricted to rates
  with ``n_eval_rate_* >= 30`` -- no nearest-neighbour infill, unlike that
  file's own ``AUC_frac_rate`` column. An instance with zero surviving rates
  is excluded and counted, never silently dropped to NaN-as-0 or infilled.

- null-cell supplementary row: flip cov=0.5 bw=0.25 at N=1000 draws
  (``null_instances.csv`` / its ``AUC_frac_usable_n1000`` column), a
  higher-precision re-measurement of the sixth flip stratum already in the
  main table -- not a tenth stratum.

Kendall tau-b, never tau-a. Bootstrap: 10,000 resamples over instances,
``np.random.default_rng(0)``, 2.5/97.5 percentiles -- see
:func:`bootstrap_tau_ci`. Gate is ``fast_gate`` exclusively (checked, not
assumed, in :func:`check_gate_is_fast_gate_only`).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

from bkrobust.core.conventions import UNREACHED
from bkrobust.robustness.survival import auc_frac as _auc_frac_fn

# ---------------------------------------------------------------------------
# Constants fixed by the brief
# ---------------------------------------------------------------------------

USABLE_MIN_N_EVAL = 30
N_BOOT = 10_000
BOOT_SEED = 0
UNDERPOWERED_N = 15
PREDICTORS = ["r_val", "shd_truth", "n_k"]

TIERED_RATE_GRID = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]

FLIP_STRATA = [(cov, bw) for cov in (0.5, 1.0) for bw in (0.00, 0.10, 0.25)]
TIERED_STRATA = [2, 3, 4]

R_ASSUMES_FOOTNOTE = (
    "Conjecture 2 (hence Anti-Exchange Case B, verified not proved); the error "
    "is one-sided, so radii can only be too large."
)
GATE_FOOTNOTE = "benchmarks.measure.fast_gate exclusively."
SHD_UNDEFINED_NOTE = "shd_truth ≡ 0 by construction there"


# ---------------------------------------------------------------------------
# Loading (read-only)
# ---------------------------------------------------------------------------


def load_inputs(root: str) -> dict[str, pd.DataFrame]:
    p = lambda name: f"{root}/{name}"
    return {
        "instances": pd.read_csv(p("survival_instances.csv")).sort_values(
            "instance_id", kind="stable"
        ).reset_index(drop=True),
        "curves": pd.read_csv(p("survival_curves.csv")).sort_values(
            ["instance_id", "d"], kind="stable"
        ).reset_index(drop=True),
        "tiered_rebinned": pd.read_csv(p("analysis_tiered_rebinned.csv")).sort_values(
            "instance_id", kind="stable"
        ).reset_index(drop=True),
        "null_instances": pd.read_csv(p("null_instances.csv")).sort_values(
            "instance_id", kind="stable"
        ).reset_index(drop=True),
        "tau_primary": pd.read_csv(p("analysis_tau_primary.csv")),
    }


def check_gate_is_fast_gate_only(instances: pd.DataFrame) -> None:
    vals = set(instances["gate"].dropna().unique())
    if vals != {"fast_gate"}:
        raise ValueError(f"gate is not fast_gate-exclusive: {vals}")


def get_r_assumes(instances: pd.DataFrame) -> str:
    vals = instances["r_assumes"].dropna().unique()
    if len(vals) != 1:
        raise ValueError(f"r_assumes not constant across dataset: {vals!r}")
    return str(vals[0])


# ---------------------------------------------------------------------------
# flip arm: AUC_frac_usable, recomputed from survival_curves.csv
# ---------------------------------------------------------------------------


def compute_auc_frac_usable(
    instances: pd.DataFrame, curves: pd.DataFrame, *, min_n_eval: int = USABLE_MIN_N_EVAL
) -> pd.DataFrame:
    """Per flip-instance AUC_frac_usable: mean S on the fixed fractional grid
    (0.1..1.0 of n_k), restricted to depths with n_eval >= min_n_eval, via the
    exact fractional-grid nearest-defined-depth algorithm
    (:func:`bkrobust.robustness.survival.auc_frac`). Instances with zero
    usable depths get NaN, never a silent 0."""
    nk_map = dict(zip(instances["instance_id"], instances["n_k"]))
    rows = []
    for iid, g in curves.groupby("instance_id", sort=True):
        nk = nk_map.get(iid)
        curve: dict[int, dict[str, Any]] = {}
        for _, r in g.iterrows():
            if r["n_eval"] >= min_n_eval and pd.notna(r["S"]):
                curve[int(r["d"])] = {"S": float(r["S"])}
        n_usable = len(curve)
        if nk is None or not (nk > 0) or n_usable == 0:
            val = float("nan")
        else:
            raw = _auc_frac_fn(curve, int(nk))
            val = float(raw) if raw != "" else float("nan")
        rows.append({"instance_id": iid, "AUC_frac_usable": val, "n_usable_depths": n_usable})
    return pd.DataFrame(rows).sort_values("instance_id", kind="stable").reset_index(drop=True)


# ---------------------------------------------------------------------------
# tiered arm: AUC_rate_usable -- the gap this task closes
# ---------------------------------------------------------------------------


@dataclass
class TieredUsableResult:
    table: pd.DataFrame  # instance_id, AUC_rate_usable, n_rates_usable
    n_excluded_no_usable_rate: int
    rates_surviving_per_instance: list[int] = field(default_factory=list)


def compute_auc_rate_usable(
    tiered_rebinned: pd.DataFrame, *, min_n_eval: int = USABLE_MIN_N_EVAL
) -> TieredUsableResult:
    """AUC_rate_usable per tiered instance: the mean of S_rate_* over only
    those corruption rates whose n_eval_rate_* >= min_n_eval. No
    nearest-neighbour infill (unlike ``AUC_frac_rate`` in the source file).
    An instance with zero surviving rates is excluded and counted."""
    rows = []
    n_excluded = 0
    counts = []
    for _, r in tiered_rebinned.iterrows():
        s_vals = []
        for rate in TIERED_RATE_GRID:
            n_eval = r[f"n_eval_rate_{rate:.2f}"]
            s = r[f"S_rate_{rate:.2f}"]
            if pd.notna(n_eval) and n_eval >= min_n_eval and pd.notna(s):
                s_vals.append(float(s))
        n_usable = len(s_vals)
        counts.append(n_usable)
        if n_usable == 0:
            n_excluded += 1
            auc = float("nan")
        else:
            auc = float(np.mean(s_vals))
        rows.append(
            {"instance_id": r["instance_id"], "AUC_rate_usable": auc, "n_rates_usable": n_usable}
        )
    table = pd.DataFrame(rows).sort_values("instance_id", kind="stable").reset_index(drop=True)
    return TieredUsableResult(
        table=table, n_excluded_no_usable_rate=n_excluded, rates_surviving_per_instance=counts
    )


# ---------------------------------------------------------------------------
# Kendall tau-b with bootstrap CI over instances
# ---------------------------------------------------------------------------


def bootstrap_tau_ci(
    x: np.ndarray, y: np.ndarray, *, n_boot: int = N_BOOT, seed: int = BOOT_SEED
) -> tuple[float, float, int]:
    """10,000 resamples over instances (rows), seed 0, explicit Generator.
    Returns (lo2.5, hi97.5, n_nan_resamples)."""
    rng = np.random.default_rng(seed)
    n = len(x)
    taus = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        t, _ = kendalltau(x[idx], y[idx], variant="b", nan_policy="propagate")
        taus[i] = t
    n_nan = int(np.isnan(taus).sum())
    valid = taus[~np.isnan(taus)]
    if len(valid) == 0:
        return float("nan"), float("nan"), n_nan
    lo, hi = np.percentile(valid, [2.5, 97.5])
    return float(lo), float(hi), n_nan


def tau_row(
    df: pd.DataFrame, predictor: str, endpoint: str, *, r_assumes: str,
    n_boot: int = N_BOOT, seed: int = BOOT_SEED,
) -> dict[str, Any]:
    """One stratum x one predictor, against the given endpoint column
    (already attached to df). Mirrors PREREGISTRATION.md section 7 exactly:
    Kendall tau-b (never tau-a), instance-level bootstrap, undefined
    (predictor constant) reported explicitly rather than as tau=0."""
    total_n = len(df)
    n_excluded_unreached = int((df["r_val"] == UNREACHED).sum()) if predictor == "r_val" else 0

    sub = df[[predictor, endpoint]].copy()
    if predictor == "r_val":
        sub = sub[sub["r_val"] != UNREACHED]
    sub = sub.dropna(subset=[predictor, endpoint])
    n = len(sub)
    n_excluded_nan = max(total_n - n_excluded_unreached - n, 0)

    row: dict[str, Any] = {
        "predictor": predictor,
        "endpoint": endpoint,
        "stratum_n_instances": total_n,
        "n": n,
        "n_excluded_unreached": n_excluded_unreached,
        "n_excluded_nan": n_excluded_nan,
        "predictor_constant": False,
        "tau_b": float("nan"),
        "p_value": float("nan"),
        "ci_lo_2p5": float("nan"),
        "ci_hi_97p5": float("nan"),
        "n_boot_nan": 0,
        "underpowered": n < UNDERPOWERED_N,
        "status": "ok",
        "r_assumes": r_assumes if predictor == "r_val" else "",
    }

    if n < 2:
        row["status"] = "insufficient_n"
        return row

    x = sub[predictor].to_numpy(dtype=float)
    y = sub[endpoint].to_numpy(dtype=float)

    if np.unique(x).size <= 1:
        row["predictor_constant"] = True
        row["status"] = "undefined (predictor constant)"
        return row
    if np.unique(y).size <= 1:
        row["status"] = "undefined (endpoint constant)"
        return row

    tau, pval = kendalltau(x, y, variant="b", nan_policy="propagate")
    row["tau_b"] = float(tau)
    row["p_value"] = float(pval)
    lo, hi, n_boot_nan = bootstrap_tau_ci(x, y, n_boot=n_boot, seed=seed)
    row["ci_lo_2p5"] = lo
    row["ci_hi_97p5"] = hi
    row["n_boot_nan"] = n_boot_nan
    return row


# ---------------------------------------------------------------------------
# Assembling the nine strata + supplementary null-cell row
# ---------------------------------------------------------------------------


def build_source_table(inputs: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Returns (source_table, diagnostics). ``source_table`` has one row per
    stratum x predictor (27 rows for the 9 strata, 3 more for the
    supplementary null-cell re-run predictors -- 30 total), matching
    ``table_tau_source.csv``. ``diagnostics`` carries the L2 filter stats and
    the raw per-stratum AUC frames for the L1 cross-check."""
    instances = inputs["instances"]
    check_gate_is_fast_gate_only(instances)
    r_assumes = get_r_assumes(instances)

    flip = instances[instances["arm"] == "flip"].copy()
    tiered = instances[instances["arm"] == "tiered"].copy()

    # --- flip: AUC_frac_usable, recomputed from source -----------------
    flip_usable = compute_auc_frac_usable(flip, inputs["curves"][inputs["curves"]["instance_id"].isin(flip["instance_id"])])
    flip = flip.merge(flip_usable, on="instance_id", how="left", validate="one_to_one")

    # --- tiered: AUC_rate_usable, the gap this task closes ---------------
    tiered_result = compute_auc_rate_usable(inputs["tiered_rebinned"])
    tiered = tiered.merge(tiered_result.table, on="instance_id", how="left", validate="one_to_one")

    rows = []
    for cov, bw in FLIP_STRATA:
        stratum_df = flip[(flip["coverage"] == cov) & (np.isclose(flip["base_wrongness"], bw))]
        for predictor in PREDICTORS:
            r = tau_row(stratum_df, predictor, "AUC_frac_usable", r_assumes=r_assumes)
            r.update({"arm": "flip", "coverage": cov, "base_wrongness": bw, "n_tiers": None})
            rows.append(r)

    for nt in TIERED_STRATA:
        stratum_df = tiered[tiered["n_tiers"] == nt]
        for predictor in PREDICTORS:
            r = tau_row(stratum_df, predictor, "AUC_rate_usable", r_assumes=r_assumes)
            r.update({"arm": "tiered", "coverage": None, "base_wrongness": None, "n_tiers": nt})
            rows.append(r)

    # --- supplementary: null-cell re-run, N=1000, flip cov=0.5 bw=0.25 ---
    null_inst = inputs["null_instances"].copy()
    null_inst = null_inst.rename(columns={"AUC_frac_usable_n1000": "AUC_frac_usable"})
    for predictor in PREDICTORS:
        r = tau_row(null_inst, predictor, "AUC_frac_usable", r_assumes=r_assumes)
        r.update(
            {
                "arm": "flip_null_rerun",
                "coverage": 0.5,
                "base_wrongness": 0.25,
                "n_tiers": None,
            }
        )
        rows.append(r)

    source = pd.DataFrame(rows)
    cols = [
        "arm", "coverage", "base_wrongness", "n_tiers", "predictor", "endpoint",
        "stratum_n_instances", "n", "n_excluded_unreached", "n_excluded_nan",
        "predictor_constant", "underpowered", "status", "tau_b", "p_value",
        "ci_lo_2p5", "ci_hi_97p5", "n_boot_nan", "r_assumes",
    ]
    source = source[cols]

    diagnostics = {
        "r_assumes": r_assumes,
        "flip_frame": flip,
        "tiered_frame": tiered,
        "tiered_result": tiered_result,
        "null_inst": null_inst,
    }
    return source, diagnostics


# ---------------------------------------------------------------------------
# Presentation layer -- stratum labels, bolding rule, markdown rendering.
# Computed once here and reused verbatim for both outputs, so the CSV and the
# markdown table can never silently disagree (L5).
# ---------------------------------------------------------------------------


def _stratum_label(row: pd.Series) -> str:
    if row["arm"] == "flip":
        return f"flip, coverage={row['coverage']:.1f}, base_wrongness={row['base_wrongness']:.2f}"
    if row["arm"] == "tiered":
        return f"tiered, n_tiers={int(row['n_tiers'])}"
    if row["arm"] == "flip_null_rerun":
        return "flip, coverage=0.5, base_wrongness=0.25 -- null-cell re-run, N=1000"
    raise ValueError(row["arm"])


def _excludes_zero(lo: float, hi: float) -> bool:
    if pd.isna(lo) or pd.isna(hi):
        return False
    return (lo > 0) or (hi < 0)


def add_presentation_columns(source: pd.DataFrame) -> pd.DataFrame:
    """Adds ``stratum_label``, ``is_supplementary`` and ``winner_bold``
    (True only for the sole predictor in a stratum whose CI excludes zero
    while the other two do not -- never for a null, never a tie) so the CSV
    carries exactly the bolding decision the markdown table renders."""
    out = source.copy()
    out["stratum_label"] = out.apply(_stratum_label, axis=1)
    out["is_supplementary"] = out["arm"] == "flip_null_rerun"
    out["excludes_zero"] = [
        _excludes_zero(lo, hi) for lo, hi in zip(out["ci_lo_2p5"], out["ci_hi_97p5"])
    ]

    winner_bold = []
    group_cols = ["arm", "coverage", "base_wrongness", "n_tiers"]
    for _, g in out.groupby(group_cols, dropna=False, sort=False):
        winners = g.loc[g["excludes_zero"], "predictor"].tolist()
        sole_winner = winners[0] if len(winners) == 1 else None
        winner_bold.extend(out.loc[g.index, "predictor"] == sole_winner)
    out["winner_bold"] = winner_bold
    return out


def _fmt_cell(row: pd.Series) -> str:
    if row["status"] == "undefined (predictor constant)":
        return "undefined (predictor constant)[^shdconst]"
    s = f"{row['tau_b']:.3f} [{row['ci_lo_2p5']:.3f}, {row['ci_hi_97p5']:.3f}]"
    return f"**{s}**" if row["winner_bold"] else s


MAIN_STRATA_ORDER = [
    ("flip", 0.5, 0.00, None),
    ("flip", 0.5, 0.10, None),
    ("flip", 0.5, 0.25, None),
    ("flip", 1.0, 0.00, None),
    ("flip", 1.0, 0.10, None),
    ("flip", 1.0, 0.25, None),
    ("tiered", None, None, 2),
    ("tiered", None, None, 3),
    ("tiered", None, None, 4),
]


def render_table_rows(annotated: pd.DataFrame) -> tuple[list[str], str | None]:
    """Builds the nine main markdown table rows plus the supplementary row.
    Returns (main_rows_markdown, supplementary_row_markdown)."""

    def row_for(mask: pd.Series) -> str:
        g = annotated[mask].set_index("predictor")
        n = int(g["n"].iloc[0])
        endpoint = str(g["endpoint"].iloc[0])
        label = str(g["stratum_label"].iloc[0])
        underpowered_flag = " **[underpowered, n<15]**" if n < UNDERPOWERED_N else ""
        cells = [_fmt_cell(g.loc[p]) for p in PREDICTORS]
        return f"| {label}{underpowered_flag} | {n} | `{endpoint}` | " + " | ".join(cells) + " |"

    main_rows = []
    for arm, cov, bw, nt in MAIN_STRATA_ORDER:
        if arm == "flip":
            mask = (
                (annotated["arm"] == arm)
                & (annotated["coverage"] == cov)
                & (np.isclose(annotated["base_wrongness"], bw))
            )
        else:
            mask = (annotated["arm"] == arm) & (annotated["n_tiers"] == nt)
        main_rows.append(row_for(mask))

    supp_mask = annotated["arm"] == "flip_null_rerun"
    supp_row = row_for(supp_mask) if supp_mask.any() else None
    return main_rows, supp_row


def render_markdown(annotated: pd.DataFrame) -> str:
    header = (
        "| stratum | n | endpoint used | τ_b(r_val) [95% CI] | τ_b(shd_truth) [95% CI] | "
        "τ_b(n_k) [95% CI] |"
    )
    sep = "|---|---|---|---|---|---|"
    main_rows, supp_row = render_table_rows(annotated)

    intro = (
        "This table reports τ_b between three predictors (`r_val`, `shd_truth`, `n_k`) and "
        "the conservative-control survival AUC (grid points with `n_eval ≥ 30`) across the "
        "nine pre-registered strata, using `AUC_frac_usable` for the flip arm (targeted-depth "
        "`d_claims` axis) and the newly computed `AUC_rate_usable` for the tiered arm "
        "(`corruption_rate` axis, per Appendix E) -- the two arms' endpoints are on different "
        "axes and are never compared row-to-row. `r_val`'s CI excludes zero in 8 of 9 strata "
        "and is the sole predictor to do so (bolded) only in flip cov=0.5 bw=0.10; elsewhere "
        "it shares significance with, or is beaten by, the SHD/`n_k` baselines, or (flip "
        "cov=1.0 bw=0.00) all three are indistinguishable from zero. `n_k`'s association is "
        "negative rather than absent: in every tiered stratum (n_tiers 2/3/4, τ ≈ −0.42/−0.38/"
        "−0.38) and in the higher-precision null-cell re-run below, its 95% CI excludes zero "
        "on the negative side, so H4's \"inert baseline\" framing is wrong. One stratum, flip "
        "cov=0.5 bw=0.25, carries a confirmed null for `r_val` at high precision: the N=1000 "
        "re-run (Appendix H) gives τ = +0.032, essentially unchanged from τ = +0.031 at N=200 "
        "unfiltered (five times the draws moved it +0.001), even though that same stratum's "
        "N=200 `AUC_frac_usable` value used in the main table row above is a markedly higher "
        "+0.279 -- a usable-depth-filter artifact at low draw count, not a draw-count "
        "attenuation effect (see the source-data note in the accompanying report). "
        "`shd_truth` is undefined by construction (predictor constant) in exactly one "
        "stratum, flip cov=1.0 bw=0.00, where full coverage and zero base wrongness make "
        "`G₀` the ground truth for every instance."
    )

    lines = [intro, "", header, sep] + main_rows
    lines += ["", "**Supplementary -- null-cell high-precision re-measurement**", "", header, sep]
    if supp_row:
        lines.append(supp_row)
    lines += [
        "",
        "[^shdconst]: `shd_truth ≡ 0` by construction there (full coverage, zero base "
        "wrongness => `G₀` *is* the ground-truth DAG for every instance in this stratum).",
        "",
        f"**Assumption (every `r_val` row):** {R_ASSUMES_FOOTNOTE}",
        "",
        f"**Gate:** {GATE_FOOTNOTE}",
    ]
    return "\n".join(lines)


def render_document(annotated: pd.DataFrame, diag: dict[str, Any]) -> str:
    title = "# Publication anchor table -- τ_b(predictor, survival AUC) by stratum\n"
    return title + "\n" + render_markdown(annotated) + "\n\n" + render_notes(diag) + "\n"


def render_notes(diag: dict[str, Any]) -> str:
    """L2 filter-statistics note (tiered AUC_rate_usable) and the source-data
    anomaly note for the flip cov=0.5 bw=0.25 null cell, both referenced from
    the intro paragraph above."""
    counts = np.array(diag["tiered_result"].rates_surviving_per_instance)
    n_excluded = diag["tiered_result"].n_excluded_no_usable_rate
    n_total = len(counts)
    lines = [
        "## Notes",
        "",
        "### Tiered `AUC_rate_usable` filter statistics (the gap this table closes)",
        "",
        f"Rates surviving the `n_eval_rate_* >= {USABLE_MIN_N_EVAL}` filter, per instance, "
        f"out of the {len(TIERED_RATE_GRID)}-point `corruption_rate` grid "
        f"({n_total} tiered instances total):",
        "",
        f"- min: {int(counts.min())}",
        f"- median: {float(np.median(counts)):.1f}",
        f"- max: {int(counts.max())}",
        f"- instances excluded for zero usable rates: {n_excluded} (of {n_total})",
        "",
        "### Source-data note: the N=200 vs. N=1000 `AUC_frac_usable` discrepancy at "
        "flip cov=0.5 bw=0.25",
        "",
        "The main table row above uses session 7's original N=200-draw survival data, "
        "filtered to `n_eval >= 30`, per the task's instruction to use the conservative "
        "usable control everywhere. For this specific stratum that gives τ_b(`r_val`, "
        "`AUC_frac_usable`) = +0.279 (CI excludes zero). The Appendix H null-cell re-run, "
        "same 220 instances, N=1000 draws, same `n_eval >= 30` filter, gives τ_b = +0.032 "
        "(CI includes zero) -- essentially the *unfiltered* N=200 value (+0.031, "
        "PREREGISTRATION.md Appendix H.1) rather than the filtered N=200 value. This is "
        "not a contradiction in the source files -- both numbers are reproduced exactly "
        "from their respective committed CSVs (see L1 in the accompanying report) -- but it "
        "is worth flagging: at N=200 draws the `n_eval >= 30` filter is aggressive in this "
        "thin, high-contradiction-rate stratum, and which depths it excludes correlates with "
        "`r_val` strongly enough to move τ from +0.03 to +0.28. At N=1000 draws the same "
        "threshold is far less exclusionary (raw counts are 5x higher), the selection effect "
        "attenuates, and τ reverts to the unfiltered value. The N=1000 result is the more "
        "trustworthy one for this stratum precisely because it is less exposed to that "
        "selection effect -- which is the whole point of the higher-precision re-run -- but "
        "the main table still reports the N=200 usable-filtered figure per the task's fixed "
        "rule of using `AUC_*_usable` for every main-table cell, with this note attached "
        "rather than silently overriding that rule.",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Entry point. Only output paths this module ever writes to:
# ``<repo_root>/table_tau_comparisons.md`` and
# ``results/axis_robustness/table_tau_source.csv`` -- both new files, never
# an existing one, and never anything under results/axis_robustness/ via
# core.resultsio.write_manifest (not imported, not called, here or anywhere
# above).
# ---------------------------------------------------------------------------


def build(results_dir: str) -> tuple[pd.DataFrame, str]:
    """Pure build: reads only, returns (annotated_source_table, markdown_doc).
    Writing is the caller's job, so this can be invoked twice for the
    determinism check without touching the filesystem."""
    inputs = load_inputs(results_dir)
    source, diag = build_source_table(inputs)
    annotated = add_presentation_columns(source)
    doc = render_document(annotated, diag)
    return annotated, doc


def main(results_dir: str, repo_root: str) -> None:
    import os

    annotated, doc = build(results_dir)
    md_path = f"{repo_root}/table_tau_comparisons.md"
    csv_path = f"{results_dir}/table_tau_source.csv"
    if os.path.exists(csv_path):
        raise FileExistsError(f"refusing to overwrite existing file: {csv_path}")
    with open(md_path, "x", encoding="utf-8") as f:
        f.write(doc)
    annotated.to_csv(csv_path, index=False)


if __name__ == "__main__":
    import sys as _sys

    _root = _sys.argv[1] if len(_sys.argv) > 1 else "."
    main(f"{_root}/results/axis_robustness", _root)
