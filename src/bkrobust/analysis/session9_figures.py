"""Figures for session 9 (Axis Robustness, real-corpus sweep).

Every figure is built exclusively from committed files under
``results/axis_robustness_real/`` (the real-network sweep; never mixed with any
synthetic-ensemble output from earlier sessions). At the time this module was
written the underlying sweep was still finishing (4 shards of one network still
outstanding), so nothing here reports, quotes, or interprets a specific
numeric value from the data — that is left to a human once the sweep is final.
This module is only the figure-drawing code; re-run it (``python -m
bkrobust.analysis.session9_figures``) after the sweep completes and every
number on every plot updates on its own.

Two arms, never merged onto one axis, same house rule as session 7. **flip**
corrupts by reversing a targeted count of asserted claims, so its native knob
is an integer depth ``d`` (``grid_kind = "d_claims"``). **tiered** corrupts by
relocating nodes between tiers, so its native knob is ``corruption_rate``
(``grid_kind = "corruption_rate"``) — a different quantity in different units.
The one deliberate exception is F6, whose cross-arm design (``results/
axis_robustness_real/xarm_*.csv``) already places both arms on a shared
*achieved-intensity* axis by construction; that is a different, purpose-built
comparison, not a violation of the flip/tiered axis rule.

Radius here is ``r_hop`` (the ``radius`` column of ``analysis_units.csv`` /
``analysis_tau.csv``), and every figure that plots it carries the one-sided
error note verbatim (see ``RADIUS_NOTE`` below). The status sentinel
``UNREACHED`` (``radius = -1``, PREREGISTRATION.md section 5.3) is a status,
never a radius: it is excluded from every numeric radius axis or bucket, and
the exclusion is counted and reported in the figure rather than silently
dropped or plotted at some numeric value.

F2 needs per-(unit, grid point) survival, which is not in any of the
pre-aggregated analysis_*.csv files (those carry only whole-curve endpoints
such as ``AUC_frac``/``AUC_rate``); ``analysis_effect_sizes.csv`` was checked
and confirmed not to carry a grid-point column, so F2 reads the shard cell
files (``results/axis_robustness_real/shards/*.cells.jsonl``) and joins them
to ``analysis_units.csv`` on the ``(shard_id, frame_row_id)`` key that both
files already carry (this is more precise than the network/x/y/stratum join
and avoids any ambiguity from replicate shard files).

Colour-blind safe Okabe-Ito palette throughout, shared with sessions 4-7;
every categorical encoding is doubled with a distinct marker shape or hatch so
nothing depends on colour alone in greyscale or print.
"""

# ruff: noqa: RUF001
# Figure text is display copy; the minus signs and en dashes are intentional.

from __future__ import annotations

import csv
import glob
import json
import re
import statistics as st
import textwrap
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

RES = Path("results/axis_robustness_real")
SHARDS = RES / "shards"
FIGDIR = Path("figures")

# Okabe-Ito, colour-blind safe; fixed assignment, never cycled per-figure.
C_BLACK = "#000000"
C_ORANGE = "#E69F00"
C_SKY = "#56B4E9"
C_GREEN = "#009E73"
C_BLUE = "#0072B2"
C_VERM = "#D55E00"
C_PINK = "#CC79A7"
C_GREY = "#888888"

N_EVAL_MIN = 30  # never plot a point backed by fewer evaluable draws than this.
N_MATCHED_MIN = 15  # xarm bin floor for "carries the comparison" (F6).

RADIUS_NOTE = (
    "r_hop assumes Conjecture 2 (proved: Anti-Exchange Case B, THEOREMS.md section 4); radii are exact."
)


def _style() -> None:
    """Apply the shared house rcParams (see ``session7_figures._style``)."""
    plt.rcParams.update(
        {
            "figure.dpi": 110,
            "savefig.bbox": "tight",
            "axes.grid": True,
            "grid.alpha": 0.3,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "font.size": 9,
        }
    )


def _save(fig: plt.Figure, name: str) -> list[Path]:
    """Write ``fig`` to ``figures/<name>.png`` (150 dpi) and ``.pdf``.

    Args:
        fig: The figure to write. It is closed afterwards.
        name: Basename without extension.

    Returns:
        The paths written, in the order png then pdf.
    """
    FIGDIR.mkdir(parents=True, exist_ok=True)
    out: list[Path] = []
    for ext in ("png", "pdf"):
        p = FIGDIR / f"{name}.{ext}"
        fig.savefig(p, dpi=150 if ext == "png" else None)
        out.append(p)
    plt.close(fig)
    return out


def _read_csv(name: str) -> list[dict]:
    """Read a CSV under ``results/axis_robustness_real`` into a list of dict rows."""
    with (RES / name).open(newline="") as fh:
        return list(csv.DictReader(fh))


def _take(fig: plt.Figure, text: str, y: float = -0.16) -> None:
    """Draw the required one-line "take from this" annotation inside the figure."""
    fig.text(0.01, y, f"Take: {text}", ha="left", va="top", fontsize=8, style="italic", color=C_BLACK)


def _marker_sizes(ns: list[int], lo: float = 18.0, hi: float = 160.0) -> list[float]:
    """Map contributing-unit counts to marker areas, sqrt-scaled, for one panel.

    Args:
        ns: Contributing-unit (or contributing-row) counts, one per point.
        lo: Marker area at the smallest ``n`` in the panel.
        hi: Marker area at the largest ``n`` in the panel.

    Returns:
        Marker areas (``s=`` for ``scatter``), same length as ``ns``.
    """
    if not ns:
        return []
    root = [n**0.5 for n in ns]
    rmin, rmax = min(root), max(root)
    if rmax == rmin:
        return [(lo + hi) / 2 for _ in ns]
    return [lo + (hi - lo) * (r - rmin) / (rmax - rmin) for r in root]


def _bucket5(value: str) -> str:
    """Bucket a non-negative integer-valued radius into ``1..4, "5+"``.

    Values ``<= 0`` (``UNREACHED = -1`` or otherwise undefined) are not a
    bucket this function will return for — callers must filter those out
    before calling it; see ``PREREGISTRATION.md`` section 5.3.
    """
    v = int(float(value))
    return str(v) if v < 5 else "5+"


_BUCKET_ORDER = ["1", "2", "3", "4", "5+"]
_BUCKET_COLOUR = {"1": C_ORANGE, "2": C_SKY, "3": C_GREEN, "4": C_VERM, "5+": C_PINK}
_BUCKET_MARKER = {"1": "o", "2": "s", "3": "^", "4": "D", "5+": "v"}


def _gapped_xy(
    points: dict[float, float], step: float, tol: float = 1e-9
) -> tuple[list[list[float]], list[list[float]]]:
    """Split ``{x: y}`` into runs of consecutive ``x`` (spacing ``step``).

    A missing x between two present ones must not be bridged by a drawn line —
    that would imply data where there is none. Returns parallel lists of x/y
    runs, sorted by x.
    """
    xs = sorted(points)
    x_runs: list[list[float]] = []
    y_runs: list[list[float]] = []
    for x in xs:
        if x_runs and abs(x - x_runs[-1][-1] - step) < tol:
            x_runs[-1].append(x)
            y_runs[-1].append(points[x])
        else:
            x_runs.append([x])
            y_runs.append([points[x]])
    return x_runs, y_runs


# ---------------------------------------------------------------------------
# Stratum-name parsing, shared by F1 and (for sort order only) elsewhere.
# ---------------------------------------------------------------------------

_FLIP_STRATUM_RE = re.compile(r"^flip_cov(\d+)_bw(\d+)$")
_TIERED_STRATUM_RE = re.compile(r"^tiered_nt(\d+)$")


def _stratum_sort_key(stratum: str) -> tuple:
    m = _FLIP_STRATUM_RE.match(stratum)
    if m:
        return (0, int(m.group(1)), int(m.group(2)))
    m2 = _TIERED_STRATUM_RE.match(stratum)
    if m2:
        return (1, int(m2.group(1)), 0)
    return (2, 0, 0)


def _stratum_display(stratum: str) -> str:
    m = _FLIP_STRATUM_RE.match(stratum)
    if m:
        cov = int(m.group(1)) / 100.0
        bw = int(m.group(2)) / 100.0
        return f"flip  cov={cov:.2f}  bw={bw:.2f}"
    m2 = _TIERED_STRATUM_RE.match(stratum)
    if m2:
        return f"tiered  n_tiers={m2.group(1)}"
    return stratum


# ---------------------------------------------------------------------------
# F1 — tau_b forest plot, nine primary strata, conservative endpoint.
# ---------------------------------------------------------------------------

_F1_PREDICTORS = ["radius", "shd_truth", "n_k"]
_F1_COLOUR = {"radius": C_BLUE, "shd_truth": C_GREEN, "n_k": C_VERM}
_F1_MARKER = {"radius": "o", "shd_truth": "^", "n_k": "s"}


def fig_tau_forest() -> list[Path]:
    """F1 — tau_b(predictor, conservative endpoint) with 95% CI, nine primary strata.

    Reads ``analysis_tau.csv``, restricted to ``stratum_class == "primary"``
    (the nine strata: six flip coverage/base-wrongness cells, three tiered
    n_tiers cells) and the three predictors ``radius``, ``shd_truth``,
    ``n_k``. The endpoint is the conservative one per arm — ``AUC_frac_usable``
    for flip strata, ``AUC_rate_usable`` for tiered — never the two arms'
    endpoints compared against each other.

    A cell whose ``status`` starts with ``"undefined"`` (predictor constant in
    that stratum, e.g. ``shd_truth`` at coverage=1.0/base_wrongness=0.0 where
    ``G0`` *is* the truth) is drawn as text, not a point at tau=0. A cell whose
    leave-one-network-out check flips sign (``loo_verdict_flips == "True"``) is
    drawn with an open marker and called out in the legend, never silently
    treated the same as a stable estimate.
    """
    _style()
    rows = [
        r
        for r in _read_csv("analysis_tau.csv")
        if r["stratum_class"] == "primary" and r["predictor"] in _F1_PREDICTORS
    ]
    conservative: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        want = "AUC_frac_usable" if r["arm"] == "flip" else "AUC_rate_usable"
        if r["endpoint"] == want:
            conservative[r["stratum"]].append(r)

    strata = sorted(conservative, key=_stratum_sort_key)
    if not strata:
        raise ValueError("no primary strata with a conservative-endpoint row found")

    fig, ax = plt.subplots(figsize=(11.5, 6.6))
    n_pred = len(_F1_PREDICTORS)
    band = 0.6
    offsets = {p: (i - (n_pred - 1) / 2) * (band / n_pred) for i, p in enumerate(_F1_PREDICTORS)}

    undefined_handle = None
    open_handle = None
    ytick_labels = []
    for si, stratum in enumerate(strata):
        by_pred = {r["predictor"]: r for r in conservative[stratum]}
        any_row = next(iter(by_pred.values()))
        ytick_labels.append(
            f"{_stratum_display(stratum)}   (n={any_row['n']}, n_networks={any_row['n_networks']})"
        )
        for p in _F1_PREDICTORS:
            row = by_pred.get(p)
            if row is None:
                continue
            y = si + offsets[p]
            colour = _F1_COLOUR[p]
            if row["status"].startswith("undefined"):
                ax.text(
                    0.02, y, "undefined (predictor constant)", fontsize=6.8, color=colour,
                    va="center", ha="left", style="italic",
                )
                undefined_handle = Line2D(
                    [], [], ls="", marker="x", color="black", markersize=8,
                    label="undefined (predictor constant) — not tau=0, not omitted",
                )
                continue
            tau = float(row["tau_b"])
            lo, hi = float(row["ci_lo_2p5"]), float(row["ci_hi_97p5"])
            ax.plot([lo, hi], [y, y], color=colour, lw=1.5, zorder=2)
            flips = row.get("loo_verdict_flips") == "True"
            if flips:
                ax.scatter(
                    [tau], [y], s=64, facecolors="none", edgecolors=colour, marker=_F1_MARKER[p],
                    linewidth=1.8, zorder=4,
                )
                open_handle = Line2D(
                    [], [], ls="", marker="o", markerfacecolor="none", markeredgecolor="black",
                    markersize=8, label="leave-one-network-out flips sign",
                )
            else:
                ax.scatter(
                    [tau], [y], s=48, color=colour, marker=_F1_MARKER[p],
                    edgecolor="white", linewidth=0.5, zorder=3,
                )

    ax.axvline(0.0, color=C_GREY, ls="--", lw=1.1, zorder=1)
    ax.set_yticks(range(len(strata)))
    ax.set_yticklabels(ytick_labels, fontsize=8.4)
    ax.invert_yaxis()
    ax.set_xlabel("tau_b(predictor, conservative endpoint), 95% CI  (flip: AUC_frac_usable; tiered: AUC_rate_usable)")
    ax.grid(axis="y", visible=False)
    ax.set_title("Primary strata: tau_b by predictor, with LOO-network sensitivity flagged", fontsize=11.5)

    handles = [
        Line2D([], [], color=_F1_COLOUR[p], marker=_F1_MARKER[p], lw=1.5, label=p) for p in _F1_PREDICTORS
    ]
    if undefined_handle is not None:
        handles.append(undefined_handle)
    if open_handle is not None:
        handles.append(open_handle)
    ax.legend(handles=handles, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0), framealpha=0.96, borderaxespad=0.0)

    footnote = "\n".join(
        textwrap.wrap(
            "analysis_tau.csv, stratum_class == primary, bootstrap CI as computed upstream. " + RADIUS_NOTE,
            width=150,
        )
    )
    fig.text(0.01, -0.04, footnote, ha="left", va="top", fontsize=7, color=C_GREY)
    _take(
        fig,
        "three predictors' rank correlation with the conservative endpoint, per primary stratum; "
        "open markers are LOO-network-unstable, text cells are predictor-constant, not tau=0.",
        y=-0.10,
    )
    return _save(fig, "s9_f1_tau_forest")


# ---------------------------------------------------------------------------
# F2 — survival by radius bucket, shard cells joined to analysis_units.csv.
# ---------------------------------------------------------------------------


def _units_index() -> dict[tuple[str, str], dict]:
    """``analysis_units.csv`` keyed by ``(shard_id, frame_row_id)``, ``status == "ok"`` only."""
    idx: dict[tuple[str, str], dict] = {}
    for r in _read_csv("analysis_units.csv"):
        if r["status"] != "ok":
            continue
        idx[(r["shard_id"], r["frame_row_id"])] = r
    return idx


def _f2_curves(arm_prefix: str, units: dict[tuple[str, str], dict]) -> tuple[
    dict[str, dict[float, tuple[float, int]]], int
]:
    """Mean ``S`` by (radius bucket, grid point) for one arm, ``n_eval >= 30`` cells only.

    Args:
        arm_prefix: ``"flip"`` or ``"tiered"`` — selects shard cell files by
            filename prefix.
        units: Output of :func:`_units_index`, providing each unit's ``radius``.

    Returns:
        ``({bucket: {grid_point: (mean_S, n_contributing_units)}}, n_unreached)``
        — the second element is the count of unit-cells excluded because their
        joined unit's radius is the ``UNREACHED`` sentinel (``<= 0``); it is
        reported in the figure's footnote rather than silently dropped.
    """
    by_bucket: dict[str, dict[float, list[float]]] = defaultdict(lambda: defaultdict(list))
    n_unreached = 0
    cell_files = sorted(glob.glob(str(SHARDS / f"{arm_prefix}__*.cells.jsonl")))
    for fp in cell_files:
        with open(fp) as fh:
            for line in fh:
                d = json.loads(line)
                if d.get("status") != "ok":
                    continue  # unresolved_all_contradictory: S undefined, never imputed.
                n_eval = d.get("n_eval", 0)
                if n_eval is None or n_eval < N_EVAL_MIN:
                    continue
                unit = units.get((d["shard_id"], d["frame_row_id"]))
                if unit is None:
                    continue
                radius_raw = unit.get("radius", "")
                if radius_raw in ("", None):
                    continue
                radius_val = float(radius_raw)
                if radius_val <= 0:
                    n_unreached += 1
                    continue
                bucket = _bucket5(radius_val)
                by_bucket[bucket][float(d["grid_point"])].append(float(d["S"]))
    out = {
        b: {gp: (st.mean(vs), len(vs)) for gp, vs in gmap.items()} for b, gmap in by_bucket.items()
    }
    return out, n_unreached


def fig_survival_by_radius_bucket() -> list[Path]:
    """F2 — mean survival S against the grid point, one line per radius bucket, two panels.

    Left panel: flip arm, x = claim-reversal depth ``d`` (``grid_kind ==
    "d_claims"``). Right panel: tiered arm, x = ``corruption_rate``
    (``grid_kind == "corruption_rate"``). The two panels have independent axes
    (different quantities in different units) and are never merged. Built by
    joining every shard's ``*.cells.jsonl`` to ``analysis_units.csv`` on
    ``(shard_id, frame_row_id)`` — pooled across every stratum within an arm,
    not restricted to the primary nine, matching the pooled-by-bucket
    convention used for the analogous session 7 figure.

    Any unit-cell whose own ``n_eval < 30`` is dropped before averaging;
    marker area encodes the number of contributing units at that (bucket,
    grid point) cell. A radius bucket's line breaks wherever the next grid
    point on that arm's native step has no surviving cell — never bridged.
    """
    _style()
    units = _units_index()
    flip_curves, flip_unreached = _f2_curves("flip", units)
    tiered_curves, tiered_unreached = _f2_curves("tiered", units)

    fig, (axl, axr) = plt.subplots(1, 2, figsize=(13.5, 5.8))

    def _panel(ax: plt.Axes, curves: dict, step: float, xlabel: str) -> None:
        for b in _BUCKET_ORDER:
            if b not in curves or not curves[b]:
                continue
            pts = {gp: sv[0] for gp, sv in curves[b].items()}
            ns = {gp: sv[1] for gp, sv in curves[b].items()}
            x_runs, y_runs = _gapped_xy(pts, step=step)
            sizes_all = _marker_sizes([ns[gp] for run in x_runs for gp in run])
            it = iter(sizes_all)
            for xs, ys in zip(x_runs, y_runs):
                sizes = [next(it) for _ in xs]
                ax.plot(xs, ys, "-", color=_BUCKET_COLOUR[b], lw=1.6, alpha=0.9, zorder=2)
                ax.scatter(
                    xs, ys, s=sizes, color=_BUCKET_COLOUR[b], marker=_BUCKET_MARKER[b],
                    edgecolor="white", linewidth=0.5, zorder=3,
                )
            n_units_total = sum(ns.values())
            ax.plot(
                [], [], "-", marker=_BUCKET_MARKER[b], color=_BUCKET_COLOUR[b],
                label=f"radius={b}  (n_cells={n_units_total})",
            )
        ax.set_xlabel(xlabel)
        ax.set_ylim(-0.04, 1.04)
        ax.grid(axis="x", visible=False)
        ax.legend(title="radius bucket", fontsize=7.6, title_fontsize=8, loc="best", framealpha=0.95)

    _panel(axl, flip_curves, step=1.0, xlabel="corruption depth  d  (claims reversed, flip arm)")
    axl.set_ylabel("mean S  (n_eval ≥ 30 cells only)")
    axl.set_title("flip arm", fontsize=11)

    _panel(axr, tiered_curves, step=0.05, xlabel="corruption rate  (fraction of nodes relocated, tiered arm)")
    axr.set_title("tiered arm", fontsize=11)

    fig.suptitle("Survival by radius bucket, each arm on its own native axis", fontsize=12.5, y=1.02)

    unreached_note = ""
    if flip_unreached or tiered_unreached:
        unreached_note = (
            f" {flip_unreached} flip and {tiered_unreached} tiered unit-cells were excluded for carrying "
            "an UNREACHED radius (never plotted on the numeric axis)."
        )
    footnote = "\n".join(
        textwrap.wrap(
            "Shard cell files joined to analysis_units.csv on (shard_id, frame_row_id); pooled across all "
            "strata within each arm. Marker area ∝ √(contributing units); a line breaks wherever the next "
            "grid point lacks a surviving cell (never bridged)." + unreached_note + " " + RADIUS_NOTE,
            width=170,
        )
    )
    fig.text(0.01, -0.10, footnote, ha="left", va="top", fontsize=7, color=C_GREY)
    _take(
        fig,
        "each arm's own native depth axis, radius-bucket lines pooled across strata; dropped cells break "
        "the line rather than being bridged.",
        y=-0.20,
    )
    return _save(fig, "s9_f2_survival_by_radius_bucket")


# ---------------------------------------------------------------------------
# F3 — per-network view, one chosen primary stratum.
# ---------------------------------------------------------------------------

_F3_STRATUM = "flip_cov100_bw010"


def fig_per_network() -> list[Path]:
    """F3 — per-network median endpoint, one horizontal dot per network, sorted.

    Stratum fixed to ``flip_cov100_bw010`` (the largest of the nine primary
    strata by unit count). Reads ``analysis_per_network.csv`` directly — no
    join, no re-derivation. Marker area encodes ``n_units``; each network's
    ``n_distinct_radius`` is annotated beside its dot. Two reference verticals
    (``median_endpoint_pair_weighted`` and ``median_endpoint_network_weighted``,
    read from ``analysis_tau.csv`` for the same stratum/endpoint — identical
    across all three predictor rows there, so any one row's values are used)
    contrast the corpus-pooled median against the mean of per-network medians,
    which is the whole point of this figure: to make network-to-network spread
    legible against a single average.
    """
    _style()
    all_rows = [r for r in _read_csv("analysis_per_network.csv") if r["stratum"] == _F3_STRATUM]
    if not all_rows:
        raise ValueError(f"no analysis_per_network.csv rows for stratum {_F3_STRATUM!r}")
    # A network can contribute units to the stratum yet have no usable endpoint at all
    # (every one of its units censored) — excluded from the dot plot, not crashed on,
    # and the exclusion is counted below rather than silently dropped.
    rows = [r for r in all_rows if r["median_endpoint_usable"] != ""]
    n_excluded_empty = len(all_rows) - len(rows)
    rows.sort(key=lambda r: float(r["median_endpoint_usable"]))

    tau_rows = [
        r for r in _read_csv("analysis_tau.csv") if r["stratum"] == _F3_STRATUM and r["endpoint"] == "AUC_frac_usable"
    ]
    ref_pair = ref_net = None
    if tau_rows:
        ref_pair = float(tau_rows[0]["median_endpoint_pair_weighted"])
        ref_net = float(tau_rows[0]["median_endpoint_network_weighted"])

    fig, ax = plt.subplots(figsize=(9.5, 8.6))
    ys = list(range(len(rows)))
    xs = [float(r["median_endpoint_usable"]) for r in rows]
    ns = [int(r["n_units"]) for r in rows]
    sizes = _marker_sizes(ns, lo=24, hi=220)

    for y, r in zip(ys, rows):
        if r["q1_endpoint_usable"] == "" or r["q3_endpoint_usable"] == "":
            continue
        q1, q3 = float(r["q1_endpoint_usable"]), float(r["q3_endpoint_usable"])
        ax.plot([q1, q3], [y, y], color=C_GREY, lw=1.0, alpha=0.5, zorder=1)

    ax.scatter(xs, ys, s=sizes, color=C_BLUE, marker="o", edgecolor="white", linewidth=0.6, zorder=3)
    for y, r in zip(ys, rows):
        ax.annotate(
            f"n_radius={r['n_distinct_radius']}",
            xy=(float(r["median_endpoint_usable"]), y),
            xytext=(8, 0),
            textcoords="offset points",
            va="center",
            fontsize=6.6,
            color=C_GREY,
        )

    if ref_pair is not None:
        ax.axvline(ref_pair, color=C_VERM, ls="--", lw=1.4, zorder=2, label="corpus median (pair-weighted)")
    if ref_net is not None:
        ax.axvline(ref_net, color=C_GREEN, ls=":", lw=1.6, zorder=2, label="mean of per-network medians")

    ax.set_yticks(ys)
    ax.set_yticklabels([r["network"] for r in rows], fontsize=7.6)
    ax.set_xlabel(f"median AUC_frac_usable  ({_F3_STRATUM})")
    ax.grid(axis="y", visible=False)
    ax.legend(fontsize=8, loc="lower right", framealpha=0.95)
    ax.set_title(f"Per-network spread, stratum {_F3_STRATUM}", fontsize=11.5)

    excl_note = (
        f" {n_excluded_empty} network(s) contributed units to this stratum but had no usable endpoint "
        "at all (every unit censored) and are excluded from the dot plot." if n_excluded_empty else ""
    )
    footnote = "\n".join(
        textwrap.wrap(
            "analysis_per_network.csv; marker area ∝ √(n_units); grey bars are each network's own "
            "q1-q3 of AUC_frac_usable; verticals from analysis_tau.csv for the same stratum/endpoint."
            + excl_note,
            width=150,
        )
    )
    fig.text(0.01, -0.05, footnote, ha="left", va="top", fontsize=7, color=C_GREY)
    _take(
        fig,
        "networks sorted by their own median; dot area is unit count, grey bars are within-network "
        "IQR, verticals contrast the corpus-pooled median against the average of per-network medians.",
        y=-0.09,
    )
    return _save(fig, "s9_f3_per_network")


# ---------------------------------------------------------------------------
# F4 — leverage vs endpoint, one point per network.
# ---------------------------------------------------------------------------

_F4_STRATUM = "flip_cov100_bw010"
_F4_HIGHLIGHT = {"paths", "hailfinder"}


def fig_leverage() -> list[Path]:
    """F4 — per-network leverage ``k_g0/n_k`` (x, log) vs median endpoint (y).

    Restricted to ``analysis_units.csv`` rows with ``stratum ==
    "flip_cov100_bw010"`` and ``status == "ok"``. Per network: leverage is the
    median of ``k_g0/n_k`` over its units (units with an empty ``k_g0`` are
    excluded from that median, never treated as leverage 0); the endpoint is
    the median ``AUC_frac_usable``. Marker area encodes the network's unit
    count. A vertical rule at leverage = 1 is annotated "K is Meek-closed"
    (leverage > 1 means the committed adjustment set draws on more of K than
    what Meek closure alone would certify). ``paths`` and ``hailfinder`` are
    highlighted with a distinct marker, per the session brief.
    """
    _style()
    rows = [r for r in _read_csv("analysis_units.csv") if r["stratum"] == _F4_STRATUM and r["status"] == "ok"]
    if not rows:
        raise ValueError(f"no ok analysis_units.csv rows for stratum {_F4_STRATUM!r}")

    by_net: dict[str, dict] = defaultdict(lambda: {"lev": [], "ep": [], "n": 0})
    for r in rows:
        by_net[r["network"]]["n"] += 1
        if r.get("k_g0", "") not in ("", None):
            n_k = float(r["n_k"])
            if n_k > 0:
                by_net[r["network"]]["lev"].append(float(r["k_g0"]) / n_k)
        if r.get("AUC_frac_usable", "") not in ("", None):
            by_net[r["network"]]["ep"].append(float(r["AUC_frac_usable"]))

    nets = sorted(n for n, d in by_net.items() if d["lev"] and d["ep"])
    lev = [st.median(by_net[n]["lev"]) for n in nets]
    ep = [st.median(by_net[n]["ep"]) for n in nets]
    ns = [by_net[n]["n"] for n in nets]
    sizes = _marker_sizes(ns, lo=30, hi=260)

    fig, ax = plt.subplots(figsize=(9.5, 7.2))
    highlight_mask = [n in _F4_HIGHLIGHT for n in nets]
    for x, y, s, n, hl in zip(lev, ep, sizes, nets, highlight_mask):
        if hl:
            ax.scatter([x], [y], s=s * 1.4, marker="*", color=C_VERM, edgecolor="black", linewidth=0.8, zorder=4)
        else:
            ax.scatter([x], [y], s=s, marker="o", color=C_BLUE, edgecolor="white", linewidth=0.5, zorder=3)
        ax.annotate(n, xy=(x, y), xytext=(6, 4), textcoords="offset points", fontsize=6.6, color=C_GREY)

    ax.set_xscale("log")
    ax.axvline(1.0, color=C_GREY, ls="--", lw=1.2, zorder=1)
    ax.text(1.0, ax.get_ylim()[1] if ep else 1.0, " K is Meek-closed", fontsize=8, color=C_GREY, va="top", ha="left")
    ax.set_xlabel("leverage  =  k_g0 / n_k  (median per network, log scale)")
    ax.set_ylabel(f"median AUC_frac_usable  ({_F4_STRATUM})")
    ax.set_title("Leverage vs endpoint, one point per network", fontsize=11.5)

    handles = [
        Line2D([], [], ls="", marker="o", color=C_BLUE, markersize=8, label="network"),
        Line2D([], [], ls="", marker="*", color=C_VERM, markersize=11, label="paths / hailfinder"),
    ]
    ax.legend(handles=handles, fontsize=8, loc="best", framealpha=0.95)

    footnote = "\n".join(
        textwrap.wrap(
            "analysis_units.csv, stratum flip_cov100_bw010, status == ok only; per-network leverage and "
            "endpoint are medians over that network's units; marker area ∝ √(n_units).",
            width=150,
        )
    )
    fig.text(0.01, -0.05, footnote, ha="left", va="top", fontsize=7, color=C_GREY)
    _take(
        fig,
        "per-network leverage k_g0/n_k against median endpoint, log-x; the leverage=1 line separates "
        "Meek-closed K from K that draws on more of the background knowledge.",
        y=-0.09,
    )
    return _save(fig, "s9_f4_leverage")


# ---------------------------------------------------------------------------
# F5 — knowledge-model paired dot plot.
# ---------------------------------------------------------------------------


def fig_knowledge_model() -> list[Path]:
    """F5 — exhaustive single-reversal contradiction rate, two knowledge models, coverage=1.0.

    Reads ``knowledge_model_diagnostic.csv`` filtered to ``coverage == "1.0"``,
    one row per (network, knowledge_model). Paired per network: the
    ``minimal_generator`` and ``all_undirected_edges`` contradiction rates are
    joined by a line, networks sorted by the difference. Each model's pooled
    rate (``knowledge_model_diagnostic.json``, keys
    ``minimal_generator_coverage_1.0`` and ``all_undirected_edges`` — the
    latter has no separate per-coverage key because it was only run at
    coverage 1.0) is drawn as a dashed vertical, read from the JSON rather
    than recomputed here.
    """
    _style()
    rows = [r for r in _read_csv("knowledge_model_diagnostic.csv") if r["coverage"] == "1.0"]
    if not rows:
        raise ValueError("no knowledge_model_diagnostic.csv rows at coverage == 1.0")
    with (RES / "knowledge_model_diagnostic.json").open() as fh:
        diag = json.load(fh)
    mg_pooled = diag.get("minimal_generator_coverage_1.0", {}).get("pooled_contradiction_rate")
    au_pooled = diag.get("all_undirected_edges", {}).get("pooled_contradiction_rate")

    mg = {r["network"]: float(r["contradiction_rate"]) for r in rows if r["knowledge_model"] == "minimal_generator"}
    au = {r["network"]: float(r["contradiction_rate"]) for r in rows if r["knowledge_model"] == "all_undirected_edges"}
    nets = sorted(set(mg) & set(au), key=lambda n: au[n] - mg[n])

    fig, ax = plt.subplots(figsize=(8.6, 8.2))
    ys = list(range(len(nets)))
    for y, n in zip(ys, nets):
        ax.plot([mg[n], au[n]], [y, y], color=C_GREY, lw=1.1, zorder=1)
    ax.scatter([mg[n] for n in nets], ys, s=46, color=C_BLUE, marker="o", edgecolor="white", linewidth=0.5, zorder=3, label="minimal_generator")
    ax.scatter([au[n] for n in nets], ys, s=46, color=C_VERM, marker="^", edgecolor="white", linewidth=0.5, zorder=3, label="all_undirected_edges")

    if mg_pooled is not None:
        ax.axvline(mg_pooled, color=C_BLUE, ls="--", lw=1.2, zorder=2, label="minimal_generator pooled rate")
    if au_pooled is not None:
        ax.axvline(au_pooled, color=C_VERM, ls=":", lw=1.4, zorder=2, label="all_undirected_edges pooled rate")

    ax.set_yticks(ys)
    ax.set_yticklabels(nets, fontsize=7.6)
    ax.set_xlabel("exhaustive single-reversal contradiction rate (coverage = 1.0)")
    ax.grid(axis="y", visible=False)
    ax.legend(fontsize=7.6, loc="lower right", framealpha=0.95)
    ax.set_title("Contradiction rate under two knowledge models, paired per network", fontsize=11)

    footnote = "\n".join(
        textwrap.wrap(
            "knowledge_model_diagnostic.csv, coverage == 1.0; networks sorted by "
            "(all_undirected_edges − minimal_generator); pooled rates from knowledge_model_diagnostic.json.",
            width=150,
        )
    )
    fig.text(0.01, -0.05, footnote, ha="left", va="top", fontsize=7, color=C_GREY)
    _take(
        fig,
        "paired per-network contradiction rate under two knowledge models, sorted by within-network "
        "difference, with each model's pooled rate marked.",
        y=-0.09,
    )
    return _save(fig, "s9_f5_knowledge_model")


# ---------------------------------------------------------------------------
# F6 — cross-arm intensity-matched comparison.
# ---------------------------------------------------------------------------

_F6_METRIC = "silent_failure_rate"


def fig_crossarm() -> list[Path]:
    """F6 — achieved-intensity histogram (left) and paired mean difference (right).

    Left: ``xarm_intensity_histogram.csv``, share of draws per intensity bin,
    one series per arm (``xarm_flip`` / ``xarm_tiered``), normalized within
    arm so the two very different totals are comparable as *shapes*, not
    absolute counts. Right: ``xarm_paired_test.csv``, metric
    ``silent_failure_rate``, the paired mean difference (tiered − flip) per
    intensity bin with its bootstrap CI computed **over networks**
    (``silent_failure_rate_boot_ci_lo_2p5`` / ``_hi_97p5``). Bins with fewer
    than 15 matched instances (``n_instances_matched``) are drawn with an open
    marker and a light shaded band, marking them as not carrying the
    comparison rather than dropping them outright.

    This is the one figure in this session where flip and tiered legitimately
    share an axis — the whole point of the xarm design is to place both arms
    on a common achieved-intensity axis, unlike the native d / corruption_rate
    axes the rest of this module keeps strictly separate.
    """
    _style()
    hist_rows = _read_csv("xarm_intensity_histogram.csv")
    paired_rows = _read_csv("xarm_paired_test.csv")
    if not hist_rows or not paired_rows:
        raise ValueError("xarm_intensity_histogram.csv or xarm_paired_test.csv is empty")

    by_arm: dict[str, dict[float, int]] = defaultdict(dict)
    for r in hist_rows:
        by_arm[r["arm"]][float(r["intensity_bin"])] = int(r["n_draws"])

    fig, (axl, axr) = plt.subplots(1, 2, figsize=(13.5, 5.6))

    arm_colour = {"xarm_flip": C_BLUE, "xarm_tiered": C_VERM}
    arm_hatch = {"xarm_flip": "//", "xarm_tiered": "\\\\"}
    width = 0.02
    for i, (arm, colour) in enumerate(arm_colour.items()):
        bins = sorted(by_arm[arm])
        total = sum(by_arm[arm].values())
        shares = [by_arm[arm][b] / total for b in bins]
        offset = (i - 0.5) * width
        axl.bar(
            [b + offset for b in bins], shares, width=width, color=colour, edgecolor="black",
            linewidth=0.3, hatch=arm_hatch[arm], alpha=0.8, label=arm,
        )
    axl.set_xlabel("achieved intensity bin")
    axl.set_ylabel("share of draws within arm")
    axl.set_title("Where each arm's mass actually lands", fontsize=10.5)
    axl.legend(fontsize=8, loc="upper right", framealpha=0.95)

    bins = [float(r["intensity_bin"]) for r in paired_rows]
    diffs = [float(r[f"{_F6_METRIC}_mean_diff_tiered_minus_flip"]) for r in paired_rows]
    lo = [float(r[f"{_F6_METRIC}_boot_ci_lo_2p5"]) for r in paired_rows]
    hi = [float(r[f"{_F6_METRIC}_boot_ci_hi_97p5"]) for r in paired_rows]
    n_matched = [int(r["n_instances_matched"]) for r in paired_rows]

    for b, d, l, h, n in zip(bins, diffs, lo, hi, n_matched):
        carries = n >= N_MATCHED_MIN
        colour = C_BLUE if carries else C_GREY
        marker = "o" if carries else "x"
        if not carries:
            axr.axvspan(b - 0.02, b + 0.02, color=C_GREY, alpha=0.12, hatch="..", zorder=0)
        axr.plot([b, b], [l, h], color=colour, lw=1.3, alpha=0.85, zorder=2)
        axr.scatter([b], [d], s=42, color=colour, marker=marker, zorder=3)

    axr.axhline(0.0, color=C_GREY, ls="--", lw=1.0, zorder=1)
    axr.set_xlabel("achieved intensity bin")
    axr.set_ylabel(f"{_F6_METRIC}: tiered − flip  (bootstrap CI over networks)")
    axr.set_title("Paired difference per intensity bin", fontsize=10.5)
    handles = [
        Line2D([], [], ls="", marker="o", color=C_BLUE, markersize=7, label=f"n_matched ≥ {N_MATCHED_MIN}"),
        Line2D([], [], ls="", marker="x", color=C_GREY, markersize=7, label=f"n_matched < {N_MATCHED_MIN} (shaded, not carrying comparison)"),
    ]
    axr.legend(handles=handles, fontsize=7.4, loc="best", framealpha=0.95)

    fig.suptitle("Cross-arm intensity-matched comparison (achieved-intensity axis, both arms)", fontsize=12, y=1.03)
    footnote = "\n".join(
        textwrap.wrap(
            "Left: xarm_intensity_histogram.csv, share of draws normalized within arm. Right: "
            "xarm_paired_test.csv, metric silent_failure_rate, bootstrap CI over networks; bins under "
            f"{N_MATCHED_MIN} matched instances are marked, not dropped.",
            width=170,
        )
    )
    fig.text(0.01, -0.10, footnote, ha="left", va="top", fontsize=7, color=C_GREY)
    _take(
        fig,
        "left shows each arm's own achieved-intensity distribution; right shows the paired tiered-minus-"
        "flip difference with across-network bootstrap CIs, low-n bins marked rather than hidden.",
        y=-0.20,
    )
    return _save(fig, "s9_f6_crossarm")


def build_all() -> list[Path]:
    """Build every session 9 figure that has data.

    Returns:
        The paths written. A figure whose required input file or column is
        missing, or whose shard join comes up empty, is skipped with a
        printed message rather than aborting the run.
    """
    out: list[Path] = []
    for fn in (
        fig_tau_forest,
        fig_survival_by_radius_bucket,
        fig_per_network,
        fig_leverage,
        fig_knowledge_model,
        fig_crossarm,
    ):
        try:
            out.extend(fn())
        except (FileNotFoundError, ValueError, IndexError, KeyError) as exc:
            print(f"skipped {fn.__name__}: {exc}")
    return out


def main() -> None:
    for p in build_all():
        print(p)


if __name__ == "__main__":
    main()
