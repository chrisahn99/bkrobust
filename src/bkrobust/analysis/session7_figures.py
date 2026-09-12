"""Figures for session 7 (Axis Robustness: survival curves and the Pareto premise).

Every figure is built from committed files under ``results/axis_robustness/`` and
must agree with ``report_fragility_and_pareto.md`` (root of the repo) and the
corrections in ``results/axis_robustness/PREREGISTRATION.md`` Appendices E and F.
Inputs are final; nothing here regenerates or re-derives a sweep.

Two arms, never merged onto one axis. **flip** corrupts by reversing a targeted
count of asserted claims, so its native knob is an integer depth ``d`` (fraction
of claims reversed = ``d / n_k``). **tiered** corrupts by relocating nodes
between tiers, so its native knob is ``corruption_rate`` (fraction of nodes
relocated) — a *different quantity in different units*, not a depth. The
tiered arm's own ``d`` column in ``survival_curves.csv`` is defective (its
``d = 0`` bin is contaminated by every corruption rate; see PREREGISTRATION
Appendix E) and is never plotted here; every tiered figure reads
``analysis_tiered_rebinned.csv`` instead.

A data-fidelity note on F5, recorded here because it changed what the figure
shows relative to a naive reading of the report. Section 2.7 of
``report_fragility_and_pareto.md`` captions a contradiction-rate-by-depth table
"flip arm" and quotes 0.491 / 0.773 / 0.896 / ... / 1.000. Recomputing that
table from ``survival_curves.csv`` restricted to ``arm == "flip"`` reproduces
none of those numbers (flip-only gives 0.627 / 0.777 / 0.876 / ... at d = 1/2/3);
recomputing it pooled over *both* arms — using the tiered arm's own defective
``d`` column — reproduces the quoted table exactly, to the third decimal. The
label in the report is wrong, not the arithmetic: the table is a both-arms pool
mislabeled as flip-only, and pooling the two arms on one ``d`` axis is exactly
what this file's house rule forbids. :func:`fig_self_revealing` therefore plots
the true flip-only curve — correct, defensible, and traceable to
``survival_curves.csv`` — rather than reproducing the mislabeled figure. The
one part of that report section that *is* flip-only, the base-wrongness
breakdown at ``d = 1`` (0.653 / 0.647 / 0.563), is reproduced exactly and is
used unchanged.

Effective sample size behind ``S(d)`` varies enormously with depth (contradiction
rates run 0.49 at d=1 to 1.00 at d=12), so every *per-depth* curve here (F2, F6)
drops any (bucket, depth) cell whose contributing rows have ``n_eval < 30`` and
further encodes the number of contributing instances as marker area, per the
pre-registered ``AUC_frac_usable`` convention (PREREGISTRATION Appendix F).
``S`` undefined (all-contradictory) is never imputed or interpolated across:
depths are only connected when they are numerically consecutive, so a dropped
depth breaks the line rather than being bridged. F1 sidesteps this floor
entirely by plotting the whole-curve endpoint ``AUC_frac`` per instance
(pre-registered, computed once per instance, no per-depth cell to threshold),
which is also why it is the primary figure and F6 (the per-depth view) is
supplementary.

Colour-blind safe Okabe-Ito palette throughout, shared with sessions 4-6;
every categorical encoding is doubled with a distinct marker shape so nothing
depends on colour alone in greyscale or print.
"""

# ruff: noqa: RUF001
# Figure text is display copy; the minus signs and en dashes are intentional.

from __future__ import annotations

import csv
import statistics as st
import textwrap
from collections import Counter, defaultdict
from pathlib import Path
from typing import Callable

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

RES = Path("results/axis_robustness")
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

N_EVAL_MIN = 30  # never plot a point backed by fewer draws than this.
RADIUS_NOTE = (
    "r_val carries Conjecture 2 (hence Anti-Exchange Case B, verified not proved); "
    "the error is one-sided, so radii can only be too large."
)


def _style() -> None:
    """Apply the shared house rcParams (see ``session4_figures._style``)."""
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
    """Read a CSV under ``results/axis_robustness`` into a list of dict rows."""
    with (RES / name).open(newline="") as fh:
        return list(csv.DictReader(fh))


def _instances() -> dict[str, dict]:
    """``survival_instances.csv`` keyed by ``instance_id`` (1,978 rows)."""
    return {r["instance_id"]: r for r in _read_csv("survival_instances.csv")}


def _curves() -> list[dict]:
    """``survival_curves.csv`` (14,066 rows), both arms, unfiltered."""
    return _read_csv("survival_curves.csv")


def _bucket5(value: str) -> str:
    """Bucket a non-negative integer-valued predictor into ``0..4, "5+"``.

    Used for both ``r_val`` (range 1-11 in this dataset) and ``shd_truth``
    (range 0-14), matching the bucket edges already committed in
    ``analysis_effect_sizes.csv``.
    """
    v = int(float(value))
    return str(v) if v < 5 else "5+"


_BUCKET_ORDER = ["0", "1", "2", "3", "4", "5+"]


def _bucket_sort_key(b: str) -> tuple[bool, int]:
    return (b == "5+", 0 if b == "5+" else int(b))


def _gapped_xy(
    points: dict[float, float], step: float, tol: float = 1e-9
) -> tuple[list[list[float]], list[list[float]]]:
    """Split ``{x: y}`` into runs of consecutive ``x`` (spacing ``step``).

    A missing x between two present ones must not be bridged by a drawn line —
    that would imply data where there is none (undefined ``S`` is a gap, never
    an interpolation). Returns parallel lists of x/y runs, sorted by x.

    Args:
        points: Mapping from x-coordinate to y-value.
        step: The expected spacing between consecutive x values.
        tol: Floating-point tolerance for "consecutive".

    Returns:
        ``(x_runs, y_runs)`` — same-length lists of runs; each run has >= 1 point.
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


def _marker_sizes(ns: list[int], lo: float = 18.0, hi: float = 130.0) -> list[float]:
    """Map contributing-instance counts to marker areas, sqrt-scaled, for one panel.

    Args:
        ns: Contributing-instance (or contributing-row) counts, one per point.
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


# ---------------------------------------------------------------------------
# F1 — Endpoint (AUC_frac) distribution by bucket, flip arm.
# ---------------------------------------------------------------------------

# Shared bucket colours/markers, fixed across F1 and F6 so the same bucket
# always reads the same way. "0" (shd_truth only) is C_BLACK, matching the
# hollow-X "predictor constant" convention used for it elsewhere (F3).
_BUCKET_COLOUR = {
    "0": C_BLACK,
    "1": C_ORANGE,
    "2": C_SKY,
    "3": C_GREEN,
    "4": C_VERM,
    "5+": C_PINK,
}
_BUCKET_MARKER = {"0": "P", "1": "o", "2": "s", "3": "^", "4": "D", "5+": "v"}


def _flip_endpoint_by_bucket(key: str, bucket_fn: Callable[[str], str]) -> dict[str, list[float]]:
    """``AUC_frac`` per instance, grouped by bucket of ``key``, flip arm only.

    Args:
        key: Column of ``survival_instances.csv`` to bucket on (``r_val`` or
            ``shd_truth``).
        bucket_fn: Maps a raw value (as a string) to a bucket label.

    Returns:
        ``{bucket: [AUC_frac, ...]}`` — one value per flip-arm instance in
        that bucket, unfiltered (``AUC_frac`` is defined for every instance;
        this is the pre-registered whole-curve endpoint, not a per-depth cell,
        so there is no ``n_eval`` floor to apply here).
    """
    out: dict[str, list[float]] = defaultdict(list)
    for rec in _instances().values():
        if rec["arm"] != "flip":
            continue
        out[bucket_fn(rec[key])].append(float(rec["AUC_frac"]))
    return dict(out)


def _box_strip_panel(
    ax: plt.Axes, data: dict[str, list[float]], buckets: list[str], rng: np.random.Generator
) -> None:
    """Draw one box+jittered-strip panel: one box per bucket, in ``buckets`` order.

    Args:
        ax: Axes to draw into.
        data: ``{bucket: [values]}`` from :func:`_flip_endpoint_by_bucket`.
        buckets: Bucket labels in display order; a bucket absent from ``data``
            is skipped (no box, no tick).
        rng: Seeded generator for the jitter — deterministic, never the global
            RNG, so re-running produces byte-identical output.
    """
    present = [b for b in buckets if b in data and data[b]]
    positions = list(range(1, len(present) + 1))
    for pos, b in zip(positions, present):
        vals = data[b]
        colour = _BUCKET_COLOUR[b]
        bp = ax.boxplot(
            [vals],
            positions=[pos],
            widths=0.5,
            patch_artist=True,
            showfliers=False,
            zorder=3,
        )
        for box in bp["boxes"]:
            box.set(facecolor=colour, alpha=0.22, edgecolor=colour, linewidth=1.4)
        for part in ("whiskers", "caps"):
            for line in bp[part]:
                line.set(color=colour, linewidth=1.2)
        for med in bp["medians"]:
            med.set(color=colour, linewidth=2.2)
        jitter = rng.uniform(-0.16, 0.16, size=len(vals))
        ax.scatter(
            [pos] * len(vals) + jitter,
            vals,
            s=10,
            marker=_BUCKET_MARKER[b],
            color=colour,
            alpha=0.28,
            linewidth=0,
            zorder=2,
        )
    ax.set_xticks(positions)
    ax.set_xticklabels([f"{b}\n(n={len(data[b])})" for b in present], fontsize=8.6)
    ax.set_xlim(0.35, len(present) + 0.65)
    ax.grid(axis="x", visible=False)


def fig_endpoint_by_bucket() -> list[Path]:
    """F1 — ``AUC_frac`` (the pre-registered endpoint) by bucket, flip arm only.

    Two panels sharing a y-axis. Left: ``AUC_frac`` distribution by ``r_val``
    bucket (1, 2, 3, 4, 5+). Right: the identical construction by
    ``shd_truth`` bucket (0, 1, 2, 3, 4, 5+). Box = quartiles + median;
    individual instances overlaid as a jittered strip (``np.random.default_rng(0)``,
    fixed seed, never the global RNG) so the reader sees the n behind each box,
    not only its summary. This is the quantity tau_b in F3 actually ranks —
    unlike the old F1 (now F6), which plotted per-depth means, a different and
    much noisier statistic than the whole-curve endpoint the analysis tests.

    Pooling every stratum (coverage, base_wrongness) onto one bucket axis, as
    done here for readability, is *not* what F3 does — F3 keeps strata
    separate and r_val's tau_b is positive in every one of them (0.03 to
    0.71). The pooled view here is expected to show the same sign but a
    weaker, less clean version of that separation, because mixing strata of
    different difficulty adds variance the stratified test does not have.
    """
    _style()
    r_data = _flip_endpoint_by_bucket("r_val", _bucket5)
    shd_data = _flip_endpoint_by_bucket("shd_truth", _bucket5)
    rng = np.random.default_rng(0)

    fig, (axl, axr) = plt.subplots(1, 2, figsize=(13.0, 5.8), sharey=True)

    _box_strip_panel(axl, r_data, ["1", "2", "3", "4", "5+"], rng)
    axl.set_xlabel("r_val bucket")
    axl.set_ylabel("AUC_frac  (flip arm, one point per instance)")
    axl.set_title("By r_val bucket: rises overall, not monotonically", fontsize=10.5)

    _box_strip_panel(axr, shd_data, _BUCKET_ORDER, rng)
    axr.set_xlabel("shd_truth bucket")
    axr.set_title("By shd_truth bucket: no visible ordering", fontsize=10.5)

    ax_ylim = (-0.06, 1.06)
    axl.set_ylim(*ax_ylim)

    fig.suptitle(
        "AUC_frac (flip arm), the pre-registered whole-curve endpoint, by bucket",
        fontsize=12.5,
        y=1.02,
    )
    footnote_text = (
        "Flip arm; AUC_frac per instance; boxes = quartiles/median, n per tick; points jittered "
        "(seed 0). shd_truth = 0 also includes every coverage=1.0/base-wrongness=0 case, where it "
        "is ≡ 0 by construction — a separate reason it cannot rank, not merely tangling. "
        + RADIUS_NOTE
        + " F3's per-stratum tau_b remains primary evidence for P1/P2."
    )
    footnote = "\n".join(textwrap.wrap(footnote_text, width=160))
    fig.text(0.5, -0.13, footnote, ha="center", va="top", fontsize=7.4, color=C_GREY)
    return _save(fig, "session7_f1_stratification")


# ---------------------------------------------------------------------------
# F6 — Supplementary: per-depth survival means by bucket (noise-dominated tails).
# ---------------------------------------------------------------------------


def _flip_bucket_curves(
    key: str, bucket_fn: Callable[[str], str]
) -> dict[str, dict[int, tuple[float, int]]]:
    """Mean ``S(d)`` by bucket of ``key``, flip arm only, ``n_eval >= 30`` cells only.

    Args:
        key: Column of ``survival_instances.csv`` to bucket on (``r_val`` or
            ``shd_truth``).
        bucket_fn: Maps a raw value (as a string) to a bucket label.

    Returns:
        ``{bucket: {depth: (mean_S, n_contributing_instances)}}``. A depth
        missing from a bucket's dict has no cell with ``n_eval >= 30`` in that
        bucket and must not be plotted, imputed, or bridged.
    """
    inst = _instances()
    by_bd: dict[str, dict[int, list[float]]] = defaultdict(lambda: defaultdict(list))
    for row in _curves():
        if row["arm"] != "flip":
            continue
        rec = inst[row["instance_id"]]
        s_raw, n_eval = row["S"], int(row["n_eval"])
        if s_raw == "" or n_eval < N_EVAL_MIN:
            continue
        by_bd[bucket_fn(rec[key])][int(row["d"])].append(float(s_raw))
    return {
        b: {d: (st.mean(vs), len(vs)) for d, vs in dmap.items()} for b, dmap in by_bd.items()
    }


def _bucket_n_instances(key: str, bucket_fn: Callable[[str], str], arm: str) -> dict[str, int]:
    """Total distinct instances per bucket of ``key``, within one arm (for legends)."""
    counts: Counter = Counter()
    for rec in _instances().values():
        if rec["arm"] != arm:
            continue
        counts[bucket_fn(rec[key])] += 1
    return dict(counts)



def fig_survival_curves_by_bucket() -> list[Path]:
    """F6 (supplementary) — per-depth mean ``S(d)`` by bucket; tails are noise-dominated.

    Two panels, same y-axis, flip arm only (corruption depth ``d`` = number of
    reversed claims; ``d / n_k`` is the fraction). Left: mean ``S(d)`` per
    ``r_val`` bucket (1, 2, 3, 4, 5+). Right: the identical construction with
    ``shd_truth`` buckets (0, 1, 2, 3, 4, 5+) in place of ``r_val``.

    This is **not** the pre-registered endpoint — that is ``AUC_frac``, a
    whole-curve summary, plotted by bucket in F1. Per-depth means are a
    different, much noisier statistic: effective sample size collapses fast
    with depth (contradiction rate 0.63 at d=1 rising to 1.00 by d=12, see the
    module docstring), so both panels here show real bucket crossings once n
    thins past d ≈ 3-4 and neither should be read as showing clean
    separation. It is kept as a supplementary, per-depth view of the same
    underlying curves; F1 is the endpoint-level summary and F3 is the
    statistical test.

    ``shd_truth`` also carries a second, independent reason it cannot rank at
    coverage 1.0 / base wrongness 0: it is *undefined there by construction*
    (``G0`` **is** the truth, so ``shd_truth ≡ 0`` — see
    ``analysis_tau_primary.csv``, ``predictor_constant = True``). Those
    instances are folded into the ``shd_truth = 0`` bucket here (a real value
    for the rest of the design), which is a second, separate failure mode from
    the tangling itself.
    """
    _style()
    r_curves = _flip_bucket_curves("r_val", _bucket5)
    shd_curves = _flip_bucket_curves("shd_truth", _bucket5)
    r_n = _bucket_n_instances("r_val", _bucket5, "flip")
    shd_n = _bucket_n_instances("shd_truth", _bucket5, "flip")

    colours = _BUCKET_COLOUR
    markers = _BUCKET_MARKER

    fig, (axl, axr) = plt.subplots(1, 2, figsize=(12.4, 5.4), sharey=True)

    def _panel(ax: plt.Axes, curves: dict, n_by_b: dict, buckets: list[str], title: str) -> None:
        all_n = [n for bmap in curves.values() for _, n in bmap.values()]
        for b in buckets:
            if b not in curves:
                continue
            pts = {d: sv[0] for d, sv in curves[b].items()}
            ns = {d: sv[1] for d, sv in curves[b].items()}
            x_runs, y_runs = _gapped_xy(pts, step=1)
            sizes_all = _marker_sizes([ns[d] for run in x_runs for d in run], hi=130.0)
            it = iter(sizes_all)
            for xs, ys in zip(x_runs, y_runs):
                sizes = [next(it) for _ in xs]
                ax.plot(xs, ys, "-", color=colours[b], lw=1.6, zorder=2, alpha=0.9)
                ax.scatter(
                    xs,
                    ys,
                    s=sizes,
                    color=colours[b],
                    marker=markers[b],
                    edgecolor="white",
                    linewidth=0.5,
                    zorder=3,
                )
            ax.plot(
                [],
                [],
                "-",
                marker=markers[b],
                color=colours[b],
                label=f"{title[:1].lower()}={b}  (n={n_by_b.get(b, 0)})",
            )
        ax.set_xlabel("corruption depth  d  (claims reversed, flip arm)")
        ax.set_ylim(-0.04, 1.04)
        ax.set_xlim(0.5, 11.5)
        ax.grid(axis="x", visible=False)

    _panel(axl, r_curves, r_n, ["1", "2", "3", "4", "5+"], "r_val")
    axl.set_ylabel("mean S(d)  (flip arm, cells with n_eval ≥ 30 only)")
    axl.set_title("r_val bucket (per-depth mean, supplementary)", fontsize=10.5)
    axl.legend(title="r_val bucket", fontsize=7.6, title_fontsize=8, loc="lower left", framealpha=0.95)

    _panel(axr, shd_curves, shd_n, _BUCKET_ORDER, "shd_truth")
    axr.set_title("shd_truth bucket (per-depth mean, supplementary)", fontsize=10.5)
    axr.legend(
        title="shd_truth bucket", fontsize=7.6, title_fontsize=8, loc="lower left", framealpha=0.95
    )

    fig.suptitle(
        "Supplementary: mean S(d) per depth, by bucket (flip arm) — noise-dominated at high d",
        fontsize=12,
        y=1.05,
    )
    footnote_text = (
        "Supplementary per-depth view: mean S(d) by bucket, flip arm, n_eval ≥ 30 cells only; marker "
        "area ∝ √(contributing instances); a line breaks wherever the next depth lacks such a cell "
        "(never bridged). n thins fast with depth, so this view is noise-dominated at high d — both "
        "panels show real crossings among the middle buckets past d ≈ 3-4. shd_truth = 0 also "
        "includes every coverage=1.0/base-wrongness=0 instance, where it is ≡ 0 by construction — a "
        "separate reason it cannot rank, not merely tangling. " + RADIUS_NOTE
        + " See F1 for the endpoint (AUC_frac) and F3 for the rank-correlation test — the primary "
        "evidence for P1/P2."
    )
    footnote = "\n".join(textwrap.wrap(footnote_text, width=150))
    fig.text(0.5, -0.14, footnote, ha="center", va="top", fontsize=7, color=C_GREY)
    return _save(fig, "session7_f6_survival_curves_by_bucket")


# ---------------------------------------------------------------------------
# F2 — Tiered stratification, re-binned axis.
# ---------------------------------------------------------------------------

_RATES = ["0.00", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40", "0.45", "0.50"]


def fig_tiered_stratification() -> list[Path]:
    """F2 — tiered survival by ``r_val`` bucket, on the re-binned ``corruption_rate`` axis.

    Reads ``analysis_tiered_rebinned.csv`` exclusively (never the tiered arm's
    ``d`` column in ``survival_curves.csv``, which Appendix E documents as
    contaminated at ``d = 0``). Each point is the mean of ``S_rate_<rate>``
    across instances in the bucket with ``n_eval_rate_<rate> >= 30``, so this
    figure has the same evidentiary floor as F1 despite reading a different
    file. ``S = 1.000`` at ``rate = 0`` for every one of the 668 instances
    individually, so every bucket's line starts exactly at 1.0 by construction,
    not by fit.
    """
    _style()
    inst = _instances()
    rows = _read_csv("analysis_tiered_rebinned.csv")

    by_b: dict[str, dict[float, list[float]]] = defaultdict(lambda: defaultdict(list))
    n_by_b: dict[str, set] = defaultdict(set)
    for row in rows:
        rec = inst.get(row["instance_id"])
        if rec is None:
            continue
        b = _bucket5(rec["r_val"])
        n_by_b[b].add(row["instance_id"])
        for rate in _RATES:
            s_raw = row.get(f"S_rate_{rate}", "")
            ne_raw = row.get(f"n_eval_rate_{rate}", "")
            if s_raw == "" or ne_raw == "" or int(ne_raw) < N_EVAL_MIN:
                continue
            by_b[b][float(rate)].append(float(s_raw))

    colours = {"1": C_ORANGE, "2": C_SKY, "3": C_GREEN, "4": C_VERM, "5+": C_PINK}
    markers = {"1": "o", "2": "s", "3": "^", "4": "D", "5+": "v"}

    fig, ax = plt.subplots(figsize=(8.0, 5.6))
    for b in ["1", "2", "3", "4", "5+"]:
        if b not in by_b:
            continue
        pts = {r: st.mean(vs) for r, vs in by_b[b].items()}
        ns = {r: len(vs) for r, vs in by_b[b].items()}
        x_runs, y_runs = _gapped_xy(pts, step=0.05)
        sizes_all = _marker_sizes([ns[r] for run in x_runs for r in run], hi=120.0)
        it = iter(sizes_all)
        for xs, ys in zip(x_runs, y_runs):
            sizes = [next(it) for _ in xs]
            ax.plot(xs, ys, "-", color=colours[b], lw=1.8, alpha=0.9, zorder=2)
            ax.scatter(
                xs, ys, s=sizes, color=colours[b], marker=markers[b],
                edgecolor="white", linewidth=0.5, zorder=3,
            )
        ax.plot(
            [], [], "-", marker=markers[b], color=colours[b],
            label=f"r_val={b}  (n={len(n_by_b[b])})",
        )

    ax.set_xlabel("corruption rate  (fraction of nodes relocated, tiered arm)")
    ax.set_ylabel("mean S(corruption_rate)  (cells with n_eval ≥ 30 only)")
    ax.set_xlim(-0.02, 0.52)
    ax.set_ylim(-0.04, 1.04)
    ax.set_xticks([float(r) for r in _RATES])
    ax.grid(axis="x", visible=False)
    ax.legend(title="r_val bucket", fontsize=8, title_fontsize=8.5, loc="lower left", framealpha=0.95)
    ax.set_title(
        "Tiered arm: r_val separates survival on the re-binned axis, same as flip",
        fontsize=11,
    )
    fig.text(
        0.5,
        -0.06,
        "analysis_tiered_rebinned.csv (668 instances); marker area ∝ √(contributing instances).\n"
        "S = 1.000 at rate 0 for all 668 instances individually (PREREGISTRATION Appendix E.2). "
        + RADIUS_NOTE,
        ha="center",
        fontsize=7,
        color=C_GREY,
    )
    return _save(fig, "session7_f2_tiered_stratification")


# ---------------------------------------------------------------------------
# F3 — Tau forest plot.
# ---------------------------------------------------------------------------

_PREDICTORS = ["r_val", "n_k", "shd_truth", "shd_cpdag"]
_PRED_COLOUR = {"r_val": C_BLUE, "n_k": C_VERM, "shd_truth": C_GREEN, "shd_cpdag": C_ORANGE}
_PRED_MARKER = {"r_val": "o", "n_k": "s", "shd_truth": "^", "shd_cpdag": "D"}


def _stratum_label(row: dict) -> str:
    """A short, ordered label for one row's stratum, plip its own endpoint's n."""
    if row["arm"] == "flip":
        return f"flip  cov={float(row['coverage']):.1f}  bw={float(row['base_wrongness']):.2f}"
    return f"tiered  n_tiers={int(float(row['n_tiers']))}"


def _stratum_sort_key(row: dict) -> tuple:
    if row["arm"] == "flip":
        return (0, float(row["coverage"]), float(row["base_wrongness"]))
    return (1, float(row["n_tiers"]), 0.0)


def fig_tau_forest() -> list[Path]:
    """F3 — tau_b(predictor, AUC) with 95% CI, grouped by stratum, coloured by predictor.

    Only ``is_primary = True`` rows of ``analysis_tau_primary.csv`` for the four
    predictors ``r_val``, ``n_k``, ``shd_truth``, ``shd_cpdag``. Flip rows score
    against ``AUC_frac``; tiered rows score against ``AUC_frac_rate`` (the
    re-binned endpoint) — different columns, stated here and in the caption,
    never conflated. Exactly one cell is ``status = "undefined (predictor
    constant)"`` (flip, coverage=1.0, base_wrongness=0.0, shd_truth — where
    ``G0 == truth`` by construction): it is drawn as its own hollow-X marker
    with no CI, not silently dropped and not plotted at tau=0.
    """
    _style()
    rows = [r for r in _read_csv("analysis_tau_primary.csv") if r["is_primary"] == "True"]
    rows = [r for r in rows if r["predictor"] in _PREDICTORS]
    strata = sorted({_stratum_label(r) for r in rows}, key=lambda lbl: lbl)
    # Re-derive ordering from rows directly, so the display order matches design
    # order (flip cov=0.5 asc bw, then cov=1.0, then tiered) rather than alnum.
    order_rows = sorted({(_stratum_sort_key(r), _stratum_label(r)) for r in rows})
    strata = [lbl for _, lbl in order_rows]

    fig, ax = plt.subplots(figsize=(11.2, 7.4))
    n_pred = len(_PREDICTORS)
    band = 0.62
    offsets = {
        p: (i - (n_pred - 1) / 2) * (band / n_pred) for i, p in enumerate(_PREDICTORS)
    }

    undefined_handle = None
    for si, stratum in enumerate(strata):
        for row in rows:
            if _stratum_label(row) != stratum:
                continue
            p = row["predictor"]
            y = si + offsets[p]
            colour = _PRED_COLOUR[p]
            if row["status"] != "ok":
                ax.scatter(
                    [0.0], [y], marker="x", s=90, color=colour, linewidth=2.2, zorder=4
                )
                ax.text(
                    0.02, y, "undefined (predictor constant)", fontsize=6.6, color=colour,
                    va="center", ha="left", style="italic",
                )
                undefined_handle = Line2D(
                    [], [], ls="", marker="x", color="black", markersize=8,
                    label="undefined (predictor constant) — not tau=0, not omitted",
                )
                continue
            tau = float(row["tau_b"])
            lo, hi = float(row["ci_lo_2p5"]), float(row["ci_hi_97p5"])
            ax.plot([lo, hi], [y, y], color=colour, lw=1.6, zorder=2)
            ax.scatter(
                [tau], [y], s=48, color=colour, marker=_PRED_MARKER[p],
                edgecolor="white", linewidth=0.5, zorder=3,
            )

    ax.axvline(0.0, color=C_GREY, ls="--", lw=1.1, zorder=1)
    ax.set_yticks(range(len(strata)))
    ax.set_yticklabels(strata, fontsize=8.6)
    ax.invert_yaxis()
    ax.set_xlabel("tau_b(predictor, survival AUC), 95% CI  (flip: AUC_frac; tiered: AUC_frac_rate)")
    ax.set_xlim(-0.65, 1.0)
    ax.grid(axis="y", visible=False)
    ax.set_title(
        "r_val is the only predictor with a consistent sign across every stratum",
        fontsize=11.5,
    )

    handles = [
        Line2D([], [], color=_PRED_COLOUR[p], marker=_PRED_MARKER[p], lw=1.6, label=p)
        for p in _PREDICTORS
    ]
    if undefined_handle is not None:
        handles.append(undefined_handle)
    ax.legend(
        handles=handles,
        fontsize=8,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        framealpha=0.96,
        borderaxespad=0.0,
    )

    fig.text(
        0.5,
        -0.02,
        "analysis_tau_primary.csv, is_primary rows only, 10,000-resample bootstrap. " + RADIUS_NOTE,
        ha="center",
        fontsize=7,
        color=C_GREY,
    )
    return _save(fig, "session7_f3_tau_forest")


# ---------------------------------------------------------------------------
# F4 — Discordance spotlight.
# ---------------------------------------------------------------------------

_SPOTLIGHT_ID = "c06_s01_cov050_seed00000007_bw000"


def fig_discordance_spotlight() -> list[Path]:
    """F4 — the one verified discordance case: SHD 0, r_val 1, immediate collapse.

    ``c06_s01_cov050_seed00000007_bw000``: ``shd_truth = 0`` (``G0`` *is* the
    ground-truth DAG — the analyst's asserted knowledge is perfectly correct),
    ``n_k = 2``, ``r_val = 1``, ``AUC_frac = 0.0``. ``S(1) = 0.0`` on 107
    non-contradictory samples, ``S(2) = 0.0`` on 200 — every corruption
    invalidates the committed adjustment set, and ``r_val = 1`` predicted
    exactly that. Plotted against the median of its own stratum (flip,
    coverage=0.5, base_wrongness=0.00, r_val bucket "1", n=80 instances),
    restricted to depths with ``n_eval >= 30`` per instance.

    This is the *only* spotlight case in this session. A second one
    (``AUC_frac = 1.0`` at ``r_val = 8``) was withdrawn after verification
    because its AUC rested on depths with ``n_eval`` of 35, 10, 1 and 1
    (PREREGISTRATION Appendix F.1) — not reintroduced here, and no other case
    is substituted without the same check, which no other case has had.
    """
    _style()
    inst = _instances()
    spot = inst[_SPOTLIGHT_ID]
    assert spot["shd_truth"] == "0" and spot["n_k"] == "2" and spot["r_val"] == "1", (
        "spotlight instance no longer matches the report's description — check "
        "results/axis_robustness/survival_instances.csv"
    )

    stratum_ids = {
        iid
        for iid, r in inst.items()
        if r["arm"] == "flip"
        and r["coverage"] == "0.5"
        and r["base_wrongness"] == "0.0"
        and _bucket5(r["r_val"]) == "1"
    }

    spot_pts: dict[int, tuple[float, int]] = {}
    strat_vals: dict[int, list[float]] = defaultdict(list)
    for row in _curves():
        if row["arm"] != "flip" or row["S"] == "":
            continue
        d = int(row["d"])
        n_eval = int(row["n_eval"])
        if row["instance_id"] == _SPOTLIGHT_ID:
            spot_pts[d] = (float(row["S"]), n_eval)
        if row["instance_id"] in stratum_ids and n_eval >= N_EVAL_MIN:
            strat_vals[d].append(float(row["S"]))

    strat_med = {d: (st.median(vs), len(vs)) for d, vs in strat_vals.items() if len(vs) >= 5}

    fig, ax = plt.subplots(figsize=(7.6, 5.6))

    xs = sorted(strat_med)
    ys = [strat_med[d][0] for d in xs]
    ns = [strat_med[d][1] for d in xs]
    sizes = _marker_sizes(ns, lo=40, hi=220)
    ax.plot(xs, ys, "-", color=C_BLUE, lw=2.0, zorder=2, label="stratum median (flip, cov=0.5, bw=0.00, r_val=1)")
    ax.scatter(xs, ys, s=sizes, color=C_BLUE, marker="o", edgecolor="white", linewidth=0.6, zorder=3)
    for d in xs:
        ax.annotate(
            f"n={strat_med[d][1]} instances",
            xy=(d, strat_med[d][0]),
            xytext=(14, 10),
            textcoords="offset points",
            ha="left",
            fontsize=7,
            color=C_BLUE,
        )

    sxs = sorted(spot_pts)
    sys_ = [spot_pts[d][0] for d in sxs]
    ax.plot(sxs, sys_, "-", color=C_VERM, lw=2.4, zorder=4, label=f"{_SPOTLIGHT_ID}  (shd_truth=0, r_val=1)")
    ax.scatter(sxs, sys_, s=140, color=C_VERM, marker="X", edgecolor="white", linewidth=0.8, zorder=5)
    for d in sxs:
        n_eval = spot_pts[d][1]
        ax.annotate(
            f"S({d})=0.0\nn_eval={n_eval}",
            xy=(d, spot_pts[d][0]),
            xytext=(10, -26),
            textcoords="offset points",
            ha="left",
            fontsize=8,
            color=C_VERM,
            fontweight="bold",
            arrowprops={"arrowstyle": "->", "color": C_VERM, "lw": 0.9},
        )

    ax.set_xlabel("corruption depth  d  (claims reversed, flip arm)")
    ax.set_ylabel("S(d)")
    ax.set_xlim(0.5, max(max(xs, default=1), max(sxs, default=1)) + 0.9)
    ax.set_ylim(-0.06, 1.06)
    ax.grid(axis="x", visible=False)
    ax.legend(fontsize=8, loc="upper right", framealpha=0.96)
    ax.set_title(
        "Zero SHD is compatible with maximal fragility — r_val=1 predicted it, SHD did not",
        fontsize=11,
    )
    fig.text(
        0.5,
        -0.08,
        "shd_truth=0 here means G0 IS the ground-truth DAG: the analyst's background knowledge is\n"
        "perfectly correct. The committed adjustment set still fails on the very first corrupted claim.\n"
        "Marker area ∝ √(contributing instances); stratum median restricted to n_eval ≥ 30 per instance "
        "and n ≥ 5 instances per depth. " + RADIUS_NOTE,
        ha="center",
        fontsize=7,
        color=C_GREY,
    )
    return _save(fig, "session7_f4_discordance_spotlight")


# ---------------------------------------------------------------------------
# F5 — Self-revealing errors.
# ---------------------------------------------------------------------------


def fig_self_revealing() -> list[Path]:
    """F5 — most orientation errors are self-revealing: contradiction rate, flip arm.

    Left: contradiction rate vs corruption depth, flip arm only, pooled sample
    counts (``n_contradictory / (n_eval + n_contradictory)`` summed over all
    flip rows at each ``d``). Right: contradiction rate at ``d = 1`` split by
    base wrongness — 0.653 / 0.647 / 0.563 at ``b = 0.00 / 0.10 / 0.25``,
    reproduced exactly from ``survival_curves.csv`` and matching
    PREREGISTRATION Appendix D.1.

    The left panel's numbers do **not** match the table printed under this
    same caption in ``report_fragility_and_pareto.md`` section 2.7 (0.491 /
    0.773 / 0.896 / ... / 1.000 at d=1/2/3/.../12): that table is labelled
    "flip arm" but is arithmetically a pool of *both* arms using the tiered
    arm's defective ``d`` column (verified by reproducing it exactly when both
    arms are pooled). Pooling the two arms onto one ``d`` axis is the exact
    violation this session's non-negotiable rules forbid, so this figure plots
    the true flip-only curve instead. See the module docstring for the
    reproduction. The right panel's numbers are unaffected — they were
    flip-only in the report and remain so here.
    """
    _style()
    inst = _instances()
    flip_curves = [r for r in _curves() if r["arm"] == "flip"]

    by_d: dict[int, list[int]] = defaultdict(lambda: [0, 0])
    for row in flip_curves:
        d = int(row["d"])
        nc, ne = int(row["n_contradictory"]), int(row["n_eval"])
        by_d[d][0] += nc
        by_d[d][1] += nc + ne

    by_bw: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for row in flip_curves:
        if int(row["d"]) != 1:
            continue
        bw = inst[row["instance_id"]]["base_wrongness"]
        nc, ne = int(row["n_contradictory"]), int(row["n_eval"])
        by_bw[bw][0] += nc
        by_bw[bw][1] += nc + ne

    fig, (axl, axr) = plt.subplots(1, 2, figsize=(11.2, 4.8), gridspec_kw={"width_ratios": [1.5, 1]})

    ds = sorted(by_d)
    rates = [by_d[d][0] / by_d[d][1] for d in ds]
    ns = [by_d[d][1] for d in ds]
    axl.plot(ds, rates, "-o", color=C_VERM, lw=1.8, ms=5, zorder=2)
    for d, rate, n in zip(ds, rates, ns):
        if d in (1, 2, 3, 12) or d == max(ds):
            axl.annotate(
                f"{rate:.3f}",
                xy=(d, rate),
                xytext=(-4, 12),
                textcoords="offset points",
                ha="right",
                fontsize=7.6,
                color=C_VERM,
            )
    axl.set_xlabel("corruption depth  d  (claims reversed, flip arm)")
    axl.set_ylabel("contradiction rate  (pooled samples, flip arm only)")
    axl.set_ylim(0, 1.08)
    axl.set_xlim(0.5, max(ds) + 0.5)
    axl.grid(axis="x", visible=False)
    axl.set_title(
        f"Contradiction rate rises fast: {rates[0]:.3f} at d=1 → "
        f"{rates[ds.index(3)]:.3f} at d=3 → 1.000 by d={ds[-1]}",
        fontsize=10,
    )

    bws = sorted(by_bw, key=float)
    bw_rates = [by_bw[bw][0] / by_bw[bw][1] for bw in bws]
    xs = list(range(len(bws)))
    axr.bar(xs, bw_rates, 0.55, color=C_BLUE, edgecolor="white", lw=0.5)
    for x, bw, rate in zip(xs, bws, bw_rates):
        axr.text(x, rate + 0.02, f"{rate:.3f}", ha="center", fontsize=9)
        axr.text(
            x, 0.03, f"n={by_bw[bw][1]}", ha="center", va="bottom", fontsize=7.4, color="white",
            transform=axr.get_xaxis_transform(),
        )
    axr.set_xticks(xs)
    axr.set_xticklabels([f"b={float(bw):.2f}" for bw in bws])
    axr.set_ylim(0, 0.85)
    axr.set_xlabel("base wrongness  (flip arm, d = 1 only)")
    axr.set_ylabel("contradiction rate at d = 1")
    axr.grid(axis="x", visible=False)
    axr.set_title("Base wrongness barely moves it at d = 1", fontsize=10)

    fig.suptitle(
        "Most orientation errors are self-revealing: Meek closure fails and the analyst finds out",
        fontsize=11.5,
        y=1.04,
    )
    fig.text(
        0.5,
        -0.09,
        "Both panels: flip arm only, pooled sample counts from survival_curves.csv (never merged with\n"
        "the tiered arm's d axis — see module docstring for a report-table mislabelling this figure\n"
        "corrects). Right panel reproduces report_fragility_and_pareto.md §2.7 exactly (0.653 / 0.647 / 0.563).",
        ha="center",
        fontsize=7,
        color=C_GREY,
    )
    return _save(fig, "session7_f5_self_revealing")


def build_all() -> list[Path]:
    """Build every session 7 figure that has data.

    Returns:
        The paths written. A figure whose input file is missing or malformed is
        skipped with a message rather than aborting the run; an assertion
        failure (the spotlight instance no longer matching its description) is
        *not* caught, because that would mean the committed data disagrees with
        the report and must be investigated, not silently skipped.
    """
    out: list[Path] = []
    for fn in (
        fig_endpoint_by_bucket,
        fig_tiered_stratification,
        fig_tau_forest,
        fig_discordance_spotlight,
        fig_self_revealing,
        fig_survival_curves_by_bucket,
    ):
        try:
            out.extend(fn())
        except (FileNotFoundError, ValueError, IndexError, KeyError) as exc:
            print(f"skipped {fn.__name__}: {exc}")
    return out


if __name__ == "__main__":
    for p in build_all():
        print(p)
