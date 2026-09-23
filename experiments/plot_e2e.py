"""Section 5.3 figures for the end-to-end speedup sweep.

Reads a flattened ``e2e_speedup.csv`` (as written by ``experiments/e2e_speedup.py``,
gated or ungated) and writes three figures -- PDF + SVG (vector) and PNG (raster
preview) -- into ``<results_dir>/figures/``:

1. ``e2e_divergence``       -- median seconds vs k_undirected, one line per leg
   (space, search_fast, hybrid), IQR band where >= MIN_N_SUMMARY runs completed,
   individual runs elsewhere, censored-at-cap markers.
2. ``e2e_rescue_scatter``   -- search_fast seconds vs hybrid seconds, log-log,
   colored by hybrid method.
3. ``e2e_censoring_bars``   -- absolute censored-run counts per leg.

Every plotted number comes straight from the CSV: no smoothing, no fitted
curves, no rescaling of ratios. Deterministic (no randomness).

Usage:
    PYTHONPATH=src .venv/bin/python experiments/plot_e2e.py results/e2e_speedup_gated/e2e_speedup.csv
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

CAP_S = 120.0
MIN_N_SUMMARY = 3

# Okabe-Ito colorblind-safe palette.
OI = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "bluish_green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "reddish_purple": "#CC79A7",
    "grey": "#999999",
}

LEG_COLOR = {
    "space": OI["blue"],
    "search_fast": OI["vermillion"],
    "hybrid": OI["bluish_green"],
}
LEG_LABEL = {
    "space": "space (build + crit BFS)",
    "search_fast": "search_fast (unbounded)",
    "hybrid": "hybrid",
}
LEG_VALUE_COL = {
    "space": "space_total_crit_s",
    "search_fast": "search_fast_s",
    "hybrid": "hybrid_total_s",
}
LEG_CENSORED_COL = {
    "space": "space_censored",
    "search_fast": "search_fast_censored",
    "hybrid": "hybrid_censored",
}
LEG_ERROR_COL = {
    "space": "space_error",
    "search_fast": "search_fast_error",
    "hybrid": "hybrid_error",
}

CENSORED_LEGS = ("space", "search_frozen", "search_fast", "hybrid")
CENSORED_COL = {
    "space": ("space_censored", "space_error"),
    "search_frozen": ("search_frozen_censored", "search_frozen_error"),
    "search_fast": ("search_fast_censored", "search_fast_error"),
    "hybrid": ("hybrid_censored", "hybrid_error"),
}

METHOD_COLOR = {
    "local_up_fast": OI["sky_blue"],
    "e1_ladder": OI["orange"],
    "degenerate": OI["grey"],
}

FONT_BASE = 9.5

plt.rcParams.update(
    {
        "font.size": FONT_BASE,
        "axes.titlesize": FONT_BASE,
        "axes.labelsize": FONT_BASE,
        "xtick.labelsize": FONT_BASE - 0.5,
        "ytick.labelsize": FONT_BASE - 0.5,
        "legend.fontsize": FONT_BASE - 1.0,
        "font.family": "sans-serif",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }
)


# --------------------------------------------------------------------------
# CSV loading
# --------------------------------------------------------------------------


def _parse_bool(s: str) -> bool | None:
    if s is None or s == "":
        return None
    return s == "True"


def _parse_float(s: str) -> float | None:
    if s is None or s == "":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _parse_int(s: str) -> int | None:
    v = _parse_float(s)
    return None if v is None else int(v)


def load_rows(csv_path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with csv_path.open(newline="") as f:
        for raw in csv.DictReader(f):
            gate_exhausted = _parse_bool(raw.get("gate_exhausted", "")) or False
            if gate_exhausted:
                # Never ran any leg -- excluded from every figure below (they
                # carry no timing data by construction).
                continue
            row = {
                "n": _parse_int(raw.get("n")),
                "edge_prob": _parse_float(raw.get("edge_prob")),
                "seed": _parse_int(raw.get("seed")),
                "k_undirected": _parse_int(raw.get("k_undirected")),
                "space_total_crit_s": _parse_float(raw.get("space_total_crit_s")),
                "space_censored": _parse_bool(raw.get("space_censored")),
                "space_error": _parse_bool(raw.get("space_error")),
                "search_frozen_censored": _parse_bool(raw.get("search_frozen_censored")),
                "search_frozen_error": _parse_bool(raw.get("search_frozen_error")),
                "search_fast_s": _parse_float(raw.get("search_fast_s")),
                "search_fast_censored": _parse_bool(raw.get("search_fast_censored")),
                "search_fast_error": _parse_bool(raw.get("search_fast_error")),
                "hybrid_total_s": _parse_float(raw.get("hybrid_total_s")),
                "hybrid_censored": _parse_bool(raw.get("hybrid_censored")),
                "hybrid_error": _parse_bool(raw.get("hybrid_error")),
                "hybrid_method": raw.get("hybrid_method") or None,
            }
            rows.append(row)
    return rows


# --------------------------------------------------------------------------
# Figure 1: e2e_divergence
# --------------------------------------------------------------------------


def fig_e2e_divergence(rows: list[dict[str, Any]], out_stub: Path) -> str:
    ks = sorted({r["k_undirected"] for r in rows if r["k_undirected"] is not None})

    fig, ax = plt.subplots(figsize=(5.2, 4.2))

    legend_handles = []
    censored_marker_used = False
    sparse_points_used = False

    for leg in ("space", "search_fast", "hybrid"):
        val_col = LEG_VALUE_COL[leg]
        cens_col = LEG_CENSORED_COL[leg]
        err_col = LEG_ERROR_COL[leg]
        color = LEG_COLOR[leg]

        xs_line, ys_med, ys_lo, ys_hi = [], [], [], []
        xs_pts, ys_pts = [], []
        xs_cap, ys_cap = [], []

        for k in ks:
            at_k = [r for r in rows if r["k_undirected"] == k]
            total = len(at_k)
            completed = [r for r in at_k if not r[cens_col] and not r[err_col] and r[val_col] is not None]
            censored = [r for r in at_k if r[cens_col]]
            vals = sorted(r[val_col] for r in completed)

            # A median and IQR are drawn only where they summarise something: at
            # least MIN_N_SUMMARY completed runs, and a completed majority (else
            # the median is of the survivors and biased low). Sparser k are shown
            # as their individual runs, never joined into a trend.
            if len(vals) >= MIN_N_SUMMARY and len(vals) >= total / 2.0:
                q1, med, q3 = statistics.quantiles(vals, n=4, method="inclusive")
                xs_line.append(k)
                ys_med.append(med)
                ys_lo.append(q1)
                ys_hi.append(q3)
            else:
                xs_line.append(k)
                ys_med.append(float("nan"))
                ys_lo.append(float("nan"))
                ys_hi.append(float("nan"))
                xs_pts.extend([k] * len(vals))
                ys_pts.extend(vals)
            if censored:
                xs_cap.append(k)
                ys_cap.append(CAP_S)

        (line,) = ax.plot(
            xs_line, ys_med, marker="o", markersize=4, linewidth=1.6, color=color, label=LEG_LABEL[leg]
        )
        ax.fill_between(xs_line, ys_lo, ys_hi, color=color, alpha=0.18, linewidth=0)
        legend_handles.append(line)
        if xs_pts:
            ax.scatter(xs_pts, ys_pts, s=14, color=color, alpha=0.75, linewidth=0, zorder=4)
            sparse_points_used = True

        if xs_cap:
            ax.scatter(
                xs_cap,
                ys_cap,
                marker="v",
                s=42,
                facecolor="none",
                edgecolor=color,
                linewidth=1.4,
                zorder=5,
            )
            censored_marker_used = True

    ax.axhline(CAP_S, color=OI["black"], linestyle=":", linewidth=1.0, alpha=0.7)
    ax.text(
        ax.get_xlim()[1] if ks else 1,
        CAP_S,
        "  cap (120 s)",
        va="bottom",
        ha="right",
        fontsize=FONT_BASE - 1.5,
        color=OI["black"],
    )

    ax.set_yscale("log")
    ax.set_xlabel("k_undirected (CPDAG undirected edges)")
    ax.set_ylabel("seconds (log scale)")
    ax.grid(True, which="major", axis="y", alpha=0.25, linewidth=0.6)

    if sparse_points_used:
        legend_handles.append(
            Line2D(
                [], [], marker="o", linestyle="none", markersize=3.5, color=OI["black"], alpha=0.75,
                label=f"individual runs (k with < {MIN_N_SUMMARY} completed)",
            )
        )
    if censored_marker_used:
        legend_handles.append(
            Line2D(
                [], [], marker="v", linestyle="none", markerfacecolor="none",
                markeredgecolor=OI["black"], label="≥ 1 run censored at the cap",
            )
        )
    ax.legend(
        handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, -0.2),
        ncol=2, frameon=False,
    )

    fig.tight_layout()
    return _save(fig, out_stub)


# --------------------------------------------------------------------------
# Figure 2: e2e_rescue_scatter
# --------------------------------------------------------------------------


def fig_e2e_rescue_scatter(rows: list[dict[str, Any]], out_stub: Path) -> str:
    fig, ax = plt.subplots(figsize=(4.6, 4.2))

    methods_present: list[str] = []
    for r in rows:
        if r["hybrid_censored"] or r["hybrid_error"] or r["hybrid_total_s"] is None:
            continue  # hybrid didn't complete -- no y value to plot
        m = r["hybrid_method"] or "unknown"
        if m not in methods_present:
            methods_present.append(m)

        y = r["hybrid_total_s"]
        sf_censored = bool(r["search_fast_censored"])
        sf_error = bool(r["search_fast_error"])
        x = CAP_S if sf_censored else r["search_fast_s"]
        if sf_error or x is None:
            continue

        color = METHOD_COLOR.get(m, OI["reddish_purple"])
        if sf_censored:
            ax.scatter(x, y, marker=">", s=34, facecolor="none", edgecolor=color, linewidth=1.3, zorder=4)
        else:
            ax.scatter(x, y, marker="o", s=22, color=color, alpha=0.75, linewidth=0, zorder=3)

    lo = 1e-5
    hi = CAP_S * 1.6
    ax.plot([lo, hi], [lo, hi], color=OI["black"], linestyle="-", linewidth=0.8, alpha=0.5, zorder=1)
    ax.axvline(CAP_S, color=OI["black"], linestyle=":", linewidth=1.0, alpha=0.6)
    ax.axhline(CAP_S, color=OI["black"], linestyle=":", linewidth=1.0, alpha=0.6)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("search_fast seconds (log)")
    ax.set_ylabel("hybrid seconds (log)")
    ax.grid(True, which="major", alpha=0.2, linewidth=0.6)

    handles = [
        Line2D(
            [], [], marker="o", linestyle="none", color=METHOD_COLOR.get(m, OI["reddish_purple"]),
            label=f"hybrid method: {m}",
        )
        for m in methods_present
    ]
    handles.append(
        Line2D(
            [], [], marker=">", linestyle="none", markerfacecolor="none", markeredgecolor=OI["black"],
            label="search_fast censored (x = cap, lower bound)",
        )
    )
    handles.append(Line2D([], [], color=OI["black"], linewidth=0.8, alpha=0.5, label="y = x"))
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=FONT_BASE - 2.0)

    fig.tight_layout()
    return _save(fig, out_stub)


# --------------------------------------------------------------------------
# Figure 3: e2e_censoring_bars
# --------------------------------------------------------------------------


def fig_e2e_censoring_bars(rows: list[dict[str, Any]], out_stub: Path) -> str:
    total = len(rows)
    fig, ax = plt.subplots(figsize=(5.0, 3.4))

    xs = list(range(len(CENSORED_LEGS)))
    censored_counts = []
    error_counts = []
    for leg in CENSORED_LEGS:
        cens_col, err_col = CENSORED_COL[leg]
        c = sum(1 for r in rows if r[cens_col])
        e = sum(1 for r in rows if r[err_col])
        censored_counts.append(c)
        error_counts.append(e)

    bars_c = ax.bar(xs, censored_counts, color=OI["vermillion"], label="censored (≥ cap)")
    bars_e = ax.bar(
        xs, error_counts, bottom=censored_counts, color=OI["reddish_purple"], label="error"
    )

    for x, c, e in zip(xs, censored_counts, error_counts):
        if e > 0:
            ax.text(x, c + e / 2.0, str(e), ha="center", va="center", fontsize=FONT_BASE - 1.5, color="white")
        top = c + e
        ax.text(x, top + max(total * 0.01, 0.6), f"{top}/{total}", ha="center", va="bottom", fontsize=FONT_BASE - 1.5)

    ax.set_xticks(xs)
    ax.set_xticklabels(list(CENSORED_LEGS), rotation=20, ha="right")
    ax.set_ylabel(f"count censored / errored (of {total} rows)")
    ax.grid(True, which="major", axis="y", alpha=0.25, linewidth=0.6)
    top_max = max((c + e) for c, e in zip(censored_counts, error_counts))
    ax.set_ylim(0, top_max * 1.25 + 1)
    if any(error_counts):
        ax.legend(loc="upper right", frameon=False)

    fig.tight_layout()
    return _save(fig, out_stub)


# --------------------------------------------------------------------------
# Saving
# --------------------------------------------------------------------------


def _save(fig, out_stub: Path) -> str:
    out_stub.parent.mkdir(parents=True, exist_ok=True)
    pdf_path = out_stub.with_suffix(".pdf")
    svg_path = out_stub.with_suffix(".svg")
    png_path = out_stub.with_suffix(".png")
    fig.savefig(pdf_path)
    fig.savefig(svg_path)
    fig.savefig(png_path, dpi=200)
    plt.close(fig)
    return f"{pdf_path.name}, {svg_path.name}, {png_path.name}"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csv_path", help="path to e2e_speedup.csv")
    p.add_argument("--out-dir", default=None, help="figures output dir (default: <csv dir>/figures)")
    args = p.parse_args()

    csv_path = Path(args.csv_path).resolve()
    out_dir = Path(args.out_dir).resolve() if args.out_dir else csv_path.parent / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = load_rows(csv_path)
    print(f"[plot_e2e] loaded {len(rows)} resolved rows (gate_exhausted rows excluded) from {csv_path}")

    written = []
    written.append(fig_e2e_divergence(rows, out_dir / "e2e_divergence"))
    written.append(fig_e2e_rescue_scatter(rows, out_dir / "e2e_rescue_scatter"))
    written.append(fig_e2e_censoring_bars(rows, out_dir / "e2e_censoring_bars"))

    for w in written:
        print(f"[plot_e2e] wrote {w}")


if __name__ == "__main__":
    main()
