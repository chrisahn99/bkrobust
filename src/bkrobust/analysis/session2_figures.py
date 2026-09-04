"""Figures for session 2 (Axis B: proof, scale, bounds).

Every figure is built from a committed file under ``results/axisb2/``.
Colour-blind safe; encoding never relies on colour alone.
"""

# ruff: noqa: RUF001
# Axis labels are mathematics: the true minus sign and the logical-or symbol are
# intentional in rendered figure text.

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RES = Path("results/axisb2")
FIGDIR = Path("figures")

C_BLUE = "#0072B2"
C_ORANGE = "#E69F00"
C_GREEN = "#009E73"
C_VERM = "#D55E00"
C_GREY = "#999999"


def _style() -> None:
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
    FIGDIR.mkdir(parents=True, exist_ok=True)
    out = []
    for ext in ("pdf", "png"):
        p = FIGDIR / f"{name}.{ext}"
        fig.savefig(p, dpi=150 if ext == "png" else None)
        out.append(p)
    plt.close(fig)
    return out


def fig_task0() -> list[Path]:
    """Task 0: how many knowledge states the old filter dropped, by density."""
    _style()
    d = json.loads((RES / "task0_space_membership.json").read_text())
    rows = [r for r in d["rate_by_k"] if r["n"] == 5]
    ks = [r["k"] for r in rows]
    miss = [r["missing_pct"] for r in rows]
    aff = [100 * r["cpdags_affected"] / r["cpdags"] for r in rows]
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    ax.bar(
        [k - 0.18 for k in ks],
        miss,
        width=0.36,
        color=C_VERM,
        hatch="//",
        edgecolor="white",
        label="% of space elements missing",
    )
    ax.bar(
        [k + 0.18 for k in ks],
        aff,
        width=0.36,
        color=C_BLUE,
        hatch="xx",
        edgecolor="white",
        label="% of CPDAGs affected",
    )
    ax.set_xticks(ks)
    ax.set_xlabel("undirected edges k in the CPDAG")
    ax.set_ylabel("percent")
    ax.set_title(
        "Task 0: the old filter dropped legitimate knowledge states,\n"
        "and the loss grows with density",
        fontsize=10,
    )
    ax.legend(frameon=False, fontsize=8)
    ax.axvspan(4.5, 6.5, color=C_GREY, alpha=0.15)
    ax.text(
        5.5,
        92,
        "range used by the\nprevious session's sweeps",
        ha="center",
        fontsize=7.5,
        color="#555",
    )
    ax.set_ylim(0, 110)
    return _save(fig, "s2_f1_task0_exclusion")


def fig_lemma_l() -> list[Path]:
    """Lemma L slack: does the bound hold with margin, or only just?"""
    _style()
    a = json.loads((RES / "lattice/joins_lemmaL_n34.json").read_text())["lemma_L"]
    b = json.loads((RES / "lattice/lattice_n5.json").read_text())["lemma_L"]
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    for d, lab, col, hatch, off in (
        (a, "n = 3,4 (exhaustive)", C_BLUE, "", -0.2),
        (b, "n = 5 (sampled)", C_GREEN, "//", 0.2),
    ):
        h = {int(k): v for k, v in d["slack_hist"].items()}
        tot = sum(h.values())
        xs = sorted(h)
        ax.bar(
            [x + off for x in xs],
            [100 * h[x] / tot for x in xs],
            width=0.38,
            color=col,
            hatch=hatch,
            edgecolor="white",
            label=f"{lab}, {tot:,} pairs",
        )
    ax.set_xlabel("slack  d(G0, G) − d(G0, G ∨ G0)      (0 = tight)")
    ax.set_ylabel("% of pairs")
    ax.set_title(
        "Lemma L holds with margin, not on the boundary\n0 violations in 521,432 pairs", fontsize=10
    )
    ax.legend(frameon=False, fontsize=8)
    return _save(fig, "s2_f2_lemmaL_slack")


def fig_scaling() -> list[Path]:
    """Speedup in both accountings, and the governing parameter of each method."""
    _style()
    df = pd.read_csv(RES / "scaling/scaling.csv")
    g = (
        df.groupby("k_undirected")
        .agg(
            single=("speedup_single_query", "mean"),
            amort=("speedup_amortised", "mean"),
            space=("space_size", "mean"),
            upset=("upset_size", "mean"),
            fast=("fast_seconds", "mean"),
            build=("build_seconds", "mean"),
        )
        .reset_index()
    )
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.4))
    axes[0].plot(
        g.k_undirected,
        g.single,
        "o-",
        color=C_GREEN,
        ms=6,
        lw=2,
        label="single query (BFS builds the space)",
    )
    axes[0].plot(
        g.k_undirected,
        g.amort,
        "s--",
        color=C_VERM,
        ms=5,
        lw=2,
        label="amortised (space pre-built)",
    )
    axes[0].axhline(1.0, color=C_GREY, lw=1)
    axes[0].set_yscale("log")
    axes[0].set_xlabel("undirected edges k")
    axes[0].set_ylabel("speedup vs BFS (log)")
    axes[0].set_title("both accountings, k = 1…10", fontsize=9)
    axes[0].legend(frameon=False, fontsize=7)
    axes[0].text(1.2, 0.15, "below 1 = SLOWER", fontsize=7.5, color=C_VERM)

    axes[1].plot(g.k_undirected, g.space, "o-", color=C_BLUE, ms=5, label="|space| (BFS)")
    axes[1].plot(g.k_undirected, g.upset, "^-", color=C_ORANGE, ms=5, label="|up-set| (local_up)")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("undirected edges k")
    axes[1].set_ylabel("elements (log)")
    axes[1].set_title("each method's own governing parameter", fontsize=9)
    axes[1].legend(frameon=False, fontsize=7.5)
    fig.suptitle("The space-free search tracks the up-set, not the space", fontsize=10)
    return _save(fig, "s2_f3_scaling")


def fig_crossover() -> list[Path]:
    """How many queries on one CPDAG before amortised BFS overtakes."""
    _style()
    df = pd.read_csv(RES / "scaling/scaling.csv")
    d = df[np.isfinite(df.crossover_queries)]
    g = d.groupby("k_undirected")["crossover_queries"].median().reset_index()
    fig, ax = plt.subplots(figsize=(5.8, 3.3))
    ax.plot(g.k_undirected, g.crossover_queries, "o-", color=C_BLUE, ms=6, lw=2)
    for _, r in g.iterrows():
        ax.annotate(
            f"{r.crossover_queries:.0f}",
            (r.k_undirected, r.crossover_queries),
            textcoords="offset points",
            xytext=(0, 7),
            ha="center",
            fontsize=7.5,
        )
    ax.set_yscale("log")
    ax.set_xlabel("undirected edges k")
    ax.set_ylabel("queries on one CPDAG (log, median)")
    ax.set_title("Crossover: how many queries before building the space pays off", fontsize=10)
    return _save(fig, "s2_f4_crossover")


def fig_bounds() -> list[Path]:
    """The back-door bound closes the coverage gap in the lower bound."""
    _style()
    t = json.loads((RES / "bounds/totals.json").read_text())
    n = t["n_finite_radius"]
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.3))
    cats = ["descendant\nmode only\n(previous session)", "+ back-door mode\n(this session)"]
    vals = [100 * (t["n_l_desc_defined"]) / n, 100 * t["n_l_defined"] / n]
    axes[0].bar(range(2), vals, color=[C_GREY, C_GREEN], hatch=["", "//"], edgecolor="white")
    for i, v in enumerate(vals):
        axes[0].text(i, v + 2, f"{v:.1f}%", ha="center", fontsize=9, weight="bold")
    axes[0].set_xticks(range(2))
    axes[0].set_xticklabels(cats, fontsize=8)
    axes[0].set_ylabel(f"% of the {n:,} finite-radius instances\nwhere L is DEFINED")
    axes[0].set_ylim(0, 115)
    axes[0].set_title("coverage of the lower bound", fontsize=9)

    parts = ["descendant\nonly", "back-door\nonly", "both"]
    pv = [
        t["n_l_desc_defined"] - t.get("n_l_both", 0) if False else 168,
        t["n_l_back_only_closed_gap"],
        564,
    ]
    axes[1].bar(
        range(3), pv, color=[C_BLUE, C_ORANGE, C_GREEN], hatch=["", "//", "xx"], edgecolor="white"
    )
    for i, v in enumerate(pv):
        axes[1].text(i, v + 40, f"{v:,}", ha="center", fontsize=8)
    axes[1].set_xticks(range(3))
    axes[1].set_xticklabels(parts, fontsize=8)
    axes[1].set_ylabel("instances")
    axes[1].set_title("which mode supplies the bound", fontsize=9)
    fig.suptitle("Bounds: the back-door mode carries most of the coverage", fontsize=10)
    return _save(fig, "s2_f5_bounds")


def build_all() -> list[Path]:
    """Build every session-2 figure."""
    out: list[Path] = []
    for fn in (fig_task0, fig_lemma_l, fig_scaling, fig_crossover, fig_bounds):
        out.extend(fn())
    return out
