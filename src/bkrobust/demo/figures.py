"""Build the six report figures from :func:`~bkrobust.demo.pipeline.run_scenario` results.

Public entrypoint, called by ``run_all.py``::

    build_all_figures(results: dict[str, dict]) -> list[Path]

``results`` maps ``"A"``/``"B"``/``"C"`` to the dict :func:`~bkrobust.demo.pipeline.run_scenario`
returns. Every figure is written twice, as a vector ``.pdf`` (primary) and a 150 dpi
``.png`` (for embedding in ``report.md``), under ``figures/`` at the repository root
(a path relative to the current working directory, matching the
``results/breakdown_radius_demo/`` convention in :mod:`bkrobust.demo.export`).

Six figures, each with a matching ``..._caption`` function that states what the
reader should *take from* the figure, not merely what is drawn:

1. :func:`fig1_example` -- ground truth, CPDAG, and each scenario's G0, one fixed layout.
2. :func:`fig2_layered_space` -- the neighbour graph laid out by shell index around G0
   (THE central figure), one call per scenario.
3. :func:`fig3_shell_profile` -- shell index against graph/failure counts and bias.
4. :func:`fig4_witness` -- G0 beside the nearest Z-invalid graph, one call per scenario.
5. :func:`fig5_baseline_scatter` -- model-oriented radius against a naive K-count radius.
   Depends on ``bkrobust.demo.baseline`` (built by a different subagent); imported
   lazily and skipped gracefully if unavailable.
6. :func:`fig6_robustness_frontier` -- r_val against asymptotic variance across every
   valid adjustment set of G0, the efficiency-vs-robustness tradeoff.

A single colour-blind-safe palette (Okabe & Ito, 2008) is used throughout, and every
figure that encodes a binary (Z-valid/invalid, optimal/not) doubles it with marker
shape so colour alone never carries the distinction.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch

from bkrobust.demo.export import elements_frame, shells_frame
from bkrobust.demo.graph import MPDAG, canon
from bkrobust.demo.pipeline import UNREACHED, robustness_frontier
from bkrobust.demo.scenario import OUTCOME, TREATMENT

logger = logging.getLogger(__name__)

#: Default output directory, relative to the current working directory -- matches
#: ``bkrobust.demo.export.RESULTS_ROOT``'s convention.
FIGURES_ROOT = Path("figures")

SCENARIO_LABELS: tuple[str, ...] = ("A", "B", "C")

# --------------------------------------------------------------------------------
# palette (Okabe & Ito 2008 -- colour-blind safe)
# --------------------------------------------------------------------------------

OKABE_ITO = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "grey": "#999999",
}

COLOR_TREATMENT = OKABE_ITO["vermillion"]
COLOR_OUTCOME = OKABE_ITO["green"]
COLOR_COVARIATE = OKABE_ITO["sky_blue"]
COLOR_DIRECTED_EDGE = "#333333"
COLOR_UNDIRECTED_EDGE = OKABE_ITO["grey"]
COLOR_VALID = OKABE_ITO["blue"]
COLOR_INVALID = OKABE_ITO["vermillion"]
COLOR_G0 = OKABE_ITO["purple"]
COLOR_WITNESS = OKABE_ITO["orange"]
COLOR_HIGHLIGHT = OKABE_ITO["vermillion"]
COLOR_OPTIMAL = OKABE_ITO["orange"]

MARKER_TREATMENT = "s"
MARKER_OUTCOME = "D"
MARKER_COVARIATE = "o"
MARKER_VALID = "o"
MARKER_INVALID = "X"

plt.rcParams.update(
    {
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
    }
)


# --------------------------------------------------------------------------------
# small shared helpers
# --------------------------------------------------------------------------------


def _save(fig: plt.Figure, name: str, out_dir: Path) -> list[Path]:
    """Write ``fig`` as both ``.pdf`` and ``.png`` under ``out_dir``; return both paths."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / f"{name}.pdf"
    png_path = out_dir / f"{name}.png"
    # bbox_inches="tight" so a figure-level legend or suptitle placed outside
    # the axes (tight_layout does not account for either) never gets clipped.
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return [pdf_path, png_path]


def _skeleton_graph(g: MPDAG) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(g.nodes)
    graph.add_edges_from(g.skeleton())
    return graph


def compute_layout(cpdag: MPDAG, seed: int = 7) -> dict[str, tuple[float, float]]:
    """A single fixed node layout for ``cpdag``'s skeleton, reused across every panel.

    All three scenarios share the same CPDAG (and hence the same skeleton), so one
    layout, computed once, lets every graph drawn in this module be compared
    position-for-position.
    """
    graph = _skeleton_graph(cpdag)
    n = max(len(graph.nodes), 1)
    pos = nx.spring_layout(graph, seed=seed, k=1.6 / n**0.5)
    return {node: (float(x), float(y)) for node, (x, y) in pos.items()}


def _edge_state(g: MPDAG, a: str, b: str) -> tuple[str, str] | None:
    """The oriented direction of the ``a``/``b`` edge in ``g``, or ``None`` if undirected."""
    if g.is_directed_edge(a, b):
        return (a, b)
    if g.is_directed_edge(b, a):
        return (b, a)
    return None


def _differing_edges(a: MPDAG, b: MPDAG) -> set[tuple[str, str]]:
    """Canonical edges whose orientation state differs between ``a`` and ``b``.

    Assumes ``a`` and ``b`` share a skeleton (true of every pair of graphs compared
    in this module -- all are refinements of the same CPDAG).
    """
    diffs: set[tuple[str, str]] = set()
    for u, v in a.skeleton() | b.skeleton():
        if _edge_state(a, u, v) != _edge_state(b, u, v):
            diffs.add(canon(u, v))
    return diffs


def _find_graph_by_edge_string(space: list[MPDAG], edge_string: str) -> MPDAG | None:
    for g in space:
        if g.edge_string() == edge_string:
            return g
    return None


def _node_style(node: str, x: str, y: str) -> tuple[str, str]:
    """(colour, marker) for a node, with X and Y distinguished from covariates."""
    if node == x:
        return COLOR_TREATMENT, MARKER_TREATMENT
    if node == y:
        return COLOR_OUTCOME, MARKER_OUTCOME
    return COLOR_COVARIATE, MARKER_COVARIATE


def _draw_graph(
    ax: plt.Axes,
    g: MPDAG,
    pos: dict[str, tuple[float, float]],
    title: str = "",
    x: str = TREATMENT,
    y: str = OUTCOME,
    highlight_edges: set[tuple[str, str]] | None = None,
) -> None:
    """Draw one MPDAG on ``ax`` using the shared ``pos`` layout.

    Directed edges are arrows in a dark neutral colour; undirected (CPDAG) edges are
    thick grey lines with no arrowhead, visually distinct from directed edges without
    relying on the reader to notice a missing arrowhead alone. X and Y are drawn with
    a different colour *and* a different marker shape from covariates. Edges in
    ``highlight_edges`` (canonical ``(a, b)`` pairs) are drawn in a strong highlight
    colour, thicker, regardless of their own directed/undirected styling.
    """
    highlight_edges = highlight_edges or set()

    for a, b in g.undirected_edges:
        is_hl = canon(a, b) in highlight_edges
        (xa, ya), (xb, yb) = pos[a], pos[b]
        ax.plot(
            [xa, xb],
            [ya, yb],
            color=COLOR_HIGHLIGHT if is_hl else COLOR_UNDIRECTED_EDGE,
            linewidth=3.4 if is_hl else 2.6,
            solid_capstyle="round",
            alpha=0.95 if is_hl else 0.75,
            zorder=1,
        )

    for a, b in g.directed_edges:
        is_hl = canon(a, b) in highlight_edges
        arrow = FancyArrowPatch(
            pos[a],
            pos[b],
            arrowstyle="-|>",
            mutation_scale=14,
            shrinkA=13,
            shrinkB=13,
            linewidth=2.2 if is_hl else 1.3,
            color=COLOR_HIGHLIGHT if is_hl else COLOR_DIRECTED_EDGE,
            zorder=2,
        )
        ax.add_patch(arrow)

    for node in g.nodes:
        px, py = pos[node]
        color, marker = _node_style(node, x, y)
        is_special = node in (x, y)
        ax.scatter(
            [px],
            [py],
            s=260 if is_special else 190,
            c=color,
            marker=marker,
            edgecolors="black",
            linewidths=0.9,
            zorder=3,
        )
        ax.annotate(
            node,
            (px, py),
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=7.5 if is_special else 6.8,
            fontweight="bold" if is_special else "normal",
            zorder=4,
        )

    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_aspect("equal")
    xs = [p[0] for p in pos.values()]
    ys = [p[1] for p in pos.values()]
    pad = 0.28
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_ylim(min(ys) - pad, max(ys) + pad)


def _graph_legend_handles(x: str = TREATMENT, y: str = OUTCOME) -> list[Any]:
    return [
        Line2D([0], [0], color=COLOR_DIRECTED_EDGE, lw=1.6, label="directed edge (oriented)"),
        Line2D(
            [0],
            [0],
            color=COLOR_UNDIRECTED_EDGE,
            lw=2.8,
            label="undirected edge (CPDAG, unresolved)",
        ),
        Line2D(
            [0],
            [0],
            marker=MARKER_TREATMENT,
            color="none",
            markerfacecolor=COLOR_TREATMENT,
            markeredgecolor="black",
            markersize=9,
            label=f"X = {x} (treatment)",
        ),
        Line2D(
            [0],
            [0],
            marker=MARKER_OUTCOME,
            color="none",
            markerfacecolor=COLOR_OUTCOME,
            markeredgecolor="black",
            markersize=9,
            label=f"Y = {y} (outcome)",
        ),
        Line2D(
            [0],
            [0],
            marker=MARKER_COVARIATE,
            color="none",
            markerfacecolor=COLOR_COVARIATE,
            markeredgecolor="black",
            markersize=9,
            label="covariate",
        ),
    ]


# --------------------------------------------------------------------------------
# figure 1 -- the example
# --------------------------------------------------------------------------------


def fig1_example(results: dict[str, dict], out_dir: Path) -> list[Path]:
    """Ground truth, CPDAG, and each scenario's G0, side by side on one fixed layout."""
    any_result = results["A"]
    cpdag = any_result["cpdag"]
    truth = any_result["truth"]
    pos = compute_layout(cpdag)

    fig, axes = plt.subplots(1, 5, figsize=(21, 4.6))
    _draw_graph(axes[0], truth, pos, title="Ground-truth DAG")
    _draw_graph(axes[1], cpdag, pos, title="Estimated CPDAG")
    for ax, label in zip(axes[2:], SCENARIO_LABELS):  # noqa: B905 - equal by construction; strict= is 3.10+
        g0 = results[label]["g0"]
        _draw_graph(ax, g0, pos, title=f"G0 -- scenario {label}")

    handles = _graph_legend_handles()
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=len(handles),
        frameon=False,
        bbox_to_anchor=(0.5, -0.04),
    )
    fig.suptitle(
        "The running example: ground truth, CPDAG, and each scenario's starting graph G0", y=1.03
    )
    fig.tight_layout()
    return _save(fig, "fig1_example", out_dir)


def fig1_example_caption(results: dict[str, dict]) -> str:
    """Caption for :func:`fig1_example`."""
    parts = [f"{lab}: {results[lab]['spec']['relation']}" for lab in SCENARIO_LABELS]
    return (
        "Panel 1 is the ground-truth DAG and panel 2 the CPDAG estimated from it; panels 3-5 "
        "are each scenario's starting graph G0, drawn on the identical fixed layout used "
        "everywhere in this module so edges can be compared position-for-position across "
        "panels. Take from this that all three scenarios share the same skeleton and differ "
        "only in how the undirected CPDAG component is resolved: "
        + "; ".join(parts)
        + " -- that resolved orientation is the single perturbation the rest of the figures "
        "trace the consequences of."
    )


# --------------------------------------------------------------------------------
# figure 2 -- THE central figure: the layered perturbation space
# --------------------------------------------------------------------------------


def _shell_layout(space: list[MPDAG], shells: dict[MPDAG, int]) -> dict[MPDAG, tuple[float, float]]:
    """Concentric layout: shell 0 at the centre, shell k at radius k."""
    by_shell: dict[int, list[MPDAG]] = {}
    for g in space:
        s = shells.get(g)
        if s is None:
            continue
        by_shell.setdefault(s, []).append(g)

    pos: dict[MPDAG, tuple[float, float]] = {}
    for s, members in by_shell.items():
        ordered = sorted(members, key=lambda g: g.edge_string())
        if s == 0:
            for g in ordered:
                pos[g] = (0.0, 0.0)
            continue
        n = len(ordered)
        stagger = 0.37 * s
        for i, g in enumerate(ordered):
            theta = 2.0 * np.pi * i / n + stagger
            pos[g] = (s * np.cos(theta), s * np.sin(theta))
    return pos


def fig2_layered_space(result: dict[str, Any], label: str, out_dir: Path) -> list[Path]:
    """The neighbour graph laid out by shell index around G0, for one scenario.

    Nodes are coloured (and shaped) by Z-validity, G0 is marked distinctly at the
    centre, and the first Z-invalid element (the r_val witness) is marked and
    annotated so the location of the first failure is immediately legible.
    """
    space = result["space"]
    shells = result["shells"]
    covers = result["covers"]
    g0 = result["g0"]
    rows_by_edge = {
        g.edge_string(): row
        for g, row in zip(space, result["rows"])  # noqa: B905 - equal by construction; strict= is 3.10+
    }
    pos = _shell_layout(space, shells)
    max_shell = max(shells.values()) if shells else 0

    witness_info = result["witnesses"]["r_val"]
    witness_g = (
        _find_graph_by_edge_string(space, witness_info["edge_string"])
        if witness_info is not None
        else None
    )

    fig, ax = plt.subplots(figsize=(8.5, 8.5))

    for s in range(max_shell + 1):
        ax.add_patch(Circle((0, 0), s, fill=False, ls="--", lw=0.7, color="#d5d5d5", zorder=0))
        if s > 0:
            ax.annotate(
                f"shell {s}",
                (s * np.cos(np.pi / 2), s * np.sin(np.pi / 2)),
                fontsize=6.5,
                color="#9a9a9a",
                ha="left",
                va="bottom",
            )

    for lower, upper in covers:
        if lower in pos and upper in pos:
            (x1, y1), (x2, y2) = pos[lower], pos[upper]
            ax.plot([x1, x2], [y1, y2], color="#e2e2e2", lw=0.6, zorder=1)

    for g in space:
        if g not in pos or g is g0 or g is witness_g:
            continue
        row = rows_by_edge[g.edge_string()]
        valid = bool(row["z_valid"])
        px, py = pos[g]
        ax.scatter(
            [px],
            [py],
            s=42,
            c=COLOR_VALID if valid else COLOR_INVALID,
            marker=MARKER_VALID if valid else MARKER_INVALID,
            alpha=0.82,
            linewidths=0.4,
            edgecolors="white" if valid else "none",
            zorder=3,
        )

    gx, gy = pos[g0]
    ax.scatter(
        [gx],
        [gy],
        s=420,
        marker="*",
        c=COLOR_G0,
        edgecolors="black",
        linewidths=1.0,
        zorder=5,
        label="G0 (analyst's graph)",
    )

    if witness_g is not None:
        wx, wy = pos[witness_g]
        w_valid = bool(rows_by_edge[witness_g.edge_string()]["z_valid"])
        ax.scatter(
            [wx],
            [wy],
            s=210,
            marker=MARKER_VALID if w_valid else MARKER_INVALID,
            c=COLOR_VALID if w_valid else COLOR_INVALID,
            edgecolors=COLOR_WITNESS,
            linewidths=2.6,
            zorder=6,
        )
        ax.annotate(
            f"first Z-failure\n(r_val = {witness_info['shell']})",
            (wx, wy),
            xytext=(14, 14),
            textcoords="offset points",
            fontsize=8.5,
            fontweight="bold",
            color=COLOR_WITNESS,
            arrowprops={"arrowstyle": "-", "color": COLOR_WITNESS, "lw": 1.0},
            zorder=6,
        )

    handles = [
        Line2D(
            [0],
            [0],
            marker="*",
            color="none",
            markerfacecolor=COLOR_G0,
            markeredgecolor="black",
            markersize=16,
            label="G0 (analyst's graph)",
        ),
        Line2D(
            [0],
            [0],
            marker=MARKER_VALID,
            color="none",
            markerfacecolor=COLOR_VALID,
            markeredgecolor="white",
            markersize=8,
            label="Z valid",
        ),
        Line2D(
            [0],
            [0],
            marker=MARKER_INVALID,
            color="none",
            markerfacecolor=COLOR_INVALID,
            markeredgecolor="none",
            markersize=9,
            label="Z invalid",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="none",
            markeredgecolor=COLOR_WITNESS,
            markersize=11,
            markeredgewidth=2.2,
            label="r_val witness",
        ),
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    n_reachable = sum(1 for g in space if g in shells)
    ax.set_title(
        f"Scenario {label}: the perturbation space, shells around G0\n"
        f"({len(space)} graphs enumerated, {n_reachable} reachable from G0)"
    )
    fig.tight_layout()
    return _save(fig, f"fig2_layered_space_{label}", out_dir)


def fig2_layered_space_caption(result: dict[str, Any], label: str) -> str:
    """Caption for :func:`fig2_layered_space`."""
    radii = result["radii"]
    r_val = radii["r_val"]
    witness_info = result["witnesses"]["r_val"]
    if r_val == UNREACHED or witness_info is None:
        fail_txt = "no element of the enumerated space breaks Z's validity -- r_val is UNREACHED"
    else:
        moves = ", ".join(witness_info["moves"])
        fail_txt = (
            f"the nearest Z-invalid graph sits at shell {r_val} of the covering-relation "
            f"neighbour graph, reached from G0 by {moves}"
        )
    return (
        f"Scenario {label}: the neighbour graph of MPDAG refinements of the CPDAG, laid out in "
        f"concentric shells by covering-relation distance from G0 (the starred centre). Blue "
        f"circles are elements where the adjustment set Z stays valid; orange X markers are "
        f"elements where it breaks. Take from this that {fail_txt} -- that shell count is "
        f"r_val, the breakdown radius, the headline robustness number reported for scenario "
        f"{label}."
    )


# --------------------------------------------------------------------------------
# figure 3 -- shell profile
# --------------------------------------------------------------------------------


def fig3_shell_profile(results: dict[str, dict], out_dir: Path) -> list[Path]:
    """Shell index against graph counts, Z-invalid counts, and mean/max bias, all scenarios."""
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.8), sharey=False)

    bar_handles: list[Any] = []
    line_handles: list[Any] = []

    for ax, label in zip(axes, SCENARIO_LABELS):  # noqa: B905 - equal by construction; strict= is 3.10+
        result = results[label]
        elements = elements_frame(result)
        shells_df = shells_frame(elements).sort_values("shell")
        x = shells_df["shell"].to_numpy()
        width = 0.38

        b1 = ax.bar(
            x - width / 2,
            shells_df["n_graphs"],
            width=width,
            color=OKABE_ITO["sky_blue"],
            label="graphs in shell",
            zorder=2,
        )
        b2 = ax.bar(
            x + width / 2,
            shells_df["n_z_invalid"],
            width=width,
            color=OKABE_ITO["vermillion"],
            label="Z invalid in shell",
            hatch="//",
            zorder=2,
        )
        ax.set_xlabel("shell index (distance from G0)")
        ax.set_ylabel("number of graphs")
        ax.set_xticks(x)

        ax2 = ax.twinx()
        (l1,) = ax2.plot(
            x,
            shells_df["mean_abs_bias"],
            color=OKABE_ITO["black"],
            marker="o",
            ms=4.5,
            lw=1.5,
            label="mean |bias|",
            zorder=3,
        )
        (l2,) = ax2.plot(
            x,
            shells_df["max_abs_bias"],
            color=OKABE_ITO["orange"],
            marker="^",
            ms=4.5,
            lw=1.3,
            ls="--",
            label="max |bias|",
            zorder=3,
        )
        ax2.set_ylabel("|bias| of the effect estimate")

        r_val = result["radii"]["r_val"]
        if r_val != UNREACHED:
            ax.axvline(r_val, color=COLOR_WITNESS, lw=1.4, ls=":", zorder=1)

        ax.set_title(f"scenario {label}  (r_val={'UNREACHED' if r_val == UNREACHED else r_val})")

        if not bar_handles:
            bar_handles = [b1, b2]
            line_handles = [l1, l2]

    handles = list(bar_handles) + list(line_handles)
    labels = [h.get_label() for h in handles]
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=len(handles),
        frameon=False,
        bbox_to_anchor=(0.5, -0.06),
    )
    fig.suptitle(
        "Shell profile: how many graphs, how many Z-invalid, and how much bias, per shell", y=1.03
    )
    fig.tight_layout()
    return _save(fig, "fig3_shell_profile", out_dir)


def fig3_shell_profile_caption(results: dict[str, dict]) -> str:
    """Caption for :func:`fig3_shell_profile`."""
    bits = []
    for label in SCENARIO_LABELS:
        r_val = results[label]["radii"]["r_val"]
        bits.append(f"{label} (r_val={'UNREACHED' if r_val == UNREACHED else r_val})")
    return (
        "Bars give the number of graphs per shell and how many have an invalid adjustment set "
        "Z; the line (twin axis) gives mean and max absolute bias per shell. Take from this "
        "that bias sits at floating-point noise until the shell where Z first turns invalid "
        "(marked by the dotted vertical line), then jumps -- it is the validity failure "
        "itself, not a gradual drift, that drives bias in the effect estimate. Breakdown "
        "radii across scenarios: " + ", ".join(bits) + "."
    )


# --------------------------------------------------------------------------------
# figure 4 -- the witness
# --------------------------------------------------------------------------------


def fig4_witness(result: dict[str, Any], label: str, out_dir: Path) -> list[Path]:
    """G0 beside the nearest Z-invalid graph, differing edges highlighted, moves annotated."""
    g0 = result["g0"]
    witness_info = result["witnesses"]["r_val"]
    pos = compute_layout(result["cpdag"])

    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
    _draw_graph(axes[0], g0, pos, title=f"G0 -- scenario {label}")

    if witness_info is None:
        axes[1].axis("off")
        axes[1].text(
            0.5,
            0.5,
            "r_val is UNREACHED:\nno Z-invalid graph found\nin the enumerated space",
            ha="center",
            va="center",
            fontsize=10,
            transform=axes[1].transAxes,
        )
        moves_text = "r_val UNREACHED -- no witness graph exists"
    else:
        witness_g = _find_graph_by_edge_string(result["space"], witness_info["edge_string"])
        diff_edges = _differing_edges(g0, witness_g) if witness_g is not None else set()
        _draw_graph(
            axes[1],
            witness_g,
            pos,
            title=f"first Z-invalid graph -- shell {witness_info['shell']}",
            highlight_edges=diff_edges,
        )
        moves_text = "moves from G0: " + "; ".join(witness_info["moves"])

    handles = [
        *_graph_legend_handles(),
        Line2D([0], [0], color=COLOR_HIGHLIGHT, lw=3, label="edge that differs from G0"),
    ]
    fig.legend(
        handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.10)
    )
    fig.suptitle(f"Witness for r_val -- scenario {label}\n{moves_text}", y=1.05, fontsize=9.5)
    fig.tight_layout()
    return _save(fig, f"fig4_witness_{label}", out_dir)


def fig4_witness_caption(result: dict[str, Any], label: str) -> str:
    """Caption for :func:`fig4_witness`."""
    witness_info = result["witnesses"]["r_val"]
    if witness_info is None:
        return (
            f"Scenario {label}: r_val is UNREACHED in the enumerated space, so no Z-invalid "
            f"witness graph exists to display -- Z stays a valid adjustment set at every graph "
            f"reachable from G0."
        )
    moves = "; ".join(witness_info["moves"])
    return (
        f"Scenario {label}: G0 (left) beside the nearest graph at which the adjustment set Z "
        f"is no longer valid (right, shell {witness_info['shell']}), with the edges that "
        f"differ from G0 drawn in a highlighted colour. Take from this the concrete, minimal "
        f"reorientation an analyst could be talked into -- {moves} -- that would already break "
        f"the adjustment set they are relying on."
    )


# --------------------------------------------------------------------------------
# figure 5 -- naive-vs-model baseline scatter
# --------------------------------------------------------------------------------


def fig5_baseline_scatter(results: dict[str, dict], out_dir: Path) -> list[Path]:
    """Model-oriented radius (shell index) against a naive K-count radius, all scenarios.

    Depends on ``bkrobust.demo.baseline.naive_vs_model_frame``, built separately.
    Imported lazily; if the module or the expected columns are unavailable, this
    figure is skipped (logged, not raised) so the rest of the module stays testable
    independent of that dependency landing.
    """
    try:
        from bkrobust.demo.baseline import naive_vs_model_frame
    except ImportError as exc:
        logger.warning(
            "bkrobust.demo.baseline unavailable (%s); skipping fig5_baseline_scatter", exc
        )
        return []

    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharex=False, sharey=False)
    rng = np.random.default_rng(0)
    any_plotted = False

    for ax, label in zip(axes, SCENARIO_LABELS):  # noqa: B905 - equal by construction; strict= is 3.10+
        result = results[label]
        try:
            df = naive_vs_model_frame(result)
        except Exception as exc:  # a broken dependency must not crash this module
            logger.warning(
                "naive_vs_model_frame(results[%r]) failed (%s); skipping panel", label, exc
            )
            ax.axis("off")
            ax.set_title(f"scenario {label} (unavailable)")
            continue

        required = {"model_distance", "naive_distance"}
        if df is None or len(df) == 0 or not required.issubset(df.columns):
            logger.warning(
                "naive_vs_model_frame(results[%r]) missing required columns %s; skipping panel",
                label,
                sorted(required),
            )
            ax.axis("off")
            ax.set_title(f"scenario {label} (unavailable)")
            continue

        # Both distance columns can in principle carry the pipeline's UNREACHED
        # sentinel (-1) for a space element outside G0's connected component.
        # Never plot that sentinel as a bare negative coordinate -- drop those
        # rows here, visibly, rather than silently.
        n_before = len(df)
        df = df[(df["model_distance"] >= 0) & (df["naive_distance"] >= 0)]
        n_dropped = n_before - len(df)
        if n_dropped:
            logger.warning(
                "scenario %r: dropped %d/%d rows with an UNREACHED (-1) distance from "
                "fig5_baseline_scatter",
                label,
                n_dropped,
                n_before,
            )
        if df.empty:
            ax.axis("off")
            ax.set_title(f"scenario {label} (all rows UNREACHED)")
            continue

        any_plotted = True
        naive = df["naive_distance"].to_numpy(dtype=float)
        model = df["model_distance"].to_numpy(dtype=float)
        jitter_x = naive + rng.uniform(-0.12, 0.12, len(df))
        jitter_y = model + rng.uniform(-0.12, 0.12, len(df))
        if "z_valid" in df.columns:
            valid = df["z_valid"].to_numpy(dtype=bool)
        else:
            valid = np.ones(len(df), dtype=bool)

        ax.scatter(
            jitter_x[valid],
            jitter_y[valid],
            c=COLOR_VALID,
            marker=MARKER_VALID,
            s=32,
            alpha=0.7,
            edgecolors="none",
            label="Z valid",
            zorder=3,
        )
        ax.scatter(
            jitter_x[~valid],
            jitter_y[~valid],
            c=COLOR_INVALID,
            marker=MARKER_INVALID,
            s=36,
            alpha=0.85,
            edgecolors="none",
            label="Z invalid",
            zorder=3,
        )
        lim = float(max(naive.max(), model.max(), 1.0)) + 1.0
        ax.plot([0, lim], [0, lim], color="#999999", ls="--", lw=1.0, zorder=1, label="y = x")
        ax.set_xlabel("naive K-count radius")
        ax.set_ylabel("model-oriented radius (shell index)")
        ax.set_title(f"scenario {label}")

    if not any_plotted:
        plt.close(fig)
        logger.warning("fig5_baseline_scatter: no scenario produced usable data; nothing written")
        return []

    handles, labels_ = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            labels_,
            loc="lower center",
            ncol=len(handles),
            frameon=False,
            bbox_to_anchor=(0.5, -0.06),
        )
    fig.suptitle("Model-oriented breakdown distance vs. a naive K-count distance", y=1.03)
    fig.tight_layout()
    return _save(fig, "fig5_baseline_scatter", out_dir)


def fig5_baseline_scatter_caption(results: dict[str, dict]) -> str:
    """Caption for :func:`fig5_baseline_scatter`."""
    return (
        "Model-oriented breakdown radius (shell index in the covering-relation neighbour graph) "
        "plotted against a naive K-count radius (a simple tally of edges changed from the "
        "analyst's asserted knowledge), one jittered point per space element per scenario, with "
        "the y=x reference line. Take from this whether the two notions of 'distance from G0' "
        "agree: points off the diagonal are graphs the naive count would rank differently -- "
        "closer or farther -- than the model-consistent covering-relation metric this "
        "demonstration argues for."
    )


# --------------------------------------------------------------------------------
# figure 6 -- the robustness frontier
# --------------------------------------------------------------------------------


def frontier_plot_frame(result: dict[str, Any], n_draws: int = 25) -> pd.DataFrame:
    """The robustness-frontier data, prepared for plotting.

    Adds a ``_y`` column that maps ``r_val`` to a plot-ready value: real (reachable)
    radii pass through unchanged, and every row where ``r_val == UNREACHED`` (``-1``)
    is instead mapped to a distinct, labelled "never fails" tier strictly above every
    real shell index -- so UNREACHED is never plotted as a bare ``-1`` on a numeric
    axis. ``never_tier`` (a scalar, repeated on every row) records that tier's value.
    """
    frontier = robustness_frontier(result, n_draws=n_draws)
    df = pd.DataFrame(frontier)
    max_shell = max(result["shells"].values()) if result["shells"] else 0
    never_tier = max_shell + 2
    df["never_fails"] = df["r_val"] == UNREACHED
    df["_y"] = df["r_val"].where(~df["never_fails"], never_tier)
    df["never_tier"] = never_tier
    df["max_shell"] = max_shell
    return df


def fig6_robustness_frontier(results: dict[str, dict], out_dir: Path) -> list[Path]:
    """r_val against asymptotic variance for every valid adjustment set of G0, all scenarios.

    UNREACHED r_val values are plotted on a separate, explicitly labelled "never
    fails" tier rather than as a numeric ``-1``. The optimal adjustment set O*(G0)
    is marked distinctly (star marker, not colour alone).
    """
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2), sharey=False)

    for ax, label in zip(axes, SCENARIO_LABELS):  # noqa: B905 - equal by construction; strict= is 3.10+
        result = results[label]
        df = frontier_plot_frame(result)
        max_shell = int(df["max_shell"].iloc[0]) if len(df) else 0
        never_tier = int(df["never_tier"].iloc[0]) if len(df) else 0

        optimal = df[df["is_optimal_at_g0"]]
        rest = df[~df["is_optimal_at_g0"]]

        ax.scatter(
            rest["asymptotic_variance"],
            rest["_y"],
            c=COLOR_COVARIATE,
            marker="o",
            s=48,
            alpha=0.8,
            edgecolors="none",
            zorder=3,
            label="valid adjustment set",
        )
        ax.scatter(
            optimal["asymptotic_variance"],
            optimal["_y"],
            c=COLOR_OPTIMAL,
            marker="*",
            s=300,
            edgecolors="black",
            linewidths=0.9,
            zorder=5,
            label="O*(G0), optimal set",
        )

        if len(df):
            ax.axhspan(never_tier - 0.5, never_tier + 0.5, color="#f0f0f0", zorder=0)
        ax.axhline(never_tier - 1 + 0.5, color="#cfcfcf", lw=0.9, ls=":", zorder=1)

        yticks = [*list(range(0, max_shell + 1)), never_tier]
        yticklabels = [str(t) for t in range(0, max_shell + 1)] + ["never fails\n(UNREACHED)"]
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)
        ax.set_xlabel("asymptotic variance of the adjusted estimator")
        ax.set_ylabel("r_val (breakdown radius for Z-validity)")
        ax.set_title(f"scenario {label}")

    handles, labels_ = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels_,
        loc="lower center",
        ncol=len(handles),
        frameon=False,
        bbox_to_anchor=(0.5, -0.06),
    )
    fig.suptitle("The robustness frontier: efficiency (variance) vs. robustness (r_val)", y=1.03)
    fig.tight_layout()
    return _save(fig, "fig6_robustness_frontier", out_dir)


def fig6_robustness_frontier_caption(results: dict[str, dict]) -> str:
    """Caption for :func:`fig6_robustness_frontier`."""
    bits = []
    for label in SCENARIO_LABELS:
        df = frontier_plot_frame(results[label])
        optimal = df[df["is_optimal_at_g0"]]
        if len(optimal):
            row = optimal.iloc[0]
            r_val_txt = "UNREACHED" if row["never_fails"] else str(int(row["r_val"]))
            bits.append(
                f"{label}: O*(G0) has r_val={r_val_txt}, variance={row['asymptotic_variance']:.3g}"
            )
    return (
        "Every valid adjustment set of G0 plotted as asymptotic variance (x) against r_val (y), "
        "with a labelled 'never fails' tier standing in for the UNREACHED sentinel rather than "
        "a bare numeric -1. The starred point is O*(G0), the asymptotically most efficient set. "
        "Take from this whether the most efficient estimator is also the least robust one -- "
        "the efficiency-robustness tradeoff, or its absence, per scenario: "
        + ("; ".join(bits) if bits else "no valid adjustment sets found")
        + "."
    )


# --------------------------------------------------------------------------------
# orchestrator
# --------------------------------------------------------------------------------


def build_all_figures(results: dict[str, dict]) -> list[Path]:
    """Build all six figures for every scenario and write them under ``figures/``.

    Args:
        results: Mapping ``"A"``/``"B"``/``"C"`` -> the dict returned by
            :func:`~bkrobust.demo.pipeline.run_scenario`.

    Returns:
        Every path written (both ``.pdf`` and ``.png`` for each figure produced).
        :func:`fig5_baseline_scatter` may contribute zero paths if its dependency
        is unavailable; every other figure always contributes exactly two.
    """
    out_dir = FIGURES_ROOT
    paths: list[Path] = []

    paths += fig1_example(results, out_dir)
    for label in SCENARIO_LABELS:
        paths += fig2_layered_space(results[label], label, out_dir)
    paths += fig3_shell_profile(results, out_dir)
    for label in SCENARIO_LABELS:
        paths += fig4_witness(results[label], label, out_dir)
    paths += fig5_baseline_scatter(results, out_dir)
    paths += fig6_robustness_frontier(results, out_dir)

    return paths
