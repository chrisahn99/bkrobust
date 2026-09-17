#!/usr/bin/env python3
"""Figures for the Panel B explainer and the T3 figure ideation.

Every number is read from ~/bkrobust/results/interval at plot time; nothing is typed in.
The L95 pipeline here is a reimplementation of the one that produced TABLE1_L95.md and it
reproduces all 27 of its Panel A cells exactly (see verify() at the bottom).

Palette: the data-viz reference instance (categorical slots 1-3, sequential blue ramp,
neutral reference gray). Marks thin, grid hairline, direct labels selective.
"""
import csv
import gzip
import json
import math
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

SRC = os.path.expanduser("~/bkrobust/results/interval")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figs")
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- palette
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8a84"
GRID, REF = "#e6e5e2", "#c9c8c4"
SEQ = ["#bcd5f2", "#6ba3e3", "#1f5aa8"]          # sequential blue, light -> dark (b = 0,1,2)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 7.5,
    "axes.edgecolor": REF,
    "axes.linewidth": 0.6,
    "axes.labelcolor": INK2,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "text.color": INK,
    "figure.dpi": 200,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,
})


def despine(ax, keep=("left", "bottom")):
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(s in keep)


# ---------------------------------------------------------------- data
def load_rows():
    return list(csv.DictReader(open(f"{SRC}/rows.csv")))


def column_population(rows):
    """(network, x, y, condition) -> panel, for informative rows where b is reachable."""
    keep = {}
    for r in rows:
        if r["informative"] in ("1", "True", "true") and r["b_realised"] != "":
            keep[(r["network"], r["x"], r["y"], r["condition"])] = r["panel"]
    return keep


PROCS = ["PLAIN", "FAR1", "NEAR1", "SCREEN", "SCREEN2", "I1", "I2", "BLANKET", "STOP1", "STOP2"]


def lam(tau, lo, hi):
    m, h = (lo + hi) / 2.0, (hi - lo) / 2.0
    if h == 0:
        return 0.0 if abs(tau - m) < 1e-12 else math.inf
    return abs(tau - m) / h


def load_replicates(keep, n_want=1000):
    """cell (panel, condition) -> proc -> list of (lambda, half-width, covered)."""
    cells = defaultdict(lambda: defaultdict(list))
    with gzip.open(f"{SRC}/replicates.csv.gz", "rt") as f:
        for r in csv.DictReader(f):
            if int(r["n"]) != n_want:
                continue
            k = (r["network"], r["x"], r["y"], r["condition"])
            if k not in keep:
                continue
            tau = float(r["tau"])
            for p in PROCS:
                lo, hi = float(r[f"{p}_lo"]), float(r[f"{p}_hi"])
                cells[(keep[k], r["condition"])][p].append(
                    (lam(tau, lo, hi), (hi - lo) / 2.0, lo - 1e-12 <= tau <= hi + 1e-12))
    return cells


def lstar(lams):
    s = sorted(lams)
    return s[min(math.ceil(0.95 * len(s)) - 1, len(s) - 1)]


def cell_L95(cell):
    out = {}
    for p, v in cell.items():
        ls = lstar([x[0] for x in v])
        cw = sum(2 * ls * x[1] for x in v) / len(v) if math.isfinite(ls) else math.inf
        out[p] = dict(lstar=ls, calw=cw, cov=sum(x[2] for x in v) / len(v), n=len(v),
                      meanh=sum(x[1] for x in v) / len(v))
    b = out["BLANKET"]["calw"]
    for p in out:
        out[p]["L95"] = out[p]["calw"] / b
    return out


# ================================================================ PART 1
def fig_lambda(cells):
    """What lambda is, on three real replicates of the one-wrong-claim condition."""
    keep = column_population(load_rows())
    picks = []
    with gzip.open(f"{SRC}/replicates.csv.gz", "rt") as f:
        for r in csv.DictReader(f):
            if int(r["n"]) != 1000 or r["condition"] != "CTRL_b1":
                continue
            k = (r["network"], r["x"], r["y"], r["condition"])
            if k not in keep or keep[k] != "A":
                continue
            tau = float(r["tau"])
            lo, hi = float(r["PLAIN_lo"]), float(r["PLAIN_hi"])
            L = lam(tau, lo, hi)
            if math.isfinite(L):
                picks.append((L, r["network"], r["x"], r["y"], tau, lo, hi))
    picks.sort()
    # three real replicates at representative factors, chosen so all three fit one axis
    def nearest(t):
        return min(picks, key=lambda q: abs(q[0] - t))
    chosen = [nearest(0.35), nearest(1.05), nearest(3.20)]

    fig, ax = plt.subplots(figsize=(6.6, 2.15))
    ticks, labels = [], []
    for i, (L, net, x, y, tau, lo, hi) in enumerate(chosen):
        yy = 2 - i
        m, h = (lo + hi) / 2, (hi - lo) / 2
        ax.plot([-1, 1], [yy, yy], color=BLUE, lw=7, solid_capstyle="butt", zorder=2)
        ax.plot([-L, L], [yy, yy], color=ORANGE, lw=3.2, solid_capstyle="butt",
                alpha=.85, zorder=3)
        ax.plot([0], [yy], marker="|", color="white", ms=9, mew=1.1, zorder=4)
        ax.plot([(tau - m) / h], [yy], marker="D", color=INK, ms=5.2, zorder=5,
                mec="white", mew=1.0)
        ax.text(max(L, 1) + .16, yy, f"$\\lambda_i$ = {L:.2f}", va="center", fontsize=7,
                color=ORANGE if L > 1 else INK2)
        ticks.append(yy)
        labels.append(f"{net}\n{x}$\\to${y}")
    ax.axvline(0, color=REF, lw=.6, zorder=0)
    ax.axvline(1, color=REF, lw=.6, ls=(0, (1, 2)), zorder=0)
    ax.axvline(-1, color=REF, lw=.6, ls=(0, (1, 2)), zorder=0)
    ax.set_xlim(-4.6, 4.9)
    ax.set_ylim(-.45, 2.6)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=6.5, linespacing=1.3)
    ax.set_xlabel("distance from the reported interval's midpoint, in half-widths")
    despine(ax, keep=("bottom",))
    ax.tick_params(axis="y", length=0)
    ax.legend(handles=[
        Line2D([], [], color=BLUE, lw=7, label="the interval as reported"),
        Line2D([], [], color=ORANGE, lw=3.2, label="scaled until it reaches the truth"),
        Line2D([], [], color=INK, marker="D", ls="none", ms=5, label="the true effect"),
    ], loc="lower center", bbox_to_anchor=(.5, 1.0), ncol=3, frameon=False, fontsize=6.8,
        handlelength=1.6, columnspacing=1.6)
    fig.savefig(f"{OUT}/f1_lambda.pdf")
    plt.close(fig)
    return chosen


def fig_calibration(cells):
    """The distribution of lambda, and lambda* at the 95th percentile."""
    c = cells[("A", "CTRL_b1")]
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.5), gridspec_kw=dict(wspace=.28))
    for ax, proc, col, name in zip(axes, ["PLAIN", "I1"], [BLUE, ORANGE],
                                   ["report the estimate as committed",
                                    "interval over one revision"]):
        v = [x[0] for x in c[proc] if math.isfinite(x[0])]
        ninf = sum(1 for x in c[proc] if not math.isfinite(x[0]))
        ls = lstar([x[0] for x in c[proc]])
        b = np.logspace(-3, math.log10(max(max(v), 2) * 1.2), 46)
        ax.hist(v, bins=b, color=col, alpha=.55, lw=0)
        ax.axvline(1, color=REF, lw=.8)
        ax.axvline(ls, color=INK, lw=1.0)
        ax.text(ls * 1.12, ax.get_ylim()[1] * .92, f"$\\lambda^*$ = {ls:.2f}",
                fontsize=7, color=INK)
        ax.text(1, -ax.get_ylim()[1] * .10, "1", ha="center", fontsize=6.4, color=MUTED)
        ax.set_xscale("log")
        ax.set_title(name, fontsize=7.2, color=INK, pad=5)
        ax.set_xlabel("$\\lambda_i$, the factor the interval needs")
        ax.grid(axis="y", color=GRID, lw=.5)
        ax.set_axisbelow(True)
        despine(ax)
        sub = f"{len(v):,} of {len(c[proc]):,} replicates finite"
        if ninf:
            sub += f"; {ninf} report a point away from the truth"
        ax.text(0, -.30, sub, transform=ax.transAxes, fontsize=6.3, color=MUTED)
    axes[0].set_ylabel("replicates")
    fig.savefig(f"{OUT}/f2_calibration.pdf")
    plt.close(fig)


def fig_ratio(cells):
    """Calibrated width, and the ratio that becomes L95."""
    c = cell_L95(cells[("A", "CTRL_b1")])
    order = ["PLAIN", "SCREEN", "I1", "I2", "BLANKET"]
    names = ["report as committed", "leave one claim out", "interval, one revision",
             "interval, two revisions", "discard the knowledge"]
    fig, ax = plt.subplots(figsize=(6.6, 2.0))
    ypos = np.arange(len(order))[::-1]
    for yy, p in zip(ypos, order):
        w = c[p]["L95"]
        col = MUTED if p == "BLANKET" else (ORANGE if w > 1 else BLUE)
        ax.barh(yy, min(w, 5.6), height=.52, color=col, alpha=.85, lw=0)
        lab = f"{w:.3f}" + ("  (off scale)" if w > 5.6 else "")
        ax.text(min(w, 5.6) + .07, yy, lab, va="center", fontsize=7,
                color=col if p != "BLANKET" else INK2)
    ax.axvline(1, color=INK, lw=.9)
    ax.text(1, len(order) - .35, "the whole class = 1", fontsize=6.6, color=INK, ha="center")
    ax.set_yticks(ypos)
    ax.set_yticklabels(names, fontsize=7)
    ax.set_xlim(0, 6.6)
    ax.set_xlabel("calibrated width, relative to discarding the knowledge  ($L_{95}$)")
    ax.grid(axis="x", color=GRID, lw=.5)
    ax.set_axisbelow(True)
    despine(ax)
    fig.savefig(f"{OUT}/f3_ratio.pdf")
    plt.close(fig)


def fig_panelB(cells):
    """Panel B as a chart: L95 and raw coverage, procedures as rows, b as a sequential ramp."""
    conds = ["CTRL_b0", "CTRL_b1", "CTRL_b2"]
    labels = ["no wrong claim", "one wrong claim", "two wrong claims"]
    order = ["PLAIN", "SCREEN", "I1", "I2", "BLANKET"]
    names = ["report the estimate as committed", "leave one claim out, re-close",
             "interval over one revision", "interval over two revisions",
             "discard the knowledge, whole class"]
    C = {k: cell_L95(cells[("A", k)]) for k in conds}
    fig, axes = plt.subplots(1, 2, figsize=(6.9, 2.6),
                             gridspec_kw=dict(width_ratios=[1.35, 1], wspace=.06))
    ypos = np.arange(len(order))[::-1]
    off = [.22, 0, -.22]
    ax = axes[0]
    for j, k in enumerate(conds):
        ax.scatter([C[k][p]["L95"] for p in order], ypos + off[j], s=30, color=SEQ[j],
                   zorder=3, lw=.7, edgecolor="white", label=labels[j])
    ax.axvline(1, color=INK, lw=.9, zorder=2)
    ax.set_xscale("log")
    ax.set_xlim(.28, 40)
    ax.set_xticks([.5, 1, 2, 5, 10, 20])
    ax.set_xticklabels([".5", "1", "2", "5", "10", "20"])
    ax.set_yticks(ypos)
    ax.set_yticklabels(names, fontsize=7)
    ax.set_ylim(-.85, len(order) - .4)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("$L_{95}$: calibrated width, class = 1")
    ax.grid(axis="x", color=GRID, lw=.5)
    ax.set_axisbelow(True)
    despine(ax)
    ax.text(1.12, -.52, "wider than throwing the knowledge away $\\rightarrow$",
            fontsize=6.2, color=MUTED, va="center")

    ax = axes[1]
    for j, k in enumerate(conds):
        ax.scatter([C[k][p]["cov"] for p in order], ypos + off[j], s=30, color=SEQ[j],
                   zorder=3, lw=.7, edgecolor="white")
    ax.axvline(.95, color=INK, lw=.9, zorder=2)
    ax.set_yticks(ypos)
    ax.set_yticklabels([])
    ax.set_ylim(-.85, len(order) - .4)
    ax.set_xlim(.42, 1.03)
    ax.set_xlabel("raw coverage of the true effect")
    ax.grid(axis="x", color=GRID, lw=.5)
    ax.set_axisbelow(True)
    despine(ax, keep=("bottom",))
    ax.tick_params(axis="y", length=0)
    ax.text(.95, -.52, "0.95", fontsize=6.4, color=INK, ha="center", va="center")
    axes[0].legend(loc="upper center", bbox_to_anchor=(.78, 1.19), ncol=3, frameon=False,
                   fontsize=6.9, handletextpad=.3, columnspacing=1.3)
    fig.savefig(f"{OUT}/f4_panelB.pdf")
    plt.close(fig)


def fig_why_above_one(cells):
    """Why a displaced midpoint makes L95 exceed 1: the scaling is about the midpoint."""
    fig, ax = plt.subplots(figsize=(6.6, 1.95))
    rows = [("a narrow interval centred right", 0.0, .55, BLUE),
            ("a narrow interval centred wrong", 2.1, .55, ORANGE),
            ("the whole class", 1.05, 2.4, MUTED)]
    for i, (name, m, h, col) in enumerate(rows):
        yy = 2 - i
        ax.plot([m - h, m + h], [yy, yy], color=col, lw=5.5, solid_capstyle="butt", zorder=3)
        ax.plot([m], [yy], marker="|", color="white", ms=8, mew=1.2, zorder=4)
        L = abs(0 - m) / h
        ax.plot([m - L * h, m + L * h], [yy, yy], color=col, lw=5.5, alpha=.22,
                solid_capstyle="butt", zorder=1)
        ax.text(-3.6, yy, name, ha="left", va="center", fontsize=7, color=INK2)
        ax.text(m + L * h + .12, yy, f"needs $\\times${L:.1f}", va="center", fontsize=6.9,
                color=col)
    ax.axvline(0, color=INK, lw=.9)
    ax.text(0, 2.62, "the true effect", ha="center", fontsize=6.8, color=INK)
    ax.set_xlim(-3.7, 6.6)
    ax.set_ylim(-.5, 2.95)
    ax.set_yticks([])
    ax.set_xticks([])
    despine(ax, keep=())
    ax.text(0, -.02, "Scaling happens about the interval's own midpoint, so an interval that is "
                     "narrow and displaced needs a larger factor than one that is wide and honest. "
                     "That is how a procedure earns a number above 1.",
            transform=ax.transAxes, fontsize=6.5, color=MUTED, va="top", wrap=True)
    fig.savefig(f"{OUT}/f5_why_above_one.pdf")
    plt.close(fig)


# ================================================================ PART 3 — T3
def load_T3():
    T = json.load(open(f"{SRC}/table2.json"))
    P = [json.loads(l) for l in open(f"{SRC}/profiles.jsonl")]
    idx = {(p["network"], p["x"], p["y"]): p for p in P if p["condition"] == "A_ORACLE"}
    out = []
    for r in T["rows"]:
        x, y = [s.strip() for s in r["query"].split("->")]
        p = idx[(r["network"], x, y)]
        b = {d["r"]: d for d in p["budgets"]}
        out.append(dict(row=r, prof=p, budgets=b, x=x, y=y))
    return out, T


def short(net):
    return {"Acid_1996": "Acid 1996", "Didelez_2010": "Didelez 2010",
            "Kampen_2014": "Kampen 2014", "Polzer_2012": "Polzer 2012",
            "Sebastiani_2005": "Sebastiani 2005", "Shrier_2008": "Shrier 2008",
            "Thoemmes_2013": "Thoemmes 2013", "mediator": "mediator",
            "paths": "paths"}.get(net, net)


def fig_T3_dumbbell(D):
    """Candidate A. Every named analysis, normalised so the reported effect sits at 1."""
    XLO, XHI = -0.55, 2.15

    def norm(d, v):
        t = d["budgets"][0]["pop"][0]
        sc = 1.0 / t if t else 0.0
        return tuple(sorted((v[0] * sc, v[1] * sc)))

    def span(d):
        a, b = norm(d, d["budgets"][1]["pop"])
        return b - a

    groups = [("published analyses, the query declared by the authors",
               sorted([x for x in D if x["row"]["source"] == "declared"],
                      key=lambda d: -span(d))),
              ("queries on the four networks with fitted coefficients",
               sorted([x for x in D if x["row"]["source"] == "frame"],
                      key=lambda d: -span(d)))]
    heights = [len(g) for _, g in groups]
    fig, axes = plt.subplots(2, 2, figsize=(6.9, .195 * sum(heights) + 1.35), sharex="col",
                             gridspec_kw=dict(height_ratios=heights, hspace=.16,
                                              width_ratios=[3.05, 1.0], wspace=.02))
    for (axp, axl), (title, g) in zip(axes, groups):
        ticks, labels = [], []
        for i, d in enumerate(g):
            yy = len(g) - 1 - i
            c0, c1 = norm(d, d["prof"]["class_pop"])
            r0, r1 = norm(d, d["budgets"][1]["pop"])
            axp.plot([max(c0, XLO), min(c1, XHI)], [yy, yy], color="#dedcd7", lw=7.5,
                     solid_capstyle="butt", zorder=1)
            axp.plot([max(r0, XLO), min(r1, XHI)], [yy, yy], color=ORANGE, lw=4.2,
                     solid_capstyle="butt", zorder=2)
            if min(r0, c0) < XLO:
                axp.plot([XLO], [yy], marker="<", color=ORANGE, ms=3.4, zorder=3)
            if max(r1, c1) > XHI:
                axp.plot([XHI], [yy], marker=">", color=ORANGE, ms=3.4, zorder=3)
            axp.plot([1], [yy], marker="o", color=BLUE, ms=4.4, zorder=4, mec="white", mew=.9)
            ticks.append(yy)
            labels.append(f"{short(d['row']['network'])}   {d['x']}$\\to${d['y']}")
            bc = d["row"].get("breaking_claims")
            txt = (" and ".join(f"{a}$\\to${b}" for a, b in bc) if bc
                   else ("nothing moves it" if not d["row"]["informative"] else ""))
            if txt:
                axl.text(.02, yy, txt, va="center", ha="left", fontsize=5.8,
                         color=INK2 if bc else MUTED)
        axp.axvline(0, color=INK, lw=.9, zorder=5)
        axp.axvline(1, color=REF, lw=.7, zorder=0)
        axp.set_yticks(ticks)
        axp.set_yticklabels(labels, fontsize=6.3)
        axp.set_ylim(-.65, len(g) - .35)
        axp.set_xlim(XLO, XHI)
        axp.grid(axis="x", color=GRID, lw=.5)
        axp.set_axisbelow(True)
        axp.tick_params(axis="y", length=0)
        despine(axp, keep=("bottom",))
        axp.set_title(title, fontsize=6.9, color=INK, style="italic", loc="left", pad=3)
        axl.set_xlim(0, 1)
        axl.set_ylim(-.65, len(g) - .35)
        axl.set_xticks([])
        axl.set_yticks([])
        despine(axl, keep=())
    axes[0][0].spines["bottom"].set_visible(False)
    axes[0][0].tick_params(axis="x", length=0)
    axes[0][1].set_title("first breaking claim", fontsize=6.6, color=INK2,
                         style="italic", loc="left", pad=3)
    axes[1][0].set_xticks([0, .5, 1, 1.5, 2])
    axes[1][0].set_xticklabels(["no effect", "half", "as reported", "", "double"], fontsize=6.5)
    axes[1][0].set_xlabel("the causal effect, in units of the effect the analyst reported")
    axes[0][0].legend(handles=[
        Line2D([], [], color=BLUE, marker="o", ls="none", ms=4.4, label="reported estimate"),
        Line2D([], [], color=ORANGE, lw=4.2, label="after retracting one claim"),
        Line2D([], [], color="#dedcd7", lw=7.5, label="with no knowledge at all"),
    ], loc="lower center", bbox_to_anchor=(.5, 1.12), ncol=3, frameon=False, fontsize=6.8,
        handlelength=1.5, columnspacing=1.8)
    fig.savefig(f"{OUT}/f6_T3_dumbbell.pdf")
    plt.close(fig)
def fig_T3_steps(D):
    """Candidate B. Interval width against retraction budget, normalised to the class width."""
    fig, ax = plt.subplots(figsize=(3.3, 2.5))
    n_step1 = 0
    n_inf = sum(1 for d in D if d["row"]["informative"])
    for d in D:
        if not d["row"]["informative"]:
            continue
        cw = d["prof"]["class_pop"][1] - d["prof"]["class_pop"][0]
        if cw <= 0:
            continue
        ws = [(d["budgets"][r]["pop"][1] - d["budgets"][r]["pop"][0]) / cw for r in range(4)]
        col = ORANGE if ws[1] > .9 else BLUE
        n_step1 += ws[1] > .9
        ax.step(range(4), ws, where="post", color=col, lw=1.1, alpha=.5)
    ax.set_xticks(range(4))
    ax.set_xlabel("claims retracted")
    ax.set_ylabel("interval width, class = 1")
    ax.set_ylim(-.04, 1.08)
    ax.grid(axis="y", color=GRID, lw=.5)
    ax.set_axisbelow(True)
    despine(ax)
    ax.text(.04, .93, f"{n_step1} of {n_inf} reach the full\nclass at one retraction",
            transform=ax.transAxes, fontsize=6.6, color=ORANGE, va="top", linespacing=1.3)
    fig.savefig(f"{OUT}/f7_T3_steps.pdf")
    plt.close(fig)


def fig_T3_strip(D):
    """Candidate C. The state of the conclusion at each budget, from the stored null radius."""
    decl = [d for d in D if d["row"]["source"] == "declared"]
    frame = [d for d in D if d["row"]["source"] == "frame"]

    def r0_of(d):
        v = str(d["row"]["r0_pop"])
        if v == "never" or v.startswith(">"):
            return math.inf          # no retraction within the search depth reaches zero
        return int(v)

    rows = sorted(decl, key=r0_of) + sorted(frame, key=r0_of)
    fig, ax = plt.subplots(figsize=(3.35, .168 * len(rows) + 1.0))
    HOLD, GONE = "#c6dcf5", "#eb6834"
    for i, d in enumerate(rows):
        yy = len(rows) - 1 - i
        for r in range(4):
            gone = r >= r0_of(d)
            ax.add_patch(Rectangle((r + .06, yy - .38), .88, .76,
                                   facecolor=GONE if gone else HOLD, lw=0))
    ax.set_xlim(0, 4)
    ax.set_ylim(-.6, len(rows) - .4)
    ax.set_xticks([r + .5 for r in range(4)])
    ax.set_xticklabels([str(r) for r in range(4)], fontsize=6.6)
    ax.set_xlabel("claims retracted")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{short(d['row']['network'])}  {d['x']}$\\to${d['y']}"
                        for d in rows][::-1], fontsize=5.9)
    ax.tick_params(length=0)
    ax.axhline(len(frame) - .5, color=INK, lw=.7)
    despine(ax, keep=())
    ax.legend(handles=[
        Line2D([], [], color=HOLD, lw=7, label="sign determined"),
        Line2D([], [], color=GONE, lw=7, label="interval contains zero"),
    ], loc="lower center", bbox_to_anchor=(.5, 1.0), ncol=2, frameon=False, fontsize=6.6,
        handlelength=1.3)
    fig.savefig(f"{OUT}/f8_T3_strip.pdf")
    plt.close(fig)
def fig_T3_card(D):
    """Candidate D. One analysis, end to end: the teaser."""
    d = next(x for x in D if x["row"]["network"] == "magic-niab" and x["x"] == "G311")
    fig, ax = plt.subplots(figsize=(6.9, 1.75))
    t = d["budgets"][0]["pop"][0]
    sc = 1.0 / t

    def nz(v):
        return tuple(sorted((v[0] * sc, v[1] * sc)))

    lanes = [("what the analyst reports", nz(d["budgets"][0]["pop"]), BLUE,
              "the estimate, on the graph their own claims produced"),
             ("if one claim is withdrawn", nz(d["budgets"][1]["pop"]), ORANGE,
              "retract  G1853 $\\to$ G311"),
             ("if no knowledge is used at all", nz(d["prof"]["class_pop"]), "#9d9c97",
              "the interval over every graph the data allows")]
    ticks, labels = [], []
    for i, (name, (a, b), col, note) in enumerate(lanes):
        yy = 2 - i
        ax.plot([a, b], [yy, yy], color=col, lw=8, solid_capstyle="butt", zorder=3)
        if b - a < .02:
            ax.plot([a], [yy], marker="o", color=col, ms=7, zorder=4)
        ticks.append(yy)
        labels.append(name)
        ax.text(1.32, yy, note, ha="left", va="center", fontsize=6.5, color=MUTED)
    ax.axvline(0, color=INK, lw=1.0, zorder=2)
    ax.axvline(1, color=REF, lw=.8, zorder=1)
    ax.set_xlim(-.30, 1.30)
    ax.set_ylim(-.55, 2.55)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=7.4)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["no effect", "the reported effect"], fontsize=7)
    ax.tick_params(length=0)
    despine(ax, keep=())
    fig.savefig(f"{OUT}/f9_T3_card.pdf")
    plt.close(fig)
def fig_T3_recommended(D):
    """Candidate E. The dumbbell, with the population it comes from beside it."""
    import csv as _csv
    led = os.path.expanduser("~/bkrobust/results/ledger/rows.csv")
    buckets = defaultdict(int)
    tot = 0
    for r in _csv.DictReader(open(led)):
        if r["status"] != "ok" or not r["arm"].startswith("D_") or r["arm"] == "D_DEGEN":
            continue
        tot += 1
        st = r["r_claim_status"]
        if st == "exact":
            v = int(r["r_claim"])
            buckets[str(v) if v <= 3 else "beyond 3"] += 1
        elif st in ("unreached", "no_retractable_edges"):
            buckets["never"] += 1
        else:
            buckets["beyond 3"] += 1
    order = ["1", "2", "3", "beyond 3", "never"]
    fig, ax = plt.subplots(figsize=(2.6, 2.0))
    vals = [buckets[k] / tot for k in order]
    cols = [ORANGE] + [MUTED] * 4
    ax.bar(range(len(order)), vals, color=cols, width=.62, lw=0)
    for i, v in enumerate(vals):
        ax.text(i, v + .012, f"{v:.0%}", ha="center", fontsize=6.6,
                color=ORANGE if i == 0 else INK2)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, fontsize=6.4)
    ax.set_ylabel("share of queries")
    ax.set_ylim(0, max(vals) * 1.22)
    ax.set_xlabel("claims whose retraction breaks the set")
    ax.grid(axis="y", color=GRID, lw=.5)
    ax.set_axisbelow(True)
    despine(ax)
    ax.text(0, 1.06, f"all {tot:,} certifiable queries", transform=ax.transAxes,
            fontsize=6.6, color=INK)
    fig.savefig(f"{OUT}/f10_T3_population.pdf")
    plt.close(fig)
    return tot, buckets



# ---------------------------------------------------------------- the full pool
PANEL_B_NETS = ("ecoli70", "arth150", "magic-niab", "magic-irri")
CENSORED_R0 = 4


def load_pool():
    """Every query on the four fitted-coefficient networks that T3 was eligible to show,
    with a flag saying whether the cap of five actually let it in."""
    P = [json.loads(l) for l in open(f"{SRC}/profiles.jsonl")]
    A = [q for q in P if q["condition"] == "A_ORACLE"]
    pool = [q for q in A if q["network"] in PANEL_B_NETS and q["source"] == "frame"
            and q["z"] is not None]
    taken = set()
    for net in PANEL_B_NETS:
        c = [q for q in pool if q["network"] == net]
        c.sort(key=lambda q: (q["r0_pop"] == CENSORED_R0, q["frame_index"]))
        for q in c[:5]:
            taken.add((q["network"], q["x"], q["y"]))
    out = []
    for q in pool:
        out.append(dict(prof=q, x=q["x"], y=q["y"],
                        printed=(q["network"], q["x"], q["y"]) in taken,
                        budgets={d["r"]: d for d in q["budgets"]},
                        row=dict(network=q["network"], informative=q["informative"],
                                 source="frame", breaking_claims=q["breaking_claims_seed0"])))
    return out


def fig_pool_dumbbell(pool):
    """The same picture as F-A, on every eligible query rather than the printed seventeen."""
    XLO, XHI = -0.55, 2.15

    def norm(d, v):
        t = d["budgets"][0]["pop"][0]
        sc = 1.0 / t if t else 0.0
        return tuple(sorted((v[0] * sc, v[1] * sc)))

    rows = sorted(pool, key=lambda d: (d["prof"]["network"],
                                       -(norm(d, d["budgets"][1]["pop"])[1]
                                         - norm(d, d["budgets"][1]["pop"])[0])))
    fig, ax = plt.subplots(figsize=(6.9, .168 * len(rows) + 1.15))
    ticks, labels = [], []
    for i, d in enumerate(rows):
        yy = len(rows) - 1 - i
        c0, c1 = norm(d, d["prof"]["class_pop"])
        r0, r1 = norm(d, d["budgets"][1]["pop"])
        ax.plot([max(c0, XLO), min(c1, XHI)], [yy, yy], color="#dedcd7", lw=6.4,
                solid_capstyle="butt", zorder=1)
        ax.plot([max(r0, XLO), min(r1, XHI)], [yy, yy], color=ORANGE, lw=3.6,
                solid_capstyle="butt", zorder=2)
        if min(r0, c0) < XLO:
            ax.plot([XLO], [yy], marker="<", color=ORANGE, ms=3.0, zorder=3)
        if max(r1, c1) > XHI:
            ax.plot([XHI], [yy], marker=">", color=ORANGE, ms=3.0, zorder=3)
        ax.plot([1], [yy], marker="o", color=BLUE, ms=3.8, zorder=4, mec="white", mew=.8)
        if d["printed"]:
            ax.plot([XLO + .03], [yy], marker="s", color=INK, ms=2.6, zorder=6)
        ticks.append(yy)
        labels.append(f"{d['prof']['network']}  {d['x']}$\\to${d['y']}")
    ax.axvline(0, color=INK, lw=.9, zorder=5)
    ax.axvline(1, color=REF, lw=.7, zorder=0)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=5.6)
    ax.set_ylim(-.65, len(rows) - .35)
    ax.set_xlim(XLO, XHI)
    ax.set_xticks([0, .5, 1, 1.5, 2])
    ax.set_xticklabels(["no effect", "half", "as reported", "", "double"], fontsize=6.5)
    ax.set_xlabel("the causal effect, in units of the effect the analyst reported")
    ax.grid(axis="x", color=GRID, lw=.5)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
    despine(ax, keep=("bottom",))
    ax.legend(handles=[
        Line2D([], [], color=BLUE, marker="o", ls="none", ms=3.8, label="reported estimate"),
        Line2D([], [], color=ORANGE, lw=3.6, label="after retracting one claim"),
        Line2D([], [], color="#dedcd7", lw=6.4, label="with no knowledge at all"),
        Line2D([], [], color=INK, marker="s", ls="none", ms=2.6, label="printed in T3"),
    ], loc="lower center", bbox_to_anchor=(.5, 1.0), ncol=4, frameon=False, fontsize=6.6,
        handlelength=1.4, columnspacing=1.4)
    fig.savefig(f"{OUT}/f11_pool_dumbbell.pdf")
    plt.close(fig)


# ================================================================ verify
def verify(cells):
    pub = {"CTRL_b0": dict(PLAIN=.464, FAR1=.655, SCREEN=.932, I1=.957, I2=.984,
                           BLANKET=1.0, STOP1=.956, STOP2=.983, NEAR1=.941),
           "CTRL_b1": dict(PLAIN=5.123, FAR1=10.771, SCREEN=.691, I1=.851, I2=.980,
                           BLANKET=1.0, STOP1=.715, STOP2=.980, NEAR1=.784),
           "CTRL_b2": dict(PLAIN=4.082, FAR1=7.539, SCREEN=5.967, I1=4.234, I2=.847,
                           BLANKET=1.0, STOP1=4.640, STOP2=.791, NEAR1=4.087)}
    bad = 0
    for cond, exp in pub.items():
        got = cell_L95(cells[("A", cond)])
        for p, v in exp.items():
            if abs(got[p]["L95"] - v) >= 0.0015:
                print(f"  MISMATCH {cond} {p}: {got[p]['L95']:.3f} vs {v}")
                bad += 1
    print(f"verify: {27 - bad}/27 Panel A cells reproduce TABLE1_L95.md")
    return bad == 0


if __name__ == "__main__":
    rows = load_rows()
    keep = column_population(rows)
    cells = load_replicates(keep)
    assert verify(cells)
    chosen = fig_lambda(cells)
    fig_calibration(cells)
    fig_ratio(cells)
    fig_panelB(cells)
    fig_why_above_one(cells)
    D, T = load_T3()
    fig_T3_dumbbell(D)
    fig_T3_steps(D)
    fig_T3_strip(D)
    fig_T3_card(D)
    tot, buckets = fig_T3_recommended(D)
    pool = load_pool()
    fig_pool_dumbbell(pool)
    # numbers the report quotes, printed so they can be checked against the text
    c1 = cell_L95(cells[("A", "CTRL_b1")])
    print(json.dumps({
        "lambda_examples": [dict(net=c[1], q=f"{c[2]}->{c[3]}", lam=round(c[0], 3)) for c in chosen],
        "b1_lstar_PLAIN": round(c1["PLAIN"]["lstar"], 3),
        "b1_lstar_I1": round(c1["I1"]["lstar"], 3),
        "b1_lstar_BLANKET": round(c1["BLANKET"]["lstar"], 3),
        "b1_meanhalfwidth_PLAIN": round(c1["PLAIN"]["meanh"], 4),
        "b1_meanhalfwidth_I1": round(c1["I1"]["meanh"], 4),
        "b1_meanhalfwidth_BLANKET": round(c1["BLANKET"]["meanh"], 4),
        "b1_calw_PLAIN": round(c1["PLAIN"]["calw"], 4),
        "b1_calw_I1": round(c1["I1"]["calw"], 4),
        "b1_calw_BLANKET": round(c1["BLANKET"]["calw"], 4),
        "b1_replicates": c1["PLAIN"]["n"],
        "T3_rows": len(D),
        "T3_class_at_r1_of_informative": sum(
            1 for d in D
            if d["row"]["informative"]
            and (d["prof"]["class_pop"][1] - d["prof"]["class_pop"][0]) > 0
            and (d["budgets"][1]["pop"][1] - d["budgets"][1]["pop"][0])
            / (d["prof"]["class_pop"][1] - d["prof"]["class_pop"][0]) > .9),
        "T3_informative": sum(1 for d in D if d["row"]["informative"]),
        "ledger_total": tot,
        "ledger_buckets": {k: round(v / tot, 4) for k, v in buckets.items()},
        "pool_size": len(pool),
        "pool_printed": sum(1 for d in pool if d["printed"]),
        "pool_informative": sum(1 for d in pool if d["row"]["informative"]),
        "pool_finite_r0": sum(1 for d in pool if d["prof"]["r0_pop"] != CENSORED_R0),
        "printed_finite_r0": sum(1 for d in pool
                                 if d["printed"] and d["prof"]["r0_pop"] != CENSORED_R0),
        "dropped_finite_r0": sum(1 for d in pool
                                 if not d["printed"] and d["prof"]["r0_pop"] != CENSORED_R0),
    }, indent=1))
