#!/usr/bin/env python3
"""Computed figures for chapters 7 to 11 (the protocol run of 11 to 13 September).

English port of the figure script of the 13 September document on the protocol run and the
Rosenbaum reframing. It reads ONLY result files under ~/bkrobust/results/ (no number typed by
hand), plots the same numbers, and writes seven vector PDFs into figs/ plus
figs/numbers_gen_figs.json, which keeps every plotted number with the file and the key it came
from. Only the words and the decimal separator changed.

    /Users/josecosta/bkrobust/.venv/bin/python gen_figs.py [--preview DIR]

--preview DIR also writes PNG proofs into DIR (outside this directory).

Style: 16 cm wide (the A4 text width with 2.4 cm margins), at most 9 cm high, fonts of at least
8 pt at print size, sans serif, Okabe-Ito colour roles (R grDevices hexes) with redundant marker,
dash or hatch, direct labels on the series and no keyed legend.

Outputs (figs/): buckets_7b.pdf, supplier_scale.pdf, cascade.pdf, pest_dangerous.pdf,
l95_crossing.pdf, r0_distribution.pdf, stop_refusal.pdf.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.transforms import blended_transform_factory  # noqa: E402

HERE = Path(__file__).resolve().parent
FIGS = HERE / "figs"
HOME = Path.home()
RES = HOME / "bkrobust" / "results"

T4_PATH = RES / "ledger" / "table4.json"
AUDIT_PATH = RES / "elicit" / "knowledge_audit.json"
CASCADE_PATH = RES / "ledger" / "cascade_test.json"
PEST_ROWS = RES / "pest" / "rows.csv"
LEDGER_ROWS = RES / "ledger" / "rows.csv"
T1_PATH = RES / "interval" / "table1.json"
INT_ROWS = RES / "interval" / "rows.csv"

CM = 1 / 2.54
W = 16 * CM

# Okabe & Ito, canonical values of R's grDevices
OI = dict(
    black="#000000", orange="#E69F00", skyblue="#56B4E9", green="#009E73",
    yellow="#F0E442", blue="#0072B2", vermilion="#D55E00", purple="#CC79A7",
)
INK = "#1A1A1A"
INK2 = "#4A4E54"   # dark grey of mugango-apostila.sty
MUTED = "#8A9097"
GRID = "#E3E4E6"
GREY = "#9AA0A6"   # medium grey of the .sty
FS = 8.0           # font floor
FS_LAB = 8.5

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica Neue", "DejaVu Sans"],
    "font.size": FS,
    "axes.labelsize": FS_LAB,
    "axes.titlesize": FS_LAB,
    "xtick.labelsize": FS,
    "ytick.labelsize": FS_LAB,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.edgecolor": MUTED,
    "axes.linewidth": 0.6,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 0,
    "text.color": INK,
    "axes.labelcolor": INK,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "hatch.linewidth": 0.8,
    "savefig.dpi": 220,
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "mathtext.sf": "Arial",
    "mathtext.tt": "Arial",
    "mathtext.cal": "Arial",
    "mathtext.bfit": "Arial:bold:italic",
})

# ------------------------------------------------------------------ number register
NUM: dict[str, list[dict]] = {}


def src(p: Path) -> str:
    return "~/" + str(p.relative_to(HOME))


def rec(fig: str, label: str, value, source: Path, key: str):
    if isinstance(value, float) and math.isinf(value):
        value_j = "inf"
    else:
        value_j = value
    NUM.setdefault(fig, []).append(
        {"label": label, "value": value_j, "source": src(source), "key": key}
    )
    return value


def fnum(x: float, nd: int = 3, sign: bool = False) -> str:
    """Number with a decimal point and a typographic minus sign."""
    s = f"{x:+.{nd}f}" if (sign and round(x, nd) != 0) else f"{abs(x) if round(x, nd) == 0 else x:.{nd}f}"
    return s.replace("-", "−")


def fpct(p: float) -> str:
    if p == 0:
        return "0"
    v = 100 * p
    return f"{v:.0f}%" if v >= 9.95 else f"{v:.1f}%"


def as_float(v) -> float:
    return float("inf") if v in ("inf", "Infinity") else float(v)


def load_json(p: Path):
    with p.open() as fh:
        return json.load(fh)


def read_csv(p: Path) -> list[dict]:
    with p.open(newline="") as fh:
        return list(csv.DictReader(fh))


def text_width_data(ax, s: str, fontsize: float) -> float:
    """Width of a text in data units of the x axis (linear axis, limits already set)."""
    r = ax.figure.canvas.get_renderer()
    t = ax.text(0, 0, s, fontsize=fontsize)
    wpx = t.get_window_extent(r).width
    t.remove()
    inv = ax.transData.inverted()
    x0 = inv.transform((ax.transData.transform((0, 0))[0], 0))[0]
    x1 = inv.transform((ax.transData.transform((0, 0))[0] + wpx, 0))[0]
    return abs(x1 - x0)


def lum(hexc: str) -> float:
    h = hexc.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))
    f = lambda c: c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4  # noqa: E731
    return 0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b)


def ink_on(hexc: str) -> str:
    return "white" if lum(hexc) < 0.30 else INK


def save(fig, name: str, preview: Path | None):
    FIGS.mkdir(exist_ok=True)
    fig.savefig(FIGS / f"{name}.pdf")
    if preview:
        preview.mkdir(parents=True, exist_ok=True)
        fig.savefig(preview / f"{name}.png", dpi=220)
    plt.close(fig)


# readable supplier names
NAME = {
    "D_LLM": "Qwen 7B",
    "D_LLM_INSTR": "Qwen 7B, second order",
    "D_LLM_SRC": "Llama 8B",
    "D_SCRAMBLED": "Qwen 7B, scrambled names",
    "D_LLM_32B": "Qwen 32B",
    "D_LLM_72B": "Qwen 72B",
    "D_LLM_72B_INSTR": "Qwen 72B, second order",
    "D_SCRAMBLED_72B": "Qwen 72B, scrambled names",
    "D_LLM_SRC_70B": "Llama 70B",
    "D_LLM_GEMMA_27B": "Gemma 27B",
    "D_RAND": "random",
    "D_DEGEN": "no knowledge",
    "A_TRUE": "true claims",
    "A_ORACLE": "oracle",
}


def family(arm: str) -> tuple[str, str, bool]:
    """(colour, marker, filled) by family; presentation variants are drawn hollow."""
    variant = arm.endswith("_INSTR") or "SCRAMBLED" in arm
    if arm.startswith("A_"):
        return INK, "o", True
    if arm == "D_RAND":
        return GREY, "D", False
    if "GEMMA" in arm:
        return OI["purple"], "^", True
    if "SRC" in arm:
        return OI["orange"], "s", True
    return OI["blue"], "o", not variant


def value_label(ax, x, y, txt, dx=5, **kw):
    kw.setdefault("fontsize", FS)
    kw.setdefault("color", INK2)
    return ax.annotate(txt, xy=(x, y), xytext=(dx, 0), textcoords="offset points",
                       ha="left" if dx >= 0 else "right", va="center", **kw)


def header(fig, ax, y, txt, x=0.012, bold=True):
    tr = blended_transform_factory(fig.transFigure, ax.transData)
    ax.text(x, y, txt, transform=tr, ha="left", va="center", fontsize=FS, color=MUTED,
            style="italic", fontweight="bold" if bold else None)


# =============================================================================================
# Figure 1: the four buckets of the primary supplier
# =============================================================================================
def fig_buckets(preview):
    F = "buckets_7b"
    t4 = load_json(T4_PATH)
    arm = t4["arms"]["D_LLM"]
    n_cert = rec(F, "certifiable rows, D_LLM", arm["certifiable"], T4_PATH,
                 "arms.D_LLM.certifiable")
    n_inv = rec(F, "rows invalid at the truth, D_LLM", arm["truth_invalid_rows"], T4_PATH,
                "arms.D_LLM.truth_invalid_rows")
    rungs = [("B2_free_rule", "free rule"), ("B6_hop", "hop radius"),
             ("B7_claim", "claim radius")]
    buckets = [("dangerous", "dangerous", OI["vermilion"], "////"),
               ("safe", "safe", OI["orange"], None),
               ("held", "held", OI["blue"], None), ("slack", "slack", OI["skyblue"], None)]

    fig = plt.figure(figsize=(W, 6.0 * CM))
    ax = fig.add_axes([0.19, 0.2, 0.60, 0.70])
    ax.set_xlim(0, n_cert)
    ax.set_ylim(-0.42, 2.78)
    H = 0.46
    ys = {d: 2 - i for i, (d, _) in enumerate(rungs)}
    tr = blended_transform_factory(ax.transAxes, ax.transData)
    for d, name in rungs:
        b = arm["ladder"][d]["buckets"]
        cnt = {k: rec(F, f"{name}: {en}", int(b.get(k, 0)), T4_PATH,
                      f"arms.D_LLM.ladder.{d}.buckets.{k}" + ("" if k in b else " (missing key = 0)"))
               for k, en, _, _ in buckets}
        assert sum(cnt.values()) == n_cert, (d, cnt)
        assert cnt["dangerous"] + cnt["safe"] == n_inv, (d, cnt)
        y = ys[d]
        x = 0
        for k, en, col, hatch in buckets:
            c = cnt[k]
            if c > 0:
                ax.barh(y, c, left=x, height=H, color=col, edgecolor="white", linewidth=1.0,
                        hatch=hatch, zorder=2)
            full, short = f"{en} {c}", f"{c}"
            if c > 0 and text_width_data(ax, full, FS_LAB) + 6 < c:
                ax.text(x + c / 2, y, full, ha="center", va="center", fontsize=FS_LAB,
                        color=ink_on(col), zorder=3)
            elif c > 0 and text_width_data(ax, short, FS_LAB) + 3 < c:
                ax.text(x + c / 2, y, short, ha="center", va="center", fontsize=FS_LAB,
                        color=ink_on(col), zorder=3)
            else:
                # empty or narrow bucket: explicit mark and label above the bar
                if c == 0:
                    ax.plot([x, x], [y - H / 2, y + H / 2], color=col, lw=2.4, zorder=4,
                            solid_capstyle="butt")
                ax.text(x + max(c, 0) + 1.2, y + H / 2 + 0.04, full, ha="left", va="bottom",
                        fontsize=FS_LAB, color=INK, zorder=3,
                        fontweight="bold" if k == "dangerous" else None)
            x += c
        det = rec(F, f"{name}: detection", arm["ladder"][d]["detection"], T4_PATH,
                  f"arms.D_LLM.ladder.{d}.detection")
        fa = rec(F, f"{name}: false alarm", arm["ladder"][d]["false_alarm"], T4_PATH,
                 f"arms.D_LLM.ladder.{d}.false_alarm")
        ax.text(1.03, y + 0.11, f"detection {fnum(det, 2)}", transform=tr, ha="left",
                va="center", fontsize=FS, color=INK2)
        ax.text(1.03, y - 0.13, f"false alarm {fnum(fa, 2)}", transform=tr, ha="left",
                va="center", fontsize=FS, color=INK2)

    ax.plot([n_inv, n_inv], [-0.36, 2.62], color=INK, lw=1.1, zorder=5, solid_capstyle="butt")
    ax.text(n_inv - 2.5, 2.66, f"invalid at the truth: {n_inv}", ha="right", va="bottom",
            fontsize=FS_LAB, color=INK, fontweight="bold")
    ax.text(n_inv + 2.5, 2.66, f"valid at the truth: {n_cert - n_inv}", ha="left", va="bottom",
            fontsize=FS_LAB, color=INK, fontweight="bold")
    ax.set_yticks([ys[d] for d, _ in rungs])
    ax.set_yticklabels([n for _, n in rungs])
    ax.tick_params(axis="y", length=0, pad=6)
    ax.spines["left"].set_visible(False)
    ax.set_xticks([0, 50, 100, 150, 200, n_cert])
    ax.set_xlabel(f"certifiable rows of the primary supplier, Qwen 7B ({n_cert})")
    save(fig, F, preview)


# =============================================================================================
# Figure 2: supplier scale
# =============================================================================================
def fig_scale(preview):
    F = "supplier_scale"
    t4 = load_json(T4_PATH)
    au = load_json(AUDIT_PATH)
    local = ["D_LLM", "D_LLM_INSTR", "D_LLM_SRC", "D_SCRAMBLED"]
    scale = ["D_LLM_32B", "D_LLM_72B", "D_LLM_72B_INSTR", "D_SCRAMBLED_72B", "D_LLM_SRC_70B",
             "D_LLM_GEMMA_27B"]
    items = ([("h", "local suppliers, 7B and 8B")] + [("r", a) for a in local]
             + [("h", "scale axis, 27B to 72B (Jean Zay)")] + [("r", a) for a in scale]
             + [("h", "reference")] + [("r", "D_RAND")])
    ypos, y = [], 0.0
    for kind, _ in items:
        y -= 0.9 if kind == "h" else 1.0
        ypos.append(y)
    fig = plt.figure(figsize=(W, 9.0 * CM))
    left, gap, right = 0.285, 0.03, 0.02
    wp = (1 - left - right - 2 * gap) / 3
    axs = [fig.add_axes([left + i * (wp + gap), 0.235, wp, 0.75]) for i in range(3)]
    ylo, yhi = min(ypos) - 0.55, max(ypos) + 0.45
    asked_all = {au["per_condition"][a]["asked"] for a in local + scale}
    assert len(asked_all) == 1
    asked = rec(F, "questions per supplier", asked_all.pop(), AUDIT_PATH,
                "per_condition.*.asked (same for all)")

    rows = [(a, yy) for (kind, a), yy in zip(items, ypos) if kind == "r"]
    for ax in axs:
        ax.set_ylim(ylo, yhi)
        ax.spines["left"].set_visible(False)
        for _, yy in rows:
            ax.axhline(yy, color=GRID, lw=0.6, zorder=0)
        ax.set_yticks([])
    axs[0].set_yticks([yy for _, yy in rows])
    axs[0].set_yticklabels([NAME[a] for a, _ in rows], fontsize=FS)
    axs[0].tick_params(axis="y", length=0, pad=4)
    for (kind, txt), yy in zip(items, ypos):
        if kind == "h":
            header(fig, axs[0], yy, txt)

    # (a) claims made
    ax = axs[0]
    ax.set_xlim(0, asked * 1.02)
    ax.axvline(asked, color=MUTED, lw=0.7, zorder=0)
    ax.text(asked - 5, ypos[0], f"of {asked} questions", ha="right", va="center", fontsize=FS,
            color=MUTED)
    for a, yy in rows:
        col, mk, filled = family(a)
        if a == "D_RAND":
            ax.text(8, yy, "random draw on\nQwen 7B's edges", ha="left", va="center",
                    fontsize=FS, color=INK2, linespacing=1.0)
            continue
        v = rec(F, f"{NAME[a]}: claims kept", au["per_condition"][a]["asserted"], AUDIT_PATH,
                f"per_condition.{a}.asserted")
        assert v == t4["elicitation"][a]["asserted"], a
        ax.plot([0, v], [yy, yy], color=col, lw=1.2, alpha=0.35, zorder=1, solid_capstyle="butt")
        ax.plot([v], [yy], marker=mk, ms=5.5, color=col, mfc=col if filled else "white", mew=1.3,
                zorder=3)
        value_label(ax, v, yy, str(v), dx=6)
    ax.set_xticks([0, 100, 200, 300])
    ax.set_xlabel("claims kept after the consistency repair")

    highlight = {"D_LLM", "D_LLM_72B", "D_RAND"}

    def ci_panel(ax, getter, xlabel):
        ax.set_xlim(0, 0.62)
        for a, yy in rows:
            col, mk, filled = family(a)
            (m, lo, hi), source, key = getter(a)
            xl1 = xlabel.replace("\n", " ")
            rec(F, f"{NAME[a]}: {xl1} (estimate)", m, source, key + "[0]")
            rec(F, f"{NAME[a]}: {xl1} (CI low)", lo, source, key + "[1]")
            rec(F, f"{NAME[a]}: {xl1} (CI high)", hi, source, key + "[2]")
            ax.plot([lo, hi], [yy, yy], color=col, lw=1.3, zorder=2, solid_capstyle="butt")
            ax.plot([m], [yy], marker=mk, ms=5.5, color=col, mfc=col if filled else "white",
                    mew=1.3, zorder=3)
            if a in highlight:
                txt = fnum(m)
                if hi + text_width_data(ax, txt, FS) * 1.15 > ax.get_xlim()[1]:
                    # no room right of the interval inside the panel: label left of it
                    value_label(ax, lo, yy, txt, dx=-4)
                else:
                    value_label(ax, hi, yy, txt, dx=4)
        ax.set_xticks([0, 0.2, 0.4, 0.6])
        ax.set_xticklabels([fnum(v, 1) for v in (0, 0.2, 0.4, 0.6)])
        ax.set_xlabel(xlabel, linespacing=1.15)

    def false_frac(a):
        if a == "D_RAND":
            return t4["arms"][a]["false_claim_fraction"], T4_PATH, f"arms.{a}.false_claim_fraction"
        return (au["per_condition"][a]["commission_network_mean"], AUDIT_PATH,
                f"per_condition.{a}.commission_network_mean")

    def invalid(a):
        return t4["arms"][a]["z_invalid_at_truth"], T4_PATH, f"arms.{a}.z_invalid_at_truth"

    ci_panel(axs[1], false_frac, "false fraction\nof the claims")
    ci_panel(axs[2], invalid, "committed set\ninvalid at the truth")
    axs[0].set_xlabel("claims kept after\nthe consistency repair")
    fig.text(0.012, 0.012, "Intervals: 95% CI, bootstrap over networks. False fraction: per-network "
             "mean over what each model asserted;\nfor random, the same fraction measured on the "
             "certifiable rows.", ha="left", va="bottom", fontsize=FS, color=INK2, linespacing=1.15)
    save(fig, F, preview)


# =============================================================================================
# Figure 3: the cascade, edges oriented per claim
# =============================================================================================
def fig_cascade(preview):
    F = "cascade"
    ct = load_json(CASCADE_PATH)
    per = ct["per_arm"]
    order = ["A_ORACLE", "A_TRUE", "D_LLM", "D_LLM_INSTR", "D_LLM_SRC", "D_SCRAMBLED", "D_LLM_32B",
             "D_LLM_72B", "D_LLM_72B_INSTR", "D_SCRAMBLED_72B", "D_LLM_SRC_70B", "D_LLM_GEMMA_27B",
             "D_RAND"]
    assert set(order) == set(per), set(per) ^ set(order)
    fig = plt.figure(figsize=(W, 9.0 * CM))
    axL = fig.add_axes([0.275, 0.155, 0.30, 0.80])
    axR = fig.add_axes([0.63, 0.155, 0.35, 0.80])

    ys, y = {}, 0.0
    for a in order:
        y -= 1.0
        if a == "A_ORACLE":
            y -= 0.35
        if a == "D_LLM":
            y -= 0.9
        ys[a] = y
    axL.set_ylim(min(ys.values()) - 0.6, 0.35)
    axL.set_xlim(0, 3.7)
    axL.spines["left"].set_visible(False)
    for a in order:
        axL.axhline(ys[a], color=GRID, lw=0.6, zorder=0)
    header(fig, axL, ys["A_ORACLE"] + 0.95, "controls that read the truth")
    header(fig, axL, ys["D_LLM"] + 0.95, "suppliers")
    ticklabels = []
    for a in order:
        s = per[a]
        mean = rec(F, f"{NAME[a]}: per-network mean", s["mean"], CASCADE_PATH, f"per_arm.{a}.mean")
        med = rec(F, f"{NAME[a]}: per-network median", s["median"], CASCADE_PATH,
                  f"per_arm.{a}.median")
        nets = rec(F, f"{NAME[a]}: networks", s["networks"], CASCADE_PATH, f"per_arm.{a}.networks")
        col, _, _ = family(a)
        yy = ys[a]
        axL.plot([med, mean], [yy, yy], color=col, lw=1.3, alpha=0.45, zorder=1)
        axL.plot([med], [yy], marker="D", ms=4.8, color=col, mfc="white", mew=1.3, zorder=3)
        axL.plot([mean], [yy], marker="o", ms=5.8, color=col, mfc=col, mew=1.0, zorder=4)
        value_label(axL, max(mean, med), yy, f"{fnum(mean, 2)} ({nets})", dx=6)
        ticklabels.append(NAME[a])
    top = order[0]
    axL.annotate("median", xy=(per[top]["median"], ys[top]), xytext=(0, 6),
                 textcoords="offset points", ha="center", va="bottom", fontsize=FS, color=INK)
    axL.annotate("mean", xy=(per[top]["mean"], ys[top]), xytext=(0, 6),
                 textcoords="offset points", ha="center", va="bottom", fontsize=FS, color=INK)
    axL.set_yticks([ys[a] for a in order])
    axL.set_yticklabels(ticklabels, fontsize=FS)
    axL.tick_params(axis="y", length=0, pad=4)
    axL.set_xticks([0, 1, 2, 3])
    axL.set_xlabel("edges oriented per claim\n(per network; networks in parentheses)",
                   linespacing=1.15)

    pairs = [("A_ORACLE - D_LLM", "oracle − Qwen 7B"),
             ("A_ORACLE - D_LLM_72B", "oracle − Qwen 72B"),
             ("A_ORACLE - D_LLM_SRC_70B", "oracle − Llama 70B"),
             ("A_TRUE - D_LLM", "true claims − Qwen 7B"),
             ("D_LLM_72B - D_LLM", "Qwen 72B − Qwen 7B")]
    VERDICT = {"separated": ("separated", OI["blue"], "o"),
               "not separated": ("not separated", INK, "o"),
               "reversed": ("reversed", OI["vermilion"], "s")}
    XL, XH = -2.0, 3.5
    axR.set_xlim(XL, XH)
    ysR = [-1.0 - 1.55 * i for i in range(len(pairs))]
    axR.set_ylim(min(ysR) - 0.75, -0.2)
    axR.axvline(0, color=MUTED, lw=0.8, zorder=0)
    axR.spines["left"].set_visible(False)
    axR.set_yticks([])
    for (k, name), yy in zip(pairs, ysR):
        p = ct["paired"][k]
        m = rec(F, f"{name}: mean difference", p["mean_diff"], CASCADE_PATH, f"paired.{k}.mean_diff")
        lo = rec(F, f"{name}: CI low", p["ci"][0], CASCADE_PATH, f"paired.{k}.ci[0]")
        hi = rec(F, f"{name}: CI high", p["ci"][1], CASCADE_PATH, f"paired.{k}.ci[1]")
        n = rec(F, f"{name}: networks", p["networks"], CASCADE_PATH, f"paired.{k}.networks")
        rec(F, f"{name}: verdict", p["verdict"], CASCADE_PATH, f"paired.{k}.verdict")
        txt, col, mk = VERDICT[p["verdict"]]
        axR.text(XL, yy + 0.3, name, ha="left", va="bottom", fontsize=FS_LAB, color=INK)
        axR.text(XH, yy + 0.3, txt, ha="right", va="bottom", fontsize=FS_LAB, color=col,
                 fontweight="bold")
        axR.plot([lo, hi], [yy, yy], color=col, lw=1.5, solid_capstyle="butt", zorder=2)
        axR.plot([m], [yy], marker=mk, ms=5.5, color=col, zorder=3)
        axR.text(XL, yy - 0.26, f"{fnum(m, 2, True)} [{fnum(lo, 2, True)}, {fnum(hi, 2, True)}], "
                 f"{n} networks", ha="left", va="top", fontsize=FS, color=INK2, zorder=4,
                 bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none"))
    axR.set_xticks([-2, -1, 0, 1, 2, 3])
    axR.set_xlabel("paired difference per network\n(95% CI, bootstrap over networks)",
                   linespacing=1.15)
    save(fig, F, preview)


# =============================================================================================
# Figure 4: the estimated CPDAG
# =============================================================================================
def fig_pest(preview):
    F = "pest_dangerous"
    pest = read_csv(PEST_ROWS)
    led = read_csv(LEDGER_ROWS)
    nets = sorted({r["network"] for r in pest})
    arms = ["D_LLM", "D_RAND", "D_DEGEN", "A_TRUE", "A_ORACLE"]

    def count(rows, source, label, filter_txt):
        ok = [r for r in rows if r["status"] == "ok"]
        c = dict(
            cert=len(ok),
            inv=sum(r["z_valid_at_truth"] == "0" for r in ok),
            b7=sum(r["bucket_B7_claim"] == "dangerous" for r in ok),
            b6=sum(r["bucket_B6_hop"] == "dangerous" for r in ok),
        )
        rec(F, f"{label}: certifiable", c["cert"], source, f"{filter_txt}, status=ok (count)")
        rec(F, f"{label}: invalid at the truth", c["inv"], source,
            f"{filter_txt}, status=ok, z_valid_at_truth=0 (count)")
        rec(F, f"{label}: claim radius dangerous", c["b7"], source,
            f"{filter_txt}, status=ok, bucket_B7_claim=dangerous (count)")
        rec(F, f"{label}: hop radius dangerous", c["b6"], source,
            f"{filter_txt}, status=ok, bucket_B6_hop=dangerous (count)")
        return c

    bars = []
    for a in arms:
        rows = [r for r in pest if r["arm"] == a]
        bars.append((a, NAME[a], count(rows, PEST_ROWS, f"{a} (estimated CPDAG)", f"arm={a}"),
                     "est"))
    ref_rows = [r for r in led if r["arm"] == "D_LLM" and r["network"] in nets]
    bars.append(("D_LLM", "Qwen 7B", count(ref_rows, LEDGER_ROWS, "D_LLM (oracle CPDAG)",
                                           f"arm=D_LLM, network in {nets}"), "orc"))
    srcs = sorted({r["cpdag_source"] for r in pest})
    rec(F, "CPDAG of the estimated panel", ";".join(srcs), PEST_ROWS,
        "cpdag_source (distinct values)")
    rec(F, "networks of the estimated panel", len(nets), PEST_ROWS, "network (distinct values)")

    fig = plt.figure(figsize=(W, 7.8 * CM))
    left, gap, right = 0.215, 0.03, 0.02
    wp = (1 - left - right - 2 * gap) / 3
    axs = [fig.add_axes([left + i * (wp + gap), 0.14, wp, 0.70]) for i in range(3)]
    ys, y = [], 0.0
    for a, _, _, kind in bars:
        y -= 1.0
        if a == "A_TRUE":
            y -= 0.8
        if kind == "orc":
            y -= 0.9
        ys.append(y)
    cols = [("inv", "set invalid\nat the truth"), ("b7", "claim radius:\ndangerous"),
            ("b6", "hop radius:\ndangerous")]
    sep1 = ys[2] - 0.55
    sep2 = (ys[4] + ys[5]) / 2 + 0.1
    for ax, (key, title) in zip(axs, cols):
        ax.set_xlim(0, 1.42)
        ax.set_ylim(min(ys) - 0.6, 0.5)
        ax.spines["left"].set_visible(False)
        ax.set_title(title, fontsize=FS_LAB, color=INK, loc="left", pad=4, linespacing=1.1)
        for (a, name, c, kind), yy in zip(bars, ys):
            col = OI["vermilion"] if kind == "est" else OI["blue"]
            v = c[key]
            share = v / c["cert"] if c["cert"] else 0.0
            if v > 0:
                ax.barh(yy, share, height=0.62, color=col, edgecolor="white", linewidth=0.8,
                        hatch=None if kind == "est" else "////", zorder=2)
            else:
                ax.plot([0, 0], [yy - 0.31, yy + 0.31], color=col, lw=2.4, zorder=3,
                        solid_capstyle="butt")
            ax.text(share + 0.03, yy, f"{v} of {c['cert']}", ha="left", va="center",
                    fontsize=FS, color=INK, fontweight="bold" if v == 0 else None)
        ax.set_xticks([0, 0.5, 1])
        ax.set_xticklabels(["0", "0.5", "1"])
        ax.axhline(sep1, color=MUTED, lw=0.7, zorder=1)
        ax.axhline(sep2, color=MUTED, lw=0.7, zorder=1)
    for ax in axs[1:]:
        ax.set_yticks([])
    axs[1].set_xlabel("fraction of the certifiable rows")
    axs[0].set_yticks(ys)
    axs[0].set_yticklabels([n for _, n, _, _ in bars])
    axs[0].tick_params(axis="y", length=0, pad=4)
    header(fig, axs[0], ys[0] + 0.72, f"estimated CPDAG ({len(nets)} networks)")
    header(fig, axs[0], ys[3] + 0.68, "controls that read the truth")
    header(fig, axs[0], sep2 - 0.36, "same networks, oracle CPDAG")
    save(fig, F, preview)


# =============================================================================================
# Figure 5: L95 against the number of wrong claims
# =============================================================================================
PROC = {
    "PLAIN": ("stated report", OI["black"], "o", "-"),
    "SCREEN": ("screen", OI["green"], "s", "-"),
    "STOP1": ("stop r=1", OI["vermilion"], "^", "-"),
    "I1": ("interval r=1", OI["orange"], "D", (0, (4, 2))),
    "STOP2": ("stop r=2", OI["blue"], "v", "-"),
    "I2": ("interval r=2", OI["skyblue"], "d", (0, (4, 2))),
    "BLANKET": ("discard K", "#6E7379", "o", "-"),
    "FAR1": ("retract far", "#B5B9BE", "x", (0, (1.2, 1.6))),
    "SCREEN2": ("screen r=2", OI["purple"], "P", "-"),
}


def spread(vals: list[float], gap: float, lo: float, hi: float) -> list[float]:
    """Push sorted positions apart to a minimum spacing, inside [lo, hi]."""
    idx = sorted(range(len(vals)), key=lambda i: vals[i])
    pos = [vals[i] for i in idx]
    for i in range(1, len(pos)):
        pos[i] = max(pos[i], pos[i - 1] + gap)
    if pos[-1] > hi:
        pos[-1] = hi
        for i in range(len(pos) - 2, -1, -1):
            pos[i] = min(pos[i], pos[i + 1] - gap)
    if pos[0] < lo:
        pos[0] = lo
        for i in range(1, len(pos)):
            pos[i] = max(pos[i], pos[i - 1] + gap)
    out = [0.0] * len(vals)
    for j, i in enumerate(idx):
        out[i] = pos[j]
    return out


def fig_l95(preview):
    F = "l95_crossing"
    t1 = load_json(T1_PATH)
    order = ["PLAIN", "SCREEN", "STOP1", "I1", "STOP2", "I2", "BLANKET", "FAR1"]
    dodge = {p: (i - (len(order) - 1) / 2) * 0.065 for i, p in enumerate(order)}
    YLO, YHI = 0.3, 30.0
    fig = plt.figure(figsize=(W, 9.0 * CM))
    axs = [fig.add_axes([0.075, 0.25, 0.285, 0.68]), fig.add_axes([0.565, 0.25, 0.285, 0.68])]
    for ax, panel in zip(axs, ["A", "B"]):
        cell = t1["controlled"][f"{panel}|1000"]
        ax.set_yscale("log")
        ax.set_ylim(YLO, YHI)
        ax.set_xlim(-0.40, 2.24)
        ax.axhline(1.0, color=INK2, lw=0.9, zorder=1)
        ends = {}
        for p in order:
            name, col, mk, ls = PROC[p]
            xs, ms = [], []
            for b in (0, 1, 2):
                c = cell[f"b={b}"]["procedures"][p]
                base = f"controlled.{panel}|1000.b={b}.procedures.{p}"
                L = rec(F, f"panel {panel}, b={b}, {p}: L95", as_float(c["L95"]), T1_PATH,
                        base + ".L95")
                lo = rec(F, f"panel {panel}, b={b}, {p}: CI low", as_float(c["L95_ci"][0]),
                         T1_PATH, base + ".L95_ci[0]")
                hi = rec(F, f"panel {panel}, b={b}, {p}: CI high", as_float(c["L95_ci"][1]),
                         T1_PATH, base + ".L95_ci[1]")
                assert not math.isinf(L), (panel, b, p)
                x = b + dodge[p]
                xs.append(x)
                ms.append(L)
                if p != "BLANKET":
                    if math.isinf(hi):
                        ax.annotate("", xy=(x, YHI), xytext=(x, max(lo, YLO)),
                                    arrowprops=dict(arrowstyle="-|>", color=col, lw=1.0,
                                                    mutation_scale=7, shrinkA=0, shrinkB=0),
                                    zorder=2, annotation_clip=False)
                    else:
                        ax.plot([x, x], [max(lo, YLO), hi], color=col, lw=1.0, zorder=2,
                                solid_capstyle="butt")
            ax.plot(xs, ms, color=col, lw=1.4 if p != "FAR1" else 1.1, ls=ls, zorder=3)
            ax.plot(xs, ms, ls="none", marker=mk, ms=4.8 if mk not in ("x", "P") else 5.2,
                    color=col, mfc=col if mk != "x" else "none", mew=1.2 if mk == "x" else 0.6,
                    mec="white" if mk != "x" else col, zorder=4)
            ends[p] = (xs[-1], ms[-1])
        logs = [math.log10(ends[p][1]) for p in order]
        pos = spread(logs, 0.115, math.log10(YLO) + 0.04, math.log10(YHI) - 0.04)
        xl = 2.42
        for p, ly in zip(order, pos):
            name, col, mk, ls = PROC[p]
            x0, y0 = ends[p]
            yl = 10 ** ly
            ax.plot([x0 + 0.05, xl - 0.04], [y0, yl], color=col, lw=0.7, zorder=2, clip_on=False)
            ax.text(xl, yl, name, ha="left", va="center", fontsize=FS, color=INK, clip_on=False)
        # label of the reference at 1, in the empty top corner, with a leader to the line
        ax.annotate("1 = width\nwith the\nknowledge\ndiscarded", xy=(-0.33, 1.0),
                    xytext=(-0.36, 24), ha="left", va="top", fontsize=FS, color=INK2,
                    linespacing=1.05,
                    arrowprops=dict(arrowstyle="-|>", color=INK2, lw=0.7, mutation_scale=6,
                                    shrinkA=1, shrinkB=0, relpos=(0.03, 0.0)),
                    zorder=5)
        ax.set_xticks([0, 1, 2])
        tl = []
        for b in (0, 1, 2):
            c = cell[f"b={b}"]
            nrows = rec(F, f"panel {panel}, b={b}: informative rows", c["rows"], T1_PATH,
                        f"controlled.{panel}|1000.b={b}.rows")
            nnets = rec(F, f"panel {panel}, b={b}: networks", c["networks"], T1_PATH,
                        f"controlled.{panel}|1000.b={b}.networks")
            tl.append(f"b = {b}\n{nrows} rows\n{nnets} networks")
        ax.set_xticklabels(tl, linespacing=1.1)
        ax.tick_params(axis="x", length=0, pad=4)
        ax.set_yticks([0.3, 0.5, 1, 2, 5, 10, 20])
        ax.set_yticklabels(["0.3", "0.5", "1", "2", "5", "10", "20"])
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.tick_params(axis="y", which="major", length=2.5, labelsize=FS)
        ax.tick_params(axis="y", which="minor", length=0)
        title = {"A": "panel A, semi-synthetic, n = 1000",
                 "B": "panel B, fitted coefficients, n = 1000"}[panel]
        ax.set_title(title, fontsize=FS_LAB, loc="left", pad=6)
    axs[0].set_ylabel("L95 (width relative to the class, log scale)")
    fig.text(0.5, 0.012, "b = wrong claims (controlled reversals). Bars: 95% CI, bootstrap over "
             "networks; arrow = infinite upper end.", ha="center", va="bottom", fontsize=FS,
             color=INK2)
    save(fig, F, preview)


# =============================================================================================
# Figure 6: distribution of the null radius r0
# =============================================================================================
def fig_r0(preview):
    F = "r0_distribution"
    rows = read_csv(INT_ROWS)
    conds = [("CTRL_b0", "control, b = 0"), ("CTRL_b1", "control, b = 1"), ("D_LLM", "Qwen 7B"),
             ("D_LLM_72B", "Qwen 72B"), ("A_ORACLE", "oracle")]
    cats = [("1", "$r_0$ = 1", OI["vermilion"]), ("2", "2", "#C9CDD2"), ("3", "3", "#8A9097"),
            ("4", "beyond 3 / never", OI["blue"])]
    fig = plt.figure(figsize=(W, 9.0 * CM))
    axA = fig.add_axes([0.16, 0.565, 0.555, 0.32])
    axB = fig.add_axes([0.16, 0.135, 0.555, 0.32])
    FILTER = "source=frame, r0_pop not in ('', '0')"
    for ax, panel in ((axA, "A"), (axB, "B")):
        nets = rec(F, f"panel {panel}: networks (source=frame)",
                   len({r["network"] for r in rows if r["panel"] == panel and r["source"] == "frame"}),
                   INT_ROWS, f"panel={panel}, source=frame, distinct networks")
        ys = [-i for i in range(len(conds))]
        ax.set_xlim(0, 1)
        ax.set_ylim(-len(conds) + 0.5, 0.5)
        ax.spines["left"].set_visible(False)
        tr = blended_transform_factory(ax.transAxes, ax.transData)
        for (cond, name), yy in zip(conds, ys):
            sel = [r for r in rows if r["condition"] == cond and r["panel"] == panel
                   and r["source"] == "frame" and r["r0_pop"] not in ("", "0")]
            n = rec(F, f"panel {panel}, {cond}: rows with a conclusion", len(sel), INT_ROWS,
                    f"condition={cond}, panel={panel}, {FILTER} (count)")
            cnt = {k: rec(F, f"panel {panel}, {cond}: r0_pop={k}",
                          sum(r["r0_pop"] == k for r in sel), INT_ROWS,
                          f"condition={cond}, panel={panel}, {FILTER}, r0_pop={k} (count)")
                   for k, _, _ in cats}
            rec(F, f"panel {panel}, {cond}: of r0_pop=4, censored at depth 3",
                sum(r["r0_pop_status"].startswith("censored") for r in sel), INT_ROWS,
                f"condition={cond}, panel={panel}, {FILTER}, r0_pop_status=censored_at_depth_3")
            assert sum(cnt.values()) == n
            x, outside = 0.0, []
            for k, lab, col in cats:
                s = cnt[k] / n
                rec(F, f"panel {panel}, {cond}: fraction r0_pop={k}", s, INT_ROWS,
                    f"derived: r0_pop={k} / rows with a conclusion")
                if s > 0:
                    ax.barh(yy, s, left=x, height=0.66, color=col, edgecolor="white",
                            linewidth=1.0, zorder=2)
                txt = fpct(s)
                if s > 0 and text_width_data(ax, txt, FS_LAB) + 0.02 < s:
                    ax.text(x + s / 2, yy, txt, ha="center", va="center", fontsize=FS_LAB,
                            color=ink_on(col), zorder=3)
                elif k in ("2", "3"):
                    outside.append(f"{k}: {txt}")
                x += s
            ax.text(1.02, yy, "   ".join([f"n = {n}"] + outside), transform=tr, ha="left",
                    va="center", fontsize=FS, color=INK2)
            if panel == "A" and cond == conds[0][0]:
                xc, x = [], 0.0
                for k, _, _ in cats:
                    s = cnt[k] / n
                    xc.append(x + s / 2)
                    x += s
                xt = [xc[0], xc[1] - 0.03, xc[2] + 0.03, xc[3]]
                for (k, lab, _), a0, a1 in zip(cats, xc, xt):
                    ax.plot([a0, a1], [yy + 0.35, yy + 0.55], color=INK2, lw=0.6, zorder=3,
                            clip_on=False)
                    ax.text(a1, yy + 0.58, lab, ha="center", va="bottom", fontsize=FS_LAB,
                            color=INK, clip_on=False)
        ax.set_yticks(ys)
        ax.set_yticklabels([nm for _, nm in conds])
        ax.tick_params(axis="y", length=0, pad=4)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
        ax.set_xticklabels(["0", "25%", "50%", "75%", "100%"])
        title = {"A": f"panel A, semi-synthetic ({nets} networks)",
                 "B": f"panel B, fitted coefficients ({nets} networks)"}[panel]
        fig.text(0.012, ax.get_position().y1 + (0.075 if panel == "A" else 0.03), title,
                 ha="left", va="bottom", fontsize=FS_LAB, fontweight="bold", color=INK)
    axB.set_xlabel("fraction of the rows in which the effect is already established with all the "
                   "knowledge ($r_0$ ≥ 1)")
    save(fig, F, preview)


# =============================================================================================
# Figure 7: the stop rule, refusal at b = 0 and the starred comparison at b = 2
# =============================================================================================
def fig_stop(preview):
    F = "stop_refusal"
    t1 = load_json(T1_PATH)
    m1b = t1["m1b"]
    fig = plt.figure(figsize=(W, 8.6 * CM))
    axL = fig.add_axes([0.13, 0.22, 0.30, 0.58])
    axR = fig.add_axes([0.575, 0.22, 0.41, 0.60])

    # --- left: refusal rate at b = 0, informative rows with a committed set
    pols = ["STOP1", "STOP2", "SCREEN", "SCREEN2"]
    ys = [-i for i in range(len(pols))]
    axL.set_xlim(0.6, 1.025)
    axL.set_ylim(-len(pols) + 0.4, 0.5)
    axL.spines["left"].set_visible(False)
    info = {}
    for yy in ys:
        axL.axhline(yy, color=GRID, lw=0.6, zorder=0)
    for panel, off, filled in (("A", 0.17, True), ("B", -0.17, False)):
        key = f"CTRL_b0|committed|{panel}"
        cell = m1b["refusal"][key]
        info[panel] = (rec(F, f"refusal b=0, panel {panel}: rows", cell["rows"], T1_PATH,
                           f"m1b.refusal.{key}.rows"),
                       rec(F, f"refusal b=0, panel {panel}: networks", cell["networks"], T1_PATH,
                           f"m1b.refusal.{key}.networks"))
        for p, yy in zip(pols, ys):
            name, col, mk, _ = PROC[p]
            r = rec(F, f"refusal b=0, panel {panel}, {p}", cell[p]["rate"], T1_PATH,
                    f"m1b.refusal.{key}.{p}.rate")
            lo = rec(F, f"refusal b=0, panel {panel}, {p}: CI low", cell[p]["ci"][0], T1_PATH,
                     f"m1b.refusal.{key}.{p}.ci[0]")
            hi = rec(F, f"refusal b=0, panel {panel}, {p}: CI high", cell[p]["ci"][1], T1_PATH,
                     f"m1b.refusal.{key}.{p}.ci[1]")
            y = yy + off
            axL.plot([lo, hi], [y, y], color=col, lw=1.3, solid_capstyle="butt", zorder=2)
            axL.plot([r], [y], marker="o" if panel == "A" else "s", ms=5.0, color=col,
                     mfc=col if filled else "white", mew=1.3, zorder=3)
            axL.annotate(f"{panel} {fnum(r)}", xy=(lo, y), xytext=(-4, 0),
                         textcoords="offset points", ha="right", va="center", fontsize=FS,
                         color=INK2)
    axL.set_yticks(ys)
    axL.set_yticklabels([PROC[p][0] for p in pols])
    axL.tick_params(axis="y", length=0, pad=4)
    axL.set_xticks([0.6, 0.7, 0.8, 0.9, 1.0])
    axL.set_xticklabels([fnum(v, 1) for v in (0.6, 0.7, 0.8, 0.9, 1.0)])
    axL.set_xlabel("refusal rate with correct knowledge (b = 0)")
    axL.text(-0.38, 1.20, "the price of the guarantee", transform=axL.transAxes, ha="left",
             va="bottom", fontsize=FS_LAB, fontweight="bold")
    axL.text(-0.38, 1.03, f"rows with a committed set\nA: {info['A'][0]} rows, "
             f"{info['A'][1]} networks · B: {info['B'][0]} rows, {info['B'][1]} networks",
             transform=axL.transAxes, ha="left", va="bottom", fontsize=FS, color=INK2,
             linespacing=1.1)

    # --- right: raw coverage at b = 2, n = 1000, rows with a committed set
    procs = ["SCREEN", "SCREEN2", "STOP2"]
    axR.set_ylim(0.6, 1.30)
    axR.set_xlim(-0.55, 6.95)
    axR.axhline(0.95, color=MUTED, lw=0.8, zorder=0)
    axR.text(-0.5, 0.955, "0.95", ha="left", va="bottom", fontsize=FS, color=MUTED)
    ticks, labs = [], []
    for g, panel in enumerate(("A", "B")):
        k = f"committed|{panel}|1000"
        st, dm = m1b["starred"][k], m1b["starred_depth_matched"][k]
        base_x = g * 4.3
        cov = {
            "SCREEN": rec(F, f"b=2 panel {panel}: coverage SCREEN", st["coverage"]["SCREEN"],
                          T1_PATH, f"m1b.starred.{k}.coverage.SCREEN"),
            "SCREEN2": rec(F, f"b=2 panel {panel}: coverage SCREEN2", dm["coverage"]["SCREEN2"],
                           T1_PATH, f"m1b.starred_depth_matched.{k}.coverage.SCREEN2"),
            "STOP2": rec(F, f"b=2 panel {panel}: coverage STOP2", st["coverage"]["STOP2"],
                         T1_PATH, f"m1b.starred.{k}.coverage.STOP2"),
        }
        assert abs(dm["coverage"]["STOP2"] - cov["STOP2"]) < 1e-12
        nrows = rec(F, f"b=2 panel {panel}: rows", st["rows"], T1_PATH, f"m1b.starred.{k}.rows")
        nnets = rec(F, f"b=2 panel {panel}: networks", st["networks"], T1_PATH,
                    f"m1b.starred.{k}.networks")
        assert nrows == dm["rows"] and nnets == dm["networks"]
        xs = {p: base_x + i * 1.05 for i, p in enumerate(procs)}
        for p in procs:
            name, col, mk, _ = PROC[p]
            axR.plot([xs[p]], [cov[p]], marker=mk, ms=7.0, color=col, mec="white", mew=0.6,
                     zorder=3)
            below = p == "SCREEN"  # below the marker, so as not to cross the 0.95 line
            axR.annotate(fnum(cov[p]), xy=(xs[p], cov[p]), xytext=(0, -6 if below else 4.5),
                         textcoords="offset points", ha="center",
                         va="top" if below else "bottom", fontsize=FS,
                         color=INK2, zorder=4,
                         bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none"))
            ticks.append(xs[p])
            labs.append(name.replace(" r=2", "\nr=2"))
        d1 = rec(F, f"b=2 panel {panel}: STOP2 − SCREEN", st["coverage_diff"], T1_PATH,
                 f"m1b.starred.{k}.coverage_diff")
        d1lo = rec(F, f"b=2 panel {panel}: STOP2 − SCREEN CI low", st["coverage_diff_ci"][0],
                   T1_PATH, f"m1b.starred.{k}.coverage_diff_ci[0]")
        d1hi = rec(F, f"b=2 panel {panel}: STOP2 − SCREEN CI high", st["coverage_diff_ci"][1],
                   T1_PATH, f"m1b.starred.{k}.coverage_diff_ci[1]")
        d2 = rec(F, f"b=2 panel {panel}: STOP2 − SCREEN2", dm["coverage_diff"], T1_PATH,
                 f"m1b.starred_depth_matched.{k}.coverage_diff")
        d2lo = rec(F, f"b=2 panel {panel}: STOP2 − SCREEN2 CI low", dm["coverage_diff_ci"][0],
                   T1_PATH, f"m1b.starred_depth_matched.{k}.coverage_diff_ci[0]")
        d2hi = rec(F, f"b=2 panel {panel}: STOP2 − SCREEN2 CI high", dm["coverage_diff_ci"][1],
                   T1_PATH, f"m1b.starred_depth_matched.{k}.coverage_diff_ci[1]")

        def bracket(x0, x1, y, txt, bold=False):
            axR.plot([x0, x0, x1, x1], [y - 0.012, y, y, y - 0.012], color=INK, lw=0.8, zorder=3)
            axR.text((x0 + x1) / 2, y + 0.007, txt, ha="center", va="bottom", fontsize=FS,
                     color=INK, fontweight="bold" if bold else None, linespacing=1.05)

        top = max(cov.values())
        bracket(xs["SCREEN2"], xs["STOP2"], top + 0.088,
                f"{fnum(d2, 3, True)}\n[{fnum(d2lo, 3, True)}, {fnum(d2hi, 3, True)}]")
        bracket(xs["SCREEN"], xs["STOP2"], top + 0.208,
                f"{fnum(d1, 3, True)}\n[{fnum(d1lo, 3, True)}, {fnum(d1hi, 3, True)}]", bold=True)
        axR.text((xs["SCREEN"] + xs["STOP2"]) / 2, -0.23,
                 f"panel {panel}\n{nrows} rows, {nnets} networks", linespacing=1.1,
                 transform=blended_transform_factory(axR.transData, axR.transAxes),
                 ha="center", va="top", fontsize=FS, color=INK)
    axR.set_xticks(ticks)
    axR.set_xticklabels(labs, linespacing=1.0)
    axR.tick_params(axis="x", length=0, pad=3)
    axR.set_yticks([0.6, 0.7, 0.8, 0.9, 1.0])
    axR.set_yticklabels([fnum(v, 1) for v in (0.6, 0.7, 0.8, 0.9, 1.0)])
    axR.tick_params(axis="y", length=2.5, labelsize=FS)
    axR.set_ylabel("raw coverage, b = 2, n = 1000")
    axR.text(-0.13, 1.20, "stop r=2 against the two screens", transform=axR.transAxes, ha="left",
             va="bottom", fontsize=FS_LAB, fontweight="bold")
    axR.text(-0.13, 1.03, "brackets: paired difference per network\n(stop minus screen), 95% CI",
             transform=axR.transAxes, ha="left", va="bottom", fontsize=FS, color=INK2,
             linespacing=1.1)
    save(fig, F, preview)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preview", type=Path, default=None)
    args = ap.parse_args()
    fig_buckets(args.preview)
    fig_scale(args.preview)
    fig_cascade(args.preview)
    fig_pest(args.preview)
    fig_l95(args.preview)
    fig_r0(args.preview)
    fig_stop(args.preview)
    out = {
        "generated_by": "gen_figs.py",
        "note": "Every plotted value, with the result file and the key (JSON) or the filter (CSV) "
                "it came from. 'inf' marks an infinite bootstrap upper end.",
        "figures": {f"figs/{k}.pdf": v for k, v in NUM.items()},
    }
    (FIGS / "numbers_gen_figs.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    for k, v in NUM.items():
        print(f"figs/{k}.pdf: {len(v)} numbers")


if __name__ == "__main__":
    main()
