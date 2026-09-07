"""Figures for session 4 (Axis B: the oracle, and the cost of enumeration).

Every figure is built from a committed file under ``results/axisb4/``.
Colour-blind safe; encoding never relies on colour alone. Censored runs are
drawn at their cap and labelled, never dropped -- omitting a timeout would
flatter the slower method.
"""

# ruff: noqa: RUF001
# Figure text is display copy; the multiplication sign is intentional.

from __future__ import annotations

import json
import statistics as st
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = Path("results/axisb4")
FIGDIR = Path("figures")

C_BLUE = "#0072B2"
C_ORANGE = "#E69F00"
C_GREEN = "#009E73"
C_VERM = "#D55E00"
C_GREY = "#999999"
UNREACHED = -1


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


def _rows(name: str) -> list[dict]:
    return [json.loads(line) for line in (RES / name).read_text().splitlines() if line.strip()]


def _profiled() -> list[dict]:
    return [
        r for r in _rows("profile_local_up.jsonl") if "total_s" in r and "timed_out_at_s" not in r
    ]


def fig_profile_split() -> list[Path]:
    """Where local_up's time goes, overall and by side.

    The point of the figure: it is enumeration-bound, but the enumeration is
    mostly cover minimality, not the validity oracle -- which is what caps the
    MPDAG criterion's speed value on its own.
    """
    _style()
    rows = _profiled()
    groups = [
        ("all\ninstances", rows),
        ("SAT side\n(a failure exists)", [r for r in rows if r["side"] == "SAT"]),
        ("UNSAT side\n(no failure)", [r for r in rows if r["side"] == "UNSAT"]),
    ]
    parts = [
        ("cover minimality\n(enumerate_dag_extensions)", "cover_ext_s", C_VERM),
        ("validity oracle\n(enumerate + check)", "oracle_s", C_ORANGE),
        ("meek_closure", "closure_s", C_GREEN),
        ("search overhead", "overhead_s", C_GREY),
    ]
    fig, ax = plt.subplots(figsize=(7.0, 3.8))
    xs = list(range(len(groups)))
    bottoms = [0.0] * len(groups)
    for label, key, col in parts:
        vals = []
        for _, g in groups:
            tot = sum(r["total_s"] for r in g)
            vals.append(100 * sum(r[key] for r in g) / tot if tot else 0.0)
        ax.bar(xs, vals, 0.55, bottom=bottoms, label=label, color=col)
        for i, v in enumerate(vals):
            if v >= 4:
                ax.text(
                    i,
                    bottoms[i] + v / 2,
                    f"{v:.1f}%",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="white" if col != C_GREY else "black",
                )
        bottoms = [bottoms[i] + vals[i] for i in range(len(groups))]
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{lab}\nn = {len(g)}" for lab, g in groups], fontsize=8)
    ax.set_ylabel("share of total wall-clock (%)")
    ax.set_ylim(0, 100)
    ax.set_title("local_up is enumeration-bound — but mostly NOT in the oracle")
    ax.legend(fontsize=7.5, loc="lower center", bbox_to_anchor=(0.5, -0.42), ncol=2)
    ax.grid(axis="x", visible=False)
    return _save(fig, "s4_f1_profile_split")


def fig_amdahl() -> list[Path]:
    """What each lever can buy, from the measured shares.

    Stated as a ceiling, not a promise: removing a component that is a fraction
    f of the time cannot beat 1/(1-f), whatever replaces it.
    """
    _style()
    rows = [r for r in _profiled() if r["side"] == "UNSAT"]
    tot = sum(r["total_s"] for r in rows)
    o = sum(r["oracle_s"] for r in rows) / tot
    c = sum(r["cover_ext_s"] for r in rows) / tot
    levers = [
        ("MPDAG criterion\nalone\n(replaces the oracle)", o, C_ORANGE),
        ("Lemma O alone\n(removes cover\nminimality)", c, C_VERM),
        ("both\ncomposed", o + c, C_BLUE),
    ]
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    xs = list(range(len(levers)))
    ceils = [1 / (1 - f) for _, f, _ in levers]
    ax.bar(xs, ceils, 0.5, color=[c3 for _, _, c3 in levers])
    for i in range(len(levers)):
        ax.text(
            i,
            ceils[i] + max(ceils) * 0.03,
            f"{ceils[i]:.1f}×\nceiling",
            ha="center",
            fontsize=9,
        )
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{lab}\n\nremoves {100 * f:.1f}%" for lab, f, _ in levers], fontsize=8)
    ax.set_ylabel("Amdahl ceiling on speedup (×)")
    ax.set_ylim(0, max(ceils) * 1.28)
    ax.set_title("What each lever can buy on the UNSAT side, from measured shares")
    ax.grid(axis="x", visible=False)
    return _save(fig, "s4_f2_amdahl")


def fig_fast_vs_frozen() -> list[Path]:
    """Measured effect of Lemma O, against the ceiling the profile predicted."""
    _style()
    rows = _rows("fast_vs_frozen.jsonl")
    both = [
        r
        for r in rows
        if "total_s" in r["slow"]
        and "total_s" in r["fast"]
        and "timed_out_at_s" not in r["slow"]
        and "timed_out_at_s" not in r["fast"]
    ]
    resc = [
        r
        for r in rows
        if "timed_out_at_s" in r["slow"]
        and "total_s" in r["fast"]
        and "timed_out_at_s" not in r["fast"]
    ]

    fig, ax = plt.subplots(figsize=(6.4, 3.9))
    for side, col, mark in (("SAT", C_BLUE, "o"), ("UNSAT", C_VERM, "^")):
        g = [r for r in both if r["fast"]["side"] == side]
        if g:
            ax.scatter(
                [r["slow"]["total_s"] for r in g],
                [r["fast"]["total_s"] for r in g],
                s=20,
                alpha=0.7,
                color=col,
                marker=mark,
                label=f"{side} side (n = {len(g)})",
            )
    if resc:
        ax.scatter(
            [120.0] * len(resc),
            [r["fast"]["total_s"] for r in resc],
            s=70,
            marker="x",
            color=C_GREEN,
            linewidths=2,
            label=f"frozen hit the 120s cap; fast finished (n = {len(resc)})",
        )
    lo = min(r["fast"]["total_s"] for r in both) / 2
    hi = 200.0
    ax.plot([lo, hi], [lo, hi], ls=":", color=C_GREY, lw=1)
    ax.annotate("equal", (hi, hi), fontsize=8, color=C_GREY, ha="right", va="bottom")
    ax.plot([lo, hi], [lo / 3.9, hi / 3.9], ls="--", color=C_GREY, lw=1)
    ax.annotate("3.9× faster", (hi, hi / 3.9), fontsize=8, color=C_GREY, ha="right", va="bottom")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("frozen local_up (seconds, log)")
    ax.set_ylabel("enumeration-free local_up (seconds, log)")
    ax.set_title("Lemma O: same radii everywhere, 3.9× less time in aggregate")
    ax.legend(fontsize=7.5, loc="upper left")
    return _save(fig, "s4_f3_fast_vs_frozen")


def fig_cost_by_component() -> list[Path]:
    """Cost against undirected-component size — the axis Axis A will need."""
    _style()
    rows = _rows("fast_vs_frozen.jsonl")
    by_slow: dict[int, list[float]] = defaultdict(list)
    by_fast: dict[int, list[float]] = defaultdict(list)
    for r in rows:
        for key, store in (("slow", by_slow), ("fast", by_fast)):
            e = r[key]
            if "total_s" in e and "timed_out_at_s" not in e:
                store[e["largest_component"]].append(e["total_s"])
    comps = sorted(set(by_slow) | set(by_fast))
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    ax.semilogy(
        comps, [st.median(by_slow[c]) for c in comps], "s-", color=C_VERM, label="frozen local_up"
    )
    ax.semilogy(
        comps,
        [st.median(by_fast[c]) for c in comps],
        "o-",
        color=C_BLUE,
        label="enumeration-free local_up",
    )
    ax.set_xlabel("largest undirected component of the CPDAG (vertices)")
    ax.set_ylabel("seconds (median, log scale)")
    ax.set_xticks(comps)
    ax.set_title("Cost against component size, not n — the axis that matters")
    ax.legend(fontsize=8)
    return _save(fig, "s4_f4_cost_by_component")


def fig_hybrid_envelope() -> list[Path]:
    """The operating envelope, on the axis Axis A will need.

    Component size, not n. Censored frozen runs are drawn at the cap and marked;
    the hybrid has none to draw.
    """
    _style()
    rows = _rows("hybrid_envelope.jsonl")
    cap = 120.0

    def ok(e: dict) -> bool:
        return "total_s" in e and "timed_out_at_s" not in e

    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    hx = [r["hybrid"]["largest_component"] for r in rows if ok(r["hybrid"])]
    hy = [r["hybrid"]["total_s"] for r in rows if ok(r["hybrid"])]
    fx = [r["frozen"]["largest_component"] for r in rows if ok(r["frozen"])]
    fy = [r["frozen"]["total_s"] for r in rows if ok(r["frozen"])]
    ax.scatter(
        fx, fy, s=18, alpha=0.55, color=C_VERM, marker="s", label=f"frozen local_up (n = {len(fx)})"
    )
    ax.scatter(hx, hy, s=18, alpha=0.7, color=C_BLUE, marker="o", label=f"hybrid (n = {len(hx)})")
    resc = [r for r in rows if "timed_out_at_s" in r["frozen"] and ok(r["hybrid"])]
    if resc:
        ax.scatter(
            [r["hybrid"]["largest_component"] for r in resc],
            [cap] * len(resc),
            s=80,
            marker="x",
            color=C_VERM,
            linewidths=2,
            label=f"frozen hit the {cap:.0f}s cap (n = {len(resc)})",
        )
        for r in resc:
            ax.annotate(
                "",
                xy=(r["hybrid"]["largest_component"], r["hybrid"]["total_s"]),
                xytext=(r["hybrid"]["largest_component"], cap),
                arrowprops={"arrowstyle": "->", "color": C_GREY, "lw": 0.8},
            )
    ax.axhline(cap, color=C_GREY, ls=":", lw=1)
    ax.set_yscale("log")
    ax.set_xlabel("largest undirected component of the CPDAG (vertices)")
    ax.set_ylabel("seconds to an exact radius (log scale)")
    ax.set_title(
        "Operating envelope: the hybrid finished every instance;\n"
        "arrows mark the ones the frozen search could not"
    )
    ax.legend(fontsize=8, loc="lower right")
    return _save(fig, "s4_f5_hybrid_envelope")


def build_all() -> list[Path]:
    """Build every session 4 figure that has data. Returns the paths written."""
    out: list[Path] = []
    for fn in (
        fig_profile_split,
        fig_amdahl,
        fig_fast_vs_frozen,
        fig_cost_by_component,
        fig_hybrid_envelope,
    ):
        try:
            out.extend(fn())
        except (FileNotFoundError, ValueError, IndexError) as exc:
            print(f"skipped {fn.__name__}: {exc}")
    return out


if __name__ == "__main__":
    for p in build_all():
        print(p)
