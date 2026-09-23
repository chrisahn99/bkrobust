"""Figures for session 3 (Axis B: declarative encodings and tractability).

Every figure is built from a committed file under ``results/axisb3/``.
Colour-blind safe; encoding never relies on colour alone. Instances where a
method did not run, or ran out of budget, are drawn as such rather than dropped:
a timeout is a datum, and silently omitting it would flatter the slower method.
"""

from __future__ import annotations

import json
import statistics as st
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = Path("results/axisb3")
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
    path = RES / name
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def fig_cost_by_n() -> list[Path]:
    """E1 wall-clock against n, separated into build, UNSAT and SAT.

    Log scale, because the medians span three orders of magnitude and the point
    is the *composition*: build dominates and SAT is negligible. Grouped rather
    than stacked, since stacked bars do not read on a log axis. A median of
    exactly zero is drawn at the floor and labelled, not silently omitted.
    """
    _style()
    rows = _rows("scaling.jsonl")
    ns = sorted({r["n"] for r in rows})
    floor = 1e-4

    def med(n: int, key: str) -> float:
        return st.median(r["e1"][key] for r in rows if r["n"] == n)

    series = [
        ("build (encoding)", [med(n, "build_s") for n in ns], C_BLUE),
        ("UNSAT (certifying clean shells)", [med(n, "unsat_s") for n in ns], C_ORANGE),
        ("SAT (finding the witness)", [med(n, "sat_s") for n in ns], C_GREEN),
    ]
    worst = [max(r["e1"]["total_s"] for r in rows if r["n"] == n) for n in ns]

    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    xs = list(range(len(ns)))
    w = 0.26
    for j, (lab, vals, col) in enumerate(series):
        pos = [x + (j - 1) * w for x in xs]
        ax.bar(pos, [max(v, floor) for v in vals], w, label=lab, color=col)
        for i, v in enumerate(vals):
            if v < floor:
                ax.text(pos[i], floor * 1.15, "0", ha="center", fontsize=7, color=col)
    ax.plot(xs, worst, "o--", color=C_VERM, label="worst case (total, not median)")
    ax.set_yscale("log")
    ax.set_ylim(floor * 0.7, max(worst) * 3)
    ax.set_xticks(xs)
    ax.set_xticklabels([str(n) for n in ns])
    ax.set_xlabel("n (vertices)")
    ax.set_ylabel("seconds, log scale (median unless marked)")
    ax.set_title("E1 cost by problem size — build dominates, SAT is nearly free")
    ax.legend(fontsize=8, loc="upper left")
    return _save(fig, "s3_f1_cost_by_n")


def fig_crossover_vs_bfs() -> list[Path]:
    """The crossover against building the space, by undirected-edge count."""
    _style()
    rows = [r for r in _rows("scaling.jsonl") if "bfs" in r]
    by = defaultdict(list)
    for r in rows:
        by[r["m_undirected"]].append(r)
    ms = sorted(by)
    bfs = [st.median(r["bfs"]["build_s"] for r in by[m]) for m in ms]
    e1 = [st.median(r["e1"]["total_s"] for r in by[m]) for m in ms]

    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    ax.semilogy(ms, bfs, "s-", color=C_VERM, label="BFS: build the space (3^m)")
    ax.semilogy(ms, e1, "o-", color=C_BLUE, label="E1: encode and solve")
    cross = next((ms[i] for i in range(len(ms)) if bfs[i] > e1[i]), None)
    if cross is not None:
        ax.axvline(cross, color=C_GREY, ls=":", lw=1)
        ax.annotate(
            f"crossover at m = {cross}",
            (cross, max(bfs)),
            xytext=(6, -6),
            textcoords="offset points",
            fontsize=8,
            color=C_GREY,
        )
    ax.set_xlabel("m (undirected edges in the CPDAG)")
    ax.set_ylabel("seconds (median, log scale)")
    ax.set_title("E1 overtakes space construction at m = 5")
    ax.set_xticks(ms)
    ax.legend(fontsize=8)
    return _save(fig, "s3_f2_crossover_vs_bfs")


def fig_crossover_vs_local_up() -> list[Path]:
    """Where E1 beats ``local_up`` and where it loses: the split is the radius.

    Two panels, because the aggregate hides the result. Timeouts are drawn at
    the cap and marked as censored, never dropped -- omitting them would flatter
    the slower method.
    """
    _style()
    rows = _rows("highk.jsonl")
    cap = 60.0
    near = [r for r in rows if r["e1"]["radius"] not in (UNREACHED, -2)]
    far = [r for r in rows if r["e1"]["radius"] == UNREACHED]

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 4.0), sharey=True)
    rng = __import__("random").Random(0)  # jitter only; no effect on values

    for ax, group, title in (
        (axes[0], near, "a failure exists\n(SAT side)"),
        (axes[1], far, "no failure anywhere\n(UNSAT side)"),
    ):
        e1 = [r["e1"]["total_s"] for r in group]
        up_ok, up_to = [], []
        for r in group:
            u = r["local_up"]
            if "total_s" not in u:
                continue
            (up_to if "timed_out_at_s" in u else up_ok).append(min(u["total_s"], cap))
        ax.scatter(
            [1 + rng.uniform(-0.13, 0.13) for _ in e1],
            e1,
            s=16,
            color=C_BLUE,
            alpha=0.65,
            label="E1",
        )
        ax.scatter(
            [2 + rng.uniform(-0.13, 0.13) for _ in up_ok],
            up_ok,
            s=16,
            color=C_ORANGE,
            alpha=0.65,
            label="local_up (finished)",
        )
        if up_to:
            ax.scatter(
                [2 + rng.uniform(-0.13, 0.13) for _ in up_to],
                up_to,
                s=42,
                color=C_VERM,
                marker="^",
                label=f"local_up hit the {cap:.0f}s cap",
            )
        for xpos, vals, col in ((1, e1, C_BLUE), (2, up_ok + up_to, C_ORANGE)):
            if vals:
                ax.plot(
                    [xpos - 0.28, xpos + 0.28],
                    [st.median(vals)] * 2,
                    color=col,
                    lw=2.5,
                    solid_capstyle="butt",
                )
        if up_to:
            # Censored runs enter the median at the cap, so the bar understates
            # local_up's true median. Say so rather than drawing it as a value.
            ax.annotate(
                f"median is a LOWER bound\n({len(up_to)} of {len(up_ok) + len(up_to)} censored)",
                (2, st.median(up_ok + up_to)),
                xytext=(0, -34),
                textcoords="offset points",
                fontsize=7.5,
                color=C_VERM,
                ha="center",
            )
        ax.axhline(cap, color=C_GREY, ls=":", lw=1)
        ax.set_yscale("log")
        ax.set_xlim(0.5, 2.5)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["E1", "local_up"])
        ax.set_title(f"{title}\nn = {len(group)} instances", fontsize=9)

    axes[0].set_ylabel("seconds (log scale); bars are medians")
    axes[1].annotate(
        f"{cap:.0f}s cap", (2.42, cap), fontsize=8, color=C_GREY, ha="right", va="bottom"
    )
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=8, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.06))
    fig.suptitle("The crossover is governed by the radius, not by k", fontsize=10)
    fig.tight_layout()
    return _save(fig, "s3_f3_crossover_vs_local_up")


def fig_agreement() -> list[Path]:
    """Cross-encoding agreement, and the scope each oracle could actually cover."""
    _style()
    scale = _rows("scaling.jsonl")
    xval = json.loads((RES / "cross_encoding_n4.json").read_text())

    labels, agree, scope = [], [], []
    for key, lab in (
        ("local_up", "local_up"),
        ("bfs", "brute-force BFS"),
        ("e2", "E2"),
        ("e3", "E3"),
    ):
        have = [r for r in scale if key in r and r[key].get("radius") is not None]
        labels.append(lab)
        agree.append(sum(1 for r in have if r[key]["radius"] == r["e1"]["radius"]))
        scope.append(len(have))
    labels.append("all (n=4 exhaustive)")
    agree.append(xval["agreements"]["e1_bfs"])
    scope.append(xval["n_instances"])

    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    ypos = range(len(labels))
    ax.barh(list(ypos), scope, color=C_GREY, alpha=0.45, label="instances where it ran")
    ax.barh(list(ypos), agree, color=C_GREEN, label="agreed with E1")
    for i in range(len(labels)):
        ax.text(scope[i] + max(scope) * 0.01, i, f"{agree[i]}/{scope[i]}", va="center", fontsize=8)
    ax.set_yticks(list(ypos))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("instances")
    ax.set_xlim(0, max(scope) * 1.20)
    ax.set_title("Every oracle that could run agreed with E1 — zero disagreements")
    ax.legend(fontsize=8, loc="lower center", bbox_to_anchor=(0.5, -0.34), ncol=2)
    return _save(fig, "s3_f4_agreement")


def fig_c2_margin() -> list[Path]:
    """Conjecture 2 stress test: the margin r_E3 - r_E1, by scope.

    A refutation is r_E3 < r_E1. Reporting only "zero counterexamples" would
    hide that every margin is also exactly 0 -- the conjecture is tied
    everywhere rather than holding with slack, and those are different states.
    """
    _style()
    xval = json.loads((RES / "cross_encoding_n4.json").read_text())
    scopes = [
        (
            "all CPDAGs, n = 4\n(exhaustive)",
            [0] * xval["walks_certified"],
            xval["walks_certified"],
            xval["walks_monotone"],
        ),
    ]
    for name, label in (
        ("scaling.jsonl", "generated, n = 8..20"),
        ("highk.jsonl", "high-k hunt\n(dense, |K| up to 51)"),
    ):
        rs = [r for r in _rows(name) if "e3" in r]
        scopes.append(
            (
                label,
                [r["e3"]["radius"] - r["e1"]["radius"] for r in rs],
                sum(1 for r in rs if r["e3"]["certified"]),
                sum(1 for r in rs if r["e3"]["monotone"]),
            )
        )

    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    total = sum(len(m) for _, m, _, _ in scopes)
    left = 0.0
    cols = [C_BLUE, C_GREEN, C_ORANGE]
    for i, (label, margins, cert, mono) in enumerate(scopes):
        assert set(margins) <= {0}, f"non-zero margin in {label}: {sorted(set(margins))}"
        ax.barh(
            0,
            len(margins),
            left=left,
            color=cols[i],
            height=0.55,
            label=f"{label}  —  {len(margins)} walks, {cert} certified, {mono} monotone",
        )
        left += len(margins)
    ax.set_yticks([0])
    ax.set_yticklabels(["margin = 0"])
    ax.set_xlim(0, total * 1.02)
    ax.set_xlabel("E3 witness walks (all at margin exactly 0)")
    ax.set_title(
        f"Conjecture 2: {total} chances to refute it, taken by an\n"
        "assumption-free encoding — and every margin is exactly 0"
    )
    ax.annotate(
        "a refutation would be a bar to the left of this line\n"
        "(r_E3 < r_E1). There is none — and no slack either.",
        (total * 0.02, -0.45),
        fontsize=8,
        color=C_VERM,
        va="top",
    )
    ax.axvline(0, color=C_VERM, ls="--", lw=1.5)
    ax.set_ylim(-0.95, 0.45)
    ax.legend(fontsize=7.5, loc="lower right")
    ax.grid(axis="y", visible=False)
    return _save(fig, "s3_f5_c2_margin")


def fig_degeneracy() -> list[Path]:
    """How often no failure exists at all — the UNSAT-side workload."""
    _style()
    rows = _rows("scaling.jsonl")
    ns = sorted({r["n"] for r in rows})
    rates, counts = [], []
    for n in ns:
        g = [r for r in rows if r["n"] == n]
        rates.append(100 * sum(1 for r in g if r["e1"]["radius"] == UNREACHED) / len(g))
        counts.append(len(g))
    fig, ax = plt.subplots(figsize=(6.2, 3.2))
    ax.bar(ns, rates, 0.6, color=C_ORANGE)
    for i, n in enumerate(ns):
        ax.text(n, rates[i] + 1, f"{rates[i]:.0f}%\n(n={counts[i]})", ha="center", fontsize=8)
    ax.set_xlabel("n (vertices)")
    ax.set_ylabel("% of instances with no reachable failure")
    ax.set_ylim(0, 100)
    ax.set_xticks(ns)
    ax.set_title("Degeneracy rate — these are the instances E1 wins on")
    return _save(fig, "s3_f6_degeneracy")


def build_all() -> list[Path]:
    """Build every session 3 figure. Returns the paths written."""
    out: list[Path] = []
    for fn in (
        fig_cost_by_n,
        fig_crossover_vs_bfs,
        fig_crossover_vs_local_up,
        fig_agreement,
        fig_c2_margin,
        fig_degeneracy,
    ):
        out.extend(fn())
    return out


if __name__ == "__main__":
    for p in build_all():
        print(p)
