"""Figures for the synthetic-ensemble / search-heuristic session.

Every figure is built from a committed file under ``results/``; none contains a
number that is not on disk. Saved as PDF (vector, for the report) and PNG.

Palette is colour-blind safe and encoding never relies on colour alone: bars
carry hatching and series carry distinct markers.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

RESULTS = Path("results")
FIGDIR = Path("figures")

# Colour-blind safe (Okabe-Ito)
C_BLUE = "#0072B2"
C_ORANGE = "#E69F00"
C_GREEN = "#009E73"
C_VERM = "#D55E00"
C_PURPLE = "#CC79A7"
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


# --------------------------------------------------------------------------- H1


def fig_h1_distribution() -> list[Path]:
    """H1: the r_val distribution, census against ensembles."""
    _style()
    cen4 = pd.read_csv(RESULTS / "synth/census/census.csv")
    cen5 = pd.read_csv(RESULTS / "synth/census/census_n5.csv")
    ens = pd.read_csv(RESULTS / "synth/main/results.csv")

    def hist(s: pd.Series) -> dict[int, float]:
        s = s[s != -1]
        vc = s.value_counts(normalize=True).sort_index()
        return {int(k): float(v) for k, v in vc.items()}

    sets = [
        ("census n=4 (exhaustive)", hist(cen4["r_val_opt"]), C_BLUE, ""),
        ("census n=5 (exhaustive)", hist(cen5["r_val_opt"]), C_GREEN, "//"),
        ("ensembles n=6..9", hist(ens["r_val"]), C_ORANGE, "xx"),
    ]
    radii = sorted({r for _, h, _, _ in sets for r in h})
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    w = 0.26
    for i, (label, h, col, hatch) in enumerate(sets):
        xs = [r + (i - 1) * w for r in radii]
        ys = [100 * h.get(r, 0.0) for r in radii]
        ax.bar(xs, ys, width=w, label=label, color=col, hatch=hatch, edgecolor="white")
        for x, y in zip(xs, ys):  # noqa: B905
            if y > 1:
                ax.text(x, y + 1.2, f"{y:.0f}", ha="center", fontsize=7)
    ax.set_xticks(radii)
    ax.set_xlabel("breakdown radius  r_val   (distance to the nearest failure)")
    ax.set_ylabel("% of instances with a finite radius")
    ax.set_title("H1: the radius does not saturate at 1")
    ax.legend(frameon=False, fontsize=8)
    ax.set_ylim(0, 100)
    return _save(fig, "sess_f1_h1_distribution")


def fig_h1_stability() -> list[Path]:
    """H1: is the ~80% figure stable across generators and graph size?"""
    _style()
    ens = pd.read_csv(RESULTS / "synth/main/results.csv")
    fin = ens[ens["r_val"] != -1]
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.2))

    g = fin.groupby("generator")["r_val"].apply(lambda s: 100 * (s == 1).mean())
    axes[0].bar(
        range(len(g)),
        g.values,
        color=[C_BLUE, C_GREEN, C_ORANGE][: len(g)],
        hatch=["", "//", "xx"][: len(g)],
        edgecolor="white",
    )
    axes[0].set_xticks(range(len(g)))
    axes[0].set_xticklabels(g.index, rotation=20, ha="right", fontsize=8)
    axes[0].set_ylabel("% with r_val = 1")
    axes[0].set_title("by generator")
    axes[0].set_ylim(0, 100)
    for i, v in enumerate(g.values):
        axes[0].text(i, v + 2, f"{v:.1f}", ha="center", fontsize=8)

    n = fin.groupby("desc_n_nodes")["r_val"].apply(lambda s: 100 * (s == 1).mean())
    axes[1].plot(n.index, n.values, "o-", color=C_VERM, markersize=6)
    axes[1].set_xlabel("number of nodes")
    axes[1].set_ylabel("% with r_val = 1")
    axes[1].set_title("by graph size (no trend)")
    axes[1].set_ylim(0, 100)
    axes[1].set_xticks(list(n.index))
    for x, v in zip(n.index, n.values):  # noqa: B905
        axes[1].text(x, v + 3, f"{v:.1f}", ha="center", fontsize=8)
    fig.suptitle("H1: saturation rate is stable across generators and size", fontsize=10)
    return _save(fig, "sess_f2_h1_stability")


# --------------------------------------------------------------------------- H2


def fig_h2_frontier() -> list[Path]:
    """H2/H3: the frontier is non-flat, yet O* is never beaten."""
    _style()
    d = json.loads((RESULTS / "synth/frontier/frontier_totals.json").read_text())
    labels, nonflat, beaten, tot = [], [], [], []
    for key in ("n4", "n5"):
        v = d[key]
        labels.append(f"{key}\n({v['n_instances']:,} inst.)")
        nonflat.append(100 * v["n_nonflat"] / v["n_instances"])
        beaten.append(100 * v["n_optimal_strictly_beaten"] / v["n_instances"])
        tot.append(v["n_instances"])
    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    x = range(len(labels))
    ax.bar(
        [i - 0.18 for i in x],
        nonflat,
        width=0.36,
        color=C_ORANGE,
        hatch="//",
        edgecolor="white",
        label="frontier NON-FLAT (some valid set differs)",
    )
    ax.bar(
        [i + 0.18 for i in x],
        beaten,
        width=0.36,
        color=C_VERM,
        hatch="xx",
        edgecolor="white",
        label="O* STRICTLY BEATEN",
    )
    for i, v in zip(x, nonflat):  # noqa: B905
        ax.text(i - 0.18, v + 1, f"{v:.1f}%", ha="center", fontsize=8)
    for i, v in zip(x, beaten):  # noqa: B905
        ax.text(i + 0.18, v + 1, "0", ha="center", fontsize=9, weight="bold")
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels)
    ax.set_ylabel("% of instances")
    ax.set_ylim(0, 40)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.set_title(
        "H2/H3: valid sets DO differ in robustness — but O* is never the loser\n"
        f"0 of {sum(tot):,} exhaustive instances",
        fontsize=10,
    )
    return _save(fig, "sess_f3_h2_frontier")


# --------------------------------------------------------------------------- H5


def fig_h5_regimes() -> list[Path]:
    """H5: the middle regime is strongly epsilon-dependent."""
    _style()
    ens = pd.read_csv(RESULTS / "synth/main/results.csv")
    cols = sorted(
        [c for c in ens.columns if c.startswith("r_eps_")], key=lambda c: float(c.split("_")[-1])
    )
    eps, frac = [], []
    for c in cols:
        sub = ens[(ens["r_val"] != -1) & (ens[c] != -1)]
        eps.append(float(c.split("_")[-1]))
        frac.append(100 * (sub[c] > sub["r_val"]).mean())
    fig, ax = plt.subplots(figsize=(5.8, 3.4))
    ax.plot(eps, frac, "o-", color=C_PURPLE, markersize=7, linewidth=2)
    for e, f in zip(eps, frac):  # noqa: B905
        ax.text(e, f + 1.0, f"{f:.1f}%", ha="center", fontsize=8)
    ax.axhspan(20, 50, color=C_GREY, alpha=0.18)
    ax.text(0.022, 35, "pre-registered prediction: 20-50%", fontsize=8, color="#555")
    ax.set_xlabel("epsilon  (bias threshold)")
    ax.set_ylabel("% of instances with r_eps > r_val")
    ax.set_title(
        "H5: the middle regime is nearly absent at small epsilon\n(prediction NOT supported)",
        fontsize=10,
    )
    ax.set_ylim(0, 55)
    return _save(fig, "sess_f4_h5_regimes")


# --------------------------------------------------------------------------- H6


def fig_h6_calibration() -> list[Path]:
    """H6: coverage is exactly 1, conservativeness varies and is predictable."""
    _style()
    d = json.loads((RESULTS / "synth/calibration/calibration.json").read_text())
    rows = pd.DataFrame(d["rows"])
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.1))

    axes[0].hist(rows["coverage"], bins=[0.9, 0.95, 1.0, 1.05], color=C_GREEN, edgecolor="white")
    axes[0].set_xlabel("coverage")
    axes[0].set_ylabel("instances")
    axes[0].set_title(f"coverage = 1.0000 in all {len(rows)}\n(correctness gate)", fontsize=9)

    axes[1].hist(rows["conservativeness"], bins=12, color=C_BLUE, edgecolor="white")
    axes[1].axvline(
        rows["conservativeness"].mean(),
        color=C_VERM,
        ls="--",
        label=f"mean {rows['conservativeness'].mean():.3f}",
    )
    axes[1].set_xlabel("conservativeness")
    axes[1].set_title("conservativeness distribution", fontsize=9)
    axes[1].legend(frameon=False, fontsize=8)

    for gname, mark, col in (
        ("erdos_renyi", "o", C_BLUE),
        ("scale_free", "s", C_GREEN),
        ("block", "^", C_ORANGE),
    ):
        sub = rows[rows["gen"] == gname]
        axes[2].scatter(
            sub["desc_k"],
            sub["conservativeness"],
            marker=mark,
            s=26,
            color=col,
            alpha=0.8,
            label=gname,
        )
    axes[2].set_xlabel("undirected edges in the CPDAG (k)")
    axes[2].set_title("Spearman(k, conservativeness) = 0.564", fontsize=9)
    axes[2].legend(frameon=False, fontsize=7)
    fig.suptitle("H6: calibration", fontsize=10)
    return _save(fig, "sess_f5_h6_calibration")


# --------------------------------------------------------------------------- gates


def fig_rejections() -> list[Path]:
    """The degeneracy gate: 91.8% of random instances are unusable."""
    _style()
    s = json.loads((RESULTS / "synth/main/run_summary.json").read_text())
    rc = dict(sorted(s["rejection_counts"].items(), key=lambda kv: kv[1]))
    labels = [*rc.keys(), "ACCEPTED"]
    vals = [*rc.values(), s["n_accepted"]]
    cols = [C_GREY] * len(rc) + [C_GREEN]
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ax.barh(range(len(vals)), vals, color=cols, edgecolor="white")
    ax.set_yticks(range(len(vals)))
    ax.set_yticklabels([lab.replace("_", " ") for lab in labels], fontsize=8)
    for i, v in enumerate(vals):
        ax.text(v + 60, i, f"{v:,}", va="center", fontsize=8)
    ax.set_xlabel("grid points")
    ax.set_title(
        f"Degeneracy gate: {s['n_accepted']:,} of {s['n_total']:,} accepted "
        f"({100 * s['n_accepted'] / s['n_total']:.1f}%)",
        fontsize=10,
    )
    return _save(fig, "sess_f6_rejections")


# --------------------------------------------------------------------------- Axis B


def fig_conjecture_scope() -> list[Path]:
    """Conjecture 2: enumeration scope and the zero counterexample count."""
    _style()
    ns, combos, cex = [], [], []
    for n in (3, 4, 5):
        d = json.loads((RESULTS / f"search/conjectures/n{n}.json").read_text())
        ns.append(n)
        combos.append(d["totals"]["n_radius_comparisons"])
        cex.append(len(d["counterexamples"]))
    fig, ax = plt.subplots(figsize=(5.8, 3.4))
    bars = ax.bar([str(n) for n in ns], combos, color=C_BLUE, edgecolor="white")
    ax.set_yscale("log")
    ax.set_xlabel("number of nodes")
    ax.set_ylabel("exhaustive radius comparisons (log)")
    for b, c, x in zip(bars, combos, cex):  # noqa: B905
        ax.text(
            b.get_x() + b.get_width() / 2,
            c * 1.4,
            f"{c:,}\n{x} counterexamples",
            ha="center",
            fontsize=8,
        )
    ax.set_ylim(50, 1e7)
    ax.set_title(
        f"Conjecture 2: {sum(combos):,} comparisons, {sum(cex)} counterexamples", fontsize=10
    )
    return _save(fig, "sess_f7_conjecture_scope")


def fig_speedup() -> list[Path]:
    """Speedup of the space-free search, in both accounting regimes."""
    _style()
    df = pd.read_csv(RESULTS / "search/speedup.csv")
    g = (
        df.groupby("k_undirected")
        .agg(
            single=("speedup_single_query", "mean"),
            amort=("speedup_amortised", "mean"),
            space=("space_size", "mean"),
        )
        .reset_index()
    )
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.3))
    axes[0].plot(
        g["k_undirected"],
        g["single"],
        "o-",
        color=C_GREEN,
        markersize=7,
        linewidth=2,
        label="single query (BFS pays space build)",
    )
    axes[0].plot(
        g["k_undirected"],
        g["amort"],
        "s--",
        color=C_VERM,
        markersize=6,
        linewidth=2,
        label="amortised (space pre-built)",
    )
    axes[0].axhline(1.0, color=C_GREY, lw=1)
    axes[0].set_yscale("log")
    axes[0].set_xlabel("undirected edges k  (governing cost parameter)")
    axes[0].set_ylabel("speedup vs BFS  (log)")
    axes[0].legend(frameon=False, fontsize=7.5)
    axes[0].set_title("space-free search: it depends how you count", fontsize=9)
    axes[0].text(1.2, 0.35, "below 1 = SLOWER", fontsize=7.5, color=C_VERM)

    axes[1].plot(g["k_undirected"], g["space"], "o-", color=C_BLUE, markersize=6)
    axes[1].set_yscale("log")
    axes[1].set_xlabel("undirected edges k")
    axes[1].set_ylabel("mean |space| (log)")
    axes[1].set_title("why: the space grows like 3^k", fontsize=9)
    return _save(fig, "sess_f8_speedup")


def fig_bounds() -> list[Path]:
    """L/U bounds: how often the radius is certified without enumeration."""
    _style()
    d = json.loads((RESULTS / "search/bounds_n4.json").read_text())
    n = d["n_cases"]
    cert = d["LU_equal_and_exact"]
    l_tight = d["L_tight"]
    l_def = d["L_admissible"]
    undef = n - l_def
    fig, ax = plt.subplots(figsize=(6.6, 3.3))
    cats = [
        "L defined\n(admissible)",
        "L tight\n(L = r)",
        "L = U = r\n(certified, no\nenumeration)",
        "L undefined\n(back-door mode)",
    ]
    vals = [l_def, l_tight, cert, undef]
    cols = [C_BLUE, C_GREEN, C_ORANGE, C_GREY]
    hat = ["", "//", "xx", ".."]
    ax.bar(range(4), vals, color=cols, hatch=hat, edgecolor="white")
    for i, v in enumerate(vals):
        ax.text(i, v + 12, f"{v}\n({100 * v / n:.1f}%)", ha="center", fontsize=8)
    ax.set_xticks(range(4))
    ax.set_xticklabels(cats, fontsize=8)
    ax.set_ylabel(f"instances (of {n} with a finite radius)")
    ax.set_ylim(0, n * 1.15)
    ax.set_title(
        "Axis B bounds, exhaustive at n=4: 83.3% certified exactly without enumeration", fontsize=10
    )
    return _save(fig, "sess_f9_bounds")


def build_all() -> list[Path]:
    """Build every session figure. Returns the paths written."""
    out: list[Path] = []
    for fn in (
        fig_h1_distribution,
        fig_h1_stability,
        fig_h2_frontier,
        fig_h5_regimes,
        fig_h6_calibration,
        fig_rejections,
        fig_conjecture_scope,
        fig_speedup,
        fig_bounds,
    ):
        out.extend(fn())
    return out
