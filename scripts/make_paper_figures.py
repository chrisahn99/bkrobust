"""Build the main-text figures of the paper from the committed results tree.

    .venv/bin/python scripts/make_paper_figures.py

Reads from the repository this script lives in (override with BKROBUST=<path>)
and writes, into figures/paper/:
  fig_running_example.pdf  -- CPDAG, the analyst's G0, and the knowledge-state
                              space in shells around G0 (running example).
  fig_survival_ranking.pdf -- (a) Kendall tau_b between each predictor and the
                              survival AUC per stratum, instance bootstrap CIs;
                              (b) survival AUC by breakdown-radius bucket.
  fig_e2e_speedup.pdf      -- (a) time per query against the number of undirected
                              edges for space construction, search and hybrid on
                              one instance family; (b) search against hybrid per
                              instance (appendix figure).

Every number is recomputed from committed files; nothing is typed in by hand.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
from scipy.stats import kendalltau

BKROBUST = Path(os.environ.get("BKROBUST", Path(__file__).resolve().parent.parent)).resolve()
sys.path.insert(0, str(BKROBUST / "src"))
OUT = BKROBUST / "figures" / "paper"

OKABE = {
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "orange": "#E69F00",
    "green": "#009E73",
    "sky": "#56B4E9",
    "purple": "#CC79A7",
    "grey": "#999999",
}

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 8,
        "axes.titlesize": 8.5,
        "axes.labelsize": 8,
        "legend.fontsize": 7,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "pdf.fonttype": 42,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    }
)


# ---------------------------------------------------------------------------
# Figure 1: the running example
# ---------------------------------------------------------------------------
def fig_running_example() -> Path:
    from bkrobust.demo.figures import _draw_graph, _shell_layout, compute_layout
    from bkrobust.demo.pipeline import run_scenario

    res = run_scenario("B", n_draws=5)
    cpdag, g0, space = res["cpdag"], res["g0"], res["space"]
    shells, covers = res["shells"], res["covers"]
    valid = {g.edge_string(): bool(r["z_valid"]) for g, r in zip(space, res["rows"])}
    wit = res["witnesses"]["r_val"]
    r_val = wit["shell"]

    fig = plt.figure(figsize=(11.0, 3.55))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.15], wspace=0.04)
    pos = compute_layout(cpdag)
    with plt.rc_context({"axes.titlesize": 11}):
        ax_c = fig.add_subplot(gs[0])
        _draw_graph(ax_c, cpdag, pos, title=r"(a) Estimated CPDAG $\widehat{\mathcal{C}}$")
        ax_g = fig.add_subplot(gs[1])
        _draw_graph(ax_g, g0, pos, title=r"(b) Analyst's state $G_0$ after Meek closure")

    ax = fig.add_subplot(gs[2])
    spos = _shell_layout(space, shells)
    max_shell = max(shells.values())
    for s in range(1, max_shell + 1):
        ax.add_patch(
            Circle((0, 0), s, fill=False, ls="-", lw=0.5,
                   color="#9fd3a8" if s < r_val else "#e3e3e3", zorder=0)
        )
    ax.add_patch(Circle((0, 0), r_val - 0.5, color="#e7f5e9", zorder=-1, lw=0))
    for lo, up in covers:
        if lo in spos and up in spos:
            (x1, y1), (x2, y2) = spos[lo], spos[up]
            ax.plot([x1, x2], [y1, y2], color="#e6e6e6", lw=0.45, zorder=1)
    for g, (px, py) in spos.items():
        if g is g0:
            continue
        ok = valid[g.edge_string()]
        ax.scatter(px, py, s=26 if ok else 30, marker="o" if ok else "X",
                   c=OKABE["blue"] if ok else OKABE["vermillion"],
                   edgecolors="white" if ok else "none", linewidths=0.4, zorder=3)
    ax.scatter(0, 0, s=230, marker="*", c=OKABE["purple"], edgecolors="black",
               linewidths=0.8, zorder=5)
    wg = next(g for g in space if g.edge_string() == wit["edge_string"])
    wx, wy = spos[wg]
    ax.scatter(wx, wy, s=120, facecolors="none", edgecolors=OKABE["orange"],
               linewidths=2.0, zorder=6)
    ax.annotate(rf"nearest failure, $r_{{\mathrm{{val}}}}={r_val}$", (wx, wy),
                xytext=(-6, -58), textcoords="offset points", fontsize=9,
                color="#9a5c00", fontweight="bold", ha="center",
                bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.9},
                arrowprops={"arrowstyle": "-", "color": OKABE["orange"], "lw": 0.8})
    handles = [
        Line2D([0], [0], marker="*", ls="", mfc=OKABE["purple"], mec="black", ms=11,
               label=r"$G_0$"),
        Line2D([0], [0], marker="o", ls="", mfc=OKABE["blue"], mec="white", ms=6,
               label=r"$Z$ valid"),
        Line2D([0], [0], marker="X", ls="", mfc=OKABE["vermillion"], mec="none", ms=7,
               label=r"$Z$ invalid"),
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.93, 1.0),
              frameon=False, fontsize=8.5, handletextpad=0.2)
    ax.set_title(rf"(c) The {len(space)} knowledge states, by distance from $G_0$",
                 fontsize=11)
    lim = max_shell + 0.4
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.axis("off")
    path = OUT / "fig_running_example.pdf"
    fig.savefig(path)
    plt.close(fig)
    print(f"running example: |space|={len(space)} r_val={r_val} "
          f"reachable={len(shells)} covers={len(covers)}")
    return path


# ---------------------------------------------------------------------------
# Figure 2: survival ranking on synthetic structure
# ---------------------------------------------------------------------------
PREDICTORS = [
    ("r_val", r"breakdown radius $r_{\mathrm{val}}$", OKABE["blue"], "o"),
    ("shd_truth", "SHD to true DAG (oracle)", OKABE["orange"], "s"),
    ("n_k", r"number of claims $|\mathcal{K}|$", OKABE["grey"], "D"),
]


def _strata(inst: pd.DataFrame) -> list[tuple[str, pd.DataFrame]]:
    out = []
    for nt in (2, 3, 4):
        out.append((f"tiered, {nt} tiers", inst[(inst.arm == "tiered") & (inst.n_tiers == nt)]))
    for cov in (0.5, 1.0):
        for bw in (0.0, 0.1, 0.25):
            s = inst[(inst.arm == "flip") & (inst.coverage == cov) & (inst.base_wrongness == bw)]
            out.append((f"flip, cov {cov:.1f}, bw {bw:.2f}", s))
    return out


def _tau(x: np.ndarray, y: np.ndarray) -> float:
    if np.all(x == x[0]):
        return float("nan")
    return float(kendalltau(x, y).statistic)


def fig_survival_ranking(n_boot: int = 2000, seed: int = 0) -> Path:
    inst = pd.read_csv(BKROBUST / "results/axis_robustness_p6/survival_instances.csv")
    inst = inst.dropna(subset=["AUC_frac"])
    rng = np.random.default_rng(seed)
    rows = []
    for label, s in _strata(inst):
        y = s.AUC_frac.to_numpy()
        idx = rng.integers(0, len(s), size=(n_boot, len(s)))
        for key, *_ in PREDICTORS:
            x = s[key].to_numpy(dtype=float)
            point = _tau(x, y)
            if np.isnan(point):
                rows.append((label, key, np.nan, np.nan, np.nan, len(s)))
                continue
            boots = np.array([_tau(x[i], y[i]) for i in idx])
            lo, hi = np.nanpercentile(boots, [2.5, 97.5])
            rows.append((label, key, point, lo, hi, len(s)))
    tab = pd.DataFrame(rows, columns=["stratum", "predictor", "tau", "lo", "hi", "n"])
    tab.to_csv(OUT / "fig_survival_ranking_tau.csv", index=False)
    print(tab.round(3).to_string())

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.2, 2.9),
                                 gridspec_kw={"width_ratios": [1.45, 1.0], "wspace": 0.32})
    labels = [lab for lab, _ in _strata(inst)]
    ypos = {lab: len(labels) - 1 - i for i, lab in enumerate(labels)}
    ax.axhspan(len(labels) - 3.5, len(labels) - 0.5, color="#eef4fa", zorder=0)
    offs = {"r_val": 0.22, "shd_truth": 0.0, "n_k": -0.22}
    for key, name, col, mk in PREDICTORS:
        sub = tab[tab.predictor == key]
        for _, r in sub.iterrows():
            yy = ypos[r.stratum] + offs[key]
            if np.isnan(r.tau):
                ax.plot(0, yy, marker="x", color=col, ms=4)
                continue
            ax.plot([r.lo, r.hi], [yy, yy], color=col, lw=1.1)
            ax.plot(r.tau, yy, marker=mk, color=col, ms=3.6, ls="")
        ax.plot([], [], marker=mk, color=col, ls="-", lw=1.1, ms=3.6, label=name)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_yticks([ypos[l] for l in labels])
    ax.set_yticklabels(labels)
    ax.set_xlim(-0.75, 1.0)
    ax.set_xlabel(r"Kendall $\tau_b$ with survival AUC (95% CI)")
    ax.set_title("(a) Ranking quality, by stratum")
    ax.legend(loc="lower center", bbox_to_anchor=(0.4, -0.40), ncol=2, frameon=False,
              columnspacing=1.0, handlelength=1.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    t = inst[inst.arm == "tiered"].copy()
    t["b"] = t.r_val.clip(upper=5).astype(int)
    groups = [t[t.b == k].AUC_frac.to_numpy() for k in range(1, 6)]
    bp = bx.boxplot(groups, positions=range(1, 6), widths=0.55, patch_artist=True,
                    showfliers=False, medianprops={"color": "black", "lw": 1.0},
                    whiskerprops={"lw": 0.7}, capprops={"lw": 0.7})
    for p in bp["boxes"]:
        p.set(facecolor="#cfe3f3", edgecolor=OKABE["blue"], lw=0.8)
    for k, g in zip(range(1, 6), groups):
        bx.text(k, 0.03, f"n={len(g)}", ha="center", fontsize=6, color="#555555")
    bx.set_xticks(range(1, 6))
    bx.set_xticklabels(["1", "2", "3", "4", "5+"])
    bx.set_xlabel(r"breakdown radius $r_{\mathrm{val}}$ at $G_0$")
    bx.set_ylabel("survival AUC under corruption")
    bx.set_ylim(0, 1.04)
    bx.set_title("(b) Tiered knowledge: survival by radius")
    for sp in ("top", "right"):
        bx.spines[sp].set_visible(False)
    path = OUT / "fig_survival_ranking.pdf"
    fig.savefig(path)
    plt.close(fig)
    med = [float(np.median(g)) for g in groups]
    print("tiered median AUC by bucket:", np.round(med, 3), [len(g) for g in groups])
    return path


# ---------------------------------------------------------------------------
# Appendix figure: end-to-end cost on one instance family
# ---------------------------------------------------------------------------
E2E_CAP = 120.0
LEGS = [
    ("space", "space_total_crit_s", "space_censored", "build the space", OKABE["blue"]),
    ("search", "search_fast_s", "search_fast_censored", "retraction search", OKABE["vermillion"]),
    ("hybrid", "hybrid_total_s", "hybrid_censored", "hybrid (Alg. 1)", OKABE["green"]),
]


def fig_e2e_speedup() -> Path:
    d = pd.read_csv(BKROBUST / "results/e2e_speedup_gated/e2e_speedup.csv")
    for c in ("space_censored", "search_fast_censored", "hybrid_censored"):
        d[c] = d[c].astype(str).eq("True")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.2, 2.8),
                                 gridspec_kw={"width_ratios": [1.35, 1.0], "wspace": 0.3})
    ks = sorted(d.k_undirected.unique())
    for key, col, cen, label, color in LEGS:
        xs, med, lo, hi = [], [], [], []
        for k in ks:
            sub = d[d.k_undirected == k]
            done = sub[~sub[cen]][col].dropna().to_numpy()
            if sub[cen].any():
                ax.plot(k, E2E_CAP, marker="v", mfc="none", mec=color, ms=4.5, ls="")
            if len(done) >= 3:
                xs.append(k)
                med.append(np.median(done))
                lo.append(np.percentile(done, 25))
                hi.append(np.percentile(done, 75))
            else:
                ax.plot([k] * len(done), done, marker="o", ms=2.2, color=color, ls="", alpha=0.8)
        # break the median line wherever consecutive k are not both summarised
        seg = [0]
        for i in range(1, len(xs)):
            if xs[i] != xs[i - 1] + 1:
                seg.append(i)
        seg.append(len(xs))
        for a, b in zip(seg[:-1], seg[1:]):
            ax.fill_between(xs[a:b], lo[a:b], hi[a:b], color=color, alpha=0.15, lw=0)
            ax.plot(xs[a:b], med[a:b], color=color, marker="o", ms=2.5, lw=1.2,
                    label=label if a == 0 else None)
    ax.axhline(E2E_CAP, color="#555555", lw=0.6, ls=":")
    ax.text(0.5, E2E_CAP * 1.3, "120 s cap", ha="left", fontsize=6.5, color="#555555")
    ax.set_yscale("log")
    ax.set_ylim(1.5e-5, 400)
    ax.set_xlabel(r"undirected edges $k$ in $\widehat{\mathcal{C}}$")
    ax.set_ylabel("seconds per query")
    ax.set_title("(a) Cost against $k$, same 256 instances")
    ax.legend(loc="lower right", frameon=False, fontsize=6.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    both = d[~d.search_fast_censored & ~d.hybrid_censored]
    for method, color, name in (("local_up_fast", OKABE["sky"], "answered by search"),
                                ("e1_ladder", OKABE["orange"], "answered by SAT leg")):
        sub = both[both.hybrid_method == method]
        bx.scatter(sub.search_fast_s, sub.hybrid_total_s, s=7, color=color, alpha=0.8,
                   lw=0, label=name)
    cen = d[d.search_fast_censored]
    bx.scatter([E2E_CAP] * len(cen), cen.hybrid_total_s, marker=">", s=22,
               facecolors="none", edgecolors=OKABE["orange"], lw=0.9,
               label="search capped")
    bx.plot([1e-5, 300], [1e-5, 300], color="#777777", lw=0.6)
    bx.set_xscale("log")
    bx.set_yscale("log")
    bx.set_xlim(3e-5, 300)
    bx.set_ylim(3e-5, 300)
    bx.set_xlabel("retraction search (s)")
    bx.set_ylabel("hybrid (s)")
    bx.set_title("(b) Per instance")
    bx.legend(loc="upper left", frameon=False, fontsize=6.5, handletextpad=0.2)
    for sp in ("top", "right"):
        bx.spines[sp].set_visible(False)
    path = OUT / "fig_e2e_speedup.pdf"
    fig.savefig(path)
    plt.close(fig)
    print("e2e: rows", len(d), "hybrid median/max",
          round(d.hybrid_total_s.median(), 4), round(d.hybrid_total_s.max(), 1))
    return path


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    print(fig_running_example())
    print(fig_survival_ranking())
    print(fig_e2e_speedup())
