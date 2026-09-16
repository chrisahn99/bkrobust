"""
Publication-quality figures.

Colour: validated with the dataviz palette validator.
  - ordinal blue ramp for rho (a magnitude, so a single hue light->dark):
      #86b6ef, #2a78d6, #104281   -> ALL CHECKS PASS (--ordinal)
  - 2-slot categorical for the two knowledge regimes:
      #2E5EAA, #D1495B            -> ALL CHECKS PASS
  - status palette (reserved roles) for the outcome composition, every segment
    directly labelled so no meaning is carried by colour alone.

Usage: python figures.py ../results
"""
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.22, "grid.linewidth": 0.5,
    "axes.axisbelow": True, "figure.dpi": 150, "axes.labelsize": 8.5,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7.2,
})

RHO_COL = {1: "#86b6ef", 2: "#2a78d6", 3: "#104281"}   # ordinal ramp
ARM_COL = {"generic": "#2E5EAA", "tiered": "#D1495B"}  # categorical
INK, MUTED = "#1a1a19", "#6b6b68"


def ecdf(v):
    v = np.sort(np.asarray(v, float))
    return v, np.arange(1, len(v) + 1) / len(v)


def fig1(res, plotdata, outdir):
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.7))

    # ---------------- (a) ECDF of ignorance-interval width / |tau|, generic K
    ax = axes[0]
    arm = "generic"
    for rho in (1, 2, 3):
        w = np.asarray(plotdata[arm]["widths"][str(rho)])
        x, y = ecdf(w)
        x = np.concatenate([[0], x])
        y = np.concatenate([[0], y])
        ax.step(x, y, where="post", color=RHO_COL[rho], lw=2.0,
                label=f"ρ = {rho}", solid_capstyle="round")
        # direct label at the point mass at zero
        p0 = float((w <= 1e-12).mean())
        ax.plot([0], [p0], "o", ms=5, color=RHO_COL[rho], mec="w", mew=1.2, zorder=5)
    p0 = float((np.asarray(plotdata[arm]["widths"]["1"]) <= 1e-12).mean())
    ax.annotate(f"{p0*100:.0f}% of SCMs:\nwidth exactly 0\n(one wrong edge\nchanges nothing)",
                xy=(0, p0), xytext=(0.42, 0.30), fontsize=7, color=INK, ha="left",
                arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.8,
                                connectionstyle="arc3,rad=-0.25"))
    ax.axvline(1.0, color=MUTED, ls=":", lw=0.9)
    ax.text(1.06, 0.60, "width = |τ|", fontsize=7, color=MUTED)
    ax.set_xlim(-0.06, 2.4)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("ignorance-interval width / |τ|")
    ax.set_ylabel("fraction of SCMs ≤ x")
    ax.set_title("(a) the ignorance interval is finite\nand usually degenerate",
                 loc="left", fontsize=8.8, color=INK)
    ax.legend(frameon=False, loc="lower right", title="knowledge-error radius",
              title_fontsize=7, borderaxespad=0.4)

    # ---------------- (b) breakdown-radius survival
    ax = axes[1]
    for arm in ("generic", "tiered"):
        if arm not in plotdata:
            continue
        for crit, ls, lab in (("sign", "-", "sign of τ reversed"),
                              ("rel50", "--", "|bias| > 50% |τ|")):
            rr = np.asarray(plotdata[arm][f"rstar_{crit}"])
            ys = [1.0] + [float((rr > r).mean()) for r in (1, 2, 3)]
            ax.plot([0, 1, 2, 3], ys, ls, color=ARM_COL[arm], marker="o", ms=4.5,
                    lw=2.0, mec="w", mew=0.9, label=f"{arm} K — {lab}")
    ax.set_ylim(0, 1.03)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xlabel("knowledge-error radius ρ")
    ax.set_ylabel("conclusion intact (fraction of SCMs)")
    ax.set_title("(b) breakdown-radius survival:\nmost conclusions never break",
                 loc="left", fontsize=8.8, color=INK)
    ax.legend(frameon=False, loc="lower left")

    # ---------------- (c) fate of ONE wrong orientation
    ax = axes[2]
    order = [("meek_inconsistent", "Meek-INCONSISTENT\ncaught for free", "#9aa0a6"),
             ("ostar_unchanged", "O* unchanged\nno bias", "#0ca30c"),
             ("not_amenable", "not amenable\nLOUD failure", "#fab219"),
             ("changed_invalid", "O* invalid\nSILENT bias", "#d03b3b")]
    arms = [a for a in ("generic", "tiered") if a in res["arms"]]
    bottom = np.zeros(len(arms))
    xs = np.arange(len(arms))
    for key, lab, col in order:
        vals = np.array([res["arms"][a]["rho1_buckets_frac"].get(key, 0.0) for a in arms])
        ax.bar(xs, vals, bottom=bottom, color=col, width=0.5, label=lab,
               edgecolor="#fcfcfb", lw=1.4)
        for i, (v, b) in enumerate(zip(vals, bottom)):
            if v > 0.05:
                ax.text(i, b + v / 2, f"{v*100:.0f}%", ha="center", va="center",
                        fontsize=8, color="w", fontweight="bold")
        bottom += vals
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{a} K" for a in arms])
    ax.set_ylim(0, 1)
    ax.set_ylabel("fraction of ρ = 1 perturbations")
    ax.set_title("(c) fate of ONE wrong orientation:\nonly the red band is dangerous",
                 loc="left", fontsize=8.8, color=INK)
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5))
    ax.grid(False)
    ax.text(0.5, -0.26, "the 5th possible outcome — O* changes but stays valid —\n"
                        "was never observed (0 / 754 O*-changing members)",
            transform=ax.transAxes, ha="center", fontsize=6.8, color=MUTED)

    fig.suptitle("How wrong would the expert have to be?  A breakdown radius for background "
                 "knowledge in causal adjustment",
                 fontsize=10.5, x=0.006, ha="left", color=INK, y=0.975)
    fig.text(0.006, 0.895,
             f"{res['n_scm']} random linear iSCMs (p = 5–8, arXiv:2406.11601)  ·  exact "
             f"population estimands, no Monte-Carlo error  ·  O* per Henckel/Perković/Maathuis "
             f"(arXiv:1907.02435)",
             fontsize=7.4, color=MUTED, ha="left")
    fig.subplots_adjust(left=0.055, right=0.855, top=0.775, bottom=0.215, wspace=0.42)
    fig.savefig(f"{outdir}/fig1_breakdown_radius.pdf")
    fig.savefig(f"{outdir}/fig1_breakdown_radius.png")
    print("wrote fig1")


def fig2(res, raw, outdir):
    fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.15))

    # (a) sortability diagnostics
    ax = axes[0]
    keys = [("varsort_naive", "var-sort\nnaive SEM"), ("varsort_iscm", "var-sort\niSCM"),
            ("r2sort_naive", "R²-sort\nnaive SEM"), ("r2sort_iscm", "R²-sort\niSCM")]
    vals = [np.array([r["diag"][k] for r in raw]) for k, _ in keys]
    bp = ax.boxplot(vals, showfliers=False, patch_artist=True, widths=0.5,
                    medianprops=dict(color="w", lw=1.4))
    for b, (k, _) in zip(bp["boxes"], keys):
        b.set(facecolor="#9aa0a6" if "naive" in k else "#2E5EAA", alpha=0.85, lw=0)
    for w in bp["whiskers"] + bp["caps"]:
        w.set(color=MUTED, lw=0.9)
    ax.axhline(0.5, color=INK, ls=":", lw=1.0)
    ax.text(4.45, 0.52, "0.5 = no leak", fontsize=7, ha="right", color=INK)
    ax.set_xticklabels([l for _, l in keys], fontsize=7.2)
    ax.set_ylabel("sortability")
    ax.set_ylim(0, 1.05)
    ax.set_title("(a) the benchmark does not leak the causal order\n"
                 "(diagnostics of Reisach et al., arXiv:2303.18211)",
                 loc="left", fontsize=8.5, color=INK)

    # (b) bias magnitude at rho=1, conditional on O* changing
    ax = axes[1]
    for arm in ("generic", "tiered"):
        b = [abs(m["bias"]) / abs(r["tau"]) for r in raw if arm in r["arms"]
             for m in r["arms"][arm]["members"]
             if m["rho"] == 1 and m.get("consistent") and m.get("amenable")
             and m.get("ostar_changed")]
        if not b:
            continue
        x, y = ecdf(b)
        ax.step(x, 1 - y, where="post", color=ARM_COL[arm], lw=2.0,
                label=f"{arm} K  (n = {len(b)})")
    ax.axvline(0.5, color=MUTED, ls=":", lw=0.9)
    ax.text(0.53, 0.9, "50% of |τ|", fontsize=7, color=MUTED)
    ax.set_xscale("log")
    ax.set_xlabel("|bias| / |τ|   (given that O* moved at all)")
    ax.set_ylabel("fraction of ρ=1 members exceeding x")
    ax.set_title("(b) when the damage lands it is LARGE:\nthe theorem cannot be a continuity bound",
                 loc="left", fontsize=8.5, color=INK)
    ax.legend(frameon=False, loc="lower left")
    ax.set_ylim(0, 1.02)

    fig.tight_layout()
    fig.savefig(f"{outdir}/fig2_diagnostics_and_bias.pdf")
    fig.savefig(f"{outdir}/fig2_diagnostics_and_bias.png")
    print("wrote fig2")


if __name__ == "__main__":
    outdir = sys.argv[1] if len(sys.argv) > 1 else "../results"
    res = json.load(open(f"{outdir}/summary.json"))
    plotdata = json.load(open(f"{outdir}/plotdata.json"))
    raw = json.load(open(f"{outdir}/linear_raw.json"))
    fig1(res, plotdata, outdir)
    fig2(res, raw, outdir)
