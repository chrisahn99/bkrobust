"""Fig 3: density/size robustness of the O*-locality mechanism."""
import json, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.22, "grid.linewidth": 0.5,
    "axes.axisbelow": True, "figure.dpi": 150,
})
P_COL = {6: "#86b6ef", 8: "#2a78d6", 10: "#104281"}   # ordinal ramp: |V|
INK, MUTED = "#1a1a19", "#6b6b68"

rows = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "../results/summary_sweep.json"))
fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.15))

def series(p, key):
    r = sorted([x for x in rows if x["p"] == p], key=lambda x: x["deg"])
    return [x["deg"] for x in r], [x[key] for x in r], r

# (a) NUMBER 1: fraction making O* invalid
ax = axes[0]
for p in (6, 8, 10):
    d, v, _ = series(p, "NUMBER1_invalid")
    ax.plot(d, v, "-o", color=P_COL[p], lw=2.0, ms=5, mec="w", mew=0.9, label=f"p = {p}")
ax.axhline(1.0, color="#d03b3b", ls="--", lw=1.2)
ax.text(1.55, 0.93, "falsifier: near-universal disruption", fontsize=7, color="#d03b3b")
ax.set_ylim(0, 1.05)
ax.set_xlabel("expected degree (graph density)")
ax.set_ylabel("frac. of consistent ρ=1 errors\nmaking O* INVALID")
ax.set_title("(a) silent-bias rate stays far\nbelow the falsifier everywhere", loc="left", fontsize=8.7, color=INK)
ax.legend(frameon=False, title="nodes", title_fontsize=7, loc="upper left",
          bbox_to_anchor=(0.0, 0.80))

# (b) ignorance width
ax = axes[1]
for p in (6, 8, 10):
    d, v, r = series(p, "width_mean")
    lo = [x["width_ci"][0] for x in r]; hi = [x["width_ci"][1] for x in r]
    ax.fill_between(d, lo, hi, color=P_COL[p], alpha=0.16, lw=0)
    ax.plot(d, v, "-o", color=P_COL[p], lw=2.0, ms=5, mec="w", mew=0.9, label=f"p = {p} (mean)")
    ax.plot(d, [x["width_median"] for x in r], ":", color=P_COL[p], lw=1.6)
ax.set_xlabel("expected degree (graph density)")
ax.set_ylabel("ρ=1 ignorance width / |τ|")
ax.set_title("(b) density erodes it gracefully;\nthe MEDIAN stays 0 in every cell", loc="left", fontsize=8.7, color=INK)
ax.set_ylim(-0.15, 3.6)
ax.plot([], [], ":", color=MUTED, lw=1.6, label="median (0.000 in every cell)")
ax.legend(frameon=False, loc="upper left")

# (c) sign informativeness
ax = axes[2]
for p in (6, 8, 10):
    d, v, r = series(p, "excl0")
    lo = [x["excl0_ci"][0] for x in r]; hi = [x["excl0_ci"][1] for x in r]
    ax.fill_between(d, lo, hi, color=P_COL[p], alpha=0.16, lw=0)
    ax.plot(d, v, "-o", color=P_COL[p], lw=2.0, ms=5, mec="w", mew=0.9, label=f"p = {p}")
ax.set_ylim(0.5, 1.02)
ax.set_xlabel("expected degree (graph density)")
ax.set_ylabel("frac. of SCMs whose ρ=1 interval\nexcludes 0")
ax.set_title("(c) the sign of the effect survives\none expert error ~90–98% of the time", loc="left", fontsize=8.7, color=INK)
ax.legend(frameon=False, loc="lower left")

fig.suptitle("Robustness of the breakdown radius to graph size and density  ·  "
             "150 linear iSCMs per cell, exact population estimands",
             fontsize=9.8, x=0.006, ha="left", color=INK, y=0.985)
fig.subplots_adjust(left=0.085, right=0.985, top=0.745, bottom=0.155, wspace=0.46)
fig.savefig("../results/fig3_density_sweep.pdf")
fig.savefig("../results/fig3_density_sweep.png")
print("wrote fig3")
