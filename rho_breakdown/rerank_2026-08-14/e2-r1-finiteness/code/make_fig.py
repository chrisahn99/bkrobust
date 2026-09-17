import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ARMS = ("G-flip", "T-flip", "T-move")
COL = {"G-flip": "#888888", "T-flip": "#c1442e", "T-move": "#1f6fb4"}

data = json.load(open("../results/e2_raw.json"))
matched = [r for r in data if all(a in r["arms"] for a in ARMS)]
M = {a: [] for a in ARMS}
for r in matched:
    for a in ARMS:
        for m in r["arms"][a]["members"]:
            m = dict(m); m["tau"] = r["tau"]; M[a].append(m)
E = {a: [m for m in M[a] if m["n_changed"] > 0] for a in ARMS}
R = {a: np.array([m["d_orient"] / m["n_changed"] for m in E[a] if m.get("consistent")])
     for a in ARMS}

fig, ax = plt.subplots(1, 4, figsize=(15.5, 3.5))

# (a) rate distribution
bins = np.arange(-0.125, 4.5, 0.25)
for a in ARMS:
    ax[0].hist(R[a], bins=bins, density=True, histtype="step", lw=2,
               color=COL[a], label=f"{a} (q90={np.quantile(R[a],.9):.1f})")
ax[0].set_xlabel(r"$\Delta$orient / $n_{changed}$")
ax[0].set_ylabel("density")
ax[0].set_title("(a) orientation change per unit perturbation")
ax[0].legend(fontsize=8, frameon=False)

# (b) R1-sufficiency
vals = []
for a in ARMS:
    c = [m for m in E[a] if m.get("consistent")]
    vals.append(np.mean([m["r1_equals_full"] for m in c]))
ax[1].bar(range(3), vals, color=[COL[a] for a in ARMS])
ax[1].axhline(0.99, ls="--", c="k", lw=1)
ax[1].set_ylim(0.7, 1.02)
ax[1].set_xticks(range(3)); ax[1].set_xticklabels(ARMS, fontsize=8)
ax[1].set_ylabel("P(R1-only closure = R1-R4 closure)")
ax[1].set_title("(b) GATE 1: is $K'$ still R1-closed?")
for i, v in enumerate(vals):
    ax[1].text(i, v + 0.004, f"{v:.4f}", ha="center", fontsize=8)

# (c) raw d_orient by n_changed
w = 0.26
for k, a in enumerate(ARMS):
    xs, ys, es = [], [], []
    for nc in (1, 2, 3):
        v = [m["d_orient"] for m in E[a] if m.get("consistent") and m["n_changed"] == nc]
        if len(v) < 5:
            continue
        xs.append(nc + (k - 1) * w); ys.append(np.median(v))
        es.append(np.quantile(v, .9) - np.median(v))
    ax[2].bar(xs, ys, width=w, color=COL[a], label=a)
    ax[2].errorbar(xs, ys, yerr=[np.zeros(len(es)), es], fmt="none", ecolor="k", lw=1)
ax[2].set_xticks([1, 2, 3]); ax[2].set_xlabel("$n_{changed}$ (statements)")
ax[2].set_ylabel(r"median $\Delta$orient (bar to q90)")
ax[2].set_title("(c) matched perturbation size")
ax[2].legend(fontsize=8, frameon=False)

# (d) damage CCDF
for a in ARMS:
    v = np.sort([abs(m["bias"]) / abs(m["tau"]) for m in E[a] if m.get("ostar_changed")])
    ax[3].plot(v, 1 - np.arange(len(v)) / len(v), color=COL[a], lw=2, label=a)
ax[3].set_xscale("log"); ax[3].set_yscale("log")
ax[3].set_xlabel(r"$|bias| / |\tau|$  given $O^*$ changed")
ax[3].set_ylabel("P(> x)")
ax[3].set_title("(d) damage tail: heavy in EVERY arm")
ax[3].legend(fontsize=8, frameon=False)

fig.tight_layout()
fig.savefig("../results/fig_e2.png", dpi=160)
fig.savefig("../results/fig_e2.pdf")
print("wrote ../results/fig_e2.png / .pdf")
