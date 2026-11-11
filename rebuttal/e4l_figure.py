"""E4L figure: hybrid time and (projected) listing time against class size."""

from __future__ import annotations

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import OUT

D = OUT / "e4l"
POPS = [("P-SYN", "#7f7f7f"), ("P-REL", "#1f77b4"), ("P-E2", "#2ca02c"), ("P-BIG", "#d62728")]


def main():
    rows = [json.loads(l) for l in open(D / "rows.jsonl")]
    ok = [r for r in rows if r.get("status") == "ok"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for ax, key, title in ((axes[0], "hyb_s", "hybrid (measured)"),
                           (axes[1], "list_projected_s", "listing the class (measured, or projected when capped)")):
        for pop, col in POPS:
            rs = [r for r in ok if r["pop"] == pop and r.get(key)]
            x = np.array([r["class_size_log10"] for r in rs])
            y = np.log10(np.maximum([r[key] for r in rs], 1e-5))
            ax.scatter(x, y, s=4, alpha=0.35, color=col, label=pop, rasterized=True)
        for h, lab in ((0, "1 s"), (np.log10(3600), "1 h"), (np.log10(3.15e7), "1 yr")):
            ax.axhline(h, color="k", lw=0.5, ls=":")
            ax.text(36, h + 0.3, lab, fontsize=7, ha="right")
        ax.set_title(title, fontsize=9)
        ax.set_xlabel(r"$\log_{10}$ |class| (exact count)")
    axes[0].set_ylabel(r"$\log_{10}$ seconds")
    axes[0].legend(fontsize=7, markerscale=3, loc="upper left")
    fig.tight_layout()
    fig.savefig(D / "fig_e4l.pdf", dpi=200)
    print("wrote", D / "fig_e4l.pdf")


if __name__ == "__main__":
    main()
