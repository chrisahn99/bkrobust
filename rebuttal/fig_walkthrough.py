"""Draft walkthrough figure (C2Bd-W1/Q2), drawn only from states E21 computed.

The retraction up-set of G0 for the walkthrough toy, as a Hasse diagram: each
node is a knowledge state drawn as a graph; each upward line is one covering
retraction, labelled with the orientation withdrawn; green = Z valid, red =
Z invalid. Reads results/rebuttal/e21/toys.json.
"""

import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch  # noqa: E402

from common import OUT, PAPER  # noqa: E402

POS = {"U": (0.0, 0.5), "X": (1.0, 0.5), "W": (2.0, 1.0), "Y": (2.0, 0.0)}
EDGES = [("U", "X"), ("X", "W"), ("X", "Y"), ("W", "Y")]


def parse(state):
    d, u = set(), set()
    for tok in state.split():
        if "->" in tok:
            a, b = tok.split("->")
            d.add((a, b))
        else:
            a, b = tok.split("-")
            u.add(frozenset((a, b)))
    return d, u


def draw_state(ax, ox, oy, state, valid, title, s=0.55):
    d, u = parse(state)
    col = "#1b7f3b" if valid else "#b3261e"
    box = plt.Rectangle((ox - 0.25, oy - 0.3), 2.5 * s + 0.5, 1.0 * s + 0.75, fill=False,
                        ec=col, lw=1.6)
    ax.add_patch(box)
    P = {k: (ox + v[0] * s, oy + v[1] * s) for k, v in POS.items()}
    for a, b in EDGES:
        pa, pb = P[a], P[b]
        if (a, b) in d or (b, a) in d:
            t, h = (a, b) if (a, b) in d else (b, a)
            ax.add_patch(FancyArrowPatch(P[t], P[h], arrowstyle="-|>", mutation_scale=10,
                                         shrinkA=7, shrinkB=7, lw=1.2, color="k"))
        elif frozenset((a, b)) in u:
            ax.plot([pa[0], pb[0]], [pa[1], pb[1]], color="0.45", lw=1.2)
    for k, (x, y) in P.items():
        ax.text(x, y, k, ha="center", va="center", fontsize=8,
                bbox=dict(boxstyle="circle,pad=0.15", fc="white", ec="k", lw=0.8))
    ax.text(ox + 1.25 * s, oy + 1.0 * s + 0.55, title, ha="center", fontsize=8, color=col)


def main():
    t = json.loads((OUT / "e21" / "toys.json").read_text())["walkthrough"]
    rows = {r["state"]: r for r in t["states_table"] if r["above_G0"]}
    g0 = t["G0"]
    g1 = "X->W X->Y U-X W-Y"
    fa, fb = "X->W U-X W-Y X-Y", "X->Y U-X W-X W-Y"
    top = "U-X W-X W-Y X-Y"
    for s in (g1, fa, fb, top):
        assert s in rows, s
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    L = {g0: (1.6, 0.0), g1: (1.6, 2.2), fa: (0.0, 4.4), fb: (3.2, 4.4), top: (1.6, 6.6)}
    lab = {g0: "$G_0$, $d=0$: Z valid", g1: "$d=1$: Z valid", fa: "$d=2$: Z fails",
           fb: "$d=2$: Z fails", top: r"$\widehat{C}$, $d=3$: Z fails"}
    for s, (x, y) in L.items():
        draw_state(ax, x, y, s, rows[s]["Z_valid"], lab[s])
    def link(a, b, text):
        (x1, y1), (x2, y2) = L[a], L[b]
        ax.annotate("", xy=(x2 + 0.7, y2 - 0.32), xytext=(x1 + 0.7, y1 + 1.25),
                    arrowprops=dict(arrowstyle="-", color="0.5", lw=1))
        mx = (x1 + x2) / 2
        ax.text(mx + 0.7 + (0.0 if x1 == x2 else (-0.75 if mx < 1.6 else 0.75)), (y1 + y2) / 2 + 0.45, text, fontsize=7.5, color="0.25",
                ha="center")
    link(g0, g1, r"retract $U\!\to\!X$")
    link(g1, fa, r"retract $X\!\to\!Y$")
    link(g1, fb, r"retract $X\!\to\!W$")
    link(fa, top, r"retract $X\!\to\!W$")
    link(fb, top, r"retract $X\!\to\!Y$")
    ax.text(5.0, 0.2, "Retracting $X\\to W$ or $X\\to Y$ first\nis not closed: Meek's R1\n"
            "re-derives both from $U\\to X$.", fontsize=7.5, color="0.25")
    ax.text(-0.3, 8.3, r"$Z=\mathbf{O}(G_0)=\emptyset$ for $(X,Y)$;  $r_{\mathrm{val}}=2$, two nearest failing states",
            fontsize=8.5)
    ax.set_xlim(-0.5, 7.0)
    ax.set_ylim(-0.6, 8.6)
    ax.axis("off")
    out = PAPER / "rebuttal" / "fig_walkthrough_draft.pdf"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(OUT / "e21" / "fig_walkthrough_draft.png", dpi=160, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    sys.exit(main())
