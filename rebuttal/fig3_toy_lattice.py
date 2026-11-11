"""Figure 3(b): every knowledge state above G0 for the Figure 3(a) toy.

Toy (chosen 10 Nov 2026 so that the drawn radius is exact and above 1):
CPDAG A-W, W-X, W-Y, X-Y, X-A; claims W->A, A->X; Meek adds W->X, W->Y, X->Y;
query (X, Y), Z = O(G0) = {W}; r_val = 2 (full space of 52 states, hybrid,
brute force over retractions; rebuttal/e21_toys.py logic). The 17 states above
G0 are enumerated as closed subsets of K_G0, drawn by number of retractions,
with covering pairs as edges.
"""

from __future__ import annotations

import sys
from itertools import combinations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import MPDAG, apply_orientations, is_closed, is_gac_valid_mpdag, k_of

BLUE, ORANGE, GREY, BAND = "#2166ac", "#e66101", "#999999", "#d1e5f0"


def states():
    cp = MPDAG(list("AWXY"), undirected=[("A", "W"), ("W", "X"), ("W", "Y"), ("X", "Y"), ("X", "A")])
    g0 = apply_orientations(cp, [("W", "A"), ("A", "X")])
    K = frozenset(k_of(cp, g0))
    out = []
    for j in range(len(K) + 1):
        for drop in combinations(sorted(K), j):
            s = frozenset(K - set(drop))
            if is_closed(cp, s):
                g = apply_orientations(cp, sorted(s))
                out.append((frozenset(drop), is_gac_valid_mpdag(g, "X", "Y", frozenset("W"))))
    return K, out


def main(path):
    K, st = states()
    rows = {}
    for drop, ok in st:
        rows.setdefault(len(drop), []).append((drop, ok))
    r = min(len(d) for d, ok in st if not ok)
    assert r == 2 and all(ok for d, ok in st if len(d) < r)
    # order within rows by barycentre of covered states below
    pos = {}
    for d in sorted(rows):
        if d == 0:
            rows[d].sort(key=lambda t: 0)
        else:
            def bary(t):
                below = [pos[b] for b in pos if len(b) == d - 1 and b < t[0]]
                return sum(below) / len(below) if below else 0
            rows[d].sort(key=lambda t: (bary(t), sorted(t[0])))
        n = len(rows[d])
        for i, (drop, ok) in enumerate(rows[d]):
            pos[drop] = (i - (n - 1) / 2) * 1.0
    fig, ax = plt.subplots(figsize=(1.15, 1.70))
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["font.family"] = "STIXGeneral"
    ax.axhspan(-0.5, r - 0.5, color=BAND, zorder=0, lw=0)
    for d in rows:
        for drop, _ in rows[d]:
            for b, _ in rows.get(d - 1, []):
                if b < drop:
                    ax.plot([pos[b], pos[drop]], [d - 1, d], color=GREY, lw=0.6, zorder=1)
    for d in rows:
        for drop, ok in rows[d]:
            x = pos[drop]
            if d == 0:
                ax.scatter([x], [d], marker="*", s=34, color="#1b1b1b", zorder=3)
            elif ok:
                ax.scatter([x], [d], marker="o", s=12, color=BLUE, zorder=3)
            else:
                ax.scatter([x], [d], marker="x", s=12, color=ORANGE, lw=1.1, zorder=3)
    lab = {("A", "X"): r"$-A{\to}X$", ("W", "A"): r"$-W{\to}A$"}
    for drop, ok in rows[1]:
        (e,) = tuple(drop)
        x = pos[drop]
        ax.text(x + (-0.22 if x < 0 else 0.22), 0.86, lab[e], fontsize=5.6, color=BLUE,
                ha="right" if x < 0 else "left", va="center")
    ax.text(pos[frozenset()] + 0.28, 0, r"$G_0$", fontsize=6.2, va="center")
    ax.axhline(r - 0.5, color=ORANGE, ls="--", lw=0.8, zorder=2)
    ax.text(2.55, r - 0.5 - 0.07, r"$r_{\mathrm{val}}=2$", fontsize=6.2, color=ORANGE, ha="right", va="top")
    for d in rows:
        ax.text(-2.85, d, str(d), fontsize=6, ha="center", va="center", color="#444444")
    ax.text(-2.85, 5.75, r"$d$", fontsize=6, ha="center", va="center", color="#444444")
    # legend
    ax.scatter([-2.6], [7.25], marker="*", s=24, color="#1b1b1b", clip_on=False)
    ax.text(-2.3, 7.25, r"$G_0$", fontsize=6, va="center")
    ax.scatter([-0.9], [7.25], marker="o", s=10, color=BLUE, clip_on=False)
    ax.text(-0.65, 7.25, r"$Z$ valid", fontsize=6, va="center")
    ax.scatter([-2.6], [6.55], marker="x", s=10, color=ORANGE, lw=1.0, clip_on=False)
    ax.text(-2.3, 6.55, r"$Z$ invalid", fontsize=6, va="center")
    ax.set_xlim(-3.1, 2.6)
    ax.set_ylim(-0.6, 7.6)
    ax.axis("off")
    fig.subplots_adjust(0, 0, 1, 1)
    fig.savefig(path)
    print("wrote", path, "r_val =", r, "states", len(st))


if __name__ == "__main__":
    main(sys.argv[1])
