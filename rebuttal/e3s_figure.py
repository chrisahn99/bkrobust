"""E3S figure: median time to an exact radius per arm, (a) against r_val on the
dense P-E3 queries, (b) against n on P-SCALE. Unsolved counts as the 300 s cap.
Families share a colour (CP-SAT blue, MILP orange, SAT/MaxSAT green, explicit
black); arms within a family differ by marker and line style."""
import json, sys, collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import OUT
from e3s_analyse import CAP, t_of

D = OUT / "e3s"
STY = {
    "cpsat_ladder": ("CP-SAT ladder (submitted)", "#2166ac", "o", "-"),
    "cpsat_min": ("CP-SAT minimise", "#2166ac", "s", "--"),
    "cpsat_core": ("CP-SAT core", "#2166ac", "^", ":"),
    "milp_scip": ("MILP SCIP", "#e66101", "o", "-"),
    "milp_highs": ("MILP HiGHS", "#e66101", "s", "--"),
    "milp_cbc": ("MILP CBC", "#e66101", "^", ":"),
    "maxsat_rc2": ("MaxSAT RC2", "#1b7837", "D", "-"),
    "sat_cadical": ("SAT ladder CaDiCaL", "#1b7837", "v", "--"),
    "explicit": ("explicit search", "#1b1b1b", "*", "-"),
}
plt.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 7})
rows = [json.loads(l) for l in open(D / "rows.jsonl")]
by = collections.defaultdict(dict)
for r in rows:
    by[r["qid"]][r["arm"]] = r
fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.5), sharey=True)
for ax, pop, key, xl in ((axs[0], "P-E3", lambda a: a["cpsat_ladder"]["r_ref"], r"$r_{\mathrm{val}}$"),
                         (axs[1], "P-SCALE", lambda a: a["cpsat_ladder"]["n_nodes"], r"$n$")):
    groups = collections.defaultdict(list)
    for q, a in by.items():
        if a["cpsat_ladder"]["pop"] == pop:
            groups[key(a)].append(a)
    xs = sorted(groups)
    for arm, (lab, c, m, ls) in STY.items():
        ys = [np.median([t_of(a[arm]) for a in groups[x]]) for x in xs]
        ax.plot(xs, ys, color=c, marker=m, ls=ls, lw=1.2, ms=4 if m != "*" else 6, label=lab,
                markeredgecolor="white", markeredgewidth=0.4)
    ax.axhline(CAP, color="#777777", lw=0.7, ls="-.")
    ax.text(xs[0], CAP * 1.15, "cap 300 s", fontsize=6, color="#555555")
    ax.set_yscale("log")
    ax.set_xlabel(xl)
    ax.set_xticks(xs)
    ax.grid(axis="y", color="#e5e5e5", lw=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    cnt = ", ".join(f"{len(groups[x])}" for x in xs)
    ax.set_title({"P-E3": "(a) dense ER, $n$ = 10..24 (132 queries)",
                  "P-SCALE": "(b) size ladder, deepest query per graph (10 per $n$)"}[pop], fontsize=7)
axs[0].set_ylabel("median time to exact $r_{\\mathrm{val}}$ (s)")
axs[0].set_ylim(5e-4, 900)
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=5, fontsize=6, frameon=False, bbox_to_anchor=(0.5, -0.02))
fig.tight_layout(rect=(0, 0.13, 1, 1))
fig.savefig(sys.argv[1] if len(sys.argv) > 1 else D / "fig_e3s.pdf")
print("ok")
