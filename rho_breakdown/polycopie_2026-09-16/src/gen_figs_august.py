#!/usr/bin/env python3
"""Computed figures for chapters 2 and 3 (the August pilot, E1', locality and X3).

English port of the figure script of the 2 September document on the August state. It reads the
same analysis files and plots the same numbers; only the words and the decimal separator changed,
plus two label placements that pointed at the wrong series. No number is typed by hand: every
plotted value is read from a result file and recorded, with its file and key, in
figs/numbers_gen_figs_august.json.

Sources (read only):
  E1'  ~/mugango/output/2026-08-19_e1prime-se-criterion/results/arm1_analysis_<ens>_K4.json
       ~/mugango/output/2026-08-19_e1prime-se-criterion/results/arm2_analysis.json
  X3   ~/mugango/output/2026-08-31_x3-skeleton-perturbation/results/x3_analysis.json
       ~/mugango/output/2026-08-31_x3-skeleton-perturbation/results/worked_example.txt

Outputs (figs/): aug_coverage.pdf, aug_degeneracy.pdf, aug_locality.pdf, x3_hops.pdf,
x3_taxonomy.pdf, aug_worked_example.pdf.

Usage:  python3 gen_figs_august.py      (a few seconds, no network, no GPU)
"""
import json
import pathlib
import re

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
HOME = pathlib.Path.home()
E1 = HOME / "mugango/output/2026-08-19_e1prime-se-criterion/results"
X3 = HOME / "mugango/output/2026-08-31_x3-skeleton-perturbation/results"
FIGS = HERE / "figs"
FIGS.mkdir(exist_ok=True)

# Okabe-Ito, darkened: the same colour roles as the mugango-apostila package.
BLUE, RED, GREEN, ORANGE, PURPLE = "#1F4E79", "#A83232", "#2E6B4F", "#C2701C", "#5B3E8E"
GREY = "#9AA0A6"

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 160, "savefig.bbox": "tight", "legend.frameon": False,
})

NUM: dict = {}


def src(p: pathlib.Path) -> str:
    return "~/" + str(p.relative_to(HOME))


def rec(fig, label, value, source, key):
    NUM.setdefault(fig, []).append(
        {"label": label, "value": value, "source": src(source), "key": key})
    return value


def f3(x):
    return f"{x:.3f}"


def label_at(ax, x, y, text, colour, dx=6, dy=0, **kw):
    """Direct label on the series."""
    ax.annotate(text, xy=(x, y), xytext=(dx, dy), textcoords="offset points",
                color=colour, fontsize=8, fontweight="bold", va="center", **kw)


ARM1_PATH = {e: E1 / f"arm1_analysis_{e}_K4.json" for e in ("original", "licensed", "large")}
ARM2_PATH = E1 / "arm2_analysis.json"
X3_PATH = X3 / "x3_analysis.json"
WORKED_PATH = X3 / "worked_example.txt"
ARM1 = {e: json.loads(p.read_text()) for e, p in ARM1_PATH.items()}
ARM2 = json.loads(ARM2_PATH.read_text())
XX = json.loads(X3_PATH.read_text())


# ------------------------------------------------------------------------- 1. locality (E1')
def fig_locality():
    """E1' locality law: censoring of rho*_any by hop distance to the query."""
    F = "aug_locality"
    ens = [("original", "original\np 5–8"), ("licensed", "licensed\np 5–10"),
           ("large", "large\np 15–25")]
    fig, ax = plt.subplots(figsize=(6.1, 3.5))
    w = 0.26
    xs = np.arange(len(ens))
    colours = [BLUE, ORANGE, RED]
    names = ["dist. 0", "dist. 1", "dist. ≥ 2"]
    last = []
    for j, d in enumerate((0, 1, 2)):
        vals, ns = [], []
        for e, _ in ens:
            loc = ARM1[e]["locality_by_min_hopdist"]["20000"]
            # key 1048576 is the "no path" sentinel; it belongs to dist >= 2
            keys = [k for k in loc if (int(k) == d if d < 2 else int(k) >= 2)]
            n = sum(loc[k]["n"] for k in keys)
            v = sum(loc[k]["censored"] * loc[k]["n"] for k in keys) / n if n else np.nan
            vals.append(v)
            ns.append(n)
            keytxt = (f"locality_by_min_hopdist.20000.{d}" if d < 2 else
                      "locality_by_min_hopdist.20000.{" + ",".join(keys) + "} (n-weighted)")
            if n:
                rec(F, f"{e}, {names[j]}: censored fraction", round(v, 4), ARM1_PATH[e],
                    keytxt + ".censored")
                rec(F, f"{e}, {names[j]}: n", n, ARM1_PATH[e], keytxt + ".n")
        ax.bar(xs + (j - 1) * w, vals, w, color=colours[j], zorder=3)
        for k, (v, n) in enumerate(zip(vals, ns)):
            if np.isnan(v):
                continue
            ax.text(xs[k] + (j - 1) * w, v + 0.022, f3(v),
                    ha="center", fontsize=7.5, color=colours[j], fontweight="bold")
            if v > 0.25:
                ax.text(xs[k] + (j - 1) * w, v / 2, f"n={n}", ha="center",
                        fontsize=6.4, color="white", rotation=90, va="center")
            else:
                ax.text(xs[k] + (j - 1) * w, v + 0.062, f"n={n}", ha="center",
                        fontsize=6.4, color=GREY)
        last.append(vals[-1])
    # direct series labels right of the last group; dist. 1 and dist. >= 2 sit close, so
    # the upper one is lifted just enough to clear the lower one
    y0, y1, y2 = last
    y2 = max(y2, y1 + 0.075)
    for yy, nm, c in ((y0, names[0], BLUE), (y1, names[1], ORANGE), (y2, names[2], RED)):
        ax.text(2.42, yy, nm, color=c, fontsize=8.5, fontweight="bold", va="center")
    ax.set_xticks(xs)
    ax.set_xticklabels([lab for _, lab in ens])
    ax.set_ylabel("fraction with $\\rho^\\star_{any}$ censored")
    ax.set_ylim(0, 1.14)
    ax.axhline(1.0, color=GREY, lw=0.6, ls=":")
    ax.grid(axis="y", color=GREY, alpha=0.25, lw=0.5, zorder=0)
    fig.savefig(FIGS / f"{F}.pdf")
    plt.close(fig)


# ---------------------------------------------------------------------- 2. degeneracy (E1')
def fig_degeneracy():
    """Censoring does not fall with n: a hard floor over five orders of magnitude."""
    F = "aug_degeneracy"
    P = ARM1_PATH["licensed"]
    sw = ARM1["licensed"]["sweep"]["z1.960"]
    ns_i = sorted(int(k) for k in sw if k != "inf")
    se = [rec(F, f"rho*_se censored fraction, n={n}", sw[str(n)]["censored"], P,
              f"sweep.z1.960.{n}.censored") for n in ns_i]
    any_ = [rec(F, f"rho*_any censored fraction, n={n}", sw[str(n)]["ANY"]["censored"], P,
                f"sweep.z1.960.{n}.ANY.censored") for n in ns_i]
    zero = [rec(F, f"inert ball (d == 0) fraction, n={n}",
                sw[str(n)]["ratio_rho1"]["frac_exactly_zero"], P,
                f"sweep.z1.960.{n}.ratio_rho1.frac_exactly_zero") for n in ns_i]
    floor = rec(F, "rho*_se censored fraction, n=inf (floor)", sw["inf"]["censored"], P,
                "sweep.z1.960.inf.censored")

    fig, ax = plt.subplots(figsize=(6.1, 3.4))
    ax.semilogx(ns_i, se, "-o", color=BLUE, ms=3.4, lw=1.6)
    ax.semilogx(ns_i, any_, "-o", color=ORANGE, ms=3.4, lw=1.6)
    ax.semilogx(ns_i, zero, "--", color=GREEN, lw=1.4)
    label_at(ax, ns_i[-1], se[-1], "$\\rho^\\star_{se}$", BLUE)
    label_at(ax, ns_i[-1], any_[-1], "$\\rho^\\star_{any}$", ORANGE)
    label_at(ax, ns_i[-1], zero[-1], "inert ball\n($d\\equiv 0$)", GREEN, dy=-2)
    ax.axhline(floor, color=GREY, lw=0.8, ls=":")
    ax.text(60, floor - 0.055, f"floor at $n=\\infty$: {f3(floor)}", fontsize=7.5, color=GREY)
    ax.set_xlabel("sample size $n$")
    ax.set_ylabel("censored fraction")
    ax.set_ylim(0, 0.78)
    ax.set_xlim(40, 3e6)
    ax.grid(color=GREY, alpha=0.22, lw=0.5)
    fig.subplots_adjust(right=0.86)
    fig.savefig(FIGS / f"{F}.pdf")
    plt.close(fig)


# ------------------------------------------------------------------------ 3. coverage (ARM 2)
def fig_coverage():
    """ARM 2: what a wrong expert costs, and why n does not help."""
    F = "aug_coverage"
    fr = ARM2["frontier"]
    rs = [0, 1, 2]
    cov, cov_rb = [], []
    for r in rs:
        c = fr[f"r{r}_n20000"]["0"]["coverage"]
        for i, part in enumerate(("estimate", "CI low", "CI high")):
            rec(F, f"naive CI coverage, r={r}, n=20000 ({part})", c[i], ARM2_PATH,
                f"frontier.r{r}_n20000.0.coverage[{i}]")
        cov.append(c)
        c = fr[f"r{r}_n20000"]["1"]["coverage"]
        for i, part in enumerate(("estimate", "CI low", "CI high")):
            rec(F, f"robust CI (rho=1) coverage, r={r}, n=20000 ({part})", c[i], ARM2_PATH,
                f"frontier.r{r}_n20000.1.coverage[{i}]")
        cov_rb.append(c)
    ns = [200, 1000, 5000, 20000]
    wr = [rec(F, f"mean width ratio, r=1, rho=1, n={n}", fr[f"r1_n{n}"]["1"]["mean_width_ratio"],
              ARM2_PATH, f"frontier.r1_n{n}.1.mean_width_ratio") for n in ns]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.6, 3.0))
    a1.errorbar(rs, [c[0] for c in cov],
                yerr=[[c[0] - c[1] for c in cov], [c[2] - c[0] for c in cov]],
                fmt="-o", color=RED, ms=4, lw=1.7, capsize=2.5)
    a1.errorbar(rs, [c[0] for c in cov_rb],
                yerr=[[c[0] - c[1] for c in cov_rb], [c[2] - c[0] for c in cov_rb]],
                fmt="-o", color=GREEN, ms=4, lw=1.7, capsize=2.5)
    a1.axhline(0.95, color=GREY, lw=0.8, ls=":")
    a1.text(0.66, 0.962, "0.95 nominal", fontsize=7.2, color=GREY)
    label_at(a1, 1.5, 0.605, "naive CI", RED, dx=0)
    label_at(a1, 0.95, 0.915, "robust CI ($\\rho=1$)", GREEN, dx=0)
    for r, c in zip(rs, cov):
        a1.annotate(f3(c[0]), xy=(r, c[0]),
                    xytext=(13 if r == 0 else 0, 8 if r == 0 else -15),
                    textcoords="offset points", ha="center", fontsize=7.6,
                    color=RED, fontweight="bold")
    a1.set_xticks(rs)
    a1.set_xlabel("reversed claims $r$")
    a1.set_ylabel("coverage of the 95% CI")
    a1.set_ylim(0.52, 1.0)
    a1.grid(color=GREY, alpha=0.22, lw=0.5)

    a2.semilogx(ns, wr, "-o", color=PURPLE, ms=4, lw=1.7)
    for n, w in zip(ns, wr):
        a2.annotate(f"{w:.2f}×", xy=(n, w), xytext=(0, 9),
                    textcoords="offset points",
                    ha="left" if n == ns[0] else ("right" if n == ns[-1] else "center"),
                    fontsize=7.6, color=PURPLE, fontweight="bold")
    a2.axhline(1.0, color=GREY, lw=0.8, ls=":")
    a2.text(210, 1.15, "no ignorance cost", fontsize=7.2, color=GREY)
    label_at(a2, 220, 6.0, "mean width divided\nby the naive width", PURPLE, dx=0)
    a2.set_xlabel("sample size $n$")
    a2.set_ylabel("width ratio, $\\rho=1$")
    a2.set_ylim(0, 8.4)
    a2.grid(color=GREY, alpha=0.22, lw=0.5)
    fig.tight_layout()
    fig.savefig(FIGS / f"{F}.pdf")
    plt.close(fig)


# --------------------------------------------------------------------------------- 4. X3 hops
def fig_x3_hops():
    """X3: silent bias by hop distance. The cut is at two hops, and it is post hoc."""
    F = "x3_hops"
    ens = ["original", "licensed", "large"]
    hops = [("hop0", "hop 0"), ("hop1", "hop 1"), ("hop2p", "hop >= 2")]
    fig, ax = plt.subplots(figsize=(6.1, 3.5))
    xs = np.arange(3)
    colours = [BLUE, ORANGE, RED]
    floor_plot = 4e-5  # a zero rate is drawn at this floor of the log axis
    series = {}
    for j, e in enumerate(ens):
        sk = XX[e]["SK"]
        raw = []
        for h, hl in hops:
            k = rec(F, f"{e}, {hl}: silent-bias count k", sk[h]["k"], X3_PATH, f"{e}.SK.{h}.k")
            n = rec(F, f"{e}, {hl}: amenable perturbations n", sk[h]["n"], X3_PATH,
                    f"{e}.SK.{h}.n")
            r = rec(F, f"{e}, {hl}: silent-bias rate", sk[h]["rate"], X3_PATH, f"{e}.SK.{h}.rate")
            if r < floor_plot:
                rec(F, f"{e}, {hl}: rate drawn at the log-axis floor (label shows k / n)",
                    floor_plot, X3_PATH, "plotting floor, not data")
            raw.append((k, n, r))
        series[e] = [max(r, floor_plot) for _, _, r in raw], raw
    for j, e in enumerate(ens):
        vals, raw = series[e]
        ax.semilogy(xs, vals, "-o", color=colours[j], ms=4.5, lw=1.7)
        dy_series = [0, -1, 8][j]
        name = ax.annotate(e, xy=(xs[-1], vals[-1]), xytext=(6, dy_series),
                           textcoords="offset points", color=colours[j], fontsize=8,
                           fontweight="bold", va="center")
        # the three series almost touch at hop 0 and hop 1: each carries its value labels
        # on a different side of the point
        dx, dy, ha = [(-6, 4, "right"), (6, 4, "left"), (0, -13, "center")][j]
        for x, v, (k, n, r) in zip(xs, vals, raw):
            txt = f"{k} / {n}" if k == 0 else f"{r:.4f}"
            if x == xs[-1]:
                # at hop >= 2 the incoming lines leave no room left of the points, so the value
                # follows the series name on the same line, on the right of its own point
                ax.annotate(txt, xy=(1, 0.5), xycoords=name, xytext=(4, 0),
                            textcoords="offset points", ha="left", va="center",
                            fontsize=7.4, color=colours[j])
            else:
                ax.annotate(txt, xy=(x, v), xytext=(dx, dy), textcoords="offset points",
                            ha=ha, fontsize=6.9, color=colours[j])
    # the pre-registered bar is a constant of the X3 pre-registration, not a result
    ax.axhline(0.02, color=GREEN, lw=1.1, ls="--")
    ax.text(2.38, 0.0235, "pre-registered bar: $r(\\mathrm{hop}\\geq 1) < 0.02$",
            fontsize=7.4, color=GREEN, ha="right")
    ax.set_xticks(xs)
    ax.set_xticklabels(["hop 0\n(touches the query)", "hop 1", "hop $\\geq$ 2"])
    ax.set_ylabel("silent-bias rate (log scale)")
    ax.set_ylim(3e-5, 0.5)
    ax.set_xlim(-0.35, 2.42)
    ax.grid(color=GREY, alpha=0.22, lw=0.5, which="both")
    fig.subplots_adjust(right=0.87)
    fig.savefig(FIGS / f"{F}.pdf")
    plt.close(fig)


# ---------------------------------------------------------------- 5. X3 taxonomy + amenability
def fig_x3_taxonomy():
    """The four outcomes, and the share the CPDAG alone cannot answer."""
    F = "x3_taxonomy"
    ens = ["original", "licensed", "large"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.7, 3.1),
                                 gridspec_kw={"width_ratios": [1.55, 1]})
    keys = [("unchanged", "unchanged", GREY, "#5f6368"),
            ("not_amenable", "loud failure", ORANGE, ORANGE),
            ("silent_bias", "silent bias", RED, RED),
            ("changed_valid", "efficiency only", GREEN, GREEN)]
    base = np.zeros(len(ens))
    xs = np.arange(len(ens))
    mids_last = []
    tot = {e: rec(F, f"{e}: all perturbations (denominator)", XX[e]["n_perturbations"],
                  X3_PATH, f"{e}.n_perturbations") for e in ens}
    for k, lab, col, _ in keys:
        vals = []
        for e in ens:
            cnt = rec(F, f"{e}: outcome {k} (count)", XX[e]["SK"]["all"]["outcomes"].get(k, 0),
                      X3_PATH, f"{e}.SK.all.outcomes.{k}"
                      + ("" if k in XX[e]["SK"]["all"]["outcomes"] else " (missing key = 0)"))
            vals.append(cnt / tot[e])
            rec(F, f"{e}: outcome {k} (fraction of all perturbations)", round(cnt / tot[e], 4),
                X3_PATH, f"derived: {e}.SK.all.outcomes.{k} / {e}.n_perturbations")
        vals = np.array(vals)
        a1.bar(xs, vals, 0.55, bottom=base, color=col, zorder=3)
        for x, v, b in zip(xs, vals, base):
            if v > 0.035:
                a1.text(x, b + v / 2, f"{100*v:.1f}%", ha="center",
                        va="center", fontsize=7.2,
                        color="white" if col != GREY else "#333333", fontweight="bold")
        mids_last.append(base[-1] + vals[-1] / 2)
        base += vals
    # direct labels right of the last column, at the segment midpoints, pushed apart
    # upwards where the thin top segments crowd each other
    ys = list(mids_last)
    for i in range(1, len(ys)):
        ys[i] = max(ys[i], ys[i - 1] + 0.052)
    for (k, lab, _, tcol), yy in zip(keys, ys):
        a1.text(2.32, yy, lab, color=tcol, fontsize=7.8, fontweight="bold", va="center")
    a1.set_xticks(xs)
    a1.set_xticklabels(ens)
    a1.set_ylabel("fraction of perturbations")
    a1.set_ylim(0, 1.08)
    a1.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    a1.set_xlim(-0.5, 3.35)

    am = []
    for e in ens:
        k = rec(F, f"{e}: CPDAG amenable without K (k)", XX[e]["cpdag_amenable"]["k"], X3_PATH,
                f"{e}.cpdag_amenable.k")
        n = rec(F, f"{e}: CPDAG amenable without K (n)", XX[e]["cpdag_amenable"]["n"], X3_PATH,
                f"{e}.cpdag_amenable.n")
        am.append(k / n)
    a2.barh(xs, am, 0.5, color=BLUE, zorder=3)
    for x, v in zip(xs, am):
        a2.text(v + 0.018, x, f"{100*v:.1f}%", va="center",
                fontsize=8, color=BLUE, fontweight="bold")
    a2.set_yticks(xs)
    a2.set_yticklabels(ens)
    a2.set_xlim(0, 0.78)
    a2.set_xlabel("CPDAG amenable WITHOUT $K$")
    a2.grid(axis="x", color=GREY, alpha=0.22, lw=0.5, zorder=0)
    fig.tight_layout()
    fig.savefig(FIGS / f"{F}.pdf")
    plt.close(fig)


# ------------------------------------------------------------------------- 6. worked example
def fig_worked_example():
    """The worked example: r_val and r_eps do not coincide, and that is what gets reported."""
    F = "aug_worked_example"
    txt = WORKED_PATH.read_text()
    rows = re.findall(r"^\s{5,}(\d)\s+\[\s*(-?\d+\.\d+),\s*(-?\d+\.\d+)\]\s+(\d+\.\d+)\s*$",
                      txt, re.M)
    r = [int(a) for a, _, _, _ in rows]
    lo = [float(b) for _, b, _, _ in rows]
    hi = [float(c) for _, _, c, _ in rows]
    wd = [float(d) for _, _, _, d in rows]
    for ri, l, h, w in zip(r, lo, hi, wd):
        rec(F, f"ignorance interval, radius {ri}: low", l, WORKED_PATH,
            "STEP 5 table, column 'ignorance interval' (low)")
        rec(F, f"ignorance interval, radius {ri}: high", h, WORKED_PATH,
            "STEP 5 table, column 'ignorance interval' (high)")
        rec(F, f"ignorance interval, radius {ri}: width", w, WORKED_PATH,
            "STEP 5 table, column 'width'")
    est = rec(F, "estimate adjusting for O*", float(
        re.search(r"estimate, adjusting for O\*\s*:\s*(-?\d+\.\d+)", txt).group(1)),
        WORKED_PATH, "STEP 3 'estimate, adjusting for O*'")
    r_val = rec(F, "r_val", int(re.search(r"r_val = (\d)", txt).group(1)), WORKED_PATH,
                "STEP 4 'r_val = '")
    treat, outc = re.search(r"Query: effect of (\w+) on (\w+)\.", txt).groups()
    # r_eps: the first radius at which the ignorance interval has positive width
    r_eps = rec(F, "r_eps (first radius with positive width)",
                next(ri for ri, w in zip(r, wd) if w > 0), WORKED_PATH,
                "derived: first row of the STEP 5 table with width > 0")
    jump = wd[r.index(r_eps)]
    pct = rec(F, "width at r_eps as % of |estimate|", round(100 * jump / abs(est)), WORKED_PATH,
              "derived: width at r_eps / |estimate|")

    fig, ax = plt.subplots(figsize=(6.1, 3.3))
    ax.fill_between(r, lo, hi, color=PURPLE, alpha=0.17, step="post", zorder=2)
    ax.step(r, lo, where="post", color=PURPLE, lw=1.8, zorder=3)
    ax.step(r, hi, where="post", color=PURPLE, lw=1.8, zorder=3)
    ax.axhline(est, color=BLUE, lw=1.4, ls="-")
    label_at(ax, r[-1], est, f"estimate\n{est:.4f}".replace("-", "−"), BLUE, dx=-70, dy=13)
    ax.axvline(r_val, color=RED, lw=1.3, ls="--")
    ax.text(r_val - 0.06, -0.93,
            f"$r_{{\\mathrm{{val}}}}={r_val}$\nthe set stops being\nprovably valid",
            fontsize=7.4, color=RED, ha="right")
    ax.axvline(r_eps, color=GREEN, lw=1.3, ls="--")
    ax.text(r_eps + 0.08, -0.93,
            f"$r_{{\\varepsilon}}={r_eps}$\nthe number finally moves:\n"
            f"{jump:.4f} = {pct}% of the effect",
            fontsize=7.4, color=GREEN)
    label_at(ax, 4.2, (lo[-1] + hi[-1]) / 2, "ignorance\ninterval", PURPLE, dx=4)
    ax.set_xlabel("radius $r$ (moves: assert or retract one claim)")
    ax.set_ylabel(f"effect of {treat} on {outc}")
    ax.set_xlim(0, 5.3)
    ax.set_ylim(-1.0, -0.28)
    ax.grid(color=GREY, alpha=0.22, lw=0.5)
    fig.subplots_adjust(right=0.84)
    fig.savefig(FIGS / f"{F}.pdf")
    plt.close(fig)


def main():
    for f in (fig_coverage, fig_degeneracy, fig_locality, fig_x3_hops, fig_x3_taxonomy,
              fig_worked_example):
        f()
        print(f"  ok {f.__name__}")
    out = {
        "generated_by": "gen_figs_august.py",
        "note": "Every plotted value, with the result file and the JSON key (or the text line) it "
                "was read from. 'derived' marks a value computed from the listed keys.",
        "figures": {f"figs/{k}.pdf": v for k, v in NUM.items()},
    }
    (FIGS / "numbers_gen_figs_august.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    for k, v in NUM.items():
        print(f"figs/{k}.pdf: {len(v)} numbers")


if __name__ == "__main__":
    main()
