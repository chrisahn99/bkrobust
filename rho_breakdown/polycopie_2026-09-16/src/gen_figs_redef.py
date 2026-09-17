#!/usr/bin/env python3
"""Computed figures for chapter 4 (the radius redefinition of 7 to 8 September).

English port of the figure script of the 8 September document on the radius redefinition. It
reads the same two result files and plots the same numbers; only the words changed, plus the
placement of two labels that sat on top of bars. No number is typed by hand: every plotted value
is recorded, with its file and key, in figs/numbers_gen_figs_redef.json.

Sources (read only), both in ~/mugango/output/2026-09-08_apostila-raio-redefinicao/:
  resultados.json   written by roda.py there (exp1 = r_val distribution, exp2 = coverage and
                    width with 0/1/2 reversed claims, exp4 = the informative stratum)
  exemplo.json      written by exemplo.py there (the five-node worked example, seed 4)
The JSON keys of those files are in Portuguese and are used as they are.

Outputs (figs/): redef_degeneracy.pdf, redef_direction.pdf, redef_hull_example.pdf,
redef_stratum.pdf. The two other figures of the original script (the coverage-width frontier and
the locality control) are not needed here and are not drawn.

Usage:  python3 gen_figs_redef.py      (about two seconds)
"""
import json
import pathlib

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
HOME = pathlib.Path.home()
SRC_DIR = HOME / "mugango/output/2026-09-08_apostila-raio-redefinicao"
RES_PATH = SRC_DIR / "resultados.json"
EX_PATH = SRC_DIR / "exemplo.json"
FIGS = HERE / "figs"
FIGS.mkdir(exist_ok=True)

R = json.loads(RES_PATH.read_text())
E = json.loads(EX_PATH.read_text())
OK = dict(az="#0072B2", lar="#E69F00", vd="#009E73", vm="#D55E00", cz="#4D4D4D", rx="#CC79A7")
plt.rcParams.update({"font.size": 9, "font.family": "sans-serif",
                     "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 160})

NUM: dict = {}


def src(p: pathlib.Path) -> str:
    return "~/" + str(p.relative_to(HOME))


def rec(fig, label, value, source, key):
    NUM.setdefault(fig, []).append(
        {"label": label, "value": value, "source": src(source), "key": key})
    return value


# --- Fig 1: the nearly binary distribution of r_val ----------------------------------------
def fig_degeneracy():
    F = "redef_degeneracy"
    e1 = R["exp1"]
    d = e1["dist"]
    n = rec(F, "problems (exp1)", e1["n"], RES_PATH, "exp1.n")
    keys = ["1", "2", "3", "UNREACHED"]
    vals = [rec(F, f"r_val = {k}: problems", d.get(k, 0), RES_PATH,
                f"exp1.dist.{k}" + ("" if k in d else " (missing key = 0)")) for k in keys]
    fb = rec(F, "fraction of mass at the two extremes (1 or not reached)", e1["frac_binary"],
             RES_PATH, "exp1.frac_binary")
    fig, ax = plt.subplots(figsize=(5.2, 2.9))
    cols = [OK["vm"], OK["cz"], OK["cz"], OK["vm"]]
    ax.bar(range(4), vals, color=cols, width=0.62)
    for i, v in enumerate(vals):
        ax.annotate(f"{v}\n{100*v/n:.1f}%", xy=(i, v), xytext=(0, 4),
                    textcoords="offset points", ha="center", fontsize=8.5,
                    color=cols[i], fontweight="bold" if v > 50 else "normal")
    ax.set_xticks(range(4))
    ax.set_xticklabels(["1", "2", "3", "not reached"])
    ax.set_xlabel("$r_{val}$")
    ax.set_ylabel("problems")
    ax.set_ylim(0, max(vals) * 1.28)
    ax.annotate(f"{100*fb:.2f}% of the mass\nat the two extremes",
                xy=(1.5, max(vals) * 0.72), ha="center", fontsize=9,
                color=OK["vm"], fontweight="bold")
    ax.set_title("The statistic is not biased, it is nearly binary", fontsize=9.5, color=OK["cz"])
    fig.tight_layout()
    fig.savefig(FIGS / f"{F}.pdf", bbox_inches="tight")
    plt.close(fig)


# --- Fig 3: how much of the hedge is a dispute about direction ------------------------------
def fig_direction():
    F = "redef_direction"
    arms = ["0", "1", "2"]
    disp, zero, shar = [], [], []
    for a in arms:
        dd = R["exp2"][a]["direction"]
        rec(F, f"{a} reversed claims: informative problems", dd["of"], RES_PATH,
            f"exp2.{a}.direction.of")
        rec(F, f"{a} reversed claims: class disagrees whether X causes Y (count)", dd["disputed"],
            RES_PATH, f"exp2.{a}.direction.disputed")
        disp.append(100 * rec(F, f"{a} reversed claims: class disagrees whether X causes Y "
                                 "(fraction)", dd["frac"], RES_PATH, f"exp2.{a}.direction.frac"))
        zero.append(100 * rec(F, f"{a} reversed claims: width among worlds where X causes Y is 0 "
                                 "(fraction)", dd["frac_zero"], RES_PATH,
                              f"exp2.{a}.direction.frac_zero"))
        shar.append(100 * rec(F, f"{a} reversed claims: that width as fraction of the blanket "
                                 "(mean)", dd["width_share"], RES_PATH,
                              f"exp2.{a}.direction.width_share"))
    fig, ax = plt.subplots(figsize=(5.6, 2.9))
    xx = np.arange(3)
    ax.bar(xx - 0.26, disp, width=0.24, color=OK["vm"])
    ax.bar(xx, zero, width=0.24, color=OK["lar"])
    ax.bar(xx + 0.26, shar, width=0.24, color=OK["az"])
    ax.annotate("the class disagrees on\nwhether $X$ causes $Y$", xy=(0 - 0.26 - 0.12, disp[0]),
                xytext=(0, 5), textcoords="offset points", fontsize=8.3, color=OK["vm"],
                fontweight="bold", ha="left", va="bottom")
    # placed above the tallest bar of its group so that it does not sit on the bars
    ax.annotate("width among the worlds where\n$X$ causes $Y$ is 0",
                xy=(1 - 0.12, max(disp[1], zero[1])), xytext=(0, 5),
                textcoords="offset points", fontsize=8.3, color=OK["lar"], fontweight="bold",
                ha="left", va="bottom")
    ax.annotate("that width, as %\nof the blanket", xy=(2 + 0.26 + 0.12, shar[2]),
                xytext=(4, 0), textcoords="offset points", fontsize=8.3, color=OK["az"],
                fontweight="bold", ha="left", va="center")
    ax.set_xticks(xx)
    ax.set_xticklabels(["0 errors", "1 error", "2 errors"])
    ax.set_ylabel("% of informative problems")
    ax.set_ylim(0, 112)
    ax.set_xlim(-0.55, 2.55)
    ax.set_title("Two thirds of the ignorance is not about confounders",
                 fontsize=9.5, color=OK["cz"])
    fig.tight_layout()
    fig.savefig(FIGS / f"{F}.pdf", bbox_inches="tight")
    plt.close(fig)


# --- Fig 4: the hull in the worked example --------------------------------------------------
def fig_hull_example():
    F = "redef_hull_example"
    effects = []
    for i, r in enumerate(E["theta_blanket"]):
        effects.append(rec(F, f"effect in extension {i} of the class", r["efeito"], EX_PATH,
                           f"theta_blanket[{i}].efeito"))
    vals = sorted(set(effects))
    bl = E["casco_blanket"]
    rh1 = E["escada_RH"][1]["casco"]
    g0 = E["casco_G0"]
    tau = rec(F, "true effect", E["verdadeiro"]["efeito_verdadeiro"], EX_PATH,
              "verdadeiro.efeito_verdadeiro")
    err = E["braco_com_erro"]["escada"][0]["casco"]
    for nm, h, key in (("blanket hull", bl, "casco_blanket"),
                       ("retraction hull RH_1", rh1, "escada_RH[1].casco"),
                       ("RH_0 with the correct expert", g0, "casco_G0"),
                       ("RH_0 with one reversed claim", err, "braco_com_erro.escada[0].casco")):
        rec(F, f"{nm}: low", h[0], EX_PATH, key + "[0]")
        rec(F, f"{nm}: high", h[1], EX_PATH, key + "[1]")
    ratio = rec(F, "RH_1 width / blanket width", E["escada_RH"][1]["razao_blanket"], EX_PATH,
                "escada_RH[1].razao_blanket")
    rec(F, "RH_0 with one reversed claim covers the true effect",
        E["braco_com_erro"]["escada"][0]["cobre_o_verdadeiro"], EX_PATH,
        "braco_com_erro.escada[0].cobre_o_verdadeiro")

    fig, ax = plt.subplots(figsize=(6.6, 3.0))
    ax.axvline(0.0, color=OK["vm"], ls=":", lw=1)
    ax.hlines(2.0, bl[0], bl[1], color=OK["cz"], lw=8, alpha=0.32)
    ax.hlines(1.0, rh1[0], rh1[1], color=OK["az"], lw=8, alpha=0.55)
    ax.plot([g0[0]], [0.0], "o", color=OK["vd"], ms=9)
    ax.plot([err[0]], [-1.0], "o", color=OK["vm"], ms=9)
    ax.plot([tau], [-1.9], "*", color=OK["lar"], ms=15, zorder=5)
    ax.vlines(tau, -1.9, 3.0, color=OK["lar"], ls="--", lw=0.9, alpha=0.6, zorder=1)
    for v in vals:
        ax.plot([v], [3.0], "o", color=OK["cz"], ms=6, zorder=4)
        ax.annotate(f"{v:g}".replace("-", "−"), xy=(v, 3.0), xytext=(0, 9),
                    textcoords="offset points", ha="center", fontsize=8, color=OK["cz"])
    ax.annotate(r"$\Theta(\widehat{C})$: the %d distinct effects the class allows" % len(vals),
                xy=(bl[1], 3.0), xytext=(9, -3), textcoords="offset points",
                fontsize=8.6, color=OK["cz"], va="center")
    ax.annotate("blanket hull, width %.3f" % (bl[1] - bl[0]), xy=(bl[1], 2.0),
                xytext=(9, 0), textcoords="offset points", fontsize=8.6, color=OK["cz"],
                va="center")
    ax.annotate("retraction hull $RH_1$, %.0f%% of the blanket" % (100 * ratio),
                xy=(rh1[1], 1.0), xytext=(9, 0), textcoords="offset points", fontsize=8.6,
                color=OK["az"], va="center")
    ax.annotate("$RH_0$ with the correct expert", xy=(g0[0], 0.0),
                xytext=(9, 0), textcoords="offset points", fontsize=8.6, color=OK["vd"],
                va="center")
    ax.annotate("$RH_0$ with one reversed claim: misses", xy=(err[0], -1.0),
                xytext=(9, 0), textcoords="offset points", fontsize=8.6, color=OK["vm"],
                va="center")
    ax.annotate("true effect", xy=(tau, -1.9), xytext=(9, 0), textcoords="offset points",
                fontsize=8.6, color=OK["lar"], va="center")
    ax.annotate("zero", xy=(0.0, -1.9), xytext=(-4, 0), textcoords="offset points",
                fontsize=8, color=OK["vm"], va="center", ha="right")
    ax.set_ylim(-2.6, 3.9)
    ax.set_xlim(-0.45, 2.25)
    ax.set_yticks([])
    ax.set_xlabel("total effect of $X$ on $Y$")
    ax.spines["left"].set_visible(False)
    ax.set_title("A hull is the smallest interval that contains the effects still possible",
                 fontsize=9.5, color=OK["cz"])
    fig.tight_layout()
    fig.savefig(FIGS / f"{F}.pdf", bbox_inches="tight")
    plt.close(fig)


# --- Fig 5: the informative stratum swaps the two statistics ---------------------------------
def fig_stratum():
    F = "redef_stratum"
    e4 = R["exp4"]
    # data keys: "todos" = all problems, "informativos" = informative stratum,
    # "nulidade" = null radius, "nunca" = never
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.0), sharey=True)
    ks = ["0", "1", "2", "3", "nunca"]
    kshow = ["0", "1", "2", "3", "never"]
    stats = (("r_val", "r_val", OK["vm"], -0.19), ("nulidade", "null radius", OK["az"], 0.19))
    for ax, (lab, key) in zip(axes, [("all problems", "todos"),
                                     ("informative only", "informativos")]):
        blk = e4[key]
        n = rec(F, f"{lab}: n", blk["n"], RES_PATH, f"exp4.{key}.n")
        xs = np.arange(len(ks))
        for name, shown, col, off in stats:
            hs = []
            for k, kl in zip(ks, kshow):
                v = blk[name]["dist"].get(k, 0)
                rec(F, f"{lab}, {shown} = {kl}: fraction", v, RES_PATH,
                    f"exp4.{key}.{name}.dist.{k}" + ("" if k in blk[name]["dist"]
                                                     else " (missing key = 0)"))
                hs.append(100 * v)
            ax.bar(xs + off, hs, width=0.36, color=col)
        ax.set_xticks(xs)
        ax.set_xticklabels(kshow)
        ax.set_title(f"{lab}  (n = {n})", fontsize=9.5, color=OK["cz"])
        ax.set_xlabel("radius")
        ev_r = rec(F, f"{lab}: r_val effective values", blk["r_val"]["valores_efetivos"],
                   RES_PATH, f"exp4.{key}.r_val.valores_efetivos")
        ev_n = rec(F, f"{lab}: null radius effective values",
                   blk["nulidade"]["valores_efetivos"], RES_PATH,
                   f"exp4.{key}.nulidade.valores_efetivos")
        ax.annotate("$r_{val}$: %.2f effective values" % ev_r,
                    xy=(1.45, 92), fontsize=8.4, color=OK["vm"], fontweight="bold")
        ax.annotate("null radius: %.2f effective values" % ev_n,
                    xy=(1.45, 83), fontsize=8.4, color=OK["az"], fontweight="bold")
    axes[0].set_ylabel("% of problems")
    axes[0].set_ylim(0, 103)
    # series names on the first panel, each above one of its own bars
    tod = e4["todos"]
    axes[0].annotate("$r_{val}$", xy=(4 - 0.19, 100 * tod["r_val"]["dist"].get("nunca", 0)),
                     xytext=(0, 4), textcoords="offset points", fontsize=9, ha="center",
                     va="bottom", color=OK["vm"], fontweight="bold")
    axes[0].annotate("null radius", xy=(0 + 0.19, 100 * tod["nulidade"]["dist"].get("0", 0)),
                     xytext=(0, 4), textcoords="offset points", fontsize=9, ha="center",
                     va="bottom", color=OK["az"], fontweight="bold")
    fig.suptitle("Conditioning on the informative stratum worsens one and improves the other",
                 y=1.03, fontsize=10)
    fig.tight_layout()
    fig.savefig(FIGS / f"{F}.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    for f in (fig_degeneracy, fig_direction, fig_hull_example, fig_stratum):
        f()
        print(f"  ok {f.__name__}")
    out = {
        "generated_by": "gen_figs_redef.py",
        "note": "Every plotted value, with the result file and the JSON key it was read from. "
                "The JSON keys of the source files are in Portuguese.",
        "figures": {f"figs/{k}.pdf": v for k, v in NUM.items()},
    }
    (FIGS / "numbers_gen_figs_redef.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    for k, v in NUM.items():
        print(f"figs/{k}.pdf: {len(v)} numbers")


if __name__ == "__main__":
    main()
