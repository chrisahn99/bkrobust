#!/usr/bin/env python3
"""A lei de localidade, medida na moeda que importa, com controle negativo."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
OK = dict(az="#0072B2", vm="#D55E00", vd="#009E73", cz="#4D4D4D", lar="#E69F00")
plt.rcParams.update({"font.size": 9, "font.family": "sans-serif",
                     "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 160})

# medido por local_hedge.py, 600 problemas informativos por braco
P = ["nenhum\nhedge", "só as\nLONGES", "só as\nPRÓXIMAS", "qualquer\numa", "blanket"]
cov1 = [0.685, 0.687, 0.998, 1.000, 1.000]
w1   = [0.014, 0.023, 0.621, 0.630, 1.000]
cov2 = [0.363, 0.370, 0.760, 0.763, 1.000]
w2   = [0.019, 0.034, 0.554, 0.563, 1.000]
cols = [OK["cz"], OK["vm"], OK["vd"], OK["az"], OK["cz"]]

fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.2), sharey=True)
for ax, cov, w, ttl in ((axes[0], cov1, w1, "1 afirmação invertida"),
                        (axes[1], cov2, w2, "2 afirmações invertidas")):
    for i, (c, ww) in enumerate(zip(cov, w)):
        ax.scatter([ww], [c], s=78, color=cols[i], zorder=3, edgecolor="white", lw=1.2)
        ax.annotate(P[i], xy=(ww, c), xytext=(0, -22), textcoords="offset points",
                    ha="center", fontsize=7.8, color=cols[i], linespacing=0.95)
    ax.plot(w, cov, "-", color=OK["cz"], lw=0.7, alpha=0.4, zorder=1)
    ax.axhline(0.95, color=OK["lar"], ls="--", lw=1)
    ax.set_xlabel("largura / blanket"); ax.set_xlim(-0.07, 1.12); ax.set_title(ttl, fontsize=9.5, color=OK["cz"])
    # a seta do controle negativo
    ax.annotate("", xy=(w[1], cov[1]), xytext=(w[0], cov[0]),
                arrowprops=dict(arrowstyle="->", color=OK["vm"], lw=1.6))
axes[0].set_ylabel("cobertura"); axes[0].set_ylim(0.30, 1.06)
axes[0].annotate("retrair as afirmações LONGE da consulta\nmove a cobertura em $+0{,}2$ pontos",
                 xy=(0.06, 0.70), xytext=(0.20, 0.53), fontsize=8.2, color=OK["vm"],
                 fontweight="bold",
                 arrowprops=dict(arrowstyle="->", color=OK["vm"], lw=0.9))
axes[1].annotate("as PRÓXIMAS recuperam\nquase tudo, por 28\\% menos\nenumerações",
                 xy=(0.554, 0.760), xytext=(0.30, 0.90), fontsize=8.2, color=OK["vd"],
                 fontweight="bold",
                 arrowprops=dict(arrowstyle="->", color=OK["vd"], lw=0.9))
fig.suptitle("A lei de localidade na moeda da cobertura, com controle negativo",
             y=1.02, fontsize=10)
fig.tight_layout(); fig.savefig("localidade.pdf", bbox_inches="tight")
print("-> localidade.pdf")
