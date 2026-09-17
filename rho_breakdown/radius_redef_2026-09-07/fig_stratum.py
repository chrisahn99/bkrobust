#!/usr/bin/env python3
"""A figura do achado: condicionar no estrato informativo troca qual estatistica funciona."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OK = dict(az="#0072B2", vm="#D55E00", vd="#009E73", cz="#4D4D4D", lar="#E69F00")
plt.rcParams.update({"font.size": 9, "font.family": "sans-serif",
                     "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 160})

# medido por stratum.py, n = 1800 sorteios / 381 informativos
lev = ["0", "1", "2", "3", "nunca"]
data = {
 ("r_val", "todos"):        [0.0, 43.9, 2.8, 0.4, 52.9],
 ("r_val", "informativos"): [0.0, 89.5, 7.9, 1.3, 1.3],
 ("nulidade", "todos"):        [47.8, 5.8, 1.9, 0.4, 44.1],
 ("nulidade", "informativos"): [14.2, 22.8, 8.4, 2.1, 52.5],
}
eff = {("r_val","todos"):2.28, ("r_val","informativos"):1.51,
       ("nulidade","todos"):2.66, ("nulidade","informativos"):3.46}

fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.1), sharey=True)
x = np.arange(len(lev))
for ax, stat, col, title in [(axes[0], "r_val", OK["vm"], "$r_{val}$: quebra de validade"),
                             (axes[1], "nulidade", OK["vd"], "raio de nulidade: o casco alcança zero")]:
    ax.bar(x - 0.19, data[(stat, "todos")], width=0.36, color=OK["cz"], alpha=0.45)
    ax.bar(x + 0.19, data[(stat, "informativos")], width=0.36, color=col)
    ax.set_xticks(x); ax.set_xticklabels(lev); ax.set_xlabel("valor da estatística")
    ax.set_title(title, fontsize=9.5, color=OK["cz"])
    ax.annotate(f"valores efetivos\n{eff[(stat,'todos')]:.2f} → "
                f"\\bf{eff[(stat,'informativos')]:.2f}".replace("\\bf",""),
                xy=(0.97, 0.93), xycoords="axes fraction", ha="right", va="top",
                fontsize=8.5, color=col, fontweight="bold")
axes[0].set_ylabel("% dos problemas")
axes[0].annotate("todos os problemas", xy=(3-0.19, data[("r_val","todos")][3]),
                 xytext=(-6, 30), textcoords="offset points", ha="right",
                 fontsize=8.3, color=OK["cz"], fontweight="bold")
axes[0].annotate("só os informativos", xy=(1+0.19, data[("r_val","informativos")][1]),
                 xytext=(10, -6), textcoords="offset points",
                 fontsize=8.3, color=OK["vm"], fontweight="bold")
fig.suptitle("No estrato onde a pergunta faz sentido, as duas estatísticas trocam de lugar",
             y=1.03, fontsize=10)
fig.tight_layout(); fig.savefig("estrato.pdf", bbox_inches="tight")
print("-> estrato.pdf")
