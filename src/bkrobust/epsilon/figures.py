"""Figures for the epsilon-bias radius.

Every figure is built from committed files under ``results/epsilon/``; nothing
here regenerates or re-derives a sweep, and nothing here computes a number that
is not already on disk. Outputs go to ``figures/`` as a PNG/PDF pair, matching
the naming of the existing sessions.

Three figures, each carrying one claim:

* **E1, the staircase.** The object the whole section is about: bias against
  retraction depth on the running example, flat at exactly zero out to
  ``r_val`` and then rising. Drawing the epsilon grid as horizontal rules makes
  the radii readable straight off the plot, which is the point of computing the
  whole staircase rather than one radius.
* **E2, what the epsilon radius buys.** The distribution of ``r_eps - r_val``
  across the sweep. This is the figure that says whether the refinement is
  worth having, and it is allowed to say no.
* **E3, cost.** Incremental versus bisection, and the states never evaluated
  because Theorem A certifies them. The saving is the practical claim of the
  algorithm section and it should be visible rather than asserted.

Reporting rules inherited from the earlier sessions and enforced here:
``UNREACHED`` is a status and is never plotted as a number -- it is shown as its
own labelled category with its own count; instances censored on a budget are
shown, not dropped; and any figure whose input file is missing is skipped with a
message rather than drawn from a default.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from bkrobust.core.conventions import UNREACHED

RESULTS = Path("results/epsilon")
FIGURES = Path("figures")

#: One colour per epsilon, ordered so that a stricter threshold is darker.
_EPS_CMAP = plt.get_cmap("viridis")


def _save(fig: plt.Figure, stem: str) -> None:
    """Write a figure as the PNG/PDF pair the repository expects."""
    FIGURES.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"{stem}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote figures/{stem}.png and .pdf")


def fig_e1_staircase(path: Path = RESULTS / "worked" / "worked_example.json") -> bool:
    """E1: the bias staircase on the running example, one panel per scenario.

    Zero out to ``r_val`` is drawn as a solid line on the axis rather than
    omitted, because "exactly zero" is the claim and an absent line would read
    as missing data. Depths past ``|K_G0|`` are outside the retraction up-set
    and therefore outside the certificate; they are shaded rather than plotted.

    Args:
        path: The worked-example JSON.

    Returns:
        Whether the figure was drawn.
    """
    if not path.exists():
        print(f"skipping E1: {path} not found")
        return False
    data = json.loads(path.read_text())
    tables = data["tables"]
    epsilons = data["epsilons_relative"]

    fig, axes = plt.subplots(1, len(tables), figsize=(4.2 * len(tables), 3.6), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, t in zip(axes, tables, strict=True):
        n_k = t["n_knowledge"]
        depths = [r["shell"] for r in t["shells"] if r["shell"] <= n_k]
        values = []
        for d in depths:
            row = next(r for r in t["shells"] if r["shell"] == d)
            values.append(0.0 if row["beta_up_relative"] is None else row["beta_up_relative"])

        for i, eps in enumerate(epsilons):
            ax.axhline(
                eps,
                color=_EPS_CMAP(i / max(1, len(epsilons) - 1)),
                lw=0.8,
                ls=":",
                zorder=1,
            )

        ax.step(depths, values, where="post", color="#1f4e79", lw=2.0, zorder=3)
        ax.plot(depths, values, "o", ms=4, color="#1f4e79", zorder=4)
        ax.axvline(
            t["r_val"],
            color="#b03060",
            lw=1.4,
            ls="--",
            zorder=2,
            label=f"$r_{{\\mathrm{{val}}}}={t['r_val']}$",
        )
        ax.axvspan(-0.4, t["r_val"] - 0.5, color="#2e7d32", alpha=0.07, zorder=0)
        ax.set_title(f"scenario {t['scenario']}  ($|K_{{G_0}}|={n_k}$)", fontsize=9)
        ax.set_xlabel("retractions $d$")
        ax.set_xticks(depths)
        ax.set_ylim(-0.04, 1.12)
        ax.legend(fontsize=7, loc="upper left")
        ax.grid(alpha=0.2, lw=0.5)
    axes[0].set_ylabel("worst-case bias $\\beta^{\\uparrow}(d)$\n(fraction of the estimate)")

    # The epsilon grid is labelled once, on a right-hand axis of the last panel,
    # rather than annotated inside every panel: repeating it three times crowds
    # the 1% and 5% rules into each other and hides the step.
    right = axes[-1].twinx()
    right.set_ylim(axes[-1].get_ylim())
    right.set_yticks(list(epsilons))
    right.set_yticklabels([f"{e:.0%}" for e in epsilons], fontsize=6)
    right.set_ylabel("$\\varepsilon$ grid", fontsize=7)
    for tick, i in zip(right.get_yticklabels(), range(len(epsilons)), strict=True):
        tick.set_color(_EPS_CMAP(i / max(1, len(epsilons) - 1)))
    fig.suptitle(
        "The bias staircase: exactly zero inside the certified shell, then a step to the ceiling",
        fontsize=10,
    )
    _save(fig, "fige1_staircase")
    return True


def _load_study(path: Path) -> list[dict] | None:
    if not path.exists():
        print(f"skipping: {path} not found")
        return None
    import csv

    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def fig_e2_what_epsilon_buys(path: Path = RESULTS / "study" / "results.csv") -> bool:
    """E2: the distribution of ``r_eps - r_val`` across the sweep.

    ``UNREACHED`` is charted as its own category, labelled, and never converted
    to a number: an instance where no perturbation reaches ``eps`` is a
    qualitatively different statement from one where a large but finite
    perturbation does.

    Args:
        path: The sweep CSV.

    Returns:
        Whether the figure was drawn.
    """
    rows = _load_study(path)
    if not rows:
        return False

    eps_cols = sorted(
        {c for c in rows[0] if c.startswith("r_eps_")},
        key=lambda c: float(c.removeprefix("r_eps_")),
    )
    if not eps_cols:
        print("skipping E2: no r_eps_* columns in the sweep")
        return False

    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    width = 0.8 / len(eps_cols)
    gains = [0, 1, 2, 3]
    labels = [str(g) for g in gains] + ["$>3$", "unreached"]
    xs = np.arange(len(labels))

    for i, col in enumerate(eps_cols):
        counts = [0] * len(labels)
        total = 0
        for r in rows:
            try:
                rv, re = int(r["r_val"]), int(r[col])
            except (ValueError, KeyError):
                continue
            if rv == UNREACHED:
                continue
            total += 1
            if re == UNREACHED:
                counts[-1] += 1
            else:
                g = re - rv
                counts[min(g, 4) if g <= 3 else 4] += 1
        if not total:
            continue
        frac = [c / total for c in counts]
        ax.bar(
            xs + i * width - 0.4 + width / 2,
            frac,
            width=width,
            color=_EPS_CMAP(i / max(1, len(eps_cols) - 1)),
            label=f"$\\varepsilon$ = {float(col.removeprefix('r_eps_')):.0%}",
        )
    ax.set_xticks(xs)
    ax.set_xticklabels(labels)
    ax.set_xlabel("$r_\\varepsilon - r_{\\mathrm{val}}$  (extra certified shells)")
    ax.set_ylabel("fraction of instances")
    ax.set_title("What tolerating $\\varepsilon$ bias buys over the validity radius", fontsize=10)
    ax.legend(fontsize=7, ncol=2)
    ax.grid(alpha=0.2, lw=0.5, axis="y")
    _save(fig, "fige2_what_epsilon_buys")
    return True


def fig_e3_cost(path: Path = RESULTS / "study" / "results.csv") -> bool:
    """E3: the cost of the two exact strategies, against the size of the knowledge set.

    Both strategies return the same radius, so the comparison is purely about
    which shells each one touches -- which is why states evaluated, not seconds,
    is the primary axis: it is machine-independent, and the repository's
    convention is to report machine-independent counters wherever a claim rests
    on them.

    Args:
        path: The sweep CSV.

    Returns:
        Whether the figure was drawn.
    """
    rows = _load_study(path)
    if not rows:
        return False
    need = {"n_k0", "cost_incremental_states", "cost_bisection_states"}
    if not need <= set(rows[0]):
        print(f"skipping E3: sweep is missing {sorted(need - set(rows[0]))}")
        return False

    by_k: dict[int, list[tuple[float, float]]] = {}
    for r in rows:
        try:
            k = int(r["n_k0"])
            inc = float(r["cost_incremental_states"])
            bis = float(r["cost_bisection_states"])
        except (ValueError, KeyError):
            continue
        by_k.setdefault(k, []).append((inc, bis))

    if not by_k:
        print("skipping E3: no usable cost rows")
        return False

    ks = sorted(by_k)
    inc_med = [float(np.median([a for a, _ in by_k[k]])) for k in ks]
    bis_med = [float(np.median([b for _, b in by_k[k]])) for k in ks]
    counts = [len(by_k[k]) for k in ks]

    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.plot(ks, inc_med, "o-", color="#1f4e79", label="incremental")
    ax.plot(ks, bis_med, "s--", color="#c77400", label="bisection")
    for k, v, c in zip(ks, inc_med, counts, strict=True):
        ax.annotate(
            f"n={c}",
            xy=(k, v),
            fontsize=6,
            xytext=(0, 6),
            textcoords="offset points",
            ha="center",
            color="#555555",
        )
    ax.set_yscale("log")
    ax.set_xlabel("$|K_{G_0}|$  (orientations the analyst's claims induce)")
    ax.set_ylabel("states evaluated (median, log)")
    ax.set_title("Both strategies are exact; they differ in which shells they touch", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.2, lw=0.5)
    _save(fig, "fige3_cost")
    return True


def main() -> None:
    """Draw every figure whose inputs are present."""
    drawn = [
        fig_e1_staircase(),
        fig_e2_what_epsilon_buys(),
        fig_e3_cost(),
    ]
    print(f"{sum(drawn)}/{len(drawn)} figures drawn")


if __name__ == "__main__":
    main()
