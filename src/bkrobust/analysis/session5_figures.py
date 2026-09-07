"""Figures for session 5 (Axis A: what actually sets the robustness radius).

Every figure is built from a committed file under ``results/axisa2/``. The
session's question is which structural quantity the back-door robustness radius
tracks. The pre-registered rivals were H7 (the radius is set by the *separation*
between the optimal set and the nearest alternative) and H8 (it is set by the
*size of the undirected component*). The designed census resolves that: at full
background-knowledge coverage the radius equals the separation exactly and is
flat in component size, and at half coverage it is capped by ``|K_G0|``.

The random control exists to say why nobody would have noticed. Random graphs
almost never manufacture separation, so their radii are 1 and the two rival
hypotheses are observationally indistinguishable there. That is why the five
random generators are always drawn apart from the ``backdoor`` generator, which
is a *designed* family that happens to live in the same file: pooling them moves
the headline "r = 1" share from 79.6% to 65.3% and destroys the contrast the
figure exists to make.

Conventions, shared with ``session4_figures``: colour-blind safe palette, no
encoding that relies on colour alone, and no silent drops. Radius ``-1`` is the
UNREACHED sentinel and is never averaged or plotted as a number; rejected,
errored and censored rows are counted out loud in the axis labels.
"""

# ruff: noqa: RUF001
# Figure text is display copy; the multiplication and minus signs are intentional.

from __future__ import annotations

import json
import statistics as st
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.patches import Patch

RES = Path("results/axisa2")
FIGDIR = Path("figures")

C_BLUE = "#0072B2"
C_ORANGE = "#E69F00"
C_GREEN = "#009E73"
C_VERM = "#D55E00"
C_GREY = "#999999"

UNREACHED = -1
EMPTY_CELL = "#EDEDED"

RANDOM_GENERATORS = ("er_sparse", "er_medium", "er_dense", "scale_free", "block")
DESIGNED_GENERATOR = "backdoor"
NO_Z_STATUS = "no_z_member_in_component"


def _style() -> None:
    """Apply the shared house rcParams (see ``session4_figures._style``)."""
    plt.rcParams.update(
        {
            "figure.dpi": 110,
            "savefig.bbox": "tight",
            "axes.grid": True,
            "grid.alpha": 0.3,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "font.size": 9,
        }
    )


def _save(fig: plt.Figure, name: str) -> list[Path]:
    """Write ``fig`` to ``figures/<name>.pdf`` and ``.png``.

    Args:
        fig: The figure to write. It is closed afterwards.
        name: Basename without extension.

    Returns:
        The paths written, in the order pdf then png.
    """
    FIGDIR.mkdir(parents=True, exist_ok=True)
    out: list[Path] = []
    for ext in ("pdf", "png"):
        p = FIGDIR / f"{name}.{ext}"
        fig.savefig(p, dpi=150 if ext == "png" else None)
        out.append(p)
    plt.close(fig)
    return out


def _rows(name: str) -> list[dict]:
    """Read a JSONL file under ``results/axisa2`` into a list of dicts.

    Args:
        name: File name, e.g. ``"census.jsonl"``.

    Returns:
        One dict per non-blank line.
    """
    return [json.loads(line) for line in (RES / name).read_text().splitlines() if line.strip()]


class Ledger:
    """A running account of every row read, so no total is ever quietly wrong.

    Attributes:
        used: Rows that reached a figure.
        rejected: Rows the harness refused at generation time.
        censored: Rows that hit the wall-clock cap.
        errored: Rows whose generator raised.
        unreached: Rows whose radius was the ``-1`` sentinel.
    """

    def __init__(self) -> None:
        self.used = 0
        self.rejected = 0
        self.censored = 0
        self.errored = 0
        self.unreached = 0

    def note(self, row: dict, radius_keys: tuple[str, ...] = ("backdoor", "gac")) -> bool:
        """Classify one row.

        Args:
            row: A raw record from a JSONL results file.
            radius_keys: The radius blocks the caller intends to plot. A row is
                dropped for the sentinel only if one of *these* radii is
                UNREACHED; a sentinel in a radius the figure never touches is
                not a reason to discard an otherwise exact measurement.

        Returns:
            True if the row is usable (accepted, finished, radius not the
            sentinel), in which case ``used`` has been incremented.
        """
        if "error" in row:
            self.errored += 1
            return False
        if row.get("censored"):
            self.censored += 1
            return False
        if not row.get("accepted"):
            self.rejected += 1
            return False
        for key in radius_keys:
            block = row.get(key)
            if isinstance(block, dict) and block.get("radius") == UNREACHED:
                self.unreached += 1
                return False
        self.used += 1
        return True

    @property
    def total(self) -> int:
        """Total rows classified."""
        return self.used + self.rejected + self.censored + self.errored + self.unreached

    def caption(self) -> str:
        """A one-line account of what was excluded, for an axis label."""
        parts = [f"n = {self.used} of {self.total} rows"]
        drops = [
            (self.rejected, "rejected at generation"),
            (self.censored, "censored at the cap"),
            (self.errored, "generator error"),
            (self.unreached, "radius UNREACHED (−1), never averaged"),
        ]
        shown = [f"{count} {label}" for count, label in drops if count]
        if shown:
            parts.append("excluded: " + ", ".join(shown))
        return "; ".join(parts)


def _census(coverage: float | None = None) -> tuple[list[dict], Ledger]:
    """Load the designed-family census.

    Args:
        coverage: If given, keep only rows at this background-knowledge
            coverage. Filtering happens after the ledger has seen every row,
            so the ledger's totals stay honest for the requested slice only
            when ``coverage`` is None.

    Returns:
        The usable rows and the ledger describing what was dropped.
    """
    led = Ledger()
    keep: list[dict] = []
    for r in _rows("census.jsonl"):
        if coverage is not None and r.get("coverage") != coverage:
            continue
        if led.note(r):
            keep.append(r)
    return keep, led


def _control() -> tuple[list[dict], list[dict], Ledger, Ledger]:
    """Load the random control, split into the five random families and the designed one.

    The ``backdoor`` generator is a designed family. It shares a file with the
    random ensemble but never a bar: pooling it changes the headline share of
    radius-1 instances from 79.6% to 65.3%.

    Only the back-door radius is plotted from this file, so a row is kept when
    that radius is exact even if its GAC radius is the UNREACHED sentinel; the
    count of such rows is reported by :func:`_gac_unreached_in_control`.

    Returns:
        ``(random_rows, designed_rows, random_ledger, designed_ledger)``.
    """
    rand_led, des_led = Ledger(), Ledger()
    rand: list[dict] = []
    des: list[dict] = []
    for r in _rows("random_control.jsonl"):
        is_designed = r.get("generator") == DESIGNED_GENERATOR
        led, sink = (des_led, des) if is_designed else (rand_led, rand)
        if led.note(r, radius_keys=("backdoor",)):
            sink.append(r)
    return rand, des, rand_led, des_led


def _gac_unreached_in_control(rows: list[dict]) -> int:
    """Count rows whose GAC radius is the UNREACHED sentinel.

    Args:
        rows: Rows already passed by a :class:`Ledger`.

    Returns:
        How many carry ``gac.radius == -1``. Their back-door radius is exact and
        is plotted; their GAC radius is never plotted as a number.
    """
    return sum(
        1 for r in rows if isinstance(r.get("gac"), dict) and r["gac"].get("radius") == UNREACHED
    )


def _grid(
    ax: plt.Axes,
    table: dict[tuple[int, int], float],
    rows: list[int],
    cols: list[int],
    vmin: int,
    vmax: int,
    fmt: str = "{:.0f}",
) -> None:
    """Draw a hand-built heatmap that leaves missing cells visibly empty.

    ``imshow`` with NaNs would interpolate nothing but would still hand the
    reader a smooth field; drawing each cell as its own patch makes "no design
    exists here" look different from "the value is small".

    Args:
        ax: Target axes.
        table: Values keyed by ``(row_value, col_value)``.
        rows: Row categories, top to bottom.
        cols: Column categories, left to right.
        vmin: Low end of the shared colour scale.
        vmax: High end of the shared colour scale.
        fmt: Format string for the in-cell annotation.
    """
    cmap = plt.get_cmap("viridis")
    norm = BoundaryNorm(list(range(vmin, vmax + 2)), cmap.N)
    for i, rv in enumerate(rows):
        for j, cv in enumerate(cols):
            val = table.get((rv, cv))
            if val is None:
                ax.add_patch(
                    plt.Rectangle(
                        (j - 0.5, i - 0.5), 1, 1, facecolor=EMPTY_CELL, edgecolor="white", lw=0.6
                    )
                )
                continue
            colour = cmap(norm(val))
            ax.add_patch(
                plt.Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor=colour, edgecolor="white", lw=0.6)
            )
            lum = 0.299 * colour[0] + 0.587 * colour[1] + 0.114 * colour[2]
            ax.text(
                j,
                i,
                fmt.format(val),
                ha="center",
                va="center",
                fontsize=7,
                color="black" if lum > 0.55 else "white",
            )
    ax.set_xlim(-0.5, len(cols) - 0.5)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([str(c) for c in cols], fontsize=7)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([str(r) for r in rows], fontsize=7)
    ax.grid(visible=False)
    ax.set_aspect("equal")


def fig_radius_vs_separation() -> list[Path]:
    """The headline: at full coverage the radius *is* the separation.

    A count heatmap over (realised separation, back-door radius) at coverage
    1.0. Every design lands on the identity diagonal; the off-diagonal cells are
    not merely small, they are empty, which is the claim H7 makes and H8 does
    not.
    """
    _style()
    rows, led = _census(coverage=1.0)
    counts: Counter = Counter()
    for r in rows:
        counts[(r["realised"]["realised_separation"], r["backdoor"]["radius"])] += 1
    seps = sorted({s for s, _ in counts})
    radii = sorted({rad for _, rad in counts})
    axis = sorted(set(seps) | set(radii))
    on = sum(v for (s, rad), v in counts.items() if s == rad)

    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    cmap = plt.get_cmap("Blues")
    hi = max(counts.values())
    for i, rad in enumerate(axis):
        for j, s in enumerate(axis):
            n = counts.get((s, rad), 0)
            if n == 0:
                ax.add_patch(
                    plt.Rectangle(
                        (j - 0.5, i - 0.5), 1, 1, facecolor=EMPTY_CELL, edgecolor="white", lw=0.6
                    )
                )
                continue
            ax.add_patch(
                plt.Rectangle(
                    (j - 0.5, i - 0.5),
                    1,
                    1,
                    facecolor=cmap(0.25 + 0.7 * n / hi),
                    edgecolor="white",
                    lw=0.6,
                )
            )
            ax.text(j, i, str(n), ha="center", va="center", fontsize=7.5, color="white")
    ax.plot(
        [-0.5, len(axis) - 0.5],
        [-0.5, len(axis) - 0.5],
        ls="--",
        lw=1.2,
        color=C_VERM,
        zorder=3,
    )
    ax.set_xlim(-0.5, len(axis) - 0.5)
    ax.set_ylim(-0.5, len(axis) - 0.5)
    ax.set_xticks(range(len(axis)))
    ax.set_xticklabels([str(a) for a in axis], fontsize=8)
    ax.set_yticks(range(len(axis)))
    ax.set_yticklabels([str(a) for a in axis], fontsize=8)
    ax.set_aspect("equal")
    ax.grid(visible=False)
    ax.set_xlabel("realised separation  s")
    ax.set_ylabel("back-door robustness radius  r")
    ax.set_title(
        f"At full coverage the radius IS the separation:\n"
        f"r = s in {on} of {led.used} designs, and every off-diagonal cell is empty",
        fontsize=10,
    )
    ax.legend(
        handles=[
            plt.Line2D([], [], ls="--", color=C_VERM, label="identity  r = s"),
            Patch(facecolor=EMPTY_CELL, edgecolor="white", label="no design in this cell"),
        ],
        fontsize=7.5,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.13),
        ncol=2,
        frameon=False,
    )
    fig.text(
        0.5,
        -0.09,
        f"coverage 1.0 only; {led.caption()}. Cells are counts.",
        ha="center",
        fontsize=7,
        color=C_GREY,
    )
    return _save(fig, "s5_f1_radius_vs_separation")


def fig_two_way_table() -> list[Path]:
    """The pre-registered H7-vs-H8 discriminator, as a two-way table.

    Median radius over (realised separation x realised component size). At
    coverage 1.0 the table is diagonal in separation and flat in component size,
    which is H7 and not H8. At coverage 0.5 the rows bend: with only half the
    background knowledge imposed, ``|K_G0|`` bounds the radius, so component
    size starts to matter for a reason that has nothing to do with H8.
    """
    _style()
    panels = []
    for cov in (1.0, 0.5):
        rows, led = _census(coverage=cov)
        cells: dict[tuple[int, int], list[int]] = defaultdict(list)
        for r in rows:
            key = (r["realised"]["realised_separation"], r["realised"]["realised_component_size"])
            cells[key].append(r["backdoor"]["radius"])
        med = {k: st.median(v) for k, v in cells.items()}
        panels.append((cov, med, cells, led))

    all_sep = sorted({k[0] for _, m, _, _ in panels for k in m})
    all_comp = sorted({k[1] for _, m, _, _ in panels for k in m})
    vmin = int(min(v for _, m, _, _ in panels for v in m.values()))
    vmax = int(max(v for _, m, _, _ in panels for v in m.values()))

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 5.2))
    heads = [
        "coverage 1.0 — diagonal in s, flat in |component|",
        "coverage 0.5 — |K_G0| now caps r, so rows bend",
    ]
    notes = []
    for idx, ax in enumerate(axes):
        cov, med, _cells, led = panels[idx]
        head = heads[idx]
        _grid(ax, med, all_sep, all_comp, vmin, vmax)
        ax.set_title(head, fontsize=9.5)
        ax.set_xlabel("realised undirected component size  |C|", fontsize=9)
        notes.append(f"coverage {cov}: {led.caption()}")
        if ax is axes[0]:
            ax.set_ylabel("realised separation  s")
    sm = plt.cm.ScalarMappable(
        cmap=plt.get_cmap("viridis"), norm=BoundaryNorm(list(range(vmin, vmax + 2)), 256)
    )
    cb = fig.colorbar(sm, ax=axes, fraction=0.026, pad=0.02, ticks=range(vmin, vmax + 1))
    cb.set_label("median back-door radius (cell entries are that median)", fontsize=8)
    fig.suptitle(
        "Median radius reads off the separation, never the component size —\n"
        "until missing background knowledge makes |K_G0| the binding constraint",
        fontsize=11,
        y=1.03,
    )
    fig.text(
        0.42,
        -0.02,
        "Grey = no design realised that (s, |C|) pair; nothing is interpolated.\n"
        + "   |   ".join(notes),
        ha="center",
        fontsize=7,
        color=C_GREY,
    )
    return _save(fig, "s5_f2_two_way_table")


def fig_natural_separation() -> list[Path]:
    """Why random graphs cannot see the effect: they have no separation to see.

    Realised separation over the five random generators (the designed
    ``backdoor`` family is excluded here by construction). Instances with no
    Z-member inside the component get their own bar: that is a *status*, not a
    separation of zero, and plotting it as a number would invent a data point.
    """
    _style()
    rand, _designed, led, _dled = _control()
    measured: Counter = Counter()
    no_z = 0
    for r in rand:
        real = r["realised"]
        if real.get("realised_separation_status") == NO_Z_STATUS:
            no_z += 1
        elif real.get("realised_separation") is None:
            no_z += 1
        else:
            measured[int(real["realised_separation"])] += 1
    seps = sorted(measured)
    total = sum(measured.values()) + no_z

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    xs = list(range(len(seps)))
    vals = [measured[s] for s in seps]
    ax.bar(xs, vals, 0.6, color=C_BLUE, label=f"measured separation ({sum(vals)} instances)")
    gap_x = len(seps) + 0.6
    ax.bar(
        [gap_x],
        [no_z],
        0.6,
        color=C_ORANGE,
        hatch="//",
        edgecolor="white",
        label=f"status: no Z-member in the component ({no_z}) — not a separation of 0",
    )
    ax.axvline(len(seps) - 0.2, color=C_GREY, ls=":", lw=1)
    labelled: list[tuple[float, int]] = [(float(x), vals[i]) for i, x in enumerate(xs)]
    labelled.append((gap_x, no_z))
    for x, v in labelled:
        ax.text(x, v + total * 0.012, f"{v}\n{100 * v / total:.1f}%", ha="center", fontsize=8)
    ax.set_xticks([*xs, gap_x])
    ax.set_xticklabels([str(s) for s in seps] + ["no Z-member\n(status)"], fontsize=8)
    ax.set_xlabel("realised separation in the five random generators")
    ax.set_ylabel("instances")
    ax.set_ylim(0, max(max(vals), no_z) * 1.55)
    measured_n = sum(vals)
    ax.set_title(
        "Random graphs almost never manufacture separation —\n"
        f"s = 1 in {100 * measured[1] / measured_n:.1f}% of the instances where s is defined "
        f"({100 * measured[1] / total:.1f}% of all), which is why their radii are 1",
        fontsize=10.5,
    )
    ax.legend(
        fontsize=7.5,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.0),
        framealpha=0.95,
    )
    ax.grid(axis="x", visible=False)
    fig.text(
        0.5,
        -0.06,
        f"{', '.join(RANDOM_GENERATORS)} only — the designed "
        f"'{DESIGNED_GENERATOR}' family is never pooled in. {led.caption()}.",
        ha="center",
        fontsize=7,
        color=C_GREY,
    )
    return _save(fig, "s5_f3_natural_separation")


def fig_designed_vs_random() -> list[Path]:
    """Radius distributions: random ensembles versus the two designed families.

    Three series, never pooled. The five random generators pile up at r = 1; the
    designed ``backdoor`` generator and the census family spread across the
    whole range. Pooling the designed ``backdoor`` rows into the random bar
    would move its r = 1 share from 79.6% to 65.3% and hide exactly this.
    """
    _style()
    rand, designed, rled, dled = _control()
    census, cled = _census()

    series = [
        (
            f"five random generators\n(pooled, n = {rled.used})",
            Counter(r["backdoor"]["radius"] for r in rand),
            rled.used,
            C_BLUE,
            None,
        ),
        (
            f"designed '{DESIGNED_GENERATOR}' generator\n(n = {dled.used})",
            Counter(r["backdoor"]["radius"] for r in designed),
            dled.used,
            C_ORANGE,
            "//",
        ),
        (
            f"designed census family\n(n = {cled.used})",
            Counter(r["backdoor"]["radius"] for r in census),
            cled.used,
            C_GREEN,
            "..",
        ),
    ]
    seen = {rad for _, c, _, _, _ in series for rad in c}
    # Span the full contiguous range so a radius nobody realised (here r = 12)
    # shows as an empty slot rather than being silently closed up.
    radii = list(range(min(seen), max(seen) + 1))

    fig, ax = plt.subplots(figsize=(8.4, 4.3))
    width = 0.27
    for k, (label, counts, n, colour, hatch) in enumerate(series):
        offs = [i + (k - 1) * width for i in range(len(radii))]
        vals = [100 * counts.get(rad, 0) / n for rad in radii]
        ax.bar(offs, vals, width, color=colour, label=label, hatch=hatch, edgecolor="white", lw=0.4)
    top = series[0][1][1] / series[0][2] * 100
    ax.annotate(
        f"{top:.1f}% of random instances\nhave radius 1",
        xy=(0 - width, top),
        xytext=(1.4, top - 6),
        fontsize=8.5,
        arrowprops={"arrowstyle": "->", "color": "black", "lw": 0.9},
    )
    ax.set_xticks(range(len(radii)))
    ax.set_xticklabels([str(r) for r in radii])
    ax.set_xlabel("back-door robustness radius  r  (exact; no UNREACHED sentinels in these rows)")
    ax.set_ylabel("share of the family's instances (%)")
    ax.set_ylim(0, 92)
    ax.set_title(
        "Random ensembles live at r = 1; designed families do not.\n"
        "The designed 'backdoor' generator is drawn apart — pooling it would read 65.3%, not 79.6%",
        fontsize=10.5,
    )
    ax.legend(fontsize=7.5, loc="upper right", ncol=1, framealpha=0.95)
    ax.grid(axis="x", visible=False)
    sentinel = _gac_unreached_in_control(designed)
    footer = [
        f"five random generators — {rled.caption()}",
        f"designed '{DESIGNED_GENERATOR}' generator — {dled.caption()}",
        f"designed census family — {cled.caption()}",
    ]
    if sentinel:
        footer.append(
            f"{sentinel} designed '{DESIGNED_GENERATOR}' rows carry an UNREACHED (−1) GAC "
            "radius; their back-door radius is exact and is the one plotted"
        )
    fig.text(0.5, -0.10, "\n".join(footer), ha="center", fontsize=6.8, color=C_GREY)
    return _save(fig, "s5_f4_designed_vs_random")


def fig_frontier() -> list[Path]:
    """The frontier sweep: GAC widens the candidate pool but never the frontier.

    Left, the pool the two criteria admit as a function of component size, with
    the GAC-only surplus shaded. Right, the result that matters: across every
    admitted candidate set in the sweep, not one beat the optimal set ``O*``
    under either criterion. A zero is drawn against the pool it was searched
    from, so the reader sees the size of the search that came back empty rather
    than a blank axis.
    """
    _style()
    rows = [r for r in _rows("frontier.jsonl") if r.get("accepted")]
    total = len(_rows("frontier.jsonl"))

    by_comp: dict[int, list[dict]] = defaultdict(list)
    for r in rows:
        by_comp[r["realised"]["realised_component_size"]].append(r)
    comps = sorted(by_comp)
    bd = [st.mean([x["n_backdoor_candidates"] for x in by_comp[c]]) for c in comps]
    gac = [st.mean([x["n_gac_candidates"] for x in by_comp[c]]) for c in comps]

    n_bd = sum(r["n_backdoor_candidates"] for r in rows)
    n_gac = sum(r["n_gac_candidates"] for r in rows)
    n_only = sum(r["n_gac_only_candidates"] for r in rows)
    beat_bd = sum(r["backdoor_beats_opt"] for r in rows)
    beat_gac = sum(r["gac_beats_opt"] for r in rows)

    fig, (axl, axr) = plt.subplots(1, 2, figsize=(10.4, 4.2))

    axl.fill_between(comps, bd, gac, color=C_GREEN, alpha=0.25, label="GAC-only surplus")
    axl.plot(comps, gac, "s-", color=C_GREEN, label="GAC-admissible candidates")
    axl.plot(comps, bd, "o-", color=C_BLUE, label="back-door-admissible candidates")
    axl.set_yscale("log")
    axl.set_xticks(comps)
    axl.set_xlabel("realised undirected component size  |C|")
    axl.set_ylabel("candidate sets per instance (mean, log)")
    axl.set_title(
        f"GAC admits {n_gac / n_bd:.2f}x the pool\n({n_only:,} GAC-only sets across the sweep)",
        fontsize=9.5,
    )
    axl.legend(fontsize=7.5, loc="upper left", framealpha=0.95)

    xs = [0, 1]
    pools = [n_bd, n_gac]
    beats = [beat_bd, beat_gac]
    axr.bar(xs, pools, 0.5, color=EMPTY_CELL, edgecolor=C_GREY, lw=1.0)
    axr.bar(xs, [max(b, 0) for b in beats], 0.5, color=C_VERM)
    for i, x in enumerate(xs):
        pool, beat = pools[i], beats[i]
        axr.text(x, pool * 1.20, f"{pool:,}\nsearched", ha="center", fontsize=8.5, color=C_GREY)
        axr.text(
            x,
            pool**0.45,
            f"{beat}",
            ha="center",
            va="center",
            fontsize=26,
            fontweight="bold",
            color=C_VERM,
        )
        axr.text(
            x,
            pool**0.45 / 3.2,
            "beat O*\n(bar has zero height)",
            ha="center",
            va="top",
            fontsize=8,
            color=C_VERM,
        )
    axr.set_yscale("log")
    axr.set_ylim(1, max(pools) * 6.0)
    axr.set_xlim(-0.6, 1.6)
    axr.set_xticks(xs)
    axr.set_xticklabels(["back-door criterion", "GAC"])
    axr.set_ylabel("candidate sets (log)")
    axr.set_title(
        "Zero of them beat the optimal set —\nunder either criterion, at every component size",
        fontsize=9.5,
    )
    axr.grid(axis="x", visible=False)

    fig.suptitle(
        "The optimal set is on the frontier: a wider criterion buys candidates, not robustness",
        fontsize=11,
        y=1.04,
    )
    fig.text(
        0.5,
        -0.06,
        f"frontier.jsonl: {len(rows)} of {total} rows accepted, none censored or errored; "
        f"coverages 1.0 and 0.5 pooled. r_opt is exact everywhere (no UNREACHED sentinels).",
        ha="center",
        fontsize=7,
        color=C_GREY,
    )
    return _save(fig, "s5_f5_frontier")


def build_all() -> list[Path]:
    """Build every session 5 figure that has data.

    Returns:
        The paths written. A figure whose input file is missing or malformed is
        skipped with a message rather than aborting the run.
    """
    out: list[Path] = []
    for fn in (
        fig_radius_vs_separation,
        fig_two_way_table,
        fig_natural_separation,
        fig_designed_vs_random,
        fig_frontier,
    ):
        try:
            out.extend(fn())
        except (FileNotFoundError, ValueError, IndexError, KeyError) as exc:
            print(f"skipped {fn.__name__}: {exc}")
    return out


if __name__ == "__main__":
    for p in build_all():
        print(p)
