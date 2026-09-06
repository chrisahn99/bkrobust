"""Figures for session 6 (Axis A, out of sample: does any of this happen in real DAGs?).

Sessions 4 and 5 were built on graphs we generated. That is the obvious
objection to the whole line of work: the back-door robustness radius is only
interesting when a CPDAG leaves something undirected *and* the optimal
adjustment set sits two or more edges away from the treatment inside that
undirected component, and a designed generator can manufacture both. Session 6
runs the same measurement over 39 published networks (bnlearn's repository plus
the dagitty examples) and asks whether the structure survives contact with
graphs nobody drew for us.

It does, and the figures here are the evidence: the median CPDAG leaves 10.7% of
its edges undirected, chain components reach 85 nodes where the synthetic
maxima were 6 and 12, and separation is at least 2 in 45.4% of the instances
where separation is defined against 14.2% under the Erdős-Rényi control from
session 5.

Conventions, shared with ``session4_figures`` and ``session5_figures``:
colour-blind safe palette, no encoding that relies on colour alone, and no
silent drops. Two sentinels govern this file and neither is ever a number:

* ``separation: null`` (status ``no_z_member_in_component``) means the instance
  has no Z-member inside X's chain component, so separation is undefined. It is
  drawn as its own hatched bar, never as a separation of 0, and never averaged.
* ``radius: -1`` is UNREACHED. Rejected screening rows carry it as filler; among
  admissible instances there are none, which :func:`_admissible` asserts rather
  than assumes.

Everything excluded is counted out loud: the M-bias ADMG (not a DAG, so no
CPDAG), the 115,143 rejected pairs, and the 3 censored per-network runs.
"""

# ruff: noqa: RUF001
# Figure text is display copy; the en dashes and minus signs are intentional.

from __future__ import annotations

import csv
import json
import statistics as st
from collections import Counter
from pathlib import Path
from typing import NamedTuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

RES = Path("results/axisa3")
FIGDIR = Path("figures")

C_BLUE = "#0072B2"
C_ORANGE = "#E69F00"
C_GREEN = "#009E73"
C_VERM = "#D55E00"
C_GREY = "#999999"

UNREACHED = -1
NO_Z_STATUS = "no_z_member_in_component"
MEASURED_STATUS = "measured"

# Session 5's random control: 85.8% of measurable separations were s = 1, so
# 14.2% were s >= 2, with a maximum of 4. See results/axisa3/preregistration.md.
ER_S1_SHARE = 0.858
ER_GE2_PCT = 14.2

# The synthetic ceilings this session is measured against.
S4_MAX_COMPONENT = 6  # session 4's largest chain component, at n = 24
S5_MAX_COMPONENT = 12  # session 5's designed generator, by construction


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


class Network(NamedTuple):
    """One row of ``descriptive_structure.csv``, for a network with a CPDAG.

    Attributes:
        name: Network identifier.
        source_format: ``bif``, ``bnjson`` or ``dagitty``.
        n_nodes: Node count.
        n_edges: Edge count of the underlying DAG.
        n_undirected: CPDAG edges left undirected.
        undirected_fraction: ``n_undirected / n_edges``.
        max_component: Size of the largest chain component, in nodes.
        n_components: Number of chain components of size >= 2.
    """

    name: str
    source_format: str
    n_nodes: int
    n_edges: int
    n_undirected: int
    undirected_fraction: float
    max_component: int
    n_components: int


class Descriptive(NamedTuple):
    """The descriptive table, split into rows that have a CPDAG and rows that do not.

    Attributes:
        networks: Networks whose input was a DAG, so a CPDAG exists.
        non_dag: Names of rows with ``is_dag=False`` and empty CPDAG columns.
    """

    networks: list[Network]
    non_dag: list[str]

    def footer(self) -> str:
        """A one-line account of what the descriptive figures exclude."""
        if not self.non_dag:
            return f"{len(self.networks)} networks, none excluded."
        return (
            f"{len(self.networks)} of {len(self.networks) + len(self.non_dag)} networks; "
            f"excluded: {', '.join(self.non_dag)} "
            "(an ADMG with bidirected edges, not a DAG, so it has no CPDAG and no "
            "undirected fraction to report)"
        )


def _descriptive() -> Descriptive:
    """Read ``descriptive_structure.csv``.

    Rows with ``is_dag=False`` carry empty CPDAG columns. They are separated out
    rather than coerced to zero: "this network has no CPDAG" and "this network's
    CPDAG is fully compelled" are different facts and the figures must not
    conflate them.

    Returns:
        The split table.
    """
    keep: list[Network] = []
    dropped: list[str] = []
    with (RES / "descriptive_structure.csv").open() as fh:
        for row in csv.DictReader(fh):
            if row["is_dag"] != "True" or not row["undirected_fraction"]:
                dropped.append(row["network"])
                continue
            keep.append(
                Network(
                    name=row["network"],
                    source_format=row["source_format"],
                    n_nodes=int(row["n_nodes"]),
                    n_edges=int(row["n_edges"]),
                    n_undirected=int(row["n_undirected_edges"]),
                    undirected_fraction=float(row["undirected_fraction"]),
                    max_component=int(row["max_component_size"]),
                    n_components=int(row["n_components"]),
                )
            )
    return Descriptive(networks=keep, non_dag=dropped)


class Screening(NamedTuple):
    """Everything ``instances.jsonl`` contains, kept apart rather than pooled.

    Attributes:
        admissible: Instances that produced a measured radius.
        rejected: Counts of screened pairs by ``reject_reason``.
        censored: The ``_meta`` rows marking a per-network run that hit the cap.
        n_meta: Total ``_meta`` header rows, which are not instances.
    """

    admissible: list[dict]
    rejected: Counter
    censored: list[dict]
    n_meta: int

    @property
    def n_screened(self) -> int:
        """Ordered (X, Y) pairs screened, admissible ones included."""
        return len(self.admissible) + sum(self.rejected.values())

    def footer(self) -> str:
        """A three-line account of what was screened out, for a figure caption."""
        top = ", ".join(f"{n:,} {reason}" for reason, n in self.rejected.most_common(3))
        nets = sorted({c["network"] for c in self.censored})
        covs = ", ".join(f"{c['coverage']:g}" for c in self.censored)
        cap = {c["wall_until_timeout_s"] for c in self.censored}
        cens = (
            f"{len(self.censored)} censored per-network runs ({'/'.join(nets)} at coverage "
            f"{covs}, each stopped at the {max(cap):.0f} s wall cap):\nthe pairs those runs "
            "had not reached were never screened and are in no total above."
            if self.censored
            else "No run was censored."
        )
        return (
            f"{len(self.admissible)} admissible instances out of {self.n_screened:,} ordered "
            f"pairs screened; {sum(self.rejected.values()):,} pairs rejected before measurement,\n"
            f"chiefly {top}.\n" + cens
        )


def _admissible() -> Screening:
    """Read ``instances.jsonl``, splitting instances from ``_meta`` bookkeeping rows.

    Returns:
        The screening account.

    Raises:
        AssertionError: If any admissible instance carries the UNREACHED radius
            sentinel. Rejected rows carry ``radius = -1`` as filler, so the
            check is meaningful only on the admissible slice; the session claims
            zero UNREACHED and this asserts it instead of trusting it.
    """
    keep: list[dict] = []
    rejected: Counter = Counter()
    censored: list[dict] = []
    n_meta = 0
    with (RES / "instances.jsonl").open() as fh:
        for line in fh:
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("_meta"):
                n_meta += 1
                if row.get("censored"):
                    censored.append(row)
                continue
            if row.get("admissible"):
                keep.append(row)
            else:
                rejected[row.get("reject_reason") or "unspecified"] += 1

    sentinels = [r for r in keep if r["radius"] == UNREACHED]
    assert not sentinels, (
        f"{len(sentinels)} admissible instances carry the UNREACHED radius sentinel "
        "(-1); they must not be plotted or averaged as numbers."
    )
    return Screening(admissible=keep, rejected=rejected, censored=censored, n_meta=n_meta)


def _separation_split(rows: list[dict]) -> tuple[Counter, int]:
    """Split separations into measured integers and the undefined status.

    ``separation`` is ``null`` exactly when ``separation_status`` is
    ``no_z_member_in_component``: X's chain component contains no member of the
    adjustment set, so there is no path inside the component to measure and the
    quantity does not exist. It is not zero and it is not small.

    Args:
        rows: Admissible instances.

    Returns:
        ``(counts_by_measured_separation, n_undefined)``.

    Raises:
        AssertionError: If a row's ``separation`` and ``separation_status``
            disagree about whether the quantity exists.
    """
    measured: Counter = Counter()
    undefined = 0
    for r in rows:
        sep, status = r["separation"], r["separation_status"]
        if sep is None:
            assert status == NO_Z_STATUS, f"null separation with unexpected status {status!r}"
            undefined += 1
        else:
            assert status == MEASURED_STATUS, (
                f"numeric separation {sep} with unexpected status {status!r}"
            )
            measured[int(sep)] += 1
    return measured, undefined


def fig_undirected_fraction() -> list[Path]:
    """The descriptive headline: how much of a real CPDAG is left undirected.

    Per network, the share of edges the CPDAG does not compel. The median is
    10.7%, so there is normally something to orient and the robustness question
    is well posed. Six networks sit at exactly zero: their CPDAG is the DAG, no
    background knowledge can change anything, and the framework has nothing to
    say about them. Those six are drawn in a different colour and named, because
    "the measurement returned zero" is a finding about the network and not a
    small number to be squinted at.
    """
    _style()
    desc = _descriptive()
    nets = sorted(desc.networks, key=lambda n: (n.undirected_fraction, n.name))
    med = st.median(n.undirected_fraction for n in nets)
    compelled = [n for n in nets if n.n_undirected == 0]

    fig, ax = plt.subplots(figsize=(7.6, 8.8))
    ys = list(range(len(nets)))
    vals = [100 * n.undirected_fraction for n in nets]
    colours = [C_VERM if n.n_undirected == 0 else C_BLUE for n in nets]
    ax.barh(ys, vals, 0.72, color=colours, edgecolor="white", lw=0.4)
    for y, n in enumerate(nets):
        if n.n_undirected == 0:
            ax.text(
                13.0,
                y,
                "0 — framework cannot apply",
                va="center",
                fontsize=6.6,
                color=C_VERM,
                fontweight="bold",
            )
        else:
            ax.text(
                100 * n.undirected_fraction + 1.2,
                y,
                f"{100 * n.undirected_fraction:.1f}%  ({n.n_undirected}/{n.n_edges})",
                va="center",
                fontsize=6.4,
                color="black",
            )
    ax.axvline(100 * med, color=C_GREEN, ls="--", lw=1.4, zorder=3)
    ax.set_yticks(ys)
    ax.set_yticklabels([n.name for n in nets], fontsize=6.8)
    ax.set_ylim(-1.0, len(nets) + 1.2)
    ax.set_xlim(0, 118)
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.set_xlabel("share of CPDAG edges left undirected (%)")
    ax.set_ylabel("network")
    ax.grid(axis="y", visible=False)
    ax.set_title(
        f"The median published network leaves {100 * med:.1f}% of its edges undirected —\n"
        f"but {len(compelled)} of {len(nets)} are fully compelled and have nothing to orient",
        fontsize=11,
    )
    ax.legend(
        handles=[
            Patch(facecolor=C_BLUE, label="something left to orient"),
            Patch(
                facecolor=C_VERM,
                label=f"fully compelled — framework cannot apply ({len(compelled)})",
            ),
            plt.Line2D([], [], ls="--", color=C_GREEN, label=f"median = {100 * med:.1f}%"),
        ],
        fontsize=7.4,
        loc="center right",
        bbox_to_anchor=(1.0, 0.42),
        framealpha=0.96,
    )
    fig.text(
        0.5,
        -0.015,
        f"Fully compelled: {', '.join(n.name for n in compelled)}.\n{desc.footer()}",
        ha="center",
        va="top",
        fontsize=7,
        color=C_GREY,
    )
    return _save(fig, "s6_f1_undirected_fraction")


def fig_component_sizes() -> list[Path]:
    """Real chain components run far past anything either synthetic study built.

    The largest chain component per network on a symmetric-log axis, so the six
    networks with no undirected edges at all can be shown at a true zero instead
    of being dropped by a log scale. Reference lines mark session 4's largest
    synthetic component (6, at n = 24) and session 5's designed generator
    maximum (12). The median real network sits at 4, below both, but the tail
    reaches 85 — seven times the larger synthetic ceiling — which is where the
    cost results of session 4 stop being reassuring.
    """
    _style()
    desc = _descriptive()
    nets = sorted(desc.networks, key=lambda n: (n.max_component, n.name))
    sizes = [n.max_component for n in nets]
    med = st.median(sizes)
    biggest = max(nets, key=lambda n: n.max_component)
    beyond = [n for n in nets if n.max_component > S5_MAX_COMPONENT]
    zeros = [n for n in nets if n.max_component == 0]

    fig, ax = plt.subplots(figsize=(7.8, 8.8))
    ys = list(range(len(nets)))
    colours = []
    for n in nets:
        if n.max_component == 0:
            colours.append(C_GREY)
        elif n.max_component > S5_MAX_COMPONENT:
            colours.append(C_VERM)
        else:
            colours.append(C_BLUE)
    ax.barh(ys, sizes, 0.72, color=colours, edgecolor="white", lw=0.4)
    ax.set_xscale("symlog", linthresh=1, linscale=0.45)
    for y, n in enumerate(nets):
        txt = "0 — no undirected edges" if n.max_component == 0 else str(n.max_component)
        ax.text(
            max(n.max_component, 0) * 1.20 + 0.07,
            y,
            txt,
            va="center",
            fontsize=6.4,
            color=C_GREY if n.max_component == 0 else "black",
        )
    ax.axvline(S4_MAX_COMPONENT, color=C_ORANGE, ls="--", lw=1.4, zorder=3)
    ax.axvline(S5_MAX_COMPONENT, color=C_GREEN, ls="-.", lw=1.4, zorder=3)
    ax.axvline(med, color="black", ls=":", lw=1.2, zorder=3)
    ax.set_yticks(ys)
    ax.set_yticklabels([n.name for n in nets], fontsize=6.8)
    ax.set_ylim(-1.0, len(nets) + 1.4)
    ax.set_xlim(0, 260)
    ax.set_xticks([0, 1, 2, 4, 6, 12, 24, 85])
    ax.set_xticklabels(["0", "1", "2", "4", "6", "12", "24", "85"])
    ax.minorticks_off()
    ax.set_xlabel("largest chain component (nodes, symmetric-log axis; 0 shown at true zero)")
    ax.set_ylabel("network")
    ax.grid(axis="y", visible=False)
    ax.annotate(
        f"{biggest.name}: {biggest.max_component} nodes,\n"
        f"{biggest.max_component / S5_MAX_COMPONENT:.0f}x the largest synthetic component",
        xy=(biggest.max_component, len(nets) - 1.35),
        xytext=(14.0, len(nets) - 20.5),
        fontsize=7.5,
        color=C_VERM,
        arrowprops={"arrowstyle": "->", "color": C_VERM, "lw": 1.0},
    )
    ax.set_title(
        f"Real chain components have a heavy tail: median {med:.0f}, maximum "
        f"{biggest.max_component}\n"
        f"({len(beyond)} networks exceed session 5's designed ceiling of {S5_MAX_COMPONENT}; "
        f"session 4 never built one past {S4_MAX_COMPONENT})",
        fontsize=11,
    )
    ax.legend(
        handles=[
            Patch(
                facecolor=C_BLUE, label=f"at or below the synthetic ceiling of {S5_MAX_COMPONENT}"
            ),
            Patch(facecolor=C_VERM, label=f"beyond it ({len(beyond)} networks)"),
            Patch(facecolor=C_GREY, label=f"no undirected edges at all ({len(zeros)})"),
            plt.Line2D(
                [], [], ls="--", color=C_ORANGE, label=f"session 4 maximum = {S4_MAX_COMPONENT}"
            ),
            plt.Line2D(
                [], [], ls="-.", color=C_GREEN, label=f"session 5 maximum = {S5_MAX_COMPONENT}"
            ),
            plt.Line2D([], [], ls=":", color="black", label=f"real median = {med:.0f}"),
        ],
        fontsize=7,
        loc="lower right",
        bbox_to_anchor=(1.0, 0.02),
        framealpha=0.96,
    )
    fig.text(0.5, -0.015, desc.footer(), ha="center", va="top", fontsize=7, color=C_GREY)
    return _save(fig, "s6_f2_component_sizes")


def fig_separation() -> list[Path]:
    """The headline: separation >= 2 is ordinary in real networks and rare in Erdős-Rényi.

    Counts of the measured separation — the graph distance inside X's undirected
    component from X to the nearest member of the optimal adjustment set — over
    the admissible instances, with the
    undefined status drawn apart under its own hatch. The dashed marker on the
    s = 1 bar is where session 5's Erdős-Rényi control would have put its mass:
    85.8% of measurable separations at s = 1, leaving 14.2% at s >= 2. Real
    networks put 45.4% at s >= 2, a factor of 3.2, and reach s = 7 where the
    control's maximum was 4.
    """
    _style()
    scr = _admissible()
    measured, undefined = _separation_split(scr.admissible)
    seps = sorted(measured)
    n_meas = sum(measured.values())
    n_ge2 = sum(v for s, v in measured.items() if s >= 2)
    pct_ge2 = 100 * n_ge2 / n_meas

    fig, ax = plt.subplots(figsize=(9.0, 5.2))
    xs = list(range(len(seps)))
    vals = [measured[s] for s in seps]
    ax.bar(xs, vals, 0.62, color=C_BLUE, label=f"measured separation s ({n_meas} instances)")
    gap_x = len(seps) + 0.9
    ax.bar(
        [gap_x],
        [undefined],
        0.62,
        color=C_ORANGE,
        hatch="///",
        edgecolor="white",
        lw=0.6,
        label=(
            f"status: no Z-member in X's component ({undefined}) —\n"
            "separation is undefined here, NOT a separation of 0, and is never averaged"
        ),
    )
    ax.axvline(len(seps) - 0.15, color=C_GREY, ls=":", lw=1.1)

    er_s1 = ER_S1_SHARE * n_meas
    ax.hlines(er_s1, -0.35, 0.35, color=C_VERM, ls="--", lw=1.8, zorder=4)
    ax.annotate(
        f"Erdős-Rényi control (session 5)\nwould put {100 * ER_S1_SHARE:.1f}% here "
        f"(≈{er_s1:.0f} instances)",
        xy=(0.36, er_s1),
        xytext=(0.95, er_s1 * 0.80),
        fontsize=8,
        color=C_VERM,
        va="center",
        arrowprops={"arrowstyle": "->", "color": C_VERM, "lw": 1.0},
    )
    top = max(max(vals), undefined)
    brace_y = top * 1.10
    ax.annotate(
        "",
        xy=(0.62, brace_y),
        xytext=(len(seps) - 1 + 0.35, brace_y),
        arrowprops={"arrowstyle": "<->", "color": C_GREEN, "lw": 1.3},
    )
    ax.text(
        (0.62 + len(seps) - 0.65) / 2,
        brace_y * 1.035,
        f"s ≥ 2 in {n_ge2} of {n_meas} = {pct_ge2:.1f}%   (Erdős-Rényi: {ER_GE2_PCT}%)",
        ha="center",
        fontsize=9,
        color=C_GREEN,
        fontweight="bold",
    )
    labelled: list[tuple[float, int]] = [(float(x), vals[i]) for i, x in enumerate(xs)]
    labelled.append((gap_x, undefined))
    for x, v in labelled:
        share = 100 * v / n_meas if x != gap_x else 100 * v / len(scr.admissible)
        of = "of measurable" if x != gap_x else "of admissible"
        # Bars below ~5% of the tallest cannot carry two lines of text without
        # colliding with their neighbours; those get the count alone.
        text = f"{v}\n{share:.1f}% {of}" if v > 0.05 * top else str(v)
        ax.text(x, v + top * 0.02, text, ha="center", fontsize=7.6)
    ax.set_xticks([*xs, gap_x])
    ax.set_xticklabels([str(s) for s in seps] + ["undefined\n(status, not a number)"], fontsize=8.5)
    ax.set_xlabel(
        "separation s: graph distance inside X's undirected component,\n"
        "from the treatment X to the nearest member of O(G0)"
    )
    ax.set_ylabel("admissible instances")
    ax.set_ylim(0, top * 1.42)
    ax.grid(axis="x", visible=False)
    ax.set_title(
        f"Separation ≥ 2 is common in real networks ({pct_ge2:.1f}% of measurable) "
        f"and rare in Erdős-Rényi ({ER_GE2_PCT}%)\n"
        f"— and the real maximum is {max(seps)}, where the random control never passed 4",
        fontsize=11,
    )
    ax.legend(
        fontsize=7.6,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.40),
        ncol=1,
        framealpha=0.96,
    )
    fig.text(0.5, -0.36, scr.footer(), ha="center", fontsize=7, color=C_GREY, wrap=True)
    return _save(fig, "s6_f3_separation")


def fig_radius_and_coverage() -> list[Path]:
    """What the radius actually is out of sample, and how background knowledge moves it.

    Left: the radius distribution over all admissible instances. Every one is
    exact and finite — there is not a single UNREACHED (-1) sentinel in the
    file, which :func:`_admissible` asserts — so the axis is a distribution and
    not a distribution plus a censoring bucket.

    Right: the share of instances at r = 1 against background-knowledge
    coverage. The share falls monotonically as knowledge is withdrawn (69.8% ->
    51.1% -> 34.0%), so with less knowledge fewer instances sit at the floor.
    The maximum radius moves the other way (14 -> 4 -> 4): withdrawing knowledge
    lifts the typical instance off r = 1 while collapsing the ceiling, which is
    session 5's ``min(s, |K_G0|)`` cap showing up out of sample. Both directions
    are drawn because reporting only the r = 1 rate would read as a single
    monotone story that the tail contradicts.
    """
    _style()
    scr = _admissible()
    rows = scr.admissible
    radii = Counter(r["radius"] for r in rows)
    span = list(range(min(radii), max(radii) + 1))

    covs = sorted({r["coverage"] for r in rows}, reverse=True)
    by_cov = {c: [r for r in rows if r["coverage"] == c] for c in covs}
    r1_pct = {c: 100 * sum(1 for r in v if r["radius"] == 1) / len(v) for c, v in by_cov.items()}
    r_max = {c: max(r["radius"] for r in v) for c, v in by_cov.items()}

    fig, (axl, axr) = plt.subplots(
        1, 2, figsize=(11.4, 4.6), gridspec_kw={"width_ratios": [1.5, 1]}
    )

    axl.bar(span, [radii.get(r, 0) for r in span], 0.66, color=C_BLUE, edgecolor="white", lw=0.4)
    for r in span:
        n = radii.get(r, 0)
        if n:
            axl.text(r, n * 1.14, str(n), ha="center", fontsize=7.4)
    axl.set_yscale("log")
    axl.set_ylim(0.7, max(radii.values()) * 3.4)
    axl.set_xticks(span)
    axl.set_xlabel("back-door robustness radius  r  (all coverages pooled)")
    axl.set_ylabel("admissible instances (log)")
    axl.grid(axis="x", visible=False)
    axl.set_title(
        f"Every one of the {len(rows)} radii is exact and finite:\n"
        "0 UNREACHED (−1) sentinels, so nothing here is censored",
        fontsize=9.8,
    )
    axl.text(
        0.98,
        0.88,
        f"median r = {st.median([r['radius'] for r in rows]):.0f}\n"
        f"maximum r = {max(radii)}\n"
        f"{100 * radii[1] / len(rows):.1f}% sit at r = 1",
        transform=axl.transAxes,
        ha="right",
        va="top",
        fontsize=7.8,
        color=C_GREY,
    )

    xs = list(range(len(covs)))
    vals = [r1_pct[c] for c in covs]
    axr.plot(xs, vals, "o-", color=C_BLUE, lw=2.0, ms=9, label="share of instances at r = 1")
    for x, c in enumerate(covs):
        axr.annotate(
            f"{r1_pct[c]:.1f}%",
            xy=(x, r1_pct[c]),
            xytext=(0, 11),
            textcoords="offset points",
            ha="center",
            fontsize=9,
            fontweight="bold",
            color=C_BLUE,
        )
        axr.text(
            x,
            0.045,
            f"n = {len(by_cov[c])}\nmax r = {r_max[c]}",
            transform=axr.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=7.8,
            color=C_VERM,
        )
    axr.set_xticks(xs)
    axr.set_xticklabels([f"{c:g}" for c in covs])
    axr.set_xlim(-0.45, len(covs) - 0.55)
    axr.set_ylim(20, 84)
    axr.set_xlabel("background-knowledge coverage")
    axr.set_ylabel("instances with r = 1 (%)")
    axr.set_title(
        "Withdrawing knowledge lifts instances off r = 1\n"
        f"({vals[0]:.1f}% → {vals[-1]:.1f}%) while the ceiling collapses "
        f"({r_max[covs[0]]} → {r_max[covs[-1]]})",
        fontsize=9.8,
    )
    axr.legend(fontsize=7.6, loc="upper right", framealpha=0.96)
    axr.grid(axis="x", visible=False)

    fig.suptitle(
        "Radius out of sample: exact everywhere, and set by how much knowledge is imposed",
        fontsize=11.5,
        y=1.03,
    )
    fig.text(0.5, -0.13, scr.footer(), ha="center", fontsize=7, color=C_GREY, wrap=True)
    return _save(fig, "s6_f4_radius_and_coverage")


def build_all() -> list[Path]:
    """Build every session 6 figure that has data.

    Returns:
        The paths written. A figure whose input file is missing or malformed is
        skipped with a message rather than aborting the run; an assertion
        failure on a sentinel is *not* caught, because that would mean the
        file's own bookkeeping is wrong.
    """
    out: list[Path] = []
    for fn in (
        fig_undirected_fraction,
        fig_component_sizes,
        fig_separation,
        fig_radius_and_coverage,
    ):
        try:
            out.extend(fn())
        except (FileNotFoundError, ValueError, IndexError, KeyError) as exc:
            print(f"skipped {fn.__name__}: {exc}")
    return out


if __name__ == "__main__":
    for p in build_all():
        print(p)
