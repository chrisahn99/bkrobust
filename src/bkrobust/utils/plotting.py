"""Shared matplotlib style for the paper figures.

One style, applied everywhere, so figures do not drift apart between
experiments. Conventions that follow from what these figures have to show:

* **Colour-blind-safe palette, and never colour alone.** The key figures
  distinguish three regimes -- benign, efficiency loss, bias -- and a reader who
  cannot separate the colours must still be able to read them, so regimes carry
  hatching or line style as well.
* **The two radii are drawn as vertical rules, always in the same style.** They
  appear in most figures and should be recognisable without the legend.
* **The band between them is shaded.** It is the paper's central object; it
  should be visible before the caption is read.
* **Vector output.** PDF for the paper, PNG only for notebooks.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    from pathlib import Path

    import matplotlib.pyplot as plt
    import pandas as pd

#: Regime colours, colour-blind safe. Keys match
#: :class:`~bkrobust.knowledge.taxonomy.Severity` values.
REGIME_COLORS: Mapping[str, str] = {
    "benign": "#4C9F70",
    "efficiency": "#E1A730",
    "validity": "#C1442E",
    "unidentifiable": "#6B5B95",
}

#: Hatch patterns, so the regimes remain distinguishable in greyscale.
REGIME_HATCH: Mapping[str, str] = {
    "benign": "",
    "efficiency": "//",
    "validity": "xx",
    "unidentifiable": "..",
}


def set_style(*, context: str = "paper", font_scale: float = 1.0) -> None:
    """Apply the shared matplotlib rcParams.

    Args:
        context: ``"paper"``, ``"notebook"`` or ``"talk"``; sets sizes.
        font_scale: Multiplier on the font sizes.
    """
    raise NotImplementedError


def figure_size(width: str = "single", aspect: float = 0.618) -> tuple[float, float]:
    """Return a figure size in inches matching the ICLR column widths.

    Args:
        width: ``"single"``, ``"double"`` or ``"full"``.
        aspect: Height as a fraction of width; the default is the golden ratio.

    Returns:
        ``(width_inches, height_inches)``.
    """
    raise NotImplementedError


def plot_breakdown_curve(
    df: pd.DataFrame,
    *,
    ax: plt.Axes | None = None,
    delta_col: str = "delta",
    bias_col: str = "abs_bias",
    variance_col: str = "variance_ratio",
    delta_opt: int | None = None,
    delta_valid: int | None = None,
) -> plt.Axes:
    """Plot bias and variance inflation against the perturbation radius.

    The paper's main figure: two curves on twin axes, vertical rules at the two
    radii, the band between them shaded. Variance lifts off at ``delta_opt`` and
    bias at ``delta_valid``; the visual claim is the gap between those two
    lift-off points.

    Args:
        df: Long-format results, one row per radius.
        ax: Axes to draw on; a new figure is created if ``None``.
        delta_col: Radius column.
        bias_col: Bias column.
        variance_col: Variance-ratio column.
        delta_opt: Where to draw the optimality rule; omitted if ``None``.
        delta_valid: Where to draw the validity rule.

    Returns:
        The axes.
    """
    raise NotImplementedError


def plot_severity_stack(
    df: pd.DataFrame,
    *,
    ax: plt.Axes | None = None,
    delta_col: str = "delta",
) -> plt.Axes:
    """Stacked area of severity classes against radius.

    Shows how the ball's composition shifts from benign through efficiency loss
    to bias -- the distributional companion to the worst-case radii, and the
    figure that says whether the failure mode is common or a corner case.
    """
    raise NotImplementedError


def plot_cascade_amplification(
    df: pd.DataFrame,
    *,
    ax: plt.Axes | None = None,
    group_col: str = "graph_family",
) -> plt.Axes:
    """Forced orientations against imposed constraints, by graph family.

    The amplification result: how far one wrong assertion travels, and whether
    hub-heavy graphs amplify harder than homogeneous ones.
    """
    raise NotImplementedError


def plot_cfm_audit(
    df: pd.DataFrame,
    *,
    ax: plt.Axes | None = None,
    checkpoint_col: str = "checkpoint",
) -> plt.Axes:
    """Bias against radius, one line per checkpoint, with the classical arm.

    The classical arm must be drawn in a visually distinct style: it is the
    reference the amortized curves are read against, not another checkpoint.
    """
    raise NotImplementedError


def save_figure(
    fig: Any,
    name: str,
    *,
    directory: str | Path = "paper/figures",
    formats: Sequence[str] = ("pdf",),
    dpi: int = 300,
) -> list[Path]:
    """Save a figure into the paper's figure directory.

    Args:
        fig: The matplotlib figure.
        name: Base filename, without extension.
        directory: Where to write.
        formats: Extensions to write.
        dpi: Raster resolution, for PNG output.

    Returns:
        The written paths.

    Note:
        Figures must carry no author or institution text, including in embedded
        font metadata. See ``docs/ANONYMITY_CHECKLIST.md``.
    """
    raise NotImplementedError
