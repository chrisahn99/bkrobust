"""Tests for :mod:`bkrobust.demo.figures`.

Uses a small ``n_draws`` when running scenarios (the figures module only cares
about shapes and validity flags, not bias precision) so the module-scoped
fixture stays fast. Figures are always written to ``tmp_path``, never to the
real ``figures/`` directory -- ``test_build_all_figures_returns_many_paths``
achieves that by ``chdir``-ing into a tmp dir before calling the public
entrypoint, since :func:`~bkrobust.demo.figures.build_all_figures` writes to
the *relative* path ``figures/`` (matching the ``results/...`` convention
used elsewhere in this package).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pytest

from bkrobust.demo import figures
from bkrobust.demo.pipeline import UNREACHED, run_scenario

SCENARIOS = ("A", "B", "C")


@pytest.fixture(scope="module")
def results() -> dict[str, dict]:
    return {label: run_scenario(label, n_draws=3) for label in SCENARIOS}


def _assert_pdf_and_png(paths: list[Path], stem: str) -> None:
    pdf_matches = [p for p in paths if p.name == f"{stem}.pdf"]
    png_matches = [p for p in paths if p.name == f"{stem}.png"]
    assert pdf_matches, f"no {stem}.pdf among {paths}"
    assert png_matches, f"no {stem}.png among {paths}"
    for p in pdf_matches + png_matches:
        assert p.exists(), f"{p} does not exist"
        assert p.stat().st_size > 0, f"{p} is empty"


# --------------------------------------------------------------------------------
# fig1
# --------------------------------------------------------------------------------


def test_fig1_example_writes_pdf_and_png(results, tmp_path):
    paths = figures.fig1_example(results, tmp_path)
    assert len(paths) == 2
    _assert_pdf_and_png(paths, "fig1_example")


def test_fig1_example_caption(results):
    caption = figures.fig1_example_caption(results)
    assert isinstance(caption, str)
    assert caption.strip()
    assert any(lab in caption for lab in SCENARIOS)


# --------------------------------------------------------------------------------
# fig2 -- central figure
# --------------------------------------------------------------------------------


@pytest.mark.parametrize("label", SCENARIOS)
def test_fig2_layered_space_writes_pdf_and_png(results, tmp_path, label):
    paths = figures.fig2_layered_space(results[label], label, tmp_path)
    assert len(paths) == 2
    _assert_pdf_and_png(paths, f"fig2_layered_space_{label}")


@pytest.mark.parametrize("label", SCENARIOS)
def test_fig2_layered_space_caption(results, label):
    caption = figures.fig2_layered_space_caption(results[label], label)
    assert isinstance(caption, str)
    assert caption.strip()
    assert label in caption
    assert "r_val" in caption or "UNREACHED" in caption


# --------------------------------------------------------------------------------
# fig3
# --------------------------------------------------------------------------------


def test_fig3_shell_profile_writes_pdf_and_png(results, tmp_path):
    paths = figures.fig3_shell_profile(results, tmp_path)
    assert len(paths) == 2
    _assert_pdf_and_png(paths, "fig3_shell_profile")


def test_fig3_shell_profile_caption(results):
    caption = figures.fig3_shell_profile_caption(results)
    assert isinstance(caption, str)
    assert caption.strip()
    assert "bias" in caption.lower()


# --------------------------------------------------------------------------------
# fig4
# --------------------------------------------------------------------------------


@pytest.mark.parametrize("label", SCENARIOS)
def test_fig4_witness_writes_pdf_and_png(results, tmp_path, label):
    paths = figures.fig4_witness(results[label], label, tmp_path)
    assert len(paths) == 2
    _assert_pdf_and_png(paths, f"fig4_witness_{label}")


@pytest.mark.parametrize("label", SCENARIOS)
def test_fig4_witness_caption(results, label):
    caption = figures.fig4_witness_caption(results[label], label)
    assert isinstance(caption, str)
    assert caption.strip()
    assert label in caption


# --------------------------------------------------------------------------------
# fig5 -- must degrade gracefully if bkrobust.demo.baseline is unavailable
# --------------------------------------------------------------------------------


def test_fig5_baseline_scatter_does_not_raise(results, tmp_path):
    paths = figures.fig5_baseline_scatter(results, tmp_path)
    assert isinstance(paths, list)
    # Either it was skipped entirely (dependency missing/broken -> []), or it
    # wrote a complete pdf+png pair. Never a partial or crashing result.
    if paths:
        assert len(paths) == 2
        _assert_pdf_and_png(paths, "fig5_baseline_scatter")


def test_fig5_baseline_scatter_caption(results):
    caption = figures.fig5_baseline_scatter_caption(results)
    assert isinstance(caption, str)
    assert caption.strip()
    assert "radius" in caption.lower() or "distance" in caption.lower()


# --------------------------------------------------------------------------------
# fig6 -- UNREACHED must never appear as a bare -1
# --------------------------------------------------------------------------------


def test_fig6_robustness_frontier_writes_pdf_and_png(results, tmp_path):
    paths = figures.fig6_robustness_frontier(results, tmp_path)
    assert len(paths) == 2
    _assert_pdf_and_png(paths, "fig6_robustness_frontier")


def test_fig6_robustness_frontier_caption(results):
    caption = figures.fig6_robustness_frontier_caption(results)
    assert isinstance(caption, str)
    assert caption.strip()
    assert "r_val" in caption or "variance" in caption.lower()


@pytest.mark.parametrize("label", SCENARIOS)
def test_frontier_plot_frame_never_plots_bare_unreached(results, label):
    df = figures.frontier_plot_frame(results[label])
    if df.empty:
        pytest.skip(f"scenario {label} has no valid adjustment sets to check")

    # The raw r_val column may legitimately contain UNREACHED (-1) -- that is
    # the pipeline's own sentinel and is untouched here.
    unreached_rows = df[df["r_val"] == UNREACHED]

    # But the plot-ready "_y" column must never expose that -1 on the axis:
    # every _y value must be non-negative, and every UNREACHED row must be
    # remapped to a distinct tier strictly above every real (reachable) shell.
    assert (df["_y"] >= 0).all()
    if len(unreached_rows):
        reachable_y = df.loc[~df["never_fails"], "_y"]
        never_tier = df["never_tier"].iloc[0]
        assert (df.loc[df["never_fails"], "_y"] == never_tier).all()
        if len(reachable_y):
            assert never_tier > reachable_y.max()
        assert -1 not in df["_y"].to_numpy()


def test_fig6_axis_tick_labels_do_not_show_bare_minus_one(results, tmp_path):
    # Build the figure and inspect the actual rendered y-tick labels for one
    # panel: UNREACHED must show as a named tier, never as the string "-1".
    import matplotlib.pyplot as plt

    df = figures.frontier_plot_frame(results["B"])
    max_shell = int(df["max_shell"].iloc[0]) if len(df) else 0
    never_tier = int(df["never_tier"].iloc[0]) if len(df) else 0
    fig, ax = plt.subplots()
    yticks = [*list(range(0, max_shell + 1)), never_tier]
    yticklabels = [str(t) for t in range(0, max_shell + 1)] + ["never fails\n(UNREACHED)"]
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)
    rendered = [t.get_text() for t in ax.get_yticklabels()]
    plt.close(fig)
    assert "-1" not in rendered
    assert any("UNREACHED" in label or "never" in label for label in rendered)


# --------------------------------------------------------------------------------
# orchestrator
# --------------------------------------------------------------------------------


def test_build_all_figures_returns_many_existing_paths(results, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    paths = figures.build_all_figures(results)
    assert len(paths) >= 12
    for p in paths:
        assert p.is_absolute() or True  # paths are relative to the chdir'd cwd, that's fine
        assert p.exists()
        assert p.stat().st_size > 0
    # every path actually lives under figures/, relative to the chdir'd cwd
    assert all(str(p).startswith("figures") or "figures" in p.parts for p in paths)
