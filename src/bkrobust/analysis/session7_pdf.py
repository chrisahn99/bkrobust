r"""Build `REPORT_FRAGILITY_AND_PARETO.pdf` — the session-7 report, same house style.

Seventh in the series, reusing the layout helpers from
:mod:`bkrobust.analysis.session_pdf` (Unicode font registration, wrapping table
cells, captioned figures) so the reports read as one sequence.

**Departure from sessions 3-5's discipline, stated up front.** Sessions 3, 4 and
5 derive every number from committed JSON/CSV at build time and never transcribe
from the markdown, on the theory that a rebuild is not a verification. This
module does the opposite: every number below is copied verbatim from
``report_fragility_and_pareto.md`` (the authority named in this session's brief)
rather than re-derived from ``results/axis_robustness/*.csv``. Three reasons,
all specific to this session rather than a rejection of the earlier discipline:

1. The brief explicitly forbids re-running any of the analysis stages
   (``run_survival``, ``run_hops``, ``run_pareto*``, ``run_analyse.py``) that
   produced these figures, and several of the reported statistics are
   10,000-resample bootstraps over a fixed seed — reproducing them here would
   mean re-implementing that analysis a second time, in a different module,
   with a real chance of silently drifting from the committed result rather
   than checking it.
2. The markdown itself already carries a documented, dated correction (the
   Section 2.7 contradiction-rate table, corrected 2026-09-10 after a pooling
   bug was caught while building the *figures* — see
   :mod:`bkrobust.analysis.session7_figures`). Re-deriving independently risks
   silently reproducing the *old* bug or a *new* one and shipping a PDF that
   once again disagrees with the markdown, which is precisely the failure mode
   sessions 3-5's discipline exists to prevent.
3. The brief's acceptance criterion (I2) is that every number in the PDF
   matches the markdown exactly — not that it is independently re-derived. A
   verbatim transcription, checked by eye against the source table-by-table, is
   the more direct way to satisfy that criterion without touching the
   forbidden pipelines.

Figures are embedded from the vector PDFs under ``figures/`` (rasterised in
memory at 300 dpi via PyMuPDF, four times the resolution of the 150 dpi PNGs
also shipped there) rather than the shared :func:`session_pdf.figure` helper,
which only reads PNGs. No new file is written to disk anywhere but the final
report.
"""

# ruff: noqa: RUF001
# This module is a DOCUMENT. Mathematical symbols and typographic punctuation in
# the rendered prose are intentional.

from __future__ import annotations

import io
import subprocess
from pathlib import Path
from typing import Any

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.platypus import Image, KeepTogether, PageBreak, SimpleDocTemplate, Spacer

from bkrobust.analysis.session_pdf import (
    GREY,
    bullets,
    para,
    styles,
    table,
)

FIGS = Path("figures")

RADIUS_NOTE = (
    "Conjecture 2 (hence Anti-Exchange Case B, verified not proved); the error "
    "is one-sided, so radii can only be too large."
)


def _figure(
    name: str, caption: str, st: dict[str, Any], max_w: float = 16.4 * cm, dpi: float = 300
) -> list[Any]:
    """Embed a figure rasterised from the vector PDF at ``dpi``, house-style caption.

    Mirrors :func:`bkrobust.analysis.session_pdf.figure` (same KeepTogether, same
    caption style, same width-capped scaling) but sources the vector PDF under
    ``figures/`` instead of the 150 dpi PNG, entirely in memory — nothing is
    written to disk. Falls back to the PNG via the shared helper if the PDF is
    missing.
    """
    pdf_path = FIGS / f"{name}.pdf"
    if not pdf_path.exists():
        from bkrobust.analysis.session_pdf import figure as _png_figure

        return _png_figure(name, caption, st, max_w)

    import fitz
    from PIL import Image as PILImage

    doc = fitz.open(str(pdf_path))
    page = doc[0]
    zoom = dpi / 72.0
    pix = page.get_pixmap(matrix=fitz.Matrix(zoom, zoom))
    png_bytes = pix.tobytes("png")
    doc.close()

    with PILImage.open(io.BytesIO(png_bytes)) as im:
        w, h = im.size
    scale = min(max_w / w, 1.0)
    img = Image(io.BytesIO(png_bytes), width=w * scale, height=h * scale)
    return [KeepTogether([img, para(caption, st["caption"])])]


def _radius_tag(st: dict[str, Any]) -> Any:
    """The standing caveat that must travel with any quoted radius."""
    return para(f"<i>Radius assumption: {RADIUS_NOTE}</i>", st["caption"])


def _session_commits() -> list[list[str]]:
    """This branch's own commits, newest first — not the shared house range.

    The reusable :func:`session_pdf.git_commits` hardcodes a commit range
    (``c201cfb..HEAD``) that is only correct for the synth/search branch
    sessions 3-5 shipped from. Session 7 lives on
    ``experiments/survival-and-pareto``, branched from ``main`` at a different
    point, so re-using that range verbatim would list unrelated history from
    the other branch's sessions. This derives the correct range instead of
    editing the shared helper.
    """
    try:
        base = subprocess.run(
            ["git", "merge-base", "main", "HEAD"], capture_output=True, text=True, check=True
        ).stdout.strip()
        out = subprocess.run(
            ["git", "log", "--format=%h\x1f%ad\x1f%s", "--date=short", f"{base}..HEAD"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
    except Exception:
        return []
    rows = []
    for line in out.split("\n"):
        if not line:
            continue
        h, d, s = line.split("\x1f")
        rows.append([h, d, s])
    return rows


def _cover(st: dict[str, Any]) -> list[Any]:
    return [
        para("Robustness survival and the Pareto frontier of causal assumptions", st["title"]),
        para(
            "Session 7 — Axis Robustness<br/>"
            "branch <font face='%FONTMONO%'>experiments/survival-and-pareto</font> · "
            "2026-09-10",
            st["subtitle"],
        ),
        table(
            [
                ["", ""],
                [
                    "Pre-registration",
                    "<i>results/axis_robustness/PREREGISTRATION.md</i> "
                    "(§1–§9 written before any survival curve or rank correlation existed; "
                    "Appendices A–F record every falsification, correction and retraction, "
                    "in order).",
                ],
                [
                    "Gate",
                    "<i>benchmarks.measure.fast_gate</i> "
                    "<b>exclusively</b>, recorded in every row. "
                    "<i>synth.runner.gate</i> was not used — it calls "
                    "the exponential <i>all_valid_adjustment_sets_mpdag</i> "
                    "internally. <i>fast_gate</i> is strictly stricter "
                    "(16/46,800 disagreements, one direction), which biases <b>against</b> the "
                    "hypothesis under test.",
                ],
                [
                    "Assumption carried on every radius",
                    f"<i>{RADIUS_NOTE}</i> All upward searches are exact iff it holds.",
                ],
            ],
            [4.3 * cm, 11.9 * cm],
            header=False,
        ),
        Spacer(1, 0.4 * cm),
        para("1. What was asked, and what the data actually supports", st["h1"]),
        para(
            "The brief asked for two demonstrations: that <font face='%FONTMONO%'>r_val</font> "
            "ranks average-case structural survival where SHD and |K| fail, and that precision "
            "trades off against robustness along a Pareto frontier (“the fragility of "
            "precision”).",
            st["body"],
        ),
        para(
            "<b>The first is supported. The second is false</b>, and the reason it is false is "
            "the more interesting result of the two.",
            st["good"],
        ),
        table(
            [
                ["Prediction", "Verdict", "Evidence"],
                [
                    "<b>P1</b> r_val positively ranks survival",
                    "<b>SUPPORTED</b>",
                    "τ_b &gt; 0 in 9/9 strata; 8/9 under the conservative endpoint",
                ],
                [
                    "<b>P2</b> r_val outranks |K|",
                    "<b>SUPPORTED</b>",
                    "9/9 strata; 7/9 under the conservative endpoint",
                ],
                [
                    "<b>P3</b> |K| is a null predictor (project hypothesis H4)",
                    "<b>NOT SUPPORTED</b>",
                    "CI covers 0 in only 3/9; |K| is <i>negatively</i> associated",
                ],
                [
                    "<b>P4</b> tiered AUC variance &gt; flip",
                    "<b>NOT COMPARABLE</b>",
                    "unit mismatch; no verdict forced (Appendix E.4)",
                ],
                [
                    "<b>P5</b> claim-axis and hop-axis agree in sign",
                    "<b>PASS</b>",
                    "6/6 strata, independent subsample",
                ],
                [
                    "<b>P6</b> a non-trivial efficiency/robustness frontier exists",
                    "<b>FALSIFIED</b>",
                    "utility exactly invariant in 32/32 strata",
                ],
                [
                    "<b>P6′</b> frontier exists for <i>achievable</i> efficiency",
                    "<b>FALSIFIED</b>",
                    "invariant again; Phase 2 abandoned per pre-registered trigger",
                ],
                [
                    "<b>P7 / P7′</b> aggressive BK trades utility for radius",
                    "<b>UNTESTABLE</b>",
                    "one axis is constant; census reported instead",
                ],
            ],
            [6.2 * cm, 3.0 * cm, 7.0 * cm],
        ),
        PageBreak(),
    ]


def _phase1_design(st: dict[str, Any]) -> list[Any]:
    return [
        para("2. Phase 1 — survival under corrupted background knowledge", st["h1"]),
        para("2.1 Design", st["h2"]),
        para(
            "An analyst holds CPDAG Ĉ, asserts knowledge K, Meek-closes to G₀, and reads "
            "Z* = optimal_adjustment_set_mpdag(G₀, x, y) <b>once, at depth 0, then holds it "
            "fixed</b>. K is then corrupted and we ask whether the set they already "
            "committed to survives — is_gac_valid_mpdag(G, x, y, Z*), polynomial, no "
            "enumeration of [G].",
            st["body"],
        ),
        para(
            "1,978 instances; 3,136,000 sampled states; 200 reps per depth; 1,656 s. "
            "Instances span component size 6–12, three separations, coverage ∈ {0.5, 1.0}, "
            "and a <b>base wrongness</b> b ∈ {0, 0.10, 0.25} applied before the sweep so that "
            "SHD(G₀, truth) is not constant — at coverage 1.0 with truthful K, G₀ <i>is</i> "
            "the ground-truth DAG and SHD ≡ 0, which would have made the baseline comparison "
            "rigged by construction.",
            st["body"],
        ),
        para(
            "Two arms, never merged: <b>flip</b> (uniform, independent, targeted depth) and "
            "<b>tiered</b> (block-correlated, re-binned on corruption_rate; see §2.5).",
            st["body"],
        ),
    ]


def _phase1_p1(st: dict[str, Any]) -> list[Any]:
    return [
        para("2.2 P1 — r_val ranks survival", st["h2"]),
        para(
            "τ_b(AUC, r_val), 10,000-resample bootstrap over instances, seed 0:",
            st["body"],
        ),
        table(
            [
                ["stratum", "τ_b", "95% CI", "n"],
                ["flip cov=0.5 bw=0.00", "+0.714", "[+0.625, +0.798]", "240"],
                ["flip cov=1.0 bw=0.00", "+0.519", "[+0.442, +0.591]", "240"],
                ["flip cov=1.0 bw=0.10", "+0.388", "[+0.296, +0.478]", "236"],
                ["flip cov=0.5 bw=0.10", "+0.380", "[+0.265, +0.488]", "240"],
                ["flip cov=1.0 bw=0.25", "+0.346", "[+0.209, +0.472]", "134"],
                ["flip cov=0.5 bw=0.25", "+0.031", "[−0.114, +0.169]", "220"],
                ["tiered n_tiers=3", "+0.811", "[+0.787, +0.833]", "238"],
                ["tiered n_tiers=2", "+0.799", "[+0.761, +0.834]", "190"],
                ["tiered n_tiers=4", "+0.777", "[+0.753, +0.798]", "240"],
            ],
            [4.6 * cm, 2.2 * cm, 4.4 * cm, 1.8 * cm],
            align_right=[1, 3],
        ),
        _radius_tag(st),
        para(
            "Positive in 9/9; CI excludes zero in 8/9. The exception is the highest-"
            "corruption, lowest-coverage flip cell, where the association is "
            "indistinguishable from zero — reported, not dropped.",
            st["body"],
        ),
        para(
            "<b>Effect sizes</b> move monotonically: for tiered n_tiers=3, median AUC rises "
            "0.67 → 0.94 → 0.99 → 0.98 → 1.00 across r_val buckets 1→5+.",
            st["body"],
        ),
        Spacer(1, 0.2 * cm),
        *_figure(
            "session7_f1_stratification",
            "<b>Figure 1.</b> AUC_frac by r_val bucket vs. by shd_truth bucket — the "
            "endpoint. The r_val panel (left) rises overall but <b>not monotonically</b> — "
            "the median dips at bucket 3 — while the shd_truth panel shows no comparable "
            "structure. " + RADIUS_NOTE,
            st,
        ),
        Spacer(1, 0.3 * cm),
        para(
            "<b>The primary stratified evidence for P1/P2.</b> Figure 2 below breaks the same "
            "association down by predictor across all nine strata — this is the figure the "
            "verdicts above rest on, not the pooled scatter above it.",
            st["callout"],
        ),
        *_figure(
            "session7_f3_tau_forest",
            "<b>Figure 2.</b> τ_b per stratum, four predictors: r_val, |K|, "
            "shd_truth = SHD(G₀, true DAG), and shd_cpdag = SHD(Ĉ, G₀), the orientations "
            "committed including Meek propagation. 95% bootstrap CIs over instances. "
            "Endpoint is AUC_frac for flip strata and AUC_frac_rate for tiered strata, as "
            "the axis states. This is "
            "the primary evidence figure for P1 and P2: r_val's interval is positive in 8 of "
            "9 strata and dominates |K| in every stratum shown. " + RADIUS_NOTE,
            st,
        ),
        PageBreak(),
    ]


def _phase1_p2p3(st: dict[str, Any]) -> list[Any]:
    return [
        para("2.3 P2/P3 — the baselines fail, but not the way the brief assumed", st["h2"]),
        para(
            "r_val outranks |K| in <b>9/9</b> strata. But |K| is <b>not inert</b>, so H4's "
            "framing is wrong: its association with survival is real and <b>negative</b> — "
            "flip −0.194/−0.205 at high corruption, tiered −0.425 to −0.470, all with CIs "
            "excluding zero.",
            st["body"],
        ),
        para("More asserted claims means <i>worse</i> survival, not better.", st["callout"]),
        para(
            "SHD(G₀, truth) fails for a different reason again: it is "
            "<font face='%FONTMONO%'>undefined (predictor constant)</font> at b = 0 <b>by "
            "construction</b>, and elsewhere its sign flips across strata (+0.240, +0.027, "
            "−0.175, −0.171, −0.134). It has no consistent direction. That is a fairer "
            "indictment than “SHD tangles”: a third of the design cannot score it "
            "at all.",
            st["body"],
        ),
    ]


def _phase1_discordance(st: dict[str, Any]) -> list[Any]:
    return [
        para("2.4 The discordance case", st["h2"]),
        para(
            "<font face='%FONTMONO%'>c06_s01_cov050_seed00000007_bw000</font>. "
            "shd_truth = 0 — G₀ <i>is</i> the ground-truth DAG, the analyst's knowledge is "
            "perfectly correct — with n_k = 2. Yet AUC_frac = 0.0: S(1) = 0.0 on 107 "
            "non-contradictory samples, S(2) = 0.0 on 200. Every corruption invalidates the "
            "committed set, and r_val = 1 says exactly that.",
            st["body"],
        ),
        para(
            "<b>Zero SHD is compatible with maximal fragility.</b> SHD and |K| describe how "
            "much was asserted; r_val describes how much can go wrong before the answer "
            "breaks.",
            st["good"],
        ),
        *_figure(
            "session7_f4_discordance_spotlight",
            "<b>Figure 3.</b> The verified discordance case — zero structural error, zero "
            "survival. " + RADIUS_NOTE,
            st,
        ),
        para(
            "A second, opposite spotlight (AUC_frac = 1.0 at r_val = 8) was "
            "<b>withdrawn on verification</b>: its AUC rested on depths with n_eval of 35, "
            "10, 1 and 1. Four surviving draws are not a survival curve.",
            st["callout"],
        ),
        PageBreak(),
    ]


def _phase1_defects(st: dict[str, Any]) -> list[Any]:
    return [
        para("2.5 Two defects found and corrected mid-analysis", st["h2"]),
        para(
            "<b>The tiered depth metric was broken (Appendix E).</b> It counted assertions "
            "added or reversed but not <i>removed</i>; tiered relocates nodes, and "
            "relocation into a shared tier silently deletes a cross-tier assertion. Its "
            "d = 0 bin was contaminated by every corruption rate and contained genuine "
            "failures. Re-binning on corruption_rate — the knob actually turned — gives a "
            "clean design: contradiction rises monotonically 0.000 → 0.852, S falls "
            "1.000 → 0.618, S = 1.000 exactly at rate 0 for all 668 instances individually, "
            "n = 133,600 per rate. Tiered τ values <b>survive</b> the correction "
            "(0.777–0.811). The flip arm was never affected.",
            st["callout"],
        ),
        *_figure(
            "session7_f2_tiered_stratification",
            "<b>Figure 4.</b> The tiered arm on its own rate axis, after re-binning on "
            "corruption_rate rather than the defective node-count depth. " + RADIUS_NOTE,
            st,
        ),
        para(
            "<b>Effective sample size varies with r_val (Appendix F).</b> τ(r_val, usable "
            "depths) runs from +0.318 to −0.438 depending on stratum, so the bias runs "
            "<i>against</i> the hypothesis in all coverage-1.0 and tiered strata and "
            "<i>for</i> it in the two strata carrying the largest τ. Under the pre-specified "
            "control AUC_frac_usable (depths with n_eval ≥ 30), <b>P1 falls 9/9 → 8/9 and P2 "
            "9/9 → 7/9</b>. Weakened, not overturned. Per Appendix F.3 the conservative "
            "figure is the one quoted in any summary claim.",
            st["body"],
        ),
        PageBreak(),
    ]


def _phase1_p5(st: dict[str, Any]) -> list[Any]:
    return [
        para("2.6 P5 — the units question", st["h2"]),
        para(
            "r_val counts BFS hops; corruption counts claims. These are not the same axis, "
            "and radius_local_up_fast cannot measure the distance because it searches "
            "<i>upward</i> only while a flip-corrupted state is not a refinement of G₀. "
            "Exact hop distance therefore required full space enumeration on a "
            "small-component subsample (1,256 instances, component size 3–7; spaces of "
            "32–607 elements, 0.007–12.2 s to build).",
            st["body"],
        ),
        para(
            "One corrupted claim ≈ <b>2.9 hops</b> (sd 1.01), converging to exactly 2·d at "
            "d ≥ 5. <b>P5 passes in 6/6 strata</b>, signs agreeing, τ positive on both axes. "
            "Cross-checked against the main dataset: 6/6 sign agreement, 4/6 rough magnitude "
            "agreement; the two divergent cells are plausibly a component-size effect (main "
            "spans 6–12, subsample 3–7) — reported, not resolved.",
            st["body"],
        ),
        _radius_tag(st),
    ]


def _phase1_unpredicted(st: dict[str, Any]) -> list[Any]:
    return [
        para("2.7 The finding that was not predicted", st["h2"]),
        para("Contradiction rates by corruption depth, <b>flip arm only</b>:", st["body"]),
        table(
            [
                ["depth", "1", "2", "3", "5", "8", "12"],
                ["rate", "0.627", "0.777", "0.876", "0.892", "0.983", "1.000"],
            ],
            [2.2 * cm] + [2.33 * cm] * 6,
            align_right=[1, 2, 3, 4, 5, 6],
        ),
        para(
            "At d = 1 and base wrongness 0.00/0.10/0.25 the rates are 0.653/0.647/0.563.",
            st["body"],
        ),
        para(
            "<i>(Corrected 2026-09-10. This table originally printed "
            "0.491/0.773/0.896/0.930/0.986/1.000 under a “flip arm” label. Those "
            "are the <b>both-arms pooled</b> figures, and the pool is computed over the "
            "tiered arm's defective d column — the exact cross-arm merge forbidden "
            "everywhere else in this report. Caught while building the session figures, "
            "which recompute every plotted number from source. Flip-only figures are shown "
            "above; the base-wrongness breakdown was already flip-only and is unchanged. "
            "See Appendix G.)</i>",
            st["callout"],
        ),
        para(
            "<b>Most orientation errors are self-revealing.</b> Reversing a single truthful "
            "claim renders the knowledge set inconsistent with the CPDAG about 6 times in "
            "10. Meek closure fails and the analyst finds out.",
            st["good"],
        ),
        para(
            "And among corrupted sets that <i>stay</i> consistent, median survival is "
            "<b>0.952</b> — most survivors are harmless. The picture is a two-stage filter: "
            "most errors announce themselves, most of the rest are benign, and the dangerous "
            "case is <b>consistent <i>and</i> invalidating</b>. That is precisely the case "
            "r_val is defined over, since its space contains only consistent MPDAGs.",
            st["body"],
        ),
        *_figure(
            "session7_f5_self_revealing",
            "<b>Figure 5.</b> Contradiction rate by depth, flip arm only, after the pooling "
            "correction above.",
            st,
        ),
        para(
            "This also dissolves an apparent tension. r_val = 1 is common while average "
            "survival at d = 1 is high. Both are true: r_val is a <b>worst-case</b> distance "
            "to the nearest failing graph, S(d) an <b>average-case</b> probability over a "
            "shell. A failure at distance 1 says nothing about how much of the shell fails. "
            "That the worst case predicts the average case is the finding, not a tautology.",
            st["body"],
        ),
        PageBreak(),
    ]


def _phase2(st: dict[str, Any]) -> list[Any]:
    return [
        para("3. Phase 2 — there is no Pareto frontier", st["h1"]),
        para("3.1 Falsified twice", st["h2"]),
        para(
            "Round 1 (24 base CPDAGs, 1,087 truthful proposals, 17,680 SEM draws): strata "
            "with more than one distinct utility_median: <b>0 of 24</b>, while r_val varied "
            "in 18. Round 2 redefined utility as <i>achievable</i> efficiency over a capped "
            "12-set menu: <b>0 of 32</b> strata varied in achievable utility, in the count of "
            "valid menu members, or in the winning member — while r_val varied in 23. The "
            "pre-registration's trigger was explicit, so Phase 2 was <b>abandoned rather "
            "than redefined a third time</b>.",
            st["callout"],
        ),
        para("3.2 Why — and this is the substantive result", st["h2"]),
        para(
            "Every truthful proposal leaves the true DAG inside [G₀]. "
            "optimal_adjustment_set_mpdag returns non-None exactly when all extensions "
            "agree; if they agree, the common answer <i>is</i> the true DAG's optimal set. "
            "Every identifying truthful proposal therefore yields the same Z, hence the same "
            "asymptotic variance.",
            st["good"],
        ),
        para(
            "<b>For truthful background knowledge, statistical efficiency is an invariant of "
            "(Ĉ, x, y) — not a function of what the analyst asserts. Knowledge moves the "
            "robustness axis and only the robustness axis.</b>",
            st["body"],
        ),
        para(
            "So the trade-off is not precision against robustness. Precision is <i>fixed</i>. "
            "What knowledge purchases is <b>identifiability</b> — 385 of 1,431 proposals "
            "(26.9%) could justify no valid set at all. The operative trade-off is "
            "<b>identifiability against robustness, and it is a step, not a frontier</b>: "
            "assert enough to identify, and every further claim is pure r_val cost at zero "
            "precision gain.",
            st["body"],
        ),
        para(
            "This converges with Phase 1's P3 from the opposite direction: there, more "
            "claims measurably <i>reduced</i> survival. Two independent designs, same "
            "conclusion — marginal truthful background knowledge is not free, and its price "
            "is paid entirely in robustness.",
            st["body"],
        ),
        para(
            "Testing the step hypothesis properly — does r_val fall monotonically with "
            "proposal size beyond the minimum identifying set? — is a separate question and "
            "is <b>not claimed here</b>.",
            st["body"],
        ),
        para("3.3 P7′ — the census stands even though the contrast cannot be run", st["h2"]),
        para(
            "Round 1 found <b>zero</b> “aggressive” proposals, an artifact of "
            "component_generator, which builds the undirected component purely from X's "
            "upstream confounders and keeps X→Y and its descendants outside it by "
            "construction.",
            st["body"],
        ),
        para(
            "On real structure the picture is different. Sampling ≤200 descendant-ordered "
            "pairs per network across the 39-network corpus (M-bias excluded as an ADMG), "
            "<b>28 networks admit at least one (x, y) pair with an undirected CPDAG edge on "
            "a directed causal path — 396 pairs</b>, densest in paths (96), Polzer_2012 (64), "
            "Kampen_2014 (33), child (22), sachs (22), insurance (21). Eight real-network "
            "bases were gated and yielded 215 aggressive and 129 conservative proposals, so "
            "the empty cell is genuinely filled — but utility is invariant there too (0 of 8 "
            "strata vary), so the intended contrast cannot be run. The residual r_val-only "
            "contrast is non-directional (2 bases lower for aggressive, 2 tied, 1 reversed) "
            "and is <b>inconclusive</b>.",
            st["body"],
        ),
        _radius_tag(st),
        PageBreak(),
    ]


def _not_claimed(st: dict[str, Any]) -> list[Any]:
    return [
        para("4. What is not claimed", st["h1"]),
        *bullets(
            [
                "<b>No cross-arm detectability comparison.</b> An earlier draft claimed "
                "correlated errors are harder to detect than uniform ones (0.354 vs 0.653 "
                "at d = 1). The tiered half came from the contaminated bin, and the "
                "comparison was never sound in principle: one arm relocates a fraction of "
                "<i>nodes</i>, the other reverses a count of <i>claims</i>. Retracted in "
                "Appendix E.3. Establishing it needs a corruption process parameterised "
                "identically across both arms — a design change, not a reanalysis.",
                "<b>No claim of perfect ranking.</b> The endpoint is a rank correlation with "
                "a bootstrap interval. One stratum is null.",
                "<b>Synthetic structure only</b> for Phase 1. Real networks appear only in "
                "the P7′ census.",
                "<b>Oracle CI throughout.</b> No finite-sample discovery; causal sufficiency "
                "assumed. These remain the project's largest scope limits and this work does "
                "not touch them.",
                "<b>AUC_frac is not quoted alone.</b> The conservative AUC_frac_usable "
                "figure governs summary claims (Appendix F.3).",
            ],
            st,
        ),
    ]


def _provenance(st: dict[str, Any]) -> list[Any]:
    return [
        para("5. Provenance", st["h1"]),
        para(
            "All results under <font face='%FONTMONO%'>results/axis_robustness/</font>: "
            "<font face='%FONTMONO%'>survival_instances.csv</font> (1,978), "
            "<font face='%FONTMONO%'>survival_curves.csv</font> (14,066), "
            "<font face='%FONTMONO%'>survival_samples.csv.gz</font> (3,136,000 rows, 48 MB "
            "compressed), <font face='%FONTMONO%'>hops_samples.csv.gz</font> (6 MB), "
            "<font face='%FONTMONO%'>analysis_tau_primary.csv</font>, "
            "<font face='%FONTMONO%'>analysis_tau_secondary.csv</font>, "
            "<font face='%FONTMONO%'>analysis_tiered_rebinned.csv</font> (668), "
            "<font face='%FONTMONO%'>analysis_decomposition.csv</font>, "
            "<font face='%FONTMONO%'>analysis_discordance.csv</font>, "
            "<font face='%FONTMONO%'>analysis_effect_sizes.csv</font>, "
            "<font face='%FONTMONO%'>pareto_proposals.csv</font> (1,087), "
            "<font face='%FONTMONO%'>pareto_sems.csv</font> (17,680), "
            "<font face='%FONTMONO%'>pareto2_proposals.csv</font> (1,431), "
            "<font face='%FONTMONO%'>pareto2_menus.csv</font> (243), "
            "<font face='%FONTMONO%'>pareto2_aggressive_census.csv</font> (39), "
            "<font face='%FONTMONO%'>hops_p5.csv</font>, "
            "<font face='%FONTMONO%'>hops_instances.csv</font>, "
            "<font face='%FONTMONO%'>hops_cost_census.csv</font>, plus manifests capturing "
            "git SHA, interpreter, platform, library versions and RADIUS_CONVENTION.",
            st["body"],
        ),
        para(
            "<font face='%FONTMONO%'>num_workers: 1</font>, "
            "<font face='%FONTMONO%'>random_seed: 0</font>, CP-SAT single-threaded, no "
            "global RNG anywhere in the live tree. Determinism verified across "
            "<font face='%FONTMONO%'>PYTHONHASHSEED</font> 0 / 12345 for every stage; "
            "bootstrap CIs bit-identical across independent runs. The decomposition identity "
            "S_contra_as_fail(d) = (1 − contradiction_rate(d)) · S(d) holds to 2.2e-16 on all "
            "8,841 rows where S is defined.",
            st["body"],
        ),
        para(
            "<b>Known incident.</b> The first survival sweep was destroyed when its worker "
            "relaunched and overwrote its own outputs in place; it was regenerated (run 2) "
            "and reproduces the original's pooled statistics exactly. The relaunch also "
            "fixed a real defect — run 1's instance_id was not unique — so several "
            "Appendix C statistics computed by grouping on it were wrong and are corrected "
            "in Appendix D, including a retracted endpoint swap. Full account in "
            "Appendices C–D.",
            st["callout"],
        ),
        para("Environment", st["h2"]),
        para(
            "Python 3.9.6, reportlab 5.0.0. Every table above is transcribed verbatim from "
            "<font face='%FONTMONO%'>report_fragility_and_pareto.md</font> (committed at "
            "<font face='%FONTMONO%'>ed9270d</font>) rather than re-derived, for the reasons "
            "given in this module's docstring; no CSV, the markdown, or any figure was "
            "modified to build this document.",
            st["body"],
        ),
    ]


def _appendix(st: dict[str, Any]) -> list[Any]:
    commits = _session_commits()
    story: list[Any] = [
        PageBreak(),
        para("Appendix — supplementary figure", st["h1"]),
        para(
            "Not referenced by any numbered section above; included for completeness as a "
            "per-depth breakdown underlying the whole-curve endpoint used throughout §2.",
            st["body"],
        ),
        *_figure(
            "session7_f6_survival_curves_by_bucket",
            "<b>Figure 6 (supplementary).</b> Survival curves S(d) by r_val bucket, per "
            "corruption depth — the per-depth view that Figure 1's whole-curve AUC_frac "
            "endpoint summarises. Supplementary: see §2.2 for the primary evidence. "
            + RADIUS_NOTE,
            st,
        ),
    ]
    if commits:
        story += [
            Spacer(1, 0.3 * cm),
            para("Commits in this series", st["h2"]),
            table(
                [["commit", "date", "subject"], *commits],
                [2.2 * cm, 2.2 * cm, 11.8 * cm],
                font_size=7,
            ),
        ]
    return story


def build(out_path: str | Path = "REPORT_FRAGILITY_AND_PARETO.pdf") -> Path:
    """Assemble the session-7 PDF. Returns the path written."""
    st = styles()
    doc = SimpleDocTemplate(
        str(out_path),
        pagesize=A4,
        leftMargin=2.2 * cm,
        rightMargin=2.2 * cm,
        topMargin=1.8 * cm,
        bottomMargin=1.8 * cm,
        title="Robustness survival and the Pareto frontier of causal assumptions",
        author="Anonymous",
    )
    story: list[Any] = []
    story += _cover(st)
    story += _phase1_design(st)
    story += _phase1_p1(st)
    story += _phase1_p2p3(st)
    story += _phase1_discordance(st)
    story += _phase1_defects(st)
    story += _phase1_p5(st)
    story += _phase1_unpredicted(st)
    story += _phase2(st)
    story += _not_claimed(st)
    story += _provenance(st)
    story += _appendix(st)

    def footer(canvas: Any, doc_: Any) -> None:
        canvas.saveState()
        canvas.setFont("Helvetica", 7.5)
        canvas.setFillColor(GREY)
        canvas.drawString(
            2.2 * cm,
            1.1 * cm,
            "Robustness survival and the Pareto frontier of causal assumptions (session 7)",
        )
        canvas.drawRightString(A4[0] - 2.2 * cm, 1.1 * cm, f"page {doc_.page}")
        canvas.restoreState()

    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    return Path(out_path)


if __name__ == "__main__":
    print(build())
