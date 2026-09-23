"""Build ``table_tau_real.md`` -- the real-structure counterpart of ``table_tau_comparisons.md``.

Read-only against ``results/axis_robustness_real/``. It never runs a sweep and
never recomputes a statistic: every number it prints is copied from
``analysis_tau.csv`` and ``analysis_units.csv``, which
:mod:`bkrobust.robustness.real_analyse` produced. The point of keeping the
formatting separate from the statistics is that the verifier can re-assert the
table's numbers against the CSVs without going through this module at all.

Schema and footnoting discipline follow ``table_tau_comparisons.md``:

* the nine pre-registered strata, in the same order as their synthetic
  counterparts, so the two tables can be read side by side -- **and never
  pooled, and never compared row to row with a synthetic row**;
* ``n`` and, because the unit of analysis here is the network rather than the
  pair, ``n_networks`` beside it;
* the endpoint actually used, named in every row;
* τ_b with a 95% interval **bootstrapped over networks**, for the radius and for
  each baseline;
* a marker on every cell where a predictor is undefined by construction, on
  every cell whose verdict does not survive dropping a single network, and on
  every cell where the raw and conservative endpoints disagree;
* the assumption string carried by every radius row.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

DEFAULT_DIR = Path("results/axis_robustness_real")
DEFAULT_TABLE = Path("table_tau_real.md")

#: The nine pre-registered strata, keyed by the tag ``real_analyse`` writes, in
#: the order their synthetic counterparts appear in
#: ``table_tau_comparisons.md`` so the two tables line up row for row.
PRIMARY_ORDER: tuple[str, ...] = (
    "flip_cov050_bw000",
    "flip_cov050_bw010",
    "flip_cov050_bw025",
    "flip_cov100_bw000",
    "flip_cov100_bw010",
    "flip_cov100_bw025",
    "tiered_nt2",
    "tiered_nt3",
    "tiered_nt4",
)

#: Human-readable labels, so the table reads like the synthetic one rather than
#: like a column of internal tags.
STRATUM_LABELS: dict[str, str] = {
    "flip_cov050_bw000": "flip, coverage=0.5, base_wrongness=0.00",
    "flip_cov050_bw010": "flip, coverage=0.5, base_wrongness=0.10",
    "flip_cov050_bw025": "flip, coverage=0.5, base_wrongness=0.25",
    "flip_cov100_bw000": "flip, coverage=1.0, base_wrongness=0.00",
    "flip_cov100_bw010": "flip, coverage=1.0, base_wrongness=0.10",
    "flip_cov100_bw025": "flip, coverage=1.0, base_wrongness=0.25",
    "tiered_nt2": "tiered, n_tiers=2",
    "tiered_nt3": "tiered, n_tiers=3",
    "tiered_nt4": "tiered, n_tiers=4",
    "flip_cov025_bw000": "flip, coverage=0.25, base_wrongness=0.00",
    "flip_cov025_bw010": "flip, coverage=0.25, base_wrongness=0.10",
    "flip_cov025_bw025": "flip, coverage=0.25, base_wrongness=0.25",
    "flip_cov025_bwa1": "flip, coverage=0.25, one claim reversed",
    "flip_cov050_bwa1": "flip, coverage=0.5, one claim reversed",
    "flip_cov100_bwa1": "flip, coverage=1.0, one claim reversed",
}


def label_for(stratum: str) -> str:
    """The human-readable label for a stratum tag.

    Args:
        stratum: The tag ``real_analyse`` writes.

    Returns:
        The label, or the tag itself if it has none.
    """
    return STRATUM_LABELS.get(stratum, stratum)

#: Column order for the predictors. ``r_hop`` is the legacy ``r_val``; the paper
#: name is used in the table and the code name in the files.
PREDICTORS: tuple[tuple[str, str], ...] = (
    ("radius", "r_hop"),
    ("shd_truth", "shd_truth"),
    ("n_k", "|K|"),
    ("k_g0", "k_g0"),
)

ASSUMPTION = "Conjecture 2 (proved: Anti-Exchange Case B, THEOREMS.md section 4)"


def read_csv(path: Path) -> list[dict[str, str]]:
    """Read a CSV into a list of dicts.

    Args:
        path: The file.

    Returns:
        The rows.

    Raises:
        FileNotFoundError: If the file is absent.
    """
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def _f(value: str | None) -> float | None:
    """Parse a CSV cell as a float, treating blanks and ``None`` as missing.

    Args:
        value: The raw cell.

    Returns:
        The float, or ``None``.
    """
    if value is None or value == "" or value.lower() in {"none", "nan"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _i(value: str | None) -> int | None:
    """Parse a CSV cell as an int.

    Args:
        value: The raw cell.

    Returns:
        The int, or ``None``.
    """
    f = _f(value)
    return None if f is None else int(f)


def format_cell(row: dict[str, str] | None, marks: list[str]) -> str:
    """Render one τ cell, or the reason there is no number in it.

    A predictor that is constant in a stratum prints
    ``undefined (predictor constant)`` and **never** ``0.000`` -- the two mean
    entirely different things and the synthetic table makes the same
    distinction.

    Args:
        row: The ``analysis_tau.csv`` row, or ``None`` if the stratum has none.
        marks: Footnote markers to append.

    Returns:
        The markdown cell.
    """
    if row is None:
        return "not run"
    status = row.get("status", "")
    if status.startswith("undefined"):
        return f"{status}[^const]"
    if status in {"insufficient_n", "bootstrap_degenerate"}:
        return f"{status}"
    tau, lo, hi = _f(row.get("tau_b")), _f(row.get("ci_lo_2p5")), _f(row.get("ci_hi_97p5"))
    if tau is None:
        return status or "—"
    body = f"{tau:+.3f}" if lo is None else f"{tau:+.3f} [{lo:+.3f}, {hi:+.3f}]"
    return body + "".join(marks)


def endpoint_for(stratum: str) -> str:
    """The conservative endpoint name for a stratum's arm.

    Args:
        stratum: The stratum label.

    Returns:
        ``AUC_rate_usable`` for the tiered arm, ``AUC_frac_usable`` otherwise.
    """
    return "AUC_rate_usable" if stratum.startswith("tiered") else "AUC_frac_usable"


def index_tau(rows: list[dict[str, str]]) -> dict[tuple[str, str, str], dict[str, str]]:
    """Index ``analysis_tau.csv`` by ``(stratum, predictor, endpoint)``.

    Args:
        rows: The τ rows.

    Returns:
        The index.
    """
    return {(r["stratum"], r["predictor"], r["endpoint"]): r for r in rows}


def build_rows(
    tau_rows: list[dict[str, str]], strata: list[str]
) -> tuple[list[list[str]], dict[str, bool]]:
    """Build the markdown body for a set of strata.

    Args:
        tau_rows: Every row of ``analysis_tau.csv``.
        strata: The strata to render, in order.

    Returns:
        ``(rows, flags)`` where ``flags`` says which footnotes were used.
    """
    idx = index_tau(tau_rows)
    flags = {"const": False, "loo": False, "gap": False, "censored": False}
    out: list[list[str]] = []
    for stratum in strata:
        endpoint = endpoint_for(stratum)
        raw_endpoint = endpoint.replace("_usable", "")
        cells: list[str] = []
        n = n_net = None
        for code, _label in PREDICTORS:
            row = idx.get((stratum, code, endpoint))
            marks: list[str] = []
            if row is not None:
                if n is None and _i(row.get("n")) is not None:
                    n, n_net = _i(row.get("n")), _i(row.get("n_networks"))
                if str(row.get("loo_verdict_flips", "")).lower() == "true":
                    marks.append("[^loo]")
                    flags["loo"] = True
                raw = idx.get((stratum, code, raw_endpoint))
                if raw is not None and _disagree(row, raw):
                    marks.append("[^gap]")
                    flags["gap"] = True
                if str(row.get("status", "")).startswith("undefined"):
                    flags["const"] = True
            cells.append(format_cell(row, marks))
        out.append([label_for(stratum), str(n if n is not None else "—"),
                    str(n_net if n_net is not None else "—"), f"`{endpoint}`", *cells])
    return out, flags


def _disagree(usable: dict[str, str], raw: dict[str, str]) -> bool:
    """Whether the conservative and raw endpoints give different verdicts.

    "Different verdict" means a different sign, or one interval excluding zero
    while the other covers it. Appendix J of the synthetic pre-registration
    established that a gap between the two endpoints is a symptom of too few
    draws rather than evidence that either is safer, so a disagreement at
    N = 1000 is worth a marker rather than a choice.

    Args:
        usable: The τ row on the conservative endpoint.
        raw: The τ row on the raw endpoint.

    Returns:
        ``True`` if the verdicts differ.
    """
    tu, tr = _f(usable.get("tau_b")), _f(raw.get("tau_b"))
    if tu is None or tr is None:
        return False
    if (tu > 0) != (tr > 0):
        return True
    def excl(row: dict[str, str]) -> bool | None:
        lo, hi = _f(row.get("ci_lo_2p5")), _f(row.get("ci_hi_97p5"))
        return None if lo is None or hi is None else (lo > 0 or hi < 0)
    eu, er = excl(usable), excl(raw)
    return eu is not None and er is not None and eu != er


def markdown(
    tau_rows: list[dict[str, str]], summary: dict[str, Any]
) -> str:
    """Render the whole table file.

    Args:
        tau_rows: Every row of ``analysis_tau.csv``.
        summary: ``analysis_summary.json``.

    Returns:
        The markdown document.
    """
    present = {r["stratum"] for r in tau_rows}
    primary = [s for s in PRIMARY_ORDER if s in present]
    supplementary = sorted(
        s for s in present
        if s not in PRIMARY_ORDER
    )
    body, flags = build_rows(tau_rows, primary)
    sup_body, sup_flags = build_rows(tau_rows, supplementary)
    for k, v in sup_flags.items():
        flags[k] = flags[k] or v

    header = ("| stratum | n | networks | endpoint used | "
              + " | ".join(f"τ_b({label}) [95% CI]" for _c, label in PREDICTORS) + " |")
    sep = "|" + "---|" * (4 + len(PREDICTORS))

    lines: list[str] = []
    lines.append("# Real-structure anchor table — τ_b(predictor, survival AUC) by stratum")
    lines.append("")
    lines.append(
        "The real-structure counterpart of `table_tau_comparisons.md`, on the **committed "
        "real-network corpus**: the 831 admissible rows on 25 networks in "
        "`results/axisa3/instances.jsonl`. Same schema, same footnoting discipline. "
        "**The two tables are never pooled and no row here may be compared with a "
        "synthetic row** — the corpora differ in structure, in how knowledge is obtained, "
        "and in what an instance is."
    )
    lines.append("")
    lines.append(
        "Three things about this corpus change how the table must be read, and all three "
        "were measured before any τ existed (`results/axis_robustness_real/PREREGISTRATION.md` "
        "§1.5, §2.1, Appendix A):"
    )
    lines.append("")
    lines.append(
        "1. **The unit of analysis is the network, not the pair.** Every interval here is a "
        "cluster bootstrap over the 25 networks, 10,000 resamples, seed 0 — not over rows. "
        "The intraclass correlation of the `r = 1` indicator is 0.40, a design effect of "
        "about 12.7, so a row-level interval would be roughly three and a half times too "
        "narrow. `networks` is printed beside `n` for that reason, and no per-row interval "
        "appears anywhere."
    )
    lines.append(
        "2. **`|K|` and `k_g0` are constants within a `(network, coverage)` cell** on this "
        "corpus: `select_knowledge` returns one claim set per network and coverage, and "
        "`flip` preserves list length. Those two baselines can therefore rank *networks* "
        "but not *queries*, which is a weaker thing than they could do synthetically. A win "
        "over them here is worth correspondingly less, and that is stated rather than "
        "banked."
    )
    lines.append(
        "3. **`SHD(G₀, truth)` needs base wrongness to exist at all.** The corpus draws the "
        "analyst's claims from the true orientations, so `G₀` *is* the ground-truth DAG at "
        "coverage 1.0 and `SHD ≡ 0`. Base wrongness is applied before the sweep, exactly as "
        "the synthetic design did. Appendix A records that the rate-based mechanism reverses "
        "nothing at all on most of this corpus because `|K|` is 1 or 2 there, and adds an "
        "absolute one-claim level, reported in the supplementary table."
    )
    lines.append("")
    lines.append("## The nine pre-registered strata")
    lines.append("")
    lines.append(header)
    lines.append(sep)
    for row in body:
        lines.append("| " + " | ".join(row) + " |")
    lines.append("")
    if sup_body:
        lines.append("## Supplementary strata")
        lines.append("")
        lines.append(
            "Not among the nine, and not comparable with the synthetic table. Coverage 0.25 "
            "has no synthetic counterpart; the `bw_abs=1` strata exist to give `SHD(G₀, truth)` "
            "a non-degenerate comparison on every network (Appendix A). **No rate is compared "
            "across coverage levels**: admissible instances collapse 543 → 182 → 106 as "
            "coverage falls and the survivors have larger separation, so such a comparison "
            "measures composition ([RE-4], [RE-6])."
        )
        lines.append("")
        lines.append(header)
        lines.append(sep)
        for row in sup_body:
            lines.append("| " + " | ".join(row) + " |")
        lines.append("")

    lines.append("### Footnotes")
    lines.append("")
    if flags["const"]:
        lines.append(
            "[^const]: **Undefined by construction, not zero.** The predictor takes a single "
            "value across the whole stratum, so no rank correlation exists. Printed as "
            "`undefined (predictor constant)` rather than `τ = 0`, and never counted as "
            "evidence that the baseline \"fails\"."
        )
    if flags["loo"]:
        lines.append(
            "[^loo]: **⚠ single-network leverage — do not quote this cell on its own.** The "
            "cell's verdict does not survive removing one network: either some "
            "leave-one-network-out τ has the opposite sign to the full-sample τ, or the "
            "full-sample interval excludes zero while the leave-one-out range spans it. "
            "Within-network variation in `r_hop` is thin on this corpus — at coverage 1.0, "
            "10 of the 25 networks have a single distinct radius across all their rows, "
            "while `paths` alone spans `r = 1…14` in 20 rows — so a single network can "
            "carry a stratum. `analysis_tau.csv` names the network at the minimum and the "
            "maximum."
        )
    if flags["gap"]:
        lines.append(
            "[^gap]: **Raw and conservative endpoints disagree at N = 1000.** The "
            "`n_eval ≥ 30` filter conditions on a quantity correlated with the radius, so it "
            "is not a conservative control (synthetic `PREREGISTRATION.md` Appendix J). "
            "Appendix J also established the diagnostic: at N = 1000 the two endpoints "
            "converged to 0.0002 on the same synthetic instances, so a disagreement at "
            "N = 1000 here means **neither** figure is settled, and the stratum is reported "
            "**unresolved** rather than resolved by choosing an endpoint. Both figures are "
            "in `analysis_tau.csv`."
        )
    lines.append("")
    lines.append(
        f"**Assumption (every `r_hop` row):** {ASSUMPTION}; the error is one-sided, so "
        "radii can only be too large, never too small (`THEOREMS.md` §4c, §6, §8)."
    )
    lines.append("")
    lines.append("**Gate:** `benchmarks.measure.fast_gate` exclusively, on every row.")
    lines.append("")
    lines.append(
        f"**Draws:** N = {summary.get('n_draws_nominal', 1000)} per grid point. "
        "**Endpoint:** the conservative `_usable` variant is tabulated; the raw variant is "
        "in `analysis_tau.csv` and every disagreement is marked."
    )
    lines.append("")
    lines.append(
        "**Provenance:** `results/axis_robustness_real/analysis_tau.csv`, built by "
        "`bkrobust.robustness.real_analyse` from the sharded sweep, which was built on the "
        "frozen frame `results/axis_robustness_real/frame.jsonl` "
        f"(sha256 `{str(summary.get('frame', {}).get('frame_sha256', ''))[:16]}…`). Regenerate with "
        "`PYTHONPATH=src python -m bkrobust.robustness.build_table_tau_real`."
    )
    lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    """Entry point.

    Args:
        argv: Argument vector, or ``None`` for ``sys.argv``.

    Returns:
        Process exit status.
    """
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in-dir", default=str(DEFAULT_DIR))
    p.add_argument("--out", default=str(DEFAULT_TABLE))
    args = p.parse_args(argv)

    in_dir = Path(args.in_dir)
    tau_rows = read_csv(in_dir / "analysis_tau.csv")
    summary = json.loads((in_dir / "analysis_summary.json").read_text())
    text = markdown(tau_rows, summary)
    Path(args.out).write_text(text)
    print(json.dumps({"written": args.out, "n_tau_rows": len(tau_rows),
                      "n_strata": len({r["stratum"] for r in tau_rows})}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
