"""The matched-population coverage contrast ([RE-6]), on the 105 triples.

Admissible instances collapse 543 -> 182 -> 106 as knowledge coverage falls, and
the survivors have systematically larger separation, so **any marginal rate
compared across coverage levels measures composition, not coverage**
(``report_real_graphs.md`` §5.3, ``docs/REMAINING_EXPERIMENTS.md`` [RE-4],
[RE-6]). This module therefore conditions on the **fixed instance set**: the
``(network, X, Y)`` triples that are admissible at *every* coverage level the
corpus carries, of which there are 105 on 8 networks.

It is **supplementary**. It produces no τ against the anchor table's strata and
it is labelled as a matched-population contrast wherever it is reported.

The sweep's other known confound applies here too and is carried through: the
coverage sweep is **not nested** ([RE-4]) -- ``select_knowledge`` uses a stride
that selects different index sets at different coverages, so ``K(0.25)`` is not
in general a subset of ``K(0.5)``. Conditioning on a fixed instance set removes
the *composition* effect, not the *nesting* one, and the second is reported
rather than fixed, because fixing it would change every committed number in
``results/axisa3/``.
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any

COLUMNS = [
    "coverage", "n_triples", "n_networks", "median_n_k", "median_k_g0",
    "median_leverage", "median_radius", "share_radius_1", "max_radius",
    "median_AUC_frac_usable", "median_AUC_frac_usable_network_weighted",
    "share_separation_measured", "median_separation",
]


def matched_triples(units: list[dict[str, Any]]) -> set[tuple[str, str, str]]:
    """The ``(network, X, Y)`` triples present at every coverage level.

    Args:
        units: Rows of ``analysis_units.csv`` for the flip arm.

    Returns:
        The matched triple set.
    """
    covs: dict[tuple[str, str, str], set[float]] = defaultdict(set)
    all_covs: set[float] = set()
    for u in units:
        cov = float(u["coverage"])
        all_covs.add(cov)
        covs[(u["network"], u["x"], u["y"])].add(cov)
    return {t for t, c in covs.items() if c == all_covs}


def _med(vals: list[float]) -> float | None:
    """Median, or ``None`` for an empty list.

    Args:
        vals: The values.

    Returns:
        The median.
    """
    return statistics.median(vals) if vals else None


def contrast(units: list[dict[str, Any]], triples: set[tuple[str, str, str]]) -> list[dict[str, Any]]:
    """One row per coverage level, over the matched triples only.

    Args:
        units: Flip-arm rows of ``analysis_units.csv`` at base wrongness 0.
        triples: The matched triple set.

    Returns:
        The contrast rows, ascending in coverage.
    """
    by_cov: dict[float, list[dict[str, Any]]] = defaultdict(list)
    for u in units:
        if (u["network"], u["x"], u["y"]) in triples:
            by_cov[float(u["coverage"])].append(u)

    rows = []
    for cov in sorted(by_cov):
        us = by_cov[cov]
        radii = [int(float(u["radius"])) for u in us if u["radius"] not in ("", "None")]
        aucs = [float(u["AUC_frac_usable"]) for u in us if u["AUC_frac_usable"] not in ("", "None")]
        per_net: dict[str, list[float]] = defaultdict(list)
        for u in us:
            if u["AUC_frac_usable"] not in ("", "None"):
                per_net[u["network"]].append(float(u["AUC_frac_usable"]))
        seps = [float(u["separation"]) for u in us if u["separation"] not in ("", "None")]
        nk = [float(u["n_k"]) for u in us]
        kg = [float(u["k_g0"]) for u in us]
        rows.append({
            "coverage": cov,
            "n_triples": len(us),
            "n_networks": len({u["network"] for u in us}),
            "median_n_k": _med(nk),
            "median_k_g0": _med(kg),
            "median_leverage": _med([k / n for k, n in zip(kg, nk) if n]),
            "median_radius": _med([float(r) for r in radii]),
            "share_radius_1": (sum(1 for r in radii if r == 1) / len(radii)) if radii else None,
            "max_radius": max(radii) if radii else None,
            "median_AUC_frac_usable": _med(aucs),
            "median_AUC_frac_usable_network_weighted": _med(
                [statistics.median(v) for v in per_net.values()]
            ),
            "share_separation_measured": (
                sum(1 for u in us if u["separation_status"] == "measured") / len(us) if us else None
            ),
            "median_separation": _med(seps),
        })
    return rows


def main(argv: list[str] | None = None) -> int:
    """Entry point.

    Args:
        argv: Argument vector, or ``None`` for ``sys.argv``.

    Returns:
        Process exit status.
    """
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in-dir", default="results/axis_robustness_real")
    args = p.parse_args(argv)
    in_dir = Path(args.in_dir)

    with (in_dir / "analysis_units.csv").open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    flip0 = [
        u for u in rows
        if u["arm"] == "flip" and u["status"] == "ok"
        and u["base_wrongness"] not in ("", "None") and float(u["base_wrongness"]) == 0.0
    ]
    triples = matched_triples(flip0)
    out = contrast(flip0, triples)

    path = in_dir / "matched_coverage_contrast.csv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        for r in out:
            w.writerow(r)

    summary = {
        "n_matched_triples": len(triples),
        "n_networks": len({t[0] for t in triples}),
        "networks": sorted({t[0] for t in triples}),
        "unmatched_note": (
            "Admissible instances collapse 543 -> 182 -> 106 across coverage, so "
            "the unmatched marginal comparison is a composition effect and is not "
            "reported. The coverage sweep is also not nested ([RE-4]); conditioning "
            "on a fixed instance set removes the composition effect, not the "
            "nesting one."
        ),
        "rows": out,
    }
    (in_dir / "matched_coverage_contrast.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
