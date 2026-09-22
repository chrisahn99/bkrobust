"""Within- versus between-network decomposition of τ_b(`r_hop`, survival).

**POST-HOC.** Not pre-registered. The pre-registered primary result is the
*pooled* τ of `analysis_tau.csv`, and this module does not replace it.

Why it exists. `PREREGISTRATION.md` §1.5 fixed the **network** as the unit of
analysis and bootstrapped over networks accordingly — but the τ itself was still
computed *pooled over rows drawn from different networks*. That statistic answers
"does a higher radius predict higher survival across all queries on all
networks", which mixes two quite different comparisons:

* **within-network** — same graph, same knowledge set `K`, same `G₀`, different
  `(X, Y)` query. This is the comparison a practitioner actually faces, and it is
  the one the radius is defined for: it prices *this* conclusion on *this*
  knowledge state.
* **between-network** — different graphs with different `|K|`, different
  component structure, and the small-`|K|` endpoint pinning of Appendix F. The
  radius is confounded with all three.

Where those two components disagree, the pooled figure is a mixture of them and
is hard to read. On this corpus they do disagree, and in the flip arm the
between-network component is the one dragging the pooled figure down.

The statistic reported here is the **mean of the per-network τ**, bootstrapped by
resampling networks with replacement — the same cluster scheme as everywhere
else. A network contributes only if **both** the radius and the endpoint vary
within it; the count of contributing networks is reported, because where it is
small the statistic is weak however many rows sit behind it.
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import kendalltau

from bkrobust.core.conventions import UNREACHED

COLUMNS = [
    "stratum", "stratum_class", "arm", "predictor", "endpoint",
    "pooled_tau_b", "n_rows", "n_networks_total",
    "n_networks_with_variation", "n_networks_positive",
    "within_mean_tau", "within_median_tau", "within_ci_lo_2p5", "within_ci_hi_97p5",
    "within_excludes_zero", "between_tau_b",
]


def decompose(
    units: list[dict[str, Any]], predictor: str, endpoint: str,
    *, n_boot: int = 10000, seed: int = 0,
) -> dict[str, Any] | None:
    """Split a stratum's rank association into within- and between-network parts.

    Args:
        units: Analysis-unit rows for one stratum.
        predictor: Predictor column.
        endpoint: Endpoint column.
        n_boot: Bootstrap resamples over networks.
        seed: Seed for ``np.random.default_rng``.

    Returns:
        The decomposition, or ``None`` if fewer than three networks carry
        variation in both the predictor and the endpoint.
    """
    sel = [
        u for u in units
        if u[predictor] not in ("", "None") and u[endpoint] not in ("", "None")
        and not (predictor == "radius" and float(u[predictor]) == UNREACHED)
    ]
    if len(sel) < 3:
        return None
    x = np.array([float(u[predictor]) for u in sel])
    y = np.array([float(u[endpoint]) for u in sel])
    pooled = float(kendalltau(x, y, variant="b")[0]) if np.unique(x).size > 1 else float("nan")

    per: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for u in sel:
        per[u["network"]].append(u)

    taus: dict[str, float] = {}
    for net, rows in per.items():
        xx = np.array([float(u[predictor]) for u in rows])
        yy = np.array([float(u[endpoint]) for u in rows])
        if np.unique(xx).size > 1 and np.unique(yy).size > 1:
            t, _ = kendalltau(xx, yy, variant="b")
            if not np.isnan(t):
                taus[net] = float(t)
    if len(taus) < 3:
        return None

    nets = sorted(taus)
    vals = np.array([taus[n] for n in nets])
    rng = np.random.default_rng(seed)
    means = np.array([
        vals[rng.integers(0, len(nets), size=len(nets))].mean() for _ in range(n_boot)
    ])
    lo, hi = np.percentile(means, [2.5, 97.5])

    med_x = [statistics.median(float(u[predictor]) for u in rows) for rows in per.values()]
    med_y = [statistics.median(float(u[endpoint]) for u in rows) for rows in per.values()]
    between = (
        float(kendalltau(np.array(med_x), np.array(med_y), variant="b")[0])
        if len(set(med_x)) > 1 else float("nan")
    )

    return {
        "pooled_tau_b": pooled, "n_rows": len(sel), "n_networks_total": len(per),
        "n_networks_with_variation": len(taus),
        "n_networks_positive": int((vals > 0).sum()),
        "within_mean_tau": float(vals.mean()),
        "within_median_tau": float(np.median(vals)),
        "within_ci_lo_2p5": float(lo), "within_ci_hi_97p5": float(hi),
        "within_excludes_zero": bool(lo > 0 or hi < 0),
        "between_tau_b": between,
    }


def main(argv: list[str] | None = None) -> int:
    """Entry point.

    Args:
        argv: Argument vector, or ``None`` for ``sys.argv``.

    Returns:
        Process exit status.
    """
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in-dir", default="results/axis_robustness_real")
    p.add_argument("--n-boot", type=int, default=10000)
    args = p.parse_args(argv)
    in_dir = Path(args.in_dir)

    with (in_dir / "analysis_units.csv").open(newline="") as fh:
        units = [u for u in csv.DictReader(fh) if u["status"] == "ok"]
    by_stratum: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for u in units:
        by_stratum[u["stratum"]].append(u)

    rows: list[dict[str, Any]] = []
    for stratum in sorted(by_stratum):
        us = by_stratum[stratum]
        arm = us[0]["arm"]
        endpoint = "AUC_rate_usable" if arm == "tiered" else "AUC_frac_usable"
        for predictor in ("radius", "shd_truth", "n_k"):
            d = decompose(us, predictor, endpoint, n_boot=args.n_boot)
            if d is None:
                continue
            rows.append({
                "stratum": stratum, "stratum_class": us[0]["stratum_class"],
                "arm": arm, "predictor": predictor, "endpoint": endpoint, **d,
            })

    path = in_dir / "analysis_within_between.csv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    prim9 = {
        "flip_cov050_bw000", "flip_cov050_bw010", "flip_cov050_bw025",
        "flip_cov100_bw000", "flip_cov100_bw010", "flip_cov100_bw025",
        "tiered_nt2", "tiered_nt3", "tiered_nt4",
    }
    rad = [r for r in rows if r["predictor"] == "radius"]
    summary = {
        "note": "POST-HOC. The pre-registered primary result is the pooled tau in analysis_tau.csv.",
        "n_rows": len(rows),
        "radius_primary9_within_excludes_zero": sorted(
            r["stratum"] for r in rad if r["stratum"] in prim9 and r["within_excludes_zero"]
        ),
        "radius_primary9_with_a_within_estimate": sorted(
            r["stratum"] for r in rad if r["stratum"] in prim9
        ),
        "radius_all_strata_within_excludes_zero": sorted(
            r["stratum"] for r in rad if r["within_excludes_zero"]
        ),
    }
    (in_dir / "analysis_within_between.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
