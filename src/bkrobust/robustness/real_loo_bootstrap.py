"""Leave-one-network-out **with the bootstrap re-run** — the stronger check.

`PREREGISTRATION.md` §6.2 specified a leave-one-network-out sensitivity on the
**point estimate only**. Appendix H records why that is too weak: it cannot see a
stratum whose interval excludes zero with every network and includes zero without
one, while every leave-one-out point estimate stays comfortably positive. All
three tiered strata pass the pre-registered flag and two of them fail this one.

This module is therefore **post-hoc**, is labelled as such everywhere it is
reported, and does not replace the pre-registered column — which is reported
unchanged beside it.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import kendalltau

from bkrobust.core.conventions import UNREACHED

COLUMNS = [
    "stratum", "stratum_class", "predictor", "endpoint", "dropped_network",
    "n", "n_networks", "tau_b", "ci_lo_2p5", "ci_hi_97p5", "excludes_zero",
    "full_tau_b", "full_excludes_zero", "loses_zero_exclusion",
]


def cluster_tau(
    rows: list[dict[str, Any]], predictor: str, endpoint: str,
    *, n_boot: int = 10000, seed: int = 0,
) -> tuple[float, float, float, int, int] | None:
    """Kendall τ-b with a cluster bootstrap over networks.

    Args:
        rows: Analysis-unit rows.
        predictor: Predictor column.
        endpoint: Endpoint column.
        n_boot: Resamples.
        seed: Seed for ``np.random.default_rng``.

    Returns:
        ``(tau, lo, hi, n, n_networks)``, or ``None`` if τ is undefined.
    """
    sel = [
        r for r in rows
        if r[predictor] not in ("", "None") and r[endpoint] not in ("", "None")
        and not (predictor == "radius" and float(r[predictor]) == UNREACHED)
    ]
    if len(sel) < 2:
        return None
    x = np.array([float(r[predictor]) for r in sel])
    y = np.array([float(r[endpoint]) for r in sel])
    if np.unique(x).size <= 1 or np.unique(y).size <= 1:
        return None
    tau, _ = kendalltau(x, y, variant="b")
    by_net: dict[str, list[int]] = defaultdict(list)
    for i, r in enumerate(sel):
        by_net[r["network"]].append(i)
    nets = sorted(by_net)
    rng = np.random.default_rng(seed)
    taus: list[float] = []
    for _ in range(n_boot):
        idx = [i for p in rng.integers(0, len(nets), size=len(nets)) for i in by_net[nets[p]]]
        t, _ = kendalltau(x[idx], y[idx], variant="b")
        if not np.isnan(t):
            taus.append(t)
    if len(taus) < 100:
        return None
    lo, hi = np.percentile(taus, [2.5, 97.5])
    return float(tau), float(lo), float(hi), len(sel), len(nets)


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
            full = cluster_tau(us, predictor, endpoint, n_boot=args.n_boot)
            if full is None:
                continue
            f_tau, f_lo, f_hi, _n, _nn = full
            f_excl = bool(f_lo > 0 or f_hi < 0)
            for net in sorted({u["network"] for u in us}):
                sub = [u for u in us if u["network"] != net]
                res = cluster_tau(sub, predictor, endpoint, n_boot=args.n_boot)
                if res is None:
                    continue
                tau, lo, hi, n, nn = res
                excl = bool(lo > 0 or hi < 0)
                rows.append({
                    "stratum": stratum, "stratum_class": us[0]["stratum_class"],
                    "predictor": predictor, "endpoint": endpoint,
                    "dropped_network": net, "n": n, "n_networks": nn,
                    "tau_b": tau, "ci_lo_2p5": lo, "ci_hi_97p5": hi,
                    "excludes_zero": excl, "full_tau_b": f_tau,
                    "full_excludes_zero": f_excl,
                    "loses_zero_exclusion": bool(f_excl and not excl),
                })

    path = in_dir / "analysis_loo_bootstrap.csv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    fragile: dict[str, list[str]] = defaultdict(list)
    for r in rows:
        if r["loses_zero_exclusion"] and r["predictor"] == "radius":
            fragile[r["stratum"]].append(r["dropped_network"])
    summary = {
        "note": "POST-HOC. The pre-registered point-estimate LOO is reported unchanged in analysis_tau.csv.",
        "n_rows": len(rows),
        "radius_strata_that_lose_zero_exclusion_when_one_network_is_dropped": dict(fragile),
        "radius_strata_whose_interval_excludes_zero_with_all_networks": sorted(
            {r["stratum"] for r in rows if r["predictor"] == "radius" and r["full_excludes_zero"]}
        ),
        "radius_strata_robust_to_dropping_any_single_network": sorted(
            {r["stratum"] for r in rows if r["predictor"] == "radius" and r["full_excludes_zero"]}
            - set(fragile)
        ),
    }
    (in_dir / "analysis_loo_bootstrap.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
