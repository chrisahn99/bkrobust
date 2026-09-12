"""Round-2 review, W1: the separation baseline at N = 1000 in the repaired cell.

`review1_unstable_stratum_n1000.py` re-measured the radius in
`flip, coverage=1.0, base_wrongness=0.00` but not the separation baseline, so the
paper's table mixed an N = 1000 radius with an N = 200 separation in one row. This
scores both predictors on the same N = 1000 endpoints, using the same bootstrap
(10,000 resamples, seed 0) as the published tau tables.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from bkrobust.core.conventions import UNREACHED  # noqa: E402

CELL = ROOT / "results" / "axis_robustness" / "review1_unstable_cell"
OUT = ROOT / "results" / "axis_robustness" / "review1_unstable_cell_separation.csv"


def boot_ci(x, y, n_boot=10_000, seed=0):
    rng = np.random.default_rng(seed)
    n = len(x)
    taus = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        taus[i] = kendalltau(x[idx], y[idx], variant="b", nan_policy="propagate")[0]
    ok = taus[~np.isnan(taus)]
    return float(np.percentile(ok, 2.5)), float(np.percentile(ok, 97.5))


def main() -> None:
    inst = pd.read_csv(CELL / "null_instances.csv")
    rows = []
    for endpoint in ("AUC_frac_n1000", "AUC_frac_usable_n1000", "AUC_frac_n200"):
        for pred in ("r_val", "separation"):
            sub = inst[[pred, endpoint]].copy()
            if pred == "r_val":
                sub = sub[sub["r_val"] != UNREACHED]
            sub = sub.dropna()
            x = sub[pred].to_numpy(float)
            y = sub[endpoint].to_numpy(float)
            t, p = kendalltau(x, y, variant="b")
            lo, hi = boot_ci(x, y)
            rows.append({"stratum": "flip, coverage=1.0, base_wrongness=0.00",
                         "endpoint": endpoint, "predictor": pred, "n": len(sub),
                         "tau_b": t, "p_value": p, "ci_lo_2p5": lo, "ci_hi_97p5": hi})
            print(f"{endpoint:24s} {pred:11s} n={len(sub):4d} tau={t:+.3f} [{lo:+.3f},{hi:+.3f}]")
    ident = int((inst["r_val"] == inst["separation"]).sum())
    print(f"\nr_val == separation on {ident}/{len(inst)} instances of this stratum")
    pd.DataFrame(rows).to_csv(OUT, index=False)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
