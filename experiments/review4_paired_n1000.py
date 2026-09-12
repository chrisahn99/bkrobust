"""All six flip strata at N = 1000, radius against separation, paired.

The decision table mixed draw counts: two cells re-measured at N = 1000 and four
still at N = 200 on an endpoint the project has shown is draw-count dependent. All
six now have an N = 1000 run, so this scores both predictors on all of them, on the
same instances, with the same paired bootstrap.
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

R = ROOT / "results" / "axis_robustness"
CELLS = {
    (0.5, 0.00): R / "review4_strata_n1000" / "cov0.5_bw0.0",
    (0.5, 0.10): R / "review4_strata_n1000" / "cov0.5_bw0.1",
    (0.5, 0.25): R,                                     # the committed null-cell run
    (1.0, 0.00): R / "review1_unstable_cell",           # earlier single-cell run
    (1.0, 0.10): R / "review4_strata_n1000" / "cov1.0_bw0.1",
    (1.0, 0.25): R / "review4_strata_n1000" / "cov1.0_bw0.25",
}
OUT = R / "review4_paired_n1000.csv"


def paired(xr, xs, y, seed=0, n_boot=10_000):
    tr = kendalltau(xr, y, variant="b")[0]
    ts = kendalltau(xs, y, variant="b")[0]
    rng = np.random.default_rng(seed)
    n = len(y)
    d = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        d[i] = (kendalltau(xr[idx], y[idx], variant="b")[0]
                - kendalltau(xs[idx], y[idx], variant="b")[0])
    ok = d[~np.isnan(d)]
    return tr, ts, tr - ts, float(np.percentile(ok, 2.5)), float(np.percentile(ok, 97.5))


def main() -> None:
    rows = []
    print(f"{'stratum':22s} {'n':>4s} {'tau(r)':>8s} {'tau(sep)':>9s} {'diff':>8s}  verdict")
    for (cov, bw), d in CELLS.items():
        f = d / "null_instances.csv"
        if not f.exists():
            print(f"  cov {cov} bw {bw}: MISSING {f}")
            continue
        inst = pd.read_csv(f)
        # the null-cell directory holds one stratum; make sure it is this one
        if "coverage" in inst and not np.isclose(inst.coverage.iloc[0], cov):
            print(f"  cov {cov} bw {bw}: directory holds coverage "
                  f"{inst.coverage.iloc[0]}, skipping")
            continue
        ep = "AUC_frac_usable_n1000"
        s = inst[["r_val", "separation", ep]].copy()
        s = s[s.r_val != UNREACHED].dropna()
        tr, ts, diff, lo, hi = paired(s.r_val.to_numpy(float),
                                      s.separation.to_numpy(float), s[ep].to_numpy(float))
        verdict = "radius" if lo > 0 else "separation" if hi < 0 else "indistinguishable"
        rows.append({"coverage": cov, "base_wrongness": bw, "n": len(s), "endpoint": ep,
                     "tau_r_hop": tr, "tau_separation": ts, "diff": diff,
                     "ci_lo": lo, "ci_hi": hi, "verdict": verdict})
        print(f"cov {cov} bw {bw:<6.2f} {len(s):>4d} {tr:+8.3f} {ts:+9.3f} "
              f"{diff:+8.3f} [{lo:+.3f},{hi:+.3f}]  {verdict}")
    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"\nwrote {OUT}")
    print("verdicts:", out.verdict.value_counts().to_dict())


if __name__ == "__main__":
    main()
