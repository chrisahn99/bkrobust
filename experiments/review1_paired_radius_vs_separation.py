"""Round-3 review, W1: test the radius-vs-separation difference instead of eyeballing it.

The decision table compared two point estimates of tau_b per stratum and declared a
winner. Both predictors are measured on the same instances, so the quantity that decides
the comparison is the paired difference tau(r_hop) - tau(separation), resampled over
instances. This computes it on both endpoints the project uses, so the endpoint
sensitivity is visible rather than hidden by showing one.
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
from bkrobust.robustness import analyse as an  # noqa: E402
from bkrobust.robustness import build_tau_comparisons as btc  # noqa: E402

R = ROOT / "results" / "axis_robustness"
OUT = R / "review1_paired_radius_vs_separation.csv"
N_BOOT = 10_000


def paired(x_r, x_s, y, seed=0, n_boot=N_BOOT):
    t_r = kendalltau(x_r, y, variant="b")[0]
    t_s = kendalltau(x_s, y, variant="b")[0]
    rng = np.random.default_rng(seed)
    n = len(y)
    d = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        d[i] = (kendalltau(x_r[idx], y[idx], variant="b")[0]
                - kendalltau(x_s[idx], y[idx], variant="b")[0])
    ok = d[~np.isnan(d)]
    return t_r, t_s, t_r - t_s, float(np.percentile(ok, 2.5)), float(np.percentile(ok, 97.5))


def main() -> None:
    inst = pd.read_csv(R / "survival_instances.csv")
    curves = pd.read_csv(R / "survival_curves.csv")
    inst = an.attach_auc_frac_usable(inst, curves)
    reb = pd.read_csv(R / "analysis_tiered_rebinned.csv")
    rate = [c for c in reb.columns if c.startswith("AUC_frac_rate")]
    if rate:
        inst = inst.merge(reb[["instance_id"] + rate[:1]], on="instance_id", how="left")
    # The anchor table scores the tiered arm on AUC_rate_usable, so the paired test must
    # use that endpoint too or it is not comparing what the paper reports.
    anchor_inputs = btc.load_inputs(str(R))
    rate_usable = btc.compute_auc_rate_usable(anchor_inputs["tiered_rebinned"]).table
    inst = inst.merge(rate_usable[["instance_id", "AUC_rate_usable"]],
                      on="instance_id", how="left")

    rows = []
    for arm, keys in (("flip", ["coverage", "base_wrongness"]), ("tiered", ["n_tiers"])):
        sub = inst[inst.arm == arm]
        endpoints = (["AUC_frac", "AUC_frac_usable"] if arm == "flip"
                     else [rate[0] if rate else "AUC_frac_rate", "AUC_rate_usable"])
        for key, g in sub.groupby(keys):
            for ep in endpoints:
                if ep not in g:
                    continue
                d = g[["r_val", "separation", ep]].copy()
                d = d[d.r_val != UNREACHED].dropna()
                if len(d) < 30:
                    continue
                t_r, t_s, diff, lo, hi = paired(
                    d.r_val.to_numpy(float), d.separation.to_numpy(float), d[ep].to_numpy(float))
                verdict = ("radius" if lo > 0 else "separation" if hi < 0 else "indistinguishable")
                rows.append({"arm": arm, "stratum": str(key), "endpoint": ep, "n": len(d),
                             "tau_r_hop": t_r, "tau_separation": t_s, "diff": diff,
                             "ci_lo": lo, "ci_hi": hi, "verdict": verdict})
                print(f"{arm:7s} {str(key):18s} {ep:18s} n={len(d):4d} "
                      f"r={t_r:+.3f} sep={t_s:+.3f} d={diff:+.3f} [{lo:+.3f},{hi:+.3f}] {verdict}")
    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"\nwrote {OUT}")
    for ep, g in out.groupby("endpoint"):
        c = g.verdict.value_counts().to_dict()
        print(f"  {ep}: {c}")


if __name__ == "__main__":
    main()
