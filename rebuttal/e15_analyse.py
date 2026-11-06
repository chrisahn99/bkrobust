"""E15 analysis: Kendall tau_b of each predictor with the survival AUC, per
stratum (as in the paper's Table 10/12), instance bootstrap (2,000 resamples),
and paired differences r_val minus AID variants on the same resamples.

Sign convention: AID counts mistakes, so larger means worse; it is reported
negated (-AID), like the negated SHD of the paper, so that a positive tau_b
means "ranks survival in the expected direction". r_val = inf (-1) is
excluded as in the paper's analysis.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

from common import OUT

NB = 2000


def stratum(r):
    if r.which == "untargeted":
        return f"untargeted, {int(r.n_tiers)} tiers"
    if r.arm == "flip":
        return f"flip, cov {r.coverage}, bw {r.base_wrongness}"
    return f"tiered, {int(r.n_tiers)} tiers"


def tau(x, y):
    if np.unique(x).size < 2 or np.unique(y).size < 2:
        return np.nan
    return kendalltau(x, y, variant="b").statistic


PRED = {"r_val": "r_val", "-aid_oracle": "aid_oracle", "-aid_analyst": "aid_analyst",
        "-aidx_oracle": "aidx_oracle", "-aidx_analyst": "aidx_analyst",
        "-shd_truth": "shd_truth", "separation": "separation"}


def main():
    d = pd.read_csv(OUT / "e15" / "aid_instances.csv")
    rep = d.reproduced.value_counts().to_dict()
    d = d[d.reproduced == "ok"].copy()
    d["r_val"] = d.r_val.astype(float)
    d = d[d.r_val >= 0]
    for k in ("aid_oracle", "aid_analyst", "aidx_oracle", "aidx_analyst", "shd_truth"):
        d[k] = -d[k].astype(float)
    d["separation"] = pd.to_numeric(d.separation, errors="coerce")
    d["stratum"] = d.apply(stratum, axis=1)
    res = {"reproduction": rep, "strata": {}}
    rng = np.random.default_rng(0)
    for s, g in sorted(d.groupby("stratum")):
        y = g.AUC_frac.values.astype(float)
        row = {"n": int(len(g))}
        boots = {p: [] for p in PRED}
        idx_all = [rng.integers(0, len(g), len(g)) for _ in range(NB)]
        for p, col in PRED.items():
            x = g[col].values.astype(float)
            m = ~np.isnan(x)
            row[p] = {"tau": tau(x[m], y[m]), "n": int(m.sum())}
            for idx in idx_all:
                xi, yi = x[idx], y[idx]
                mm = ~np.isnan(xi)
                boots[p].append(tau(xi[mm], yi[mm]))
            b = np.array(boots[p])
            if np.isfinite(b).sum() > 100:
                row[p]["lo"], row[p]["hi"] = map(float, np.nanpercentile(b, [2.5, 97.5]))
        for p in ("-aid_oracle", "-aid_analyst", "-aidx_oracle", "-aidx_analyst", "-shd_truth"):
            diff = np.array(boots["r_val"]) - np.array(boots[p])
            if np.isfinite(diff).sum() > 100:
                row[f"r_val minus {p}"] = {
                    "delta": row["r_val"]["tau"] - row[p]["tau"],
                    "lo": float(np.nanpercentile(diff, 2.5)),
                    "hi": float(np.nanpercentile(diff, 97.5))}
        # raw sign: the comparison against +AID (= -(-AID))
        for p in ("-aid_oracle", "-aid_analyst", "-aidx_oracle", "-aidx_analyst", "-shd_truth"):
            diff = np.array(boots["r_val"]) + np.array(boots[p])
            if np.isfinite(diff).sum() > 100:
                row[f"r_val minus raw {p[1:]}"] = {
                    "delta": row["r_val"]["tau"] + row[p]["tau"],
                    "lo": float(np.nanpercentile(diff, 2.5)),
                    "hi": float(np.nanpercentile(diff, 97.5))}
        res["strata"][s] = row
    (OUT / "e15" / "aid_tau.json").write_text(json.dumps(res, indent=1, default=float))
    print("reproduction:", rep)
    for s, row in res["strata"].items():
        def f(p):
            v = row[p]
            return (f"{v['tau']:+.2f}[{v.get('lo', float('nan')):+.2f},{v.get('hi', float('nan')):+.2f}]"
                    if v["tau"] == v["tau"] else "  n/a")
        print(f"{s:28s} n={row['n']:4d} r={f('r_val')} -AIDo={f('-aid_oracle')} "
              f"-AIDa={f('-aid_analyst')} -AIDxa={f('-aidx_analyst')} -SHD={f('-shd_truth')} "
              f"s={f('separation')}")
        for k in [k for k in row if k.startswith("r_val minus raw")]:
            if k in row:
                print(f"{'':30s}{k}: {row[k]['delta']:+.2f} [{row[k]['lo']:+.2f},{row[k]['hi']:+.2f}]")


if __name__ == "__main__":
    main()
