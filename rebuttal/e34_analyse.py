"""E34 analysis: Kendall tau_b between the radius on C_hat (inf highest) and the survival
AUC, pooled per (arm, learner, n), graph-cluster bootstrap (2,000 draws, seed 20261111);
the true CPDAG is the reference. Writes results/rebuttal/e34/analysis.json."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

from common import OUT

INF = 10 ** 6


def tau(df):
    if len(df) < 3 or df.rnum.nunique() < 2 or df.auc.nunique() < 2:
        return np.nan
    return kendalltau(df.rnum, df.auc, variant="b").statistic


def main():
    d = pd.DataFrame([json.loads(x) for x in (OUT / "e34" / "units.jsonl").open()])
    d = d[(d.status == "ok") & d.auc.notna()].copy()
    d["rnum"] = d.r.map(lambda v: INF if str(v) == "inf" else float(v))
    d["n"] = d.n.astype(str)
    rng = np.random.default_rng(20261111)
    out = []
    for (arm, lr, n), c in d.groupby(["arm", "learner", "n"]):
        t = tau(c)
        groups = [g.index.values for _, g in c.groupby("graph_id")]
        bs = []
        for _ in range(2000):
            ix = np.concatenate([groups[i] for i in rng.integers(0, len(groups), len(groups))])
            bs.append(tau(c.loc[ix]))
        bs = np.array([b for b in bs if not np.isnan(b)])
        out.append({"arm": arm, "learner": lr, "n": n, "units": len(c), "graphs": len(groups),
                    "tau_b": None if np.isnan(t) else round(float(t), 3),
                    "ci": [round(float(np.quantile(bs, q)), 3) for q in (0.025, 0.975)] if len(bs) else None,
                    "share_r_inf": round(float((c.rnum >= INF).mean()), 3),
                    "mean_auc": round(float(c.auc.mean()), 3)})
    t = pd.DataFrame(out)
    t.to_csv(OUT / "e34" / "e34_table.csv", index=False)
    (OUT / "e34" / "analysis.json").write_text(json.dumps(out, indent=1))
    print(t.to_string())


if __name__ == "__main__":
    main()
