"""E17: worked applied case, and the share invalid at truth by radius bin.

From release/per_query_records.csv (status ok rows): (1) the Polzer_2012 and
Acid_1996 states of the primary condition (D_LLM), every query with its radius,
witness and validity at truth; (2) the share of committed sets invalid in the
data-generating DAG by radius bin (1, 2, >=3, inf) for the eight real-name
elicited conditions pooled, for D_LLM alone and for the matched random states
(D_RAND, three replicates), with network-cluster bootstrap intervals
(10,000 resamples). Also the AUC of 1/r for invalidity at truth.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from common import OUT, RELEASE

REAL = ["D_LLM", "D_LLM_INSTR", "D_LLM_SRC", "D_LLM_72B", "D_LLM_72B_INSTR", "D_LLM_32B",
        "D_LLM_SRC_70B", "D_LLM_GEMMA_27B"]


def rbin(r):
    if r == -1:
        return "inf"
    return "1" if r == 1 else ("2" if r == 2 else ">=3")


def auc(score, y):
    pos, neg = score[y == 1], score[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    s = 0.0
    for p in pos:
        s += np.sum(p > neg) + 0.5 * np.sum(p == neg)
    return s / (len(pos) * len(neg))


def table(df, nb=10_000, seed=0):
    df = df.copy()
    df["bin"] = df.r_hop.map(rbin)
    df["inv"] = 1 - df.z_valid_at_truth
    nets = sorted(df.network.unique())
    rng = np.random.default_rng(seed)
    g = {n: d for n, d in df.groupby("network")}
    out = {}
    for b in ["1", "2", ">=3", "inf"]:
        sub = df[df.bin == b]
        out[b] = {"n": int(len(sub)), "share_invalid": float(sub.inv.mean()) if len(sub) else None}
    boots = {b: [] for b in out}
    aucs = []
    for _ in range(nb):
        pick = [nets[i] for i in rng.integers(0, len(nets), len(nets))]
        bd = pd.concat([g[n] for n in pick])
        for b in out:
            s = bd[bd.bin == b]
            if len(s):
                boots[b].append(s.inv.mean())
        sc = np.where(bd.r_hop == -1, 0.0, 1.0 / bd.r_hop.clip(lower=1))
        aucs.append(auc(sc, bd.inv.values))
    for b in out:
        if boots[b]:
            out[b]["lo"], out[b]["hi"] = map(float, np.percentile(boots[b], [2.5, 97.5]))
    sc = np.where(df.r_hop == -1, 0.0, 1.0 / df.r_hop.clip(lower=1))
    out["auc_1_over_r_for_invalid"] = float(auc(sc, df.inv.values))
    out["auc_ci"] = list(map(float, np.nanpercentile(aucs, [2.5, 97.5])))
    out["n_networks"] = len(nets)
    return out


if __name__ == "__main__":
    d = pd.read_csv(RELEASE / "per_query_records.csv")
    ok = d[(d.status == "ok") & d.r_hop.notna() & d.z_valid_at_truth.notna()]
    res = {"real_conditions": table(ok[ok.arm.isin(REAL)]),
           "D_LLM": table(ok[ok.arm == "D_LLM"]),
           "D_RAND": table(ok[ok.arm == "D_RAND"]),
           "A_TRUE": table(ok[ok.arm == "A_TRUE"])}
    cols = ["x", "y", "z_size", "r_hop", "r_claim", "witness", "z_valid_at_truth"]
    for net in ("Polzer_2012", "Acid_1996"):
        s = ok[(ok.arm == "D_LLM") & (ok.network == net)]
        res[f"worked_{net}"] = s[cols].to_dict(orient="records")
        res[f"worked_{net}_state"] = s[["n_k", "k_g0", "n_wrong_claims"]].iloc[0].to_dict()
    p = OUT / "e17" / "worked_case.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(res, indent=1, default=float))
    for k in ("real_conditions", "D_LLM", "D_RAND", "A_TRUE"):
        print(k, json.dumps(res[k]))
