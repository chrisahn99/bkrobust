"""E5/E6/E7 analysis.

E5: by knowledge arm and sample size (oracle CPDAG as reference): share of
frame queries where the pipeline commits a set (status ok), where it wrongly
concludes no causal path (y_not_descendant, although every frame pair has a
true causal path), and, among committed sets, the share invalid in the true
DAG. Also invalidity by radius bin and by skeleton errors next to the query.
E6: AUC of (1 - bootstrap stability) for invalidity at truth, and the share
invalid among fully stable and among unstable sets.
E7: invalidity by where the hidden variable sits.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from common import OUT


def auc(score, y):
    pos, neg = score[y == 1], score[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return None
    return float(np.mean([(p > neg).mean() + 0.5 * (p == neg).mean() for p in pos]))


def rb(r):
    return "inf" if str(r) == "inf" else ("1" if str(r) in ("1", "1.0") else
                                          ("2" if str(r) in ("2", "2.0") else ">=3"))


def main():
    out = {}
    p = OUT / "e5" / "learned_rows_boot.csv"
    if not p.exists():
        p = OUT / "e5" / "learned_rows.csv"
    d = pd.read_csv(p, low_memory=False)
    d["inv"] = 1 - d.valid_truth.map({True: 1, False: 0, "True": 1, "False": 0})
    summ = []
    for (arm, n), g in d.groupby(["arm", "n"]):
        ok = g[g.status == "ok"]
        summ.append({"arm": arm, "n": n, "queries": len(g),
                     "share_committed": round((g.status == "ok").mean(), 3),
                     "share_no_path_claimed": round((g.status == "y_not_descendant").mean(), 3),
                     "share_not_amenable": round((g.status == "not_amenable").mean(), 3),
                     "n_committed": len(ok),
                     "share_invalid_at_truth": round(ok.inv.mean(), 3) if len(ok) else None,
                     "by_network_invalid": {k: round(v, 3) for k, v in ok.groupby("network").inv.mean().items()}})
    out["e5_by_arm_n"] = summ
    lr = d[(d.status == "ok") & (d.n != "oracle")].copy()
    lr["bin"] = lr.r.map(rb)
    out["e5_invalid_by_radius_bin"] = {
        arm: {b: {"n": int(len(s)), "share_invalid": round(s.inv.mean(), 3)}
              for b, s in g.groupby("bin")} for arm, g in lr.groupby("arm")}
    lr["near"] = (pd.to_numeric(lr.near_skel_err) > 0).map({True: "near_error", False: "no_near_error"})
    out["e5_invalid_by_near_skeleton_error"] = {
        arm: {k: {"n": int(len(s)), "share_invalid": round(s.inv.mean(), 3)}
              for k, s in g.groupby("near")} for arm, g in lr.groupby("arm")}
    sc = np.where(lr.r.astype(str) == "inf", 0.0, 1.0 / pd.to_numeric(lr.r, errors="coerce").fillna(1e9).clip(lower=1))
    out["e5_auc_1_over_r"] = auc(sc, lr.inv.values)
    out["e5_auc_near_error"] = auc(pd.to_numeric(lr.near_skel_err).values.astype(float), lr.inv.values)
    if "boot_stability" in d:
        b = d[d.boot_stability.notna() & (d.status == "ok")].copy()
        out["e6"] = {"n": int(len(b)), "n_invalid": int(b.inv.sum()),
                     "auc_instability": auc(1 - b.boot_stability.values, b.inv.values),
                     "share_invalid_if_stable_1": round(b[b.boot_stability >= 0.999].inv.mean(), 3),
                     "n_stable_1": int((b.boot_stability >= 0.999).sum()),
                     "share_invalid_if_unstable": round(b[b.boot_stability < 0.999].inv.mean(), 3),
                     "n_unstable": int((b.boot_stability < 0.999).sum()),
                     "share_invalid_if_stable_ge_0.9": round(b[b.boot_stability >= 0.9].inv.mean(), 3),
                     "auc_combined_instability_plus_near": auc(
                         (1 - b.boot_stability.values) + (pd.to_numeric(b.near_skel_err).values > 0),
                         b.inv.values),
                     "by_arm": {arm: {"n": int(len(g)),
                                      "auc": auc(1 - g.boot_stability.values, g.inv.values),
                                      "inv_if_stable": round(g[g.boot_stability >= 0.999].inv.mean(), 3)
                                      if (g.boot_stability >= 0.999).any() else None}
                                for arm, g in b.groupby("arm")}}
    h = OUT / "e5" / "hidden_rows.csv"
    if h.exists():
        e = pd.read_csv(h)
        e["inv"] = 1 - e.valid_truth.map({True: 1, False: 0, "True": 1, "False": 0})
        out["e7"] = {w: {"n": int(len(g)), "committed": int((g.status == "ok").sum()),
                         "no_path_claimed": int((g.status == "y_not_descendant").sum()),
                         "share_invalid": round(g[g.status == "ok"].inv.mean(), 3),
                         "share_invalid_hidden_confounder": round(
                             g[(g.status == "ok") & (g.hidden_is_confounder)].inv.mean(), 3)
                         if ((g.status == "ok") & (g.hidden_is_confounder)).any() else None}
                     for w, g in e.groupby("where")}
    (OUT / "e5" / "analysis.json").write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
