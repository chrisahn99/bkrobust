"""E11 analysis: the Table 1 statistics on the CPDAG-only frame.

Units are built by the submitted analysis (bkrobust.robustness.llm_analyse,
run on results/rebuttal/e11). Statistics are the paper's: pooled tau_b with a
network-cluster bootstrap (real_analyse.tau_for_stratum) and the within-state
stratified tau_b (pooled within-state pair counts, states where r_val varies,
network-cluster bootstrap), 10,000 resamples, seed 0. Reported for the whole
frame and for its bridge stratum (pairs that were also in the paper's Table 1
admissible set) versus the rest.
"""

from __future__ import annotations

import csv
import json

from common import OUT
from w4_corrected_tables import paired_within, states, strat, tfs


def load():
    rows = list(csv.DictReader(open(OUT / "e11" / "analysis_units.csv")))
    frame = {(r["network"], r["x"], r["y"]): r for r in
             map(json.loads, open(OUT / "e11_frame" / "frame.jsonl"))}
    num = ("radius", "k_g0", "n_k", "shd_truth", "AUC_frac", "k_accuracy")
    for u in rows:
        for c in num:
            v = u.get(c)
            try:
                u[c] = float(v) if v not in ("", None) else None
            except ValueError:
                u[c] = None
        u["separation"] = None
        f = frame.get((u["network"], u["x"], u["y"]), {})
        u["bridge"] = bool(f.get("in_committed_admissible"))
    return rows


def block(units, label):
    pooled = [u for u in units if u.get("naming") == "real"]
    out = {"label": label, "n_rows": len(pooled),
           "n_scored": sum(1 for u in pooled if u.get("AUC_frac") is not None and u.get("radius") is not None)}
    out["pooled"] = {p: tfs(pooled, p) for p in ("radius", "k_g0", "n_k", "shd_truth", "k_accuracy")}
    st = states(units)
    out["within"] = strat(st, "radius") if st else None
    if st:
        loo = {}
        for net in sorted({k[0] for k in st}):
            sub = {k: v for k, v in st.items() if k[0] != net}
            if sub:
                loo[net] = round(strat_point(sub), 3)
        out["within_loo_min"] = min(loo.values())
        out["within_loo_max"] = max(loo.values())
        out["within_loo_argmin"] = min(loo, key=loo.get)
    return out


def strat_point(st):
    import llm_paired_v2 as base
    return base._stratified_tau_b(st, "radius", "AUC_frac")["tau_b"]


if __name__ == "__main__":
    import sys
    from common import ROOT
    sys.path.insert(0, str(ROOT / "experiments"))
    u = load()
    res = {"all": block(u, "CPDAG-only frame, 528 queries"),
           "bridge": block([x for x in u if x["bridge"]], "bridge: also in Table 1 frame"),
           "non_bridge": block([x for x in u if not x["bridge"]], "not in Table 1 frame")}
    st_all = states(u)
    res["status_counts"] = {}
    for x in u:
        if x.get("naming") == "real":
            res["status_counts"][x["status"]] = res["status_counts"].get(x["status"], 0) + 1
    (OUT / "e11" / "e11_tables.json").write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps(res, indent=1, default=str)[:6000])
