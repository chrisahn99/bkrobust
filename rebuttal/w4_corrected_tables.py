"""C2Bd-W4 and E16: every published-network table on the corrected cohort.

Corrected cohort = the original analysis units with the 94 real-name units
that the pipeline wrongly excluded (valid empty optimal set) re-scored. The
earlier corrected-cohort script appended the re-scored rows but kept the
superseded (excluded, non-scored) rows of the same units; here the superseded
rows are dropped, so each (network, x, y, condition) appears once.

Outputs (results/rebuttal/w4/):
  table1.json      pooled / balanced 8-of-8 / within-state, all predictors
  table11.json     balanced levels >= k/8
  table13.json     within-state tau_b per condition (+ scrambled controls)
  e16.json         radius where separation is undefined (pooled, within state),
                   and the within-state paired lead radius - separation
Same statistics, bootstrap (network clusters, 10,000 resamples, seed 0) and
code paths as the paper (experiments/llm_paired_v2.py, real_analyse.py).
"""

from __future__ import annotations

import sys
from collections import Counter, defaultdict

import numpy as np

from common import OUT, ROOT, dump

sys.path.insert(0, str(ROOT / "experiments"))
import llm_paired_v2 as base  # noqa: E402
import table1_corrected_cohort as t1  # noqa: E402

from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.robustness.real_analyse import tau_for_stratum  # noqa: E402

E = "AUC_frac"
NB, SEED = 10_000, 0


def key(u):
    return (u["network"], u["x"], u["y"], u["condition"])


def cohort():
    old, new, _ = t1.build_cohorts()
    newkeys = {key(u) for u in new}
    superseded = [u for u in old if key(u) in newkeys]
    kept = [u for u in old if key(u) not in newkeys]
    comb = kept + new
    c = Counter(key(u) for u in comb)
    assert max(c.values()) == 1
    return old, comb, len(superseded)


def states(units, cond_filter=None, require_var=True):
    by = defaultdict(list)
    for u in base.real_naming(units) if cond_filter is None else units:
        if cond_filter is not None and u["condition"] != cond_filter:
            continue
        rv, ev = u.get("radius"), u.get(E)
        if rv is None or ev is None or rv == UNREACHED:
            continue
        by[(u["network"], u["condition"])].append(u)
    if require_var:
        by = {k: v for k, v in by.items() if len({r["radius"] for r in v}) > 1}
    return by


def strat(st, x):
    r = base._stratified_tau_b(st, x, E)
    lo, hi, nb = base._within_state_bootstrap(st, x, E, n_boot=NB, seed=SEED)
    return {"tau_b": r["tau_b"], "lo": lo, "hi": hi, "n_units": sum(map(len, st.values())),
            "n_states": len(st), "n_networks": len({k[0] for k in st}), "C": r["C"], "D": r["D"]}


def paired_within(st, a, b):
    """Network-cluster bootstrap of the within-state tau_b difference a - b on
    the same resampled states (same resampling as _within_state_bootstrap)."""
    keys = list(st)
    nets = sorted({k[0] for k in keys})
    by_net = defaultdict(list)
    for k in keys:
        by_net[k[0]].append(k)
    rng = np.random.default_rng(SEED)
    arr = np.array(nets, dtype=object)
    d = []
    for _ in range(NB):
        res = {}
        for i, n in enumerate(arr[rng.integers(0, len(nets), size=len(nets))]):
            for k in by_net[n]:
                res[(i, k)] = st[k]
        ta = base._stratified_tau_b(res, a, E)["tau_b"]
        tb = base._stratified_tau_b(res, b, E)["tau_b"]
        if ta is not None and tb is not None:
            d.append(ta - tb)
    obs = base._stratified_tau_b(st, a, E)["tau_b"] - base._stratified_tau_b(st, b, E)["tau_b"]
    lo, hi = np.percentile(d, [2.5, 97.5])
    return {"delta": obs, "lo": float(lo), "hi": float(hi), "n_valid_boot": len(d),
            "p_delta_gt_0": float(np.mean(np.array(d) > 0)),
            "n_networks": len(nets), "n_states": len(keys)}


def tfs(rows, p):
    r = tau_for_stratum(rows, p, E, n_boot=NB, seed=SEED)
    return {k: r[k] for k in ("tau_b", "ci_lo_2p5", "ci_hi_97p5", "n", "n_networks", "status",
                              "loo_tau_min", "loo_tau_max")}


def run(units, label):
    pooled = base.real_naming(units)
    okt = base.ok_conditions_by_triple(units)
    out = {"cohort": label}
    scored = [u for u in pooled if u.get(E) is not None and u.get("radius") is not None]
    out["pooled_scored_units"] = len(scored)
    bal = base.subset_at_least(units, okt, 8)
    bal_s = [u for u in bal if u.get(E) is not None]
    out["balanced_scored_units"] = len(bal_s)
    out["balanced_queries"] = len({(u["network"], u["x"], u["y"]) for u in bal_s})
    t = {"pooled": {}, "balanced": {}}
    for p in t1.PREDICTORS:
        t["pooled"][p] = tfs(pooled, p)
        t["balanced"][p] = tfs(bal, p)
    st = states(units)
    t["within"] = {"radius": strat(st, "radius")}
    sep_st = {k: [r for r in v if r.get("separation") is not None] for k, v in st.items()}
    sep_st = {k: v for k, v in sep_st.items() if len(v) >= 2}
    t["within"]["separation_where_defined"] = strat(sep_st, "separation")
    t["within"]["radius_where_sep_defined"] = strat(sep_st, "radius")
    t["within"]["paired_radius_minus_sep"] = paired_within(sep_st, "radius", "separation")
    # E16: units where separation is undefined
    und_pooled = [u for u in pooled if u.get("separation") is None]
    t["e16_pooled_radius_sep_undefined"] = tfs(und_pooled, "radius")
    und_st = {k: [r for r in v if r.get("separation") is None] for k, v in st.items()}
    und_st = {k: v for k, v in und_st.items() if len(v) >= 2 and len({r["radius"] for r in v}) > 1}
    t["e16_within_radius_sep_undefined"] = strat(und_st, "radius")
    t["share_within_units_sep_undefined"] = 1 - (
        sum(map(len, sep_st.values())) / max(1, sum(map(len, st.values()))))
    out["table1"] = t
    out["table11"] = []
    for k in (8, 7, 6):
        sub = base.subset_at_least(units, okt, k)
        out["table11"].append({"k": k, "n_scored": len([u for u in sub if u.get(E) is not None]),
                               "n_networks": len({u["network"] for u in sub if u.get(E) is not None}),
                               **{p: tfs(sub, p) for p in ("radius", "k_g0", "n_k", "shd_truth")}})
    out["table13"] = {}
    conds = sorted({u["condition"] for u in units})
    for c in conds:
        stc = states(units, cond_filter=c)
        if stc:
            out["table13"][c] = strat(stc, "radius")
    out["paired_pooled"] = {}
    from bkrobust.robustness.paired_resample import paired_comparison
    for other in ("shd_truth", "k_g0", "n_k"):
        r = paired_comparison(pooled, "radius", other, E, "network", n_boot=NB, seed=SEED)
        out["paired_pooled"][other] = {k: r.get(k) for k in ("delta", "ci_lo", "ci_hi", "ci_lo_2p5",
                                                              "ci_hi_97p5", "p_delta_gt_0")}
    return out


if __name__ == "__main__":
    old, comb, n_sup = cohort()
    print("superseded rows dropped:", n_sup)
    res = {"n_superseded_dropped": n_sup, "original": run(old, "original"),
           "corrected": run(comb, "corrected")}
    dump(res, OUT / "w4" / "tables.json")
    for lab in ("original", "corrected"):
        r = res[lab]
        t = r["table1"]
        print(f"== {lab}: pooled scored {r['pooled_scored_units']}, balanced "
              f"{r['balanced_scored_units']} ({r['balanced_queries']} queries)")
        for p in t1.PREDICTORS:
            a, b = t["pooled"][p], t["balanced"][p]
            print(f"  {p:12s} pooled {a['tau_b']} [{a['ci_lo_2p5']},{a['ci_hi_97p5']}] n={a['n']}"
                  f" | balanced {b['tau_b']} [{b['ci_lo_2p5']},{b['ci_hi_97p5']}] n={b['n']}")
        for k, v in t["within"].items():
            print("  within", k, v)
        print("  E16 pooled", t["e16_pooled_radius_sep_undefined"])
        print("  E16 within", t["e16_within_radius_sep_undefined"])
        print("  share sep undefined (within)", t["share_within_units_sep_undefined"])
        print("  table13", {c: (round(v["tau_b"], 3), v["lo"], v["hi"], v["n_units"])
                            for c, v in r["table13"].items()})
        print("  paired", r["paired_pooled"])
