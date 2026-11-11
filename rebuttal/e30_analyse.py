"""E30-E33 analysis (registered statistics of PREREGISTRATION.md, section E30-E34).

Writes results/rebuttal/e30/analysis.json and a few CSV tables.
"""

from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

from common import OUT

ODIR = OUT / "e30"
B = 2000
SEED = 20261111
INF = 10 ** 6


def load(path=ODIR / "rows.jsonl"):
    G, R, errs, tos = [], [], [], []
    for line in path.open():
        r = json.loads(line)
        if r.get("learn_timeout"):
            tos.append({k: r[k] for k in ("network", "n", "seed", "learner")})
            continue
        if "error" in r:
            errs.append({k: r[k] for k in ("id", "error")})
            continue
        rows = r.pop("rows")
        r.pop("cpdag", None)
        G.append(r)
        for x in rows:
            R.append({"network": r["network"], "n": r["n"], "seed": r["seed"], "learner": r["learner"], **x})
    G, R = pd.DataFrame(G), pd.DataFrame(R)
    R["ok"] = R.status == "ok"
    R["inv"] = R.ok & (R.valid_truth == False)  # noqa: E712
    R["valid"] = R.ok & (R.valid_truth == True)  # noqa: E712
    R["rnum"] = R.r.map(lambda v: INF if str(v) == "inf" else (np.nan if pd.isna(v) else float(v)))
    if "class_ok" in R:
        R["class_ok"] = R.class_ok.fillna(False).astype(bool)
        R["g0_ok"] = R.g0_ok.fillna(False).astype(bool)
        R["guarantee"] = R.class_ok & R.k_err.notna() & ((R.rnum >= INF) | (R.k_err < R.rnum))
        R["violation"] = R.guarantee & R.inv
    load.timeouts = pd.DataFrame(tos)
    return G, R, errs


def boot_ratio(cells, num, den, rng, draws=B):
    """cells: DataFrame with columns network, seed, num, den (per graph sums).
    Resample seeds within each network; return the pooled ratio's percentile interval."""
    nets = [g for _, g in cells.groupby("network")]
    vals = []
    for _ in range(draws):
        n_s = d_s = 0.0
        for g in nets:
            idx = rng.integers(0, len(g), len(g))
            n_s += g[num].values[idx].sum()
            d_s += g[den].values[idx].sum()
        vals.append(n_s / d_s if d_s else np.nan)
    v = np.array(vals)
    v = v[~np.isnan(v)]
    if not len(v):
        return [None, None], v
    return [round(float(np.quantile(v, 0.025)), 4), round(float(np.quantile(v, 0.975)), 4)], v


def per_graph(R, arm="TRUE"):
    d = R[(R.arm == arm) & (R.n != "oracle")]
    agg = d.groupby(["learner", "n", "network", "seed"]).agg(
        frame=("status", "size"), ok=("ok", "sum"), inv=("inv", "sum"), valid=("valid", "sum"),
        lc=("class_ok", "sum"), g0ok=("g0_ok", "sum"), guar=("guarantee", "sum"),
        viol=("violation", "sum"),
        inv_lc=("inv", lambda s: int((s & d.loc[s.index, "class_ok"]).sum())),
        inv_nlc=("inv", lambda s: int((s & ~d.loc[s.index, "class_ok"]).sum())),
    ).reset_index()
    agg["nlc"] = agg.ok - agg.lc
    return agg


def e30_tables(R, rng):
    agg = per_graph(R)
    out = []
    for (lr, n), c in agg.groupby(["learner", "n"]):
        rec = {"learner": lr, "n": int(n), "graphs": len(c), "frame": int(c.frame.sum()),
               "committed": int(c.ok.sum())}
        for name, num, den in (("I", "inv", "ok"), ("Cm", "ok", "frame"), ("V", "valid", "frame"),
                               ("P_LC", "lc", "ok"), ("P_guar", "guar", "ok"),
                               ("inv_given_LC", "inv_lc", "lc"), ("inv_given_notLC", "inv_nlc", "nlc"),
                               ("P_g0ok", "g0ok", "ok")):
            den_s = c[den].sum()
            rec[name] = round(float(c[num].sum() / den_s), 4) if den_s else None
            rec[name + "_ci"] = boot_ratio(c, num, den, rng)[0] if den_s else [None, None]
        rec["violations"] = int(c.viol.sum())
        rec["n_LC"] = int(c.lc.sum())
        rec["n_notLC"] = int(c.nlc.sum())
        out.append(rec)
    return pd.DataFrame(out).sort_values(["learner", "n"])


def threshold(t, lr, i_max, cm_min):
    s = t[t.learner == lr].sort_values("n")
    ns = list(s.n)
    ok = [(r.I is not None and r.I <= i_max and r.Cm >= cm_min) for r in s.itertuples()]
    for k in range(len(ns)):
        if all(ok[k:]):
            return int(ns[k])
    return None


def diff_ci(agg, lr, n1, n2, num, den, rng):
    a = agg[(agg.learner == lr) & (agg.n == n1)]
    b = agg[(agg.learner == lr) & (agg.n == n2)]
    if not len(a) or not len(b):
        return None
    _, va = boot_ratio(a, num, den, np.random.default_rng(SEED))
    _, vb = boot_ratio(b, num, den, np.random.default_rng(SEED + 1))
    m = min(len(va), len(vb))
    d = va[:m] - vb[:m]
    pt = a[num].sum() / a[den].sum() - b[num].sum() / b[den].sum()
    return {"point": round(float(pt), 4), "ci": [round(float(np.quantile(d, 0.025)), 4),
                                                  round(float(np.quantile(d, 0.975)), 4)]}


def e33(R, arm="TRUE"):
    orc = R[(R.n == "oracle") & (R.arm == arm) & R.ok][["network", "x", "y", "rnum", "z"]].rename(
        columns={"rnum": "r_or", "z": "z_or"})
    d = R[(R.n != "oracle") & (R.arm == arm) & R.ok].merge(orc, on=["network", "x", "y"], how="inner")
    out = []
    for (lr, n), c in d.groupby(["learner", "n"]):
        tau = kendalltau(c.rnum, c.r_or, variant="b").statistic if len(c) > 2 else np.nan
        rec = {"learner": lr, "n": int(n), "pairs": len(c),
               "P_r_equal": round(float((c.rnum == c.r_or).mean()), 4),
               "tau_b": None if np.isnan(tau) else round(float(tau), 4),
               "P_z_equal": round(float((c.z.fillna("") == c.z_or.fillna("")).mean()), 4)}
        # graph-cluster bootstrap of P_r_equal and tau_b
        rng = np.random.default_rng(SEED)
        keys = c.groupby(["network", "seed"]).indices
        bynet = {}
        for (net, s), ix in keys.items():
            bynet.setdefault(net, []).append(ix)
        pe, tb = [], []
        for _ in range(500):
            ix = np.concatenate([lst[j] for lst in bynet.values()
                                 for j in rng.integers(0, len(lst), len(lst))])
            s = c.iloc[ix]
            pe.append((s.rnum == s.r_or).mean())
            if len(s) > 2:
                tb.append(kendalltau(s.rnum, s.r_or, variant="b").statistic)
        rec["P_r_equal_ci"] = [round(float(np.quantile(pe, q)), 4) for q in (0.025, 0.975)]
        tb = np.array([t for t in tb if not np.isnan(t)])
        rec["tau_b_ci"] = [round(float(np.quantile(tb, q)), 4) for q in (0.025, 0.975)] if len(tb) else None
        out.append(rec)
    if not out:
        return pd.DataFrame(columns=["learner", "n"])
    return pd.DataFrame(out).sort_values(["learner", "n"])


def e5_repro(R):
    p = OUT / "e5" / "learned_rows_boot.csv"
    e5 = pd.read_csv(p, low_memory=False)
    e5 = e5[(e5.arm == "TRUE") & (e5.n != "oracle")].copy()
    e5["n"] = e5.n.astype(int)
    e5["seed"] = e5.seed.astype(int)
    mine = R[(R.learner == "pc01d3") & (R.arm == "TRUE") & (R.n.isin([200, 1000, 5000]))
             & (R.seed.isin(range(5)))].copy()
    mine["n"] = mine.n.astype(int)
    m = e5.merge(mine, on=["network", "n", "seed", "x", "y"], suffixes=("_e5", ""), how="outer",
                 indicator=True)
    both = m[m._merge == "both"]
    z5 = both.z_e5.fillna("").astype(str)
    zm = both.z.fillna("").astype(str)
    st = both.status_e5 == both.status
    zeq = (z5 == zm) | ~both.ok
    v5 = both.valid_truth_e5.astype(str)
    vm = both.valid_truth.astype(str)
    veq = (v5 == vm) | ~both.ok
    return {"rows_e5": int(len(e5)), "rows_here": int(len(mine)), "matched": int(len(both)),
            "status_equal": int(st.sum()), "z_equal": int(zeq.sum()), "validity_equal": int(veq.sum())}


def graph_table(G):
    d = G[G.n != "oracle"]
    return d.groupby(["learner", "n"]).agg(
        graphs=("id", "size"), shd=("shd_cpdag", "mean"), skel_rec=("skel_rec", "mean"),
        skel_prec=("skel_prec", "mean"), vs_rec=("vs_rec", "mean"), vs_prec=("vs_prec", "mean"),
        missing=("missing", "mean"), extra=("extra", "mean"), is_cpdag=("is_cpdag", "mean"),
        learn_s=("learn_secs_laptop", "median"), score_s=("score_secs_laptop", "median"),
    ).round(3).reset_index()


def main():
    G, R, errs = load()
    rng = np.random.default_rng(SEED)
    res = {"errors": errs, "graphs": int(len(G)), "rows": int(len(R))}
    to = load.timeouts
    res["learn_timeouts"] = (to.groupby(["learner", "network", "n"]).size().reset_index(name="count")
                             .to_dict("records") if len(to) else [])
    t = e30_tables(R, rng)
    t.to_csv(ODIR / "e30_table.csv", index=False)
    res["e30"] = t.to_dict("records")
    res["thresholds"] = {lr: {"0.05/0.95": threshold(t, lr, 0.05, 0.95),
                              "0.10/0.90": threshold(t, lr, 0.10, 0.90)}
                         for lr in sorted(t.learner.unique())}
    agg = per_graph(R)
    res["delta_I_100k_minus_1k"] = {lr: diff_ci(agg, lr, 100000, 1000, "inv", "ok", rng)
                                    for lr in sorted(agg.learner.unique())}
    res["delta_I_pc01d3_minus_pc01"] = {
        int(n): None for n in (20000, 50000, 100000)}
    for n in (20000, 50000, 100000):
        a = agg[(agg.learner == "pc01d3") & (agg.n == n)]
        b = agg[(agg.learner == "pc01") & (agg.n == n)]
        if len(a) and len(b):
            _, va = boot_ratio(a, "inv", "ok", np.random.default_rng(SEED))
            _, vb = boot_ratio(b, "inv", "ok", np.random.default_rng(SEED + 1))
            m = min(len(va), len(vb))
            res["delta_I_pc01d3_minus_pc01"][n] = {
                "point": round(float(a.inv.sum() / a.ok.sum() - b.inv.sum() / b.ok.sum()), 4),
                "ci": [round(float(np.quantile(va[:m] - vb[:m], q)), 4) for q in (0.025, 0.975)]}
    # theorem violations, every arm
    res["violations_all_arms"] = R.groupby("arm").violation.sum().astype(int).to_dict()
    res["violations_total"] = int(R.violation.sum())
    # per network
    per_net = []
    for (lr, n, net), c in agg.groupby(["learner", "n", "network"]):
        per_net.append({"learner": lr, "n": int(n), "network": net,
                        "I": round(c.inv.sum() / c.ok.sum(), 3) if c.ok.sum() else None,
                        "Cm": round(c.ok.sum() / c.frame.sum(), 3),
                        "P_LC": round(c.lc.sum() / c.ok.sum(), 3) if c.ok.sum() else None,
                        "committed": int(c.ok.sum())})
    pd.DataFrame(per_net).to_csv(ODIR / "e30_per_network.csv", index=False)
    e = e33(R)
    e.to_csv(ODIR / "e33_table.csv", index=False)
    res["e33_TRUE"] = e.to_dict("records")
    e2 = e33(R, "LLM")
    res["e33_LLM"] = e2.to_dict("records")
    gt = graph_table(G)
    gt.to_csv(ODIR / "graph_table.csv", index=False)
    res["graph_table"] = gt.to_dict("records")
    try:
        res["e5_reproduction"] = e5_repro(R)
    except Exception as ex:  # noqa: BLE001
        res["e5_reproduction"] = f"failed: {ex}"
    # oracle reference
    o = R[R.n == "oracle"]
    res["oracle"] = {arm: {"committed": int(c.ok.sum()), "frame": int(len(c)),
                           "invalid": int(c.inv.sum())} for arm, c in o.groupby("arm")}
    (ODIR / "analysis.json").write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps({k: res[k] for k in ("graphs", "rows", "thresholds", "violations_total",
                                          "e5_reproduction", "oracle")}, indent=1, default=str))


if __name__ == "__main__":
    sys.exit(main())
