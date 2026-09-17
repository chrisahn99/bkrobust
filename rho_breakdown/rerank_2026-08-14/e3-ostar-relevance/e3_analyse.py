"""
E3 step 2 - head-to-head. O*-relevance vs edit distance vs the trivial
query-incidence feature, on the pre-registered populations and thresholds.

All CIs are CLUSTER bootstraps resampling the 600 SCMs (members are nested in
SCMs and are strongly dependent). Delta-AUC CIs are PAIRED on the same resample.

Usage: python e3_analyse.py
"""
import json
import os
from collections import Counter, defaultdict

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
ROWS = os.path.join(HERE, "results", "e3_rows.json")
OUT = os.path.join(HERE, "results", "e3_analysis.json")

NBOOT = 1000
RNG = np.random.default_rng(20260814)

BASELINES = ["d_edit", "d_xy_inc", "n_xy_inc"]
CANDIDATES = ["d_ostar", "d_ostar_graded", "d_ostar_any", "d_ostar_dil"]
NEGDIST = ["d_gdist"]          # smaller = more dangerous -> score is -d_gdist
ALL_FEATS = BASELINES + CANDIDATES + NEGDIST


def auc(score, label):
    """Mann-Whitney AUC with ties (midranks). label in {0,1}."""
    score = np.asarray(score, dtype=float)
    label = np.asarray(label, dtype=bool)
    n1, n0 = int(label.sum()), int((~label).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    order = np.argsort(score, kind="mergesort")
    s = score[order]
    ranks = np.empty(len(s), dtype=float)
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        ranks[i:j + 1] = (i + j) / 2.0 + 1.0
        i = j + 1
    r = np.empty(len(s), dtype=float)
    r[order] = ranks
    return float((r[label].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def score_of(rows_idx, feat, data):
    v = data[feat][rows_idx]
    return -v if feat in NEGDIST else v


def cluster_boot_indices(seeds_arr, uniq_seeds, seed_to_idx, rng):
    pick = rng.choice(len(uniq_seeds), size=len(uniq_seeds), replace=True)
    return np.concatenate([seed_to_idx[uniq_seeds[k]] for k in pick])


def analyse_population(data, mask, label_key, name, feats=ALL_FEATS, extra_oracle=None):
    idx = np.flatnonzero(mask)
    y = data[label_key][idx].astype(bool)
    seeds = data["seed"][idx]
    uniq = np.unique(seeds)
    s2i = {s: np.flatnonzero(seeds == s) for s in uniq}

    point = {f: auc(score_of(idx, f, data), y) for f in feats}
    if extra_oracle is not None:
        point["ORACLE_" + extra_oracle] = auc(score_of(idx, extra_oracle, data), y)

    boots = defaultdict(list)
    for _ in range(NBOOT):
        bi = cluster_boot_indices(seeds, uniq, s2i, RNG)
        yb = y[bi]
        if yb.sum() == 0 or (~yb).sum() == 0:
            continue
        cur = {}
        for f in feats:
            cur[f] = auc(score_of(idx, f, data)[bi], yb)
            boots[f].append(cur[f])
        # paired deltas against the two pre-registered comparators
        for c in CANDIDATES:
            boots[f"D_{c}_vs_d_edit"].append(cur[c] - cur["d_edit"])
            boots[f"D_{c}_vs_d_xy_inc"].append(cur[c] - cur["d_xy_inc"])
            boots[f"D_{c}_vs_d_gdist"].append(cur[c] - (-1) * 0 + 0)  # placeholder
            boots[f"D_{c}_vs_d_gdist"][-1] = cur[c] - cur["d_gdist"]

    ci = {k: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
          for k, v in boots.items()}
    deltas = {}
    for c in CANDIDATES:
        for comp in ("d_edit", "d_xy_inc", "d_gdist"):
            key = f"D_{c}_vs_{comp}"
            deltas[key] = dict(point=point[c] - point[comp], ci=ci[key])

    return dict(name=name, n=int(len(idx)), n_pos=int(y.sum()),
                n_scm=int(len(uniq)), auc=point,
                auc_ci={f: ci[f] for f in feats if f in ci}, deltas=deltas)


def main():
    with open(ROWS) as f:
        blob = json.load(f)
    rows = blob["rows"]
    assert not blob["gate_fail"], "integrity gates failed upstream"

    keys_num = ["seed", "rho", "d_edit", "d_xy_inc", "n_xy_inc", "d_gdist",
                "d_ostar", "d_ostar_graded", "d_ostar_any", "d_ostar_dil"]
    data = {k: np.array([r[k] for r in rows], dtype=float) for k in keys_num}
    data["seed"] = np.array([r["seed"] for r in rows], dtype=int)
    for k in ["consistent", "amenable", "silent", "cpdag_amenable"]:
        data[k] = np.array([bool(r[k]) for r in rows])
    data["ostar_changed"] = np.array([1.0 if r["ostar_changed"] else 0.0 for r in rows])
    data["arm"] = np.array([r["arm"] for r in rows])
    data["abs_rel_bias"] = np.array(
        [np.nan if r["abs_rel_bias"] is None else r["abs_rel_bias"] for r in rows])

    out = {"n_rows": len(rows), "n_arms": blob["n_arms"]}

    # ---------------------------------------------------------------- degeneracy
    stmt = {tuple(k.split("|")): v for k, v in blob["stmt_rel"].items()}
    per_arm = defaultdict(Counter)
    for (arm, r, r2, xy, d), n in stmt.items():
        per_arm[arm]["total"] += n
        per_arm[arm][f"r={r}"] += n
        per_arm[arm][f"r2={r2}"] += n
        per_arm[arm][f"xy={xy}"] += n
        per_arm[arm][f"dist={d}"] += n
        per_arm[arm][f"r={r},xy={xy}"] += n
    out["statement_marginals"] = {a: dict(c) for a, c in per_arm.items()}

    # ---------------------------------------------------------------- populations
    gen = data["arm"] == "generic"
    tie = data["arm"] == "tiered"
    pops = {
        "P1_all_generic": gen,
        "P2_consistent_generic": gen & data["consistent"],
        "P3_undetectable_generic": gen & data["consistent"] & data["amenable"],
        "P3_undetectable_tiered": tie & data["consistent"] & data["amenable"],
        "P3_undetectable_both": data["consistent"] & data["amenable"],
    }
    out["populations"] = {}
    for name, m in pops.items():
        out["populations"][name] = analyse_population(
            data, m, "silent", name, extra_oracle="ostar_changed")

    # ---------------------------------------------------------------- per-rho
    out["by_rho"] = {}
    base = gen & data["consistent"] & data["amenable"]
    for rho in (1, 2, 3):
        m = base & (data["rho"] == rho)
        if m.sum() > 0 and data["silent"][m].sum() > 0:
            out["by_rho"][f"rho={rho}"] = analyse_population(
                data, m, "silent", f"generic P3 rho={rho}")

    # ---------------------------------------------------------------- incremental
    # does d_ostar say anything INSIDE the stratum where no statement touches X or Y?
    inc = {}
    for nm, m in [("xy_inc=0", base & (data["d_xy_inc"] == 0)),
                  ("xy_inc=1", base & (data["d_xy_inc"] == 1))]:
        y = data["silent"][m]
        inc[nm] = dict(n=int(m.sum()), n_silent=int(y.sum()),
                       rate=float(y.mean()) if m.sum() else None)
        if y.sum() > 0 and (~y).sum() > 0:
            for f in CANDIDATES + ["d_edit"]:
                inc[nm][f"auc_{f}"] = auc(score_of(np.flatnonzero(m), f, data), y)
    out["stratified_by_xy_incidence"] = inc

    # cross-tab r(k) vs xy-incidence at the STATEMENT level (generic)
    ct = Counter()
    for (arm, r, r2, xy, d), n in stmt.items():
        if arm == "generic":
            ct[(r, xy)] += n
    out["statement_crosstab_r_by_xy_generic"] = {f"r={r},xy={xy}": n for (r, xy), n in sorted(ct.items())}

    # ---------------------------------------------------------------- distance table
    dt = {}
    for d in range(0, 6):
        m = base & (data["d_gdist"] == d)
        if m.sum() == 0:
            continue
        dt[f"min_gdist={d}"] = dict(n=int(m.sum()), n_silent=int(data["silent"][m].sum()),
                                    rate=float(data["silent"][m].mean()))
    out["silent_rate_by_min_graph_distance"] = dt

    # per-STATEMENT distance table at rho=1 (this is the campaign's 0/791 claim)
    r1 = gen & (data["rho"] == 1)
    st = {}
    for d in range(0, 6):
        m = r1 & (data["d_gdist"] == d)
        if m.sum() == 0:
            continue
        st[f"dist={d}"] = dict(n=int(m.sum()), n_silent=int(data["silent"][m].sum()),
                               rate=float(data["silent"][m].mean()))
    out["rho1_all_statements_by_distance_generic"] = st
    r1c = r1 & data["consistent"] & data["amenable"]
    st2 = {}
    for d in range(0, 6):
        m = r1c & (data["d_gdist"] == d)
        if m.sum() == 0:
            continue
        st2[f"dist={d}"] = dict(n=int(m.sum()), n_silent=int(data["silent"][m].sum()),
                                rate=float(data["silent"][m].mean()))
    out["rho1_undetectable_by_distance_generic"] = st2

    # ---------------------------------------------------------------- continuous outcome
    m = base & np.isfinite(data["abs_rel_bias"])
    sp = {}
    for f in ALL_FEATS:
        s = score_of(np.flatnonzero(m), f, data)
        rho_s, pv = spearmanr(s, data["abs_rel_bias"][m])
        sp[f] = dict(spearman=float(rho_s), p=float(pv))
    # cluster-bootstrap CI on the two headline features
    seeds = data["seed"][m]
    uniq = np.unique(seeds)
    s2i = {s: np.flatnonzero(seeds == s) for s in uniq}
    yv = data["abs_rel_bias"][m]
    bs = defaultdict(list)
    idx_m = np.flatnonzero(m)
    for _ in range(NBOOT):
        bi = cluster_boot_indices(seeds, uniq, s2i, RNG)
        for f in ["d_edit", "d_xy_inc", "d_ostar_graded", "d_ostar", "d_gdist"]:
            s = score_of(idx_m, f, data)[bi]
            if np.all(s == s[0]):
                continue
            bs[f].append(spearmanr(s, yv[bi])[0])
    sp_ci = {f: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
             for f, v in bs.items() if len(v) > 50}
    out["spearman_abs_rel_bias_generic_P3"] = dict(
        n=int(m.sum()), point=sp, ci=sp_ci)

    # ---------------------------------------------------------------- verdict
    P = out["populations"]["P3_undetectable_generic"]
    A_ost = P["auc"]["d_ostar_graded"]
    A_edit = P["auc"]["d_edit"]
    A_xy = P["auc"]["d_xy_inc"]
    d1 = P["deltas"]["D_d_ostar_graded_vs_d_edit"]
    d2 = P["deltas"]["D_d_ostar_graded_vs_d_xy_inc"]
    c1 = A_ost >= 0.70
    c2 = (d1["point"] >= 0.10) and (d1["ci"][0] > 0)
    c3 = (d2["point"] >= 0.05) and (d2["ci"][0] > 0)
    verdict = ("SUPPORTED" if (c1 and c2 and c3)
               else "PARTIAL_ABSORBED" if (c1 and c2)
               else "REFUTED")
    out["verdict"] = dict(A_ost=A_ost, A_edit=A_edit, A_xy=A_xy,
                          cond1_auc_ge_070=bool(c1),
                          cond2_beats_edit_by_010=bool(c2),
                          cond3_beats_xy_inc_by_005=bool(c3),
                          delta_vs_edit=d1, delta_vs_xy_inc=d2,
                          oracle_auc=P["auc"].get("ORACLE_ostar_changed"),
                          verdict=verdict)

    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out["verdict"], indent=1))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
