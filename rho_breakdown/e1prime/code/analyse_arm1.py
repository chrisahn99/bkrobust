"""
ARM 1 analysis: the pre-registered decision rule, the gates, the placebos.

Usage: python analyse_arm1.py ../results/arm1_licensed.json licensed
"""
import json
import sys
from collections import Counter

import numpy as np
from scipy.stats import spearmanr

from se import (Z95, check_monotone_n, check_monotone_rho, deviation_radius,
                ident_radius, ratio_curve, rho_star_oracle, rho_star_se)

# the registered falsifier's own grid, + diagnostic extensions
N_REG = [200, 500, 1000, 2000, 5000, 10000, 20000, np.inf]
N_DIAG = [50, 100, 50000, 100000, 1000000]
N_GRID = sorted(set(N_REG) | set(N_DIAG))
ZS = [1.0, Z95, 2.5758293035489004]
TAU_FLOOR = 1e-3
RNG = np.random.default_rng(20260819)


def boot_ci(v, B=2000):
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return (float("nan"),) * 3
    idx = RNG.integers(0, len(v), size=(B, len(v)))
    bs = v[idx].mean(axis=1)
    return (float(v.mean()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)))


def se_at(se_factor, k_reg, n):
    if np.isinf(n):
        return 0.0
    if n <= k_reg + 1:
        return float("inf")
    return se_factor / np.sqrt(n - k_reg - 1)


def profile(rhos, max_rho):
    """bins 1..max_rho plus censored (max_rho+1)."""
    rhos = np.asarray(rhos)
    n = len(rhos)
    bins = {str(b): round(float((rhos == b).sum()) / n, 4)
            for b in range(1, max_rho + 2)}
    cens = float((rhos == max_rho + 1).mean())
    return dict(n=n, bins=bins, censored=round(cens, 4),
                censored_ci=[round(x, 4) for x in boot_ci((rhos == max_rho + 1).astype(float))],
                max_bin=round(max(bins.values()), 4),
                n_bins_ge_5pct=sum(1 for v in bins.values() if v >= 0.05),
                median=float(np.median(rhos)))


def main():
    path, label = sys.argv[1], sys.argv[2]
    key = sys.argv[3] if len(sys.argv) > 3 else "K4"
    res = [r for r in json.load(open(path)) if key in r["arms"]
           and r["arms"][key].get("mpdag_amenable")]
    res_full = [r for r in res if abs(r["tau"]) >= TAU_FLOOR]

    S = {"source": path, "label": label, "arm_key": key,
         "n_scm_amenable": len(res), "n_scm_analysed": len(res_full),
         "tau_floor": TAU_FLOOR,
         "p_dist": dict(sorted(Counter(r["p"] for r in res_full).items())),
         "n_und": {"mean": round(float(np.mean([r["n_undirected"] for r in res_full])), 3),
                   "min": int(min(r["n_undirected"] for r in res_full)),
                   "max": int(max(r["n_undirected"] for r in res_full))}}

    max_rho = res_full[0]["max_rho"]
    D, DF, SEF, KREG, EST0, TAU, IDENT = [], [], [], [], [], [], []
    NCONS, NAMEN, NMEMB = [], [], []
    ORACLE = {c: [] for c in ("sign", "rel50", "ident")}
    HOPMIN = []
    for r in res_full:
        a = r["arms"][key]
        top = min(max_rho, a["n_K"])
        D.append(deviation_radius(a["members"], a["est0"], top))
        DF.append(deviation_radius(a["members_frozen"], a["est0"], top))
        SEF.append(a["se_factor"]); KREG.append(a["k_reg"])
        EST0.append(a["est0"]); TAU.append(r["tau"])
        IDENT.append(ident_radius(a["members"], top))
        for c in ORACLE:
            ORACLE[c].append(rho_star_oracle(a["members"], r["tau"], a["est0"], top, c))
        HOPMIN.append(min(a["stmt_hopdist"]))
        ms = [m for m in a["members"] if m["rho"] <= top]
        NMEMB.append(len(ms))
        NCONS.append(sum(1 for m in ms if m.get("consistent")))
        NAMEN.append(sum(1 for m in ms if m.get("consistent") and m.get("amenable")))
    D = np.array(D); DF = np.array(DF)
    SEF = np.array(SEF); KREG = np.array(KREG)
    EST0 = np.array(EST0); TAU = np.array(TAU); HOPMIN = np.array(HOPMIN)
    IDENT = np.array(IDENT)
    NCONS = np.array(NCONS); NAMEN = np.array(NAMEN); NMEMB = np.array(NMEMB)
    DLAST = D[:, -1]

    # ---------------- gates
    G = {}
    G["G4_monotone_rho"] = bool(all(check_monotone_rho(d) for d in D))
    # G1: silent-bias rate at rho=1
    # NOTE 2026-08-19: the denominator is the PILOT'S -- consistent members only
    # (E1 RESULTS.md 6: "on the pilot's own denominator (consistent members)").
    # The first cut of this gate divided by ALL rho=1 members and read low by
    # exactly the Meek-inconsistency rate; both denominators are reported.
    n1 = s1 = c1 = 0
    for r in res_full:
        for m in r["arms"][key]["members"]:
            if m["rho"] != 1:
                continue
            n1 += 1
            if m.get("consistent"):
                c1 += 1
            if m.get("consistent") and m.get("amenable") \
               and m.get("ostar_changed") and not m.get("ostar_valid"):
                s1 += 1
    G["G1_silent_bias_rho1"] = round(s1 / c1, 4) if c1 else None
    G["G1_silent_bias_all_members_denom"] = round(s1 / n1, 4) if n1 else None
    G["G1_meek_inconsistent_rho1"] = round(1 - c1 / n1, 4) if n1 else None
    G["G1_pass"] = bool(c1 and 0.12 <= s1 / c1 <= 0.18)
    G["changed_still_valid"] = int(sum(
        1 for r in res_full for m in r["arms"][key]["members"]
        if m.get("consistent") and m.get("amenable")
        and m.get("ostar_changed") and m.get("ostar_valid")))

    # ---------------- the sweep
    sweep, sweep_frozen, mono = {}, {}, []
    per_scm_rho, per_scm_any = {}, {}
    for z in ZS:
        zk = f"z{z:.3f}"
        sweep[zk], sweep_frozen[zk] = {}, {}
        for n in N_GRID:
            ses = np.array([se_at(SEF[i], KREG[i], n) for i in range(len(D))])
            rhos = np.array([rho_star_se(D[i], ses[i], z, min(max_rho, len(D[i]) - 1))
                             for i in range(len(D))])
            rf = np.array([rho_star_se(DF[i], ses[i], z, min(max_rho, len(DF[i]) - 1))
                           for i in range(len(DF))])
            nk = "inf" if np.isinf(n) else str(int(n))
            anyr = np.minimum(rhos, IDENT)
            sweep[zk][nk] = profile(rhos, max_rho)
            sweep[zk][nk]["ANY"] = profile(anyr, max_rho)
            sweep_frozen[zk][nk] = profile(rf, max_rho)
            # censored-bin decomposition (ADDENDUM 1)
            for tag, rr in (("se", rhos), ("any", anyr)):
                cens = rr == max_rho + 1
                nc = int(cens.sum())
                sweep[zk][nk][f"censored_decomp_{tag}"] = dict(
                    n_censored=nc,
                    all_caught_free=round(float((cens & (NCONS == 0)).sum() / nc), 4) if nc else None,
                    all_loud=round(float((cens & (NCONS > 0) & (NAMEN == 0)).sum() / nc), 4) if nc else None,
                    genuinely_robust=round(float((cens & (NAMEN > 0)).sum() / nc), 4) if nc else None,
                    of_which_d_exactly_zero=round(float(
                        (cens & (NAMEN > 0) & (DLAST == 0)).sum() / nc), 4) if nc else None)
            if abs(z - Z95) < 1e-9:
                per_scm_rho[nk] = rhos
                per_scm_any[nk] = anyr
                R = np.array([ratio_curve(D[i], ses[i], z)[min(1, len(D[i]) - 1)]
                              for i in range(len(D))])
                Rf = R[np.isfinite(R)]
                sweep[zk][nk]["ratio_rho1"] = dict(
                    mean=round(float(Rf.mean()), 4) if len(Rf) else None,
                    median=round(float(np.median(Rf)), 4) if len(Rf) else None,
                    q10=round(float(np.percentile(Rf, 10)), 4) if len(Rf) else None,
                    q90=round(float(np.percentile(Rf, 90)), 4) if len(Rf) else None,
                    frac_gt1=round(float((R > 1).mean()), 4),
                    frac_inf=round(float(np.isinf(R).mean()), 4),
                    frac_exactly_zero=round(float((R == 0).mean()), 4))
    # G3 monotone in n
    ns = [n for n in N_GRID]
    keys = ["inf" if np.isinf(n) else str(int(n)) for n in ns]
    M = np.array([per_scm_rho[k] for k in keys])          # (n_grid, n_scm)
    ordn = np.argsort([1e18 if np.isinf(n) else n for n in ns])
    M = M[ordn]
    G["G3_monotone_n"] = bool(np.all(np.diff(M, axis=0) <= 0))
    G["G3_violations"] = int((np.diff(M, axis=0) > 0).sum())
    # G5 placebo: frozen ball must be 100% censored at every n
    G["G5_placebo_frozen_censored"] = sorted({sweep_frozen[f"z{Z95:.3f}"][k]["censored"]
                                              for k in sweep_frozen[f"z{Z95:.3f}"]})
    G["G5_pass"] = bool(G["G5_placebo_frozen_censored"] == [1.0])
    # G6 floor: SE = 0
    G["G6_floor_censored_at_inf"] = sweep[f"z{Z95:.3f}"]["inf"]["censored"]
    G["G6_frac_d_exactly_zero"] = round(float((D[:, -1] == 0).mean()), 4)
    G["G6_pass"] = bool(abs(G["G6_floor_censored_at_inf"] - G["G6_frac_d_exactly_zero"]) < 1e-9)

    # ---------------- oracle criteria, for continuity with E1
    S["oracle"] = {c: profile(np.array(v), max_rho) for c, v in ORACLE.items()}

    # ---------------- honest-risk comparator (PREREG 3c)
    n2k = per_scm_rho["2000"]
    ses2k = np.array([se_at(SEF[i], KREG[i], 2000) for i in range(len(D))])
    tstat = np.abs(EST0) / np.where(ses2k > 0, ses2k, np.nan)
    ok = np.isfinite(tstat)
    rs = spearmanr(n2k[ok], tstat[ok])
    S["honest_risk_vs_tstat_n2000"] = dict(
        spearman=round(float(rs.statistic), 4), n=int(ok.sum()),
        collapses_to_one_sentence=bool(abs(rs.statistic) > 0.90))

    a2k = per_scm_any["2000"]
    rs2 = spearmanr(a2k, IDENT)
    S["honest_risk_any_vs_ident_n2000"] = dict(
        spearman=round(float(rs2.statistic), 4),
        se_half_contributes_nothing=bool(abs(rs2.statistic) > 0.95),
        frac_any_strictly_below_ident=round(float((a2k < IDENT).mean()), 4))

    # ---------------- locality stratification (the D6 mechanism)
    loc = {}
    for nk in ("2000", "20000", "inf"):
        rr = per_scm_any[nk]
        loc[nk] = {}
        for h in sorted(set(HOPMIN.tolist())):
            sel = HOPMIN == h
            if sel.sum() < 10:
                continue
            loc[nk][str(h)] = dict(n=int(sel.sum()),
                                   censored=round(float((rr[sel] == max_rho + 1).mean()), 4))
    S["locality_by_min_hopdist"] = loc

    # ---------------- pre-registered decision
    z95k = f"z{Z95:.3f}"
    def first_n(pred, sub=None):
        for n in sorted([x for x in N_GRID if not np.isinf(x)]) + [np.inf]:
            nk = "inf" if np.isinf(n) else str(int(n))
            p = sweep[z95k][nk]
            if sub:
                p = p[sub]
            if pred(p):
                return nk
        return None

    def le(a, b):
        return a is not None and (np.inf if a == "inf" else float(a)) <= b

    def rule(sub):
        n50 = first_n(lambda p: p["censored"] <= 0.50, sub)
        nND = first_n(lambda p: p["max_bin"] <= 0.60 and p["n_bins_ge_5pct"] >= 3, sub)
        if le(n50, 20000) and le(nND, 20000):
            v = "SUPPORTED"
        elif (not le(n50, 1e6)) or (not le(nND, 1e6)):
            v = "REFUTED"
        else:
            v = "INCONCLUSIVE"
        g = lambda nk: (sweep[z95k][nk][sub] if sub else sweep[z95k][nk])
        return dict(n50=n50, nND=nND, verdict=v,
                    censored_at_50=g("50")["censored"],
                    censored_at_200=g("200")["censored"],
                    censored_at_20000=g("20000")["censored"],
                    censored_at_inf=g("inf")["censored"],
                    bar_was_free=bool(g("200")["censored"] <= 0.50),
                    registered_falsifier_fires=bool(all(
                        g("inf" if np.isinf(n) else str(int(n)))["censored"] > 0.755
                        for n in N_REG)))

    S["PRIMARY_rho_any"] = rule("ANY")          # ADDENDUM 1 primary
    S["SECONDARY_rho_se"] = rule(None)          # the isolating statistic
    S["rho_ident_alone"] = dict(
        censored=S["oracle"]["ident"]["censored"],
        max_bin=S["oracle"]["ident"]["max_bin"],
        n_bins_ge_5pct=S["oracle"]["ident"]["n_bins_ge_5pct"])
    S["verdict"] = S["PRIMARY_rho_any"]["verdict"]
    S["verdicts_disagree"] = bool(
        S["PRIMARY_rho_any"]["verdict"] != S["SECONDARY_rho_se"]["verdict"])
    S["registered_falsifier_fires"] = S["PRIMARY_rho_any"]["registered_falsifier_fires"]

    S["gates"] = G
    S["sweep"] = sweep
    S["sweep_placebo_frozen"] = {k: {kk: dict(censored=vv["censored"])
                                     for kk, vv in v.items()}
                                 for k, v in sweep_frozen.items()}
    print(json.dumps({k: v for k, v in S.items() if k != "sweep"}, indent=1))
    with open(f"../results/arm1_analysis_{label}_{key}.json", "w") as f:
        json.dump(S, f, indent=1)


if __name__ == "__main__":
    main()
