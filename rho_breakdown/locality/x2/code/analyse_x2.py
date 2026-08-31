"""
X2 ANALYSIS -- every fraction with its denominator, every stratum with its n,
every bar as pre-registered (PREREG section 4, as amended by the binding
mitigations M1-M19 of AUDIT.md).

Usage: python analyse_x2.py <results_dir> <out.json>
"""
import json
import os
import sys
from collections import defaultdict

import numpy as np

import x2lib as X

RD = sys.argv[1] if len(sys.argv) > 1 else "../results"
OUT = sys.argv[2] if len(sys.argv) > 2 else "../results/x2_analysis.json"

SCAN_ENS = ["licensed", "large", "k8", "deep", "deep8"]
FAR_BUCKETS = ("1", "2", "3", "ge4")

# M1: per-ensemble G4a bands, fixed in AUDIT.md from archived E1' data BEFORE
# any X2 code existed.  G4a is a REPRODUCTION check, not evidence.
G4A_BANDS = {"original/K4": (0.32, 0.42), "licensed/K4": (0.42, 0.52),
             "large/K4": (0.15, 0.25), "k8/K4": (0.71, 0.81),
             "k8/K8": (0.51, 0.61)}
G4B_MIN = 0.05          # pre-committed before ARM E ran: dynamic range floor
# M7: S5 thresholds, per ensemble, set from the AUDIT's realised yields
S5_SUPPORTED = 20000
S5_WEAK = 3000
CELL_MIN_D = 300        # M19: cells below this are printed, never used in a verdict


def load(name):
    p = os.path.join(RD, name)
    return json.load(open(p)) if os.path.exists(p) else None


def agg_scan(d, arm, dist="dmin"):
    """-> (by_bucket, by_bucket_true, by_bucket_false), summed over SCMs."""
    tot = defaultdict(lambda: defaultdict(int))
    bysplit = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    byop = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for r in d["recs"]:
        for kk, v in r["acc"][arm][dist].items():
            parts = kk.split("|")
            b, t = parts[0], parts[1]
            op = parts[2] if len(parts) > 2 else "all"
            for c, n in v.items():
                tot[b][c] += n
                bysplit[t][b][c] += n
                if t == "F":
                    byop[op][b][c] += n
    return ({b: dict(v) for b, v in tot.items()},
            {t: {b: dict(v) for b, v in bysplit[t].items()} for t in ("T", "F")},
            {o: {b: dict(v) for b, v in byop[o].items()} for o in byop})


def roll(tab, buckets):
    out = defaultdict(int)
    for b in buckets:
        for c, n in tab.get(b, {}).items():
            out[c] += n
    return dict(out)


def derived(c):
    D_all = c.get("D_all", 0)
    D_con = c.get("D_con", 0)
    D_am = c.get("D_am", 0)
    D_s5 = c.get("D_s5", 0)
    N_ident = c.get("N_ident", 0)
    N_gain = c.get("N_gain", 0)
    N_set = c.get("N_set", 0)
    return dict(D_all=D_all, D_con=D_con, D_am=D_am, D_s5=D_s5,
                D_ident=D_am + N_ident, D_gain=D_con - D_am - N_ident,
                N_set=N_set, N_bias=c.get("N_bias", 0), N_ident=N_ident,
                N_ident_nopath=c.get("N_ident_nopath", 0),
                N_ident_unamen=c.get("N_ident_unamen", 0),
                N_gain=N_gain, N_gain_paths=c.get("N_gain_paths", 0),
                N_gain_orient=c.get("N_gain_orient", 0),
                N_meek=c.get("N_meek", 0), N_null=c.get("N_null", 0),
                N_est_inf=c.get("N_est_inf", 0), N_est_n=c.get("N_est_n", 0),
                cascade_pos=c.get("cascade_pos", 0))


def rate(n, d):
    return (n / d) if d else float("nan")


def main():
    A = {}
    step0 = load("step0_archive.json")
    ver = load("verify_armB.json")
    mach = load("machinery_tests.json")

    # ================================================================ GATES
    gates = {}
    gates["machinery_tests"] = dict(
        n=len(mach or []), n_pass=sum(1 for m in (mach or []) if m["ok"]),
        pass_=all(m["ok"] for m in (mach or [])))
    gates["G3_replication"] = {k: v["mismatches"] for k, v in (ver or {}).items()}
    gates["G3_pass"] = all(not v["mismatches"] for v in (ver or {}).values())
    gates["G10_placebo"] = {k: v["armD_placebo_G10"] for k, v in (step0 or {}).items()}
    gates["G10_pass"] = all(all(x["Q"] == 0 for x in v["armD_placebo_G10"].values())
                            for v in (step0 or {}).values())

    # G4a -- reproduction of the archived positive control (M1)
    g4a = {}
    for cellk, band in G4A_BANDS.items():
        c = (step0 or {}).get(cellk)
        if not c:
            continue
        near = c["per_statement_rho1_by_dmin"].get("0", {})
        L = rate(near.get("N_set", 0), near.get("D_am", 0))
        g4a[cellk] = dict(L_set_0=L, band=list(band), n=near.get("D_am", 0),
                          pass_=bool(band[0] <= L <= band[1]))
    gates["G4a_reproduction"] = g4a
    gates["G4a_pass"] = all(v["pass_"] for v in g4a.values())

    # ================================================================ SCAN ARMS
    arms = {}
    for ens in SCAN_ENS:
        d = load(f"scan_{ens}.json")
        if d is None:
            continue
        cell = dict(n_scm=d["n_scm"], n_abort=d["n_abort"], n_draws=d["n_draws"],
                    abort_rate=d["n_abort"] / max(d["n_draws"], 1))
        for arm in ("A", "E"):
            tot, split, byop = agg_scan(d, arm, "dmin")
            totx, _, _ = agg_scan(d, arm, "dmax")
            near = derived(roll(tot, ("0",)))
            far = derived(roll(tot, FAR_BUCKETS))
            disc = derived(roll(tot, ("disc",)))
            far_F = derived(roll(split["F"], FAR_BUCKETS))     # misstatements only
            far_T = derived(roll(split["T"], FAR_BUCKETS))
            nearx = derived(roll(totx, ("0",)))
            farx = derived(roll(totx, FAR_BUCKETS))
            # per-SCM cluster counts for the far D_s5 denominator (M9)
            n_clust = 0
            n_clust_ev = 0
            for r in d["recs"]:
                s = 0
                ev = 0
                for kk, v in r["acc"][arm]["dmin"].items():
                    pp = kk.split("|")
                    b, tt = pp[0], pp[1]
                    if b in FAR_BUCKETS and tt == "F":
                        s += v.get("D_s5", 0)
                        ev += v.get("N_set", 0)
                if s > 0:
                    n_clust += 1
                    n_clust_ev += (ev > 0)
            # radii
            rstar = dict(set=-1.0, ident=-1.0, gain=-1.0)
            for b in ("0",) + FAR_BUCKETS:
                v = tot.get(b, {})
                bv = 4.0 if b == "ge4" else float(b)
                if v.get("N_set", 0):
                    rstar["set"] = max(rstar["set"], bv)
                if v.get("N_ident", 0):
                    rstar["ident"] = max(rstar["ident"], bv)
                if v.get("N_gain", 0):
                    rstar["gain"] = max(rstar["gain"], bv)
            rstar["elicit"] = max(rstar.values())
            cell[f"arm{arm}"] = dict(
                by_dmin={b: derived(v) for b, v in sorted(tot.items())},
                by_dmax={b: derived(v) for b, v in sorted(totx.items())},
                by_dmin_true={b: derived(v) for b, v in sorted(split["T"].items())},
                by_dmin_false={b: derived(v) for b, v in sorted(split["F"].items())},
                by_op={o: dict(near=derived(roll(byop[o], ("0",))),
                               far=derived(roll(byop[o], FAR_BUCKETS)),
                               L_set_near=rate(roll(byop[o], ("0",)).get("N_set", 0),
                                               roll(byop[o], ("0",)).get("D_am", 0)))
                       for o in byop},
                near=near, far=far, disc=disc, far_misstmt=far_F, far_true=far_T,
                near_dmax=nearx, far_dmax=farx,
                L_set_near=rate(near["N_set"], near["D_am"]),
                L_set_far=rate(far["N_set"], far["D_am"]),
                D_s5_far=far_F["D_s5"], N_set_far_misstmt=far_F["N_set"],
                n_clusters_far=n_clust, n_clusters_with_event=n_clust_ev,
                bound_naive=X.rule_of_three(far_F["D_s5"]),
                bound_cluster=X.rule_of_three(n_clust),
                design_effect=rate(far_F["D_s5"], n_clust),
                r_star=rstar,
                composite_near=rate(near["N_set"] + near["N_ident"] + near["N_gain"], near["D_con"]),
                composite_far=rate(far["N_set"] + far["N_ident"] + far["N_gain"], far["D_con"]),
                cascade_near=rate(near["cascade_pos"], near["D_con"]),
                cascade_far=rate(far["cascade_pos"], far["D_con"]),
            )
            # k_base decomposition of ARM E's positive control
        if arm == "E":
            kb = defaultdict(lambda: defaultdict(int))
            for r in d["recs"]:
                for kk, v in r["acc"]["E"]["dmin"].items():
                    if not kk.startswith("0|"):
                        continue
                    for c, n in v.items():
                        kb[r["k_base"]][c] += n
            cell["armE_near_by_k_base"] = {
                str(k): dict(D_am=v.get("D_am", 0), N_set=v.get("N_set", 0),
                             L=rate(v.get("N_set", 0), v.get("D_am", 0)))
                for k, v in sorted(kb.items())}
    # M13 matched within-SCM contrast, ARM E
        m = dict(n_scm=0, D_am_near=0, N_set_near=0, D_am_far=0, N_set_far=0)
        for r in d["recs"]:
            ne, fa = r.get("near_E"), r.get("far_E")
            if ne and fa and ne["D_am"] and fa["D_am"]:
                m["n_scm"] += 1
                m["D_am_near"] += ne["D_am"]
                m["N_set_near"] += ne["N_set"]
                m["D_am_far"] += fa["D_am"]
                m["N_set_far"] += fa["N_set"]
        cell["matched_within_scm_armE"] = m
        # ARM A mechanism statistics (AUDIT A2)
        recs = d["recs"]
        amenC = [r["amen_C"] for r in recs]
        cell["mechanism"] = dict(
            amen_C_rate=float(np.mean(amenC)),
            amen0_E_rate=float(np.mean([bool(r["amen0_E"]) for r in recs])),
            mean_nU=float(np.mean([r["nU"] for r in recs])),
            mean_ncn_given_amenC=float(np.mean([r["ncn_C"] for r in recs if r["amen_C"]]) if any(amenC) else 0),
            mean_u_inc_cn_given_amenC=float(np.mean([r["u_inc_cn_C"] for r in recs if r["amen_C"]]) if any(amenC) else 0),
            frac_amenC_with_any_u_inc_cn=float(np.mean([r["u_inc_cn_C"] > 0 for r in recs if r["amen_C"]]) if any(amenC) else 0),
        )
        # far events actually recorded
        fr = [rr for r in recs for rr in r["far_recs"]]
        cell["far_event_counts"] = dict(
            n=len(fr),
            by_event=dict(sorted({e: sum(1 for z in fr if z["event"] == e)
                                  for e in set(z["event"] for z in fr)}.items())))
        cell["far_set_events"] = [z for z in fr if z["event"] == "set"][:40]
        arms[ens] = cell
    A["scan"] = arms

    # G4b: ARM E dynamic range, pre-committed floor
    g4b = {}
    for ens, cell in arms.items():
        L = cell["armE"]["L_set_near"]
        n = cell["armE"]["near"]["D_am"]
        g4b[ens] = dict(L_set_0=L, n=n, floor=G4B_MIN,
                        pass_=bool(n >= CELL_MIN_D and L >= G4B_MIN))
    gates["G4b_armE_dynamic_range"] = g4b
    gates["G4b_pass"] = all(v["pass_"] for v in g4b.values() if v["n"] >= CELL_MIN_D)
    # G7 sentinel: `disc` is its own row everywhere -- structural (T2)
    gates["G7_sentinel"] = "disc is a separate bucket in every table; T2 proves an edge is wholly in or out"
    # G9 aborts
    gates["G9_aborts"] = {e: arms[e]["abort_rate"] for e in arms}
    gates["G9_pass"] = all(v <= 0.05 for v in gates["G9_aborts"].values())
    # G11 null-trial placebo (ARM E only): a swept TRUE statement already in the
    # base reproduces the base exactly, so it must fire nothing.
    g11 = {}
    for ens, cell in arms.items():
        nnull = sum(v["N_null"] for v in cell["armE"]["by_dmin"].values())
        g11[ens] = dict(n_null=nnull)
    gates["G11_null_trials"] = g11
    A["gates"] = gates

    # ================================================================ S5
    s5 = dict(rho1_by_ensemble={}, rho_ge2_from_archive={}, dmax_headline={})
    for ens, cell in arms.items():
        e = cell["armE"]
        s5["rho1_by_ensemble"][ens] = dict(
            N_set=e["far_misstmt"]["N_set"], D_s5=e["D_s5_far"],
            D_am=e["far_misstmt"]["D_am"], D_con=e["far_misstmt"]["D_con"],
            D_all=e["far_misstmt"]["D_all"],
            n_clusters=e["n_clusters_far"],
            bound_naive=e["bound_naive"], bound_cluster=e["bound_cluster"],
            design_effect=e["design_effect"],
            L_set_near=e["L_set_near"], D_am_near=e["near"]["D_am"],
            verdict=("REFUTED" if e["far_misstmt"]["N_set"] >= 1 else
                     "SUPPORTED" if e["D_s5_far"] >= S5_SUPPORTED else
                     "SUPPORTED-WEAK" if e["D_s5_far"] >= S5_WEAK else "INCONCLUSIVE"))
        s5["dmax_headline"][ens] = dict(
            near_dmax_D_am=e["near_dmax"]["D_am"], near_dmax_N_set=e["near_dmax"]["N_set"],
            far_dmax_D_am=e["far_dmax"]["D_am"], far_dmax_N_set=e["far_dmax"]["N_set"])
    for cellk, c in (step0 or {}).items():
        s5["rho_ge2_from_archive"][cellk] = c["rho_ge2"]
        s5.setdefault("dmax_archive", {})[cellk] = c["per_statement_rho1_by_dmax"]
    n_ref = sum(v["N_set"] for v in s5["rho1_by_ensemble"].values())
    s5["overall_rho1"] = "REFUTED" if n_ref else "SEE PER-ENSEMBLE"
    # rho>=2, all flips far: does the falsifier fire?
    tot2 = defaultdict(int)
    for cellk, c in (step0 or {}).items():
        for kk, v in c["rho_ge2"].items():
            if kk.endswith("_all_far"):
                tot2["N_set"] += int(v.get("N_set", 0))
                tot2["D_am"] += int(v.get("D_am", 0))
                tot2["N_bias"] += int(v.get("N_bias", 0))
    s5["rho_ge2_all_far_pooled"] = dict(tot2)
    s5["S5_as_written_verdict"] = ("REFUTED" if tot2["N_set"] >= 1 else "NOT REFUTED")
    A["S5"] = s5

    # ================================================================ S8
    s8 = dict(armD={}, Qset={})
    for cellk, c in (step0 or {}).items():
        s8["armD"][cellk] = c["armD"]
    for cellk, v in (ver or {}).items():
        s8["Qset"][cellk] = v["Qset"]
        s8.setdefault("mean_n_Ostar", {})[cellk] = v["mean_n_Ostar"]
    verdicts = {}
    for r in ("-1", "0", "1", "2", "3"):
        sup = ref = inc = 0
        used = []
        for cellk, c in s8["armD"].items():
            v = c[r]
            if v["n_scm"] < CELL_MIN_D:
                continue
            used.append(cellk)
            lo, hi = v["Q_ci"]
            cost_ok = v["cost_mean"] <= 0.5
            if lo > 0.01:
                ref += 1
            elif hi < 0.001 and cost_ok:
                sup += 1
            else:
                inc += 1
        verdicts[r] = dict(cells_used=used, n_supported=sup, n_refuted=ref,
                           n_inconclusive=inc,
                           verdict=("REFUTED" if ref else
                                    "SUPPORTED" if sup == len(used) and used else
                                    "INCONCLUSIVE"))
    s8["verdict_by_radius"] = verdicts
    s8["overall"] = ("SUPPORTED at some radius"
                     if any(v["verdict"] == "SUPPORTED" for v in verdicts.values())
                     else "NOT SUPPORTED AT ANY RADIUS")
    A["S8"] = s8

    # ================================================================ S9
    s9 = {}
    for ens, cell in arms.items():
        for arm in ("A", "E"):
            e = cell[f"arm{arm}"]
            rows = {}
            for b in ("0",) + FAR_BUCKETS + ("disc",):
                v = e["by_dmin"].get(b)
                if not v:
                    continue
                rows[b] = dict(D_con=v["D_con"], D_ident=v["D_ident"], D_gain=v["D_gain"],
                               N_ident=v["N_ident"], nopath=v["N_ident_nopath"],
                               unamen=v["N_ident_unamen"],
                               N_gain=v["N_gain"], gain_paths=v["N_gain_paths"],
                               gain_orient=v["N_gain_orient"],
                               rate_ident=rate(v["N_ident"], v["D_ident"]),
                               rate_gain=rate(v["N_gain"], v["D_gain"]))
            s9[f"{ens}/arm{arm}"] = dict(rows=rows, r_star=e["r_star"],
                                         composite_near=e["composite_near"],
                                         composite_far=e["composite_far"])
    rmax = max((v["r_star"]["elicit"] for v in s9.values()), default=-1)
    s9["r_elicit_max_over_ensembles"] = rmax
    s9["verdict"] = ("SUPPORTED" if rmax <= 0 else
                     "PARTIAL (radius 1)" if rmax <= 1 else
                     f"REFUTED (r*_elicit = {rmax:g})")
    A["S9"] = s9

    # ================================================================ M12  S3/S4/S5
    m12 = {}
    for cellk, c in (step0 or {}).items():
        t = c["per_statement_rho1_by_dmin"]
        D_all = sum(v.get("D_all", 0) for v in t.values())
        D_con = sum(v.get("D_con", 0) for v in t.values())
        N_bias = sum(v.get("N_bias", 0) for v in t.values())
        near = t.get("0", {})
        m12[cellk] = dict(
            S3_frac_consistent=rate(D_con, D_all), D_all=D_all, D_con=D_con,
            S4_frac_biasing=rate(N_bias, D_con), N_bias=N_bias,
            P_near_given_consistent=rate(near.get("D_con", 0), D_con),
            L_bias_near=rate(near.get("N_bias", 0), near.get("D_con", 0)),
            product=rate(near.get("D_con", 0), D_con) * rate(near.get("N_bias", 0), near.get("D_con", 0)),
            N_bias_far=N_bias - near.get("N_bias", 0),
            E_bias_ne_E_set=c["E_bias_ne_E_set"], E_est_ne_E_set=c["E_est_inf_ne_E_set"],
            matched=c["matched_within_scm"])
    A["M12_S3_S4_S5_contingency"] = m12

    # ================================================================ lemma
    lem = {}
    for p in (4, 5, 6):
        d = load(f"lemma_p{p}.json")
        if d:
            lem[str(p)] = d
    A["lemma"] = lem

    A["step0_contradiction"] = {k: dict(
        E1_locality=v["E1_locality_by_min_hopdist_reproduced"],
        decomposition=v["noncensored_hopmin_ge1_decomposition"])
        for k, v in (step0 or {}).items()}
    A["armB_ident_split"] = {k: v["ident_split"] for k, v in (ver or {}).items()}

    with open(OUT, "w") as f:
        json.dump(A, f, indent=1, default=float)
    print(f"[analyse] -> {OUT}")

    # ---------------------------------------------------------------- report
    print("\n" + "=" * 78)
    print("GATES")
    print("=" * 78)
    print(f"  machinery tests : {gates['machinery_tests']['n_pass']}/{gates['machinery_tests']['n']} pass")
    print(f"  G3 replication  : {'PASS' if gates['G3_pass'] else 'FAIL'}  "
          f"({sum(v['n_members'] for v in (ver or {}).values())} archived members re-derived)")
    print(f"  G4a reproduction: {'PASS' if gates['G4a_pass'] else 'FAIL'}")
    for k, v in g4a.items():
        print(f"      {k:<14} L_set(0)={v['L_set_0']:.4f} band={v['band']} n={v['n']}  "
              f"{'ok' if v['pass_'] else 'OUT OF BAND'}")
    print(f"  G4b ARM E range : {'PASS' if gates['G4b_pass'] else 'FAIL'}")
    for k, v in g4b.items():
        print(f"      {k:<10} L_set(0)={v['L_set_0']:.4f} n={v['n']} floor={v['floor']}  "
              f"{'ok' if v['pass_'] else 'BELOW FLOOR/UNDERPOWERED'}")
    print(f"  G9 aborts       : {'PASS' if gates['G9_pass'] else 'FAIL'}  {gates['G9_aborts']}")
    print(f"  G10 placebo     : {'PASS' if gates['G10_pass'] else 'FAIL'}")
    print(f"  G11 null trials : {gates['G11_null_trials']}")

    print("\n" + "=" * 78)
    print("S5  rho=1, ARM E (exhaustive statement sweep), MISSTATEMENTS, FINITE dmin>=1")
    print("=" * 78)
    print(f"{'ensemble':<10}{'N_set':>7}{'D_s5':>9}{'D_am':>9}{'D_con':>9}{'D_all':>9}"
          f"{'clusters':>10}{'3/D_s5':>11}{'3/nSCM':>11}  verdict")
    for ens, v in s5["rho1_by_ensemble"].items():
        print(f"{ens:<10}{v['N_set']:>7}{v['D_s5']:>9}{v['D_am']:>9}{v['D_con']:>9}"
              f"{v['D_all']:>9}{v['n_clusters']:>10}{v['bound_naive']:>11.2e}"
              f"{v['bound_cluster']:>11.2e}  {v['verdict']}")
    print("\n  dmax reading (M10) -- the OTHER reading of 'graph-distance >= 1':")
    print(f"{'ensemble':<10}{'N_set(dmax=0)':>15}{'D_am(dmax=0)':>14}"
          f"{'N_set(dmax>=1)':>16}{'D_am(dmax>=1)':>15}")
    for ens, v in s5["dmax_headline"].items():
        print(f"{ens:<10}{v['near_dmax_N_set']:>15}{v['near_dmax_D_am']:>14}"
              f"{v['far_dmax_N_set']:>16}{v['far_dmax_D_am']:>15}")
    print(f"\n  rho>=2 (ARM B archive, ALL flips far, finite): "
          f"N_set={tot2['N_set']} N_bias={tot2['N_bias']} D_am={tot2['D_am']}  "
          f"-> S5 as written: {s5['S5_as_written_verdict']}")

    print("\n" + "=" * 78)
    print("S8  ARM D: Q_r (interval) and Qset_r (reachable O* set)")
    print("=" * 78)
    for cellk, c in s8["armD"].items():
        qs = s8["Qset"].get(cellk, {})
        print(f"  {cellk}")
        for r in ("-1", "0", "1", "2", "3"):
            v = c[r]
            q2 = qs.get(r, {})
            print(f"    r={r:>2} n={int(v['n_scm']):5d} Q={v['Q']:.4f} "
                  f"CI=[{v['Q_ci'][0]:.5f},{v['Q_ci'][1]:.5f}] cost={v['cost_mean']:.3f} "
                  f"| Qset={q2.get('Q', float('nan')):.4f} "
                  f"CI=[{q2.get('Q_ci',[float('nan')]*2)[0]:.5f},{q2.get('Q_ci',[float('nan')]*2)[1]:.5f}]")
    for r, v in s8["verdict_by_radius"].items():
        print(f"  radius {r:>2}: {v['verdict']}  (cells used: {len(v['cells_used'])})")
    print(f"  OVERALL S8: {s8['overall']}")

    print("\n" + "=" * 78)
    print("S9  elicitation radius")
    print("=" * 78)
    for k, v in s9.items():
        if not isinstance(v, dict) or "r_star" not in v:
            continue
        print(f"  {k}: r*_set={v['r_star']['set']:g} r*_ident={v['r_star']['ident']:g} "
              f"r*_gain={v['r_star']['gain']:g} r*_elicit={v['r_star']['elicit']:g}")
        for b, row in v["rows"].items():
            print(f"      dmin={b:>4} D_con={row['D_con']:>7} "
                  f"N_ident={row['N_ident']:>5} (nopath {row['nopath']}/unamen {row['unamen']}) "
                  f"N_gain={row['N_gain']:>5} (paths {row['gain_paths']}/orient {row['gain_orient']}) "
                  f"D_gain={row['D_gain']:>7}")
    print(f"  r*_elicit over all ensembles = {rmax:g}  -> S9: {s9['verdict']}")

    print("\n" + "=" * 78)
    print("M12  S3 / S4 / S5 are ONE contingency table (ARM B archive)")
    print("=" * 78)
    for cellk, v in m12.items():
        print(f"  {cellk}: S3={v['S3_frac_consistent']:.4f} ({v['D_con']}/{v['D_all']})  "
              f"S4={v['S4_frac_biasing']:.4f} ({v['N_bias']}/{v['D_con']})  "
              f"= P(near|cons)={v['P_near_given_consistent']:.4f} x "
              f"L_bias(near)={v['L_bias_near']:.4f} = {v['product']:.4f}  "
              f"N_bias(far)={v['N_bias_far']}")

    print("\n" + "=" * 78)
    print("LEMMA  exhaustive small-p existence search")
    print("=" * 78)
    for p, d in lem.items():
        print(f"  p={p} ({d['mode']}): {d['n_cpdags']} CPDAGs, {d['n_mpdags']} MPDAGs, "
              f"far E_set hits = {d['n_hits_far']}")
        for b, v in d["by_dmin"].items():
            print(f"      dmin={b:>4}: " + " ".join(f"{a}={int(c)}" for a, c in sorted(v.items())))


if __name__ == "__main__":
    main()
