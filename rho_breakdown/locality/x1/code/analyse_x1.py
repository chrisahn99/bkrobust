"""
X1 ANALYSIS.  Every fraction carries its denominator, its n, and its conditioning.
No stratum is suppressed for small n.  The sentinel is a labelled row, never a number.

Usage: python analyse_x1.py <ensemble>
Writes ../results/x1_analysis_<ensemble>.json
"""
import json
import sys
from collections import defaultdict

import numpy as np

from stats_x1 import (Z95, cliffs_delta, cluster_bootstrap_diff,
                      cluster_bootstrap_rr, design_effect, rate, standardise,
                      wilson)
from x1_ops import UNREACH

N_GRID = [200, 500, 1000, 2000, 5000, 10000, 20000, float("inf")]
Z_GRID = [1.0, Z95, 2.5758293035489004]
BOOT = 10000
LEGACY = "754/8085 = 0.0933 [0.0871,0.0998] -- LEGACY: pilot grid, generic+tiered arms, rho<=3, |K| in {3,4}, member unit, DISJOINT SCM draw (default_rng(20260717)). Quoted once; never the comparator."


def hoplab(h):
    """AUDIT M11 / gate G10: the sentinel is a LABEL, never a number."""
    return "unreachable" if h >= UNREACH else str(int(h))


def silent(m):
    """E1' definition kept VERBATIM (AUDIT A13#5: the ostar_changed clause is
    empirically null because O0 is always valid; equality is stated, not assumed)."""
    return bool(m.get("consistent") and m.get("amenable")
                and m.get("ostar_valid") is False)


def silent_with_changed(m):
    return bool(m.get("consistent") and m.get("amenable")
                and m.get("ostar_changed") and m.get("ostar_valid") is False)


def loud(m):
    return bool(m.get("consistent") and not m.get("amenable"))


def mpdag_valid(m):
    return bool(m.get("ext", True) and not m.get("cycle", False))


# --------------------------------------------------------------------- helpers
def per_scm(scms, arm, pred_num, pred_den, rho=1, extra=None):
    """Returns (num_per_scm, den_per_scm) -- the unit the cluster bootstrap needs."""
    a = np.zeros(len(scms)); b = np.zeros(len(scms))
    for i, s in enumerate(scms):
        for m in s[arm]:
            if rho is not None and m["rho"] != rho:
                continue
            if extra is not None and not extra(s, m):
                continue
            if not pred_den(m):
                continue
            b[i] += 1
            if pred_num(m):
                a[i] += 1
    return a, b


def tot(a, b):
    return rate(int(a.sum()), int(b.sum()))


# ------------------------------------------------------------------------ main
def main():
    ens = sys.argv[1]
    raw = json.load(open(f"../results/x1_{ens}.json"))
    scms = raw["scms"]
    meta = raw["meta"]
    out = dict(ensemble=ens, meta=meta)
    n = len(scms)

    # ============================ A17: the header numbers, BEFORE any verdict
    hop1_uni = sum(1 for s in scms for m in s["Suni"]
                   if m["rho"] == 1 and 1 <= m["hopC"] < UNREACH)
    hop1_unreach = sum(1 for s in scms for m in s["Suni"]
                       if m["rho"] == 1 and m["hopC"] >= UNREACH)
    out["header"] = dict(
        n_analysed=n,
        verdict_readable=(n >= 1800) if ens == "original" else True,
        n_stmt_rho1_per_arm=sum(len([m for m in s["Suni"] if m["rho"] == 1])
                                for s in scms),
        n_hop_ge1_reachable_Suni_rho1=hop1_uni,
        n_hop_unreachable_Suni_rho1=hop1_unreach,
        P_O0_empty=rate(sum(1 for s in scms if len(s["O0"]) == 0), n),
        mean_O0=float(np.mean([len(s["O0"]) for s in scms])),
        drops=meta["drops"])

    # ============================ GATES
    g = {}
    # G1 -- structural, not measured: O0/est0/se_factor0/k_reg0/sigma2 are
    # computed ONCE per SCM and passed to every ball.  There is no per-arm anchor
    # that could diverge.  Declared, not claimed as evidence.
    g["G1_anchor_shared"] = dict(verdict="STRUCTURAL",
                                 note="one anchor computed per SCM, passed to all balls; "
                                      "not a measurement")
    # G3 -- definitional gates (M7: S-uni ONLY)
    n_uni = sum(len(s["Suni"]) for s in scms)
    g["G3_S2_catch_rate_Suni"] = dict(value=0.0, n=n_uni, verdict="PASS",
                                      note="definitional: bk_assert never raises")
    det_uni = sum(1 for s in scms for t in range(4)
                  if s["lab_uni"][t][0] == "pure_spurious")
    g["G3_skeleton_detect_Suni"] = dict(**rate(det_uni, 4 * n),
                                        verdict="PASS" if det_uni == 4 * n else "FAIL")
    det_loc = sum(1 for s in scms for t in range(4)
                  if s["lab_loc"][t][0] == "pure_spurious")
    g["G3_skeleton_detect_Sloc_MEASURED"] = dict(
        **rate(det_loc, 4 * n), verdict="MEASURED (M7: not a gate)")
    # G4 P-frozen
    g["G4_Pfrozen"] = dict(max_abs_dev=max(s["Pfrozen"]["max_abs_dev"] for s in scms),
                           ostar_changed=sum(s["Pfrozen"]["ostar_changed"] for s in scms),
                           silent=sum(s["Pfrozen"]["silent"] for s in scms),
                           n=sum(s["Pfrozen"]["n"] for s in scms))
    g["G4_Pfrozen"]["verdict"] = "PASS" if (g["G4_Pfrozen"]["max_abs_dev"] == 0.0
                                           and g["G4_Pfrozen"]["ostar_changed"] == 0
                                           and g["G4_Pfrozen"]["silent"] == 0) else "FAIL"
    # G5a / G5b P-null (M6)
    ida = [m for s in scms for m in s["Pnull_id"]]
    g["G5a_Pnull_identity"] = dict(n=len(ida),
                                   G_eq_G0=sum(m["G_eq_G0"] for m in ida),
                                   ostar_changed=sum(m["ostar_changed"] for m in ida),
                                   silent=sum(1 for m in ida
                                              if m["amenable"] and not m["ostar_valid"]))
    g["G5a_Pnull_identity"]["verdict"] = "PASS" if (
        g["G5a_Pnull_identity"]["G_eq_G0"] == len(ida)
        and g["G5a_Pnull_identity"]["ostar_changed"] == 0
        and g["G5a_Pnull_identity"]["silent"] == 0) else "FAIL"
    ap = [s["Pnull_app"] for s in scms if s["Pnull_app"] is not None]
    g["G5b_Pnull_append"] = dict(n=len(ap), n_scm_without_directed_edge=n - len(ap),
                                 G_eq_G0=sum(m["G_eq_G0"] for m in ap),
                                 ostar_changed=sum(m["ostar_changed"] for m in ap))
    g["G5b_Pnull_append"]["verdict"] = "PASS" if (
        g["G5b_Pnull_append"]["G_eq_G0"] == len(ap)
        and g["G5b_Pnull_append"]["ostar_changed"] == 0) else "FAIL"
    # G6 monotone d(rho), per SCM per arm
    viol = 0; nchk = 0
    for s in scms:
        for arm in ("R", "Suni", "Sloc", "Drop"):
            d = np.zeros(s["max_rho"] + 1)
            for r in range(1, s["max_rho"] + 1):
                best = d[r - 1]
                for m in s[arm]:
                    if m["rho"] == r and m.get("consistent") and m.get("amenable"):
                        best = max(best, abs(m["est"] - s["est0"]))
                d[r] = best
            nchk += 1
            if np.any(np.diff(d) < -1e-12):
                viol += 1
    g["G6_monotone_d_rho"] = dict(n=nchk, violations=viol,
                                  verdict="PASS" if viol == 0 else "FAIL")
    # G7 re-stamp no-op  /  G8 closure order
    pre = cur = sa = nS = 0
    for s in scms:
        for arm in ("Suni", "Sloc"):
            for m in s[arm]:
                nS += 1
                pre += int(m.get("restamp_pre", False))
                cur += int(m.get("restamp_cur", False))
                sa += int(m.get("stampall_differs", False))
    g["G7_restamp_noop"] = dict(n=nS, changed_prefix=pre, changed_current=cur,
                                verdict="PASS" if pre == 0 and cur == 0 else "FAIL")
    g["G8_closure_order_MEASUREMENT"] = dict(**rate(sa, nS),
                                             note="close-after-each vs stamp-all-then-close; "
                                                  "no bar, this is a measurement")
    # G10 sentinel discipline: every hop stratum key is a label
    g["G10_sentinel"] = dict(verdict="PASS",
                             note="hop strata keyed by hoplab(); 'unreachable' is its own row")
    # G12 -- the added edge alone (needs ARM D, which the PREREG does not have)
    ge = gr = gn = 0
    for s in scms:
        dmap = {(m["rho"], tuple(m["flip"])): m for m in s["Drop"]}
        for m in s["Suni"]:
            if m["hopC"] < UNREACH:
                continue
            d = dmap[(m["rho"], tuple(m["flip"]))]
            gn += 1
            gr += int(bool(m.get("ostar_changed", False)))
            if m.get("O") != d.get("O") or m.get("amenable") != d.get("amenable"):
                ge += 1
    g["G12_unreachable_structural_zero"] = dict(
        n=gn, edge_attributable_Ostar_moves=ge, raw_moves_vs_G0=gr,
        verdict="PASS" if ge == 0 else "FAIL",
        note="edge-attributable = O*(S) != O*(D) at matched (rho, flip). "
             "raw = O* != O0, which also contains the WITHDRAWN true statement.")
    out["gates"] = g

    # ============================ M13 design effects
    de = {}
    for arm in ("R", "Suni", "Sloc", "Drop"):
        a, b = per_scm(scms, arm, silent, lambda m: True, rho=1)
        de[arm] = design_effect(a, 4)
    out["design_effects_rho1_silent"] = de
    out["design_effect_note"] = ("AUDIT M13: > 1.25 forces every CI to the "
                                 "cluster bootstrap. All primary contrasts use "
                                 "the paired SCM-level cluster bootstrap regardless.")

    # ============================ X1-P1 catchability
    p1 = {}
    p1["R_meek_inconsistent_rho1"] = tot(*per_scm(
        scms, "R", lambda m: not m["consistent"], lambda m: True, rho=1))
    p1["R_meek_inconsistent_rho_le4"] = tot(*per_scm(
        scms, "R", lambda m: not m["consistent"], lambda m: True, rho=None))
    for arm in ("Suni", "Sloc"):
        d = {}
        d["S2_catch"] = dict(rate=0.0, n=sum(len(s[arm]) for s in scms),
                             note="DEFINITIONAL (PREREG 5.1)")
        d["S1_caught_intrinsic_rho1"] = tot(*per_scm(
            scms, arm, lambda m: not m["consistent_S1"], lambda m: True, rho=1))
        # M4.3 four-way decomposition, each with its n
        conf = cyc = nonext_nocyc = 0; nn = 0
        vsC = vs_shield_only = vs_gained = vs_lost_other = 0
        for s in scms:
            for m in s[arm]:
                if m["rho"] != 1:
                    continue
                nn += 1
                if m["conflict"]:
                    conf += 1
                elif m["cycle"]:
                    cyc += 1
                elif not m["ext"]:
                    nonext_nocyc += 1
                if m["vsC"]:
                    vsC += 1
                    if m["vs_gained"] == 0 and m["vs_lost_other"] == 0 \
                       and m["vs_lost_shield"] > 0:
                        vs_shield_only += 1
                    if m["vs_gained"] > 0:
                        vs_gained += 1
                    if m["vs_lost_other"] > 0:
                        vs_lost_other += 1
        d["S1_decomposition_rho1"] = dict(
            n=nn, conflict=rate(conf, nn), cycle=rate(cyc, nn),
            non_extendable_acyclic=rate(nonext_nocyc, nn))
        d["vstruct_vs_C_DIAGNOSTIC_rho1"] = dict(
            n=nn, fires=rate(vsC, nn),
            purely_collider_of_C_shielded_by_added_edge=rate(vs_shield_only, nn),
            any_new_unshielded_collider=rate(vs_gained, nn),
            lost_collider_not_explained_by_added_edge=rate(vs_lost_other, nn),
            note="AUDIT M4: this predicate is contaminated by a purely skeletal "
                 "artefact and is NOT the S1 predicate. Diagnostic only.")
        d["mpdag_valid_rho1"] = tot(*per_scm(scms, arm, mpdag_valid,
                                             lambda m: True, rho=1))
        d["mpdag_valid_by_rho"] = {
            str(r): tot(*per_scm(scms, arm, mpdag_valid, lambda m: True, rho=r))
            for r in range(1, 5)}
        d["path_blowup"] = tot(*per_scm(scms, arm, lambda m: m.get("blowup", False),
                                        lambda m: True, rho=None))
        p1[arm] = d
    out["P1_catchability"] = p1

    # ============================ X1-P2 the headline
    def damage_block(arm_s, arm_r, restrict=None, label=""):
        ns, ds = per_scm(scms, arm_s, silent, lambda m: True, rho=1, extra=restrict)
        nr, dr = per_scm(scms, arm_r, silent, lambda m: m["consistent"], rho=1,
                         extra=restrict)
        nru, dru = per_scm(scms, arm_r, silent, lambda m: True, rho=1, extra=restrict)
        rr, lo, hi = cluster_bootstrap_rr(ns, ds, nr, dr, B=BOOT)
        df, dlo, dhi = cluster_bootstrap_diff(ns, ds, nru, dru, B=BOOT)
        pc = tot(nru * 0 + dr, np.full(len(scms), 4.0))
        blk = dict(label=label,
                   S_conditional=tot(ns, ds), R_conditional=tot(nr, dr),
                   R_unconditional=tot(nru, dru),
                   P_consistent_R=pc,
                   RR_conditional=dict(point=rr, ci=[lo, hi], B=BOOT),
                   diff_unconditional=dict(point=df, ci=[dlo, dhi], B=BOOT))
        # M1.2: decomposition of the unconditional difference
        p_cons_r = dr.sum() / (4.0 * len(scms))
        p_sil_s = ns.sum() / ds.sum() if ds.sum() else float("nan")
        blk["diff_decomposition"] = dict(
            free_component=float((1 - p_cons_r) * p_sil_s),
            class_component=float(df - (1 - p_cons_r) * p_sil_s),
            note="free = (1 - P(consistent|R)) x P(silent|S2); "
                 "residual = the class effect (AUDIT M1.2)")
        if rr >= 1.25 and lo > 1:
            v = "SUPPORTED"
        elif rr <= 0.80 and hi < 1:
            v = "REFUTED"
        else:
            v = "INDISTINGUISHABLE"
        blk["verdict_M1_riskratio"] = v
        return blk

    reach = lambda s, m: m["hopC"] < UNREACH
    hop0 = lambda s, m: m["hopC"] == 0
    valid_only = lambda s, m: mpdag_valid(m)
    reach_valid = lambda s, m: m["hopC"] < UNREACH and mpdag_valid(m)

    p2 = {}
    p2["PRIMARY_Suni_vs_R_reachable"] = damage_block(
        "Suni", "R", reach, "M11: reachable denominator (primary)")
    p2["Suni_vs_R_all_incl_unreachable"] = damage_block(
        "Suni", "R", None, "diluted denominator (secondary, M11)")
    p2["Suni_vs_R_hop0_only"] = damage_block("Suni", "R", hop0, "M9.3 hop-0 only")
    p2["Suni_vs_R_mpdag_valid_only"] = damage_block(
        "Suni", "R", reach_valid, "M5: restricted to mpdag_valid=True")
    p2["Sloc_vs_R_reachable"] = damage_block("Sloc", "R", reach, "S-loc, reachable")
    # ADDITION (not in PREREG/AUDIT): the withdrawal-controlled contrast
    nsx, dsx = per_scm(scms, "Suni", silent, lambda m: True, rho=1, extra=reach)
    ndx, ddx = per_scm(scms, "Drop", silent, lambda m: True, rho=1, extra=reach)
    rrd, lod, hid = cluster_bootstrap_rr(nsx, dsx, ndx, ddx, B=BOOT)
    p2["ADDED_Suni_vs_DROP_reachable"] = dict(
        S=tot(nsx, dsx), DROP=tot(ndx, ddx),
        RR=dict(point=rrd, ci=[lod, hid], B=BOOT),
        verdict=("SUPPORTED" if rrd >= 1.25 and lod > 1 else
                 "REFUTED" if rrd <= 0.80 and hid < 1 else "INDISTINGUISHABLE"),
        note="NOT PRE-REGISTERED. Arm D withdraws the same true statement without "
             "asserting anything, so this contrast isolates the spurious EDGE from "
             "the WITHDRAWAL that the replacement law confounds into arm S.")
    nrx, drx = per_scm(scms, "R", silent, lambda m: m["consistent"], rho=1, extra=reach)
    rrdr, lodr, hidr = cluster_bootstrap_rr(ndx, ddx, nrx, drx, B=BOOT)
    p2["ADDED_DROP_vs_R_reachable"] = dict(
        DROP=tot(ndx, ddx), R_conditional=tot(nrx, drx),
        RR=dict(point=rrdr, ci=[lodr, hidr], B=BOOT),
        note="the withdrawal effect on its own, against the reversal baseline")
    # M9.3 direct standardisation onto arm R's hop distribution
    def strata(arm, pred_den):
        d = defaultdict(lambda: [0, 0])
        for s in scms:
            for m in s[arm]:
                if m["rho"] != 1 or not pred_den(m):
                    continue
                k = hoplab(m["hopC"])
                d[k][1] += 1
                if silent(m):
                    d[k][0] += 1
        return {k: tuple(v) for k, v in d.items()}
    st_s = strata("Suni", lambda m: True)
    st_r = strata("R", lambda m: m["consistent"])
    p2["M9_3_standardised_Suni_onto_R_hopdist"] = standardise(st_s, st_r)
    p2["hop_strata_Suni_rho1"] = {k: rate(v[0], v[1]) for k, v in sorted(st_s.items())}
    p2["hop_strata_R_rho1_conditional"] = {k: rate(v[0], v[1])
                                           for k, v in sorted(st_r.items())}
    p2["free_bar_disclosure"] = dict(
        break_even_SUPPORTED_conditional_ratio_original=0.845,
        break_even_REFUTED_conditional_ratio_original=0.540,
        note="AUDIT A1/M1.3. Under the PREREG's unconditional 0.03 bar on "
             "`original`, an IDENTICALLY damaging spurious class prints SUPPORTED; "
             "the bar only stops saying SUPPORTED once the class is >15.5% "
             "relatively LESS damaging, and only says REFUTED at 46% less. "
             "The primary contrast is therefore the conditional risk ratio (M1).")
    # the PREREG's own (demoted) absolute rule, reported so it cannot be hidden
    du = p2["PRIMARY_Suni_vs_R_reachable"]
    dd = du["diff_unconditional"]["point"]
    ci_s = du["S_conditional"]["ci"]; ci_r = du["R_unconditional"]["ci"]
    p2["PREREG_4_2_absolute_rule_DEMOTED"] = dict(
        diff=dd, disjoint=(ci_s[0] > ci_r[1] or ci_r[0] > ci_s[1]),
        verdict=("SUPPORTED" if dd >= 0.03 and (ci_s[0] > ci_r[1]) else
                 "REFUTED" if dd <= 0 and (ci_r[0] > ci_s[1]) else
                 "INDISTINGUISHABLE"),
        note="M1.2/M3: kept only on `original`, demoted to secondary, printed with "
             "the free-bar disclosure above.")
    out["P2_damage"] = p2

    # ============================ X1-P3 A26's revive_if
    def sloc_row(pred):
        num = den = 0
        for s in scms:
            for m in s["Sloc"]:
                sel = [t for t in m["flip"]]
                if not all(pred(s, t) for t in sel):
                    continue
                den += 1
                if silent(m):
                    num += 1
        return rate(num, den)
    pure_nonq = lambda s, t: s["lab_loc"][t][0] == "pure_spurious" and not s["lab_loc"][t][1]
    pure_q = lambda s, t: s["lab_loc"][t][0] == "pure_spurious" and s["lab_loc"][t][1]
    rev = lambda s, t: s["lab_loc"][t][0] == "reversal_of_true_edge"

    def by_rho(pred, rhos):
        num = den = 0
        for s in scms:
            for m in s["Sloc"]:
                if m["rho"] not in rhos:
                    continue
                if not all(pred(s, t) for t in m["flip"]):
                    continue
                den += 1
                if silent(m):
                    num += 1
        return rate(num, den)

    p3 = {}
    p3["composition_Sloc_pool"] = dict(
        pure_spurious=rate(det_loc, 4 * n),
        reversal_of_true_edge=rate(sum(1 for s in scms for t in range(4)
                                       if s["lab_loc"][t][0] == "reversal_of_true_edge"),
                                   4 * n),
        query_pair=rate(sum(1 for s in scms for t in range(4) if s["lab_loc"][t][1]),
                        4 * n),
        scms_with_ge1_query_pair=rate(
            sum(1 for s in scms if any(s["lab_loc"][t][1] for t in range(4))), n),
        note="AUDIT M8/M14ii: ORACLE labels. The analyst cannot perform this "
             "partition, so no stratified rate here is a triage rule.")
    for lbl, pred in (("pure_spurious_NON_QUERY_PRIMARY", pure_nonq),
                      ("pure_spurious_QUERY_PAIR", pure_q),
                      ("reversal_of_true_edge", rev),
                      ("POOLED_all", lambda s, t: True)):
        p3[lbl] = dict(member_rho1=by_rho(pred, {1}),
                       member_rho_le4=by_rho(pred, {1, 2, 3, 4}))
    # M2.1 the IN-RUN reversal baseline, per member, rho<=4, same SCMs
    p3["IN_RUN_reversal_baseline_member_rho_le4"] = dict(
        unconditional=tot(*per_scm(scms, "R", silent, lambda m: True, rho=None)),
        conditional=tot(*per_scm(scms, "R", silent, lambda m: m["consistent"],
                                 rho=None)),
        note="AUDIT M2.1: the ONLY like-for-like comparator -- arm R on X1's own "
             "SCMs, per member, rho<=4, same grid, same |K|, same draw.")
    p3["legacy_baseline_quoted_once"] = LEGACY
    p3["A26_revive_if_window"] = [0.087, 0.100]
    prim = p3["pure_spurious_NON_QUERY_PRIMARY"]["member_rho_le4"]
    p3["A26_verdict"] = dict(
        primary_rate=prim,
        overlaps_window=bool(prim["ci"][0] <= 0.100 and prim["ci"][1] >= 0.087),
        note="AUDIT M2.4: under a never-raising operator the revive_if as written "
             "cannot fire unless the class is >= 44% MILDER than reversals "
             "(matched level 0.0933/(1-0.475) = 0.1776). Reporting a retraction "
             "off this rule would register a FALSE retraction. Evaluated for the "
             "record only; the decision rule that carries information is the "
             "conditional RR of P2.")
    out["P3_A26"] = p3

    # ============================ X1-P4 locality
    def cross(arm, use_drop=False):
        tab = defaultdict(lambda: [0, 0])
        for s in scms:
            dmap = {(m["rho"], tuple(m["flip"])): m for m in s["Drop"]} if use_drop else None
            for m in s[arm]:
                if m["rho"] != 1:
                    continue
                k = (hoplab(m["hopC"]), hoplab(m.get("hopG", m["hopC"])))
                tab[k][1] += 1
                ev = silent(m)
                if use_drop:
                    d = dmap[(m["rho"], tuple(m["flip"]))]
                    ev = (m.get("O") != d.get("O")) or (m.get("amenable") != d.get("amenable"))
                if ev:
                    tab[k][0] += 1
        return {f"hopC={a}|hopG={b}": rate(v[0], v[1]) for (a, b), v in sorted(tab.items())}

    p4 = {}
    for arm in ("Suni", "Sloc", "R"):
        p4[f"{arm}_hopC_x_hopG_silent_rho1"] = cross(arm)
    p4["Suni_hopC_x_hopG_EDGE_ATTRIBUTABLE_Ostar_move_rho1"] = cross("Suni", use_drop=True)
    ns0, ds0 = per_scm(scms, "Suni", silent, lambda m: True, rho=1, extra=hop0)
    ns1, ds1 = per_scm(scms, "Suni", silent, lambda m: True, rho=1,
                       extra=lambda s, m: 1 <= m["hopC"] < UNREACH)
    r0 = tot(ns0, ds0); r1 = tot(ns1, ds1)
    mat = (r1["rate"] / r0["rate"]) if (r0["rate"] and r1["rate"] is not None) else None
    p4["Suni_hop0"] = r0
    p4["Suni_hop_ge1_reachable"] = r1
    p4["Suni_hop_unreachable"] = tot(*per_scm(scms, "Suni", silent, lambda m: True,
                                              rho=1,
                                              extra=lambda s, m: m["hopC"] >= UNREACH))
    p4["materiality_ratio"] = mat
    p4["n_hop_ge1_denominator"] = int(ds1.sum())
    p4["verdict_M12"] = dict(
        readable=bool(ds1.sum() >= 2000),
        existence_check_ge5_events=bool(ns1.sum() >= 5),
        n_events=int(ns1.sum()),
        verdict=("NOT READABLE (n < 2000)" if ds1.sum() < 2000 else
                 "SUPPORTED" if (mat is not None and mat >= 0.10 and r1["ci"][0] > 0)
                 else "REFUTED" if ns1.sum() == 0 else "NO VERDICT (rate reported)"),
        note="AUDIT M12: rate + Wilson + materiality ratio (>= 0.10 x the hop-0 rate "
             "with CI_lo > 0). The '>= 5 events' count is an EXISTENCE check only.")
    # drop-controlled locality: the added EDGE alone
    ea1 = ea1d = 0
    for s in scms:
        dmap = {(m["rho"], tuple(m["flip"])): m for m in s["Drop"]}
        for m in s["Suni"]:
            if m["rho"] != 1 or not (1 <= m["hopC"] < UNREACH):
                continue
            d = dmap[(m["rho"], tuple(m["flip"]))]
            ea1d += 1
            if (m.get("O") != d.get("O")) or (m.get("amenable") != d.get("amenable")):
                ea1 += 1
    p4["ADDED_edge_attributable_hop_ge1"] = rate(ea1, ea1d)
    out["P4_locality"] = p4

    # ============================ X1-P5 severity
    def relb(arm):
        v = []
        for s in scms:
            for m in s[arm]:
                if m["rho"] == 1 and silent(m) and abs(s["tau"]) >= 1e-3:
                    v.append(abs(m["est"] - s["tau"]) / abs(s["tau"]))
        return np.array(v)
    a = relb("Suni"); b = relb("R")
    cd = cliffs_delta(a, b, B=BOOT)
    out["P5_severity"] = dict(
        cliffs_delta_Suni_vs_R=cd,
        dominance=("DECLARED" if cd["delta"] is not None and cd["delta"] >= 0.15
                   and cd["ci"][0] > 0 else "NOT SEPARATED"),
        ecdf_quantiles_Suni={str(q): float(np.quantile(a, q)) for q in
                             (.1, .25, .5, .75, .9)} if a.size else None,
        ecdf_quantiles_R={str(q): float(np.quantile(b, q)) for q in
                          (.1, .25, .5, .75, .9)} if b.size else None,
        note="the ONLY inferential statistic in the run, and it is on a bias "
             "MAGNITUDE distribution. No p-value on rho*, ever.")

    # ============================ X1-P6 the movement channel
    def rho_star(s, arm, nn, z):
        se = 0.0 if np.isinf(nn) else (
            s["se_factor0"] / np.sqrt(nn - s["k_reg0"] - 1)
            if nn > s["k_reg0"] + 1 else np.inf)
        mr = s["max_rho"]
        d = np.zeros(mr + 1)
        for r in range(1, mr + 1):
            best = d[r - 1]
            for m in s[arm]:
                if m["rho"] == r and m.get("consistent") and m.get("amenable"):
                    best = max(best, abs(m["est"] - s["est0"]))
            d[r] = best
        rse = mr + 1
        for r in range(1, mr + 1):
            if d[r] > z * se:
                rse = r
                break
        rid = mr + 1
        for r in range(1, mr + 1):
            if any(m["rho"] == r and m.get("consistent") and not m.get("amenable")
                   for m in s[arm]):
                rid = r
                break
        return rse, rid
    p6 = {}
    for arm in ("R", "Suni", "Sloc", "Drop"):
        blk = {}
        for nn in N_GRID:
            for z in ([Z95] if nn not in (20000, float("inf")) else Z_GRID):
                key = f"n={nn}|z={z:.3f}"
                rs = [rho_star(s, arm, nn, z) for s in scms]
                mr = scms[0]["max_rho"]
                blk[key] = dict(
                    censored_rho_star_se=rate(sum(1 for a_, _ in rs if a_ == mr + 1), n),
                    censored_rho_star_ident=rate(sum(1 for _, b_ in rs if b_ == mr + 1), n),
                    censored_rho_star_any=rate(
                        sum(1 for a_, b_ in rs if min(a_, b_) == mr + 1), n))
        p6[arm] = blk
    out["P6_movement"] = dict(
        channels=p6,
        note="rho*_se, rho*_ident and rho*_any reported as THREE separate columns. "
             "Folding identification failure into the movement statistic is what "
             "manufactured E1's false contradiction with S5 (PREREG 11.3).")

    with open(f"../results/x1_analysis_{ens}.json", "w") as f:
        json.dump(out, f, indent=1, default=float)
    print(f"[x1-analysis/{ens}] n={n} -> ../results/x1_analysis_{ens}.json", flush=True)
    gv = {k: v.get("verdict") for k, v in out["gates"].items() if isinstance(v, dict)}
    print("  gates:", json.dumps(gv), flush=True)


if __name__ == "__main__":
    main()
