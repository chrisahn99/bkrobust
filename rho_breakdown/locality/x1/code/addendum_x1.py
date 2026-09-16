"""
X1 ADDENDUM -- items the main analysis' own numbers made necessary.

1. M13 is now BINDING: design effects on S-loc are 1.64/2.10/1.76 (> 1.25), so
   every headline rate gets an SCM-level cluster-bootstrap CI beside its Wilson.
2. The DROP arm's ostar_changed rate -- quantifies the withdrawal confound on
   outcome (b) and shows it is null on outcome (c).
3. The PREREG's ORIGINAL (contaminated) S1 predicate, conflict OR cycle OR
   vstruct_vs_C, reported for the record beside M4's intrinsic one.
4. est0 == tau whenever O0 is valid -- the identity that explains why d(rho) > 0
   and `silent` are the same event, and why the DROP arm censors at 1.000.
"""
import json
import sys

import numpy as np

from stats_x1 import rate, wilson
from x1_ops import UNREACH


def boot_rate(num, den, B=10000, seed=20260819):
    num = np.asarray(num, float); den = np.asarray(den, float)
    m = len(num)
    pt = num.sum() / den.sum() if den.sum() else float("nan")
    rng = np.random.default_rng(seed)
    out = []
    done = 0
    while done < B:
        b = min(500, B - done)
        idx = rng.integers(0, m, size=(b, m))
        with np.errstate(divide="ignore", invalid="ignore"):
            out.append(num[idx].sum(1) / den[idx].sum(1))
        done += b
    o = np.concatenate(out); o = o[np.isfinite(o)]
    lo, hi = np.quantile(o, [0.025, 0.975])
    return dict(rate=float(pt), k=int(num.sum()), n=int(den.sum()),
                ci_wilson=list(wilson(int(num.sum()), int(den.sum()))),
                ci_cluster=[float(lo), float(hi)])


def silent(m):
    return bool(m.get("consistent") and m.get("amenable")
                and m.get("ostar_valid") is False)


out = {}
for e in ("original", "licensed", "large"):
    raw = json.load(open(f"../results/x1_{e}.json"))["scms"]
    n = len(raw)
    blk = {}
    # --- 1. cluster-bootstrap CIs on the headline rates ----------------------
    def pv(arm, den_pred, num_pred=silent, rho=1, reach=True):
        a = np.zeros(n); b = np.zeros(n)
        for i, s in enumerate(raw):
            for m in s[arm]:
                if rho is not None and m["rho"] != rho:
                    continue
                if reach and m["hopC"] >= UNREACH:
                    continue
                if not den_pred(m):
                    continue
                b[i] += 1
                if num_pred(m):
                    a[i] += 1
        return a, b
    blk["M13_cluster_CIs_rho1_reachable"] = dict(
        Suni_S2=boot_rate(*pv("Suni", lambda m: True)),
        Sloc_S2=boot_rate(*pv("Sloc", lambda m: True)),
        R_conditional=boot_rate(*pv("R", lambda m: m["consistent"])),
        R_unconditional=boot_rate(*pv("R", lambda m: True)),
        Suni_S1_conditional=boot_rate(*pv("Suni", lambda m: m["consistent_S1"])),
    )
    blk["M13_cluster_CIs_Sloc_rho_le4_member"] = boot_rate(
        *pv("Sloc", lambda m: True, rho=None, reach=False))
    # --- 2. the DROP arm: outcome (b) vs outcome (c) -------------------------
    oc = si = am = tot_ = 0
    for s in raw:
        for m in s["Drop"]:
            if m["rho"] != 1:
                continue
            tot_ += 1
            am += int(bool(m.get("amenable")))
            oc += int(bool(m.get("ostar_changed", False)))
            si += int(silent(m))
    blk["DROP_arm_rho1"] = dict(
        ostar_changed=rate(oc, tot_), silent=rate(si, tot_),
        amenable=rate(am, tot_),
        note="the WITHDRAWAL confound, quantified: it moves O* (outcome b) but "
             "NEVER produces silent bias (outcome c). So the P2 contrast on "
             "`silent` is clean; a P2 contrast on `ostar_changed` would not be.")
    # matched: does the added edge ever move O* beyond what the drop already did?
    ea = ead = 0
    for s in raw:
        dm = {(m["rho"], tuple(m["flip"])): m for m in s["Drop"]}
        for m in s["Suni"]:
            if m["rho"] != 1:
                continue
            d = dm[(m["rho"], tuple(m["flip"]))]
            ead += 1
            if m.get("O") != d.get("O") or m.get("amenable") != d.get("amenable"):
                ea += 1
    blk["edge_attributable_Ostar_move_rho1_all_hops"] = rate(ea, ead)
    # --- 3. the PREREG's original contaminated S1 predicate ------------------
    for arm in ("Suni", "Sloc"):
        c = t = 0
        for s in raw:
            for m in s[arm]:
                if m["rho"] != 1:
                    continue
                t += 1
                if m["conflict"] or m["cycle"] or m["vsC"]:
                    c += 1
        blk[f"{arm}_S1_PREREG_contaminated_predicate_rho1"] = dict(
            **rate(c, t),
            note="conflict OR cycle OR (v_structures(G) != v_structures(C)). "
                 "AUDIT M4 replaced this with the intrinsic Dor-Tarsi predicate "
                 "because it fires for a purely skeletal reason. Reported for "
                 "the record only.")
    # --- 4. est0 == tau whenever O0 is valid ---------------------------------
    bad = 0; chk = 0; mx = 0.0
    for s in raw:
        if not s["O0_valid"]:
            continue
        chk += 1
        d = abs(s["est0"] - s["tau"])
        mx = max(mx, d)
        if d > 1e-9:
            bad += 1
    blk["est0_equals_tau_when_O0_valid"] = dict(
        n=chk, violations=bad, max_abs_diff=float(mx),
        n_O0_invalid=sum(1 for s in raw if not s["O0_valid"]),
        note="ALL valid adjustment sets give the same POPULATION estimate (= tau). "
             "Hence |est_m - est0| > 0  <=>  O_m invalid  <=>  `silent`, which is "
             "why E1's changed_still_valid was 0 everywhere and why the DROP arm "
             "censors rho*_se at exactly 1.000.")
    out[e] = blk

json.dump(out, open("../results/x1_addendum.json", "w"), indent=1, default=float)
for e, blk in out.items():
    print("=" * 96)
    print("ENSEMBLE", e)
    for k, v in blk["M13_cluster_CIs_rho1_reachable"].items():
        print("  %-22s %.4f  wilson[%.4f,%.4f]  CLUSTER[%.4f,%.4f]  k=%d n=%d"
              % (k, v["rate"], v["ci_wilson"][0], v["ci_wilson"][1],
                 v["ci_cluster"][0], v["ci_cluster"][1], v["k"], v["n"]))
    v = blk["M13_cluster_CIs_Sloc_rho_le4_member"]
    print("  %-22s %.4f  wilson[%.4f,%.4f]  CLUSTER[%.4f,%.4f]  k=%d n=%d"
          % ("Sloc rho<=4 (member)", v["rate"], v["ci_wilson"][0], v["ci_wilson"][1],
             v["ci_cluster"][0], v["ci_cluster"][1], v["k"], v["n"]))
    d = blk["DROP_arm_rho1"]
    print("  DROP rho=1: ostar_changed k=%d n=%d rate=%.4f | silent k=%d rate=%.4f "
          "| amenable rate=%.4f" % (d["ostar_changed"]["k"], d["ostar_changed"]["n"],
                                    d["ostar_changed"]["rate"], d["silent"]["k"],
                                    d["silent"]["rate"], d["amenable"]["rate"]))
    ea = blk["edge_attributable_Ostar_move_rho1_all_hops"]
    print("  edge-attributable O* move (all hops, rho=1): %.4f [%.4f,%.4f] k=%d n=%d"
          % (ea["rate"], ea["ci"][0], ea["ci"][1], ea["k"], ea["n"]))
    for arm in ("Suni", "Sloc"):
        v = blk[f"{arm}_S1_PREREG_contaminated_predicate_rho1"]
        print("  %s S1 under the PREREG's contaminated predicate: %.4f [%.4f,%.4f] k=%d n=%d"
              % (arm, v["rate"], v["ci"][0], v["ci"][1], v["k"], v["n"]))
    v = blk["est0_equals_tau_when_O0_valid"]
    print("  est0 == tau when O0 valid: %d/%d violations, max|diff|=%.3e ; "
          "O0 invalid in %d SCMs" % (v["violations"], v["n"], v["max_abs_diff"],
                                     v["n_O0_invalid"]))
