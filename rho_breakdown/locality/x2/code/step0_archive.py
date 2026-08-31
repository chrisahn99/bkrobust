"""
X2 STEP 0 -- dissolve the E1' / S5 contradiction in code, and recover ARM D.

Reads ONLY ~/latent-causal/e1prime-se/results/arm1_*.json (read-only; nothing
under e1prime-se/ is written).  Produces:

  (1) E1's own per-SCM locality table reproduced from raw, then EVERY
      non-censored dmin>=1 SCM decomposed into {ident-only, rho*_se=1,
      rho*_se>=2} at n=inf and n=20000.  This is the caller's step 1:
      hypotheses (a) genuine O*-change / (b) min-vs-max endpoint artefact /
      (c) rho*_any folding in the identification criterion / (d) a bug.
  (2) The per-STATEMENT rho=1 table with all five denominators (M8), by dmin
      bucket AND by dmax bucket (M10), plus `disconnected` as its own row (G7).
  (3) ARM D: Q_r and Cost_r for r in {-1,0,1,2,3} with Clopper-Pearson
      intervals (M19), mean Cost + per-SCM fraction below 0.5 (M5).
  (4) The S3/S4/S5 contingency table and its decomposition (M12).
  (5) The within-SCM matched near/far contrast (M13).
  (6) The competing-risk table D_am/D_con by stratum and the composite change
      rate on D_con (M14).
  (7) P(cascade > 0) per stratum (M18) and the E_bias == E_set identity (P9).
  (8) rho>=2: "all flips far" and ">=1 flip far" strata (M18 #3).

Usage: python step0_archive.py <arch_dir> <out.json>
"""
import json
import sys
from collections import defaultdict
from itertools import combinations

import numpy as np

import x2lib as X
from se import deviation_radius, ident_radius, rho_star_se, Z95

ARCH = sys.argv[1] if len(sys.argv) > 1 else "/home/costaj/latent-causal/e1prime-se/results"
OUT = sys.argv[2] if len(sys.argv) > 2 else "../results/step0_archive.json"

CELLS = [("original", "K4", 4), ("licensed", "K4", 4), ("large", "K4", 4),
         ("k8", "K4", 3), ("k8", "K6", 3), ("k8", "K8", 3)]
TAU_FLOOR = 1e-3          # E1's analysis-time filter, reproduced verbatim
N_GRID = [("inf", np.inf), ("20000", 20000.0)]


def flips_of(k, max_rho):
    out = []
    for rho in range(1, min(max_rho, k) + 1):
        for fl in combinations(range(k), rho):
            out.append(fl)
    return out


def se_at(arm, n):
    if not np.isfinite(n):
        return 0.0
    df = n - arm["k_reg"] - 1
    return arm["se_factor"] / np.sqrt(df) if df > 0 else np.inf


def load(ens):
    return json.load(open(f"{ARCH}/arm1_{ens}.json"))


def main():
    out = {}
    for ens, key, max_rho in CELLS:
        recs = load(ens)
        cell = dict(ens=ens, arm=key, max_rho=max_rho)

        # ---------------------------------------------------- per-SCM (E1's unit)
        per_scm = []
        for rec in recs:
            if abs(rec["tau"]) < TAU_FLOOR:
                continue
            arm = rec["arms"].get(key)
            if arm is None or not arm.get("mpdag_amenable"):
                continue
            k = arm["n_K"]
            fl = flips_of(k, max_rho)
            ms = arm["members"]
            assert len(ms) == len(fl)
            hop = [float(h) if h < X.SENTINEL else np.inf for h in arm["stmt_hopdist"]]
            d = deviation_radius(ms, arm["est0"], max_rho)
            ident = ident_radius(ms, max_rho)
            row = dict(seed=rec["seed"], p=rec["p"], deg=rec["deg"], k=k,
                       hopmin=min(hop), hop=hop, ident=int(ident))
            for lab, n in N_GRID:
                rse = rho_star_se(d, se_at(arm, n), Z95, max_rho)
                row[f"rse_{lab}"] = int(rse)
                row[f"any_{lab}"] = int(min(rse, ident))
                row[f"cens_{lab}"] = bool(min(rse, ident) == max_rho + 1)
            per_scm.append(row)
        cell["n_scm_analysed"] = len(per_scm)

        # E1's locality table, reproduced (per-SCM HOPMIN strata)
        loc = {}
        for lab, _ in N_GRID:
            byh = defaultdict(list)
            for r in per_scm:
                byh[X.dbucket(r["hopmin"])].append(r)
            loc[lab] = {h: dict(n=len(v),
                                censored=float(np.mean([r[f"cens_{lab}"] for r in v])))
                        for h, v in sorted(byh.items())}
        cell["E1_locality_by_min_hopdist_reproduced"] = loc

        # the decomposition the caller asked for: non-censored, hopmin >= 1
        dec = {}
        for lab, _ in N_GRID:
            sel = [r for r in per_scm if r["hopmin"] >= 1 and not r[f"cens_{lab}"]]
            tot = [r for r in per_scm if r["hopmin"] >= 1]
            dec[lab] = dict(
                n_hopmin_ge1=len(tot),
                n_not_censored=len(sel),
                ident_only=sum(1 for r in sel if r["ident"] <= max_rho and r[f"rse_{lab}"] > max_rho),
                rse_eq_1=sum(1 for r in sel if r[f"rse_{lab}"] == 1),
                rse_ge_2=sum(1 for r in sel if 2 <= r[f"rse_{lab}"] <= max_rho),
                both=sum(1 for r in sel if r["ident"] <= max_rho and r[f"rse_{lab}"] <= max_rho),
                disconnected_only=sum(1 for r in sel if not np.isfinite(r["hopmin"])),
            )
        cell["noncensored_hopmin_ge1_decomposition"] = dec

        # ---------------------------------------------------- per-STATEMENT rho=1
        st = defaultdict(lambda: defaultdict(int))       # bucket -> counters
        stx = defaultdict(lambda: defaultdict(int))      # dmax bucket -> counters
        matched = []                                     # M13
        rstar = dict(set=-1.0, ident=-1.0)
        ebias_ne_eset = 0
        est_ne_set = 0
        for rec in recs:
            if abs(rec["tau"]) < TAU_FLOOR:
                continue
            arm = rec["arms"].get(key)
            if arm is None or not arm.get("mpdag_amenable"):
                continue
            k = arm["n_K"]
            fl = flips_of(k, max_rho)
            ms = arm["members"]
            hop = [float(h) if h < X.SENTINEL else np.inf for h in arm["stmt_hopdist"]]
            Kst = arm["K"]
            x, y = rec["x"], rec["y"]
            # dmax per statement: an endpoint is at distance 0 iff it is x or y,
            # so dmax = dmin unless exactly one endpoint is the query.
            dmaxs = []
            for i, (u, v) in enumerate(Kst):
                inc_u = u in (x, y)
                inc_v = v in (x, y)
                if inc_u and inc_v:
                    dmaxs.append(0.0)
                elif inc_u or inc_v:
                    dmaxs.append(max(hop[i], 1.0) if np.isfinite(hop[i]) else np.inf)
                else:
                    dmaxs.append(hop[i])       # >=1 anyway; exact value not stored
            loc_near = dict(n=0, N=0)
            loc_far = dict(n=0, N=0)
            for i, (flip, m) in enumerate(zip(fl, ms)):
                if len(flip) != 1:
                    continue
                j = flip[0]
                b = X.dbucket(hop[j])
                bx = X.dbucket(dmaxs[j])
                cons = bool(m.get("consistent"))
                amen = bool(m.get("amenable"))
                st[b]["D_all"] += 1
                stx[bx]["D_all"] += 1
                if cons:
                    st[b]["D_con"] += 1
                    stx[bx]["D_con"] += 1
                    if m.get("cascade", 0) and m["cascade"] > 0:
                        st[b]["cascade_pos"] += 1
                    if not amen:
                        st[b]["N_ident"] += 1
                        stx[bx]["N_ident"] += 1
                        rstar["ident"] = max(rstar["ident"], hop[j] if np.isfinite(hop[j]) else -1.0)
                    else:
                        st[b]["D_am"] += 1
                        stx[bx]["D_am"] += 1
                        if np.isfinite(hop[j]):
                            st[b]["D_s5"] += 1     # ARM B: the flipped stmt is always false
                        chg = bool(m.get("ostar_changed"))
                        moved = abs(m["est"] - arm["est0"]) > X.EPS_BETA
                        if chg:
                            st[b]["N_set"] += 1
                            stx[bx]["N_set"] += 1
                            rstar["set"] = max(rstar["set"], hop[j] if np.isfinite(hop[j]) else -1.0)
                            if m.get("ostar_valid"):
                                ebias_ne_eset += 1      # changed but STILL VALID
                            else:
                                st[b]["N_bias"] += 1
                                stx[bx]["N_bias"] += 1
                        if moved:
                            st[b]["N_est_inf"] += 1
                        if moved != chg:
                            est_ne_set += 1
                        if hop[j] == 0:
                            loc_near["n"] += 1
                            loc_near["N"] += int(chg)
                        elif np.isfinite(hop[j]):
                            loc_far["n"] += 1
                            loc_far["N"] += int(chg)
            if loc_near["n"] and loc_far["n"]:
                matched.append((loc_near["n"], loc_near["N"], loc_far["n"], loc_far["N"]))
        cell["per_statement_rho1_by_dmin"] = {b: dict(v) for b, v in sorted(st.items())}
        cell["per_statement_rho1_by_dmax"] = {b: dict(v) for b, v in sorted(stx.items())}
        cell["r_star"] = rstar
        cell["E_bias_ne_E_set"] = ebias_ne_eset
        cell["E_est_inf_ne_E_set"] = est_ne_set
        if matched:
            a = np.array(matched, dtype=float)
            cell["matched_within_scm"] = dict(
                n_scm=len(matched),
                D_am_near=int(a[:, 0].sum()), N_set_near=int(a[:, 1].sum()),
                D_am_far=int(a[:, 2].sum()), N_set_far=int(a[:, 3].sum()))
        else:
            cell["matched_within_scm"] = None

        # ---------------------------------------------------- rho >= 2 strata (M18 #3)
        r2 = defaultdict(lambda: defaultdict(int))
        for rec in recs:
            if abs(rec["tau"]) < TAU_FLOOR:
                continue
            arm = rec["arms"].get(key)
            if arm is None or not arm.get("mpdag_amenable"):
                continue
            k = arm["n_K"]
            hop = [float(h) if h < X.SENTINEL else np.inf for h in arm["stmt_hopdist"]]
            for flip, m in zip(flips_of(k, max_rho), arm["members"]):
                if len(flip) < 2:
                    continue
                hs = [hop[j] for j in flip]
                allfar = all(np.isfinite(h) and h >= 1 for h in hs)
                anyfar = any(np.isfinite(h) and h >= 1 for h in hs)
                for lab, ok in (("all_far", allfar), ("any_far", anyfar), ("all", True)):
                    if not ok:
                        continue
                    kk = f"rho{len(flip)}_{lab}"
                    r2[kk]["D_all"] += 1
                    if m.get("consistent"):
                        r2[kk]["D_con"] += 1
                        if m.get("amenable"):
                            r2[kk]["D_am"] += 1
                            if m.get("ostar_changed"):
                                r2[kk]["N_set"] += 1
                                if not m.get("ostar_valid"):
                                    r2[kk]["N_bias"] += 1
                        else:
                            r2[kk]["N_ident"] += 1
        cell["rho_ge2"] = {k2: dict(v) for k2, v in sorted(r2.items())}

        # ---------------------------------------------------- ARM D
        armd = {}
        for r in (-1, 0, 1, 2, 3):
            nq = nchg = 0
            costs = []
            for rec in recs:
                if abs(rec["tau"]) < TAU_FLOOR:
                    continue
                arm = rec["arms"].get(key)
                if arm is None or not arm.get("mpdag_amenable"):
                    continue
                hop = [float(h) if h < X.SENTINEL else np.inf for h in arm["stmt_hopdist"]]
                vf, vp, nf, npr = X.prune_partition(arm["K"], hop, arm["members"], max_rho, r)
                If = X.interval(vf, arm["est0"])
                Ip = X.interval(vp, arm["est0"])
                nq += 1
                if abs(If[0] - Ip[0]) > 1e-9 or abs(If[1] - Ip[1]) > 1e-9:
                    nchg += 1
                costs.append(npr / nf)
            costs = np.array(costs)
            lo, hi = X.clopper_pearson(nchg, nq)
            armd[str(r)] = dict(n_scm=nq, n_changed=nchg, Q=nchg / max(nq, 1),
                                Q_ci=[lo, hi], cost_mean=float(costs.mean()),
                                cost_median=float(np.median(costs)),
                                cost_frac_below_half=float((costs <= 0.5).mean()))
        cell["armD"] = armd

        # ---------------------------------------------------- G10 placebo
        pl = {}
        for r in (-1, 0, 1, 2, 3):
            nchg = nq = 0
            for rec in recs:
                if abs(rec["tau"]) < TAU_FLOOR:
                    continue
                arm = rec["arms"].get(key)
                if arm is None or not arm.get("mpdag_amenable"):
                    continue
                hop = [float(h) if h < X.SENTINEL else np.inf for h in arm["stmt_hopdist"]]
                vf, vp, _, _ = X.prune_partition(arm["K"], hop, arm["members_frozen"], max_rho, r)
                If = X.interval(vf, arm["est0"])
                Ip = X.interval(vp, arm["est0"])
                nq += 1
                if abs(If[0] - Ip[0]) > 1e-9 or abs(If[1] - Ip[1]) > 1e-9:
                    nchg += 1
            pl[str(r)] = dict(n=nq, Q=nchg / max(nq, 1))
        cell["armD_placebo_G10"] = pl

        out[f"{ens}/{key}"] = cell
        print(f"[step0] {ens}/{key}: {cell['n_scm_analysed']} SCMs, "
              f"per-stmt rho=1 D_all={sum(v.get('D_all',0) for v in st.values())}", flush=True)

    with open(OUT, "w") as f:
        json.dump(out, f, indent=1, default=float)
    print(f"[step0] -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
