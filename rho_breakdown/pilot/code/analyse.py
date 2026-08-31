"""
Turn linear_raw.json into the two headline numbers, the breakdown radius, and
the figures.

Usage: python analyse.py ../results/linear_raw.json ../results
"""
import json
import sys
from collections import defaultdict

import numpy as np

RNG = np.random.default_rng(7)


def boot_ci(v, f=np.mean, B=2000):
    v = np.asarray(v, dtype=float)
    v = v[~np.isnan(v)]
    if len(v) == 0:
        return (np.nan, np.nan, np.nan)
    stat = f(v)
    idx = RNG.integers(0, len(v), size=(B, len(v)))
    bs = f(v[idx], axis=1)
    return (float(stat), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)))


def fmt(t):
    return f"{t[0]:.3f} [{t[1]:.3f}, {t[2]:.3f}]"


# ------------------------------------------------------------------ per-SCM summaries
def interval_at(arm, rho, est0):
    """Ignorance interval over all consistent+amenable members at distance <= rho,
    including the rho=0 point (the asserted K itself)."""
    ests = [est0] + [m["est"] for m in arm["members"]
                     if m["rho"] <= rho and m.get("consistent") and m.get("amenable")]
    n_nonamen = sum(1 for m in arm["members"]
                    if m["rho"] <= rho and m.get("consistent") and not m.get("amenable"))
    return min(ests), max(ests), n_nonamen, len(ests)


def breakdown_radius(arm, tau, est0, max_rho, crit):
    """Smallest rho at which SOME consistent member overturns the conclusion.
    crit='sign'  : estimate crosses 0 / flips sign vs tau
    crit='rel50' : |est - tau| > 0.5*|tau|
    crit='ident' : some member is non-amenable (effect not identified at all)
    Returns max_rho+1 if never overturned inside the ball."""
    for rho in range(1, max_rho + 1):
        for m in arm["members"]:
            if m["rho"] != rho or not m.get("consistent"):
                continue
            if crit == "ident":
                if not m.get("amenable"):
                    return rho
                continue
            if not m.get("amenable"):
                continue
            e = m["est"]
            if crit == "sign" and np.sign(e) != np.sign(tau):
                return rho
            if crit == "rel50" and abs(e - tau) > 0.5 * abs(tau):
                return rho
    return max_rho + 1


def main():
    raw_path = sys.argv[1] if len(sys.argv) > 1 else "../results/linear_raw.json"
    outdir = sys.argv[2] if len(sys.argv) > 2 else "../results"
    res = json.load(open(raw_path))
    print(f"{len(res)} SCMs")

    S = {"n_scm": len(res), "arms": {}}

    # ---------------- sortability diagnostics
    S["sortability"] = {
        k: fmt(boot_ci([r["diag"][k] for r in res]))
        for k in ("varsort_iscm", "r2sort_iscm", "varsort_naive", "r2sort_naive")
    }

    S["cpdag_amenable_frac"] = fmt(boot_ci([float(r["cpdag_amenable"]) for r in res]))

    per_arm_plot = {}
    for arm_name in ("generic", "tiered"):
        rs = [r for r in res if arm_name in r["arms"]]
        if not rs:
            continue
        A = {"n_scm": len(rs)}
        A["n_K"] = fmt(boot_ci([r["arms"][arm_name]["n_K"] for r in rs]))

        # ---------- NUMBER 1: what happens to O* under a single orientation error
        buckets = defaultdict(int)
        per_scm_changed, per_scm_invalid = [], []
        expels = []
        cascades = []
        for r in rs:
            arm = r["arms"][arm_name]
            ch = iv = n = 0
            for m in arm["members"]:
                if m["rho"] != 1:
                    continue
                if not m["consistent"]:
                    buckets["meek_inconsistent"] += 1
                    continue
                expels.append(float(m["expels_true_dag"]))
                cascades.append(m["cascade"])
                if not m["amenable"]:
                    buckets["not_amenable"] += 1
                    n += 1
                    continue
                n += 1
                if not m["ostar_changed"]:
                    buckets["ostar_unchanged"] += 1
                elif m["ostar_valid"]:
                    buckets["changed_still_valid"] += 1
                    ch += 1
                else:
                    buckets["changed_invalid"] += 1
                    ch += 1
                    iv += 1
            if n:
                per_scm_changed.append(ch / n)
                per_scm_invalid.append(iv / n)
        tot = sum(buckets.values())
        A["rho1_members_total"] = tot
        A["rho1_buckets"] = {k: v for k, v in sorted(buckets.items())}
        A["rho1_buckets_frac"] = {k: round(v / tot, 4) for k, v in sorted(buckets.items())}
        cons = tot - buckets["meek_inconsistent"]
        A["rho1_frac_meek_inconsistent"] = round(buckets["meek_inconsistent"] / tot, 4)
        A["NUMBER1_frac_consistent_flips_changing_Ostar"] = round(
            (buckets["changed_still_valid"] + buckets["changed_invalid"] +
             buckets["not_amenable"]) / cons, 4)
        A["NUMBER1_frac_consistent_flips_making_Ostar_INVALID"] = round(
            buckets["changed_invalid"] / cons, 4)
        A["NUMBER1_per_scm_frac_changing_Ostar"] = fmt(boot_ci(per_scm_changed))
        A["NUMBER1_per_scm_frac_invalid_Ostar"] = fmt(boot_ci(per_scm_invalid))
        A["consistency_correctness_gap_frac_consistent_flips_expelling_true_DAG"] = fmt(
            boot_ci(expels))
        A["cascade_extra_orientations_rho1"] = fmt(boot_ci(cascades))
        A["cascade_frac_nonzero_rho1"] = fmt(boot_ci([float(c > 0) for c in cascades]))

        # ---------- NUMBER 2: ignorance interval
        A["ignorance"] = {}
        widths_by_rho = {}
        for rho in (1, 2, 3):
            w, wn, excl0, cov, nonamen_any = [], [], [], [], []
            for r in rs:
                arm = r["arms"][arm_name]
                if arm["n_K"] < rho:
                    continue
                lo, hi, n_na, k = interval_at(arm, rho, arm["est0"])
                tau = r["tau"]
                w.append(hi - lo)
                wn.append((hi - lo) / abs(tau))
                excl0.append(float(lo > 0 or hi < 0))
                cov.append(float(lo <= tau <= hi))
                nonamen_any.append(float(n_na > 0))
            A["ignorance"][f"rho{rho}"] = {
                "width_abs": fmt(boot_ci(w)),
                "width_over_|tau|_mean": fmt(boot_ci(wn)),
                "width_over_|tau|_median": fmt(boot_ci(wn, f=lambda a, axis=None: np.median(a, axis=axis))),
                "width_over_|tau|_q90": float(np.percentile(wn, 90)),
                "frac_interval_excludes_0": fmt(boot_ci(excl0)),
                "frac_interval_covers_tau": fmt(boot_ci(cov)),
                "frac_scm_with_a_nonidentified_member": fmt(boot_ci(nonamen_any)),
            }
            widths_by_rho[rho] = wn
        per_arm_plot[arm_name] = {"widths": widths_by_rho}

        # ---------- breakdown radius
        A["breakdown_radius"] = {}
        for crit in ("sign", "rel50", "ident"):
            rr = []
            for r in rs:
                arm = r["arms"][arm_name]
                mr = min(3, arm["n_K"])
                rr.append(breakdown_radius(arm, r["tau"], arm["est0"], mr, crit))
            rr = np.array(rr)
            A["breakdown_radius"][crit] = {
                "frac_rho_star_eq_1": fmt(boot_ci((rr == 1).astype(float))),
                "frac_rho_star_le_2": fmt(boot_ci((rr <= 2).astype(float))),
                "frac_rho_star_le_3": fmt(boot_ci((rr <= 3).astype(float))),
                "frac_never_overturned_in_ball": fmt(boot_ci((rr > 3).astype(float))),
                "median": float(np.median(rr)),
            }
            per_arm_plot[arm_name][f"rstar_{crit}"] = rr.tolist()

        # ---------- stratified by whether the CPDAG alone is already amenable
        A["strata"] = {}
        for lab, sel in (("cpdag_amenable", True), ("cpdag_NOT_amenable", False)):
            sub = [r for r in rs if r["cpdag_amenable"] == sel]
            if not sub:
                continue
            wn, inval = [], []
            for r in sub:
                arm = r["arms"][arm_name]
                lo, hi, _, _ = interval_at(arm, 1, arm["est0"])
                wn.append((hi - lo) / abs(r["tau"]))
                mm = [m for m in arm["members"]
                      if m["rho"] == 1 and m["consistent"] and m.get("amenable")]
                if mm:
                    inval.append(np.mean([float(not m["ostar_valid"]) for m in mm]))
            A["strata"][lab] = {
                "n": len(sub),
                "rho1_width_over_|tau|_median": float(np.median(wn)),
                "rho1_width_over_|tau|_mean": fmt(boot_ci(wn)),
                "rho1_frac_Ostar_invalid": fmt(boot_ci(inval)),
            }
        S["arms"][arm_name] = A

    with open(f"{outdir}/summary.json", "w") as f:
        json.dump(S, f, indent=1)
    print(json.dumps(S, indent=1))

    with open(f"{outdir}/plotdata.json", "w") as f:
        json.dump(per_arm_plot, f)


if __name__ == "__main__":
    main()
