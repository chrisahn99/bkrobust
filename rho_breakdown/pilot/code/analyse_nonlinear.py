"""
Nonlinear-arm summary. Everything here is finite-sample: read the rho>=1 numbers
against the rho=0 estimator-noise baseline printed first.

Usage: python analyse_nonlinear.py ../results/nonlinear_raw.json ../results
"""
import json
import sys
from collections import defaultdict

import numpy as np

RNG = np.random.default_rng(11)


def boot_ci(v, B=2000):
    v = np.asarray(v, float)
    v = v[~np.isnan(v)]
    if not len(v):
        return "n/a"
    bs = v[RNG.integers(0, len(v), size=(B, len(v)))].mean(axis=1)
    return f"{v.mean():.3f} [{np.percentile(bs,2.5):.3f}, {np.percentile(bs,97.5):.3f}]"


def main():
    raw = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "../results/nonlinear_raw.json"))
    outdir = sys.argv[2] if len(sys.argv) > 2 else "../results"
    S = {"n_scm": len(raw)}
    S["sortability"] = {k: boot_ci([r["diag"][k] for r in raw])
                        for k in ("varsort_iscm", "r2sort_iscm")}
    S["arms"] = {}
    for arm in ("generic", "tiered"):
        rs = [r for r in raw if arm in r["arms"]]
        if not rs:
            continue
        A = {"n_scm": len(rs)}
        # rho=0 estimator-noise baseline: correct O*, so deviation = estimation error
        A["BASELINE_rho0_|est-tau|/|tau|"] = boot_ci(
            [abs(r["arms"][arm]["est0"] - r["tau"]) / abs(r["tau"]) for r in rs])
        buckets = defaultdict(int)
        for r in rs:
            for m in r["arms"][arm]["members"]:
                if m["rho"] != 1:
                    continue
                if not m["consistent"]:
                    buckets["meek_inconsistent"] += 1
                elif not m["amenable"]:
                    buckets["not_amenable"] += 1
                elif not m["ostar_changed"]:
                    buckets["ostar_unchanged"] += 1
                elif m["ostar_valid"]:
                    buckets["changed_still_valid"] += 1
                else:
                    buckets["changed_invalid"] += 1
        tot = sum(buckets.values())
        cons = tot - buckets["meek_inconsistent"]
        A["rho1_buckets_frac"] = {k: round(v / tot, 4) for k, v in sorted(buckets.items())}
        A["NUMBER1_frac_consistent_flips_changing_Ostar"] = round(
            (buckets["changed_still_valid"] + buckets["changed_invalid"] +
             buckets["not_amenable"]) / cons, 4)
        A["NUMBER1_frac_consistent_flips_making_Ostar_INVALID"] = round(
            buckets["changed_invalid"] / cons, 4)
        for rho in (1, 2):
            wn, ex = [], []
            for r in rs:
                ests = [r["arms"][arm]["est0"]] + [
                    m["est"] for m in r["arms"][arm]["members"]
                    if m["rho"] <= rho and m["consistent"] and m.get("amenable")]
                lo, hi = min(ests), max(ests)
                wn.append((hi - lo) / abs(r["tau"]))
                ex.append(float(lo > 0 or hi < 0))
            A[f"rho{rho}_width_over_|tau|_mean"] = boot_ci(wn)
            A[f"rho{rho}_width_over_|tau|_median"] = float(np.median(wn))
            A[f"rho{rho}_frac_excludes_0"] = boot_ci(ex)
        S["arms"][arm] = A
    json.dump(S, open(f"{outdir}/summary_nonlinear.json", "w"), indent=1)
    print(json.dumps(S, indent=1))


if __name__ == "__main__":
    main()
