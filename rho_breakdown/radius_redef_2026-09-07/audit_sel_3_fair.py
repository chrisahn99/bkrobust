"""SEL_A audit 3 -- the fair version.

Two corrections to audit 2, both in SEL_A's favour:
  * separate the harness defect (tau(D) computed through O*(D), which is not the
    causal effect when D is non-amenable) from the definition's own behaviour, by
    restricting to problems where O*(D_true) IS a valid adjustment set;
  * score the retrieval with random tie-breaking and report the full ranking, not
    just argmax, plus the recall of "Delta_f > 0 at all".
"""
import sys, itertools, time
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from bkrobust.demo.evaluate import (optimal_adjustment_set_dag, is_valid_adjustment_set_dag,
                                    adjusted_estimand)
from audit_sel_2_coverage import build, methods, Z196, calibrated_width

NOBS = 10000


def amenable(pr):
    d = pr["dag"]
    z = optimal_adjustment_set_dag(d, pr["x"], pr["y"])
    return is_valid_adjustment_set_dag(d, pr["x"], pr["y"], z)


def run(n_false, nprob, seed):
    rng = np.random.default_rng(seed)
    rows, tried, t0 = [], 0, time.time()
    while len(rows) < nprob and tried < 200000 and time.time()-t0 < 300:
        tried += 1
        pr = build(rng, n_false)
        if pr is None:
            continue
        cache = {}
        try:
            iv, se, delta, D1, D2, A1, A2 = methods(pr, NOBS, cache)
        except Exception:
            continue
        rows.append(dict(pr=pr, iv=iv, se=se, delta=delta, D1=D1, A1=A1, A2=A2,
                         amen=amenable(pr)))
    return rows


def table(rows, tag):
    print(f"\n--- {tag}  (n={len(rows)}) ---")
    print(f"  {'method':10s} {'raw cov':>8s} {'raw width':>10s} {'lambda':>8s} "
          f"{'cal cov':>8s} {'HONEST WIDTH':>13s}")
    for name in ("point", "SEL_A1", "SEL_A2", "RHIG_1", "RHIG_2", "blanket"):
        ivs = [r["iv"][name] for r in rows]
        ses = [r["se"] for r in rows]
        tau = np.array([r["pr"]["tau_true"] for r in rows])
        th0 = np.array([r["pr"]["theta0"] for r in rows])
        L = np.array([t - (t-iv[0]) - Z196*s for iv, s, t in zip(ivs, ses, th0)])
        H = np.array([t + (iv[1]-t) + Z196*s for iv, s, t in zip(ivs, ses, th0)])
        raw_cov = np.mean((tau >= L) & (tau <= H))
        lam, cov, wid = calibrated_width(ivs, ses, tau, th0)
        print(f"  {name:10s} {raw_cov:8.3f} {np.mean(H-L):10.4f} {lam:8.2f} {cov:8.3f} {wid:13.4f}")


if __name__ == "__main__":
    for n_false in (1, 2):
        rows = run(n_false, 400, 777 + n_false)
        am = [r for r in rows if r["amen"]]
        print(f"\n{'='*78}\n{n_false} FALSE claim(s), problems={len(rows)}, "
              f"of which O*(true DAG) is a valid adjustment set: {len(am)} "
              f"({100*len(am)/len(rows):.1f}%)")
        table(rows, "ALL problems")
        table(am, "AMENABLE only -- the sound-case test")

        # soundness check, order 1, amenable subset, exactly n_false false claims
        bad = 0
        for r in am:
            lo, hi = r["iv"]["SEL_A1"]
            t = r["pr"]["tau_true"]
            if not (lo - r["delta"] <= t <= hi + r["delta"]):
                bad += 1
        print(f"  [soundness] SEL_A1 (+/- delta) misses tau_true on the amenable subset: "
              f"{bad}/{len(am)} = {100*bad/max(len(am),1):.1f}%")
        badr = 0
        for r in am:
            lo, hi = r["iv"]["RHIG_1"]
            t = r["pr"]["tau_true"]
            if not (lo - r["delta"] <= t <= hi + r["delta"]):
                badr += 1
        print(f"  [soundness] RHIG_1  (+/- delta) misses tau_true on the amenable subset: "
              f"{badr}/{len(am)} = {100*badr/max(len(am),1):.1f}%")

        # retrieval, done properly
        rng = np.random.default_rng(0)
        p1_d = p1_b = p1_r = 0.0
        recall_any = 0.0
        n_sc = 0
        rank_sum = 0.0
        for r in am:
            pr = r["pr"]
            if abs(pr["tau_true"] - pr["theta0"]) <= 1e-9:
                continue                      # nothing to retrieve: no damage
            n_sc += 1
            D = np.asarray(r["D1"], float)
            m = len(D)
            F = pr["false_idx"]
            noise = rng.random(m) * 1e-12
            top = int(np.argmax(D + noise))
            p1_d += top in F
            touch = np.array([1.0 if (pr["x"] in e or pr["y"] in e) else 0.0
                              for e in pr["K"]])
            p1_b += int(np.argmax(touch + rng.random(m)*1e-9)) in F
            p1_r += len(F) / m
            recall_any += any(D[i] > 0 for i in F)
            order = np.argsort(-(D + noise))
            rank_sum += min(int(np.where(order == i)[0][0]) for i in F) + 1
        if n_sc:
            print(f"  [retrieval on {n_sc} DAMAGING amenable problems] "
                  f"precision@1: Delta {p1_d/n_sc:.3f} | 'touches X or Y' {p1_b/n_sc:.3f} "
                  f"| random {p1_r/n_sc:.3f}")
            print(f"     Delta_f > 0 for at least one false claim: {recall_any/n_sc:.3f}; "
                  f"mean rank of the best-ranked false claim: {rank_sum/n_sc:.2f} "
                  f"of mean {np.mean([len(r['pr']['K']) for r in am]):.2f} claims")
            # how often is the TOP-ranked claim a TRUE one with positive Delta?
            tp = 0
            for r in am:
                pr = r["pr"]
                D = np.asarray(r["D1"], float)
                if D.max() <= 0:
                    continue
                top = int(np.argmax(D))
                tp += top not in pr["false_idx"]
            print(f"     top-ranked claim is a TRUE claim (load-bearing but correct) in "
                  f"{tp}/{sum(1 for r in am if np.asarray(r['D1']).max()>0)} of the problems "
                  f"where the screen fires")
