"""
ARM 2 analysis: the coverage-width frontier and the vacuity falsifier.

Usage: python analyse_arm2.py ../results/arm2.json
"""
import json
import sys
from collections import defaultdict

import numpy as np

RNG = np.random.default_rng(20260819)


def boot(v, B=2000):
    v = np.asarray(v, dtype=float)
    if len(v) == 0:
        return (float("nan"),) * 3
    idx = RNG.integers(0, len(v), size=(B, len(v)))
    bs = v[idx].mean(axis=1)
    return (round(float(v.mean()), 4), round(float(np.percentile(bs, 2.5)), 4),
            round(float(np.percentile(bs, 97.5)), 4))


def main():
    d = json.load(open(sys.argv[1]))
    scms = [s for s in d["scms"] if not s.get("gate_fail")]
    rhos = [str(r) for r in range(0, 5)]

    S = dict(source=sys.argv[1], n_scm=len(scms), reps=d["reps"], z=d["z"],
             gate_regen_failures=d["gate_regen_failures"])

    # how often a wrong expert is caught for free / loudly, per r
    cf, ld = defaultdict(list), defaultdict(list)
    for s in scms:
        for r, v in s["caught_free"].items():
            cf[r].append(float(v))
        for r, v in s["loud"].items():
            ld[r].append(float(v))
    S["caught_free_rate"] = {r: boot(v) for r, v in sorted(cf.items())}
    S["loud_rate"] = {r: boot(v) for r, v in sorted(ld.items())}

    cov = defaultdict(lambda: defaultdict(list))
    wid = defaultdict(lambda: defaultdict(list))
    wr = defaultdict(lambda: defaultdict(list))
    for s in scms:
        for c in s["cells"]:
            key = (c["r"], c["n"])
            w0 = c["w"]["0"]
            for rho in rhos:
                if rho not in c["cov"]:
                    continue
                cov[key][rho].append(float(c["cov"][rho]))
                wid[key][rho].append(c["w"][rho])
                if w0 > 0:
                    wr[key][rho].append(c["w"][rho] / w0)

    F = {}
    for key in sorted(cov):
        r, n = key
        F[f"r{r}_n{n}"] = {rho: dict(
            coverage=boot(cov[key][rho]),
            median_width=round(float(np.median(wid[key][rho])), 4),
            median_width_ratio=round(float(np.median(wr[key][rho])), 4),
            mean_width_ratio=round(float(np.mean(wr[key][rho])), 4),
            q75_width_ratio=round(float(np.percentile(wr[key][rho], 75)), 4),
            q90_width_ratio=round(float(np.percentile(wr[key][rho], 90)), 4),
            q99_width_ratio=round(float(np.percentile(wr[key][rho], 99)), 4),
            frac_width_ratio_is_1=round(float(np.mean(
                np.abs(np.asarray(wr[key][rho]) - 1.0) < 1e-9)), 4),
            n_cells=len(cov[key][rho])) for rho in rhos if rho in cov[key]}
    S["frontier"] = F

    def g(r, n, rho, f):
        return F.get(f"r{r}_n{n}", {}).get(str(rho), {}).get(f)

    # ---- gates and the pre-registered rule
    cap = g(1, 20000, 0, "coverage")
    S["capability_check"] = dict(
        naive_coverage_r1_n20000=cap,
        pass_=bool(cap is not None and cap[0] < 0.90))
    g7 = g(0, 20000, 0, "coverage")
    S["G7_naive_coverage_r0_n20000"] = dict(value=g7,
        pass_=bool(g7 is not None and abs(g7[0] - 0.95) <= 0.03))

    c11 = g(1, 20000, 1, "coverage")
    w11 = g(1, 20000, 1, "median_width_ratio")
    if c11 is None or w11 is None:
        S["verdict"] = "NO DATA"
    elif c11[0] >= 0.90 and w11 < 3.0:
        S["verdict"] = "SUPPORTED"
    elif c11[0] < 0.80 or w11 > 10.0:
        S["verdict"] = "REFUTED"
    else:
        S["verdict"] = "INCONCLUSIVE"
    S["decision_inputs"] = dict(coverage_r1_rho1_n20000=c11,
                                median_width_ratio_r1_rho1_n20000=w11)
    # the non-definitional cell: ball misses the truth (rho < r)
    S["ball_misses_truth_r2_rho1"] = {
        f"n{n}": dict(coverage=g(2, n, 1, "coverage"),
                      width_ratio=g(2, n, 1, "median_width_ratio"))
        for n in d["n_grid"]}

    print("\n=== COVERAGE / WIDTH FRONTIER (z=1.96) ===")
    print(f"{'cell':>12} {'rho':>4} {'coverage':>22} {'w_ratio med/mean/q90':>26} {'inert%':>7} {'cells':>8}")
    for key in sorted(F, key=lambda k: (int(k.split('_')[0][1:]), int(k.split('_n')[1]))):
        for rho in rhos:
            if rho not in F[key]:
                continue
            v = F[key][rho]
            c = v["coverage"]
            print(f"{key:>12} {rho:>4}  {c[0]:.4f} [{c[1]:.4f},{c[2]:.4f}]"
                  f"  {v['median_width_ratio']:>7.3f}/{v['mean_width_ratio']:>7.3f}/{v['q90_width_ratio']:>7.3f}"
                  f"  {100*v['frac_width_ratio_is_1']:>6.1f} {v['n_cells']:>8}")
    print()
    print(json.dumps({k: v for k, v in S.items() if k != "frontier"}, indent=1))
    with open("../results/arm2_analysis.json", "w") as f:
        json.dump(S, f, indent=1)


if __name__ == "__main__":
    main()
