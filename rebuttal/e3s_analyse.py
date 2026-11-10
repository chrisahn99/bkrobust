"""E3S analysis: the registered statistics of PREREGISTRATION.md (section E3S)."""

from __future__ import annotations

import json
import math
from collections import defaultdict

import numpy as np

from common import OUT

D = OUT / "e3s"
CAP = 300.0
ARMS = ["cpsat_ladder", "cpsat_min", "cpsat_core", "milp_scip", "milp_highs", "milp_cbc",
        "maxsat_rc2", "sat_cadical", "explicit"]
RCLS = [("1-3", 1, 3), ("4-6", 4, 6), (">=7", 7, 10**9)]


def wilson(k, n):
    if n == 0:
        return (None, None)
    z = 1.96
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (round(max(0.0, c - h), 3), round(min(1.0, c + h), 3))


def med_ci(xs):
    xs = sorted(xs)
    n = len(xs)
    if n == 0:
        return None, (None, None)
    lo = max(0, int(math.floor(n / 2 - 1.96 * math.sqrt(n) / 2)) - 1)
    hi = min(n - 1, int(math.ceil(n / 2 + 1.96 * math.sqrt(n) / 2)))
    return float(np.median(xs)), (xs[lo], xs[hi])


def solved(r):
    return r.get("status") == "ok" and r.get("exact")


def t_of(r):
    """Time to an exact answer; unsolved counts as the cap (censored)."""
    return float(r["total_s"]) if solved(r) else CAP


def graph_id(qid):
    parts = qid.split(":")
    if parts[0] in ("e3",):
        return ":".join(parts[:5])
    if parts[0] == "scale":
        return ":".join(parts[:6])
    return qid


def rcls(r):
    for lab, lo, hi in RCLS:
        if lo <= r <= hi:
            return lab
    return "other"


def main():
    rows = [json.loads(l) for l in open(D / "rows.jsonl")]
    by_q = defaultdict(dict)
    for r in rows:
        by_q[r["qid"]][r["arm"]] = r
    qinfo = {q: next(iter(a.values())) for q, a in by_q.items()}
    out = {"n_rows": len(rows), "n_queries": len(by_q),
           "status": {}, "mismatch": [], "witness_fail": [], "errors": []}
    for r in rows:
        s = r.get("status", "?").split(":")[0]
        out["status"][s] = out["status"].get(s, 0) + 1
        if solved(r) and r.get("r") is not None and r["r"] != r["r_ref"]:
            out["mismatch"].append((r["qid"], r["arm"], r["r"], r["r_ref"]))
        if solved(r) and r.get("unreached"):
            out["mismatch"].append((r["qid"], r["arm"], "unreached", r["r_ref"]))
        if r.get("witness_ok") is False:
            out["witness_fail"].append((r["qid"], r["arm"]))
        if r.get("status", "").startswith(("error", "crashed")):
            out["errors"].append((r["qid"], r["arm"], r["status"]))
    complete = [q for q, a in by_q.items() if all(x in a for x in ARMS)]
    out["n_queries_all_arms"] = len(complete)

    def stats(qs, label):
        res = {}
        rng = np.random.default_rng(0)
        for arm in ARMS:
            rs = [by_q[q][arm] for q in qs]
            n = len(rs)
            k = sum(solved(r) for r in rs)
            ts = [t_of(r) for r in rs]
            m, ci = med_ci(ts)
            # paired ratio to the ladder, cluster bootstrap over graphs
            ratios = {}
            for q in qs:
                a1 = t_of(by_q[q]["cpsat_ladder"])
                ratios.setdefault(graph_id(q), []).append(t_of(by_q[q][arm]) / max(a1, 1e-6))
            allr = [x for v in ratios.values() for x in v]
            keys = list(ratios)
            boots = []
            for _ in range(1000):
                pick = [keys[i] for i in rng.integers(0, len(keys), len(keys))]
                boots.append(np.median([x for kk in pick for x in ratios[kk]]))
            fastest = 0
            for q in qs:
                tt = {a: t_of(by_q[q][a]) for a in ARMS}
                if solved(by_q[q][arm]) and tt[arm] <= min(tt.values()):
                    fastest += 1
            solve_only = [float(r["solve_s"]) for r in rs if solved(r) and r.get("solve_s") is not None]
            res[arm] = {
                "n": n, "solved": k, "solved_share": round(k / n, 3) if n else None,
                "solved_ci": wilson(k, n),
                "median_s": m, "median_ci": ci, "max_s_solved": max([t for t, r in zip(ts, rs) if solved(r)], default=None),
                "ratio_to_ladder_median": float(np.median(allr)) if allr else None,
                "ratio_ci": (float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))) if allr else None,
                "fastest_share": round(fastest / n, 3) if n else None,
                "solve_only_median_s": float(np.median(solve_only)) if solve_only else None,
                "build_median_s": float(np.median([r["build_s"] for r in rs if r.get("build_s") is not None]))
                if any(r.get("build_s") is not None for r in rs) else None,
                "trans_median_s": float(np.median([r["trans_s"] for r in rs if r.get("trans_s") is not None]))
                if any(r.get("trans_s") is not None for r in rs) else None,
                "size_median": {kk: float(np.median([r[kk] for r in rs if r.get(kk) is not None]))
                                for kk in ("n_vars", "n_cons", "milp_rows", "clauses", "cnf_vars")
                                if any(r.get(kk) is not None for r in rs)},
            }
        return res

    tabs = {}
    for pop in ("P-E3", "P-SATR", "P-SCALE", "ALL"):
        qp = [q for q in complete if pop == "ALL" or qinfo[q]["pop"] == pop]
        tabs[pop] = {"all": {"n": len(qp), "arms": stats(qp, pop)}}
        for lab, lo, hi in RCLS:
            qq = [q for q in qp if lo <= qinfo[q]["r_ref"] <= hi]
            if qq:
                tabs[pop][lab] = {"n": len(qq), "arms": stats(qq, lab)}
    out["tables"] = tabs
    # P-SCALE by n
    sc = {}
    for q in complete:
        if qinfo[q]["pop"] != "P-SCALE":
            continue
        n = by_q[q]["cpsat_ladder"].get("n_nodes")
        sc.setdefault(n, []).append(q)
    out["scale_by_n"] = {n: {"n": len(qs), "arms": {a: {"solved": sum(solved(by_q[q][a]) for q in qs),
                                                        "median_s": med_ci([t_of(by_q[q][a]) for q in qs])[0]}
                                                    for a in ARMS},
                             "radius": sorted(qinfo[q]["r_ref"] for q in qs)} for n, qs in sorted(sc.items())}
    # ladder certificate time
    cert = {}
    for arm in ("cpsat_ladder", "sat_cadical"):
        unsat, sat = [], []
        for q in complete:
            r = by_q[q][arm]
            if solved(r) and r.get("rungs"):
                ok_code = 4 if arm == "cpsat_ladder" else 1
                feas = [t for k, t, st in r["rungs"] if (st in (2, 4) if arm == "cpsat_ladder" else st == 1)]
                inf = [t for k, t, st in r["rungs"] if not (st in (2, 4) if arm == "cpsat_ladder" else st == 1)]
                sat.append(sum(feas))
                unsat.append(sum(inf))
        cert[arm] = {"n": len(sat), "unsat_rungs_median_s": float(np.median(unsat)) if unsat else None,
                     "sat_rung_median_s": float(np.median(sat)) if sat else None,
                     "unsat_share_of_solve_median": float(np.median([u / max(u + s, 1e-9) for u, s in zip(unsat, sat)]))
                     if sat else None}
    out["ladder_certificate"] = cert
    (D / "analysis.json").write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps({k: v for k, v in out.items() if k not in ("tables",)}, indent=1, default=str)[:4000])
    for pop, t in tabs.items():
        for cl, v in t.items():
            print(f"== {pop} {cl} n={v['n']}")
            for arm, s in v["arms"].items():
                print(f"  {arm:13s} solved {s['solved']:5d}/{s['n']:<5d} {s['solved_ci']} med {s['median_s']:.4g} "
                      f"{tuple(round(x, 4) for x in s['median_ci'])} ratio {s['ratio_to_ladder_median']:.3g} "
                      f"fastest {s['fastest_share']} solve_only {s['solve_only_median_s']}")


if __name__ == "__main__":
    main()
