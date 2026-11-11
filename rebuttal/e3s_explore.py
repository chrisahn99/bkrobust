"""E3S exploratory (not registered) descriptives: fastest among solver arms only,
solve-only paired ratios, explicit-search time by radius on P-E3."""
import json, collections
import numpy as np
from common import OUT
from e3s_analyse import ARMS, CAP, solved, t_of
D = OUT / "e3s"
rows = [json.loads(l) for l in open(D / "rows.jsonl")]
by = collections.defaultdict(dict)
for r in rows:
    by[r["qid"]][r["arm"]] = r
SOLV = [a for a in ARMS if a != "explicit"]
out = {}
for pop in ("P-E3", "P-SATR", "P-SCALE", "ALL"):
    for lab, lo, hi in (("1-3", 1, 3), (">=4", 4, 99), ("4-6", 4, 6), (">=7", 7, 99)):
        qs = [q for q, a in by.items() if (pop == "ALL" or a["cpsat_ladder"]["pop"] == pop) and lo <= a["cpsat_ladder"]["r_ref"] <= hi]
        if not qs:
            continue
        fast = collections.Counter(min(SOLV, key=lambda a: t_of(by[q][a])) for q in qs)
        fast_all = collections.Counter(min(ARMS, key=lambda a: t_of(by[q][a])) for q in qs)
        so = {}
        for a in ("cpsat_min", "maxsat_rc2", "sat_cadical", "milp_scip"):
            rat = [by[q][a]["solve_s"] / max(by[q]["cpsat_min"]["solve_s"], 1e-6) for q in qs
                   if solved(by[q][a]) and solved(by[q]["cpsat_min"])]
            so[a] = float(np.median(rat)) if rat else None
        tr = {a: float(np.median([by[q][a]["trans_s"] for q in qs if by[q][a].get("trans_s") is not None]))
              for a in ("milp_scip", "maxsat_rc2", "sat_cadical")}
        bd = float(np.median([by[q]["cpsat_min"]["build_s"] for q in qs]))
        out[f"{pop} {lab}"] = {"n": len(qs), "fastest_among_solvers": dict(fast), "fastest_all": dict(fast_all),
                               "solve_only_ratio_to_cpsat_min": so, "trans_median_s": tr, "e1_build_median_s": bd}
ex = collections.defaultdict(list)
for q, a in by.items():
    if a["explicit"]["pop"] == "P-E3":
        ex[a["explicit"]["r_ref"]].append((round(t_of(a["explicit"]), 3), round(t_of(a["cpsat_min"]), 3)))
out["P-E3 explicit vs cpsat_min by r"] = {r: {"n": len(v), "explicit_med": float(np.median([x for x, _ in v])),
                                             "explicit_max": max(x for x, _ in v), "cpsat_min_med": float(np.median([y for _, y in v])),
                                             "cpsat_min_max": max(y for _, y in v),
                                             "explicit_unsolved": sum(x >= CAP for x, _ in v)} for r, v in sorted(ex.items())}
(D / "explore.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
