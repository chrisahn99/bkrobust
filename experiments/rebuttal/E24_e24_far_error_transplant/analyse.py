#!/usr/bin/env python3
"""E24 analysis: the registered statistic (invalid share at d = 2, clustered Wilson), a table, results.json.

    .venv/bin/python experiments/rebuttal/E24_e24_far_error_transplant/analyse.py --out results/rebuttal/E24

Exit 0 when every consistency check holds, 2 otherwise.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
BAR = 0.05
Z = 1.959964
LEVELS = ("none", "1", "2", "3", "all")
E5_POP = {"pc": 408, "ges": 579}  # certifiable A_ORACLE_R learned rows in E5 (results/rebuttal/E5)
E5_FULL_PC = (0.5385, 0.6593)  # E5 registered statistic, PC n = 200 and n = 5000
ANCHOR = (25, 0)  # E5 true-CPDAG anchor: certifiable rows, invalid


def load_rows(out: Path) -> list[dict]:
    rows = {}
    for p in sorted(out.glob("rows.*.jsonl")):
        for line in p.open():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            rows[r["unit"]] = r
    return list(rows.values())


def wilson(p: float, n: float) -> tuple[float | None, float | None]:
    if n <= 0:
        return None, None
    d = 1 + Z * Z / n
    c = (p + Z * Z / (2 * n)) / d
    h = Z * math.sqrt(max(p * (1 - p) / n + Z * Z / (4 * n * n), 0.0)) / d
    return max(0.0, c - h), min(1.0, c + h)


def clustered_share(units: list[dict]) -> dict:
    """Pooled share invalid, Wilson interval on n / design effect, clusters = learned graphs."""
    by_g: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for u in units:
        by_g[u["gid"]][0] += 1 - int(u["z_valid_at_truth"])
        by_g[u["gid"]][1] += 1
    k = sum(v[0] for v in by_g.values())
    n = sum(v[1] for v in by_g.values())
    if n == 0:
        return {"k": 0, "n": 0, "share": None, "ci95": [None, None], "deff": None, "n_graphs": 0}
    p = k / n
    G = len(by_g)
    deff = 1.0
    if 0 < p < 1 and G > 1:
        var = G / (G - 1) * sum((a - p * m) ** 2 for a, m in by_g.values()) / n ** 2
        deff = max(1.0, var / (p * (1 - p) / n))
    lo, hi = wilson(p, n / deff)
    return {"k": k, "n": n, "share": round(p, 4), "ci95": [round(lo, 4), round(hi, 4)],
            "deff": round(deff, 3), "n_eff": round(n / deff, 1), "n_graphs": G}


def far_only(u: dict) -> bool:
    return u.get("err_min_hop") is None or u["err_min_hop"] >= 2


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/rebuttal/E24")
    args = ap.parse_args(argv)
    out = Path(args.out)
    out = out if out.is_absolute() else REPO / out
    rows = load_rows(out)
    cfgs = [json.loads(p.read_text()) for p in sorted(out.glob("config.*.json"))]
    smoke = any(c.get("smoke") for c in cfgs)
    fails = [(r["unit"], r["consistency_fail"]) for r in rows if r.get("consistency_fail")]

    g0 = [r for r in rows if r["level"] == "g0"]
    anchor = [r for r in rows if r["level"] == "anchor"]
    hyb = [r for r in rows if r["unit"].startswith("h|")]
    pop = Counter(r["method"] for r in g0 if r["status"] == "ok")
    a_ok = [r for r in anchor if r["status"] == "ok"]
    a_inv = sum(1 - r["z_valid_at_truth"] for r in a_ok)
    checks = {
        "row_failures": len(fails),
        "criteria": sorted({r["criterion"] for r in rows}),
        "repo_commits": sorted({r["repo_commit"] for r in rows}),
        "anchor": [len(a_ok), a_inv],
        "population": dict(pop),
        "repair_crosschecked_rows": sum(1 for r in hyb if r.get("repair_crosschecked")),
    }
    ok = not fails and checks["criteria"] == ["gac"] and len(checks["repo_commits"]) == 1 and a_inv == 0
    if not smoke:
        ok = ok and tuple(checks["anchor"]) == ANCHOR and dict(pop) == E5_POP
    checks["consistency_ok"] = bool(ok)

    table = []
    for level in LEVELS:
        for method in ("pc", "ges", "pooled"):
            us = [r for r in hyb if r["level"] == level and (method == "pooled" or r["method"] == method)]
            st = Counter(r.get("hstatus") for r in us)
            com = [r for r in us if r.get("hstatus") == "ok"]
            s = clustered_share(com)
            sf = clustered_share([r for r in com if far_only(r)])
            table.append({"level": level, "method": method, "units": len(us),
                          "no_extension": st.get("no_extension", 0), "cyclic": st.get("cyclic", 0),
                          "no_committed_set": len(us) - len(com) - st.get("no_extension", 0) - st.get("cyclic", 0),
                          "committed": len(com), "invalid": s["k"], "share": s["share"], "ci95": s["ci95"],
                          "deff": s["deff"], "far_only_n": sf["n"], "far_only_share": sf["share"],
                          "far_only_ci95": sf["ci95"]})
    per_net = []
    for net in sorted({r["network"] for r in hyb}):
        com = [r for r in hyb if r["level"] == "2" and r["network"] == net and r.get("hstatus") == "ok"]
        s = clustered_share(com)
        per_net.append({"network": net, "committed": s["n"], "invalid": s["k"], "share": s["share"],
                        "ci95": s["ci95"]})
    g0_ok = [r for r in g0 if r["status"] == "ok"]
    g0_share = clustered_share(g0_ok)

    reg_units = [r for r in hyb if r["level"] == "2" and r.get("hstatus") == "ok"]
    reg = clustered_share(reg_units)
    if reg["share"] is None:
        verdict = "no_data"
    elif reg["share"] >= BAR:
        verdict = "against_authors"
    else:
        verdict = "as_predicted"
    reg.update({"level": 2, "methods": "pc+ges pooled", "bar": BAR, "verdict": verdict,
                "anchor_true_class": {"certifiable": len(a_ok), "invalid": a_inv},
                "full_learned_pc_e5": E5_FULL_PC,
                "discarded_per_level": {t["level"]: {"cyclic": t["cyclic"], "no_extension": t["no_extension"]}
                                        for t in table if t["method"] == "pooled"}})

    print(f"E24 analysis of {out} ({'smoke' if smoke else 'full'}): {len(rows)} rows, "
          f"{len(g0_ok)} certifiable units ({dict(pop)}), consistency {'ok' if ok else 'FAILED'}")
    print(f"true-class anchor: {len(a_ok)} certifiable, {a_inv} invalid (E5: 25, 0)")
    print(f"full learned graph (G0, this run): invalid {g0_share['k']}/{g0_share['n']} = {g0_share['share']} "
          f"{g0_share['ci95']}")
    print()
    print("d    | method | units | no ext | cyclic | no set | committed | invalid | share  | 95% CI (clustered) | far-only n, share")
    for t in table:
        print(f"{t['level']:<4} | {t['method']:<6} | {t['units']:>5} | {t['no_extension']:>6} | {t['cyclic']:>6} | "
              f"{t['no_committed_set']:>6} | {t['committed']:>9} | {t['invalid']:>7} | {str(t['share']):<6} | "
              f"{str(t['ci95']):<18} | {t['far_only_n']}, {t['far_only_share']}")
    print()
    print("d = 2 by network:", "; ".join(f"{p['network']} {p['invalid']}/{p['committed']} = {p['share']}"
                                        for p in per_net))
    print()
    print(f"REGISTERED: invalid-at-truth share at d = 2 (PC and GES pooled) = {reg['k']}/{reg['n']} = "
          f"{reg['share']}, 95% clustered Wilson {reg['ci95']} (design effect {reg['deff']}, "
          f"{reg['n_graphs']} graphs); bar {BAR}; verdict {verdict}")
    print(f"  against: true-class anchor {len(a_ok)} certifiable, {a_inv} invalid; full learned PC {E5_FULL_PC}")
    print(f"  discarded per d: {reg['discarded_per_level']}")
    if fails:
        for u, f in fails[:20]:
            print(f"  FAIL {u}: {f}")
    res = {"experiment": "E24_e24_far_error_transplant", "smoke": smoke, "registered_statistic": reg,
           "table": table, "d2_by_network": per_net, "g0_full_learned": g0_share,
           "consistency": checks, "consistency_ok": bool(ok)}
    (out / "results.json").write_text(json.dumps(res, indent=1))
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())
