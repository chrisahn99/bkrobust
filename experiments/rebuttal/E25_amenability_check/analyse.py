"""E25 analysis: the registered shares (a) and (b), the comparison with the first author's
amenability radius, timings, and every exception with its graph. Writes results.json.

    .venv/bin/python experiments/rebuttal/E25_amenability_check/analyse.py --out results/rebuttal/E25
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
UNREACHED = -1
CHRIS_CSV = ROOT / "experiments/rebuttal/_data/chris_e10/amenability_radius.csv"


def cp_interval(k: int, n: int, alpha: float = 0.05):
    """95% Clopper-Pearson interval (scipy's beta quantiles; exact closed form at k = 0 or n)."""
    if n == 0:
        return None
    if k == n:
        return [round((alpha / 2) ** (1 / n), 5), 1.0]
    if k == 0:
        return [0.0, round(1 - (alpha / 2) ** (1 / n), 5)]
    from scipy.stats import beta
    return [round(float(beta.ppf(alpha / 2, k, n - k + 1)), 5),
            round(float(beta.ppf(1 - alpha / 2, k + 1, n - k)), 5)]


def share(k: int, n: int) -> dict:
    return {"k": k, "n": n, "share": round(k / n, 5) if n else None, "ci95": cp_interval(k, n)}


def finals(r: dict):
    sc, hb, am = r["scan"], r["hybrid"], r["amen_search"]
    rv = sc["r_val"] if sc["r_val_decided"] else (hb["r"] if hb["exact"] else None)
    ra = sc["r_amen"] if sc["r_amen_decided"] else (am["r"] if am["exact"] else None)
    return rv, ra


def med(v):
    v = [x for x in v if x is not None]
    return round(statistics.median(v), 6) if v else None


def analyse_group(rows: list[dict]) -> dict:
    ok = [r for r in rows if r.get("status") == "ok"]
    a_all = [r for r in ok if r["scan"]["r_val_decided"] and r["scan"]["r_val"] not in (None, UNREACHED)
             and r["scan"].get("fail_level", {}).get("complete")]
    a_k = sum(r["scan"]["fail_level"]["n_failing_amenable"] == 0 for r in a_all)
    a_wit = [r for r in ok if r["scan"].get("fail_level")]
    a_wit_k = sum(r["scan"]["fail_level"]["n_failing_amenable"] == 0 for r in a_wit)
    b_rows = [(r, *finals(r)) for r in ok]
    b_dec = [(r, rv, ra) for r, rv, ra in b_rows if rv is not None and ra is not None]
    b_k = sum(rv == ra for _, rv, ra in b_dec)
    b_fin = [(r, rv, ra) for r, rv, ra in b_dec if rv != UNREACHED or ra != UNREACHED]
    b_fin_k = sum(rv == ra for _, rv, ra in b_fin)
    lit_agree = sum(r["scan"]["fail_level"]["literal_agree"] for r in a_wit)
    lit_und = sum(r["scan"]["fail_level"]["literal_undecided"] for r in a_wit)
    lit_n = sum(r["scan"]["fail_level"]["n_failing"] for r in a_wit)
    p42 = Counter()
    for r in ok:
        for k2 in ("amenable_states", "pa_invalid", "y_in_pa"):
            p42[k2] += r["scan"]["prop42"][k2]
    reasons = Counter()
    for r in a_wit:
        reasons.update(r["scan"]["fail_level"]["reasons"])
    t_call = [r["t_amen_call"] for r in ok]
    t_hyb = [r["hybrid"]["seconds"] for r in ok]
    t_am = [r["amen_search"]["seconds"] for r in ok]
    t_scan = [r["scan"]["seconds"] for r in ok]
    ratio_hyb_call = [h / c for h, c in zip(t_hyb, t_call) if c > 0]
    ratio_hyb_am = [h / a for h, a in zip(t_hyb, t_am) if a > 0]
    return {
        "n_rows": len(rows), "status": dict(Counter(r.get("status") for r in rows)), "n_ok": len(ok),
        "registered_a_every_nearest_failing_state_non_amenable": share(a_k, len(a_all)),
        "a_found_failing_states_non_amenable_incl_incomplete_levels": share(a_wit_k, len(a_wit)),
        "a_undecided_units": len([r for r in ok if r["scan"].get("fail_level")
                                  and not r["scan"]["fail_level"]["complete"]]),
        "registered_b_radius_equals_amenability_radius": share(b_k, len(b_dec)),
        "b_finite_only": share(b_fin_k, len(b_fin)),
        "b_undecided_units": len(b_rows) - len(b_dec),
        "infinite_radius_decided_by_top_state_rule_not_exhausted": sum(
            bool(r["scan"].get("decided_by_top_state")) for r in ok),
        "radius_distribution": dict(Counter(str(rv) for _, rv, _ in b_dec).most_common()),
        "nearest_failing_states_total": lit_n, "gac_reason_at_nearest_failing_states": dict(reasons),
        "literal_b_amenability_agrees": {"agree": lit_agree, "undecided": lit_und,
                                         "of": lit_n},
        "prop42_pa_valid_when_amenable": dict(p42),
        "hard_failures_units": sum(bool(r.get("hard_failures")) for r in rows),
        "timing_seconds_median": {"one_call_amenability_at_g0": med(t_call), "hybrid_radius": med(t_hyb),
                                  "amenability_radius_search": med(t_am), "level_scan": med(t_scan),
                                  "ratio_hybrid_over_one_call": med(ratio_hyb_call),
                                  "ratio_hybrid_over_amenability_search": med(ratio_hyb_am)},
        "timing_seconds_total": {"one_call": round(sum(t_call), 3), "hybrid": round(sum(t_hyb), 3),
                                 "amenability_search": round(sum(t_am), 3), "level_scan": round(sum(t_scan), 3)},
    }


def chris_comparison(rows: list[dict]) -> dict:
    if not CHRIS_CSV.exists():
        return {"missing": str(CHRIS_CSV)}
    his = list(csv.DictReader(CHRIS_CSV.open()))
    ours_rel = {}
    for r in rows:
        if r.get("source") == "released" and r.get("status") == "ok":
            ours_rel.setdefault((r["network"], r["arm"], r["x"], r["y"]), []).append(r)
    ours_syn = {(r["which"], r["instance_id"]): r for r in rows
                if r.get("source") == "synthetic" and r.get("status") == "ok"}
    out = {"his_rows": len(his), "his_equal_column_true": sum(h["equal"] == "True" for h in his)}
    matched = agree_amen = agree_val = 0
    disagreements = []
    unmatched = Counter()
    for h in his:
        if h["src"] == "release":
            cands = ours_rel.get((h["network"], h["arm"], h["x"], h["y"]), [])
            r = cands[0] if cands else None
        else:
            r = ours_syn.get((h["src"], h["instance_id"]))
        if r is None:
            unmatched[h["src"]] += 1
            continue
        rv, ra = finals(r)
        if rv is None or ra is None:
            unmatched[h["src"] + "_undecided_ours"] += 1
            continue
        matched += 1
        ok_a = int(h["r_amen"]) == ra
        ok_v = int(h["r_val"]) == rv
        agree_amen += ok_a
        agree_val += ok_v
        if not (ok_a and ok_v) and len(disagreements) < 25:
            disagreements.append({"key": r["key"], "his_r_val": h["r_val"], "his_r_amen": h["r_amen"],
                                  "our_r_val": rv, "our_r_amen": ra})
    out.update({"matched": matched, "unmatched": dict(unmatched),
                "our_r_amen_equals_his": share(agree_amen, matched),
                "our_r_val_equals_his": share(agree_val, matched),
                "registered_b_on_his_queries_our_radius_equals_his_amenability_radius": None,
                "disagreements": disagreements})
    # (b) on his queries: our radius against HIS amenability radius
    k = n = 0
    for h in his:
        r = (ours_rel.get((h["network"], h["arm"], h["x"], h["y"]), [None])[0] if h["src"] == "release"
             else ours_syn.get((h["src"], h["instance_id"])))
        if r is None:
            continue
        rv, _ = finals(r)
        if rv is None:
            continue
        n += 1
        k += rv == int(h["r_amen"])
    out["registered_b_on_his_queries_our_radius_equals_his_amenability_radius"] = share(k, n)
    return out


def exceptions(rows: list[dict]) -> list[dict]:
    exc = []
    for r in rows:
        if r.get("status") != "ok":
            continue
        rv, ra = finals(r)
        fl = r["scan"].get("fail_level", {})
        if fl.get("n_failing_amenable", 0) or (rv is not None and ra is not None and rv != ra) \
                or r.get("hard_failures"):
            exc.append({"key": r["key"], "r_val": rv, "r_amen": ra,
                        "n_failing_amenable": fl.get("n_failing_amenable"),
                        "examples": fl.get("exceptions"), "hard_failures": r.get("hard_failures"),
                        "graphs": r.get("exception_graphs")})
    return exc


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=ROOT / "results/rebuttal/E25")
    args = ap.parse_args(argv)
    rows = {}
    for line in (args.out / "rows.jsonl").open():
        try:
            r = json.loads(line)
        except json.JSONDecodeError:
            continue
        rows[r["key"]] = r  # last record per unit wins (a resumed error is replaced)
    rows = list(rows.values())
    cfg = json.loads((args.out / "config.json").read_text())
    groups = {"all": rows}
    for s in ("released", "e11", "synthetic"):
        groups[s] = [r for r in rows if r.get("source") == s]
    res = {"experiment": "E25", "config": cfg, "n_rows_written": len(rows),
           "n_units_selected": cfg.get("n_units"),
           "groups": {g: analyse_group(v) for g, v in groups.items() if v},
           "first_author_comparison": chris_comparison(rows),
           "exceptions": exceptions(rows)}
    a = res["groups"]["all"]["registered_a_every_nearest_failing_state_non_amenable"]
    b = res["groups"]["all"]["registered_b_radius_equals_amenability_radius"]
    bh = res["first_author_comparison"].get(
        "registered_b_on_his_queries_our_radius_equals_his_amenability_radius") or {}
    res["prediction_met"] = (a["share"] == 1.0 and b["share"] == 1.0 and bh.get("share") in (1.0, None))
    (args.out / "results.json").write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps({"a": a, "b": b, "b_his": bh, "n_exceptions": len(res["exceptions"]),
                      "prediction_met": res["prediction_met"],
                      "timing_all": res["groups"]["all"]["timing_seconds_median"]}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
