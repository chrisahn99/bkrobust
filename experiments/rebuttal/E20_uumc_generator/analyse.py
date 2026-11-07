"""E20 analysis: the registered statistic of PREREG.md, from the units written by run.py.

    .venv/bin/python experiments/rebuttal/E20_uumc_generator/analyse.py --out results/rebuttal/E20

Reads every units.*.jsonl in --out, writes results.json there and prints a plain-text table.
Exit 2 if any unit carries a consistency failure (then no number is to be reported).
"""
from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
EPS = ["0.02", "0.1", "0.25", "0.4", "0.6", "0.9"]
ARMS = ("naive", "iscm", "uumc")
CONTRAST = {"barley": ("barley", "dg25", "s2225"), "schipf": ("Schipf_2010", "A", "U"),
            "kampen": ("Kampen_2014", "AFF", "CDR")}
UNR = -1
PAPER_NAIVE_DRAW0 = {"S": 22, "C": 32, "R": 22, "B": 21, "n": 49}


def load(out: Path) -> list[dict]:
    rows = {}
    for path in sorted(out.glob("units.*.jsonl")):
        for line in path.open():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            rows[(r["network"], r["x"], r["y"], r["arm"], r["draw"])] = r  # last write wins
    return list(rows.values())


def sig(r):
    return tuple(r["r_mu_rel"][e] for e in EPS)


def counts(rs: list[dict]) -> dict:
    b = [round(r["b_top_rel"], 3) for r in rs]
    s_ = {sig(r) for r in rs}
    c_ = {(sig(r), bb) for r, bb in zip(rs, b, strict=True)}
    r_ = {(r["r_val"], bb) for r, bb in zip(rs, b, strict=True)}
    return {"n": len(rs), "S": len(s_), "C": len(c_), "R": len(r_), "B": len(set(b)),
            "added": len(c_) - len(r_), "n_radius_values": len({r["r_val"] for r in rs})}


def contrast(rs: list[dict], eps: str = "0.6"):
    """Barley crosses eps strictly earlier than both r_val = 1 queries; UNR = infinity."""
    by = {(r["network"], r["x"], r["y"]): r for r in rs}
    if any(CONTRAST[k] not in by for k in CONTRAST):
        return None
    val = {k: by[q]["r_mu_rel"][eps] for k, q in CONTRAST.items()}
    inf = lambda v: float("inf") if v == UNR else v  # noqa: E731
    if val["barley"] == UNR:
        return False
    return inf(val["barley"]) < min(inf(val["schipf"]), inf(val["kampen"]))


def boot_median(xs: list[float], n_boot: int = 4000, seed: int = 0):
    import numpy as np
    if not xs:
        return None, None
    rng = np.random.default_rng(seed)
    a = np.asarray(xs, dtype=float)
    meds = np.median(rng.choice(a, size=(n_boot, len(a)), replace=True), axis=1)
    return float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))


def _num(v):
    return float("nan") if v is None else v


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", default="results/rebuttal/E20")
    args = p.parse_args()
    out = Path(args.out)
    out = out if out.is_absolute() else ROOT / out
    rows = load(out)
    if not rows:
        print(f"no units in {out}")
        return 1
    failures = [r for r in rows if r.get("consistency_failures")]
    status = {}
    for r in rows:
        status[r["status"]] = status.get(r["status"], 0) + 1

    per_draw: dict[str, dict[int, dict]] = {a: {} for a in ARMS}
    sort_: dict[str, dict[str, list[float]]] = {a: {} for a in ARMS}
    for a in ARMS:
        draws = sorted({r["draw"] for r in rows if r["arm"] == a})
        for k in draws:
            rs = [r for r in rows if r["arm"] == a and r["draw"] == k and r["status"] == "ok"]
            if not rs:
                continue
            c = counts(rs)
            c["contrast_60"] = contrast(rs)
            c["contrast_by_eps"] = {e: contrast(rs, e) for e in EPS}
            per_draw[a][k] = c
        for key in ("var_sortability_pair", "r2_sortability_pair",
                    "var_sortability_reisach", "r2_sortability_reisach"):
            sort_[a][key] = [r[key] for r in rows if r["arm"] == a and key in r]

    summary = {}
    for a in ARMS:
        d = per_draw[a]
        if not d:
            continue
        s_vals = [c["S"] for c in d.values()]
        added = [c["added"] for c in d.values()]
        cons = [c["contrast_60"] for c in d.values() if c["contrast_60"] is not None]
        lo, hi = boot_median(s_vals)
        summary[a] = {
            "n_draws": len(d),
            "S_median": statistics.median(s_vals), "S_ci95": [lo, hi],
            "S_min": min(s_vals), "S_max": max(s_vals),
            "added_median": statistics.median(added),
            "draws_added_positive": sum(v > 0 for v in added),
            "contrast_draws_holding": sum(cons), "contrast_draws_defined": len(cons),
            "sortability_mean": {k: (statistics.fmean(v) if v else None) for k, v in sort_[a].items()},
        }

    # Prereg decision, read mechanically; None where an arm is missing.
    verdict = {}
    if "iscm" in summary and "uumc" in summary:
        si, su = summary["iscm"], summary["uumc"]
        def share_hi(s):  # noqa: E306
            return s["contrast_draws_holding"] >= 6
        def share_lo(s):  # noqa: E306
            return s["contrast_draws_holding"] <= 4
        rev_a = (share_hi(si) and share_lo(su)) or (share_hi(su) and share_lo(si))
        rev_b = (si["added_median"] > 0) != (su["added_median"] > 0)
        verdict = {"reversal_a_contrast": rev_a, "reversal_b_added": rev_b,
                   "reversal": rev_a or rev_b,
                   "uumc_profile_adds_nothing": su["S_median"] <= 3,
                   "prediction_S_ge_10": su["S_median"] >= 10,
                   "prediction_added_in_8_of_10": su["draws_added_positive"] >= 8}
    repro = None
    if 0 in per_draw["naive"]:
        c0 = per_draw["naive"][0]
        repro = {k: (c0[k], v, c0[k] == v) for k, v in PAPER_NAIVE_DRAW0.items()}

    res = {"out": str(out), "n_units": len(rows), "status": status,
           "n_consistency_failures": len(failures),
           "consistency_failures": [{k: r[k] for k in ("network", "x", "y", "arm", "draw",
                                                       "consistency_failures")} for r in failures[:50]],
           "naive_draw0_vs_paper": repro, "per_arm": summary, "verdict": verdict,
           "per_draw": {a: {str(k): v for k, v in d.items()} for a, d in per_draw.items()}}
    (out / "results.json").write_text(json.dumps(res, indent=2) + "\n")

    print(f"E20  {out}  units {len(rows)}  status {status}  consistency failures {len(failures)}")
    if repro:
        print("naive draw 0 vs paper (mine, paper, equal): "
              + ", ".join(f"{k} {m}/{p_} {'ok' if e else 'DIFF'}" for k, (m, p_, e) in repro.items()))
    print(f"{'arm':6} {'draws':>5} {'S med':>6} {'95% CI':>13} {'min-max':>8} {'added':>6} "
          f"{'contrast':>9} {'var-sort':>8} {'R2-sort':>8}")
    for a, s in summary.items():
        ci = s["S_ci95"]
        print(f"{a:6} {s['n_draws']:>5} {s['S_median']:>6} {f'[{ci[0]:.1f}, {ci[1]:.1f}]':>13} "
              f"{str(s['S_min']) + '-' + str(s['S_max']):>8} {s['added_median']:>6} "
              f"{str(s['contrast_draws_holding']) + '/' + str(s['contrast_draws_defined']):>9} "
              f"{_num(s['sortability_mean']['var_sortability_pair']):>8.3f} "
              f"{_num(s['sortability_mean']['r2_sortability_pair']):>8.3f}")
    print("radius alone: " + ", ".join(f"{a} {next(iter(per_draw[a].values()))['n_radius_values']}"
                                       for a in summary))
    if verdict:
        print("verdict: " + ", ".join(f"{k}={v}" for k, v in verdict.items()))
    if failures:
        print("CONSISTENCY FAILURES: no number above is to be reported")
    return 2 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
