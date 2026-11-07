"""E19 analysis: the registered statistic of PREREG.md, with intervals, from run.py's output.

    .venv/bin/python experiments/rebuttal/E19_bias_radius_internal_standardisation/analyse.py [--out DIR] [--smoke]
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
UNREACHED = -1  # bkrobust.core.conventions.UNREACHED; restated so analysis needs no imports
ARMS = ("naive", "posthoc", "internal")
EPS_KEYS = ["0.01", "0.02", "0.05", "0.1", "0.25", "0.5", "1"]
PAPER = {"never_reached_eps1": 0.701, "step_share": 0.841, "ramp_share": 0.151}  # App. E.4, naive
AGAINST = {"never_reached_eps1": 0.40, "step_share": 0.60}
N_BOOT = 2000


def wilson(k: int, n: int, z: float = 1.959964) -> tuple[float | None, float | None]:
    if n == 0:
        return None, None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def share(flags: list[bool]) -> dict:
    k, n = int(sum(flags)), len(flags)
    lo, hi = wilson(k, n)
    return {"k": k, "n": n, "share": k / n if n else None, "ci95": [lo, hi]}


def per_arm(rows: list[dict], arm: str) -> dict:
    a = [r["arms"][arm] for r in rows]
    out: dict = {"n": len(a)}
    out["never_reached"] = {e: share([x["r_eps"][e] == UNREACHED for x in a]) for e in EPS_KEYS}
    comp = [x for x in a if x["r_val"] != UNREACHED]
    out["r_eps_gt_r_val"] = {
        e: share([x["r_eps"][e] == UNREACHED or x["r_eps"][e] > x["r_val"] for x in comp]) for e in EPS_KEYS}
    nontriv = [x for x in a if x["staircase_shape"] != "n/a"]
    for s in ("step", "ramp", "flat"):
        out[f"{s}_share"] = share([x["staircase_shape"] == s for x in nontriv])
    out["r0_equals_r_val"] = share([x["r0_d"] == x["r_val"] for x in a])
    full = [x for x in a if not x["beta_up_at_kfalse_is_lower_bound"]]
    out["certificate_holds"] = share([bool(x["certificate_holds"]) for x in full])
    # DEVIATIONS.md: both threshold conventions, tolerance 1e-9 at the threshold
    for conv in ("ge", "gt"):
        out[f"never_reached_{conv}"] = {e: share([x[f"r_eps_{conv}"][e] == UNREACHED for x in a]) for e in EPS_KEYS}
        nt = [x for x in a if x[f"staircase_shape_{conv}"] != "n/a"]
        out[f"step_share_{conv}"] = share([x[f"staircase_shape_{conv}"] == "step" for x in nt])
    out["n_with_threshold_tie"] = sum(1 for x in a if x.get("threshold_ties"))
    out["var_sortability_mean"] = float(np.mean([x["var_sortability"] for x in a])) if a else None
    out["r2_sortability_mean"] = float(np.mean([x["r2_sortability"] for x in a])) if a else None
    return out


def primary_flags(rows: list[dict], arm: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    never = np.array([r["arms"][arm]["r_eps"]["1"] == UNREACHED for r in rows], dtype=float)
    step = np.array([r["arms"][arm]["staircase_shape"] == "step" for r in rows], dtype=float)
    nontriv = np.array([r["arms"][arm]["staircase_shape"] != "n/a" for r in rows], dtype=bool)
    return never, step, nontriv


def paired_diff(rows: list[dict]) -> dict:
    if not rows:
        return {}
    ni, si, ti = primary_flags(rows, "internal")
    nn, sn, tn = primary_flags(rows, "naive")

    def stat(idx: np.ndarray) -> tuple[float, float]:
        d_never = ni[idx].mean() - nn[idx].mean()
        mi, mn = ti[idx], tn[idx]
        si_, sn_ = (si[idx][mi].mean() if mi.any() else np.nan), (sn[idx][mn].mean() if mn.any() else np.nan)
        return d_never, si_ - sn_

    point = stat(np.arange(len(rows)))
    rng = np.random.default_rng(0)
    boots = np.array([stat(rng.integers(0, len(rows), len(rows))) for _ in range(N_BOOT)])
    ci = [np.nanpercentile(boots[:, j], [2.5, 97.5]).tolist() for j in range(2)]
    return {"never_reached_eps1_internal_minus_naive": {"diff": point[0], "ci95": ci[0]},
            "step_share_internal_minus_naive": {"diff": point[1], "ci95": ci[1]},
            "n_boot": N_BOOT, "boot_seed": 0}


def fmt(s: dict) -> str:
    if not s or s.get("share") is None:
        return "   n/a"
    lo, hi = s["ci95"]
    return f"{s['share']:6.1%} [{lo:5.1%},{hi:5.1%}] n={s['n']}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "results" / "rebuttal" / "E19"))
    ap.add_argument("--smoke", action="store_true", help="read OUT/smoke")
    args = ap.parse_args()
    out_dir = Path(args.out) / ("smoke" if args.smoke else "")
    files = sorted(out_dir.glob("units.shard*of*.jsonl"))
    if not files:
        print(f"no units.shard*.jsonl in {out_dir}", file=sys.stderr)
        return 1
    rows: dict[str, dict] = {}
    for f in files:
        for line in f.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                rows[r["instance_id"]] = r  # a unit rewritten after a resume keeps its last line
    allrows = list(rows.values())
    failed = [r for r in allrows if r.get("consistency_failures")]
    excluded = [r for r in allrows if not r.get("status_all_ok")]
    paired = [r for r in allrows if r.get("status_all_ok") and not r.get("consistency_failures")]

    res: dict = {
        "experiment": "E19", "criterion": sorted({r.get("criterion") for r in allrows}),
        "repo_commits": sorted({r.get("repo_commit") for r in allrows}),
        "any_code_dirty": any(r.get("repo_code_dirty") for r in allrows),
        "n_units": len(allrows), "n_paired": len(paired), "n_excluded_timeout_or_error": len(excluded),
        "excluded_fraction": len(excluded) / len(allrows) if allrows else None,
        "n_consistency_failures": len(failed),
        "consistency_failure_examples": [{"instance_id": r["instance_id"], "failures": r["consistency_failures"]}
                                         for r in failed[:10]],
        "arms": {a: per_arm(paired, a) for a in ARMS},
        "paired": paired_diff(paired),
    }
    by_gen: dict[str, list] = defaultdict(list)
    for r in paired:
        by_gen[r["generator"]].append(r)
    res["per_generator"] = {
        g: {a: {"never_reached_eps1": per_arm(rs, a)["never_reached"]["1"],
                "step_share": per_arm(rs, a)["step_share"],
                "var_sortability_mean": per_arm(rs, a)["var_sortability_mean"],
                "r2_sortability_mean": per_arm(rs, a)["r2_sortability_mean"]} for a in ("naive", "internal")}
        for g, rs in sorted(by_gen.items())}
    ia = res["arms"]["internal"]
    never, step = ia["never_reached"]["1"]["share"], ia["step_share"]["share"]
    res["n_tie_only_rows"] = sum(1 for r in allrows if r.get("tie_only_differences"))
    res["tie_only_examples"] = [{"instance_id": r["instance_id"], "differences": r["tie_only_differences"]}
                                for r in allrows if r.get("tie_only_differences")][:20]
    res["registered_both_conventions"] = {
        conv: {"never_reached_eps1_internal": ia[f"never_reached_{conv}"]["1"],
               "step_share_internal": ia[f"step_share_{conv}"],
               "never_reached_eps1_naive": res["arms"]["naive"][f"never_reached_{conv}"]["1"],
               "step_share_naive": res["arms"]["naive"][f"step_share_{conv}"]}
        for conv in ("ge", "gt")}
    res["registered"] = {
        "never_reached_eps1_internal": ia["never_reached"]["1"],
        "step_share_internal": ia["step_share"],
        "against_authors": (never is not None and never < AGAINST["never_reached_eps1"])
                           or (step is not None and step < AGAINST["step_share"]),
        "thresholds": AGAINST, "paper_naive_values": PAPER,
        "reduced_set": bool(res["excluded_fraction"] and res["excluded_fraction"] > 0.05),
    }

    print(f"E19  units {len(allrows)}  paired {len(paired)}  excluded {len(excluded)}  "
          f"consistency failures {len(failed)}  criterion {res['criterion']}")
    print(f"{'statistic':28s} " + " ".join(f"{a:>34s}" for a in ARMS))
    lines = [("never reached, eps=" + e, lambda x, e=e: x["never_reached"][e]) for e in EPS_KEYS]
    lines += [("r_eps > r_val, eps=" + e, lambda x, e=e: x["r_eps_gt_r_val"][e]) for e in ("0.01", "0.25", "1")]
    lines += [("step share", lambda x: x["step_share"]), ("ramp share", lambda x: x["ramp_share"]),
              ("r_0 = r_val", lambda x: x["r0_equals_r_val"]), ("certificate holds", lambda x: x["certificate_holds"])]
    for name, get in lines:
        print(f"{name:28s} " + " ".join(f"{fmt(get(res['arms'][a])):>34s}" for a in ARMS))
    for a in ARMS:
        print(f"{a:9s} var-sortability {res['arms'][a]['var_sortability_mean']}  "
              f"R2-sortability {res['arms'][a]['r2_sortability_mean']}")
    pd_ = res["paired"]
    if pd_:
        for k in ("never_reached_eps1_internal_minus_naive", "step_share_internal_minus_naive"):
            print(f"{k}: {pd_[k]['diff']:+.3f} [{pd_[k]['ci95'][0]:+.3f}, {pd_[k]['ci95'][1]:+.3f}]")
    print(f"against the authors (never<40% or step<60% under internal): {res['registered']['against_authors']}")
    print(f"rows whose only differences are threshold ties: {res['n_tie_only_rows']}")
    for conv, label in (("ge", "bias >= tolerance counts as exceeding"), ("gt", "bias > tolerance counts as exceeding")):
        print(f"convention {conv} ({label}), tie tolerance 1e-9:")
        for name, get in ([("  never reached, eps=" + e, lambda x, e=e, c=conv: x[f"never_reached_{c}"][e]) for e in EPS_KEYS]
                          + [("  step share", lambda x, c=conv: x[f"step_share_{c}"])]):
            print(f"{name:28s} " + " ".join(f"{fmt(get(res['arms'][a])):>34s}" for a in ARMS))
    if failed:
        print(f"WARNING: {len(failed)} units failed a consistency check; see results.json")

    (out_dir / "results.json").write_text(json.dumps(res, indent=2, default=float))
    print(f"wrote {out_dir / 'results.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
