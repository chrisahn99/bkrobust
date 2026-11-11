"""E7 analysis: the registered paired difference Δ_near, Δ_far with bootstrap intervals.

Reads every units_*.jsonl in --out, pairs each hidden unit's rows with the
nothing-hidden unit of the same network and seed, prints a plain-text table and
writes results.json beside the units. Works on the smoke output too.

    python experiments/rebuttal/E7_one_hidden_variable/analyse.py --out results/rebuttal/E7
"""

from __future__ import annotations

import argparse
import collections
import glob
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
NONE = "-"
UNREACHED = -1
BOOT_SEED = 20261106
REGISTERED_ARM = "A_ORACLE"
AGAINST = 0.10


def load_units(out: Path) -> dict[str, dict]:
    units: dict[str, dict] = {}
    for f in sorted(glob.glob(str(out / "units_*.jsonl"))):
        for line in open(f):
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue  # a line cut by a kill
            units.setdefault(r["unit"], r)
    return units


def stratum(dist) -> str:
    return "near" if dist == 1 else "far"


def build_pairs(units: dict[str, dict]) -> tuple[list[dict], dict]:
    base: dict[tuple, dict] = {}
    for u in units.values():
        if u["h"] == NONE and u["status"] == "ok":
            for r in u["rows"]:
                base[(u["network"], u["seed"], r["q"], r["arm"])] = r
    pairs, missing = [], collections.Counter()
    for u in units.values():
        if u["h"] == NONE or u["status"] != "ok":
            continue
        for r in u["rows"]:
            b = base.get((u["network"], u["seed"], r["q"], r["arm"]))
            if b is None:
                missing[u["network"]] += 1
                continue
            pairs.append({
                "network": u["network"], "seed": u["seed"], "h": u["h"], "q": r["q"], "arm": r["arm"],
                "stratum": stratum(r["dist"]), "dist": r["dist"],
                "common_ancestor": r.get("h_common_ancestor", 0),
                "bad_h": r["invalid_truth"], "bad_0": b["invalid_truth"],
                "commit_h": r.get("committed", 0), "commit_0": b.get("committed", 0),
                "radius_h": r.get("radius_radius"), "radius_exact_h": r.get("radius_exact"),
                "radius_lb_h": r.get("radius_lb"), "radius_status_h": r.get("radius_status"),
                "radius_0": b.get("radius_radius"), "radius_exact_0": b.get("radius_exact"),
            })
    return pairs, dict(missing)


def delta(ps: list[dict]) -> float | None:
    return float(np.mean([p["bad_h"] - p["bad_0"] for p in ps])) if ps else None


def boot_delta(ps: list[dict], B: int, seed: int) -> dict[str, list[float | None]]:
    """Percentile interval of Δ per stratum, resampling hidden variables within each network."""
    # per (network, h): the stratum sums, so a bootstrap draw is a sum over drawn h
    sums: dict[str, dict[str, np.ndarray]] = collections.defaultdict(dict)
    for p in ps:
        s = sums[p["network"]].setdefault(p["h"], np.zeros(4))
        k = 0 if p["stratum"] == "near" else 2
        s[k] += p["bad_h"] - p["bad_0"]
        s[k + 1] += 1
    rng = np.random.default_rng(BOOT_SEED if seed is None else seed)
    nets = sorted(sums)
    mats = [np.array([sums[n][h] for h in sorted(sums[n])]) for n in nets]
    draws = {"near": [], "far": []}
    for _ in range(B):
        tot = np.zeros(4)
        for m in mats:
            tot += m[rng.integers(0, len(m), size=len(m))].sum(axis=0)
        draws["near"].append(tot[0] / tot[1] if tot[1] else np.nan)
        draws["far"].append(tot[2] / tot[3] if tot[3] else np.nan)
    out = {}
    for k, v in draws.items():
        v = np.array(v, dtype=float)
        v = v[~np.isnan(v)]
        out[k] = [round(float(np.percentile(v, 2.5)), 4), round(float(np.percentile(v, 97.5)), 4)] if len(v) else [None, None]
        out[k + "_draws_used"] = int(len(v))
    return out


def no_warning(p: dict) -> int | None:
    """1 when the radius on the hidden graph was 2 or more, or unreached; None when unknown."""
    if p["radius_status_h"] != "ok":
        return None
    if p["radius_exact_h"]:
        return int(p["radius_h"] == UNREACHED or p["radius_h"] >= 2)
    return int(p["radius_lb_h"] is not None and p["radius_lb_h"] >= 2)


def summarise(ps: list[dict]) -> dict:
    by = {s: [p for p in ps if p["stratum"] == s] for s in ("near", "far")}
    out = {}
    for s, q in by.items():
        valid0 = [p for p in q if p["commit_0"] and not p["bad_0"]]
        out[s] = {
            "pairs": len(q),
            "delta": None if delta(q) is None else round(delta(q), 4),
            "invalid_hidden": round(float(np.mean([p["bad_h"] for p in q])), 4) if q else None,
            "invalid_full": round(float(np.mean([p["bad_0"] for p in q])), 4) if q else None,
            "commit_hidden": round(float(np.mean([p["commit_h"] for p in q])), 4) if q else None,
            "commit_full": round(float(np.mean([p["commit_0"] for p in q])), 4) if q else None,
            "became_invalid": round(float(np.mean([p["bad_h"] for p in valid0])), 4) if valid0 else None,
            "became_invalid_n": len(valid0),
        }
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=str(ROOT / "results" / "rebuttal" / "E7"))
    ap.add_argument("--B", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=BOOT_SEED)
    args = ap.parse_args(argv)
    out = Path(args.out)
    cfg = json.loads((out / "config.json").read_text()) if (out / "config.json").exists() else {}
    units = load_units(out)
    if not units:
        print(f"no units in {out}", file=sys.stderr)
        return 2
    problems = [f"{u['unit']}: {p}" for u in units.values() for p in u.get("problems", [])]

    # PC timeouts per network, and the flag the prereg sets at 5%
    tcount = collections.defaultdict(lambda: [0, 0])
    for u in units.values():
        tcount[u["network"]][1] += 1
        tcount[u["network"]][0] += int(u["status"] == "timeout")
    timeouts = {n: {"timeout": a, "units": b, "flag": a / b > 0.05} for n, (a, b) in sorted(tcount.items())}

    pairs, missing_base = build_pairs(units)
    reg = [p for p in pairs if p["arm"] == REGISTERED_ARM]
    main_s = summarise(reg)
    ci = boot_delta(reg, args.B, args.seed)
    for s in ("near", "far"):
        main_s[s]["ci95"] = ci[s]
        main_s[s]["boot_draws_used"] = ci[s + "_draws_used"]

    d_far, d_near = main_s["far"]["delta"], main_s["near"]["delta"]
    verdict = {
        "against_authors": None if d_far is None else bool(d_far > AGAINST),
        "prediction_met": None if d_far is None or d_near is None else bool(
            d_near > d_far and d_far < 0.03 and ci["far"][1] is not None and ci["far"][1] < AGAINST
        ),
    }

    per_network = {n: summarise([p for p in reg if p["network"] == n]) for n in sorted({p["network"] for p in reg})}
    per_seed = {str(s): summarise([p for p in reg if p["seed"] == s]) for s in sorted({p["seed"] for p in reg})}
    degen = summarise([p for p in pairs if p["arm"] == "D_DEGEN"])
    common = summarise([p for p in reg if p["common_ancestor"]])
    became = [p for p in reg if p["commit_0"] and not p["bad_0"] and p["bad_h"]]
    warn = [no_warning(p) for p in became]
    known = [w for w in warn if w is not None]
    radius = {
        "became_invalid_pairs": len(became),
        "radius_known": len(known),
        "no_warning_share": round(float(np.mean(known)), 4) if known else None,
        "by_stratum": {
            s: {"pairs": sum(1 for p in became if p["stratum"] == s),
                "no_warning": sum(1 for p, w in zip(became, warn) if p["stratum"] == s and w == 1)}
            for s in ("near", "far")
        },
        "hidden_radius_status": dict(collections.Counter(str(p["radius_status_h"]) for p in reg if p["commit_h"])),
    }
    base_rows = [r for u in units.values() if u["h"] == NONE and u["status"] == "ok" for r in u["rows"] if r["arm"] == REGISTERED_ARM]
    results = {
        "experiment": "E7",
        "registered_statistic": "Δ_s = mean over pairs of [Z invalid at truth, h hidden] − [same, nothing hidden], arm A_ORACLE, near (hop 1) vs far (hop ≥ 2)",
        "criterion": cfg.get("criterion"),
        "repo_commit": cfg.get("repo_commit"),
        "smoke": cfg.get("smoke"),
        "units": len(units),
        "pairs_registered_arm": len(reg),
        "registered": main_s,
        "verdict": verdict,
        "full_graph_invalid_at_truth": round(float(np.mean([r["invalid_truth"] for r in base_rows])), 4) if base_rows else None,
        "per_network": per_network,
        "per_seed": per_seed,
        "D_DEGEN": degen,
        "h_common_ancestor_of_x_and_y": common,
        "radius": radius,
        "pc_timeouts": timeouts,
        "pairs_without_base_unit": missing_base,
        "consistency_problems": problems,
        "bootstrap": {"B": args.B, "seed": args.seed, "resample": "hidden variables within network"},
    }
    (out / "results.json").write_text(json.dumps(results, indent=1))

    print(f"E7  {len(units)} units, {len(reg)} paired rows (arm {REGISTERED_ARM}), criterion {cfg.get('criterion')}, smoke={cfg.get('smoke')}")
    print(f"{'stratum':8} {'pairs':>7} {'delta':>8} {'95% CI':>18} {'inv|h':>7} {'inv|0':>7} {'became':>7} {'commit|h':>9}")
    for s in ("near", "far"):
        m = main_s[s]
        lo, hi = m["ci95"]
        cis = f"[{lo}, {hi}]" if lo is not None else "-"
        print(f"{s:8} {m['pairs']:>7} {str(m['delta']):>8} {cis:>18} {str(m['invalid_hidden']):>7} "
              f"{str(m['invalid_full']):>7} {str(m['became_invalid']):>7} {str(m['commit_hidden']):>9}")
    print("per network (delta near / far):")
    for n, m in per_network.items():
        print(f"  {n:12} {m['near']['delta']} / {m['far']['delta']}   (pairs {m['near']['pairs']} / {m['far']['pairs']})")
    print(f"D_DEGEN delta near / far: {degen['near']['delta']} / {degen['far']['delta']}")
    print(f"became invalid: {radius['became_invalid_pairs']} pairs, radius gave no warning in {radius['no_warning_share']} of {radius['radius_known']} known")
    print(f"against the authors (delta_far > {AGAINST}): {verdict['against_authors']}; prediction met: {verdict['prediction_met']}")
    flagged = [n for n, t in timeouts.items() if t["flag"]]
    if flagged:
        print(f"FLAG: more than 5% PC timeouts on {flagged}")
    if problems:
        print(f"CONSISTENCY PROBLEMS: {len(problems)} (first: {problems[0]})")
        return 3
    print(f"wrote {out / 'results.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
