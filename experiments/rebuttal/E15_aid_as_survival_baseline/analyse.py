"""E15 analysis: does the AID rank survival as well as the radius?

    PYTHONPATH=src .venv/bin/python experiments/rebuttal/E15_aid_as_survival_baseline/analyse.py \
        --out results/rebuttal/E15

Reads every units_shard*.jsonl in --out, writes results.json there, prints a table.
The registered statistic is PREREG.md's: on the 480 untargeted instances,
Delta = tau_b(radius, AUC_frac) - tau_b(-AID_cpdag graph-level, AUC_frac).
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.stats import kendalltau

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "src"))

from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.robustness.paired_resample import (  # noqa: E402
    paired_cluster_bootstrap,
    paired_comparison,
)

ENDPOINT = "AUC_frac"
RADIUS = "radius_gac"
REGISTERED_BASELINE = "neg_aid_cpdag_graph_norm"
MIN_N = 30
SENTINELS = {RADIUS: (UNREACHED,)}

OURS = [f"aid_{ref}_{scope}" for ref in ("cpdag", "oracle") for scope in ("graph_norm", "row_norm", "pair")]
GADJID = [f"gadjid_{ref}_{name}_aid_norm" for ref in ("cpdag", "oracle") for name in ("oset", "parent", "ancestor")]


def load(out_dir: Path) -> list[dict]:
    recs = {}
    for f in sorted(out_dir.glob("units_shard*.jsonl")):
        for line in f.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                if r["unit_id"] in recs:
                    raise SystemExit(f"duplicate unit {r['unit_id']} in {f}")
                recs[r["unit_id"]] = r
    return list(recs.values())


def flat(rec: dict) -> dict:
    """One analysis row; AID columns negated so every predictor reads 'higher = safer'."""
    def num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    row = {
        "instance_id": rec["unit_id"],
        "dataset": rec["dataset"],
        "arm": rec.get("stored_arm"),
        "component_size": rec.get("stored_component_size"),
        "n_tiers": rec.get("stored_n_tiers"),
        ENDPOINT: num(rec.get(f"stored_{ENDPOINT}")),
        RADIUS: rec.get(RADIUS),
        "separation": num(rec.get("stored_separation")),
        "shd_truth": num(rec.get("stored_shd_truth")),
    }
    for col in OURS + GADJID:
        v = num(rec.get(col))
        row[f"neg_{col}"] = None if v is None else -v
    return row


def marginal(rows, col, n_boot, seed=0):
    pts = [(r[col], r[ENDPOINT]) for r in rows
           if r.get(col) is not None and r[ENDPOINT] is not None and r[col] not in SENTINELS.get(col, ())]
    if len(pts) < MIN_N:
        return {"n": len(pts), "status": "undetermined"}
    x, y = np.array(pts, dtype=float).T
    tau = kendalltau(x, y).statistic
    if not np.isfinite(tau):
        return {"n": len(pts), "status": "constant predictor"}
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        i = rng.integers(0, len(x), len(x))
        t = kendalltau(x[i], y[i]).statistic
        if np.isfinite(t):
            boots.append(t)
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return {"n": len(pts), "tau_b": float(tau), "ci": [float(lo), float(hi)], "status": "ok"}


def paired(rows, baseline, n_boot, with_loo=False):
    fn = paired_comparison if with_loo else paired_cluster_bootstrap
    if sum(1 for r in rows if r.get(baseline) is not None) < MIN_N:
        return {"status": "undetermined", "n": sum(1 for r in rows if r.get(baseline) is not None)}
    res = fn(rows, RADIUS, baseline, ENDPOINT, "instance_id", n_boot=n_boot, seed=0, sentinels=SENTINELS)
    res.pop("loo_by_cluster", None)
    return res


def strata(rows):
    un = [r for r in rows if r["dataset"] == "untargeted"]
    tg = [r for r in rows if r["dataset"] == "targeted"]
    out = {"untargeted|pooled": un}
    for nt in sorted({r["n_tiers"] for r in un}):
        out[f"untargeted|n_tiers={nt}"] = [r for r in un if r["n_tiers"] == nt]
    for arm in sorted({r["arm"] for r in tg}):
        sub = [r for r in tg if r["arm"] == arm]
        out[f"targeted|{arm}|pooled"] = sub
        for c in sorted({r["component_size"] for r in sub}, key=lambda s: int(s)):
            out[f"targeted|{arm}|c={c}"] = [r for r in sub if r["component_size"] == c]
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/rebuttal/E15")
    ap.add_argument("--n-boot", type=int, default=10_000, help="untargeted pooled stratum (registered)")
    ap.add_argument("--n-boot-secondary", type=int, default=2_000, help="every other stratum")
    args = ap.parse_args(argv)
    out_dir = Path(args.out) if Path(args.out).is_absolute() else REPO / args.out

    recs = load(out_dir)
    status = Counter(r["status"] for r in recs)
    failures = Counter(f for r in recs for f in r.get("consistency_failures", []))
    ok = [flat(r) for r in recs if r["status"] == "ok"]
    agree = {ref: Counter(r.get(f"gadjid_{ref}_oset_agrees") for r in recs) for ref in ("cpdag", "oracle")}
    gadjid_err = Counter(k for r in recs for k in r if k.startswith("gadjid_") and k.endswith("_error"))

    baselines = [f"neg_{c}" for c in OURS + GADJID] + ["separation", "shd_truth"]
    groups = strata(ok)
    table = {}
    for name, rows in groups.items():
        entry = {"n": len(rows), "marginal": {}, "paired": {}}
        pooled = name.endswith("pooled")
        nb = args.n_boot if name == "untargeted|pooled" else args.n_boot_secondary
        entry["n_boot"] = nb
        for col in [RADIUS] + baselines:
            if pooled:
                entry["marginal"][col] = marginal(rows, col, nb)
        for b in baselines:
            registered = name == "untargeted|pooled" and b == REGISTERED_BASELINE
            entry["paired"][b] = paired(rows, b, nb, with_loo=registered)
        table[name] = entry

    reg = table["untargeted|pooled"]["paired"][REGISTERED_BASELINE]
    opp = table["untargeted|pooled"]["marginal"].get(REGISTERED_BASELINE, {})
    rad = table["untargeted|pooled"]["marginal"].get(RADIUS, {})
    if reg.get("status") != "ok" or opp.get("status") != "ok" or rad.get("status") != "ok":
        verdict = "undetermined"
    else:
        # AID with the opposite sign ranking at least as well also counts against us.
        opposite_sign_wins = -opp["tau_b"] >= rad["tau_b"]
        verdict = ("prediction holds" if reg["ci_lo_2p5"] > 0 and not opposite_sign_wins
                   else "against the authors")
    results = {
        "experiment": "E15", "criterion": "gac",
        "repo_commits": sorted({r.get("repo_commit") for r in recs}),
        "n_units": len(recs), "status_counts": dict(status), "consistency_failures": dict(failures),
        "gadjid_oset_agreement": {k: {str(a): b for a, b in v.items()} for k, v in agree.items()},
        "gadjid_errors": dict(gadjid_err),
        "registered": {"statistic": f"tau_b({RADIUS}) - tau_b({REGISTERED_BASELINE}) on untargeted pooled",
                       "result": reg, "verdict": verdict},
        "strata": table,
    }
    (out_dir / "results.json").write_text(json.dumps(results, indent=2, default=str))

    print(f"E15: {len(recs)} units, status {dict(status)}, failures {dict(failures) or 'none'}")
    print(f"gadjid oset agreement with ours: {results['gadjid_oset_agreement']}; errors {dict(gadjid_err) or 'none'}")
    print(f"\n{'stratum':30s} {'baseline':38s} {'n':>5s} {'tau_r':>7s} {'tau_b':>7s} {'delta':>7s}  95% CI")
    for name, entry in table.items():
        for b, res in entry["paired"].items():
            if res.get("status") != "ok":
                print(f"{name:30s} {b:38s} {res.get('n', ''):>5}  {res.get('status')}")
                continue
            print(f"{name:30s} {b:38s} {res['n']:>5d} {res['tau_a_point']:7.3f} {res['tau_b_point']:7.3f} "
                  f"{res['delta_point']:7.3f}  [{res['ci_lo_2p5']:.3f}, {res['ci_hi_97p5']:.3f}]")
    print(f"\nREGISTERED ({results['registered']['statistic']}): {verdict}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
