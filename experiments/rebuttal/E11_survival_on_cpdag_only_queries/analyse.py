"""E11 analysis: within-state tau_b of radius vs AUC_frac on each frame (see PREREG.md).

The statistic is the pooled-pair within-state tau_b of Appendix G, the same one
experiments/llm_paired_v2.py computes. Pair counts are taken once per state and the
network-cluster bootstrap re-weights them by how often each network is drawn, with
the same random draws as llm_paired_v2._within_state_bootstrap; both functions are
called on the data as a check, and a disagreement exits with status 2.

    .venv/bin/python experiments/rebuttal/E11_survival_on_cpdag_only_queries/analyse.py --out results/rebuttal/E11
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "experiments"))

import llm_paired_v2 as base  # noqa: E402

UNREACHED = -1
ENDPOINT = "AUC_frac"
MIN_STATES = 5
PUBLISHED = {"tau_b": 0.43, "ci": [0.07, 0.69], "source": "Section 6.3, corrected cohort, back-door radius"}


def load_units(out: Path) -> list[dict[str, Any]]:
    by_key: dict[str, dict[str, Any]] = {}
    for f in sorted(out.glob("units*.jsonl")):
        for line in f.read_text().splitlines():
            if not line.strip():
                continue
            r = json.loads(line)
            old = by_key.get(r["key"])
            if old is not None and (old.get("radius"), old.get(ENDPOINT)) != (r.get("radius"), r.get(ENDPOINT)):
                raise SystemExit(f"conflicting duplicate unit {r['key']}")
            by_key[r["key"]] = r
    return list(by_key.values())


def states(units: list[dict[str, Any]], rcol: str) -> dict[tuple, list[dict[str, Any]]]:
    """(network, condition) -> scorable units, keeping only states where the radius varies."""
    by: dict[tuple, list[dict[str, Any]]] = defaultdict(list)
    for u in units:
        rv, ev = u.get(rcol), u.get(ENDPOINT)
        if u.get("status") != "ok" or rv is None or ev is None or rv == UNREACHED:
            continue
        by[(u["network"], u["condition"])].append({**u, "radius": float(rv), ENDPOINT: float(ev)})
    return {k: v for k, v in by.items() if len({r["radius"] for r in v}) > 1}


def pair_counts(xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
    i, j = np.triu_indices(len(xs), 1)
    sx, sy = np.sign(xs[i] - xs[j]), np.sign(ys[i] - ys[j])
    prod = sx * sy
    return np.array([(prod > 0).sum(), (prod < 0).sum(),
                     ((sx == 0) & (sy != 0)).sum(), ((sy == 0) & (sx != 0)).sum()], dtype=float)


def tau_from(c: np.ndarray) -> float | None:
    C, D, Tx, Ty = c
    den = math.sqrt((C + D + Tx) * (C + D + Ty))
    return (C - D) / den if den > 0 else None


def within_state(st: dict[tuple, list[dict[str, Any]]], n_boot: int, n_perm: int, seed: int = 0) -> dict[str, Any]:
    if not st:
        return {"tau_b": None, "n_states": 0, "n_units": 0, "n_networks": 0, "verdict": "undetermined"}
    keys = list(st)
    arrs = {k: (np.array([r["radius"] for r in st[k]]), np.array([r[ENDPOINT] for r in st[k]])) for k in keys}
    counts = {k: pair_counts(*arrs[k]) for k in keys}
    total = sum(counts.values())
    tau = tau_from(total)
    networks = sorted({k[0] for k in keys})
    net_counts = np.array([sum(counts[k] for k in keys if k[0] == n) for n in networks])
    # same generator and draw sequence as base._within_state_bootstrap
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        drawn = rng.integers(0, len(networks), size=len(networks))
        t = tau_from(np.bincount(drawn, minlength=len(networks)) @ net_counts)
        if t is not None:
            vals.append(t)
    lo = hi = None
    if len(vals) >= 100:
        lo, hi = (float(v) for v in np.percentile(vals, [2.5, 97.5]))
    p = None
    if n_perm:
        prng = np.random.default_rng(seed)
        hits = 0
        for _ in range(n_perm):
            c = sum(pair_counts(prng.permutation(arrs[k][0]), arrs[k][1]) for k in keys)
            t = tau_from(c)
            hits += t is not None and tau is not None and t >= tau
        p = (1 + hits) / (1 + n_perm)
    return {"tau_b": tau, "ci_lo_2p5": lo, "ci_hi_97p5": hi, "n_boot": n_boot, "n_boot_valid": len(vals),
            "perm_p_one_sided": p, "n_perm": n_perm, "C_D_Tx_Ty": [int(v) for v in total],
            "n_states": len(keys), "n_units": sum(len(v) for v in st.values()), "n_networks": len(networks),
            "n_states_more_concordant": int(sum(counts[k][0] > counts[k][1] for k in keys)),
            "verdict": verdict(tau, lo, hi, len(keys))}


def verdict(tau, lo, hi, n_states) -> str:
    if n_states < MIN_STATES or tau is None or lo is None:
        return "undetermined"
    if lo > 0:
        return "holds: interval excludes zero"
    if hi < 0:
        return "reversed: interval below zero"
    return "against: interval includes zero"


def check_against_base(st: dict[tuple, list[dict[str, Any]]], res: dict[str, Any]) -> list[str]:
    """The vectorised statistic must equal llm_paired_v2's, point and (short) bootstrap."""
    if not st:
        return []
    errs = []
    ref = base._stratified_tau_b(st, "radius", ENDPOINT)
    if ref["tau_b"] is None or res["tau_b"] is None or abs(ref["tau_b"] - res["tau_b"]) > 1e-12:
        errs.append(f"point tau_b {res['tau_b']} != llm_paired_v2 {ref['tau_b']}")
    lo, hi, _ = base._within_state_bootstrap(st, "radius", ENDPOINT, n_boot=200, seed=0)
    mine = within_state(st, n_boot=200, n_perm=0)
    if lo is not None and (abs(lo - mine["ci_lo_2p5"]) > 1e-12 or abs(hi - mine["ci_hi_97p5"]) > 1e-12):
        errs.append(f"bootstrap ({mine['ci_lo_2p5']}, {mine['ci_hi_97p5']}) != llm_paired_v2 ({lo}, {hi})")
    return errs


def reproduction(units: list[dict[str, Any]], committed: Path, corrected: Path) -> dict[str, Any]:
    """Table 1 frame rescored with the back-door radius vs the committed panel (+ the 94 corrected units)."""
    if not committed.is_file():
        return {"available": False}
    ref: dict[tuple, dict[str, str]] = {}
    for path in (committed, corrected):
        if path.is_file():
            for r in csv.DictReader(path.open(encoding="utf-8")):
                if r.get("status") == "ok" and r.get("naming") == "real":
                    ref[(r["condition"], r["network"], r["x"], r["y"])] = r
    n = n_rad = n_auc = 0
    examples = []
    for u in units:
        if u["frame"] != "table1" or u.get("status") != "ok":
            continue
        r = ref.get((u["condition"], u["network"], u["x"], u["y"]))
        if r is None:
            continue
        n += 1
        rad_ok = r["radius"] != "" and u.get("bd_radius") is not None and float(r["radius"]) == float(u["bd_radius"])
        auc_ok = (r[ENDPOINT] != "" and u.get(ENDPOINT) is not None
                  and abs(float(r[ENDPOINT]) - float(u[ENDPOINT])) <= 1e-9)
        n_rad += rad_ok
        n_auc += auc_ok
        if (not rad_ok or not auc_ok) and len(examples) < 10:
            examples.append({"key": u["key"], "committed_radius": r["radius"], "bd_radius": u.get("bd_radius"),
                             "committed_auc": r[ENDPOINT], "auc": u.get(ENDPOINT)})
    return {"available": True, "n_matched_ok_units": n, "n_committed_ok_units": len(ref),
            "radius_agree": n_rad, "auc_agree": n_auc,
            "frac_radius_agree": n_rad / n if n else None, "frac_auc_agree": n_auc / n if n else None,
            "examples_disagree": examples}


def fmt(r: dict[str, Any]) -> str:
    if r.get("tau_b") is None:
        return "undetermined"
    ci = f"[{r['ci_lo_2p5']:.3f}, {r['ci_hi_97p5']:.3f}]" if r.get("ci_lo_2p5") is not None else "[n/a]"
    return f"{r['tau_b']:+.3f} {ci}"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(REPO / "results/rebuttal/E11"))
    ap.add_argument("--n-boot", type=int, default=10_000)
    ap.add_argument("--n-perm", type=int, default=2_000)
    ap.add_argument("--committed-units", default=str(REPO / "results/axis_robustness_llm/analysis_units.csv"))
    ap.add_argument("--corrected-units", default=str(REPO / "results/axis_robustness_llm_v2/empty_valid_units.csv"))
    args = ap.parse_args(argv)
    out = Path(args.out)
    units = load_units(out)
    if not units:
        raise SystemExit(f"no units*.jsonl under {out}")
    real = [u for u in units if u.get("naming") == "real"]

    def sel(frame: str, extra=lambda u: True) -> list[dict[str, Any]]:
        return [u for u in real if u["frame"] == frame and extra(u)]

    analyses = {
        "registered__cpdag_only__gac": (sel("cpdag_only"), "radius", args.n_perm),
        "cpdag_only__backdoor": (sel("cpdag_only"), "bd_radius", 0),
        "cpdag_only_without_table1_pairs__gac": (sel("cpdag_only", lambda u: not u.get("in_table1_frame")), "radius", 0),
        "table1__gac": (sel("table1"), "radius", 0),
        "table1__backdoor": (sel("table1"), "bd_radius", 0),
    }
    results: dict[str, Any] = {}
    errors: list[str] = []
    for name, (us, rcol, n_perm) in analyses.items():
        st = states(us, rcol)
        res = within_state(st, args.n_boot, n_perm)
        if name == "registered__cpdag_only__gac":
            errors += check_against_base(st, res)
        res["radius_column"] = rcol
        res["n_units_in_population"] = len(us)
        res["status_counts"] = dict(sorted(_count(us, "status").items()))
        res["n_ok_unreached_excluded"] = sum(1 for u in us if u.get("status") == "ok" and u.get(rcol) == UNREACHED)
        results[name] = res
    per_cond = {}
    for c in sorted({u["condition"] for u in real}):
        st = states([u for u in sel("cpdag_only") if u["condition"] == c], "radius")
        per_cond[c] = within_state(st, args.n_boot, 0)
    xc = [u["xcheck_enum"] for u in units if u.get("xcheck_enum")]
    summary = {
        "experiment": "E11", "criterion": "gac", "published_reference": PUBLISHED,
        "registered": results["registered__cpdag_only__gac"],
        "analyses": results, "per_condition_cpdag_only_gac": per_cond,
        "reproduction_table1_backdoor_vs_committed": reproduction(
            units, Path(args.committed_units), Path(args.corrected_units)),
        "n_units_total": len(units), "smoke": any(u.get("smoke") for u in units),
        "commits": sorted({u.get("commit") for u in units if u.get("commit")}),
        "criteria_recorded": sorted({u.get("criterion") for u in units if u.get("criterion")}),
        "n_xchecked": len(xc), "n_xcheck_exact": sum(1 for x in xc if x.get("exact")),
        "consistency_failures_on_disk": sum(len(f.read_text().splitlines()) for f in out.glob("consistency_failures*.jsonl")),
        "analysis_errors": errors,
    }
    (out / "results.json").write_text(json.dumps(summary, indent=1, default=str))

    print(f"E11 within-state tau_b (radius vs {ENDPOINT}), network-cluster bootstrap n={args.n_boot}")
    print(f"{'analysis':42s} {'tau_b [95% CI]':26s} {'states':>6s} {'units':>6s} {'nets':>5s}  verdict")
    for name, r in results.items():
        print(f"{name:42s} {fmt(r):26s} {r['n_states']:6d} {r['n_units']:6d} {r['n_networks']:5d}  {r['verdict']}")
    print(f"{'published (corrected, back-door)':42s} +0.430 [0.070, 0.690]")
    print(f"registered permutation p (one-sided): {results['registered__cpdag_only__gac'].get('perm_p_one_sided')}")
    rep = summary["reproduction_table1_backdoor_vs_committed"]
    if rep.get("available"):
        print(f"reproduction vs committed: {rep['n_matched_ok_units']} matched, radius agree {rep['radius_agree']}, "
              f"AUC agree {rep['auc_agree']}")
    for c, r in per_cond.items():
        print(f"  {c:20s} {fmt(r):26s} states={r['n_states']}")
    print(f"wrote {out / 'results.json'}")
    if errors:
        print("ANALYSIS CHECK FAILED: " + "; ".join(errors), file=sys.stderr)
        return 2
    return 0


def _count(us: list[dict[str, Any]], col: str) -> dict[str, int]:
    c: dict[str, int] = defaultdict(int)
    for u in us:
        c[str(u.get(col))] += 1
    return c


if __name__ == "__main__":
    raise SystemExit(main())
