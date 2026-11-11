#!/usr/bin/env python3
"""E28 analysis: the registered statistics of PREREG.md, the comparison with the orientation-only radius, results.json.

    .venv/bin/python experiments/rebuttal/E28_joint_radius_edges/analyse.py --out results/rebuttal/E28

Exit 2 when a row carries a consistency failure or the population is incomplete (full run only), 0 otherwise.
"""

from __future__ import annotations

import argparse
import collections
import json
import statistics
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
B, BOOT_SEED = 2000, 20261107
BAR = 0.8
N_POPULATION, N_INVALID, N_INVALID_MISSING = 984, 534, 452


def load_rows(out: Path) -> tuple[list[dict], list[dict]]:
    rows: dict[str, dict] = {}
    errors = []
    for p in sorted(Path(out).glob("rows.*.jsonl")):
        for line in p.open():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            rows.setdefault(r["unit"], r)
    for p in sorted(Path(out).glob("errors.*.jsonl")):
        errors += [json.loads(line) for line in p.open() if line.strip()]
    done = set(rows)
    errors = [e for e in errors if e["unit"] not in done]
    return list(rows.values()), errors


def cluster_boot(rows, pred, n_boot=B, seed=BOOT_SEED) -> dict:
    by = collections.defaultdict(list)
    for r in rows:
        by[r["gid"]].append(1 if pred(r) else 0)
    gids = sorted(by)
    if not gids:
        return {"count": 0, "n": 0, "share": None, "lo": None, "hi": None, "n_graphs": 0}
    num = np.array([sum(by[g]) for g in gids], float)
    den = np.array([len(by[g]) for g in gids], float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(gids), size=(n_boot, len(gids)))
    s = num[idx].sum(1) / den[idx].sum(1)
    return {"count": int(num.sum()), "n": int(den.sum()), "share": float(num.sum() / den.sum()),
            "lo": float(np.percentile(s, 2.5)), "hi": float(np.percentile(s, 97.5)), "n_graphs": len(gids)}


def fmt(c: dict) -> str:
    if c["share"] is None:
        return "n/a"
    return f"{c['count']:>4}/{c['n']:<4} {100 * c['share']:5.1f} % [{100 * c['lo']:5.1f}, {100 * c['hi']:5.1f}]"


def auc(scores_pos, scores_neg) -> float | None:
    """P(score of an invalid row > score of a valid row), ties counted one half."""
    if not scores_pos or not scores_neg:
        return None
    pos, neg = np.array(scores_pos, float), np.array(scores_neg, float)
    gt = (pos[:, None] > neg[None, :]).sum()
    eq = (pos[:, None] == neg[None, :]).sum()
    return float((gt + 0.5 * eq) / (len(pos) * len(neg)))


def inv_r(r) -> float:
    return 0.0 if r is None else 1.0 / r


def bucket_orient(r: dict) -> str:
    if r.get("r_orient") is not None:
        return str(min(r["r_orient"], 3))
    st = r.get("r_orient_status") or ""
    if st in ("unreached", "no_retractable_edges"):
        return "inf"
    if st.startswith("gt_"):
        return ">depth"
    return "censored"


def bucket_joint(r: dict) -> str:
    if r.get("r_joint") is not None:
        return str(r["r_joint"])
    return ">2" if r.get("r_joint_status") == "gt_2" else "censored"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="results/rebuttal/E28")
    ap.add_argument("--n-boot", type=int, default=B)
    args = ap.parse_args(argv)
    out = Path(args.out)
    if not out.is_absolute():
        out = REPO / out
    rows, errors = load_rows(out)
    cfgs = [json.loads(p.read_text()) for p in sorted(out.glob("config.*.json"))]
    smoke = any(c.get("smoke") for c in cfgs)
    fails = [r for r in rows if r.get("consistency_fail")]
    scored = [r for r in rows if not r.get("consistency_fail") and r.get("reconstruct_ok") == 1]
    inv = [r for r in scored if r["z_valid_at_truth"] == 0]
    val = [r for r in scored if r["z_valid_at_truth"] == 1]
    invm = [r for r in inv if r.get("true_missing_in_A")]
    complete = (len(scored) == N_POPULATION and len(invm) == N_INVALID_MISSING and len(inv) == N_INVALID)
    nb = args.n_boot

    def boot(rs, pred):
        return cluster_boot(rs, pred, n_boot=nb)

    L = []
    L.append(f"E28 rows {len(rows)} (population {N_POPULATION}), scored {len(scored)}: {len(inv)} invalid "
             f"({len(invm)} with a true edge missing in A), {len(val)} valid; smoke={smoke}")
    L.append(f"consistency failures: {len(fails)}; rows lost to an error: {len(errors)}")
    inc1 = [r for r in scored if not r.get("depth1_complete")]
    cen2 = [r for r in scored if r.get("r_joint_status") == "censored_at_2"]
    L.append(f"depth 1 incomplete (time limit): {len(inc1)}; depth 2 censored: {len(cen2)}")

    le2 = lambda r: r.get("r_joint") is not None and r["r_joint"] <= 2  # noqa: E731
    cov = lambda r: r.get("depth1_complete") == 1 and r.get("covered_pair") == 1  # noqa: E731
    s1 = boot(invm, le2)
    s2 = boot(invm, cov)
    med_w = statistics.median([r["n_witness_pairs"] for r in invm]) if invm else None
    med_k = statistics.median([r["n_k"] for r in invm]) if invm else None
    L.append("")
    L.append(f"REGISTERED on the {len(invm)} invalid rows with a true edge missing in An_true({{x,y}} u Z):")
    L.append(f"  S1 joint radius 1 or 2 (censored counted as not):          {fmt(s1)}   bar {BAR}")
    L.append(f"  S2 a true missing edge is a distance-1 addition witness:   {fmt(s2)}   bar {BAR}")
    L.append(f"  S3 burden: median witness pairs {med_w} against median orientations in K {med_k}"
             f"   (dies if the first is larger)")
    verdict = {
        "S1_met": None if s1["share"] is None else bool(s1["share"] >= BAR),
        "S2_met": None if s2["share"] is None else bool(s2["share"] >= BAR),
        "S3_met": None if med_w is None else bool(med_w <= med_k),
    }
    L.append(f"  verdicts: {verdict}")

    L.append("")
    L.append("reported with it, not registered")
    L.append(f"  S1 on all invalid rows:                {fmt(boot(inv, le2))}")
    L.append(f"  S1 on valid rows (no ranking power if as high): {fmt(boot(val, le2))}")
    L.append(f"  joint radius exactly 1, invalid | valid: {fmt(boot(inv, lambda r: r.get('r_joint') == 1))} | "
             f"{fmt(boot(val, lambda r: r.get('r_joint') == 1))}")
    L.append(f"  S2 with the true direction only:      {fmt(boot(invm, lambda r: cov(r) and r.get('covered_oriented') == 1))}")
    L.append(f"  a true missing edge lies in the search region: {fmt(boot(invm, lambda r: bool(r.get('true_missing_in_region'))))}")
    L.append(f"  adding all true missing edges in A breaks Z:   "
             f"{fmt(boot([r for r in invm if r.get('all_true_missing_added_breaks') is not None], lambda r: r['all_true_missing_added_breaks'] == 1))}")
    rt = collections.Counter(str(r.get("r_true_missing")) for r in invm)
    L.append(f"  fewest true missing edges (true direction, up to 3) that break Z: {dict(sorted(rt.items()))}")
    for name, rs in (("invalid", inv), ("valid", val)):
        if rs:
            q = np.percentile([r["n_witness_pairs"] for r in rs], [10, 50, 90])
            qk = np.percentile([r["n_k"] for r in rs], [10, 50, 90])
            L.append(f"  witness pairs p10/p50/p90 on {name} rows: {q.tolist()}; orientations in K: {qk.tolist()}")

    L.append("")
    L.append("per row: orientation-only radius (rows) against joint radius (columns)")
    cols = ["1", "2", ">2", "censored"]
    for name, rs in (("invalid", inv), ("valid", val)):
        tab = collections.Counter((bucket_orient(r), bucket_joint(r)) for r in rs)
        L.append(f"  {name}: r_orient \\ r_joint  " + "  ".join(f"{c:>8}" for c in cols))
        for ro in ["1", "2", "3", ">depth", "inf", "censored"]:
            L.append(f"    {ro:>8}               " + "  ".join(f"{tab[(ro, c)]:>8}" for c in cols))
    lower = sum(1 for r in scored if r.get("r_joint") is not None and
                (r.get("r_orient") is None or r["r_joint"] < r["r_orient"]))
    L.append(f"  rows where the joint radius is below the orientation-only radius: {lower}/{len(scored)}")
    a_or = auc([inv_r(r.get("r_orient")) for r in inv], [inv_r(r.get("r_orient")) for r in val])
    a_jo = auc([inv_r(r.get("r_joint")) for r in inv], [inv_r(r.get("r_joint")) for r in val])
    a_w = auc([r["n_witness_pairs"] for r in inv], [r["n_witness_pairs"] for r in val])
    L.append(f"  AUC for invalidity: 1/r_orient {a_or}, 1/r_joint {a_jo}, witness pairs {a_w}")

    L.append("")
    L.append("per network: n invalid-missing, S1, S2, median witness pairs, median |K|")
    per_net = {}
    for net in sorted({r["network"] for r in scored}):
        rs = [r for r in invm if r["network"] == net]
        if not rs:
            continue
        a, b_ = sum(map(le2, rs)), sum(map(cov, rs))
        mw = statistics.median([r["n_witness_pairs"] for r in rs])
        mk = statistics.median([r["n_k"] for r in rs])
        per_net[net] = {"n": len(rs), "S1": a, "S2": b_, "median_witness_pairs": mw, "median_n_k": mk}
        L.append(f"  {net:<11} {len(rs):>4}  {a:>4}  {b_:>4}  {mw:>6}  {mk:>5}")
    secs = [r["seconds"] for r in scored if "seconds" in r]
    if secs:
        L.append(f"\nseconds per row: median {statistics.median(secs):.1f}, max {max(secs):.0f}, total {sum(secs) / 3600:.2f} h")

    text = "\n".join(L)
    print(text)
    res = {"n_rows": len(rows), "n_scored": len(scored), "n_invalid": len(inv), "n_invalid_missing": len(invm),
           "n_valid": len(val), "smoke": smoke, "complete": complete, "n_consistency_fail": len(fails),
           "n_errors": len(errors), "depth1_incomplete": len(inc1), "depth2_censored": len(cen2),
           "S1": s1, "S2": s2, "S3": {"median_witness_pairs": med_w, "median_n_k": med_k}, "verdict": verdict,
           "auc": {"inv_r_orient": a_or, "inv_r_joint": a_jo, "witness_pairs": a_w},
           "joint_below_orient": lower, "per_network": per_net,
           "consistency_ok": (not fails) and (smoke or complete)}
    (out / "results.json").write_text(json.dumps(res, indent=1))
    if fails or (not smoke and not complete):
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
