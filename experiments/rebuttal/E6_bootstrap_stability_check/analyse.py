"""E6 analysis: AUC of bootstrap instability for "Z invalid at truth", network-clustered.

    python experiments/rebuttal/E6_bootstrap_stability_check/analyse.py --out results/rebuttal/E6
    python experiments/rebuttal/E6_bootstrap_stability_check/analyse.py --out results/rebuttal/E6_smoke

Reads reference.jsonl and every units_*.jsonl in the folder, writes results.json there.
The registered statistic and its interval are fixed in PREREG.md.
"""

from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

import numpy as np

N_DRAWS = 2000
CI_SEED = 20261106
AGAINST_US = 0.6


def auc(scores: list[float], labels: list[int]) -> float | None:
    """P(score of a positive > score of a negative), ties one half; None without both classes."""
    s = np.asarray(scores, float)
    lab = np.asarray(labels, int)
    pos, neg = s[lab == 1], s[lab == 0]
    if len(pos) == 0 or len(neg) == 0:
        return None
    gt = (pos[:, None] > neg[None, :]).sum()
    eq = (pos[:, None] == neg[None, :]).sum()
    return float((gt + 0.5 * eq) / (len(pos) * len(neg)))


def cluster_ci(rows: list[dict], key: str, rng: np.random.Generator) -> dict:
    """Two-stage cluster bootstrap: networks, then query pairs inside each drawn network."""
    by_net: dict[str, dict[str, list[dict]]] = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rows:
        by_net[r["network"]][f"{r['x']}|{r['y']}"].append(r)
    nets = sorted(by_net)
    vals, skipped = [], 0
    for _ in range(N_DRAWS):
        drawn: list[dict] = []
        for net in rng.choice(nets, size=len(nets), replace=True):
            pairs = sorted(by_net[net])
            for p in rng.choice(pairs, size=len(pairs), replace=True):
                drawn += by_net[net][p]
        a = auc([d[key] for d in drawn], [d["invalid"] for d in drawn])
        if a is None:
            skipped += 1
        else:
            vals.append(a)
    if not vals:
        return {"lo": None, "hi": None, "draws_used": 0, "draws_skipped": skipped}
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return {"lo": round(float(lo), 4), "hi": round(float(hi), 4), "draws_used": len(vals), "draws_skipped": skipped}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/rebuttal/E6", help="the run's output folder")
    args = ap.parse_args(argv)
    out = Path(args.out)
    ref = [json.loads(line) for line in (out / "reference.jsonl").open()]
    units: dict[str, dict] = {}
    for path in sorted(out.glob("units_*.jsonl")):
        for line in path.read_text().splitlines():
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue  # a line cut by a kill
            units.setdefault(rec["unit"], rec)
    config = json.loads((out / "config.json").read_text())

    replicates = [u for u in units.values() if u["b"] >= 0]
    identity = [u for u in units.values() if u["b"] < 0]
    finished = [u for u in replicates if u["status"] == "ok"]
    per_net = collections.Counter(u["network"] for u in finished)
    timeouts = collections.Counter(u["network"] for u in replicates if u["status"] == "timeout")
    criteria = {u["criterion"] for u in units.values()}
    commits = {u["repo_commit"] for u in units.values()}

    agg: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    for u in finished:
        for r in u["rows"]:
            c = agg[r["rid"]]
            c["n"] += 1
            c["valid"] += r["z_valid"]
            c["state_ok"] += r["state_ok"]
            c["z_same"] += r["z_same"]
            c["skel"] += r["skel_nb_stable"]
            c["nb"] += r["state_nb_stable"]
            c["err"] += int(bool(r["err"]))
    rows = []
    for r in ref:
        c = agg.get(r["rid"])
        if not c or not c["n"]:
            continue
        n = c["n"]
        rows.append(
            {
                "rid": r["rid"],
                "network": r["network"],
                "arm": r["arm"],
                "x": r["x"],
                "y": r["y"],
                "invalid": 1 - r["z_valid_at_truth"],
                "n_rep": n,
                "stability": c["valid"] / n,
                "score": 1 - c["valid"] / n,
                "score_state": 1 - c["state_ok"] / n,
                "score_zsame": 1 - c["z_same"] / n,
                "score_skel": 1 - c["skel"] / n,
                "score_nb": 1 - c["nb"] / n,
                "errors": c["err"],
            }
        )
    if not rows:
        print("no finished replicates yet", file=sys.stderr)
        return 1

    rng = np.random.default_rng(CI_SEED)
    head = auc([r["score"] for r in rows], [r["invalid"] for r in rows])
    ci = cluster_ci(rows, "score", rng)
    inv = [r for r in rows if r["invalid"]]
    val = [r for r in rows if not r["invalid"]]
    secondary = {
        key: auc([r[key] for r in rows], [r["invalid"] for r in rows])
        for key in ("score_state", "score_zsame", "score_skel", "score_nb")
    }
    by = {}
    for field in ("network", "arm"):
        by[field] = {}
        for g in sorted({r[field] for r in rows}):
            sub = [r for r in rows if r[field] == g]
            by[field][g] = {
                "rows": len(sub),
                "invalid": sum(r["invalid"] for r in sub),
                "auc": auc([r["score"] for r in sub], [r["invalid"] for r in sub]),
                "mean_stability_invalid": round(float(np.mean([r["stability"] for r in sub if r["invalid"]])), 4) if any(r["invalid"] for r in sub) else None,
                "mean_stability_valid": round(float(np.mean([r["stability"] for r in sub if not r["invalid"]])), 4) if any(not r["invalid"] for r in sub) else None,
            }
    results = {
        "experiment": "E6",
        "statistic": "AUC of (1 - bootstrap share with Z GAC-valid) for Z invalid at truth",
        "auc": None if head is None else round(head, 4),
        "ci95_cluster_network_then_pair": ci,
        "against_authors": None if head is None else bool(head <= AGAINST_US),
        "rows": len(rows),
        "rows_reference": len(ref),
        "rows_invalid": len(inv),
        "rows_valid": len(val),
        "stable_but_wrong_share": round(sum(r["stability"] >= 0.95 for r in inv) / len(inv), 4) if inv else None,
        "valid_flagged_share": round(sum(r["stability"] < 0.5 for r in val) / len(val), 4) if val else None,
        "secondary_auc": {k: None if v is None else round(v, 4) for k, v in secondary.items()},
        "by": by,
        "replicates_finished": dict(per_net),
        "replicates_timeout": dict(timeouts),
        "networks_under_80_replicates": sorted(n for n in config["networks"] if per_net.get(n, 0) < 80),
        "identity_reproduces_reference": {u["network"]: u.get("is_reference_cpdag") for u in identity},
        "rows_with_check_errors": sum(1 for r in rows if r["errors"]),
        "criterion": sorted(criteria),
        "repo_commit": sorted(commits),
        "smoke": config["smoke"],
        "B_planned": config["B"],
    }
    (out / "results.json").write_text(json.dumps(results, indent=1))
    (out / "per_row.json").write_text(json.dumps(rows, indent=1))

    print(f"E6  {'SMOKE  ' if config['smoke'] else ''}criterion {','.join(sorted(criteria))}  commit {','.join(c[:10] for c in sorted(commits))}")
    print(f"rows {len(rows)} of {len(ref)} (invalid at truth {len(inv)}, valid {len(val)}); replicates finished {dict(per_net)}, timeouts {dict(timeouts)}")
    print(f"AUC {results['auc']}  95% CI [{ci['lo']}, {ci['hi']}]  ({ci['draws_used']} draws, {ci['draws_skipped']} skipped)  against authors (<= {AGAINST_US}): {results['against_authors']}")
    print(f"stable but wrong (invalid, s >= 0.95): {results['stable_but_wrong_share']}   valid flagged (s < 0.5): {results['valid_flagged_share']}")
    print(f"secondary AUC: {results['secondary_auc']}")
    print(f"{'group':<22}{'rows':>6}{'inval':>7}{'AUC':>8}{'s|inv':>8}{'s|val':>8}")
    for field in ("network", "arm"):
        for g, d in by[field].items():
            a = "-" if d["auc"] is None else f"{d['auc']:.3f}"
            si = "-" if d["mean_stability_invalid"] is None else f"{d['mean_stability_invalid']:.3f}"
            sv = "-" if d["mean_stability_valid"] is None else f"{d['mean_stability_valid']:.3f}"
            print(f"{field[:3]}:{g:<18}{d['rows']:>6}{d['invalid']:>7}{a:>8}{si:>8}{sv:>8}")
    print(f"identity replicate = reference CPDAG: {results['identity_reproduces_reference']}")
    print(f"-> {out / 'results.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
