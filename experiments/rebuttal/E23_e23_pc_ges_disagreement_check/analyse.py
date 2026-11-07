"""E23 analysis: AUC of PC/GES disagreement for "Z invalid at truth", network-then-query clustered.

    python experiments/rebuttal/E23_e23_pc_ges_disagreement_check/analyse.py --out results/rebuttal/E23
    python experiments/rebuttal/E23_e23_pc_ges_disagreement_check/analyse.py --out results/rebuttal/E23_smoke

Reads every rows.*.jsonl in the folder and writes results.json there. Statistic fixed in PREREG.md.
Exit 0, or 2 when a consistency check failed.
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
BAR = 0.65


def auc(scores, labels) -> float | None:
    """P(score of a positive > score of a negative), ties one half; None without both classes."""
    s, lab = np.asarray(scores, float), np.asarray(labels, int)
    pos, neg = s[lab == 1], s[lab == 0]
    if len(pos) == 0 or len(neg) == 0:
        return None
    gt = (pos[:, None] > neg[None, :]).sum()
    eq = (pos[:, None] == neg[None, :]).sum()
    return float((gt + 0.5 * eq) / (len(pos) * len(neg)))


def cluster_ci(rows: list[dict], key: str, seed: int = CI_SEED) -> dict:
    """Two-stage cluster bootstrap as in E6: networks, then query pairs inside each drawn network."""
    rng = np.random.default_rng(seed)
    by_net: dict = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rows:
        by_net[r["network"]][f"{r['x']}|{r['y']}"].append(r)
    nets = sorted(by_net)
    vals, skipped = [], 0
    for _ in range(N_DRAWS):
        drawn = []
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


def r4(v):
    return None if v is None else round(v, 4)


def summary(sub: list[dict]) -> dict:
    inv = [r for r in sub if r["invalid"]]
    val = [r for r in sub if not r["invalid"]]
    return {
        "rows": len(sub), "invalid": len(inv),
        "auc": r4(auc([r["flag"] for r in sub], [r["invalid"] for r in sub])),
        "auc_or_diff": r4(auc([r["flag_or_diff"] for r in sub], [r["invalid"] for r in sub])),
        "auc_moved": r4(auc([r["flag_moved"] for r in sub], [r["invalid"] for r in sub])),
        "shared_error_share": r4(sum(1 - r["flag"] for r in inv) / len(inv)) if inv else None,
        "false_alarm_share": r4(sum(r["flag"] for r in val) / len(val)) if val else None,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/rebuttal/E23")
    args = ap.parse_args(argv)
    out = Path(args.out)
    rows: dict = {}
    for p in sorted(out.glob("rows.*.jsonl")):
        for line in p.open():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            rows.setdefault(r["unit"], r)
    configs = [json.loads(p.read_text()) for p in sorted(out.glob("config.*.json"))]
    if not rows or not configs:
        print("no rows yet", file=sys.stderr)
        return 1
    smoke = any(c["smoke"] for c in configs)
    fails = [(u, r["consistency_fail"]) for u, r in rows.items() if r["consistency_fail"]]
    e5 = [json.loads(p.read_text()) for p in sorted(out.glob("e5_check.*.json"))]
    sharded = any(c["n_shards"] > 1 for c in configs)
    if sharded and not smoke:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from run import e5_check  # the rebuild check needs every shard's rows
        gids = {r["gid"] for r in rows.values()}
        e5_fails, n_e5 = e5_check(rows, Path(configs[0]["e5_dir"]), gids)
        e5 = [{"e5_rows": n_e5, "fails": e5_fails, "e5_present": n_e5 > 0}]
    for e in e5:
        fails += e["fails"]

    paired = []
    no_partner = collections.Counter()
    for r in rows.values():
        if r["partner"] is None:
            no_partner[r["network"]] += 1
            continue
        paired.append({**r, "invalid": 1 - r["z_valid_at_truth"]})
    head = auc([r["flag"] for r in paired], [r["invalid"] for r in paired])
    ci = cluster_ci(paired, "flag") if head is not None else {"lo": None, "hi": None, "draws_used": 0, "draws_skipped": 0}
    against = None if head is None else bool(head <= BAR or (ci["lo"] is not None and ci["lo"] < 0.5))
    by = {f: {str(g): summary([r for r in paired if r[f] == g]) for g in sorted({r[f] for r in paired})}
          for f in ("network", "n", "direction")}
    results = {
        "experiment": "E23",
        "statistic": "AUC of (Z not GAC-valid in the other method's state) for Z invalid at truth, pooled over n and both directions",
        "auc": r4(head), "ci95_cluster_network_then_pair": ci, "bar": BAR, "against_authors": against,
        **{k: v for k, v in summary(paired).items() if k != "auc"},
        "rows_no_partner_by_network": dict(no_partner),
        "by": by,
        "e5_rebuild_check": e5,
        "consistency_ok": not fails, "consistency_failures": len(fails),
        "criterion": sorted({r["criterion"] for r in rows.values()}),
        "repo_commit": sorted({r["repo_commit"] for r in rows.values()}),
        "smoke": smoke,
    }
    (out / "results.json").write_text(json.dumps(results, indent=1))

    print(f"E23 {'SMOKE  ' if smoke else ''}criterion {','.join(results['criterion'])}  "
          f"commit {','.join(c[:10] for c in results['repo_commit'])}")
    print(f"paired rows {len(paired)} (invalid at truth {results['invalid']}); rows without a partner graph {dict(no_partner)}")
    print(f"REGISTERED  AUC {results['auc']}  95% CI [{ci['lo']}, {ci['hi']}]  ({ci['draws_used']} draws, "
          f"{ci['draws_skipped']} skipped)  bar {BAR}  against authors: {against}")
    print(f"shared error (invalid, other state agrees Z valid): {results['shared_error_share']}   "
          f"false alarm (valid, flagged): {results['false_alarm_share']}")
    print(f"secondary AUC: union flag {results['auc_or_diff']}, moved knowledge {results['auc_moved']}")
    print(f"{'group':<26}{'rows':>6}{'inval':>7}{'AUC':>8}{'AUCor':>8}{'AUCmv':>8}{'shared':>8}{'falarm':>8}")
    for f, d in by.items():
        for g, s in d.items():
            cells = ["-" if s[k] is None else f"{s[k]:.3f}" for k in
                     ("auc", "auc_or_diff", "auc_moved", "shared_error_share", "false_alarm_share")]
            print(f"{f[:3]}:{g:<22}{s['rows']:>6}{s['invalid']:>7}" + "".join(f"{c:>8}" for c in cells))
    print(f"consistency: {'ok' if not fails else f'{len(fails)} FAILURES'}   -> {out / 'results.json'}")
    return 2 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
