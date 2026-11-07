"""E9 analysis: AUC of parametric invalidity for "Z invalid at truth", and its gain over E6.

    python experiments/rebuttal/E9_parametric_bootstrap_validity/analyse.py --out results/rebuttal/E9
    python experiments/rebuttal/E9_parametric_bootstrap_validity/analyse.py --out results/rebuttal/E9_smoke --e6 results/rebuttal/E6_smoke

Reads reference.jsonl and every units_*.jsonl in --out, and E6's per_row.json in --e6
when it exists; writes results.json and per_row.json in --out. The registered
statistic and its interval are fixed in PREREG.md.
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
MIN_GAIN = 0.05
MIN_FINISHED = 160


def auc(scores, labels) -> float | None:
    """P(score of a positive > score of a negative), ties one half; None without both classes."""
    s = np.asarray(scores, float)
    lab = np.asarray(labels, int)
    pos, neg = s[lab == 1], s[lab == 0]
    if len(pos) == 0 or len(neg) == 0:
        return None
    gt = (pos[:, None] > neg[None, :]).sum()
    eq = (pos[:, None] == neg[None, :]).sum()
    return float((gt + 0.5 * eq) / (len(pos) * len(neg)))


def cluster_draws(rows: list[dict], rng: np.random.Generator) -> list[list[int]]:
    """Index lists of the two-stage cluster bootstrap: networks, then query pairs in each."""
    by_net: dict = collections.defaultdict(lambda: collections.defaultdict(list))
    for i, r in enumerate(rows):
        by_net[r["network"]][f"{r['x']}|{r['y']}"].append(i)
    nets = sorted(by_net)
    draws = []
    for _ in range(N_DRAWS):
        idx: list[int] = []
        for net in rng.choice(nets, size=len(nets), replace=True):
            pairs = sorted(by_net[net])
            for p in rng.choice(pairs, size=len(pairs), replace=True):
                idx += by_net[net][p]
        draws.append(idx)
    return draws


def interval(vals: list[float], skipped: int) -> dict:
    if not vals:
        return {"lo": None, "hi": None, "draws_used": 0, "draws_skipped": skipped}
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return {"lo": round(float(lo), 4), "hi": round(float(hi), 4), "draws_used": len(vals), "draws_skipped": skipped}


def r4(v):
    return None if v is None else round(v, 4)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/rebuttal/E9", help="the run's output folder")
    ap.add_argument("--e6", default="results/rebuttal/E6", help="E6's output folder (per_row.json)")
    args = ap.parse_args(argv)
    out = Path(args.out)
    ref_all = [json.loads(line) for line in (out / "reference.jsonl").open()]
    # rows whose state has no simulated truth are listed with their reason, never scored
    excluded = [r for r in ref_all if r.get("excluded")]
    ref = [r for r in ref_all if not r.get("excluded")]
    excluded_reasons = dict(collections.Counter(r["excluded"] for r in excluded))
    excluded_networks = sorted({r["network"] for r in excluded} - {r["network"] for r in ref})
    units: dict[str, dict] = {}
    for path in sorted(out.glob("units_*.jsonl")):
        for line in path.read_text().splitlines():
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue  # a line cut by a kill
            units.setdefault(rec["unit"], rec)
    config = json.loads((out / "config.json").read_text())

    sims = [u for u in units.values() if u["b"] >= 0]
    identity = [u for u in units.values() if u["b"] < 0]
    finished = [u for u in sims if u["status"] == "ok"]
    per_net = collections.Counter(u["network"] for u in finished)
    timeouts = collections.Counter(u["network"] for u in sims if u["status"] == "timeout")
    criteria = {u["criterion"] for u in units.values()}
    commits = {u["repo_commit"] for u in units.values()}

    agg: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    for u in finished:
        for r in u["rows"]:
            c = agg[r["rid"]]
            c["n"] += 1
            c["valid"] += r["zb_valid"]
            c["exists"] += r["zb_exists"]
            c["invalid_given_exists"] += r["zb_exists"] - r["zb_valid"]
            c["same"] += r["zb_same"]
            c["fix"] += r["zfix_valid_state"]
            c["ext_valid"] += r["zb_valid_ext"]
            c["ext_n"] += r["n_ext"]
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
                "n_sim": n,
                "validity": c["valid"] / n,
                "score": 1 - c["valid"] / n,
                "score_fixed_state": 1 - c["fix"] / n,
                "score_zchange": 1 - c["same"] / n,
                # no committed set anywhere: the most suspicious value, counted below
                "score_conditional": c["invalid_given_exists"] / c["exists"] if c["exists"] else 1.0,
                "score_ext": 1 - c["ext_valid"] / c["ext_n"] if c["ext_n"] else None,
                "never_committed": int(c["exists"] == 0),
                "errors": c["err"],
            }
        )
    if not rows:
        print("no finished simulations yet", file=sys.stderr)
        return 1
    labels = [r["invalid"] for r in rows]

    # E6 joined by row id; the paired difference uses the same cluster draws
    e6_path = Path(args.e6) / "per_row.json"
    e6 = {}
    e6_smoke = None
    if e6_path.exists():
        e6 = {r["rid"]: r["score"] for r in json.loads(e6_path.read_text())}
        e6_cfg = Path(args.e6) / "config.json"
        e6_smoke = json.loads(e6_cfg.read_text()).get("smoke") if e6_cfg.exists() else None
    for r in rows:
        r["score_e6"] = e6.get(r["rid"])
    paired = [i for i, r in enumerate(rows) if r["score_e6"] is not None]

    rng = np.random.default_rng(CI_SEED)
    draws = cluster_draws(rows, rng)
    head = auc([r["score"] for r in rows], labels)
    v9, vd, sk9, skd = [], [], 0, 0
    pset = set(paired)
    for idx in draws:
        a = auc([rows[i]["score"] for i in idx], [labels[i] for i in idx])
        if a is None:
            sk9 += 1
        else:
            v9.append(a)
        if paired:
            pi = [i for i in idx if i in pset]
            a9 = auc([rows[i]["score"] for i in pi], [labels[i] for i in pi])
            a6 = auc([rows[i]["score_e6"] for i in pi], [labels[i] for i in pi])
            if a9 is None or a6 is None:
                skd += 1
            else:
                vd.append(a9 - a6)
    auc9_p = auc([rows[i]["score"] for i in paired], [labels[i] for i in paired]) if paired else None
    auc6_p = auc([rows[i]["score_e6"] for i in paired], [labels[i] for i in paired]) if paired else None
    delta = None if auc9_p is None or auc6_p is None else auc9_p - auc6_p
    ci_d = interval(vd, skd) if paired else None

    inv = [r for r in rows if r["invalid"]]
    val = [r for r in rows if not r["invalid"]]
    secondary = {
        key: r4(auc([r[key] for r in rows if r[key] is not None], [r["invalid"] for r in rows if r[key] is not None]))
        for key in ("score_fixed_state", "score_zchange", "score_conditional", "score_ext")
    }
    by: dict = {}
    for field in ("network", "arm"):
        by[field] = {}
        for g in sorted({r[field] for r in rows}):
            sub = [r for r in rows if r[field] == g]
            si = [r["validity"] for r in sub if r["invalid"]]
            sv = [r["validity"] for r in sub if not r["invalid"]]
            by[field][g] = {
                "rows": len(sub),
                "invalid": len(si),
                "auc": r4(auc([r["score"] for r in sub], [r["invalid"] for r in sub])),
                "mean_validity_invalid": round(float(np.mean(si)), 4) if si else None,
                "mean_validity_valid": round(float(np.mean(sv)), 4) if sv else None,
            }
    e6_final = bool(paired) and e6_smoke is False and not config["smoke"]
    results = {
        "experiment": "E9",
        "statistic": "AUC of (1 - share of simulations with recomputed Z_b GAC-valid in the fitted truth) for Z invalid at truth; Delta = AUC(E9) - AUC(E6) on the paired rows",
        "auc": r4(head),
        "ci95_cluster_network_then_pair": interval(v9, sk9),
        "paired_rows": len(paired),
        "auc_e9_paired": r4(auc9_p),
        "auc_e6_paired": r4(auc6_p),
        "delta": r4(delta),
        "delta_ci95": ci_d,
        "against_authors": None if delta is None else bool(delta < MIN_GAIN),
        "provisional": not e6_final,
        "e6_folder": str(Path(args.e6)),
        "e6_smoke": e6_smoke,
        "rows": len(rows),
        "rows_reference": len(ref_all),
        "rows_simulated_reference": len(ref),
        "rows_excluded": len(excluded),
        "rows_excluded_invalid": sum(1 - r["z_valid_at_truth"] for r in excluded),
        "excluded_reasons": excluded_reasons,
        "networks_excluded": excluded_networks,
        "rows_invalid": len(inv),
        "rows_valid": len(val),
        "valid_but_wrong_share": round(sum(r["validity"] >= 0.95 for r in inv) / len(inv), 4) if inv else None,
        "rows_never_committed": sum(r["never_committed"] for r in rows),
        "secondary_auc": secondary,
        "by": by,
        "simulations_finished": dict(per_net),
        "simulations_timeout": dict(timeouts),
        "networks_under_160_simulations": sorted(
            n for n in config["networks"] if n not in excluded_networks and per_net.get(n, 0) < MIN_FINISHED
        ),
        "identity_reproduces_reference": {u["network"]: u.get("is_reference_cpdag") for u in identity},
        "rows_with_check_errors": sum(1 for r in rows if r["errors"]),
        "criterion": sorted(criteria),
        "repo_commit": sorted(commits),
        "smoke": config["smoke"],
        "S_planned": config["S"],
    }
    (out / "results.json").write_text(json.dumps(results, indent=1))
    (out / "per_row.json").write_text(json.dumps(rows, indent=1))

    ci = results["ci95_cluster_network_then_pair"]
    print(f"E9  {'SMOKE  ' if config['smoke'] else ''}criterion {','.join(sorted(criteria))}  commit {','.join(c[:10] for c in sorted(commits))}")
    print(f"rows {len(rows)} of {len(ref_all)} (invalid at truth {len(inv)}, valid {len(val)}); simulations finished {dict(per_net)}, timeouts {dict(timeouts)}")
    if excluded:
        print(f"excluded {len(excluded)} rows ({results['rows_excluded_invalid']} invalid at truth), networks {excluded_networks}:")
        for why, c in sorted(excluded_reasons.items()):
            print(f"  {c}: {why}")
    print(f"AUC(E9) {results['auc']}  95% CI [{ci['lo']}, {ci['hi']}]  ({ci['draws_used']} draws, {ci['draws_skipped']} skipped)")
    if paired:
        print(f"paired rows {len(paired)}: AUC(E9) {results['auc_e9_paired']}  AUC(E6) {results['auc_e6_paired']}  Delta {results['delta']}  95% CI [{ci_d['lo']}, {ci_d['hi']}]")
        print(f"against authors (Delta < {MIN_GAIN}): {results['against_authors']}{'   PROVISIONAL (E6 or E9 is a smoke run)' if results['provisional'] else ''}")
    else:
        print(f"no E6 rows at {e6_path}: Delta not computed, result PROVISIONAL")
    print(f"valid but wrong (invalid, v >= 0.95): {results['valid_but_wrong_share']}   rows never committed: {results['rows_never_committed']}")
    print(f"secondary AUC: {secondary}")
    print(f"{'group':<22}{'rows':>6}{'inval':>7}{'AUC':>8}{'v|inv':>8}{'v|val':>8}")
    for field in ("network", "arm"):
        for g, d in by[field].items():
            a = "-" if d["auc"] is None else f"{d['auc']:.3f}"
            si = "-" if d["mean_validity_invalid"] is None else f"{d['mean_validity_invalid']:.3f}"
            sv = "-" if d["mean_validity_valid"] is None else f"{d['mean_validity_valid']:.3f}"
            print(f"{field[:3]}:{g:<18}{d['rows']:>6}{d['invalid']:>7}{a:>8}{si:>8}{sv:>8}")
    print(f"identity unit = reference CPDAG: {results['identity_reproduces_reference']}")
    print(f"-> {out / 'results.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
