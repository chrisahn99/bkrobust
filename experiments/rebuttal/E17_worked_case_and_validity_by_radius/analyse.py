"""E17 analysis: the worked case, and share invalid at truth by radius bin, network-clustered.

    .venv/bin/python experiments/rebuttal/E17_worked_case_and_validity_by_radius/analyse.py --out results/rebuttal/E17

Reads <out>/rows.*.jsonl (and the released records for the secondary table), prints a
plain-text table and writes <out>/results.json. Registered statistic: PREREG.md.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
BINS = ("1", "2", ">=3", "unreached")
ARMS = ("D_LLM", "D_RANDOM_MATCHED")


def bin_of(radius: int) -> str:
    if radius == -1:
        return "unreached"
    return "1" if radius == 1 else "2" if radius == 2 else ">=3"


def load_rows(out: Path) -> list[dict]:
    rows: dict[str, dict] = {}
    for f in sorted(glob.glob(str(out / "rows.*.jsonl"))):
        for line in open(f):
            if line.strip():
                r = json.loads(line)
                rows[r["key"]] = r  # a resumed shard can only repeat a key with the same unit
    return list(rows.values())


def counts_by_network(units: list[tuple[str, str, bool]]) -> dict[str, np.ndarray]:
    """units = (network, bin, invalid) -> per network an array [bin][n, invalid]."""
    out: dict[str, np.ndarray] = {}
    for net, b, inv in units:
        a = out.setdefault(net, np.zeros((len(BINS), 2)))
        a[BINS.index(b), 0] += 1
        a[BINS.index(b), 1] += inv
    return out


def gap(a: np.ndarray) -> float | None:
    """share invalid at radius 1 minus share invalid at radius >= 2 or unreached."""
    n1, i1 = a[..., 0, 0], a[..., 0, 1]
    nl, il = a[..., 1:, 0].sum(-1), a[..., 1:, 1].sum(-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where((n1 > 0) & (nl > 0), i1 / n1 - il / nl, np.nan)


def arm_table(per_net: dict[str, np.ndarray], weights: np.ndarray | None, nets: list[str]) -> dict:
    tot = sum(per_net[n] for n in nets) if nets else np.zeros((len(BINS), 2))
    stack = np.stack([per_net[n] for n in nets]) if nets else None
    res = {"n_networks": len(nets), "bins": {}}
    for i, b in enumerate(BINS):
        n, inv = int(tot[i, 0]), int(tot[i, 1])
        cell = {"n": n, "invalid": inv, "share": inv / n if n else None, "ci95": None}
        if weights is not None and n:
            bn = weights @ stack[:, i, 0]
            bi = weights @ stack[:, i, 1]
            ok = bn > 0
            if ok.any():
                cell["ci95"] = [float(v) for v in np.percentile(bi[ok] / bn[ok], [2.5, 97.5])]
        res["bins"][b] = cell
    g = gap(tot)
    res["gap"] = None if np.isnan(g) else float(g)
    if weights is not None and nets:
        bg = gap(np.einsum("bn,nkc->bkc", weights, stack))
        ok = ~np.isnan(bg)
        res["gap_ci95"] = [float(v) for v in np.percentile(bg[ok], [2.5, 97.5])] if ok.any() else None
        res["gap_resamples_dropped"] = int((~ok).sum())
    return res


def compare(per_net: dict[str, dict[str, np.ndarray]], arms: tuple[str, str], n_res: int,
            seed: int) -> dict:
    """Both arms' tables and Delta = gap(arm0) - gap(arm1), resampling the shared networks."""
    common = sorted(set(per_net[arms[0]]) & set(per_net[arms[1]]))
    rng = np.random.default_rng(seed)
    w = None
    if common:
        idx = rng.integers(0, len(common), size=(n_res, len(common)))
        w = np.stack([np.bincount(r, minlength=len(common)) for r in idx]).astype(float)
    out = {"networks_common": common, "arms": {}}
    for a in arms:
        out["arms"][a] = arm_table(per_net[a], w, common)
        out["arms"][a]["all_networks"] = arm_table(per_net[a], None, sorted(per_net[a]))
    g0, g1 = out["arms"][arms[0]]["gap"], out["arms"][arms[1]]["gap"]
    out["delta"] = None if g0 is None or g1 is None else g0 - g1
    out["delta_ci95"], out["delta_resamples_dropped"] = None, None
    if w is not None:
        s = [np.stack([per_net[a][n] for n in common]) for a in arms]
        d = gap(np.einsum("bn,nkc->bkc", w, s[0])) - gap(np.einsum("bn,nkc->bkc", w, s[1]))
        ok = ~np.isnan(d)
        out["delta_ci95"] = [float(v) for v in np.percentile(d[ok], [2.5, 97.5])] if ok.any() else None
        out["delta_resamples_dropped"] = int((~ok).sum())
    return out


def worked_case(rows: list[dict], order: list[str]) -> dict:
    for net in order:
        ok = sorted((r for r in rows if r["kind"] == "worked" and r["network"] == net
                     and r.get("status") == "ok" and r.get("radius_exact")),
                    key=lambda r: (r["x"], r["y"]))
        one = [r for r in ok if r["radius"] == 1 and r.get("claim_witness")]
        far = [r for r in ok if r["radius"] == -1] or sorted(
            (r for r in ok if r["radius"] >= 2), key=lambda r: -r["radius"])
        if one and far:
            keep = ("x", "y", "k", "z", "radius", "radius_method", "claim_radius",
                    "claim_witness", "claim_witness_wrong_at_truth", "wrong_claims",
                    "n_wrong_closure", "valid_gac_truth", "valid_backdoor_truth", "release",
                    "checks", "notes")
            return {"network": net, "found": True,
                    "radius_one": {c: one[0].get(c) for c in keep},
                    "large_radius": {c: far[0].get(c) for c in keep},
                    "n_queries_ok": len(ok), "n_radius_one": len(one),
                    "n_unreached": sum(r["radius"] == -1 for r in ok),
                    "share_invalid_at_radius_one": sum(not r["valid_gac_truth"] for r in one) / len(one)}
    return {"found": False}


def release_units(path: Path) -> dict[str, list]:
    out = {"D_LLM": [], "D_RAND": []}
    if not path.exists():
        return out
    with path.open(newline="") as fh:
        for r in csv.DictReader(fh):
            arm = r["arm"]
            if arm not in out or (arm == "D_LLM" and r["rep"] != "0"):
                continue
            if r["status"] != "ok" or r["r_hop"] in ("", "0") or r["hop_exact"] != "1":
                continue
            out[arm].append((r["network"], bin_of(int(r["r_hop"])), r["z_valid_at_truth"] == "0"))
    return out


def fmt(cell: dict) -> str:
    if not cell["n"]:
        return f"{'-':>22}"
    ci = cell["ci95"]
    ci_s = f" [{ci[0]:.2f},{ci[1]:.2f}]" if ci else ""
    return f"{cell['invalid']:>4}/{cell['n']:<5}{cell['share']:.2f}{ci_s}".rjust(22)


def print_comparison(title: str, c: dict, arms: tuple[str, str]) -> None:
    print(f"\n{title}  (networks in both arms: {len(c['networks_common'])})")
    print(f"{'arm':<18}" + "".join(f"{'r=' + b:>22}" for b in BINS) + f"{'gap':>22}")
    for a in arms:
        t = c["arms"][a]
        g, ci = t["gap"], t.get("gap_ci95")
        gs = "-" if g is None else f"{g:+.3f}" + (f" [{ci[0]:+.2f},{ci[1]:+.2f}]" if ci else "")
        print(f"{a:<18}" + "".join(fmt(t["bins"][b]) for b in BINS) + f"{gs:>22}")
    d, ci = c["delta"], c["delta_ci95"]
    print("Delta gap = " + ("undefined" if d is None else f"{d:+.3f}")
          + (f"  95% CI [{ci[0]:+.3f}, {ci[1]:+.3f}]" if ci else "")
          + f"  (resamples dropped: {c['delta_resamples_dropped']})")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/rebuttal/E17")
    ap.add_argument("--release-records", default="experiments/rebuttal/_data/per_query_records.csv")
    ap.add_argument("--worked-networks", default="Acid_1996,Polzer_2012")
    ap.add_argument("--n-resamples", type=int, default=10000)
    ap.add_argument("--boot-seed", type=int, default=20261106)
    args = ap.parse_args(argv)
    rp = lambda s: Path(s) if Path(s).is_absolute() else ROOT / s  # noqa: E731
    out = rp(args.out)
    rows = load_rows(out)
    if not rows:
        print(f"no rows under {out}", file=sys.stderr)
        return 1

    res: dict = {"experiment": "E17", "criterion": sorted({r.get("criterion") for r in rows}),
                 "commits": sorted({r.get("commit") for r in rows}),
                 "smoke": any(r.get("smoke") for r in rows), "n_rows": len(rows)}
    res["failed_checks"] = [{"key": r["key"], "checks": r["checks"]} for r in rows if r.get("checks")]
    draws = [r for r in rows if r["kind"] == "draw"]
    res["draws"] = {"n": len(draws), "exhausted": sum(bool(r.get("exhausted")) for r in draws),
                    "not_reproduced": sum("draw_not_reproduced" in (r.get("checks") or []) for r in draws)}

    res["worked_case"] = worked_case(rows, [n for n in args.worked_networks.split(",") if n])

    table = [r for r in rows if r["kind"] == "table"]
    denom: dict = {}
    per_net: dict = {"gac": {a: {} for a in ARMS}, "backdoor": {a: {} for a in ARMS}}
    for a in ARMS:
        ra = [r for r in table if r["arm"] == a]
        d = denom.setdefault(a, {})
        for r in ra:
            s = r.get("status")
            if s == "ok" and r["radius"] == 0:
                s = "degenerate_radius_0"
            elif s == "ok" and not r.get("radius_exact"):
                s = "radius_not_exact"
            d[s] = d.get(s, 0) + 1
        inf = [r for r in ra if r.get("status") == "ok" and r["radius"] != 0 and r.get("radius_exact")]
        d["informative"] = len(inf)
        d["share_at_radius_1"] = (sum(r["radius"] == 1 for r in inf) / len(inf)) if inf else None
        for crit, col in (("gac", "valid_gac_truth"), ("backdoor", "valid_backdoor_truth")):
            per_net[crit][a] = counts_by_network([(r["network"], bin_of(r["radius"]), not r[col])
                                                  for r in inf])
    res["denominators"] = denom
    res["registered"] = compare(per_net["gac"], ARMS, args.n_resamples, args.boot_seed)
    res["secondary_backdoor"] = compare(per_net["backdoor"], ARMS, args.n_resamples, args.boot_seed)
    rel = release_units(rp(args.release_records))
    rel_net = {a: counts_by_network(u) for a, u in rel.items()}
    res["secondary_release_records"] = compare(rel_net, ("D_LLM", "D_RAND"), args.n_resamples,
                                               args.boot_seed)

    wc = res["worked_case"]
    print(f"E17  rows={len(rows)}  smoke={res['smoke']}  criterion={res['criterion']}  "
          f"failed checks={len(res['failed_checks'])}  draws not reproduced="
          f"{res['draws']['not_reproduced']}/{res['draws']['n']}")
    if wc["found"]:
        print(f"\nWorked case: {wc['network']}  K={wc['radius_one']['k']}  "
              f"wrong claims at truth={wc['radius_one']['wrong_claims']}")
        for lab in ("radius_one", "large_radius"):
            q = wc[lab]
            print(f"  {q['x']} -> {q['y']}: Z={q['z']} radius={q['radius']} "
                  f"claim radius={q['claim_radius']} witness={q['claim_witness']} "
                  f"valid at truth (GAC/back-door)={q['valid_gac_truth']}/{q['valid_backdoor_truth']}")
        print(f"  {wc['n_radius_one']} of {wc['n_queries_ok']} queries at radius 1, "
              f"{wc['n_unreached']} unreached; share invalid at radius 1 = "
              f"{wc['share_invalid_at_radius_one']:.2f}")
    else:
        print("\nWorked case: none found in the worked-case networks")
    print("\nDenominators:", json.dumps(denom))
    print_comparison("REGISTERED: share GAC-invalid at truth by GAC radius", res["registered"], ARMS)
    print_comparison("secondary: back-door validity at truth", res["secondary_backdoor"], ARMS)
    print_comparison("secondary: released records (r_hop, z_valid_at_truth)",
                     res["secondary_release_records"], ("D_LLM", "D_RAND"))
    (out / "results.json").write_text(json.dumps(res, indent=2, default=str))
    print(f"\nwrote {out / 'results.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
