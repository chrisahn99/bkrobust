"""E29 analysis: registered statistics, secondary cells, results.json. Works on full and smoke output.

    PYTHONPATH=src .venv/bin/python experiments/rebuttal/E29_radius_vs_degree_sep_undefined/analyse.py --out results/rebuttal/E29

The within-state helpers (pair_counts, tau_from_counts, the network-draw bootstrap) are copied from
experiments/rebuttal/E16_radius_where_separation_undefined/run.py at commit 8469bbf, unchanged except that
they take the endpoint as an argument, so the E16 within-state point (0.0827 on 757 units, 48 states) is
reproduced here as a check before the paired cell is read.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import beta

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "src"))

from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.robustness.real_analyse import tau_for_stratum  # noqa: E402

ENDPOINT = "AUC_frac"
SEED = 0
N_BOOT = 10_000
E16_WITHIN = {"point": 0.08271901100062913, "n_units": 757, "n_states": 48}
E16_POOLED = {"tau_b": 0.1337765673492525, "n": 1033}


# ------------------------------------------------------------------ copied from E16 run.py (8469bbf)

def pair_counts(xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
    i, j = np.triu_indices(len(xs), 1)
    dx, dy = xs[i] - xs[j], ys[i] - ys[j]
    keep = ~((dx == 0) & (dy == 0))
    dx, dy = dx[keep], dy[keep]
    tx = np.sum(dx == 0)
    ty = np.sum((dx != 0) & (dy == 0))
    both = (dx != 0) & (dy != 0)
    c = np.sum(both & (dx * dy > 0))
    d = np.sum(both & (dx * dy < 0))
    return np.array([c, d, tx, ty], dtype=np.int64)


def tau_from_counts(cnt: np.ndarray) -> float | None:
    c, d, tx, ty = (int(v) for v in cnt)
    denom = math.sqrt((c + d + tx) * (c + d + ty))
    return (c - d) / denom if denom > 0 else None


def network_count_matrix(states, col, endpoint):
    networks = sorted({k[0] for k in states})
    m = np.zeros((len(networks), 4), dtype=np.int64)
    pos = {n: i for i, n in enumerate(networks)}
    for k, rows in states.items():
        xs = np.array([float(r[col]) for r in rows])
        ys = np.array([float(r[endpoint]) for r in rows])
        m[pos[k[0]]] += pair_counts(xs, ys)
    return networks, m


def within_bootstrap(states, cols, endpoint, n_boot, seed):
    mats = []
    networks = None
    for col in cols:
        networks, m = network_count_matrix(states, col, endpoint)
        mats.append(m)
    rng = np.random.default_rng(seed)
    out = np.full((n_boot, len(cols)), np.nan)
    for b in range(n_boot):
        idx = rng.integers(0, len(networks), size=len(networks))
        for c, m in enumerate(mats):
            t = tau_from_counts(m[idx].sum(axis=0))
            if t is not None:
                out[b, c] = t
    return out


def percentile_ci(vals: np.ndarray):
    valid = vals[~np.isnan(vals)]
    if len(valid) < 100:
        return None, None, int(len(vals) - len(valid))
    lo, hi = np.percentile(valid, [2.5, 97.5])
    return float(lo), float(hi), int(len(vals) - len(valid))


# ------------------------------------------------------------------ E29 cells

def load_rows(out: Path) -> tuple[list[dict], list[dict], list[dict]]:
    rows = {}
    for f in sorted(glob.glob(str(out / "units*.jsonl"))):
        for line in open(f, encoding="utf-8"):
            if line.strip():
                r = json.loads(line)
                rows[r["key"]] = r  # last line for a key wins (resumed runs)
    fails = [json.loads(l) for f in sorted(glob.glob(str(out / "consistency_failures*.jsonl")))
             for l in open(f, encoding="utf-8") if l.strip()]
    errs = [json.loads(l) for f in sorted(glob.glob(str(out / "errors*.jsonl")))
            for l in open(f, encoding="utf-8") if l.strip()]
    return list(rows.values()), fails, errs


def tie_ceiling(x: np.ndarray, y: np.ndarray) -> float | None:
    """Largest tau_b attainable with these two marginals: pair the sorted values (comonotone)."""
    return tau_from_counts(pair_counts(np.sort(x), np.sort(y)))


def share_with_ci(rows, flag, n_boot, seed):
    k = sum(1 for r in rows if r[flag])
    n = len(rows)
    out = {"k": k, "n": n, "share": k / n if n else None}
    if n:
        out["clopper_pearson"] = [float(beta.ppf(0.025, k, n - k + 1)) if k > 0 else 0.0,
                                  float(beta.ppf(0.975, k + 1, n - k)) if k < n else 1.0]
        by_net = defaultdict(list)
        for r in rows:
            by_net[r["network"]].append(bool(r[flag]))
        nets = sorted(by_net)
        kk = np.array([sum(by_net[s]) for s in nets])
        nn = np.array([len(by_net[s]) for s in nets])
        rng = np.random.default_rng(seed)
        vals = []
        for _ in range(n_boot):
            idx = rng.integers(0, len(nets), size=len(nets))
            vals.append(kk[idx].sum() / nn[idx].sum())
        lo, hi = np.percentile(vals, [2.5, 97.5])
        out["network_cluster_ci"] = [float(lo), float(hi)]
    return out


def e16_states(rows: list[dict]) -> dict[tuple, list[dict]]:
    """E16's within-state population: (network, condition) states where the radius varies over all
    scorable real-name units, restricted to the separation-undefined units, kept when 2 or more remain."""
    radii = defaultdict(set)
    with (REPO / "results/axis_robustness_llm/analysis_units.csv").open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("naming") != "real" or r["radius"] == "" or r[ENDPOINT] == "":
                continue
            rad = float(r["radius"])
            if rad == UNREACHED:
                continue
            radii[(r["network"], r["condition"])].add(rad)
    varying = {k for k, v in radii.items() if len(v) > 1}
    states = defaultdict(list)
    for r in rows:
        k = (r["network"], r["condition"])
        if k in varying:
            states[k].append(r)
    return {k: v for k, v in states.items() if len(v) >= 2}


def fmt(p, lo, hi):
    if p is None:
        return "n/a"
    ci = f"[{lo:+.3f}, {hi:+.3f}]" if lo is not None and hi is not None else "[n/a]"
    return f"{p:+.3f} {ci}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(REPO / "results/rebuttal/E29"))
    ap.add_argument("--n-boot", type=int, default=None, help="default 10000 (200 when the run is a smoke)")
    args = ap.parse_args(argv)
    out = Path(args.out)
    all_rows, fails, errs = load_rows(out)
    if not all_rows:
        print(f"no rows in {out}", file=sys.stderr)
        return 1
    smoke = any(r.get("smoke") for r in all_rows)
    n_boot = args.n_boot or (200 if smoke else N_BOOT)
    ok = [r for r in all_rows if r.get("status") == "ok"]
    status = Counter(r.get("status") for r in all_rows)
    res: dict = {"experiment": "E29", "smoke": smoke, "n_boot": n_boot, "seed": SEED,
                 "n_rows": len(all_rows), "n_ok": len(ok), "status_counts": dict(status),
                 "consistency_failures": len(fails), "consistency_failure_kinds": dict(Counter(f["kind"] for f in fails)),
                 "errors": len(errs)}
    lines = [f"E29 analysis ({'SMOKE' if smoke else 'full'}), {len(ok)} of {len(all_rows)} units ok, "
             f"{len(fails)} consistency failures, {len(errs)} item errors, {n_boot} resamples"]

    # registered 1: tau_b(radius, undirected degree of x), pooled, network-cluster bootstrap
    t1 = tau_for_stratum(ok, "radius", "deg_u_x", n_boot=n_boot, seed=SEED)
    res["registered_tau_radius_degree"] = t1
    lines.append(f"REGISTERED 1  tau_b(radius, undirected degree of X) = "
                 f"{fmt(t1['tau_b'], t1['ci_lo_2p5'], t1['ci_hi_97p5'])}, n = {t1['n']}, "
                 f"{t1['n_networks']} networks, LOO [{t1['loo_tau_min']}, {t1['loo_tau_max']}]")

    # registered 2: share of radius-1 units whose only analyst orientation at x is the one retracted
    r1 = [r for r in ok if r["radius"] == 1]
    s = share_with_ci(r1, "only_kx_retracted", n_boot, SEED)
    res["registered_share_only_kx_retracted"] = s
    if s["n"]:
        lines.append(f"REGISTERED 2  radius-1 units whose only orientation at X is the retracted one: "
                     f"{s['k']} / {s['n']} = {s['share']:.3f}, Clopper-Pearson "
                     f"[{s['clopper_pearson'][0]:.3f}, {s['clopper_pearson'][1]:.3f}], network-cluster "
                     f"[{s['network_cluster_ci'][0]:.3f}, {s['network_cluster_ci'][1]:.3f}]")
    s_any = share_with_ci(r1, "some_d1_retracts_at_x", n_boot, SEED)
    res["share_some_d1_retracts_at_x"] = s_any
    if s_any["n"]:
        lines.append(f"  beside it: radius-1 units with some distance-1 witness retracting an orientation at X: "
                     f"{s_any['k']} / {s_any['n']} = {s_any['share']:.3f}")

    # sensitivity: the registered statistics without units whose GAC-recompute Z has a member in X's component
    clean = [r for r in ok if r["z_in_component"] == 0]
    t1c = tau_for_stratum(clean, "radius", "deg_u_x", n_boot=n_boot, seed=SEED)
    sc = share_with_ci([r for r in clean if r["radius"] == 1], "only_kx_retracted", n_boot, SEED)
    res["sensitivity_z_outside_component"] = {"n_units": len(clean), "n_dropped": len(ok) - len(clean),
                                              "tau_radius_degree": t1c, "share_only_kx_retracted": sc}
    lines.append(f"  sensitivity, {len(ok) - len(clean)} units with a Z member in X's component dropped: "
                 f"tau_b {fmt(t1c['tau_b'], t1c['ci_lo_2p5'], t1c['ci_hi_97p5'])}, share "
                 f"{sc['k']} / {sc['n']}")

    # secondary (a): tie ceiling of tau_b(radius, AUC)
    x = np.array([float(r["radius"]) for r in ok])
    y = np.array([float(r[ENDPOINT]) for r in ok])
    if len(ok) >= 2:
        obs = tau_from_counts(pair_counts(x, y))
        ceil = tie_ceiling(x, y)
        res["tie_ceiling"] = {"observed_tau_b": obs, "ceiling_tau_b": ceil,
                              "share_of_ceiling": obs / ceil if obs is not None and ceil else None,
                              "e16_pooled_tau_b": E16_POOLED["tau_b"]}
        lines.append(f"(a) tau_b(radius, AUC) observed {obs:+.4f} (E16: {E16_POOLED['tau_b']:+.4f}), "
                     f"ceiling under these ties {ceil:+.4f}, share of ceiling "
                     f"{(obs / ceil if ceil else float('nan')):.3f}")

    # secondary (b): failure mechanism of the witness, by radius
    mech = defaultdict(Counter)
    for r in ok:
        mech[int(r["radius"])][r.get("witness_mechanism")] += 1
    d1_mech = Counter(w["mechanism"] for r in r1 for w in r["d1_witnesses"])
    res["mechanism_by_radius"] = {str(k): dict(v) for k, v in sorted(mech.items())}
    res["mechanism_all_d1_witnesses"] = dict(d1_mech)
    lines.append("(b) witness mechanism by radius: " + "; ".join(
        f"r={k}: " + ", ".join(f"{m} {c}" for m, c in sorted(v.items())) for k, v in sorted(mech.items())))
    lines.append(f"    all distance-1 witnesses of radius-1 units: {dict(d1_mech)}")

    # secondary (c): within-state paired tau_b(radius, AUC) - tau_b(degree, AUC), E16's states
    states = e16_states(ok)
    c_res: dict = {"n_units": sum(len(v) for v in states.values()), "n_states": len(states),
                   "n_networks": len({k[0] for k in states})}
    if states:
        pts = [tau_from_counts(network_count_matrix(states, col, ENDPOINT)[1].sum(axis=0))
               for col in ("radius", "deg_u_x")]
        boot = within_bootstrap(states, ["radius", "deg_u_x"], ENDPOINT, n_boot, SEED)
        delta = boot[:, 0] - boot[:, 1]
        lo, hi, nan = percentile_ci(delta)
        c_res.update(tau_radius=pts[0], tau_degree=pts[1],
                     point=None if None in pts else pts[0] - pts[1], ci_lo_2p5=lo, ci_hi_97p5=hi,
                     n_boot_nan=nan, ci_radius=list(percentile_ci(boot[:, 0])[:2]),
                     ci_degree=list(percentile_ci(boot[:, 1])[:2]),
                     n_states_degree_varies=sum(1 for v in states.values() if len({r["deg_u_x"] for r in v}) > 1),
                     e16_reproduced=(not smoke and c_res["n_units"] == E16_WITHIN["n_units"]
                                     and c_res["n_states"] == E16_WITHIN["n_states"]
                                     and pts[0] is not None and abs(pts[0] - E16_WITHIN["point"]) < 1e-9))
        lines.append(f"(c) within-state, {c_res['n_units']} units, {c_res['n_states']} states: tau_b(radius, AUC) "
                     f"{fmt(pts[0], *c_res['ci_radius'])}, tau_b(degree, AUC) {fmt(pts[1], *c_res['ci_degree'])}, "
                     f"paired radius minus degree {fmt(c_res['point'], lo, hi)}; E16 point reproduced: "
                     f"{c_res['e16_reproduced']}")
    res["within_paired_radius_minus_degree"] = c_res

    # secondary (d): pooled degree and component size against AUC and radius; cross-tabs
    for col, ep in (("deg_u_x", ENDPOINT), ("comp_size", "radius"), ("comp_size", ENDPOINT), ("n_k_x", "radius")):
        t = tau_for_stratum(ok, col, ep, n_boot=n_boot, seed=SEED)
        res[f"pooled_tau_{col}_vs_{ep}"] = t
        lines.append(f"(d) tau_b({col}, {ep}) = {fmt(t['tau_b'], t['ci_lo_2p5'], t['ci_hi_97p5'])}, n = {t['n']}")
    xt = defaultdict(Counter)
    for r in ok:
        xt[int(r["radius"])][int(r["deg_u_x"])] += 1
    res["crosstab_radius_by_degree"] = {str(k): {str(d): c for d, c in sorted(v.items())} for k, v in sorted(xt.items())}
    lines.append("    radius x undirected degree of X: " + "; ".join(
        f"r={k}: " + ", ".join(f"d{d}:{c}" for d, c in sorted(v.items())) for k, v in sorted(xt.items())))
    xk = defaultdict(Counter)
    for r in ok:
        xk[int(r["radius"])][int(r["n_k_x"])] += 1
    res["crosstab_radius_by_n_k_x"] = {str(k): {str(d): c for d, c in sorted(v.items())} for k, v in sorted(xk.items())}
    lines.append("    radius x analyst orientations at X: " + "; ".join(
        f"r={k}: " + ", ".join(f"k{d}:{c}" for d, c in sorted(v.items())) for k, v in sorted(xk.items())))
    res["radius_counts"] = dict(Counter(str(int(r["radius"])) for r in ok))
    res["median_t_s"] = float(np.median([r["t_s"] for r in ok])) if ok else None
    lines.append(f"radius counts {res['radius_counts']}; median seconds per unit {res['median_t_s']}")

    (out / "results.json").write_text(json.dumps(res, indent=1, default=float))
    (out / "analyse.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
