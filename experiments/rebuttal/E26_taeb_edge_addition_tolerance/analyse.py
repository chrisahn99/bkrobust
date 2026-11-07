"""E26 analysis: the registered statistic of PREREG.md from results/rebuttal/E26/rows*.jsonl.

    .venv/bin/python experiments/rebuttal/E26_taeb_edge_addition_tolerance/analyse.py --out results/rebuttal/E26

Writes <out>/results.json and <out>/analyse.txt. Exit 0 always once the rows parse; the verdict
is a field, not an exit code.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import beta, kendalltau

INF = 10 ** 6
CAP = 4  # tolerance above the search depth (3) is the category "4 or more"
N_BOOT = 2000
BOOT_SEED = 20261107
E10_ARMS = {"D_LLM", "D_LLM_INSTR", "D_LLM_SRC", "D_SCRAMBLED", "D_LLM_72B", "D_LLM_72B_INSTR",
            "D_SCRAMBLED_72B", "D_LLM_32B", "D_LLM_SRC_70B", "D_LLM_GEMMA_27B", "D_RAND"}
DIES_TAU, DIES_LOWER = 0.8, 0.7


def load(out: Path) -> list[dict]:
    rows: dict[str, dict] = {}
    for f in sorted(out.glob("rows*.jsonl")):
        for line in f.open():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            if r.get("status") == "error" and r["key"] in rows:
                continue
            rows[r["key"]] = r
    return list(rows.values())


def cp(k: int, n: int) -> list:
    if n == 0:
        return [None, None]
    lo = 0.0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(0.975, k + 1, n - k))
    return [round(lo, 4), round(hi, 4)]


def radius_value(r: dict):
    """Released radius: an int, INF for unreached (-1), None when censored (hop_exact 0)."""
    if r["released_hop_exact"] != "1" or r["released_r_hop"] == "":
        return None
    v = int(r["released_r_hop"])
    return INF if v == -1 else v


def tol_value(t: dict):
    """Tolerance: exact value, CAP when the full depth-3 search found nothing, None if stopped by time."""
    if t.get("exact"):
        return t["t"]
    if t.get("stop") == "depth" or (t.get("lb") or 0) > 3:
        return CAP
    return None


def pair_shares(r: np.ndarray, t: np.ndarray) -> dict:
    """Over pairs the radius orders strictly: share ordered the same way by the tolerance,
    share reversed, share tied by the tolerance."""
    n = len(r)
    same = rev = tied = untied_r = 0
    for i in range(0, n, 512):
        dr = np.sign(r[i:i + 512, None] - r[None, :])
        dt = np.sign(t[i:i + 512, None] - t[None, :])
        idx = np.arange(i, min(i + 512, n))[:, None] < np.arange(n)[None, :]
        m = idx & (dr != 0)
        untied_r += int(m.sum())
        same += int((m & (dr == dt)).sum())
        rev += int((m & (dt == -dr)).sum())
        tied += int((m & (dt == 0)).sum())
    if untied_r == 0:
        return {"pairs_ordered_by_radius": 0}
    return {"pairs_ordered_by_radius": untied_r, "share_same": round(same / untied_r, 4),
            "share_reversed": round(rev / untied_r, 4), "share_tied_by_tolerance": round(tied / untied_r, 4),
            "share_not_ordered_alike": round(1 - same / untied_r, 4)}


def tau_b(r, t) -> float | None:
    if len(set(r)) < 2 or len(set(t)) < 2:
        return None
    v = kendalltau(r, t, variant="b").statistic
    return None if v != v else float(v)


def boot_tau(nets: list[str], r: np.ndarray, t: np.ndarray) -> list:
    by = defaultdict(list)
    for i, n in enumerate(nets):
        by[n].append(i)
    keys = sorted(by)
    rng = np.random.default_rng(BOOT_SEED)
    vals = []
    for _ in range(N_BOOT):
        pick = rng.choice(len(keys), size=len(keys), replace=True)
        idx = np.concatenate([by[keys[j]] for j in pick])
        v = tau_b(r[idx], t[idx])
        if v is not None:
            vals.append(v)
    if not vals:
        return [None, None, 0]
    return [round(float(np.percentile(vals, 2.5)), 4), round(float(np.percentile(vals, 97.5)), 4), len(vals)]


def auc(score: np.ndarray, label: np.ndarray) -> float | None:
    """Probability that an invalid row (label 1) scores higher than a valid one, ties half."""
    pos, neg = score[label == 1], score[label == 0]
    if len(pos) == 0 or len(neg) == 0:
        return None
    gt = (pos[:, None] > neg[None, :]).sum()
    eq = (pos[:, None] == neg[None, :]).sum()
    return round(float((gt + 0.5 * eq) / (len(pos) * len(neg))), 4)


def cat_r(v):
    if v is None:
        return "censored"
    if v == INF:
        return "inf"
    return str(v) if v < CAP else "4+"


def cat_t(v):
    if v is None:
        return "undecided"
    return str(v) if v < CAP else "4+ (none found to depth 3)"


def block(rows: list[dict], label: str) -> dict:
    ok = [r for r in rows if r.get("status") == "ok"]
    out: dict = {"population": label, "n_rows_ok": len(ok)}
    t_vals = [tol_value(r["tol_order"]) for r in ok]
    decided = [v for v in t_vals if v is not None]
    k1 = sum(1 for v in decided if v == 1)
    out["tolerance_order"] = {
        "n_decided": len(decided), "n_undecided_time": len(t_vals) - len(decided),
        "distribution": dict(sorted(Counter(cat_t(v) for v in t_vals).items())),
        "share_t1": round(k1 / len(decided), 4) if decided else None, "share_t1_cp95": cp(k1, len(decided)),
        "share_t1_worst_case_undecided_not_1": round(k1 / len(t_vals), 4) if t_vals else None,
    }
    a_vals = [tol_value(r["tol_acyclic"]) for r in ok]
    a_dec = [v for v in a_vals if v is not None]
    ka = sum(1 for v in a_dec if v == 1)
    out["tolerance_acyclic"] = {"n_decided": len(a_dec), "distribution": dict(sorted(Counter(cat_t(v) for v in a_vals).items())),
                                "share_t1": round(ka / len(a_dec), 4) if a_dec else None, "share_t1_cp95": cp(ka, len(a_dec))}

    pairs = [(radius_value(r), tol_value(r["tol_order"]), r["network"], int(r["z_valid_at_truth"])) for r in ok]
    both = [p for p in pairs if p[0] is not None and p[1] is not None]
    out["n_in_tau"] = len(both)
    out["n_radius_censored"] = sum(1 for p in pairs if p[0] is None)
    if len(both) >= 2:
        rr = np.array([p[0] for p in both], dtype=float)
        tt = np.array([p[1] for p in both], dtype=float)
        nets = [p[2] for p in both]
        tb = tau_b(rr, tt)
        out["tau_b"] = None if tb is None else round(tb, 4)
        out["tau_b_network_bootstrap95"] = boot_tau(nets, rr, tt) if tb is not None else None
        out["n_networks"] = len(set(nets))
        out["pairs"] = pair_shares(rr, tt)
        rc = np.minimum(rr, CAP)
        out["query_level"] = {
            "share_equal_after_cap4": round(float(np.mean(tt == rc)), 4),
            "share_tolerance_below_radius": round(float(np.mean(tt < rc)), 4),
            "share_tolerance_above_radius": round(float(np.mean(tt > rc)), 4),
        }
        cross = Counter((cat_r(p[0] if p[0] is not None else None), cat_t(p[1])) for p in both)
        out["crosstab_radius_x_tolerance"] = {f"r={a} | t={b}": c for (a, b), c in sorted(cross.items())}
        inval = np.array([1 - p[3] for p in both])
        out["auc_invalid_at_truth"] = {"n_invalid": int(inval.sum()), "n": len(both),
                                       "radius_1_over_r": auc(1.0 / rr, inval),
                                       "tolerance_1_over_t": auc(1.0 / tt, inval)}
    ext = [e for r in ok for e in r.get("extensions", []) if e.get("exact") is not None]
    d0_t1 = {r["key"]: tol_value(r["tol_order"]) == 1 for r in ok}
    agree = [(d0_t1[r["key"]] == e["t1"]) for r in ok for e in r.get("extensions", []) if e.get("decided", e["exact"])]
    out["extension_sensitivity"] = {"n_extension_checks": len(ext), "n_decided": len(agree),
                                    "share_same_depth1_verdict_as_d0": round(float(np.mean(agree)), 4) if agree else None,
                                    "share_extension_identical_to_d0": round(float(np.mean([e["same_as_d0"] for e in ext])), 4) if ext else None}
    out["by_arm_share_t1"] = {}
    for arm in sorted({r["arm"] for r in ok}):
        v = [tol_value(r["tol_order"]) for r in ok if r["arm"] == arm]
        v = [x for x in v if x is not None]
        out["by_arm_share_t1"][arm] = [sum(1 for x in v if x == 1), len(v)]
    sec = [r["tol_order"]["seconds"] for r in ok]
    out["seconds_order_median_max"] = [round(float(np.median(sec)), 3), round(float(max(sec)), 3)] if sec else None
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    rows = load(args.out)
    status = Counter(r.get("status") for r in rows)
    hard = [r["key"] for r in rows if r.get("hard_failures")]
    res = {"experiment": "E26", "n_rows": len(rows), "status": dict(status), "n_hard_failures": len(hard),
           "hard_failure_keys_first10": hard[:10],
           "primary": block(rows, "all released status-ok queries, all arms"),
           "secondary_e10_arms": block([r for r in rows if r["arm"] in E10_ARMS], "ten elicited arms and D_RAND")}
    p = res["primary"]
    tb, ci = p.get("tau_b"), p.get("tau_b_network_bootstrap95") or [None, None]
    if tb is None:
        verdict = "undefined (a measure is constant on the decided rows)"
        against = False
    else:
        against = tb >= DIES_TAU and ci[0] is not None and ci[0] > DIES_LOWER
        verdict = ("AGAINST the novelty sentence: tolerance and radius order queries alike"
                   if against else "novelty sentence survives: the tolerance does not reproduce the radius ordering")
    s1 = p["tolerance_order"]["share_t1"]
    res["verdict"] = {
        "registered_dies_if": f"tau_b >= {DIES_TAU} with bootstrap lower limit > {DIES_LOWER}",
        "against_method": against, "text": verdict,
        "prediction_share_t1_at_least_0.8": None if s1 is None else s1 >= 0.8,
        "prediction_tau_b_at_most_0.3": None if tb is None else tb <= 0.3,
    }
    (args.out / "results.json").write_text(json.dumps(res, indent=1))
    lines = [f"E26 rows {len(rows)} status {dict(status)} hard failures {len(hard)}",
             f"tolerance (order) distribution {p['tolerance_order']['distribution']}",
             f"share t=1 {s1} CP95 {p['tolerance_order']['share_t1_cp95']} (decided {p['tolerance_order']['n_decided']})",
             f"tau_b {tb} network bootstrap 95 {ci} on {p.get('n_in_tau')} queries",
             f"pairs {p.get('pairs')}", f"query level {p.get('query_level')}",
             f"AUC invalid at truth {p.get('auc_invalid_at_truth')}",
             f"extension sensitivity {p.get('extension_sensitivity')}",
             f"VERDICT {res['verdict']}"]
    (args.out / "analyse.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
