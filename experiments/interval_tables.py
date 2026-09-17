"""Tables 1 and 2, the Figure 1 candidates and the prediction scores for measurement M1.

Reads what ``experiments/interval_panel.py`` wrote under ``results/interval/`` and
applies the metric fixed in ``results/interval/PREREGISTRATION.md``: per cell,
each procedure's interval is scaled about its midpoint by the smallest factor that
reaches 95 % coverage of the true effect, and L95 is the mean scaled width over
the mean scaled width of the class-wide band on the same replicates, with a
network cluster bootstrap. Raw coverage and the raw width ratio sit beneath.

M1b adds the policy rows (the radius stop rule and the leave-one-out screens), the
starred comparison of STOP2 against SCREEN at two reversed claims, the depth-one
identity check, the refusal rates at correct knowledge and the ceiling of the starred
comparison, computed after every M1 cell with a separate bootstrap generator so no M1
number moves.

Writes ``TABLE1_L95.md``, ``table1.json``, ``TABLE2_REPORTS.md``, ``table2.json``,
``fig1_candidates.json`` and ``summary.json``.

    python experiments/interval_tables.py
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "interval"
DECLARED = ROOT / "results" / "frame" / "declared_rows.csv"

B_BOOT = 2000
BOOT_SEED = 20260912
TABLE1_PROCS = ("PLAIN", "SCREEN", "STOP1", "STOP2", "FAR1", "I1", "I2", "BLANKET")
POLICY_PROCS = ("SCREEN", "SCREEN2", "STOP1", "STOP2")
ALL_PROCS = ("PLAIN", "FAR1", "NEAR1", "I1", "I2", "BLANKET", "CH", *POLICY_PROCS)
M1B_SEED = 20260913
PANEL_B = ("ecoli70", "arth150", "magic-niab", "magic-irri")
CENSORED = 4
MIN_ROWS = 20


def clean(obj: object) -> object:
    """Replace non-finite floats by strings so every JSON output is standard JSON."""
    if isinstance(obj, float) and not math.isfinite(obj):
        return "nan" if math.isnan(obj) else ("inf" if obj > 0 else "-inf")
    if isinstance(obj, dict):
        return {k: clean(v) for k, v in obj.items()}
    if isinstance(obj, list | tuple):
        return [clean(v) for v in obj]
    return obj


def dump(path: Path, obj: object) -> None:
    """Write ``obj`` as indented, standard JSON."""
    path.write_text(json.dumps(clean(obj), indent=1))


def fnum(s: str | None) -> float:
    """A CSV cell as a float, NaN when empty."""
    return float(s) if s not in ("", None, "None") else float("nan")


# ----------------------------------------------------------------------------- metric


def lambdas(lo: np.ndarray, hi: np.ndarray, tau: np.ndarray, tol: np.ndarray) -> np.ndarray:
    """Per replicate, the scaling of the half-width about the midpoint needed to cover."""
    mid = (lo + hi) / 2
    half = (hi - lo) / 2
    gap = np.abs(tau - mid)
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = np.where(half > 0, gap / half, np.where(gap <= tol, 0.0, np.inf))
    return lam


def calibrated(lam: np.ndarray, half: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Calibrated mean width for each weight vector (rows of ``weights``)."""
    order = np.argsort(lam, kind="stable")
    lam_s = lam[order]
    w_s = weights[:, order]
    cum = np.cumsum(w_s, axis=1)
    tot = cum[:, -1]
    k = np.argmax(cum >= 0.95 * tot[:, None] - 1e-9, axis=1)
    lam_star = lam_s[k]
    mean_half = (weights * half[None, :]).sum(axis=1) / tot
    with np.errstate(invalid="ignore"):
        width = np.where(np.isinf(lam_star), np.inf, 2 * lam_star * mean_half)
    return width


def cell_metrics(
    reps: dict[str, np.ndarray],
    procs: tuple[str, ...],
    rng: np.random.Generator,
    finite: bool,
) -> dict:
    """L95 with bootstrap intervals, raw coverage and raw width ratio for one cell."""
    tau, tol, net = reps["tau"], reps["tol"], reps["net"]
    n = len(tau)
    out: dict = {"replicates": n, "rows": len(set(reps["row"])), "networks": len(set(net))}
    if n == 0:
        return out
    blo, bhi = reps["BLANKET_lo"], reps["BLANKET_hi"]
    b_width = float(np.mean(bhi - blo))
    nets, inv = np.unique(net, return_inverse=True)
    boots = rng.integers(0, len(nets), size=(B_BOOT, len(nets)))
    counts = np.stack([np.bincount(b, minlength=len(nets)) for b in boots])
    w_boot = counts[:, inv].astype(float)
    w_one = np.ones((1, n))
    lam_b = lambdas(blo, bhi, tau, tol)
    half_b = (bhi - blo) / 2
    cal_b = calibrated(lam_b, half_b, w_one)[0] if finite else None
    cal_b_boot = calibrated(lam_b, half_b, w_boot) if finite else None
    per: dict = {}
    for p in procs:
        lo, hi = reps[f"{p}_lo"], reps[f"{p}_hi"]
        ok = ~np.isnan(lo)
        if not ok.all():
            per[p] = {"available": int(ok.sum())}
            continue
        cover = (lo - tol <= tau) & (tau <= hi + tol)
        rec = {
            "coverage": float(cover.mean()),
            "width_ratio": float(np.mean(hi - lo) / b_width) if b_width > 0 else None,
        }
        if finite:
            lam = lambdas(lo, hi, tau, tol)
            half = (hi - lo) / 2
            cal = calibrated(lam, half, w_one)[0]
            boot = calibrated(lam, half, w_boot) / cal_b_boot
            rec["L95"] = float(cal / cal_b) if math.isfinite(cal) else None
            rec["L95_ci"] = [
                float(np.quantile(boot, 0.025, method="lower")),
                float(np.quantile(boot, 0.975, method="higher")),
            ]
            rec["L95_boot_infinite_share"] = float(np.mean(~np.isfinite(boot)))
            rec["_boot"] = boot
        per[p] = rec
    if finite and "_boot" in per.get("I1", {}) and "_boot" in per.get("PLAIN", {}):
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio = per["I1"]["_boot"] / per["PLAIN"]["_boot"]
        ratio = np.where(np.isnan(ratio), np.inf, ratio)
        out["ratio_I1_over_PLAIN_ci"] = [
            float(np.quantile(ratio, 0.025, method="lower")),
            float(np.quantile(ratio, 0.975, method="higher")),
        ]
    for rec in per.values():
        rec.pop("_boot", None)
    out["procedures"] = per
    return out


# ----------------------------------------------------------------------------- loading


def load(out: Path) -> tuple[list[dict], dict[tuple, dict], list[dict], list[dict], dict]:
    """Rows, rows keyed by (network, x, y, condition), replicates, profiles, panel summary."""
    rows = list(csv.DictReader((out / "rows.csv").open()))
    keyed = {
        (r["network"], r["x"], r["y"], r["condition"]): r for r in rows if r["source"] == "frame"
    }
    with gzip.open(out / "replicates.csv.gz", "rt") as fh:
        reps = list(csv.DictReader(fh))
    profiles = [json.loads(line) for line in (out / "profiles.jsonl").open()]
    panel = json.loads((out / "panel_summary.json").read_text())
    return rows, keyed, reps, profiles, panel


def select(
    reps: list[dict], keyed: dict[tuple, dict], cond: set[str], n: int, keep: object
) -> dict[str, np.ndarray]:
    """Replicates of the given conditions and ``n`` whose row passes ``keep``, as arrays."""
    chosen = []
    for r in reps:
        if r["condition"] not in cond or int(r["n"]) != n:
            continue
        row = keyed[(r["network"], r["x"], r["y"], r["condition"])]
        if keep(row):
            chosen.append((r, row))
    arr: dict[str, np.ndarray] = {
        "tau": np.array([fnum(r["tau"]) for r, _ in chosen]),
        "tol": np.array([fnum(row["tol"]) for _, row in chosen]),
        "net": np.array([r["network"] for r, _ in chosen]),
        "row": np.array([f"{r['network']}|{r['x']}|{r['y']}|{r['condition']}" for r, _ in chosen]),
        "dispute": np.array([fnum(row["dispute_share"]) for _, row in chosen]),
    }
    for p in ALL_PROCS:
        arr[f"{p}_lo"] = np.array([fnum(r[f"{p}_lo"]) for r, _ in chosen])
        arr[f"{p}_hi"] = np.array([fnum(r[f"{p}_hi"]) for r, _ in chosen])
    return arr


def subset(arr: dict[str, np.ndarray], mask: np.ndarray) -> dict[str, np.ndarray]:
    """The replicates selected by ``mask``."""
    return {k: v[mask] for k, v in arr.items()}


def dispute_summary(arr: dict[str, np.ndarray]) -> dict:
    """Median and mean direction-dispute share over the distinct rows of a cell."""
    seen: dict[str, float] = {}
    for key, d in zip(arr["row"], arr["dispute"], strict=True):
        seen[str(key)] = float(d)
    vals = np.array(list(seen.values()))
    if not len(vals):
        return {"rows": 0}
    return {"rows": len(vals), "median": float(np.median(vals)), "mean": float(np.mean(vals))}


# ----------------------------------------------------------------------------- rendering


def fmt_cell(rec: dict | None, finite: bool) -> str:
    """``L95 [lo, hi] / cov`` or ``width / cov`` for one procedure in one cell."""
    if rec is None or "coverage" not in rec:
        return "n/a"
    cov = f"{rec['coverage']:.3f}"
    if not finite:
        wr = rec.get("width_ratio")
        return f"{wr:.3f} / {cov}" if wr is not None else f"- / {cov}"
    if rec.get("L95") is None:
        return f"inf / {cov}"
    ci = rec.get("L95_ci")
    ci_s = f" [{ci[0]:.2f}, {ci[1]:.2f}]".replace("inf]", "inf)") if ci else ""
    return f"{rec['L95']:.3f}{ci_s} / {cov}"


def table_md(
    title: str, cols: list[str], cells: dict[str, dict], procs: tuple[str, ...], finite: bool
) -> list[str]:
    """One markdown table: procedures as rows, cells as columns."""
    lines = [
        f"**{title}**",
        "",
        "| procedure | " + " | ".join(cols) + " |",
        "|---|" + "---|" * len(cols),
    ]
    for p in procs:
        vals = [fmt_cell(cells[c].get("procedures", {}).get(p), finite) for c in cols]
        lines.append(f"| {p} | " + " | ".join(vals) + " |")
    lines.append(
        "| informative rows (replicates) | "
        + " | ".join(f"{cells[c].get('rows', 0)} ({cells[c].get('replicates', 0)})" for c in cols)
        + " |"
    )
    lines.append(
        "| networks | " + " | ".join(str(cells[c].get("networks", 0)) for c in cols) + " |"
    )
    lines.append(
        "| direction-dispute share, median | "
        + " | ".join(
            (
                f"{cells[c]['dispute']['median']:.3f}"
                if cells[c].get("dispute", {}).get("rows")
                else "-"
            )
            for c in cols
        )
        + " |"
    )
    lines.append("")
    return lines


# ----------------------------------------------------------------------------- M1b


def boot_indices(nets: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Network-cluster bootstrap weights, shape ``(B_BOOT, len(nets))``."""
    uniq, inv = np.unique(nets, return_inverse=True)
    draws = rng.integers(0, len(uniq), size=(B_BOOT, len(uniq)))
    counts = np.stack([np.bincount(d, minlength=len(uniq)) for d in draws])
    return counts[:, inv].astype(float)


def ci(values: np.ndarray) -> list[float]:
    """Percentile 95 % interval as order statistics, so infinities are kept."""
    return [
        float(np.quantile(values, 0.025, method="lower")),
        float(np.quantile(values, 0.975, method="higher")),
    ]


def paired(
    arr: dict[str, np.ndarray], a: str, b: str, rng: np.random.Generator, finite: bool
) -> dict:
    """Policy ``a`` minus policy ``b`` on the same replicates, with one bootstrap for both."""
    tau, tol, nets = arr["tau"], arr["tol"], arr["net"]
    out: dict = {"replicates": len(tau), "rows": len(set(arr["row"])), "networks": len(set(nets))}
    if not len(tau):
        return out
    cov = {p: (arr[f"{p}_lo"] - tol <= tau) & (tau <= arr[f"{p}_hi"] + tol) for p in (a, b)}
    w = boot_indices(nets, rng)
    d = cov[a].astype(float) - cov[b].astype(float)
    out["coverage"] = {a: float(cov[a].mean()), b: float(cov[b].mean())}
    out["coverage_diff"] = float(d.mean())
    out["coverage_diff_ci"] = ci((w * d).sum(axis=1) / w.sum(axis=1))
    per_net = {}
    for net in sorted(set(nets)):
        m = nets == net
        per_net[str(net)] = {"replicates": int(m.sum()), "diff": float(d[m].mean())}
    out["per_network"] = per_net
    if finite:
        one = np.ones((1, len(tau)))
        cal = {}
        for p in (a, b, "BLANKET"):
            lo, hi = arr[f"{p}_lo"], arr[f"{p}_hi"]
            lam, half = lambdas(lo, hi, tau, tol), (hi - lo) / 2
            cal[p] = (calibrated(lam, half, one)[0], calibrated(lam, half, w))
        l95 = {p: cal[p][0] / cal["BLANKET"][0] for p in (a, b)}
        boot = {p: cal[p][1] / cal["BLANKET"][1] for p in (a, b)}
        with np.errstate(invalid="ignore"):
            diff_b = np.where(np.isnan(boot[a] - boot[b]), 0.0, boot[a] - boot[b])
        out["L95"] = {p: float(v) for p, v in l95.items()}
        out["L95_diff"] = float(l95[a] - l95[b])
        out["L95_diff_ci"] = ci(diff_b)
    return out


def m1b_section(
    reps: list[dict], keyed: dict[tuple, dict], rows: list[dict]
) -> tuple[dict, list[str], dict]:
    """The starred comparison, the identity check, the refusal rates and the ceiling."""
    rng = np.random.default_rng(M1B_SEED)
    panels = ("A", "B")
    frame = [r for r in rows if r["source"] == "frame" and r.get("stop1_certified", "") != ""]

    def committed(r: dict) -> bool:
        return r.get("plain_kind") == "committed"

    def rval_is(r: dict, v: int) -> bool:
        return r.get("r_val", "") not in ("", None) and int(r["r_val"]) == v

    out: dict = {
        "starred": {},
        "starred_depth_matched": {},
        "decomposition": {},
        "identity": {},
        "refusal": {},
        "ceiling": {},
    }
    md = [
        "## M1b: policy rows, the starred comparison and the price of the guarantee",
        "",
        "SCREEN: leave one claim out, re-close, widen to I1 if the committed result changes. "
        "STOP1, STOP2: report PLAIN when the committed set survives every retraction of at "
        "most r claims, otherwise the band I_r. SCREEN2: the screen over one or two claims, "
        "widening to I2. Policy comparisons are on informative rows with a committed set; "
        "the all-row numbers use the literal reading (no set: the stop rule widens, the "
        "screen compares verdicts). Network cluster bootstrap, 2000 resamples, seed "
        f"{M1B_SEED}.",
        "",
    ]

    # ---- the starred comparison at b = 2
    md += [
        "**Starred comparison: STOP2 minus SCREEN at b = 2**",
        "",
        "| panel | rows | n | coverage STOP2 / SCREEN | difference [95 %] | L95 STOP2 / "
        "SCREEN | difference [95 %] | networks |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for pop_name, keep_pop in (("committed", committed), ("all", lambda r: True)):
        for pan in panels:
            for n in (0, 1000, 20000):
                arr = select(
                    reps,
                    keyed,
                    {"CTRL_b2"},
                    n,
                    lambda r, pan=pan, kp=keep_pop: r["panel"] == pan and kp(r),
                )
                res = paired(arr, "STOP2", "SCREEN", rng, n > 0)
                out["starred"][f"{pop_name}|{pan}|{n}"] = res
                dm = paired(arr, "STOP2", "SCREEN2", rng, n > 0)
                out["starred_depth_matched"][f"{pop_name}|{pan}|{n}"] = dm
                if not res.get("replicates"):
                    continue
                cov = res["coverage"]
                l95 = (
                    (
                        f"{res['L95']['STOP2']:.3f} / {res['L95']['SCREEN']:.3f} | "
                        f"{res['L95_diff']:+.3f} [{res['L95_diff_ci'][0]:+.3f}, "
                        f"{res['L95_diff_ci'][1]:+.3f}]"
                    )
                    if n
                    else "- | -"
                )
                md.append(
                    f"| {pan} | {pop_name} | {'population' if n == 0 else n} | "
                    f"{cov['STOP2']:.3f} / {cov['SCREEN']:.3f} | {res['coverage_diff']:+.3f} "
                    f"[{res['coverage_diff_ci'][0]:+.3f}, {res['coverage_diff_ci'][1]:+.3f}] | "
                    f"{l95} | {res['networks']} |"
                )
                if pop_name == "committed" and n:
                    strata = {}
                    for label, fn in (
                        ("interaction", lambda r: r.get("interaction") == "1"),
                        ("r_val = 1", lambda r: rval_is(r, 1)),
                        ("rest", lambda r: r.get("interaction") != "1" and not rval_is(r, 1)),
                    ):
                        mask = np.array(
                            [fn(keyed[(*k.split("|")[:3], "CTRL_b2")]) for k in arr["row"]],
                            dtype=bool,
                        )
                        sub = subset(arr, mask)
                        if not len(sub["tau"]):
                            strata[label] = {"rows": 0}
                            continue
                        ca = (sub["STOP2_lo"] - sub["tol"] <= sub["tau"]) & (
                            sub["tau"] <= sub["STOP2_hi"] + sub["tol"]
                        )
                        cb = (sub["SCREEN_lo"] - sub["tol"] <= sub["tau"]) & (
                            sub["tau"] <= sub["SCREEN_hi"] + sub["tol"]
                        )
                        diff = float(ca.mean() - cb.mean())
                        share = len(sub["tau"]) / len(arr["tau"])
                        strata[label] = {
                            "rows": len(set(sub["row"])),
                            "replicate_share": share,
                            "coverage_diff": diff,
                            "contribution": share * diff,
                        }
                    out["decomposition"][f"{pan}|{n}"] = strata
    md += [
        "",
        "Depth-matched: STOP2 minus SCREEN2 (both widen to I2), informative rows with "
        "a committed set:",
        "",
    ]
    for pan in panels:
        for n in (1000, 20000):
            dm = out["starred_depth_matched"][f"committed|{pan}|{n}"]
            if dm.get("replicates"):
                md.append(
                    f"- Panel {pan}, n = {n}: coverage difference "
                    f"{dm['coverage_diff']:+.3f} [{dm['coverage_diff_ci'][0]:+.3f}, "
                    f"{dm['coverage_diff_ci'][1]:+.3f}], L95 difference "
                    f"{dm['L95_diff']:+.3f} [{dm['L95_diff_ci'][0]:+.3f}, "
                    f"{dm['L95_diff_ci'][1]:+.3f}]"
                )
    md += [
        "",
        "Where the STOP2 minus SCREEN coverage gap comes from (committed rows, b = 2): "
        "stratum, rows, share of replicates, coverage difference inside it, contribution to "
        "the pooled difference.",
        "",
    ]
    for key, strata in out["decomposition"].items():
        parts = [
            f"{lab} {v['rows']} rows, {v.get('replicate_share', 0):.2f}, "
            f"{v.get('coverage_diff', 0):+.3f}, {v.get('contribution', 0):+.3f}"
            for lab, v in strata.items()
        ]
        md.append(f"- Panel {key.replace('|', ', n = ')}: " + "; ".join(parts))

    # ---- identity check, refusal and ceiling, per condition
    def decisions(r: dict) -> dict:
        return {
            "STOP1": r["stop1_certified"] == "0",
            "STOP2": r["stop2_certified"] == "0",
            "SCREEN": r["screen_triggered"] == "1",
            "SCREEN2": r["screen2_triggered"] == "1",
        }

    def reports_differ(r: dict) -> bool:
        stop2 = "PLAIN" if r["stop2_certified"] == "1" else "I2"
        screen = "I1" if r["screen_triggered"] == "1" else "PLAIN"
        if stop2 == screen:
            return False
        return not (stop2 == "I2" and screen == "I1" and r["sets_b1_equal_b2"] == "1")

    for cond in ("CTRL_b0", "CTRL_b1", "CTRL_b2", "D_LLM", "D_LLM_72B"):
        for pop_name, keep_pop in (("committed", committed), ("all", lambda r: True)):
            for pan in (*panels, "pooled"):
                sel = [
                    r
                    for r in frame
                    if r["condition"] == cond
                    and r.get("status") == "ok"
                    and r["informative"] == "1"
                    and keep_pop(r)
                    and (pan == "pooled" or r["panel"] == pan)
                ]
                key = f"{cond}|{pop_name}|{pan}"
                classes: collections.Counter = collections.Counter()
                for r in sel:
                    dec = decisions(r)
                    if dec["STOP1"] == dec["SCREEN"]:
                        continue
                    if not committed(r):
                        classes["no committed set"] += 1
                    elif dec["SCREEN"] and r["screen_change_z_valid"] == "1":
                        classes["set changes while Z stays valid"] += 1
                    elif not dec["SCREEN"]:
                        classes["screen silent while Z breaks"] += 1
                    else:
                        classes["other"] += 1
                out["identity"][key] = {
                    "rows": len(sel),
                    "differ": sum(classes.values()),
                    "classes": dict(classes),
                }
                if not sel:
                    continue
                nets = np.array([r["network"] for r in sel])
                w = boot_indices(nets, rng)
                ref = {}
                for p in POLICY_PROCS:
                    v = np.array([decisions(r)[p] for r in sel], dtype=float)
                    ref[p] = {"rate": float(v.mean()), "ci": ci((w * v).sum(1) / w.sum(1))}
                out["refusal"][key] = {"rows": len(sel), "networks": len(set(nets)), **ref}
                inter = np.array([r.get("interaction") == "1" for r in sel], dtype=float)
                lit = np.array([reports_differ(r) for r in sel], dtype=float)
                out["ceiling"][key] = {
                    "rows": len(sel),
                    "interaction": float(inter.mean()),
                    "interaction_ci": ci((w * inter).sum(1) / w.sum(1)),
                    "reports_differ": float(lit.mean()),
                    "r_val_1": float(np.mean([rval_is(r, 1) for r in sel])),
                    "r_val_2": float(np.mean([rval_is(r, 2) for r in sel])),
                }

    md += [
        "",
        "**Price of the guarantee: refusal rate at b = 0 (knowledge correct)**, share of "
        "informative rows on which the policy widens [network bootstrap 95 %]",
        "",
        "| panel | rows | n rows (networks) | STOP1 | STOP2 | SCREEN | SCREEN2 |",
        "|---|---|---|---|---|---|---|",
    ]
    for pop_name in ("committed", "all"):
        for pan in panels:
            rr = out["refusal"].get(f"CTRL_b0|{pop_name}|{pan}")
            if not rr:
                continue
            cells = [
                f"{rr[p]['rate']:.3f} [{rr[p]['ci'][0]:.2f}, {rr[p]['ci'][1]:.2f}]"
                for p in ("STOP1", "STOP2", "SCREEN", "SCREEN2")
            ]
            md.append(
                f"| {pan} | {pop_name} | {rr['rows']} ({rr['networks']}) | "
                + " | ".join(cells)
                + " |"
            )
    md += [
        "",
        "**Identity check and ceiling** (informative rows; `differ` = STOP1 and SCREEN "
        "take different decisions; interaction = r_val = 2, or a pair retraction changes the "
        "committed result while no single one does; reports differ = STOP2 and SCREEN "
        "report different intervals)",
        "",
        "| condition | rows | panel | n rows | STOP1 vs SCREEN differ (classes) | interaction "
        "[95 %] | reports differ | r_val = 1 | r_val = 2 |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for cond in ("CTRL_b0", "CTRL_b1", "CTRL_b2", "D_LLM", "D_LLM_72B"):
        for pop_name in ("committed", "all"):
            for pan in (*panels, "pooled"):
                key = f"{cond}|{pop_name}|{pan}"
                idn, cl = out["identity"].get(key), out["ceiling"].get(key)
                if not idn or not cl:
                    continue
                cls = ", ".join(f"{k}: {v}" for k, v in idn["classes"].items()) or "-"
                md.append(
                    f"| {cond} | {pop_name} | {pan} | {idn['rows']} | {idn['differ']} "
                    f"({cls}) | {cl['interaction']:.3f} [{cl['interaction_ci'][0]:.2f}, "
                    f"{cl['interaction_ci'][1]:.2f}] | {cl['reports_differ']:.3f} | "
                    f"{cl['r_val_1']:.3f} | {cl['r_val_2']:.3f} |"
                )
    md.append("")

    # ---- SCREEN2 as a row, with the policy rows, committed population
    policy_cells: dict = {}
    md += [
        "**Policy rows with SCREEN2, informative rows with a committed set** "
        "(L95 [95 %] / raw coverage)",
        "",
    ]
    for pan in panels:
        for n in (1000, 20000):
            cells = {}
            for b in (0, 1, 2):
                arr = select(
                    reps,
                    keyed,
                    {f"CTRL_b{b}"},
                    n,
                    lambda r, pan=pan: r["panel"] == pan and committed(r),
                )
                cells[f"b={b}"] = cell_metrics(
                    arr,
                    ("PLAIN", "SCREEN", "SCREEN2", "STOP1", "STOP2", "I1", "I2", "BLANKET"),
                    rng,
                    True,
                )
                cells[f"b={b}"]["dispute"] = dispute_summary(arr)
            policy_cells[f"{pan}|{n}"] = cells
            md += table_md(
                f"Panel {pan}, n = {n}",
                list(cells),
                cells,
                ("PLAIN", "SCREEN", "SCREEN2", "STOP1", "STOP2", "I1", "I2", "BLANKET"),
                True,
            )
    out["policy_cells"] = policy_cells

    # ---- the predictions
    t1 = {}
    for pan in panels:
        for n in (1000, 20000):
            for b in (0, 1):
                arr = select(reps, keyed, {f"CTRL_b{b}"}, n, lambda r, pan=pan: r["panel"] == pan)
                t1[f"{pan}|{n}|{b}"] = cell_metrics(arr, ("STOP1", "I1"), rng, True)
    preds: dict = {}
    idn = out["identity"]["CTRL_b1|committed|pooled"]
    share = idn["differ"] / idn["rows"] if idn["rows"] else float("nan")
    preds["Q1"] = {
        "differ": idn["differ"],
        "rows": idn["rows"],
        "share": share,
        "classes": idn["classes"],
        "all_rows": out["identity"]["CTRL_b1|all|pooled"],
        "verdict": "supported" if share <= 0.01 else "falsified",
    }
    parts = {}
    for n in (1000, 20000):
        res = out["starred"][f"committed|A|{n}"]
        parts[f"coverage_A_{n}"] = bool(res["coverage_diff_ci"][0] > 0)
    for pan in panels:
        parts[f"ceiling_{pan}"] = bool(
            out["ceiling"][f"CTRL_b2|committed|{pan}"]["interaction"] <= 0.31
        )
    preds["Q2"] = {
        "parts": parts,
        "verdict": "supported"
        if all(parts.values())
        else ("partial" if any(parts.values()) else "falsified"),
    }
    q3 = {}
    q5 = {}
    for pan in panels:
        for n in (1000, 20000):
            c0 = t1[f"{pan}|{n}|0"]["procedures"]
            q3[f"{pan}|{n}"] = bool(
                (c0["STOP1"]["L95"] or math.inf) < (c0["I1"]["L95"] or math.inf)
            )
            c1 = t1[f"{pan}|{n}|1"]["procedures"]
            q5[f"{pan}|{n}"] = bool(c1["STOP1"]["coverage"] >= 0.95)

    def verdict(cells: dict) -> str:
        v = list(cells.values())
        return "supported" if all(v) else ("partial" if any(v) else "falsified")

    preds["Q3"] = {"cells": q3, "verdict": verdict(q3)}
    q4 = {
        pan: bool(out["refusal"][f"CTRL_b0|committed|{pan}"]["STOP1"]["rate"] < 0.5)
        for pan in panels
    }
    preds["Q4"] = {
        "cells": q4,
        "rates": {pan: out["refusal"][f"CTRL_b0|committed|{pan}"] for pan in panels},
        "verdict": verdict(q4),
    }
    preds["Q5"] = {"cells": q5, "verdict": verdict(q5)}
    out["q_cells"] = t1
    return out, md, preds


# ----------------------------------------------------------------------------- main


def main(argv: list[str] | None = None) -> None:
    """Build every table and score the predictions."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="")
    args = ap.parse_args(argv)
    out = Path(args.out) if args.out else OUT
    rows, keyed, reps, profiles, panel = load(out)
    rng = np.random.default_rng(BOOT_SEED)
    panels = ("A", "B")
    ns = (0, 1000, 20000)

    def is_panel(p: str) -> object:
        return lambda row: row["panel"] == p

    reach2 = {
        (r["network"], r["x"], r["y"])
        for r in rows
        if r["condition"] == "CTRL_b2" and r.get("status") == "ok"
    }

    table1: dict = {
        "controlled": {},
        "controlled_matched": {},
        "controlled_committed": {},
        "elicited": {},
        "elicited_by_arm": {},
        "ch_defined": {},
        "ch_descendant": {},
    }
    md: list[str] = [
        "# Table 1 — L95 and raw coverage (measurement M1)",
        "",
        "Each cell: L95 [network-bootstrap 95 % interval] / raw coverage. L95 is "
        "the mean width after scaling each interval about its midpoint by the "
        "smallest factor giving 95 % coverage in the cell, over the same for "
        "BLANKET. Population cells print raw width ratio / raw coverage. "
        "Informative rows only. Panel A: semi-synthetic parameters on the swept "
        "structures. Panel B: fitted coefficients, four networks, so its bootstrap "
        "has four clusters.",
        "",
    ]
    for pan in panels:
        for n in ns:
            finite = n > 0
            cells, cells_m, cells_c, cells_ch = {}, {}, {}, {}
            for b in (0, 1, 2):
                cond = f"CTRL_b{b}"
                arr = select(reps, keyed, {cond}, n, is_panel(pan))
                cell = cell_metrics(arr, (*TABLE1_PROCS, "NEAR1"), rng, finite)
                cell["dispute"] = dispute_summary(arr)
                cells[f"b={b}"] = cell
                mask = np.array([tuple(k.split("|")[:3]) in reach2 for k in arr["row"]], dtype=bool)
                cm = cell_metrics(subset(arr, mask), TABLE1_PROCS, rng, finite)
                cm["dispute"] = dispute_summary(subset(arr, mask))
                cells_m[f"b={b}"] = cm
                committed = np.array(
                    [
                        keyed[(*k.split("|")[:3], cond)]["plain_kind"] == "committed"
                        for k in arr["row"]
                    ],
                    dtype=bool,
                )
                cc = cell_metrics(subset(arr, committed), TABLE1_PROCS, rng, finite)
                cc["dispute"] = dispute_summary(subset(arr, committed))
                cells_c[f"b={b}"] = cc
                ch_def = np.array(
                    [keyed[(*k.split("|")[:3], cond)].get("ch_defined") == "1" for k in arr["row"]],
                    dtype=bool,
                )
                chd = cell_metrics(
                    subset(arr, ch_def), ("CH", "PLAIN", "I1", "BLANKET"), rng, finite
                )
                chd["dispute"] = dispute_summary(subset(arr, ch_def))
                cells_ch[f"b={b}"] = chd
            label = "population" if n == 0 else f"n = {n}"
            table1["controlled"][f"{pan}|{n}"] = cells
            table1["controlled_matched"][f"{pan}|{n}"] = cells_m
            table1["controlled_committed"][f"{pan}|{n}"] = cells_c
            table1["ch_defined"][f"{pan}|{n}"] = cells_ch
            md += table_md(
                f"Panel {pan}, controlled reversals, {label}",
                list(cells),
                cells,
                (*TABLE1_PROCS, "NEAR1"),
                finite,
            )
        # elicited, stratified by realised false claims in the treatment's component
        for n in ns:
            finite = n > 0
            cells = {}
            for stratum in ("0", "1", "2+"):

                def keep(row: dict, stratum: str = stratum, pan: str = pan) -> bool:
                    if row["panel"] != pan or row.get("n_false_local", "") == "":
                        return False
                    k = int(row["n_false_local"])
                    return str(k) == stratum if stratum != "2+" else k >= 2

                arr = select(reps, keyed, {"D_LLM", "D_LLM_72B"}, n, keep)
                cell = cell_metrics(arr, TABLE1_PROCS, rng, finite)
                cell["dispute"] = dispute_summary(arr)
                cells[f"false={stratum}"] = cell
                for arm in ("D_LLM", "D_LLM_72B"):
                    a_arr = select(reps, keyed, {arm}, n, keep)
                    table1["elicited_by_arm"][f"{pan}|{n}|{arm}|{stratum}"] = cell_metrics(
                        a_arr, TABLE1_PROCS, rng, finite
                    )
            table1["elicited"][f"{pan}|{n}"] = cells
            label = "population" if n == 0 else f"n = {n}"
            md += table_md(
                f"Panel {pan}, elicited knowledge (both suppliers pooled) by false "
                f"claims in the treatment's component, {label}",
                list(cells),
                cells,
                TABLE1_PROCS,
                finite,
            )
    md += ["## Secondary: one population for all three columns (rows where b = 2 is reachable)", ""]
    for pan in panels:
        for n in (1000, 20000):
            c = table1["controlled_matched"][f"{pan}|{n}"]
            md += table_md(f"Panel {pan}, n = {n}", list(c), c, TABLE1_PROCS, True)
    md += ["## Secondary: rows where the analyst graph commits to a set", ""]
    for pan in panels:
        for n in (1000, 20000):
            c = table1["controlled_committed"][f"{pan}|{n}"]
            md += table_md(f"Panel {pan}, n = {n}", list(c), c, TABLE1_PROCS, True)
    md += [
        "## Feasibility M2: the oracle confounding bound on rows where its parameter is defined",
        "",
    ]
    for pan in panels:
        for n in (1000, 20000):
            c = table1["ch_defined"][f"{pan}|{n}"]
            md += table_md(
                f"Panel {pan}, n = {n}", list(c), c, ("CH", "PLAIN", "I1", "BLANKET"), True
            )

    # ---- P9: the CH formula on descendant rows
    def desc_keep(row: dict) -> bool:
        return (
            row.get("z_has_true_descendant") == "1"
            and row.get("z_valid_at_truth") == "0"
            and row.get("pop_CH_lo", "") != ""
        )

    p9 = {}
    for n in (0, 20000):
        arr = select(reps, keyed, {"CTRL_b1", "CTRL_b2", "D_LLM", "D_LLM_72B"}, n, desc_keep)
        p9[n] = cell_metrics(arr, ("CH", "PLAIN", "I1", "BLANKET"), rng, n > 0)
    table1["ch_descendant"] = p9
    ch_def_pop = select(
        reps,
        keyed,
        set(("CTRL_b0", "CTRL_b1", "CTRL_b2", "D_LLM", "D_LLM_72B")),
        0,
        lambda row: row.get("ch_defined") == "1",
    )
    ch_pop = cell_metrics(ch_def_pop, ("CH",), rng, False)
    m1b, m1b_md, m1b_preds = m1b_section(reps, keyed, rows)
    table1["m1b"] = m1b
    md += m1b_md
    dump(out / "table1.json", table1)
    (out / "TABLE1_L95.md").write_text("\n".join(md) + "\n")

    # ----------------------------------------------------------------- stratum, r0, P7
    frame_rows = {
        (r["network"], r["x"], r["y"]): r
        for r in rows
        if r["source"] == "frame" and r.get("class_width", "") != ""
    }
    inf_all = [r for r in frame_rows.values() if r["informative"] == "1"]
    strat = {
        "frame_rows": len(frame_rows),
        "informative": len(inf_all),
        "share": len(inf_all) / len(frame_rows),
        "by_panel": {
            p: {
                "frame_rows": sum(1 for r in frame_rows.values() if r["panel"] == p),
                "informative": sum(1 for r in inf_all if r["panel"] == p),
            }
            for p in panels
        },
        "class_censored_rows": sum(
            1 for r in rows if r.get("status") == "class_parent_sets_censored"
        )
        // 6,
    }
    disp = np.array([fnum(r["dispute_share"]) for r in inf_all])
    rshare = np.array([fnum(r["restricted_share"]) for r in inf_all])
    dispute = {
        "rows": len(disp),
        "dispute_share_median": float(np.median(disp)),
        "dispute_share_mean": float(np.mean(disp)),
        "restricted_share_mean": float(np.mean(rshare)),
        "restricted_share_exactly_zero": float(np.mean(rshare <= 1e-9)),
        "by_panel": {
            p: {
                "median": float(
                    np.median([fnum(r["dispute_share"]) for r in inf_all if r["panel"] == p])
                ),
                "mean": float(
                    np.mean([fnum(r["dispute_share"]) for r in inf_all if r["panel"] == p])
                ),
            }
            for p in panels
            if any(r["panel"] == p for r in inf_all)
        },
    }

    def r0_dist(sel: list[dict]) -> dict:
        c = collections.Counter(
            r["r0_pop_status"] if r["r0_pop_status"] != "exact" else r["r0_pop"] for r in sel
        )
        return {str(k): v for k, v in sorted(c.items(), key=lambda kv: str(kv[0]))}

    ok_rows = [r for r in rows if r.get("status") == "ok" and r["source"] == "frame"]
    r0 = {"population_all_rows": {}, "population_informative": {}, "band_vs_plugin": {}}
    for cond in ("CTRL_b0", "CTRL_b1", "CTRL_b2", "D_LLM", "D_LLM_72B", "A_ORACLE"):
        for pan in panels:
            sel = [r for r in ok_rows if r["condition"] == cond and r["panel"] == pan]
            r0["population_all_rows"][f"{cond}|{pan}"] = r0_dist(sel)
            r0["population_informative"][f"{cond}|{pan}"] = r0_dist(
                [r for r in sel if r["informative"] == "1"]
            )

    def above(sel: list[dict], key: str) -> dict:
        tot = up = down = eq = 0
        for r in sel:
            pop = int(r["r0_pop"])
            for v in r[key].split(";"):
                v = int(v)
                tot += 1
                up += v > pop
                down += v < pop
                eq += v == pop
        return {
            "replicates": tot,
            "above_population": up / tot if tot else None,
            "below_population": down / tot if tot else None,
            "equal": eq / tot if tot else None,
        }

    ctrl = [r for r in ok_rows if r["condition"] in ("CTRL_b0", "CTRL_b1", "CTRL_b2")]
    for scope, sel in (
        ("controlled_all_rows", ctrl),
        ("controlled_informative", [r for r in ctrl if r["informative"] == "1"]),
        ("elicited_all_rows", [r for r in ok_rows if r["condition"] in ("D_LLM", "D_LLM_72B")]),
    ):
        for n in (1000, 20000):
            r0["band_vs_plugin"][f"{scope}|n{n}"] = {
                "plugin": above(sel, f"r0_plug_n{n}"),
                "band": above(sel, f"r0_band_n{n}"),
                "by_panel": {
                    p: {
                        "plugin": above([r for r in sel if r["panel"] == p], f"r0_plug_n{n}"),
                        "band": above([r for r in sel if r["panel"] == p], f"r0_band_n{n}"),
                    }
                    for p in panels
                },
            }

    # ----------------------------------------------------------------- censoring
    cens: dict = collections.Counter()
    for r in rows:
        st = r.get("status", "")
        if st != "ok":
            cens[f"{r['condition']}|{r['panel']}|{st}"] += 1
        if r.get("ball_status") not in ("", None, "exhaustive", "no_component"):
            cens[f"{r['condition']}|{r['panel']}|ball_{r['ball_status']}"] += 1
    rval_src = collections.Counter(f"{r['condition']}|{r.get('r_val_source', '')}" for r in ok_rows)
    mismatch = sum(1 for r in ok_rows if r.get("ledger_verdict_match") == "0")
    ch_counts = collections.Counter(
        f"{r['condition']}|{r['panel']}|{r.get('ch_reason', '')}"
        for r in ok_rows
        if r.get("plain_kind") == "committed"
    )
    plain_kinds = collections.Counter(
        f"{r['condition']}|{r['panel']}|{r.get('plain_kind')}"
        for r in ok_rows
        if r["informative"] == "1"
    )

    # ----------------------------------------------------------------- predictions
    def cellp(kind: str, pan: str, n: int, b: str) -> dict:
        return table1[kind][f"{pan}|{n}"][b].get("procedures", {})

    def l95(c: dict, p: str) -> float:
        v = c.get(p, {}).get("L95")
        return math.inf if v is None else v

    preds: dict = {}
    cells4 = [(pan, n) for pan in panels for n in (1000, 20000)]

    def verdict(results: dict) -> str:
        scored = [v for v in results.values() if v is not None]
        if not scored:
            return "not scored"
        if all(scored):
            return "supported"
        return "partial" if any(scored) else "falsified"

    res = {}
    for pan, n in cells4:
        c = cellp("controlled", pan, n, "b=0")
        others = [l95(c, p) for p in ("FAR1", "I1", "I2", "BLANKET")]
        res[f"{pan}|{n}"] = bool(l95(c, "PLAIN") < min(others) and l95(c, "I1") > l95(c, "PLAIN"))
    preds["P1"] = {"cells": res, "verdict": verdict(res)}
    res = {}
    for pan, n in cells4:
        c = cellp("controlled", pan, n, "b=1")
        res[f"{pan}|{n}"] = bool(
            c["PLAIN"]["coverage"] < 0.95 and l95(c, "PLAIN") > l95(c, "I1") and l95(c, "I1") < 1.0
        )
    preds["P2"] = {"cells": res, "verdict": verdict(res)}
    res = {}
    for pan, n in cells4:
        ok = True
        for b in ("b=1", "b=2"):
            c = cellp("controlled", pan, n, b)
            ok &= abs(c["FAR1"]["coverage"] - c["PLAIN"]["coverage"]) <= 0.02
        res[f"{pan}|{n}"] = bool(ok)
    preds["P3"] = {"cells": res, "verdict": verdict(res)}
    res = {}
    for pan, n in cells4:
        c = cellp("controlled", pan, n, "b=2")
        res[f"{pan}|{n}"] = bool(c["I2"]["coverage"] >= 0.95 and c["I1"]["coverage"] < 0.95)
    preds["P4"] = {"cells": res, "verdict": verdict(res)}
    preds["P5"] = {
        "median_dispute_share": dispute["dispute_share_median"],
        "verdict": "supported" if dispute["dispute_share_median"] >= 0.5 else "falsified",
    }
    preds["P6"] = {
        "share": strat["share"],
        "verdict": "supported" if 0.15 <= strat["share"] <= 0.40 else "falsified",
    }
    bvp = r0["band_vs_plugin"]["controlled_all_rows|n1000"]
    preds["P7"] = {
        "plugin_above": bvp["plugin"]["above_population"],
        "band_above": bvp["band"]["above_population"],
        "verdict": "supported"
        if bvp["plugin"]["above_population"] >= 0.15 and bvp["band"]["above_population"] <= 0.05
        else "falsified",
    }
    res = {}
    for pan, n in cells4:
        e = table1["elicited"][f"{pan}|{n}"]
        c0, c1 = e["false=0"], e["false=1"]
        if c0.get("rows", 0) < MIN_ROWS or c1.get("rows", 0) < MIN_ROWS:
            res[f"{pan}|{n}"] = None
            continue
        p0, p1 = c0["procedures"], c1["procedures"]
        others = [l95(p0, p) for p in ("FAR1", "I1", "I2", "BLANKET")]
        first = l95(p0, "PLAIN") < min(others) and l95(p0, "I1") > l95(p0, "PLAIN")
        second = (
            p1["PLAIN"]["coverage"] < 0.95
            and l95(p1, "PLAIN") > l95(p1, "I1")
            and l95(p1, "I1") < 1.0
        )
        res[f"{pan}|{n}"] = bool(first and second)
    preds["P8"] = {
        "cells": res,
        "verdict": verdict(res),
        "rows": {
            f"{pan}|{n}": {
                s: table1["elicited"][f"{pan}|{n}"][s].get("rows", 0)
                for s in ("false=0", "false=1", "false=2+")
            }
            for pan, n in cells4
        },
    }
    p9c = p9[20000].get("procedures", {}).get("CH", {})
    p9_rows = p9[20000].get("rows", 0)
    preds["P9"] = {
        "rows": p9_rows,
        "coverage_n20000": p9c.get("coverage"),
        "population_coverage_on_descendant_rows": p9[0]
        .get("procedures", {})
        .get("CH", {})
        .get("coverage"),
        "check_population_coverage_on_defined_rows": ch_pop.get("procedures", {})
        .get("CH", {})
        .get("coverage"),
        "verdict": (
            "underpowered"
            if p9_rows < MIN_ROWS
            else "supported"
            if p9c.get("coverage", 1.0) < 0.95
            else "falsified"
        ),
    }

    # ----------------------------------------------------------------- Table 2
    declared_all = list(csv.DictReader(DECLARED.open()))
    t2: list[dict] = []

    def r0_seed_range(p: dict) -> str:
        vals = p["r0_band_seeds"]
        lab = [">3" if v == CENSORED else str(v) for v in (min(vals), max(vals))]
        return lab[0] if lab[0] == lab[1] else f"{lab[0]}..{lab[1]}"

    def r0_label(v: int, status: str = "exact") -> str:
        return str(v) if v != CENSORED else ("never" if status == "never" else ">3")

    def rval_label(v: object, status: str) -> str:
        if v not in ("", None):
            return str(v)
        if status in ("unreached", "no_retractable_edges"):
            return "never"
        if status.startswith("gt_"):
            return ">" + status.split("_")[1]
        if status.startswith("censored_at_depth_"):
            return ">=" + status.rsplit("_", 1)[1]
        return status

    def t2_row(p: dict) -> dict:
        b = {q["r"]: q for q in p["budgets"]}
        return {
            "network": p["network"],
            "provenance": (
                "structure and fitted coefficients from the network file"
                if p["params"] == "fitted"
                else "structure from the published DAG, semi-synthetic parameters"
            ),
            "source": p["source"],
            "query": f"{p['x']} -> {p['y']}",
            "n_k": p["n_k"],
            "n_k_loc": p["n_k_loc"],
            "tau_true": p["tau"],
            "estimate": p.get("estimate"),
            "estimate_ci": p.get("estimate_ci"),
            "committed_verdict": p["committed_verdict"],
            "I1_band": b.get(1, {}).get("band_seed0"),
            "I1_pop": b.get(1, {}).get("pop"),
            "I2_band": b.get(2, {}).get("band_seed0"),
            "I2_pop": b.get(2, {}).get("pop"),
            "r0_pop": r0_label(p["r0_pop"], p["r0_pop_status"]),
            "r0_band_seed0": r0_label(p["r0_band_seeds"][0]),
            "r0_band_range": r0_seed_range(p),
            "r_val": rval_label(p["r_val"], p["r_val_status"]),
            "breaking_claims": p["breaking_claims_seed0"],
            "claims_local": p["claims_local"],
            "informative": p["informative"],
            "frame_index": p["frame_index"],
        }

    for p in sorted((q for q in profiles if q["source"] == "declared"), key=lambda q: q["network"]):
        t2.append(t2_row(p))
    for net in PANEL_B:
        cand = [
            q
            for q in profiles
            if q["network"] == net
            and q["source"] == "frame"
            and q["condition"] == "A_ORACLE"
            and q["z"] is not None
        ]
        cand.sort(key=lambda q: (q["r0_pop"] == CENSORED, q["frame_index"]))
        t2 += [t2_row(q) for q in cand[:5]]
    excluded = [
        {
            "network": r["network"],
            "exposure": r["exposure"],
            "outcome": r["outcome"],
            "status": r["status"],
        }
        for r in declared_all
        if r["status"] != "ok"
    ]
    dump(
        out / "table2.json",
        {
            "rows": t2,
            "excluded_declared": excluded,
            "knowledge": "A_ORACLE (recovering set)",
            "n": 20000,
            "seed": 0,
        },
    )

    def iv(v: list | None) -> str:
        return f"[{v[0]:+.3f}, {v[1]:+.3f}]" if v else "-"

    t2md = [
        "# Table 2 — sensitivity reports on named analyses (measurement M1)",
        "",
        "Knowledge: the recovering set. Estimate and bands at n = 20000, seed 0; population "
        "interval in parentheses. r0: population / band at seed 0 / range of the band over "
        "20 seeds. r_val: the claim radius, depth 3, 300-subset budget. Declared queries use "
        "semi-synthetic parameters on the published structure; the real-coefficient rows "
        "are the first five frame rows with a committed set, finite population r0 first, "
        "then frame order.",
        "",
        "| network | params | query | K (local) | true effect | estimate [95 % CI] | I1 band "
        "(pop) | I2 band (pop) | r0 pop / band / range | r_val | first breaking claims |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in t2:
        est = (
            f"{r['estimate']:+.3f} {iv(r['estimate_ci'])}"
            if r["estimate"] is not None
            else r["committed_verdict"]
        )
        brk = (
            " and ".join(f"{a}->{b}" for a, b in r["breaking_claims"])
            if r["breaking_claims"]
            else "-"
        )
        t2md.append(
            f"| {r['network']} | {'fitted' if 'fitted' in r['provenance'] else 'semi-synth.'} | "
            f"{r['query']} | {r['n_k']} ({r['n_k_loc']}) | {r['tau_true']:+.3f} | {est} | "
            f"{iv(r['I1_band'])} ({iv(r['I1_pop'])}) | {iv(r['I2_band'])} ({iv(r['I2_pop'])}) | "
            f"{r['r0_pop']} / {r['r0_band_seed0']} / {r['r0_band_range']} | {r['r_val']} | "
            f"{brk} |"
        )
    t2md += ["", "Declared queries excluded:", ""]
    t2md += [
        f"- {e['network']}: {e['exposure']} on {e['outcome']}, {e['status']}" for e in excluded
    ]
    (out / "TABLE2_REPORTS.md").write_text("\n".join(t2md) + "\n")

    # ----------------------------------------------------------------- Figure 1
    applied = {
        "Acid_1996",
        "Didelez_2010",
        "Kampen_2014",
        "Polzer_2012",
        "Schipf_2010",
        "Sebastiani_2005",
        "Shrier_2008",
        "Thoemmes_2013",
    }

    def fig_check(p: dict) -> tuple[list[str], float]:
        fails = []
        r0p, r0b = p["r0_pop"], p["r0_band_seeds"][0]
        if r0p not in (2, 3):
            fails.append(f"population r0 = {r0_label(r0p, p['r0_pop_status'])}")
        if r0b not in (2, 3):
            fails.append(f"band r0 = {r0_label(r0b)}")
        jump = float("nan")
        b = {q["r"]: q for q in p["budgets"] if q.get("enumerated")}
        if r0b in (2, 3) and all(r in b for r in range(r0b + 1)):
            w = {r: b[r]["band_seed0"][1] - b[r]["band_seed0"][0] for r in b}
            if w[r0b - 1] > 1.5 * w[0]:
                fails.append(
                    f"no plateau: width at {r0b - 1} is {w[r0b - 1] / w[0]:.2f} x width at 0"
                )
            jump = w[r0b] / w[r0b - 1] if w[r0b - 1] > 0 else float("inf")
        return fails, jump

    order_groups = []
    for cond in ("A_ORACLE", "CTRL_b0"):
        order_groups.append(
            [q for q in profiles if q["condition"] == cond and q["network"] in PANEL_B]
        )
        order_groups.append(
            [q for q in profiles if q["condition"] == cond and q["source"] == "declared"]
        )
        order_groups.append(
            [
                q
                for q in profiles
                if q["condition"] == cond and q["network"] in applied and q["source"] == "frame"
            ]
        )
    chosen, misses = [], []
    for group in order_groups:
        scored = []
        for p in group:
            fails, jump = fig_check(p)
            (scored if not fails else misses).append((p, fails, jump))
        scored.sort(key=lambda t: -t[2] if math.isfinite(t[2]) else -1e9)
        for p, _, jump in scored:
            if len(chosen) < 6:
                chosen.append({**p, "jump": jump})
    near = sorted(
        (m for m in misses if len(m[1]) == 1),
        key=lambda m: (m[0]["network"] not in PANEL_B, m[0]["condition"] != "A_ORACLE"),
    )
    tier_rank = {"B": 0, "declared": 1, "applied": 2}
    cond_rank = {
        c: i
        for i, c in enumerate(("A_ORACLE", "CTRL_b0", "CTRL_b1", "CTRL_b2", "D_LLM", "D_LLM_72B"))
    }
    relaxed = []
    for p in profiles:
        group = (
            "B"
            if p["network"] in PANEL_B
            else "declared"
            if p["source"] == "declared"
            else "applied"
            if p["network"] in applied
            else None
        )
        if group is None or p["r0_pop"] not in (2, 3) or p["r0_band_seeds"][0] not in (2, 3):
            continue
        r0b = p["r0_band_seeds"][0]
        b = {q["r"]: q for q in p["budgets"] if q.get("enumerated")}
        if not all(r in b for r in range(r0b + 1)):
            continue
        w = {r: b[r]["band_seed0"][1] - b[r]["band_seed0"][0] for r in b}
        relaxed.append(
            {
                **p,
                "group": group,
                "plateau_ratio": w[r0b - 1] / w[0] if w[0] > 0 else float("inf"),
                "jump": w[r0b] / w[r0b - 1] if w[r0b - 1] > 0 else float("inf"),
                "r0_band_stable_share": float(np.mean(np.array(p["r0_band_seeds"]) == r0b)),
            }
        )
    relaxed.sort(
        key=lambda q: (
            tier_rank[q["group"]],
            q["n_false_local"] > 0,
            q["plateau_ratio"],
            cond_rank[q["condition"]],
            -q["jump"],
        )
    )
    fig = {
        "rule": "A_ORACLE then CTRL_b0; real-coefficient networks, then declared queries, then "
        "applied-DAG frame rows; population r0 and seed-0 band r0 in {2, 3}; band width "
        "at r0 - 1 at most 1.5 x width at 0; ranked by width(r0) / width(r0 - 1).",
        "n": 20000,
        "seed": 0,
        "candidates": chosen,
        "near_misses": [{**m[0], "failed": m[1]} for m in near[: max(0, 6 - len(chosen)) + 6]],
        "relaxed_rule": "NOT the pre-registered rule, added after it returned no candidate: any "
        "condition, population r0 and seed-0 band r0 in {2, 3}, plateau "
        "condition dropped; real-coefficient networks, then declared queries, "
        "then applied-DAG frame rows; rows whose local claims are all true first; "
        "flattest band before r0 first; then truthful, controlled, elicited.",
        "relaxed_candidates": relaxed[:6],
        "relaxed_total": len(relaxed),
    }
    dump(out / "fig1_candidates.json", fig)

    summary = {
        "informative_stratum": strat,
        "direction_dispute": dispute,
        "r0": r0,
        "predictions": preds,
        "censoring": dict(sorted(cens.items())),
        "r_val_sources": dict(sorted(rval_src.items())),
        "ledger_verdict_mismatches": mismatch,
        "ch_rows_by_reason": dict(sorted(ch_counts.items())),
        "plain_kind_informative": dict(sorted(plain_kinds.items())),
        "ch_population_coverage_defined_rows": ch_pop.get("procedures", {}).get("CH"),
        "panel": panel,
        "bootstrap": {"B": B_BOOT, "seed": BOOT_SEED, "cluster": "network"},
        "table2_rows": len(t2),
        "fig1_candidates": len(chosen),
        "m1b_predictions": m1b_preds,
        "m1b": {k: v for k, v in m1b.items() if k != "policy_cells"},
    }
    dump(out / "summary.json", summary)
    print(json.dumps(clean(preds), indent=1))
    print(json.dumps(clean(m1b_preds), indent=1))


if __name__ == "__main__":
    main()
