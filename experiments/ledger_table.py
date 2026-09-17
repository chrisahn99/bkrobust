"""The certificate ledger: one row per knowledge arm, deployment above the rule, controls below.

This is the table the paper is built around. Every number in it resolves to a
row of ``results/ledger/rows.csv`` and every rate carries a network cluster
bootstrap, because the network is the unit of evidence and the design effect
on this corpus is near fourteen.

Block 1, what the analyst sees, computed from observables: how many frame rows
the arm can certify at all and why the rest cannot be, the size of the asserted
set against the size of its Meek closure, the certificate in claims and its
share at one, the hop radius beside it, the fraction of single retractions that
break the set, and how much the two replicates disagree.

Block 2, the truth revealed: the fraction of asserted claims that were false,
the fraction of orientable claims never made, how often the committed set was
invalid at the truth, and the four-way audit bucket for each instrument on the
same rows. The dangerous rate is the over-certification rate; the rule-of-three
floor is printed beside it so no cell certifies below what its size permits.

Writes ``results/ledger/TABLE4.md`` and ``results/ledger/table4.json``.

    python experiments/ledger_table.py
"""

from __future__ import annotations

import collections
import csv
import json
import math
from collections.abc import Callable
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ROWS = ROOT / "results" / "ledger" / "rows.csv"
KNOWLEDGE = ROOT / "results" / "elicit" / "knowledge.json"
OUT = ROOT / "results" / "ledger"

ARM_LABEL = {
    "D_LLM": "language model, real names",
    "D_LLM_INSTR": "same model, second presentation order",
    "D_LLM_SRC": "second model family, real names",
    "D_SCRAMBLED": "same model, scrambled names",
    "D_LLM_72B": "same family, 72B, real names",
    "D_LLM_72B_INSTR": "same family, 72B, second presentation order",
    "D_SCRAMBLED_72B": "same family, 72B, scrambled names",
    "D_LLM_32B": "same family, 32B, real names",
    "D_LLM_SRC_70B": "second family, 70B, real names",
    "D_LLM_GEMMA_27B": "third family, 27B, real names",
    "D_LLM_PLUS_ANC": "language model, real names, plus reduced ancestral claims",
    "D_RAND": "chance level, matched on the model's edges",
    "D_DEGEN": "no knowledge asserted",
    "A_TRUE": "truthful orientations on the model's edges",
    "A_ORACLE": "greedy recovering set (fully oriented graph)",
}
DEPLOYMENT = (
    "D_LLM",
    "D_LLM_INSTR",
    "D_LLM_SRC",
    "D_SCRAMBLED",
    "D_LLM_72B",
    "D_LLM_72B_INSTR",
    "D_SCRAMBLED_72B",
    "D_LLM_32B",
    "D_LLM_SRC_70B",
    "D_LLM_GEMMA_27B",
    "D_LLM_PLUS_ANC",
    "D_RAND",
    "D_DEGEN",
)
AUDIT = ("A_TRUE", "A_ORACLE")
RUNGS = ("B2_free_rule", "B4_kg0", "B5_k", "B6_hop", "B7_claim")
MAX_DEPTH = 3  # the claim instrument's search depth; certified moves are capped at it
B = 4000
SEED = 20260913


def boot(
    by_net: dict[str, list[float]], reps: int = B, seed: int = SEED
) -> tuple[float, float, float]:
    """Mean and a percentile interval, resampling networks."""
    flat = [v for vs in by_net.values() for v in vs]
    if not flat:
        return math.nan, math.nan, math.nan
    names = sorted(by_net)
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(reps):
        pick = rng.integers(0, len(names), size=len(names))
        vals = [v for i in pick for v in by_net[names[int(i)]]]
        if vals:
            draws.append(sum(vals) / len(vals))
    draws.sort()
    return sum(flat) / len(flat), draws[int(0.025 * len(draws))], draws[int(0.975 * len(draws)) - 1]


def certified_moves(r: dict) -> int:
    """Safe moves the claim instrument licenses on a row, capped at the search depth.

    An exact radius ``r`` licenses ``r - 1``. A search censored at depth ``d``
    checked every smaller subset and licenses ``d - 1``; one that finished its
    depth licenses ``depth``; one that exhausted every subset without a failure
    licenses everything and is capped. The cap keeps the mean comparable across
    arms whose censoring differs.
    """
    if r["r_claim"] not in ("", "None"):
        return min(int(r["r_claim"]), MAX_DEPTH + 1) - 1
    st = r["r_claim_status"]
    if st in ("unreached", "no_retractable_edges"):
        return MAX_DEPTH
    if st.startswith("censored_at_depth_"):
        return int(st.rsplit("_", 1)[1]) - 1
    if st.startswith("gt_"):
        return min(int(st.rsplit("_", 1)[1]), MAX_DEPTH)
    raise ValueError(st)


def leverage(r: dict) -> float | None:
    """Hop radius over claim radius where both are exact and finite."""
    if r["hop_exact"] != "1" or r["r_hop"] in ("", "-1") or r["r_claim"] in ("", "None"):
        return None
    return int(r["r_hop"]) / int(r["r_claim"])


def fmt_rate(m: float, lo: float, hi: float) -> str:
    """A rate with its interval."""
    if math.isnan(m):
        return "—"
    return f"{m:.3f} [{lo:.3f}, {hi:.3f}]"


def summarise(rows: list[dict]) -> dict:
    """Everything the ledger prints for one arm, from its rows."""
    total = len(rows)
    ok = [r for r in rows if r["status"] == "ok"]
    nets = {r["network"] for r in rows}
    status = collections.Counter(r["status"] for r in rows)

    def by_net(
        pred: Callable[[dict], object], subset: list[dict] | None = None
    ) -> dict[str, list[float]]:
        subset = ok if subset is None else subset
        d: dict[str, list[float]] = collections.defaultdict(list)
        for r in subset:
            d[r["network"]].append(float(pred(r)))
        return d

    out: dict = {
        "rows": total,
        "networks": len(nets),
        "certifiable": len(ok),
        "status": dict(status),
        "n_k_median": int(np.median([int(r["n_k"]) for r in rows])) if rows else None,
        "k_g0_median": int(np.median([int(r["k_g0"]) for r in ok])) if ok else None,
    }
    if not ok:
        return out
    claim_int = [int(r["r_claim"]) for r in ok if r["r_claim"] not in ("", "None")]
    out["r_claim_share_at_1"] = boot(by_net(lambda r: r["r_claim"] == "1"))
    out["r_claim_distribution"] = dict(
        sorted(
            collections.Counter(
                r["r_claim"] if r["r_claim"] not in ("", "None") else r["r_claim_status"]
                for r in ok
            ).items()
        )
    )
    out["safe_moves_mean"] = float(np.mean([c - 1 for c in claim_int])) if claim_int else None
    out["certified_moves"] = boot(by_net(certified_moves))
    levs = [v for v in (leverage(r) for r in ok) if v is not None]
    out["leverage"] = {
        "rows": len(levs),
        "median": float(np.median(levs)) if levs else None,
        "max": float(max(levs)) if levs else None,
        "share_above_1": sum(1 for v in levs if v > 1) / len(levs) if levs else None,
        "share_above_5": sum(1 for v in levs if v > 5) / len(levs) if levs else None,
    }
    out["hop_unreached"] = sum(1 for r in ok if r["r_hop"] == "-1")
    out["hop_no_retractable"] = sum(1 for r in ok if r["hop_method"] == "no_retractable_edges")
    out["hop_inexact"] = sum(1 for r in ok if r["hop_exact"] == "0")
    phis = [float(r["phi_1"]) for r in ok if r["phi_1"] not in ("", "None")]
    out["phi_1_median"] = float(np.median(phis)) if phis else None
    # the truth
    out["false_claim_fraction"] = boot(
        by_net(
            lambda r: int(r["n_wrong_claims"]) / max(int(r["n_k"]), 1),
            [r for r in ok if int(r["n_k"]) > 0],
        )
    )
    out["omitted_mean"] = float(np.mean([int(r["n_omitted"]) for r in ok]))
    out["z_invalid_at_truth"] = boot(by_net(lambda r: r["z_valid_at_truth"] == "0"))
    out["severity"] = dict(collections.Counter(r["severity"] for r in ok))
    n_invalid = sum(1 for r in ok if r["z_valid_at_truth"] == "0")
    out["rule_of_three_floor"] = round(3 / len(ok), 4) if ok else None
    ladder = {}
    for rung in RUNGS:
        bk = collections.Counter(r[f"bucket_{rung}"] for r in ok)
        inv = bk["dangerous"] + bk["safe"]
        val = bk["held"] + bk["slack"]
        m, lo, hi = boot(by_net(lambda r, rung=rung: r[f"bucket_{rung}"] == "dangerous"))
        ladder[rung] = {
            "dangerous": (m, lo, hi),
            "detection": bk["safe"] / inv if inv else math.nan,
            "false_alarm": bk["slack"] / val if val else math.nan,
            "buckets": dict(bk),
        }
    out["ladder"] = ladder
    out["truth_invalid_rows"] = n_invalid
    return out


def main() -> None:
    """Assemble Table 4 from the ledger rows."""
    rows = list(csv.DictReader(ROWS.open()))
    know = json.loads(KNOWLEDGE.read_text()) if KNOWLEDGE.exists() else {}
    by_arm: dict[str, list[dict]] = collections.defaultdict(list)
    for r in rows:
        by_arm[r["arm"]].append(r)
    table = {a: summarise(by_arm[a]) for a in [*DEPLOYMENT, *AUDIT] if a in by_arm}

    # replicate disagreement, per network, between D_LLM and each replicate
    def disagreement(a: str, b: str) -> tuple[float, float, float] | None:
        if a not in know or b not in know:
            return None
        d: dict[str, list[float]] = {}
        for net, ka in know[a]["networks"].items():
            kb = know[b]["networks"].get(net)
            if not kb:
                continue
            sa = {tuple(e) for e in ka["k"]}
            sb = {tuple(e) for e in kb["k"]}
            both = {frozenset(e) for e in sa} & {frozenset(e) for e in sb}
            if not both:
                continue
            flipped = sum(
                1
                for e in both
                if (tuple(sorted(e)) in {tuple(sorted(x)) for x in sa})
                and (
                    (next(x for x in sa if frozenset(x) == e))
                    != (next(x for x in sb if frozenset(x) == e))
                )
            )
            d[net] = [flipped / len(both)]
        return boot(d) if d else None

    deltas = {
        "delta_instrument": disagreement("D_LLM", "D_LLM_INSTR"),
        "delta_source": disagreement("D_LLM", "D_LLM_SRC"),
        "delta_scrambled": disagreement("D_LLM", "D_SCRAMBLED"),
        "delta_scale_7B_72B": disagreement("D_LLM", "D_LLM_72B"),
        "delta_instrument_72B": disagreement("D_LLM_72B", "D_LLM_72B_INSTR"),
        "delta_scrambled_72B": disagreement("D_LLM_72B", "D_SCRAMBLED_72B"),
        "delta_source_72B_70B": disagreement("D_LLM_72B", "D_LLM_SRC_70B"),
        "delta_source_72B_gemma": disagreement("D_LLM_72B", "D_LLM_GEMMA_27B"),
    }
    # elicitation-level statistics from the knowledge file
    elicit = {}
    for a, entry in know.items():
        nets = entry["networks"]
        elicit[a] = {
            "model": entry["model"],
            "family": entry.get("family"),
            "asked": sum(v["n_asked"] for v in nets.values()),
            "asserted": sum(v["n_asserted"] for v in nets.values()),
            "declined": sum(v["declined"] for v in nets.values()),
            "not_reached": sum(v["not_reached"] for v in nets.values()),
            "meek_inconsistent_networks": sum(
                1 for v in nets.values() if not v["meek_consistent_pre"]
            ),
            "dropped_for_consistency": sum(v["dropped_for_consistency"] for v in nets.values()),
            "control_right": sum(
                (v["compelled_control"] or {}).get("right", 0) for v in nets.values()
            ),
            "control_wrong": sum(
                (v["compelled_control"] or {}).get("wrong", 0) for v in nets.values()
            ),
        }

    # paired contrasts: the same networks under two arms, difference of network means
    def paired(a: str, b: str, stat: Callable[[dict], float | None]) -> dict | None:
        if a not in by_arm or b not in by_arm:
            return None
        means: dict[str, dict[str, float]] = {a: {}, b: {}}
        for arm in (a, b):
            per: dict[str, list[float]] = collections.defaultdict(list)
            for r in by_arm[arm]:
                if r["status"] == "ok" and (v := stat(r)) is not None:
                    per[r["network"]].append(v)
            means[arm] = {n: sum(v) / len(v) for n, v in per.items()}
        shared = sorted(set(means[a]) & set(means[b]))
        if not shared:
            return None
        diff = {n: [means[a][n] - means[b][n]] for n in shared}
        m, lo, hi = boot(diff)
        return {"networks": len(shared), "mean_diff": m, "ci": [lo, hi]}

    def false_fraction(r: dict) -> float | None:
        return int(r["n_wrong_claims"]) / int(r["n_k"]) if int(r["n_k"]) > 0 else None

    contrasts = {
        "false_claims D_RAND - D_LLM": paired("D_RAND", "D_LLM", false_fraction),
        "invalid_at_truth D_RAND - D_LLM": paired(
            "D_RAND", "D_LLM", lambda r: float(r["z_valid_at_truth"] == "0")
        ),
        "certified_moves A_TRUE - D_LLM": paired(
            "A_TRUE", "D_LLM", lambda r: float(certified_moves(r))
        ),
        "certified_moves D_LLM - D_RAND": paired(
            "D_LLM", "D_RAND", lambda r: float(certified_moves(r))
        ),
    }

    # the paired instrument ladder: claim certificate against the free rule and the
    # hop radius on the same rows, network-paired, with and without the textbook tier
    def instrument_rates(rows_: list[dict], rung: str) -> dict[str, dict[str, float]]:
        per: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
        for r in rows_:
            per[r["network"]][r[f"bucket_{rung}"]] += 1
        out: dict[str, dict[str, float]] = {}
        for net, c in per.items():
            n = sum(c.values())
            valid = c["held"] + c["slack"]
            out[net] = {"dangerous": c["dangerous"] / n}
            if valid:
                out[net]["false_alarm"] = c["slack"] / valid
        return out

    def verdict(lo: float, hi: float, margin: float = 0.05) -> str:
        """Sharpness verdict on a paired interval; lower false alarm is better."""
        if hi < 0:
            return "radius wins"
        if lo > -margin and hi < margin:
            return "equivalent"
        if lo > 0:
            return "radius loses"
        return "UNDERPOWERED"

    def validity_verdict(claim_rows: int, other_rows: int) -> str:
        """Validity is absolute: a certificate is safe when it over-certifies on no row."""
        if claim_rows == 0 and other_rows == 0:
            return "both safe"
        if claim_rows == 0:
            return f"claim safe, other dangerous on {other_rows}"
        return f"claim dangerous on {claim_rows}"

    paired_ladder: dict[str, dict] = {}
    for arm in DEPLOYMENT:
        ok_rows = [r for r in by_arm.get(arm, []) if r["status"] == "ok"]
        if not ok_rows:
            continue
        for tier_label, subset in (
            ("all", ok_rows),
            ("without T0", [r for r in ok_rows if r["tier"] != "T0"]),
        ):
            claim = instrument_rates(subset, "B7_claim")
            for other_rung, other_label in (("B2_free_rule", "free rule"), ("B6_hop", "hop")):
                other = instrument_rates(subset, other_rung)
                entry: dict = {"networks": 0}
                for key in ("dangerous", "false_alarm"):
                    shared = [
                        n for n in claim if n in other and key in claim[n] and key in other[n]
                    ]
                    if not shared:
                        entry[key] = None
                        continue
                    m, lo, hi = boot({n: [claim[n][key] - other[n][key]] for n in shared})
                    if key == "dangerous":
                        v = validity_verdict(
                            sum(1 for r in subset if r["bucket_B7_claim"] == "dangerous"),
                            sum(1 for r in subset if r[f"bucket_{other_rung}"] == "dangerous"),
                        )
                    else:
                        v = verdict(lo, hi)
                    entry[key] = {"mean_diff": m, "ci": [lo, hi], "verdict": v}
                    entry["networks"] = max(entry["networks"], len(shared))
                paired_ladder[f"{arm} | {tier_label} | claim - {other_label}"] = entry

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "table4.json").write_text(
        json.dumps(
            {
                "arms": table,
                "deltas": deltas,
                "elicitation": elicit,
                "contrasts": contrasts,
                "paired_ladder": paired_ladder,
            },
            indent=1,
            default=str,
        )
    )

    lines: list[str] = []
    add = lines.append

    def row(cells: list[object]) -> None:
        add("| " + " | ".join(str(c) for c in cells) + " |")

    def sep(n: int) -> None:
        add("|" + "---|" * n)

    def num(v: object, nd: int = 3) -> str:
        return "—" if v is None else (f"{v:.{nd}f}" if isinstance(v, float) else str(v))

    add("# Table 4 — the certificate ledger")
    add("")
    add("Generated by `experiments/ledger_table.py` from `results/ledger/rows.csv`. Every rate")
    add(f"carries a network cluster bootstrap (B = {B}); the network is the unit. Deployment arms")
    add("above the rule assert nothing read off the truth. Their CPDAG is the true DAG's, so every")
    add("row is stamped `via_Chat_only`, and none carries a deployment claim until an estimated")
    add("CPDAG replaces it. `query reversed` counts rows where the analyst's own graph says the")
    add("treatment cannot cause the outcome; the frame drew the pair from a CPDAG in which it was")
    add("possible and the knowledge oriented it away.")
    add("")
    add("## Block 1 — what the analyst sees")
    add("")
    head1 = [
        "arm",
        "rows",
        "nets",
        "certifiable",
        "query reversed",
        "not amenable",
        "inconsistent",
        "median claims",
        "median closure",
        "share `r_claim = 1`",
        "certified moves (cap 3)",
        "leverage median / max",
        "hop unreached",
        "`phi_1` median",
    ]
    row(head1)
    sep(len(head1))
    for group in (DEPLOYMENT, AUDIT):
        for a in group:
            if a not in table:
                continue
            t = table[a]
            st = t["status"]
            share = fmt_rate(*t["r_claim_share_at_1"]) if "r_claim_share_at_1" in t else "—"
            row(
                [
                    f"`{a}` {ARM_LABEL.get(a, '')}",
                    t["rows"],
                    t["networks"],
                    t["certifiable"],
                    st.get("y_not_possible_descendant", 0),
                    st.get("not_amenable", 0),
                    st.get("knowledge_inconsistent", 0),
                    t["n_k_median"],
                    num(t.get("k_g0_median")),
                    share,
                    fmt_rate(*t["certified_moves"]) if "certified_moves" in t else "—",
                    (
                        f"{num(t['leverage']['median'], 2)} / {num(t['leverage']['max'], 1)}"
                        if t.get("leverage", {}).get("rows")
                        else "—"
                    ),
                    num(t.get("hop_unreached")),
                    num(t.get("phi_1_median")),
                ]
            )
        if group is DEPLOYMENT:
            row(
                ["**— audit controls below this rule read the truth to build their knowledge —**"]
                + [""] * (len(head1) - 1)
            )
    add("")
    add("## Replicate disagreement, fraction of shared edges oriented differently")
    add("")
    for k, v in deltas.items():
        add(f"- `{k}`: {fmt_rate(*v) if v else 'ABSENT'}")
    add("")
    add("## Paired contrasts, difference of network means on shared networks")
    add("")
    for k, v in contrasts.items():
        if v:
            ci = f"[{v['ci'][0]:+.3f}, {v['ci'][1]:+.3f}]"
            add(f"- `{k}`: {v['mean_diff']:+.3f} {ci} on {v['networks']} networks")
        else:
            add(f"- `{k}`: ABSENT")
    add("")
    add("## Elicitation, per supplier")
    add("")
    head_e = [
        "condition",
        "model",
        "asked",
        "asserted",
        "declined",
        "not reached",
        "Meek-inconsistent nets",
        "dropped",
        "control right",
        "control wrong",
    ]
    row(head_e)
    sep(len(head_e))
    for a, e in elicit.items():
        row(
            [
                f"`{a}`",
                e["model"],
                e["asked"],
                e["asserted"],
                e["declined"],
                e["not_reached"],
                e["meek_inconsistent_networks"],
                e["dropped_for_consistency"],
                e["control_right"],
                e["control_wrong"],
            ]
        )
    add("")
    add("## Block 2 — the truth revealed, on the certifiable rows")
    add("")
    head2 = [
        "arm",
        "false claims / asserted",
        "mean omitted",
        "set invalid at truth",
        "severity",
        "rule-of-three floor",
    ]
    row(head2)
    sep(len(head2))
    for a in [*DEPLOYMENT, *AUDIT]:
        if a not in table or "false_claim_fraction" not in table[a]:
            continue
        t = table[a]
        row(
            [
                f"`{a}`",
                fmt_rate(*t["false_claim_fraction"]),
                f"{t['omitted_mean']:.2f}",
                fmt_rate(*t["z_invalid_at_truth"]),
                t["severity"],
                t["rule_of_three_floor"],
            ]
        )
    add("")
    add("## The ladder, per arm: over-certification (dangerous), detection, false alarm")
    add("")
    head3 = ["arm", "instrument", "dangerous [CI]", "detection", "false alarm", "buckets"]
    row(head3)
    sep(len(head3))
    for a in [*DEPLOYMENT, *AUDIT]:
        if a not in table or "ladder" not in table[a]:
            continue
        for rung, v in table[a]["ladder"].items():
            m, lo, hi = v["dangerous"]
            det = "—" if math.isnan(v["detection"]) else f"{v['detection']:.3f}"
            fa = "—" if math.isnan(v["false_alarm"]) else f"{v['false_alarm']:.3f}"
            row([f"`{a}`", f"`{rung}`", fmt_rate(m, lo, hi), det, fa, v["buckets"]])
    add("")
    add("## The paired instrument ladder: the claim certificate against the free rule and the hop")
    add("")
    add("Difference of network means, claim minus the other instrument, on the same certifiable")
    add("rows. Lower is better on both columns. Validity is absolute: a certificate is safe when")
    add("it over-certifies on no row, and the column says which of the two does. Sharpness is a")
    add("paired verdict: `radius wins` when the interval is below zero, `equivalent` when it lies")
    add("inside ±0.05, `radius loses` when above zero, and `UNDERPOWERED` otherwise. `without T0`")
    add("drops the textbook diagrams.")
    add("")
    head4 = [
        "arm",
        "rows",
        "contrast",
        "dangerous, claim minus other",
        "false alarm, claim minus other",
        "nets",
    ]
    row(head4)
    sep(len(head4))
    for key, e in paired_ladder.items():
        arm, tier_label, contrast = (k.strip() for k in key.split("|"))

        def cell(v: dict | None) -> str:
            if not v:
                return "—"
            return f"{v['mean_diff']:+.3f} [{v['ci'][0]:+.3f}, {v['ci'][1]:+.3f}] {v['verdict']}"

        row(
            [
                f"`{arm}`",
                tier_label,
                contrast,
                cell(e["dangerous"]),
                cell(e["false_alarm"]),
                e["networks"],
            ]
        )
    (OUT / "TABLE4.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
