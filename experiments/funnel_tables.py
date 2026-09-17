"""Table 0, the applicability funnel, and Table 3, the status composition.

Table 0 says what fraction of a real network the framework speaks about, per
substrate tier: networks parsed against networks yielding a frame, nodes inside
a chain component, orientable edges, candidate pairs and the sampled rows, and,
beside every certifiable rate, the number of networks, the design effect and the
effective sample size of the network cluster. Table 3 is the fixed denominator:
one row per tier and knowledge arm, every status a column, so an arm's
population loss is a reported outcome rather than a silent shrinkage. Both also
print the two things the protocol asks for twice: the maximum radius with and
without the textbook tier, and the number of certifiable rows whose analyst
graph keeps at least one undirected edge, which is the Lever 0 admission count.

Writes ``results/ledger/TABLE0_FUNNEL.md``, ``TABLE3_STATUS.md`` and ``funnel.json``.

    python experiments/funnel_tables.py
"""

from __future__ import annotations

import collections
import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "results" / "ledger"
FRAME = ROOT / "results" / "frame"
TIER = {
    **{n: "T0" for n in ("confounding", "mediator", "paths", "M-bias")},
    **{
        n: "T2"
        for n in (
            "Acid_1996",
            "Didelez_2010",
            "Kampen_2014",
            "Polzer_2012",
            "Schipf_2010",
            "Sebastiani_2005",
            "Shrier_2008",
            "Thoemmes_2013",
        )
    },
    **{n: "T3" for n in ("arth150", "ecoli70", "magic-irri", "magic-niab")},
}
TIERS = ("T0", "T1", "T2", "T3")
STATUSES = (
    "ok",
    "y_not_possible_descendant",
    "not_amenable",
    "no_valid_set",
    "knowledge_inconsistent",
)


def design_effect(groups: list[list[float]]) -> tuple[float, float, float]:
    """One-way ANOVA intraclass correlation, design effect and effective size."""
    sizes = [len(g) for g in groups if g]
    if len(sizes) < 2:
        return 0.0, 1.0, float(sum(sizes))
    allv = np.concatenate([np.array(g) for g in groups if g])
    grand = allv.mean()
    n = len(allv)
    k = len(sizes)
    ssb = sum(len(g) * (np.mean(g) - grand) ** 2 for g in groups if g)
    ssw = sum(((np.array(g) - np.mean(g)) ** 2).sum() for g in groups if g)
    msb = ssb / (k - 1)
    msw = ssw / max(n - k, 1)
    n0 = (n - sum(s * s for s in sizes) / n) / (k - 1)
    icc = (msb - msw) / (msb + (n0 - 1) * msw) if (msb + (n0 - 1) * msw) > 0 else 0.0
    icc = max(0.0, min(1.0, icc))
    mbar = n / k
    de = 1 + (mbar - 1) * icc
    return float(icc), float(de), float(n / de)


def main() -> None:
    """Assemble the funnel and the status composition from the frame and the ledger."""
    per_net = list(csv.DictReader((FRAME / "per_network.csv").open()))
    rows = list(csv.DictReader((LEDGER / "rows.csv").open()))
    arms = sorted({r["arm"] for r in rows}, key=lambda a: (a.startswith("A_"), a))
    tier_of = {r["network"]: TIER.get(r["network"], "T1") for r in per_net}
    for r in rows:
        tier_of.setdefault(r["network"], TIER.get(r["network"], "T1"))

    # ---- Table 0
    funnel: dict = {}
    swept = {r["network"] for r in rows}
    for t in TIERS:
        nets = [p for p in per_net if tier_of[p["network"]] == t]
        yielding = [p for p in nets if int(p["frame_pairs"]) > 0]
        funnel[t] = {
            "networks_parsed": len(nets),
            "networks_with_frame": len(yielding),
            "networks_swept": len([p for p in nets if p["network"] in swept]),
            "nodes": sum(int(p["nodes"]) for p in nets),
            "undirected_edges": sum(int(p["undirected_edges"]) for p in nets),
            "components_ge2": sum(int(p["components_ge2"]) for p in nets),
            "candidate_treatments": sum(int(p["candidate_treatments"]) for p in nets),
            "frame_pairs": sum(int(p["frame_pairs"]) for p in nets),
            "sampled_pairs": sum(int(p["sampled_pairs"]) for p in nets),
            "committed_admissible": sum(int(p["committed_admissible"]) for p in nets),
        }
    # certifiable rate per tier under the primary supplier, with the cluster design effect
    for t in TIERS:
        by_net: dict[str, list[float]] = collections.defaultdict(list)
        for r in rows:
            if r["arm"] == "D_LLM" and tier_of[r["network"]] == t:
                by_net[r["network"]].append(float(r["status"] == "ok"))
        groups = list(by_net.values())
        icc, de, neff = design_effect(groups)
        n = sum(len(g) for g in groups)
        funnel[t].update(
            {
                "rows_primary": n,
                "certifiable_primary": int(sum(sum(g) for g in groups)),
                "icc": round(icc, 4),
                "design_effect": round(de, 2),
                "n_eff": round(neff, 1),
            }
        )

    # ---- Table 3, statuses per tier x arm; max radius twice; Lever 0 count
    status: dict[str, dict[str, collections.Counter]] = collections.defaultdict(
        lambda: collections.defaultdict(collections.Counter)
    )
    for r in rows:
        status[r["arm"]][tier_of[r["network"]]][r["status"]] += 1
    radius: dict = {}
    for a in arms:
        ok = [r for r in rows if r["arm"] == a and r["status"] == "ok"]
        okx = [r for r in ok if tier_of[r["network"]] != "T0"]

        def mx(rs: list[dict], col: str) -> int | None:
            vals = [int(r[col]) for r in rs if r[col] not in ("", "None", "-1")]
            return max(vals) if vals else None

        radius[a] = {
            "max_r_claim_all": mx(ok, "r_claim"),
            "max_r_claim_without_T0": mx(okx, "r_claim"),
            "max_r_hop_all": mx(ok, "r_hop"),
            "max_r_hop_without_T0": mx(okx, "r_hop"),
            "r_claim_ge3_all": sum(
                1 for r in ok if r["r_claim"] not in ("", "None") and int(r["r_claim"]) >= 3
            ),
            "r_claim_ge3_without_T0": sum(
                1 for r in okx if r["r_claim"] not in ("", "None") and int(r["r_claim"]) >= 3
            ),
            "certifiable": len(ok),
            "g0_keeps_undirected": sum(
                1
                for r in ok
                if r["g0_undirected"] not in ("", "None") and int(r["g0_undirected"]) >= 1
            ),
        }
    (LEDGER / "funnel.json").write_text(
        json.dumps(
            {
                "table0": funnel,
                "table3": {a: {t: dict(c) for t, c in v.items()} for a, v in status.items()},
                "radius": radius,
            },
            indent=1,
        )
    )

    lines = ["# Table 0 — the applicability funnel", ""]
    lines.append(
        "Per substrate tier: what the frame contains and what the primary supplier can certify. The"
        " intraclass correlation, design effect and effective sample size are for the certifiable"
        " indicator with the network as the cluster. Networks above 500 nodes and the one whose"
        " chain components all exceed the questionnaire cap are in the frame but not swept."
    )
    lines.append("")
    head = [
        "tier",
        "networks parsed / with frame / swept",
        "nodes",
        "undirected edges",
        "components ≥ 2",
        "candidate treatments",
        "candidate pairs",
        "sampled rows",
        "committed pairs (old gate)",
        "certifiable (primary)",
        "ICC / design effect / n_eff",
    ]
    lines.append("| " + " | ".join(head) + " |")
    lines.append("|" + "---|" * len(head))
    for t in TIERS:
        f = funnel[t]
        cells = [
            t,
            f"{f['networks_parsed']} / {f['networks_with_frame']} / {f['networks_swept']}",
            f["nodes"],
            f["undirected_edges"],
            f["components_ge2"],
            f["candidate_treatments"],
            f["frame_pairs"],
            f["sampled_pairs"],
            f["committed_admissible"],
            f"{f['certifiable_primary']} of {f['rows_primary']}",
            f"{f['icc']} / {f['design_effect']} / {f['n_eff']}",
        ]
        lines.append("| " + " | ".join(str(c) for c in cells) + " |")
    lines += ["", "## The radius printed twice, and the Lever 0 admission count", ""]
    head2 = [
        "arm",
        "certifiable",
        "analyst graph keeps an undirected edge",
        "max r_claim all / without T0",
        "r_claim ≥ 3 all / without T0",
        "max r_hop all / without T0",
    ]
    lines.append("| " + " | ".join(head2) + " |")
    lines.append("|" + "---|" * len(head2))
    for a in arms:
        v = radius[a]
        cells = [
            f"`{a}`",
            v["certifiable"],
            v["g0_keeps_undirected"],
            f"{v['max_r_claim_all']} / {v['max_r_claim_without_T0']}",
            f"{v['r_claim_ge3_all']} / {v['r_claim_ge3_without_T0']}",
            f"{v['max_r_hop_all']} / {v['max_r_hop_without_T0']}",
        ]
        lines.append("| " + " | ".join(str(c) for c in cells) + " |")
    (LEDGER / "TABLE0_FUNNEL.md").write_text("\n".join(lines) + "\n")

    lines = ["# Table 3 — the status composition, per tier and knowledge arm", ""]
    lines.append(
        "The fixed denominator. Every arm reports the same rows; the columns say why a row could"
        " not be certified. `ok` rows carry a certificate."
    )
    lines.append("")
    head3 = ["arm", "tier", "rows", *STATUSES]
    lines.append("| " + " | ".join(head3) + " |")
    lines.append("|" + "---|" * len(head3))
    for a in arms:
        for t in TIERS:
            c = status[a].get(t)
            if not c:
                continue
            lines.append(
                f"| `{a}` | {t} | {sum(c.values())} | "
                + " | ".join(str(c.get(s, 0)) for s in STATUSES)
                + " |"
            )
    (LEDGER / "TABLE3_STATUS.md").write_text("\n".join(lines) + "\n")
    print((LEDGER / "TABLE0_FUNNEL.md").read_text())


if __name__ == "__main__":
    main()
