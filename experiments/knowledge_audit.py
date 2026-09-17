"""Stage 7, the knowledge-level audit: how wrong each supplier was, and where.

This reads the frozen knowledge and the true DAGs and nothing from the sweep,
so it answers the questions about the suppliers themselves. Commission is the
fraction of asserted claims whose orientation is false at the truth; omission
is the fraction of asked pairs the supplier did not assert, whether by
declining or by never reaching them. The two are separate axes and are never
summed. Both are reported per network and pooled with a network cluster
bootstrap, and split by substrate tier.

The compelled-edge control block is scored against what the data compel, which
needs no truth. Whether its error rate lies inside the interval of the
supplier's error rate on open edges, which does need the truth, is what decides
whether the truth-free bridge the protocol wants for the no-truth tier holds on
this corpus.

Replicate disagreement is the fraction of edges two conditions both asserted
that they oriented differently: the same model on a second presentation order
is instrument variance, a second family is source variance, and scrambled
names against real names is the recitation probe.

Writes ``results/elicit/knowledge_audit.json`` and prints the table.

    python experiments/knowledge_audit.py
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from bkrobust.benchmarks.describe import parse_file, to_mpdag  # noqa: E402

MODELS = ROOT / "results" / "axisa3" / "networks" / "example_models"
KNOWLEDGE = ROOT / "results" / "elicit" / "knowledge.json"
OUT = ROOT / "results" / "elicit"
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
B = 4000


def boot(
    by_net: dict[str, float], reps: int = B, seed: int = 20260913
) -> tuple[float, float, float]:
    """Mean of per-network values with a percentile interval over networks."""
    names = sorted(by_net)
    if not names:
        return math.nan, math.nan, math.nan
    vals = np.array([by_net[n] for n in names])
    rng = np.random.default_rng(seed)
    draws = np.sort(
        [vals[rng.integers(0, len(names), size=len(names))].mean() for _ in range(reps)]
    )
    return float(vals.mean()), float(draws[int(0.025 * reps)]), float(draws[int(0.975 * reps) - 1])


def main() -> None:
    """Score every supplier's knowledge against the truth, per network and pooled."""
    know = json.loads(KNOWLEDGE.read_text())
    dags = {}
    for path in sorted(MODELS.iterdir()):
        parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
        if not parsed.bidirected:
            dags[parsed.name] = to_mpdag(parsed)

    per: dict[str, dict[str, dict]] = {}
    for cond, entry in know.items():
        per[cond] = {}
        for net, v in entry["networks"].items():
            dag = dags[net]
            k = [tuple(e) for e in v["k"]]
            false = sum(1 for a, b in k if not dag.is_directed_edge(a, b))
            ctl = v.get("compelled_control") or {}
            per[cond][net] = {
                "tier": TIER.get(net, "T1"),
                "asked": v["n_asked"],
                "asserted": len(k),
                "false": false,
                "commission": false / len(k) if k else None,
                "omission": (v["n_asked"] - len(k)) / v["n_asked"] if v["n_asked"] else None,
                "declined": v["declined"],
                "not_reached": v["not_reached"],
                "dropped_for_consistency": v["dropped_for_consistency"],
                "meek_consistent_pre": v["meek_consistent_pre"],
                "control_wrong_rate": (
                    ctl["wrong"] / (ctl["right"] + ctl["wrong"])
                    if ctl and (ctl["right"] + ctl["wrong"])
                    else None
                ),
                "control_answered": (ctl.get("right", 0) + ctl.get("wrong", 0)) if ctl else 0,
            }

    def pooled(cond: str, key: str, tier: str | None = None) -> tuple[float, float, float]:
        return boot(
            {
                n: r[key]
                for n, r in per[cond].items()
                if r[key] is not None and (tier is None or r["tier"] == tier)
            }
        )

    def disagreement(a: str, b: str) -> tuple[float, float, float, int]:
        d: dict[str, float] = {}
        shared = 0
        for net, va in know[a]["networks"].items():
            vb = know[b]["networks"].get(net)
            if not vb:
                continue
            oa = {frozenset(e): tuple(e) for e in va["k"]}
            ob = {frozenset(e): tuple(e) for e in vb["k"]}
            both = oa.keys() & ob.keys()
            if not both:
                continue
            shared += len(both)
            d[net] = sum(1 for e in both if oa[e] != ob[e]) / len(both)
        return (*boot(d), shared)

    summary: dict = {"per_condition": {}, "deltas": {}, "per_network": per}
    for cond in know:
        rows = per[cond].values()
        tot_asserted = sum(r["asserted"] for r in rows)
        tot_false = sum(r["false"] for r in rows)
        summary["per_condition"][cond] = {
            "model": know[cond]["model"],
            "networks": len(per[cond]),
            "asked": sum(r["asked"] for r in rows),
            "asserted": tot_asserted,
            "false": tot_false,
            "commission_pooled_ratio": tot_false / tot_asserted if tot_asserted else None,
            "commission_network_mean": pooled(cond, "commission"),
            "omission_network_mean": pooled(cond, "omission"),
            "commission_by_tier": {
                t: pooled(cond, "commission", t) for t in ("T0", "T1", "T2", "T3")
            },
            "meek_inconsistent_networks": sum(1 for r in rows if not r["meek_consistent_pre"]),
            "dropped_for_consistency": sum(r["dropped_for_consistency"] for r in rows),
            "control_wrong_rate_network_mean": pooled(cond, "control_wrong_rate"),
            "control_answered": sum(r["control_answered"] for r in rows),
        }
    for base, other, name in (
        ("D_LLM", "D_LLM_INSTR", "delta_instrument"),
        ("D_LLM", "D_LLM_SRC", "delta_source"),
        ("D_LLM", "D_SCRAMBLED", "delta_scrambled"),
        ("D_LLM", "D_LLM_72B", "delta_scale_7B_72B"),
        ("D_LLM", "D_LLM_32B", "delta_scale_7B_32B"),
        ("D_LLM_72B", "D_LLM_72B_INSTR", "delta_instrument_72B"),
        ("D_LLM_72B", "D_SCRAMBLED_72B", "delta_scrambled_72B"),
        ("D_LLM_72B", "D_LLM_SRC_70B", "delta_source_72B_70B"),
        ("D_LLM_72B", "D_LLM_GEMMA_27B", "delta_source_72B_gemma"),
        ("D_LLM_SRC", "D_LLM_SRC_70B", "delta_scale_8B_70B"),
    ):
        if base in know and other in know:
            m, lo, hi, shared = disagreement(base, other)
            summary["deltas"][name] = {"mean": m, "ci": [lo, hi], "shared_edges": shared}

    # paired contrasts, difference of per-network values on shared networks
    def paired(a: str, b: str, key: str) -> dict | None:
        if a not in per or b not in per:
            return None
        shared = [
            n
            for n in per[a]
            if n in per[b] and per[a][n][key] is not None and per[b][n][key] is not None
        ]
        if not shared:
            return None
        m, lo, hi = boot({n: per[a][n][key] - per[b][n][key] for n in shared})
        return {"mean_diff": m, "ci": [lo, hi], "networks": len(shared)}

    def control_vs_open(cond: str) -> dict | None:
        if cond not in per:
            return None
        shared = [
            n
            for n, r in per[cond].items()
            if r["control_wrong_rate"] is not None and r["commission"] is not None
        ]
        if not shared:
            return None
        m, lo, hi = boot(
            {n: per[cond][n]["control_wrong_rate"] - per[cond][n]["commission"] for n in shared}
        )
        return {"mean_diff": m, "ci": [lo, hi], "networks": len(shared)}

    summary["paired"] = {
        "commission 72B - 7B": paired("D_LLM_72B", "D_LLM", "commission"),
        "commission 32B - 7B": paired("D_LLM_32B", "D_LLM", "commission"),
        "commission 72B - 32B": paired("D_LLM_72B", "D_LLM_32B", "commission"),
        "omission 72B - 7B": paired("D_LLM_72B", "D_LLM", "omission"),
        "commission scrambled - real, 7B": paired("D_SCRAMBLED", "D_LLM", "commission"),
        "commission scrambled - real, 72B": paired("D_SCRAMBLED_72B", "D_LLM_72B", "commission"),
        "commission llama70B - llama8B": paired("D_LLM_SRC_70B", "D_LLM_SRC", "commission"),
        "commission gemma27B - qwen72B": paired("D_LLM_GEMMA_27B", "D_LLM_72B", "commission"),
        "commission gemma27B - llama70B": paired("D_LLM_GEMMA_27B", "D_LLM_SRC_70B", "commission"),
        "control - open, 7B": control_vs_open("D_LLM"),
        "control - open, 72B": control_vs_open("D_LLM_72B"),
    }

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "knowledge_audit.json").write_text(json.dumps(summary, indent=1))

    def f(t: tuple[float, float, float]) -> str:
        return "—" if math.isnan(t[0]) else f"{t[0]:.3f} [{t[1]:.3f}, {t[2]:.3f}]"

    head = (
        f"{'condition':12s} {'asserted':>8} {'false':>5} {'commission (net mean)':>26} "
        f"{'omission':>26} {'control wrong':>26} {'inconsistent':>12}"
    )
    print(head)
    for cond, s in summary["per_condition"].items():
        print(
            f"{cond:12s} {s['asserted']:8d} {s['false']:5d} {f(s['commission_network_mean']):>26} "
            f"{f(s['omission_network_mean']):>26} {f(s['control_wrong_rate_network_mean']):>26} "
            f"{s['meek_inconsistent_networks']:>12}"
        )
    print("\ncommission by tier (network mean):")
    for cond, s in summary["per_condition"].items():
        print(f"  {cond:12s}", {t: f(v) for t, v in s["commission_by_tier"].items()})
    print("\npaired contrasts (difference of per-network values):")
    for k, v in summary["paired"].items():
        if v:
            ci = f"[{v['ci'][0]:+.3f}, {v['ci'][1]:+.3f}]"
            print(f"  {k:36s} {v['mean_diff']:+.3f} {ci}  on {v['networks']} networks")
    print("\nreplicate disagreement on shared edges:")
    for k, v in summary["deltas"].items():
        ci = f"[{v['ci'][0]:.3f}, {v['ci'][1]:.3f}]"
        print(f"  {k:18s} {v['mean']:.3f} {ci}  on {v['shared_edges']} shared edges")


if __name__ == "__main__":
    main()
