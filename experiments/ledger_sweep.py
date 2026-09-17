"""Stage 5 and 7: the certificate on every frozen frame row, then the truth revealed.

For every row of the frozen frame and every knowledge arm, the analyst graph is
built from the CPDAG and that arm's knowledge, the committed set is read off it
in polynomial time, and every instrument answers from observables alone: the
retraction certificate in claims, the hop radius with its dispatch leg, the
fraction of single retractions that break the set, the free hop-count rule, the
closure size and the claim count. The row is then frozen.

Only after that does the generating graph enter, and only to score: how many
asserted claims were false, how many orientable claims were never made, whether
the committed set is valid at the truth under the generalised adjustment
criterion, written on networkx so the audit shares no code with the
instrument, and which of the four audit buckets each instrument's answer falls
in. That ordering is the invariant, and each row
carries a provenance stamp saying where the true DAG touched its inputs.

Arms above the rule assert nothing that came from the truth. The CPDAG they
start from is the true DAG's, which is the idealisation panel, and every such
row is stamped ``via_Chat_only`` rather than ``false``; only an estimated CPDAG
could lift that stamp. Arms below the rule are read off the truth and are
stamped ``via_K``.

Writes ``results/ledger/rows.csv`` and ``results/ledger/sweep_summary.json``.
Networks can be split across processes with ``--networks`` and ``--out`` and
joined back by ``experiments/ledger_merge.py``.

    python experiments/ledger_sweep.py --rand-reps 3
    python experiments/ledger_sweep.py --networks andes,diabetes --out results/ledger/parts/big
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import itertools
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from bkrobust.benchmarks.audit import is_valid_adjustment_set_gac_dag  # noqa: E402
from bkrobust.benchmarks.describe import parse_file, to_mpdag  # noqa: E402
from bkrobust.benchmarks.measure import component_of, separation  # noqa: E402
from bkrobust.core.conventions import NO_RETRACTABLE_EDGES, UNREACHED  # noqa: E402
from bkrobust.demo.evaluate import (  # noqa: E402
    is_valid_adjustment_set_dag,
    optimal_adjustment_set_dag,
)
from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag  # noqa: E402
from bkrobust.hybrid import DEFAULT_SEARCH_BUDGET, breakdown_radius  # noqa: E402
from bkrobust.mpdag_criterion.criterion import is_amenable  # noqa: E402
from bkrobust.mpdag_criterion.optimal import committed_adjustment_set  # noqa: E402
from bkrobust.search.exact import retractable_edges  # noqa: E402

MODELS = ROOT / "results" / "axisa3" / "networks" / "example_models"
FRAME = ROOT / "results" / "frame" / "frame.jsonl"
KNOWLEDGE = ROOT / "results" / "elicit" / "knowledge.json"
OUT = ROOT / "results" / "ledger"

MAX_DEPTH = 3
SUBSET_BUDGET = 300
# The hop radius: exhaustive when the retractable set is small enough for the
# whole space to be walked, otherwise the default depth; the ladder only where
# its quartic model build is affordable, and a closure budget on the walk.
EXHAUSTIVE_MAX_RETRACTABLE = 12
LADDER_MAX_NODES = 60
CLOSURE_BUDGET = 4000
DEPLOYMENT_ARMS = (
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
    # the primary supplier plus the soundly reduced ancestral claims; lives in a second
    # knowledge file (results/elicit/knowledge_ancestral.json) read via --knowledge-extra
    "D_LLM_PLUS_ANC",
)
TIER = {
    **{n: "T0" for n in ("confounding", "mediator", "paths", "M-bias", "M-structure")},
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
RUNGS = ("B2_free_rule", "B4_kg0", "B5_k", "B6_hop", "B7_claim")


def load(name: str) -> tuple[MPDAG, MPDAG, dict[str, str]]:
    """Parse one network and return its DAG, CPDAG and raw-name map."""
    path = next(p for p in sorted(MODELS.iterdir()) if p.name.startswith(name))
    parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
    dag = to_mpdag(parsed)
    return dag, dag_to_cpdag(dag), dict(parsed.name_map)


def sha(g: MPDAG) -> str:
    """A hash of a graph's edge string, for the provenance stamp."""
    return hashlib.sha256(g.edge_string().encode()).hexdigest()[:16]


def breaking_depths(
    cpdag: MPDAG, k: list, x: str, y: str, z: frozenset
) -> tuple[int | None, str, int | None, tuple]:
    """Smallest retraction breaking validity (with its witness) and identifiability."""
    r_claim = r_id = None
    witness: tuple = ()
    status = f"gt_{min(len(k), MAX_DEPTH)}" if k else NO_RETRACTABLE_EDGES
    used = 0
    for depth in range(1, MAX_DEPTH + 1):
        if depth > len(k):
            break
        if used + math.comb(len(k), depth) > SUBSET_BUDGET:
            if r_claim is None:
                status = f"censored_at_depth_{depth}"
            break
        for subset in itertools.combinations(range(len(k)), depth):
            drop = set(subset)
            g = apply_orientations(cpdag, [e for i, e in enumerate(k) if i not in drop])
            used += 1
            if g is None:
                continue
            if r_claim is None and not is_gac_valid_mpdag(g, x, y, z):
                r_claim, status, witness = depth, "exact", tuple(k[i] for i in subset)
            if r_id is None and not is_amenable(g, x, y):
                r_id = depth
            if r_claim is not None and r_id is not None:
                break
        if r_claim is not None and r_id is not None:
            break
    if r_claim is None and status.startswith("gt_") and len(k) <= MAX_DEPTH:
        status = "unreached"  # every non-empty subset was tried
    return r_claim, status, r_id, witness


def claim_lower_bound(r_claim: int | None, status: str) -> int | None:
    """What the claim instrument certifies: the radius, or its censored floor.

    ``None`` means no retraction breaks the set, exhaustively, which the audit
    treats as an unbounded radius. A search censored at depth ``d`` checked
    every subset smaller than ``d`` and so certifies ``r_claim >= d``; one that
    finished its depth without a failure certifies ``r_claim >= depth + 1``.
    """
    if r_claim is not None:
        return r_claim
    if status in ("unreached", NO_RETRACTABLE_EDGES):
        return None
    if status.startswith("censored_at_depth_"):
        return int(status.rsplit("_", 1)[1])
    if status.startswith("gt_"):
        return int(status.rsplit("_", 1)[1]) + 1
    raise ValueError(status)


def bucket(r: int | None, w_c: int, z_valid: bool) -> str:
    """The four-way audit verdict for one instrument's certified radius."""
    covered = True if r is None else w_c <= r - 1
    if covered:
        return "dangerous" if not z_valid else "held"
    return "safe" if not z_valid else "slack"


def chance_level(
    cpdag: MPDAG, base: list, rng: np.random.Generator, retries: int = 20
) -> list | None:
    """Same edges as ``base``, orientations uniform among Meek-consistent choices."""
    for _ in range(retries):
        drawn = [(a, b) if rng.integers(2) else (b, a) for (a, b) in base]
        if apply_orientations(cpdag, drawn) is not None:
            return drawn
    return None


def instrument_row(
    cpdag: MPDAG,
    g0: MPDAG,
    k: list,
    x: str,
    y: str,
    names: dict[str, str],
    radius_limit: float,
) -> tuple[dict, frozenset | None, dict | None]:
    """Everything the analyst can compute on one row, from observables alone.

    Returns the row's fields, the committed set, and the radius each instrument
    certifies keyed by rung. The set and the certificates are ``None`` when the
    row is not certifiable, in which case ``status`` says why. The generating
    DAG is not an input, so nothing here can read it.
    """
    rec: dict = {
        "g0_undirected": len(g0.undirected_edges),
        "k_g0": sum(
            1 for (a, b) in g0.directed_edges if tuple(sorted((a, b))) in cpdag.undirected_edges
        ),
    }
    cs = committed_adjustment_set(g0, x, y)
    rec["committed_verdict"] = cs.verdict
    if cs.z is None:
        rec["status"] = cs.verdict
        return rec, None, None
    z = cs.z
    rec["status"] = "ok"
    rec["z_size"] = len(z)
    sep, sep_status = separation(cpdag, x, z)
    comp = component_of(cpdag, x)
    rec.update(
        {
            "separation": sep,
            "separation_status": sep_status,
            "component_size": len(comp) if comp else 0,
        }
    )
    n_retractable = len(retractable_edges(cpdag, g0))
    if not n_retractable:
        rec.update({"r_hop": UNREACHED, "hop_method": NO_RETRACTABLE_EDGES, "hop_exact": 1})
        hop_lower = None
    else:
        hop = breakdown_radius(
            cpdag,
            None,
            x,
            y,
            z,
            g0=g0,
            search_budget=(
                n_retractable
                if n_retractable <= EXHAUSTIVE_MAX_RETRACTABLE
                else DEFAULT_SEARCH_BUDGET
            ),
            time_limit_s=radius_limit,
            use_ladder=len(cpdag.nodes) <= LADDER_MAX_NODES,
            search_closure_budget=CLOSURE_BUDGET,
        )
        # an inexact UNREACHED is the anytime statement r >= depth + 1
        hop_lower = None if hop.exact or hop.radius != UNREACHED else hop.stats.depth_checked + 1
        rec.update(
            {
                "r_hop": hop.radius,
                "hop_method": hop.method,
                "hop_exact": int(hop.exact),
                "hop_lower_bound": hop_lower,
            }
        )
    r_claim, r_claim_status, r_id, witness = breaking_depths(cpdag, k, x, y, z)
    n_break = sum(
        1
        for i in range(len(k))
        if (gi := apply_orientations(cpdag, [e for j, e in enumerate(k) if j != i])) is not None
        and not is_gac_valid_mpdag(gi, x, y, z)
    )
    rec.update(
        {
            "r_claim": r_claim,
            "r_claim_status": r_claim_status,
            "r_id": r_id,
            "phi_1": round(n_break / len(k), 4) if k else None,
            "witness": " and ".join(f"{names.get(a, a)}->{names.get(b, b)}" for a, b in witness),
        }
    )
    cert = {
        "B2_free_rule": sep if sep_status == "measured" else 1,
        "B4_kg0": rec["k_g0"],
        "B5_k": len(k),
        "B6_hop": (rec["r_hop"] if rec["r_hop"] != UNREACHED else hop_lower),
        "B7_claim": claim_lower_bound(r_claim, r_claim_status),
    }
    return rec, z, cert


def audit_row(dag: MPDAG, cpdag: MPDAG, k: list, x: str, y: str, z: frozenset, cert: dict) -> dict:
    """The truth, revealed after the freeze: commission, omission, validity, buckets."""
    rec: dict = {}
    comp = component_of(cpdag, x)
    w_c = sum(1 for (a, b) in k if not dag.is_directed_edge(a, b))
    comp_edges = {tuple(sorted(e)) for e in cpdag.undirected_edges if comp and e[0] in comp}
    asserted = {tuple(sorted(e)) for e in k}
    rec["n_omitted"] = len(comp_edges - asserted)
    rec["n_wrong_claims"] = w_c
    z_ok = is_valid_adjustment_set_gac_dag(dag, x, y, z)
    rec["z_valid_at_truth"] = int(z_ok)
    rec["z_backdoor_at_truth"] = int(is_valid_adjustment_set_dag(dag, x, y, z))
    # validity first: a set equal to the DAG-level optimal set can still be
    # invalid when the truth orients the query away, since that formula
    # returns the empty set where there is no causal path
    if not z_ok:
        rec["severity"] = "validity"
    elif set(z) == set(optimal_adjustment_set_dag(dag, x, y)):
        rec["severity"] = "benign"
    else:
        rec["severity"] = "efficiency"
    for rung in RUNGS:
        rec[f"bucket_{rung}"] = bucket(cert[rung], w_c, z_ok)
    return rec


def main(argv: list[str] | None = None) -> None:
    """Run every arm on every frozen frame row and write the ledger rows."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--rand-reps", type=int, default=3)
    ap.add_argument("--seed", type=int, default=20260913)
    ap.add_argument("--radius-limit", type=float, default=5.0)
    ap.add_argument("--arms", default="")
    ap.add_argument("--max-nodes", type=int, default=500)
    ap.add_argument("--networks", default="", help="comma list; empty means every network")
    ap.add_argument("--out", default="", help="output directory; default results/ledger")
    ap.add_argument(
        "--knowledge-extra",
        default="",
        help="a second knowledge file; its arm-shaped keys are added, never overriding the first",
    )
    args = ap.parse_args(argv)
    out = Path(args.out) if args.out else OUT
    only = {n for n in args.networks.split(",") if n}

    know = json.loads(KNOWLEDGE.read_text())
    if args.knowledge_extra:
        extra = json.loads(Path(args.knowledge_extra).read_text())
        for a, v in extra.items():
            if a in DEPLOYMENT_ARMS and a not in know and "networks" in v:
                know[a] = v
    present = [a for a in DEPLOYMENT_ARMS if a in know]
    arms = [a for a in args.arms.split(",") if a] or [
        *present,
        "D_RAND",
        "D_DEGEN",
        "A_TRUE",
        "A_ORACLE",
    ]
    frame = [json.loads(line) for line in FRAME.open()]
    nets_with_k = set(know[present[0]]["networks"]) if present else set()
    frame = [r for r in frame if r["network"] in nets_with_k]
    if only:
        frame = [r for r in frame if r["network"] in only]
    by_net: dict[str, list] = collections.defaultdict(list)
    for r in frame:
        by_net[r["network"]].append(r)

    graphs: dict[str, tuple[MPDAG, MPDAG, dict]] = {}
    recovering: dict[str, list] = {}
    excluded: list[dict] = []
    for n in sorted(by_net):
        g = load(n)
        if len(g[0].nodes) > args.max_nodes:
            excluded.append({"network": n, "nodes": len(g[0].nodes)})
            del by_net[n]
            continue
        graphs[n] = g
        recovering[n] = sorted(knowledge_to_recover(g[0], g[1]))
    frame = [r for r in frame if r["network"] in graphs]
    print(f"{len(frame)} frame rows over {len(by_net)} networks; arms {arms}", flush=True)
    print(f"excluded for size: {excluded}", flush=True)

    out.mkdir(parents=True, exist_ok=True)
    fields = [
        "network",
        "tier",
        "x",
        "y",
        "arm",
        "rep",
        "true_dag_on_path",
        "chat_source_sha256",
        "generating_dag_sha256",
        "k_sha256",
        "n_k",
        "status",
        "committed_verdict",
        "z_size",
        "g0_undirected",
        "k_g0",
        "separation",
        "separation_status",
        "component_size",
        "r_claim",
        "r_claim_status",
        "r_id",
        "r_hop",
        "hop_method",
        "hop_exact",
        "hop_lower_bound",
        "phi_1",
        "witness",
        "n_wrong_claims",
        "n_omitted",
        "z_valid_at_truth",
        "z_backdoor_at_truth",
        "severity",
        *[f"bucket_{r}" for r in RUNGS],
    ]
    fh = (out / "rows.csv").open("w", newline="")
    w = csv.DictWriter(fh, fieldnames=fields)
    w.writeheader()

    counts: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    t0 = time.perf_counter()
    done = 0
    for net in sorted(by_net):
        dag, cpdag, names = graphs[net]
        dag_sha, chat_sha = sha(dag), sha(cpdag)
        tier = TIER.get(net, "T1")
        k_true = recovering[net]
        # knowledge per arm for this network
        plans: list[tuple[str, int, list, str]] = []
        base_llm: list = []
        for a in present:
            entry = know[a]["networks"].get(net)
            kk = [tuple(e) for e in (entry or {}).get("k", [])]
            if a == "D_LLM":
                base_llm = kk
            plans.append((a, 0, kk, "via_Chat_only"))
        if "D_RAND" in arms and base_llm:
            for rep in range(args.rand_reps):
                net_seed = int(hashlib.sha256(net.encode()).hexdigest()[:8], 16)
                rng = np.random.default_rng([args.seed, net_seed, rep])
                kk = chance_level(cpdag, base_llm, rng)
                if kk is not None:
                    plans.append(("D_RAND", rep, kk, "via_Chat_only"))
        if "D_DEGEN" in arms:
            plans.append(("D_DEGEN", 0, [], "via_Chat_only"))
        if "A_TRUE" in arms and base_llm:
            truthful = [(a, b) if dag.is_directed_edge(a, b) else (b, a) for (a, b) in base_llm]
            plans.append(("A_TRUE", 0, truthful, "via_K"))
        if "A_ORACLE" in arms:
            plans.append(("A_ORACLE", 0, list(k_true), "via_K"))

        for arm, rep, k, provenance in plans:
            if arm not in arms:
                continue
            k_sha = hashlib.sha256(json.dumps(sorted(k)).encode()).hexdigest()[:16]
            g0 = apply_orientations(cpdag, k)
            for row in by_net[net]:
                x, y = row["x"], row["y"]
                done += 1
                rec: dict = {
                    "network": net,
                    "tier": tier,
                    "x": x,
                    "y": y,
                    "arm": arm,
                    "rep": rep,
                    "true_dag_on_path": provenance,
                    "chat_source_sha256": chat_sha,
                    "generating_dag_sha256": dag_sha,
                    "k_sha256": k_sha,
                    "n_k": len(k),
                }
                if g0 is None:
                    rec["status"] = "knowledge_inconsistent"
                    counts[arm][rec["status"]] += 1
                    w.writerow(rec)
                    continue
                # observables only, from here to the freeze
                irec, z, cert = instrument_row(cpdag, g0, k, x, y, names, args.radius_limit)
                rec.update(irec)
                if z is None or cert is None:
                    counts[arm][rec["status"]] += 1
                    w.writerow(rec)
                    continue
                # ---- the truth, revealed only now
                rec.update(audit_row(dag, cpdag, k, x, y, z, cert))
                counts[arm]["ok"] += 1
                w.writerow(rec)
        print(f"  {net}: {done} evaluations, {time.perf_counter() - t0:.0f} s", flush=True)
    fh.close()

    summary = {
        "frame_rows": len(frame),
        "networks": len(by_net),
        "arms": arms,
        "evaluations": done,
        "seconds": round(time.perf_counter() - t0, 1),
        "status_by_arm": {a: dict(c) for a, c in sorted(counts.items())},
        "settings": {
            "rand_reps": args.rand_reps,
            "seed": args.seed,
            "max_depth": MAX_DEPTH,
            "subset_budget": SUBSET_BUDGET,
            "radius_limit_s": args.radius_limit,
            "exhaustive_max_retractable": EXHAUSTIVE_MAX_RETRACTABLE,
            "ladder_max_nodes": LADDER_MAX_NODES,
            "closure_budget": CLOSURE_BUDGET,
        },
        "excluded_for_size": excluded,
    }
    if args.knowledge_extra:
        # outside "settings" so ledger_merge.py still joins parts run without it
        summary["knowledge_extra"] = str(Path(args.knowledge_extra))
    (out / "sweep_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary["status_by_arm"], indent=1))


if __name__ == "__main__":
    main()
