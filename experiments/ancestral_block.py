"""The ancestral block: the one knowledge granularity the causal-order questionnaire leaves out.

The frozen questionnaire asks a causal order per chain component and reads it
back per adjacent pair. That already yields tiers (the order is one), required
edges (an orientation on a skeleton edge is one) and forbidden edges (on a
skeleton edge, forbidding one direction is requiring the other). What it cannot
express is an *ancestral* claim between two variables that are not adjacent:
"A causes B, possibly through other variables". This block asks that question
and converts the answers into orientation claims without guessing.

**Pairs.** For each chain component of the frozen questionnaire, every pair of
component nodes that is not adjacent and sits at distance two in the component's
undirected skeleton, so that at least one possibly directed path joins them in
each direction and the question is never empty. At most eight pairs per network,
drawn by a seeded permutation. Rendered like the open items: real names, the
same domain sentence and variable states, node and pair order permuted under
order seed 0. The frozen questionnaire is not modified; this block is a separate
file with its own hash.

**Grammar.** Per pair: ``a causes b``, ``b causes a``, ``neither``, ``DECLINE``,
with a confidence that is retained and quarantined. ``neither`` is recorded and
never reduced; a pair the model does not mention is ``not_reached``.

**The reduction rule** (:mod:`bkrobust.knowledge.ancestral`, tested in
``tests/knowledge/``). An ancestral claim ``a`` causes ``b`` becomes orientation
claims only when it entails them on the CPDAG. Let ``P`` be the possibly
directed paths from ``a`` to ``b``. If every path in ``P`` shares its first
edge and that edge is undirected, orient it away from ``a``; if every path
shares its last edge and that edge is undirected, orient it into ``b``. Any
other claim is ``not_reducible`` and is counted, never converted: the paths
diverge, the unique end edges are already compelled, or no path exists and the
claim contradicts the data. The rule is sound and incomplete by design.

**Freeze, then audit.** The reduced orientation claims are written and hashed
before the true DAG is read. Only then is each reduced claim scored at the
truth (commission) and each ancestral claim checked for ancestry in the true
DAG. The merged deployment arm ``D_LLM_PLUS_ANC`` is the primary supplier's
frozen orientation set plus the reduced claims: a reduced claim on an edge the
order already oriented is dropped as ``already_asserted`` when it agrees and as
``conflicts_with_orientation`` when it does not (the direct question wins); the
union is Meek-checked and, if inconsistent, repaired by dropping reduced claims
in increasing confidence, and the count is printed. ``ledger_sweep.py`` reads
the arm through ``--knowledge-extra``.

    python experiments/ancestral_block.py build              # questionnaire_ancestral.json, .sha256
    python experiments/ancestral_block.py run                # knowledge_ancestral.json
    python experiments/ancestral_block.py run --export-only  # ancestral_prompts.jsonl, no call
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT / "src", ROOT):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from bkrobust.benchmarks.describe import parse_file, to_mpdag  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.knowledge.ancestral import REDUCED, reduce_ancestral_claim  # noqa: E402
from experiments.elicit_run import CONDITIONS, FAMILY, HOST, ask  # noqa: E402
from experiments.questionnaire_build import DOMAIN, MODELS, states_from_bif  # noqa: E402

ELICIT = ROOT / "results" / "elicit"
QUESTIONNAIRE = ELICIT / "questionnaire_ancestral.json"
KNOWLEDGE = ELICIT / "knowledge_ancestral.json"
PROMPTS = ELICIT / "ancestral_prompts.jsonl"

PAIRS_PER_NETWORK = 8
ORDER_SEED = 0
BASE_CONDITION = "D_LLM"
ARM = "D_LLM_PLUS_ANC"

_FENCE = re.compile(r"^```(?:json)?\s*|\s*```$", re.M)


def net_int(name: str) -> int:
    """A stable per-network integer for seeding; the built-in hash is salted per process."""
    return int(hashlib.sha256(name.encode()).hexdigest()[:8], 16)


def graphs() -> dict[str, tuple[MPDAG, MPDAG, dict[str, str], dict[str, list[str]]]]:
    """Every network's DAG, CPDAG, raw-name map and variable states."""
    out = {}
    for path in sorted(MODELS.iterdir()):
        parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
        if parsed.bidirected:
            continue
        dag = to_mpdag(parsed)
        states = states_from_bif(path) if parsed.source_format == "bif" else {}
        out[parsed.name] = (dag, dag_to_cpdag(dag), dict(parsed.name_map), states)
    return out


def render(
    domain: str,
    names: dict[str, str],
    nodes: list[str],
    pairs: list[tuple[str, str]],
    states: dict[str, list[str]],
) -> str:
    """The prompt. Order of ``nodes`` and ``pairs`` is the caller's, already permuted."""
    lines = [
        "You are asked about causal influence among variables from " + domain + ".",
        "Some of these variables are directly related, but the data alone cannot tell which "
        "of each related pair causes the other.",
        "",
        "Variables in this group:",
    ]
    for n in nodes:
        st = states.get(n)
        lines.append(f"  {names[n]}" + (f" (values: {', '.join(st)})" if st else ""))
    lines += [
        "",
        "Each question below is about two variables that are NOT directly related; one may "
        "still be a cause of the other through other variables in the group.",
    ]
    for i, (a, b) in enumerate(pairs, 1):
        lines.append(
            f"  {i}. Is {names[a]} a cause of {names[b]}, possibly through other variables?"
        )
    lines += [
        "",
        'For each question answer "a causes b" if the first variable named is a cause of the '
        'second, possibly through other variables; "b causes a" if instead the second is a cause '
        'of the first, possibly through other variables; "neither" if neither is a cause of the '
        "other, even indirectly; or DECLINE if you cannot judge it. Give a confidence from 0 to 1.",
        "",
        'Answer ONLY with JSON: {"pairs": [{"a": "<first name as asked>", "b": "<second name as '
        'asked>", "answer": "a causes b" | "b causes a" | "neither" | "DECLINE", '
        '"confidence": 0.0}, ...]}',
    ]
    return "\n".join(lines)


def candidate_pairs(cpdag: MPDAG, comp: list[str]) -> list[tuple[str, str]]:
    """Non-adjacent pairs of the component at distance two in its undirected skeleton."""
    inside = set(comp)
    out = []
    for i, a in enumerate(comp):
        na = cpdag.neighbors(a) & inside
        for b in comp[i + 1 :]:
            if cpdag.has_edge(a, b):
                continue
            if na & cpdag.neighbors(b):
                out.append((a, b))
    return sorted(out)


def build(seed: int) -> None:
    """Derive the ancestral questionnaire from the frozen one and hash it."""
    q = json.loads((ELICIT / "questionnaire.json").read_text())
    qhash = (ELICIT / "questionnaire.sha256").read_text().strip()
    comps: dict[str, list[tuple[int, list[str]]]] = {}
    for it in q["items"]:
        if it["kind"] == "open" and it["naming"] == "real" and it["order_seed"] == ORDER_SEED:
            comps.setdefault(it["network"], []).append((it["component"], it["nodes"]))

    gs = graphs()
    items: list[dict] = []
    per_network: list[dict] = []
    for net in sorted(comps):
        _dag, cpdag, raw, states_raw = gs[net]
        real = {n: raw.get(n, n) for n in cpdag.nodes}
        states = {n: states_raw[raw.get(n, n)] for n in cpdag.nodes if raw.get(n, n) in states_raw}
        cands = [
            (ci, a, b) for ci, nodes in sorted(comps[net]) for a, b in candidate_pairs(cpdag, nodes)
        ]
        n_cand = len(cands)
        if len(cands) > PAIRS_PER_NETWORK:
            rng = np.random.default_rng([seed, 13, net_int(net)])
            idx = sorted(rng.choice(len(cands), PAIRS_PER_NETWORK, replace=False))
            cands = [cands[int(i)] for i in idx]
        by_comp: dict[int, list[tuple[str, str]]] = {}
        for ci, a, b in cands:
            by_comp.setdefault(ci, []).append((a, b))
        nodes_of = dict(comps[net])
        for ci in sorted(by_comp):
            nodes = sorted(nodes_of[ci])
            pairs = sorted(by_comp[ci])
            rng = np.random.default_rng([seed, ORDER_SEED, ci, net_int(net)])
            node_order = [nodes[int(i)] for i in rng.permutation(len(nodes))]
            pair_order = [pairs[int(i)] for i in rng.permutation(len(pairs))]
            prompt = render(
                DOMAIN.get(net, "an applied study"), real, node_order, pair_order, states
            )
            items.append(
                {
                    "network": net,
                    "component": ci,
                    "kind": "ancestral",
                    "naming": "real",
                    "order_seed": ORDER_SEED,
                    "nodes": nodes,
                    "pairs": pairs,
                    "label_of": {n: real[n] for n in nodes},
                    "prompt": prompt,
                    "prompt_sha256": hashlib.sha256(prompt.encode()).hexdigest(),
                }
            )
        per_network.append(
            {
                "network": net,
                "components": len(comps[net]),
                "candidate_pairs": n_cand,
                "asked_pairs": len(cands),
                "items": len(by_comp),
            }
        )

    payload = {
        "seed": seed,
        "order_seed": ORDER_SEED,
        "pairs_per_network": PAIRS_PER_NETWORK,
        "pair_rule": "non-adjacent, distance two in the component's undirected skeleton",
        "questionnaire_sha256": qhash,
        "cpdag_source": "dag_to_cpdag(true DAG); true_dag_on_path = via_Chat_only",
        "per_network": per_network,
        "items": items,
    }
    text = json.dumps(payload, indent=1, sort_keys=True)
    QUESTIONNAIRE.write_text(text)
    digest = hashlib.sha256(text.encode()).hexdigest()
    QUESTIONNAIRE.with_suffix(".sha256").write_text(digest + "\n")
    print(
        f"{len(items)} ancestral items over {len(per_network)} networks, "
        f"{sum(p['asked_pairs'] for p in per_network)} pairs asked of "
        f"{sum(p['candidate_pairs'] for p in per_network)} candidates"
    )
    print("questionnaire_ancestral sha256", digest)


def parse(response: str, labels: dict[str, str], pairs: list[list[str]]) -> dict:
    """Read the model's JSON into ancestral claims over the asked pairs."""
    inv = {v: k for k, v in labels.items()}
    text = _FENCE.sub("", response.strip())
    start, end = text.find("{"), text.rfind("}")
    out: dict = {
        "claims": {},
        "confidence": {},
        "neither": [],
        "declined": [],
        "not_reached": [],
        "invented": 0,
        "unparsed_answer": 0,
        "parse_ok": False,
    }
    data = None
    if start >= 0 and end >= 0:
        try:
            data = json.loads(text[start : end + 1])
        except json.JSONDecodeError:
            data = None
    if not isinstance(data, dict):
        out["not_reached"] = [f"{a}|{b}" for a, b in pairs]
        return out
    out["parse_ok"] = True
    asked = {frozenset(p) for p in pairs}
    answered: dict[frozenset, tuple[str, tuple[str, str] | None, float | None]] = {}
    for p in data.get("pairs") or []:
        if not isinstance(p, dict):
            out["invented"] += 1
            continue
        a, b = inv.get(str(p.get("a", "")).strip()), inv.get(str(p.get("b", "")).strip())
        if not a or not b or a == b or frozenset((a, b)) not in asked:
            out["invented"] += 1
            continue
        ans = re.sub(r"\s+", " ", str(p.get("answer", "")).strip().lower())
        conf = p.get("confidence")
        try:
            conf = float(conf) if conf is not None else None
        except (TypeError, ValueError):
            conf = None
        la, lb = labels[a].lower(), labels[b].lower()
        if ans in ("a causes b", f"{la} causes {lb}", "a->b", "a -> b"):
            answered[frozenset((a, b))] = ("asserted", (a, b), conf)
        elif ans in ("b causes a", f"{lb} causes {la}", "b->a", "b -> a"):
            answered[frozenset((a, b))] = ("asserted", (b, a), conf)
        elif ans in ("neither", "none", "no", "neither causes the other"):
            answered[frozenset((a, b))] = ("neither", None, conf)
        elif ans == "decline":
            answered[frozenset((a, b))] = ("declined", None, conf)
        else:
            out["unparsed_answer"] += 1
            answered[frozenset((a, b))] = ("declined", None, conf)
    for a, b in pairs:
        key = f"{a}|{b}"
        kind, claim, conf = answered.get(frozenset((a, b)), (None, None, None))
        if kind is None:
            out["not_reached"].append(key)
        elif kind == "asserted" and claim is not None:
            out["claims"][key] = list(claim)
            out["confidence"][key] = conf
        else:
            out[kind].append(key)
    return out


def host_reachable() -> bool:
    """Whether the local GPU host answers on its model API, within a short wall."""
    try:
        proc = subprocess.run(
            [
                "ssh",
                "-o",
                "BatchMode=yes",
                "-o",
                "ConnectTimeout=10",
                HOST,
                "curl",
                "-s",
                "-m",
                "8",
                "http://localhost:11434/api/tags",
            ],
            capture_output=True,
            text=True,
            timeout=30,
        )
    except (subprocess.SubprocessError, OSError):
        return False
    return proc.returncode == 0 and "models" in proc.stdout


def export(items: list[dict]) -> None:
    """Write the prompts in the format of ``elicit_export.py``, for a supplier elsewhere."""
    with PROMPTS.open("w") as fh:
        for it in items:
            fh.write(
                json.dumps(
                    {
                        "prompt_sha256": it["prompt_sha256"],
                        "network": it["network"],
                        "kind": it["kind"],
                        "naming": it["naming"],
                        "order_seed": it["order_seed"],
                        "prompt": it["prompt"],
                    }
                )
                + "\n"
            )
    print(f"{len(items)} prompts -> {PROMPTS}")


def merge_arm(
    cpdag: MPDAG, base: list[list[str]], reduced: dict[tuple[str, str], float | None]
) -> dict:
    """The primary supplier's frozen orientations plus the reduced claims, made consistent."""
    base_t = [tuple(e) for e in base]
    have = set(base_t)
    added: dict[tuple[str, str], float | None] = {}
    already = conflicts = 0
    for (t, h), conf in sorted(reduced.items(), key=lambda kv: kv[0]):
        if (t, h) in have:
            already += 1
        elif (h, t) in have:
            conflicts += 1
        else:
            added[(t, h)] = conf
    consistent = apply_orientations(cpdag, [*base_t, *added]) is not None
    dropped = 0
    if not consistent:
        order = sorted(added, key=lambda e: (added[e] is None, added[e] or 0.0, e))
        while added and apply_orientations(cpdag, [*base_t, *added]) is None:
            added.pop(order[dropped])
            dropped += 1
    k = sorted({*base_t, *added})
    return {
        "k": [list(e) for e in k],
        "k_sha256": hashlib.sha256(json.dumps([list(e) for e in k]).encode()).hexdigest(),
        "n_asserted": len(k),
        "n_from_order": len(base_t),
        "n_from_ancestral": len(added),
        "already_asserted": already,
        "conflicts_with_orientation": conflicts,
        "meek_consistent_pre": consistent,
        "dropped_for_consistency": dropped,
    }


def run(export_only: bool) -> None:
    """Ask the primary supplier, reduce soundly, freeze, then audit."""
    q = json.loads(QUESTIONNAIRE.read_text())
    qhash = QUESTIONNAIRE.with_suffix(".sha256").read_text().strip()
    assert hashlib.sha256(QUESTIONNAIRE.read_text().encode()).hexdigest() == qhash, "hash drift"
    spec = CONDITIONS[BASE_CONDITION]
    model = spec["model"]
    items = q["items"]
    if export_only or not host_reachable():
        export(items)
        if not export_only:
            print(f"host {HOST} unreachable: prompts exported, nothing asked, no knowledge written")
        return

    gs = graphs()
    base_know = json.loads((ELICIT / "knowledge.json").read_text())[BASE_CONDITION]
    answers = []
    per_net: dict[str, dict] = {}
    hits = 0
    t0 = time.perf_counter()
    for i, it in enumerate(items):
        rec = ask(model, it["prompt"])
        hits += rec["cache_hit"]
        parsed = parse(rec["response"], it["label_of"], it["pairs"])
        answers.append(
            {
                **{k: it[k] for k in ("network", "component", "kind", "naming", "order_seed")},
                "prompt_sha256": it["prompt_sha256"],
                "model": model,
                "cache_key": rec["key"],
                **parsed,
            }
        )
        pn = per_net.setdefault(
            it["network"],
            {
                "asked": 0,
                "claims": {},
                "confidence": {},
                "neither": 0,
                "declined": 0,
                "not_reached": 0,
                "invented": 0,
                "unparsed_answer": 0,
                "parse_fail": 0,
            },
        )
        pn["asked"] += len(it["pairs"])
        pn["claims"].update(parsed["claims"])
        pn["confidence"].update(parsed["confidence"])
        for key in ("neither", "declined", "not_reached"):
            pn[key] += len(parsed[key])
        pn["invented"] += parsed["invented"]
        pn["unparsed_answer"] += parsed["unparsed_answer"]
        pn["parse_fail"] += 0 if parsed["parse_ok"] else 1
        if i % 10 == 0:
            print(f"  {i}/{len(items)} items, {hits} cache hits, {time.perf_counter() - t0:.0f} s")

    # reduce, then freeze: no truth has been read yet
    networks: dict[str, dict] = {}
    arm_networks: dict[str, dict] = {}
    for net, pn in sorted(per_net.items()):
        cpdag = gs[net][1]
        claims = []
        reduced: dict[tuple[str, str], float | None] = {}
        for key, (tail, head) in sorted(pn["claims"].items()):
            red = reduce_ancestral_claim(cpdag, tail, head)
            conf = pn["confidence"].get(key)
            claims.append(
                {
                    "pair": key,
                    "tail": tail,
                    "head": head,
                    "confidence": conf,
                    "status": red.status,
                    "reason": red.reason,
                    "n_first_edges": red.n_first_edges,
                    "n_last_edges": red.n_last_edges,
                    "orientations": [list(e) for e in red.orientations],
                }
            )
            for e in red.orientations:
                prev = reduced.get(e)
                reduced[e] = conf if prev is None else (prev if conf is None else min(prev, conf))
        k_anc = sorted(reduced)
        networks[net] = {
            "asked": pn["asked"],
            "asserted": len(claims),
            "reduced": sum(1 for c in claims if c["status"] == REDUCED),
            "not_reducible": sum(1 for c in claims if c["status"] != REDUCED),
            "neither": pn["neither"],
            "declined": pn["declined"],
            "not_reached": pn["not_reached"],
            "invented": pn["invented"],
            "unparsed_answer": pn["unparsed_answer"],
            "parse_fail": pn["parse_fail"],
            "not_reducible_reasons": sorted(
                {c["reason"] for c in claims if c["status"] != REDUCED}
            ),
            "ancestral_claims": claims,
            "k_anc": [list(e) for e in k_anc],
            "k_anc_confidence": {f"{t}|{h}": reduced[(t, h)] for t, h in k_anc},
            "k_anc_sha256": hashlib.sha256(
                json.dumps([list(e) for e in k_anc]).encode()
            ).hexdigest(),
            "confidence_quarantined": {k: v for k, v in pn["confidence"].items()},
        }
    for net, entry in sorted(base_know["networks"].items()):
        cpdag = gs[net][1]
        anc = networks.get(net, {"k_anc": [], "k_anc_confidence": {}})
        reduced = {tuple(e): anc["k_anc_confidence"][f"{e[0]}|{e[1]}"] for e in anc["k_anc"]}
        arm_networks[net] = merge_arm(cpdag, entry["k"], reduced)

    def totals(nets: dict[str, dict], keys: tuple[str, ...]) -> dict[str, int]:
        return {k: sum(v[k] for v in nets.values()) for k in keys}

    payload = {
        "ANCESTRAL": {
            "model": model,
            "naming": "real",
            "order_seed": ORDER_SEED,
            "family": FAMILY.get(model, "?"),
            "questionnaire_sha256": q["questionnaire_sha256"],
            "questionnaire_ancestral_sha256": qhash,
            "pairs_per_network": PAIRS_PER_NETWORK,
            "reduction_rule": (
                "a causes b -> orient the unique undirected first edge of every possibly directed "
                "a..b path away from a, and the unique undirected last edge into b; else "
                "not_reducible"
            ),
            "cache_hits": hits,
            "networks": networks,
            "totals": totals(
                networks,
                (
                    "asked",
                    "asserted",
                    "reduced",
                    "not_reducible",
                    "neither",
                    "declined",
                    "not_reached",
                    "invented",
                    "unparsed_answer",
                    "parse_fail",
                ),
            )
            | {"reduced_orientations": sum(len(v["k_anc"]) for v in networks.values())},
        },
        ARM: {
            **spec,
            "family": FAMILY.get(model, "?"),
            "questionnaire_sha256": q["questionnaire_sha256"],
            "questionnaire_ancestral_sha256": qhash,
            "base_condition": BASE_CONDITION,
            "networks": arm_networks,
            "true_dag_on_path": "false",
            "arm": "deployment",
            "totals": totals(
                arm_networks,
                (
                    "n_asserted",
                    "n_from_order",
                    "n_from_ancestral",
                    "already_asserted",
                    "conflicts_with_orientation",
                    "dropped_for_consistency",
                ),
            ),
        },
    }
    KNOWLEDGE.write_text(json.dumps(payload, indent=1, sort_keys=True))
    (ELICIT / f"answers_{ARM}.json").write_text(
        json.dumps(
            {
                "condition": "ANCESTRAL",
                "model": model,
                "questionnaire_ancestral_sha256": qhash,
                "cache_hits": hits,
                "items": answers,
            },
            indent=1,
        )
    )
    frozen_hashes = {n: v["k_anc_sha256"] for n, v in networks.items()}
    print("frozen: " + json.dumps(payload["ANCESTRAL"]["totals"]))
    print("merged: " + json.dumps(payload[ARM]["totals"]))

    # ---- the truth, read only now, and only to score
    for net, v in networks.items():
        dag = gs[net][0]
        assert v["k_anc_sha256"] == frozen_hashes[net]
        k_anc = [tuple(e) for e in v["k_anc"]]
        false_o = sum(1 for t, h in k_anc if not dag.is_directed_edge(t, h))
        false_a = sum(1 for c in v["ancestral_claims"] if c["tail"] not in dag.ancestors(c["head"]))
        false_a_red = sum(
            1
            for c in v["ancestral_claims"]
            if c["status"] == REDUCED and c["tail"] not in dag.ancestors(c["head"])
        )
        v["audit"] = {
            "true_dag_on_path": "via_K",
            "reduced_orientations": len(k_anc),
            "reduced_false": false_o,
            "reduced_commission": false_o / len(k_anc) if k_anc else None,
            "ancestral_asserted": len(v["ancestral_claims"]),
            "ancestral_false": false_a,
            "ancestral_commission": (
                false_a / len(v["ancestral_claims"]) if v["ancestral_claims"] else None
            ),
            "ancestral_false_among_reduced": false_a_red,
        }
    for net, v in arm_networks.items():
        dag = gs[net][0]
        k = [tuple(e) for e in v["k"]]
        v["audit"] = {
            "true_dag_on_path": "via_K",
            "false": sum(1 for t, h in k if not dag.is_directed_edge(t, h)),
            "commission": (
                sum(1 for t, h in k if not dag.is_directed_edge(t, h)) / len(k) if k else None
            ),
        }
    n_red = sum(v["audit"]["reduced_orientations"] for v in networks.values())
    n_false = sum(v["audit"]["reduced_false"] for v in networks.values())
    n_anc = sum(v["audit"]["ancestral_asserted"] for v in networks.values())
    n_anc_false = sum(v["audit"]["ancestral_false"] for v in networks.values())
    per_net_comm = [
        v["audit"]["reduced_commission"]
        for v in networks.values()
        if v["audit"]["reduced_commission"] is not None
    ]
    payload["ANCESTRAL"]["audit"] = {
        "true_dag_on_path": "via_K",
        "reduced_orientations": n_red,
        "reduced_false": n_false,
        "reduced_commission_pooled": n_false / n_red if n_red else None,
        "reduced_commission_network_mean": (float(np.mean(per_net_comm)) if per_net_comm else None),
        "networks_with_reduced": len(per_net_comm),
        "ancestral_asserted": n_anc,
        "ancestral_false": n_anc_false,
        "ancestral_commission_pooled": n_anc_false / n_anc if n_anc else None,
    }
    KNOWLEDGE.write_text(json.dumps(payload, indent=1, sort_keys=True))
    print("audit: " + json.dumps(payload["ANCESTRAL"]["audit"]))
    print(f"{KNOWLEDGE} written, {time.perf_counter() - t0:.0f} s")


def main(argv: list[str] | None = None) -> None:
    """Build the block, or run it."""
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=("build", "run"))
    ap.add_argument("--seed", type=int, default=20260913)
    ap.add_argument("--export-only", action="store_true")
    args = ap.parse_args(argv)
    if args.stage == "build":
        build(args.seed)
    else:
        run(args.export_only)


if __name__ == "__main__":
    main()
