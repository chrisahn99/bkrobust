"""Stage 4: the suppliers answer the frozen questionnaire, and the answers become knowledge.

Every call is temperature zero, cached against the SHA-256 of ``model | prompt``,
and never repeated once cached, so the results reproduce without re-querying a
model that may have changed underneath them. The model runs locally on the GPU
host through its HTTP API; nothing leaves the two machines.

Conditions, each a column of the ledger:

``D_LLM``       the primary supplier, real names, order seed 0
``D_LLM_INSTR`` the same model, same names, order seed 1 -- instrument variance
``D_LLM_SRC``   a second model family, real names, order seed 0 -- source variance
``D_SCRAMBLED`` the primary model on scrambled names -- the recitation probe

The order the model returns is the primary instrument; the per-pair directions
are read against it. A pair whose stated direction contradicts the model's own
order is counted and the order wins, because an order cannot contain a cycle and
a list of directions can. ``DECLINE`` is an answer and is kept as one; a pair the
model never mentioned is ``not_reached``, which is a different thing.

Per network the asserted orientations are Meek-closed against the CPDAG. A set
admitting no consistent MPDAG is recorded as such, that is the Meek rejection
rate the audit reports, and is then repaired the way a practitioner repairs it:
drop claims in increasing order of confidence until the closure succeeds, and
record how many were dropped. Both the pre-filter and post-filter sets are kept,
and the post-filter set is frozen and hashed before any radius is computed.

The compelled-edge control block is scored here, against what the data compel,
and never enters the knowledge set.

Writes ``results/elicit/answers_<condition>.json`` and ``results/elicit/knowledge.json``.

    python experiments/elicit_run.py --conditions D_LLM
    python experiments/elicit_run.py --conditions D_LLM,D_LLM_INSTR,D_LLM_SRC,D_SCRAMBLED
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import json
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from bkrobust.benchmarks.describe import parse_file, to_mpdag  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402

MODELS = ROOT / "results" / "axisa3" / "networks" / "example_models"
ELICIT = ROOT / "results" / "elicit"
CACHE = ELICIT / "cache"
HOST = "betelgeuse"

CONDITIONS: dict[str, dict] = {
    "D_LLM": {"model": "qwen2.5:7b-instruct", "naming": "real", "order_seed": 0},
    "D_LLM_INSTR": {"model": "qwen2.5:7b-instruct", "naming": "real", "order_seed": 1},
    "D_LLM_SRC": {"model": "llama3.1:8b", "naming": "real", "order_seed": 0},
    "D_SCRAMBLED": {"model": "qwen2.5:7b-instruct", "naming": "scrambled", "order_seed": 0},
    # cluster suppliers, answered through vLLM on the shared checkpoints and imported
    # into the same cache by ``elicit_import.py``; never called from here
    "D_LLM_72B": {"model": "Qwen2.5-72B-Instruct", "naming": "real", "order_seed": 0},
    "D_LLM_72B_INSTR": {"model": "Qwen2.5-72B-Instruct", "naming": "real", "order_seed": 1},
    "D_SCRAMBLED_72B": {"model": "Qwen2.5-72B-Instruct", "naming": "scrambled", "order_seed": 0},
    "D_LLM_32B": {"model": "Qwen2.5-32B-Instruct", "naming": "real", "order_seed": 0},
    "D_LLM_SRC_70B": {"model": "Llama-3.3-70B-Instruct", "naming": "real", "order_seed": 0},
    "D_LLM_GEMMA_27B": {"model": "gemma-3-27b-it", "naming": "real", "order_seed": 0},
}
FAMILY = {
    "qwen2.5:7b-instruct": "qwen2",
    "llama3.1:8b": "llama",
    "Qwen2.5-72B-Instruct": "qwen2",
    "Qwen2.5-32B-Instruct": "qwen2",
    "Llama-3.3-70B-Instruct": "llama",
    "gemma-3-27b-it": "gemma",
}
CLUSTER_MODELS = {
    "Qwen2.5-72B-Instruct",
    "Qwen2.5-32B-Instruct",
    "Llama-3.3-70B-Instruct",
    "gemma-3-27b-it",
}

_FENCE = re.compile(r"^```(?:json)?\s*|\s*```$", re.M)


def ask(model: str, prompt: str, num_predict: int = 700, timeout: int = 240) -> dict:
    """One temperature-zero completion, from the cache if it is there."""
    key = hashlib.sha256(f"{model}|{prompt}".encode()).hexdigest()
    path = CACHE / f"{key}.json"
    if path.exists():
        rec = json.loads(path.read_text())
        rec["cache_hit"] = True
        return rec
    if model in CLUSTER_MODELS:
        raise SystemExit(f"{model}: answer not in the cache; import it with elicit_import.py")
    payload = json.dumps(
        {
            "model": model,
            "prompt": prompt,
            "stream": False,
            "think": False,
            "options": {"temperature": 0.0, "seed": 0, "num_predict": num_predict},
        }
    )
    t = time.perf_counter()
    proc = subprocess.run(
        [
            "ssh",
            "-o",
            "BatchMode=yes",
            HOST,
            "curl",
            "-s",
            "-m",
            str(timeout),
            "http://localhost:11434/api/generate",
            "-d",
            "@-",
        ],
        input=payload,
        capture_output=True,
        text=True,
        timeout=timeout + 30,
    )
    raw = json.loads(proc.stdout or "{}")
    rec = {
        "key": key,
        "model": model,
        "prompt_sha256": hashlib.sha256(prompt.encode()).hexdigest(),
        "response": raw.get("response", ""),
        "eval_count": raw.get("eval_count"),
        "eval_ms": (raw.get("eval_duration") or 0) // 1_000_000,
        "wall_s": round(time.perf_counter() - t, 2),
        "cache_hit": False,
    }
    CACHE.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rec, indent=1))
    return rec


def parse(response: str, labels: dict[str, str], pairs: list[list[str]]) -> dict:
    """Read the model's JSON into orientations over the asked pairs.

    Returns the orientations, the confidences, and the bookkeeping the audit
    needs: pairs declined, pairs never reached, pairs invented, and pairs whose
    stated direction contradicts the model's own order.
    """
    inv = {v: k for k, v in labels.items()}  # label -> node
    text = _FENCE.sub("", response.strip())
    start, end = text.find("{"), text.rfind("}")
    out = {
        "orientations": {},
        "confidence": {},
        "declined": [],
        "not_reached": [],
        "invented": 0,
        "contradicts_order": 0,
        "parse_ok": False,
        "order_len": 0,
    }
    if start < 0 or end < 0:
        out["not_reached"] = [f"{a}|{b}" for a, b in pairs]
        return out
    try:
        data = json.loads(text[start : end + 1])
    except json.JSONDecodeError:
        out["not_reached"] = [f"{a}|{b}" for a, b in pairs]
        return out
    out["parse_ok"] = True
    order = [inv.get(str(v).strip()) for v in data.get("order") or []]
    order = [v for v in order if v]
    rank = {v: i for i, v in enumerate(order)}
    out["order_len"] = len(order)

    asked = {frozenset(p) for p in pairs}
    answered: dict[frozenset, tuple[str, float | None]] = {}
    for p in data.get("pairs") or []:
        a, b = inv.get(str(p.get("a", "")).strip()), inv.get(str(p.get("b", "")).strip())
        if not a or not b or frozenset((a, b)) not in asked:
            out["invented"] += 1
            continue
        d = str(p.get("direction", "")).strip()
        conf = p.get("confidence")
        try:
            conf = float(conf) if conf is not None else None
        except (TypeError, ValueError):
            conf = None
        la, lb = labels[a], labels[b]
        if d.upper() == "DECLINE":
            answered[frozenset((a, b))] = ("DECLINE", conf)
        elif d in ("a->b", f"{la}->{lb}"):
            answered[frozenset((a, b))] = (f"{a}->{b}", conf)
        elif d in ("b->a", f"{lb}->{la}"):
            answered[frozenset((a, b))] = (f"{b}->{a}", conf)
        else:
            answered[frozenset((a, b))] = ("DECLINE", conf)

    for a, b in pairs:
        k = frozenset((a, b))
        stated, conf = answered.get(k, (None, None))
        from_order = None
        if a in rank and b in rank and rank[a] != rank[b]:
            from_order = f"{a}->{b}" if rank[a] < rank[b] else f"{b}->{a}"
        if stated == "DECLINE":
            out["declined"].append(f"{a}|{b}")
            continue
        if stated is None and from_order is None:
            out["not_reached"].append(f"{a}|{b}")
            continue
        chosen = from_order or stated
        if stated and from_order and stated != from_order:
            out["contradicts_order"] += 1
        tail, head = chosen.split("->")
        out["orientations"][f"{a}|{b}"] = [tail, head]
        out["confidence"][f"{a}|{b}"] = conf
    return out


def repair(
    cpdag: MPDAG, claims: dict[str, list[str]], confidence: dict[str, float | None]
) -> tuple[list, int]:
    """Drop claims in increasing confidence until the Meek closure succeeds."""
    ordered = sorted(claims, key=lambda k: (confidence.get(k) is None, confidence.get(k) or 0.0))
    kept = dict(claims)
    dropped = 0
    while kept and apply_orientations(cpdag, [tuple(v) for v in kept.values()]) is None:
        kept.pop(ordered[dropped])
        dropped += 1
    return [tuple(v) for v in kept.values()], dropped


def run_condition(
    cond: str, spec: dict, items: list[dict], graphs: dict[str, MPDAG], qhash: str
) -> tuple[dict, dict, dict]:
    """Answer ``items`` under one condition and freeze the knowledge it becomes.

    Returns the answers payload, the frozen per-network knowledge, and a
    summary. ``graphs`` holds the CPDAG each network's orientations are
    Meek-closed against; it is the true DAG's CPDAG on the ledger and the
    estimated one on the estimated panel, and nothing here reads which.
    """
    answers = []
    per_net: dict[str, dict] = collections.defaultdict(
        lambda: {
            "claims": {},
            "confidence": {},
            "declined": 0,
            "not_reached": 0,
            "invented": 0,
            "contradicts_order": 0,
            "parse_fail": 0,
            "asked": 0,
            "control": None,
        }
    )
    hits = 0
    t0 = time.perf_counter()
    for i, it in enumerate(items):
        rec = ask(spec["model"], it["prompt"])
        hits += rec["cache_hit"]
        parsed = parse(rec["response"], it["label_of"], it["pairs"])
        answers.append(
            {
                **{
                    k: it[k]
                    for k in (
                        "network",
                        "component",
                        "kind",
                        "naming",
                        "order_seed",
                        "prompt_sha256",
                    )
                },
                "model": spec["model"],
                "cache_key": rec["key"],
                **parsed,
            }
        )
        pn = per_net[it["network"]]
        if it["kind"] == "open":
            pn["claims"].update(parsed["orientations"])
            pn["confidence"].update(parsed["confidence"])
            pn["declined"] += len(parsed["declined"])
            pn["not_reached"] += len(parsed["not_reached"])
            pn["invented"] += parsed["invented"]
            pn["contradicts_order"] += parsed["contradicts_order"]
            pn["parse_fail"] += 0 if parsed["parse_ok"] else 1
            pn["asked"] += len(it["pairs"])
        else:
            truth = it["compelled_direction"]
            right = wrong = 0
            for k, (tail, head) in parsed["orientations"].items():
                if k in truth:
                    if [tail, head] == truth[k]:
                        right += 1
                    else:
                        wrong += 1
            pn["control"] = {
                "asked": len(it["pairs"]),
                "right": right,
                "wrong": wrong,
                "declined": len(parsed["declined"]),
                "not_reached": len(parsed["not_reached"]),
            }
        if i % 20 == 0:
            print(
                f"  {i}/{len(items)} items, {hits} cache hits, {time.perf_counter() - t0:.0f} s",
                flush=True,
            )

    # Meek closure, repair, freeze
    frozen = {}
    for net, pn in sorted(per_net.items()):
        cpdag = graphs[net]
        pre = [tuple(v) for v in pn["claims"].values()]
        consistent = apply_orientations(cpdag, pre) is not None
        post, dropped = (pre, 0) if consistent else repair(cpdag, pn["claims"], pn["confidence"])
        post = sorted(post)
        frozen[net] = {
            "k_pre": sorted(pre),
            "k": post,
            "k_sha256": hashlib.sha256(json.dumps(post).encode()).hexdigest(),
            "n_asked": pn["asked"],
            "n_asserted_pre": len(pre),
            "n_asserted": len(post),
            "meek_consistent_pre": consistent,
            "dropped_for_consistency": dropped,
            "declined": pn["declined"],
            "not_reached": pn["not_reached"],
            "invented": pn["invented"],
            "contradicts_order": pn["contradicts_order"],
            "parse_fail": pn["parse_fail"],
            "compelled_control": pn["control"],
            "confidence_quarantined": {k: v for k, v in pn["confidence"].items()},
        }
    payload = {
        "condition": cond,
        **spec,
        "family": FAMILY.get(spec["model"], "?"),
        "questionnaire_sha256": qhash,
        "cache_hits": hits,
        "items": answers,
    }
    summary = {
        "networks": len(frozen),
        "asked": sum(v["n_asked"] for v in frozen.values()),
        "asserted": sum(v["n_asserted"] for v in frozen.values()),
        "declined": sum(v["declined"] for v in frozen.values()),
        "not_reached": sum(v["not_reached"] for v in frozen.values()),
        "invented": sum(v["invented"] for v in frozen.values()),
        "contradicts_own_order": sum(v["contradicts_order"] for v in frozen.values()),
        "meek_inconsistent_networks": sum(
            1 for v in frozen.values() if not v["meek_consistent_pre"]
        ),
        "dropped_for_consistency": sum(v["dropped_for_consistency"] for v in frozen.values()),
        "control_right": sum(
            (v["compelled_control"] or {}).get("right", 0) for v in frozen.values()
        ),
        "control_wrong": sum(
            (v["compelled_control"] or {}).get("wrong", 0) for v in frozen.values()
        ),
        "cache_hits": hits,
        "seconds": round(time.perf_counter() - t0, 1),
    }
    return payload, frozen, summary


def main(argv: list[str] | None = None) -> None:
    """Answer the questionnaire under the requested conditions and freeze the knowledge."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--conditions", default="D_LLM")
    ap.add_argument("--networks", default="", help="comma-separated subset, for a pilot")
    ap.add_argument(
        "--questionnaire",
        default=str(ELICIT / "questionnaire.json"),
        help="the frozen questionnaire; its .sha256 sits beside it",
    )
    ap.add_argument(
        "--cpdag-source",
        default="",
        help="JSON of estimated CPDAGs to close against; by default the true DAG's CPDAG",
    )
    ap.add_argument("--out", default=str(ELICIT), help="output directory")
    args = ap.parse_args(argv)
    conds = [c for c in args.conditions.split(",") if c]
    only = {n for n in args.networks.split(",") if n}
    qpath = Path(args.questionnaire)
    out = Path(args.out)

    q = json.loads(qpath.read_text())
    qhash = qpath.with_suffix(".sha256").read_text().strip()
    graphs: dict[str, MPDAG] = {}
    if args.cpdag_source:
        from questionnaire_build import load_cpdag_source

        graphs = load_cpdag_source(Path(args.cpdag_source))
    else:
        for path in sorted(MODELS.iterdir()):
            parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
            if parsed.bidirected:
                continue
            graphs[parsed.name] = dag_to_cpdag(to_mpdag(parsed))

    out.mkdir(parents=True, exist_ok=True)
    knowledge_all = (
        json.loads((out / "knowledge.json").read_text())
        if (out / "knowledge.json").exists()
        else {}
    )
    for cond in conds:
        spec = CONDITIONS[cond]
        items = [
            it
            for it in q["items"]
            if it["naming"] == spec["naming"]
            and (it["kind"] == "compelled_control" or it["order_seed"] == spec["order_seed"])
            and (not only or it["network"] in only)
        ]
        print(
            f"\n== {cond}: {spec['model']}, {spec['naming']} names, "
            f"order seed {spec['order_seed']}: {len(items)} items",
            flush=True,
        )
        payload, frozen, summary = run_condition(cond, spec, items, graphs, qhash)
        (out / f"answers_{cond}.json").write_text(json.dumps(payload, indent=1))
        knowledge_all[cond] = {
            **spec,
            "family": FAMILY.get(spec["model"], "?"),
            "questionnaire_sha256": qhash,
            "networks": frozen,
            "true_dag_on_path": "false",
            "arm": "deployment",
        }
        (out / "knowledge.json").write_text(json.dumps(knowledge_all, indent=1, sort_keys=True))
        print("  " + json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
