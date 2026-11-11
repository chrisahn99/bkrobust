"""Functions copied verbatim from experiments/rebuttal/E5_learned_cpdags_by_sample_size/run.py
(sha256 0e8d66763a8434f13ae2fe7f34c9e93ed821fc307145bd9088e4f08125d99df6, copied 2026-11-06):
constants NETWORKS, DRAND_SEED; stable_seed, canon, file_sha, repo_commit, load_inputs,
greedy_consistent, truthful, plans_transfer. Only the imports were trimmed to what these use.
"""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

from bkrobust.benchmarks.describe import parse_file, to_mpdag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from vendor_bkrobust import chance_level, graph_from_json, sha

REPO = Path(__file__).resolve().parents[3]
NETWORKS = ("ecoli70", "arth150", "magic-niab", "magic-irri")
DRAND_SEED = 20260913  # the stored sweep's --seed, so D_RAND reproduces


def stable_seed(name: str) -> int:
    """experiments/questionnaire_build.py::stable_seed, copied to avoid importing that module."""
    return int(hashlib.sha256(name.encode()).hexdigest()[:8], 16)


def canon(e) -> tuple[str, str]:
    return tuple(sorted(e))


def k_hash(k: list) -> str:
    return hashlib.sha256(json.dumps(sorted(k)).encode()).hexdigest()[:16]


def file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repo_commit() -> str:
    try:
        out = subprocess.run(
            ["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True, timeout=20
        )
        return out.stdout.strip() or "unknown"
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"



def load_inputs(data_dir: Path, models_dir: Path, networks: tuple[str, ...]) -> dict:
    """Everything a worker needs, read once. Fails loudly on a missing or mismatched file."""
    chat = json.loads((data_dir / "chat_e.json").read_text())
    know = json.loads((data_dir / "knowledge.json").read_text())
    frame = [json.loads(line) for line in (data_dir / "frame.jsonl").open() if line.strip()]
    learned = []
    for name in ("graphs.jsonl", "graphs_retry.jsonl"):
        p = data_dir / name
        if p.exists():
            learned += [json.loads(line) for line in p.open() if line.strip()]
    stored: dict[tuple, dict] = {}
    with (data_dir / "rows.csv").open() as fh:
        for r in csv.DictReader(fh):
            stored[(r["network"], r["arm"], int(r["rep"]), r["x"], r["y"])] = r
    dags, names = {}, {}
    for net in networks:
        path = models_dir / f"{net}.json"
        parsed = parse_file(path, file_sha(path))
        dags[net] = to_mpdag(parsed)
        names[net] = dict(parsed.name_map)
    graphs: dict[str, dict] = {}
    for net in networks:
        e = chat["networks"][net]
        g = graph_from_json(e)
        if sha(g) != e["sha256"]:
            raise SystemExit(f"{net}: chat_e.json reference CPDAG hash mismatch")
        graphs[f"{net}|reference"] = {"net": net, "method": "pc", "n": chat["meta"]["n"],
                                      "seed": e["seed"], "kind": "reference", "g": g}
        cp = dag_to_cpdag(dags[net])
        graphs[f"{net}|true_cpdag|true|-"] = {"net": net, "method": "true_cpdag", "n": None,
                                              "seed": None, "kind": "true_cpdag", "g": cp}
    for e in learned:
        if e["network"] not in networks:
            continue
        gid = f"{e['network']}|{e['method']}|{e['n']}|{e['seed']}"
        if gid in graphs:
            raise SystemExit(f"duplicate learned graph {gid}")
        graphs[gid] = {"net": e["network"], "method": e["method"], "n": e["n"], "seed": e["seed"],
                       "kind": "learned", "g": graph_from_json(e)}
    return {"chat": chat, "know": know, "frame": frame, "stored": stored, "dags": dags,
            "names": names, "graphs": graphs}



def greedy_consistent(g: MPDAG, claims: list) -> tuple[list, int]:
    """Keep claims in order, skipping any whose addition makes the Meek closure fail."""
    kept, dropped = [], 0
    for e in claims:
        if apply_orientations(g, [*kept, e]) is None:
            dropped += 1
        else:
            kept.append(e)
    return kept, dropped


def truthful(dag: MPDAG, edges) -> list:
    return [(a, b) if dag.is_directed_edge(a, b) else (b, a) for (a, b) in edges if dag.has_edge(a, b)]


def plans_stored(inp: dict, net: str, rand_reps: int) -> list[tuple]:
    """The stored sweep's arms on the reference CPDAG, exactly as pest_panel.stage_sweep built them."""
    dag, g = inp["dags"][net], inp["graphs"][f"{net}|reference"]["g"]
    base = [tuple(e) for e in inp["know"]["D_LLM"]["networks"].get(net, {}).get("k", [])]
    plans = [("D_LLM", 0, base, "false", {})]
    if base:
        for rep in range(rand_reps):
            kk = chance_level(g, base, np.random.default_rng([DRAND_SEED, stable_seed(net), rep]))
            if kk is not None:
                plans.append(("D_RAND", rep, kk, "false", {}))
    plans.append(("D_DEGEN", 0, [], "false", {}))
    plans.append(("A_TRUE", 0, truthful(dag, base), "via_K", {}))
    plans.append(("A_ORACLE", 0, truthful(dag, sorted(g.undirected_edges)), "via_K", {}))
    return plans


def plans_transfer(inp: dict, gid: str, rand_reps: int) -> list[tuple]:
    """The same arms on another graph: the supplier's claims restricted to the graph's open edges."""
    info = inp["graphs"][gid]
    net, g, dag = info["net"], info["g"], inp["dags"][info["net"]]
    stored = [tuple(e) for e in inp["know"]["D_LLM"]["networks"].get(net, {}).get("k", [])]
    on_graph = [e for e in stored if canon(e) in g.undirected_edges]
    base, dropped = greedy_consistent(g, on_graph)
    note = {"n_llm_stored": len(stored), "n_llm_off_graph": len(stored) - len(on_graph),
            "n_llm_dropped": dropped}
    plans = [("D_LLM", 0, base, "false", note)]
    if base:
        for rep in range(rand_reps):
            kk = chance_level(g, base, np.random.default_rng([DRAND_SEED, stable_seed(net), rep]))
            if kk is not None:
                plans.append(("D_RAND", rep, kk, "false", note))
    plans.append(("D_DEGEN", 0, [], "false", {}))
    plans.append(("A_TRUE", 0, truthful(dag, base), "via_K", note))
    oracle = truthful(dag, sorted(g.undirected_edges))
    plans.append(("A_ORACLE", 0, oracle, "via_K", {}))
    # the stored A_ORACLE is unusable whenever one truthful claim contradicts a wrongly learned
    # arrow; the repaired arm keeps the truthful claims that close, so every graph contributes
    kept, dropped_o = greedy_consistent(g, oracle)
    plans.append(("A_ORACLE_R", 0, kept, "via_K", {"n_oracle_dropped": dropped_o}))
    return plans
