"""Functions E22 needs from E5, copied verbatim on 2026-11-06.

Source: experiments/rebuttal/E5_learned_cpdags_by_sample_size/run.py
    sha256 at copy time: 0e8d66763a8434f13ae2fe7f34c9e93ed821fc307145bd9088e4f08125d99df6
Copied: canon, k_hash, file_sha, repo_commit, greedy_consistent, truthful, and the learned-graph part of
load_inputs (as load_learned, cut to what E22 reads). The A_ORACLE_R plan of plans_transfer is repeated in
oracle_r_knowledge with the same two calls.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from vendor_bkrobust import graph_from_json


def canon(e) -> tuple[str, str]:
    return tuple(sorted(e))


def k_hash(k: list) -> str:
    return hashlib.sha256(json.dumps(sorted(k)).encode()).hexdigest()[:16]


def file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repo_commit(repo: Path) -> str:
    try:
        out = subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "HEAD"], capture_output=True, text=True, timeout=20
        )
        return out.stdout.strip() or "unknown"
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"


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


def oracle_r_knowledge(dag: MPDAG, g: MPDAG) -> tuple[list, int]:
    """plans_transfer's A_ORACLE_R arm: truthful claims on every open edge, kept while they close."""
    oracle = truthful(dag, sorted(g.undirected_edges))
    return greedy_consistent(g, oracle)


def load_learned(data_dir: Path, networks) -> dict:
    """load_inputs' learned-graph loop: gid -> record with the MPDAG under 'g'."""
    learned = []
    for name in ("graphs.jsonl", "graphs_retry.jsonl"):
        p = data_dir / name
        if p.exists():
            learned += [json.loads(line) for line in p.open() if line.strip()]
    graphs: dict[str, dict] = {}
    for e in learned:
        if e["network"] not in networks:
            continue
        gid = f"{e['network']}|{e['method']}|{e['n']}|{e['seed']}"
        if gid in graphs:
            raise SystemExit(f"duplicate learned graph {gid}")
        graphs[gid] = {"net": e["network"], "method": e["method"], "n": e["n"], "seed": e["seed"],
                       "kind": "learned", "g": graph_from_json(e)}
    return graphs
