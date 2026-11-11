"""Vendored from experiments/ledger_sweep.py (the sweep that wrote the
released per_query_records.csv).

commit ba7ddc12c446f1087e58dc428720c84d33b1b622
git status of that file at copy time: 'A  experiments/ledger_sweep.py' (staged, not committed).

``sha`` and ``chance_level`` are copied unchanged. ``load`` is ledger_sweep.load with the models
directory passed in instead of read from a module constant. ``rebuild_knowledge`` restates, per
row, the per-arm knowledge that ledger_sweep.main builds in its ``plans`` loop (same seeds, same
order); the release holds no A_TRUE or A_ORACLE rows, so those arms are not rebuilt here.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

from bkrobust.benchmarks.describe import parse_file, to_mpdag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations

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
)


def load(models: Path, name: str) -> tuple[MPDAG, MPDAG, dict[str, str]]:
    """Parse one network and return its DAG, CPDAG and raw-name map."""
    path = next(p for p in sorted(models.iterdir()) if p.name.startswith(name))
    parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
    dag = to_mpdag(parsed)
    return dag, dag_to_cpdag(dag), dict(parsed.name_map)


def sha(g: MPDAG) -> str:
    """A hash of a graph's edge string, for the provenance stamp."""
    return hashlib.sha256(g.edge_string().encode()).hexdigest()[:16]


def chance_level(
    cpdag: MPDAG, base: list, rng: np.random.Generator, retries: int = 20
) -> list | None:
    """Same edges as ``base``, orientations uniform among Meek-consistent choices."""
    for _ in range(retries):
        drawn = [(a, b) if rng.integers(2) else (b, a) for (a, b) in base]
        if apply_orientations(cpdag, drawn) is not None:
            return drawn
    return None


def rebuild_knowledge(net: str, arm: str, rep: int, cpdag: MPDAG, know: dict, seed: int) -> list | None:
    """The knowledge ledger_sweep.main gave this (network, arm, rep)."""
    base_llm = [tuple(e) for e in (know["D_LLM"]["networks"].get(net) or {}).get("k", [])]
    if arm == "D_RAND":
        if not base_llm:
            return None
        net_seed = int(hashlib.sha256(net.encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng([seed, net_seed, rep])
        return chance_level(cpdag, base_llm, rng)
    if arm == "D_DEGEN":
        return []
    if arm in DEPLOYMENT_ARMS:
        return [tuple(e) for e in (know[arm]["networks"].get(net) or {}).get("k", [])]
    raise ValueError(f"arm {arm!r} is not rebuilt by E13")
