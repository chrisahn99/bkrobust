"""State rebuild of the released sweep, vendored from experiments/ledger_sweep.py.

commit ba7ddc12c446f1087e58dc428720c84d33b1b622
git status of that file at copy time: 'A  experiments/ledger_sweep.py' (staged; working tree of 2026-11-06).

``load`` and ``chance_level`` are copied from that file; ``knowledge_for_row`` restates the
knowledge plan of its ``main`` (lines 431-451 there) for one released row. The same rebuild is
used by experiments/rebuttal/E1_cpdag_check_on_lower_bound_rows/run.py, which verified it
against the released hashes.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

from bkrobust.benchmarks.describe import parse_file, to_mpdag
from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations

RELEASE_SEED = 20260913


def load(models: Path, name: str) -> tuple[MPDAG, MPDAG, str]:
    """Parse one network and return its DAG, CPDAG and file name (ledger_sweep.load)."""
    path = next(p for p in sorted(models.iterdir()) if p.name.startswith(name))
    parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
    dag = to_mpdag(parsed)
    return dag, dag_to_cpdag(dag), path.name


def chance_level(cpdag: MPDAG, base: list, rng: np.random.Generator, retries: int = 20) -> list | None:
    """Same edges as ``base``, orientations uniform among Meek-consistent choices."""
    for _ in range(retries):
        drawn = [(a, b) if rng.integers(2) else (b, a) for (a, b) in base]
        if apply_orientations(cpdag, drawn) is not None:
            return drawn
    return None


def knowledge_for_row(row: dict, dag: MPDAG, cpdag: MPDAG, know: dict, seed: int = RELEASE_SEED) -> list | None:
    """The knowledge the released sweep asserted for this row's arm on its network."""
    net, arm = row["network"], row["arm"]
    base_llm = [tuple(e) for e in (know["D_LLM"]["networks"].get(net) or {}).get("k", [])]
    if arm == "D_RAND":
        net_seed = int(hashlib.sha256(net.encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng([seed, net_seed, int(row["rep"])])
        return chance_level(cpdag, base_llm, rng)
    if arm == "A_TRUE":
        return [(a, b) if dag.is_directed_edge(a, b) else (b, a) for (a, b) in base_llm]
    if arm == "A_ORACLE":
        return sorted(knowledge_to_recover(dag, cpdag))
    if arm == "D_DEGEN":
        return []
    entry = know[arm]["networks"].get(net) or {}
    return [tuple(e) for e in entry.get("k", [])]
