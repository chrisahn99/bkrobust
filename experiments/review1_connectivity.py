"""Review round 1, weakness 6: is the covering graph of the knowledge-state space
connected?

"Hop distance on the covering graph is a metric" is trivially true on any
connected graph and false (as a metric) when it is not: a disconnected space has
undefined distances, and the paper's hedge "every triple whose three distances
are defined" concedes as much. This settles it by construction.

Sketch of the argument the check is testing: the CPDAG is the top of the order,
every element is a Meek closure of some consistent orientation set, and
retracting orientations walks upward, so every element is joined to the top.
Connectivity therefore follows from gradedness -- but gradedness rests on
Anti-Exchange Case B, which is verified and not proved, so the connectivity of
the *covering graph* is worth checking directly rather than inheriting.
"""
from __future__ import annotations

import itertools
import json
import pathlib
import sys
from collections import deque

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.space import bfs_distances, covering_pairs  # noqa: E402
from bkrobust.core.spacelib import enumerate_space_correct  # noqa: E402

OUT = pathlib.Path(__file__).resolve().parent.parent / "results" / "axisb2"


def all_cpdags(n: int) -> list[MPDAG]:
    """Every CPDAG on ``n`` labelled nodes, via the project's own space code."""
    from bkrobust.demo.example import dag_to_cpdag

    nodes = [f"V{i}" for i in range(n)]
    seen: dict[str, MPDAG] = {}
    pairs = list(itertools.combinations(nodes, 2))
    for mask in range(1 << len(pairs)):
        chosen = [pairs[i] for i in range(len(pairs)) if mask >> i & 1]
        for orient in itertools.product([0, 1], repeat=len(chosen)):
            edges = {(a, b) if o == 0 else (b, a) for (a, b), o in zip(chosen, orient)}
            d = MPDAG(nodes=set(nodes), directed=edges, undirected=set())
            if not d.is_acyclic():
                continue
            c = dag_to_cpdag(d)
            seen.setdefault(c.edge_string(), c)
    return list(seen.values())


def is_connected(space: list[MPDAG], covers: list[tuple[MPDAG, MPDAG]]) -> bool:
    adj: dict[str, set[str]] = {g.edge_string(): set() for g in space}
    for a, b in covers:
        adj[a.edge_string()].add(b.edge_string())
        adj[b.edge_string()].add(a.edge_string())
    start = next(iter(adj))
    seen = {start}
    q = deque([start])
    while q:
        for nb in adj[q.popleft()]:
            if nb not in seen:
                seen.add(nb)
                q.append(nb)
    return len(seen) == len(adj)


def main() -> None:
    report = {"scope": "all CPDAGs on 3 and 4 labelled nodes, exhaustive", "spaces": 0,
              "elements": 0, "disconnected_spaces": 0, "undefined_pairs": 0, "examples": []}
    for n in (3, 4):
        for c in all_cpdags(n):
            space = enumerate_space_correct(c)
            if len(space) < 2:
                report["spaces"] += 1
                report["elements"] += len(space)
                continue
            covers = covering_pairs(space)
            conn = is_connected(space, covers)
            dists = bfs_distances(space, covers) if not conn else None
            report["spaces"] += 1
            report["elements"] += len(space)
            if not conn:
                report["disconnected_spaces"] += 1
                if len(report["examples"]) < 3:
                    report["examples"].append(c.edge_string())
            if dists is not None:
                missing = sum(
                    1
                    for a in space
                    for b in space
                    if a != b and (a, b) not in dists
                )
                report["undefined_pairs"] += missing
    path = OUT / "review1_connectivity.json"
    path.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
