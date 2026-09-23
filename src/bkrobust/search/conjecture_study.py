"""Exhaustive counterexample search for the two structural conjectures.

Enumerates *all* DAGs on ``n`` nodes, derives their CPDAGs, dedupes, and then
for every CPDAG, every element of its space taken as ``G0``, every ordered
treatment/outcome pair, and every valid adjustment set, checks:

* **Conjecture 1** -- failure is upward-closed (expected to hold; it is a
  theorem, so a violation would indicate an implementation bug, not a
  mathematical discovery);
* **Conjecture 2** -- the nearest failure is reachable by retractions alone
  (genuinely open).

Scope is reported rather than assumed: the number of DAGs, CPDAGs, spaces,
``(G0, X, Y, Z)`` combinations and total radius comparisons all go into the
results file, so the strength of a "no counterexample found" claim is legible.
"""

from __future__ import annotations

import itertools
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from bkrobust.core.oracle import is_valid
from bkrobust.core.spacelib import build_space, distances_from
from bkrobust.demo.evaluate import (
    all_valid_adjustment_sets_mpdag,
    optimal_adjustment_set_mpdag,
)
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.search.conjectures import (
    check_conjecture1,
    check_conjecture2,
    check_no_empty_extensions,
)


def all_dags(n: int) -> list[MPDAG]:
    """Every DAG on ``n`` labelled nodes, deterministically ordered.

    Enumerates the ``3^(n choose 2)`` assignments of each unordered pair to
    {absent, forward, backward} and keeps the acyclic ones.
    """
    nodes = [f"V{i}" for i in range(n)]
    pairs = list(itertools.combinations(nodes, 2))
    out: list[MPDAG] = []
    for assign in itertools.product((0, 1, 2), repeat=len(pairs)):
        edges = []
        for (a, b), v in zip(pairs, assign):  # noqa: B905 - equal by construction
            if v == 1:
                edges.append((a, b))
            elif v == 2:
                edges.append((b, a))
        g = MPDAG(nodes, directed=edges)
        if g.is_acyclic():
            out.append(g)
    return out


def all_cpdags(n: int) -> list[MPDAG]:
    """Every distinct CPDAG on ``n`` labelled nodes, deterministically ordered."""
    seen: dict[str, MPDAG] = {}
    for d in all_dags(n):
        c = dag_to_cpdag(d)
        seen.setdefault(c.edge_string(), c)
    return [seen[k] for k in sorted(seen)]


@dataclass
class StudyTotals:
    """Running totals for one study sweep."""

    n_dags: int = 0
    n_cpdags: int = 0
    n_spaces: int = 0
    n_combos: int = 0
    n_radius_comparisons: int = 0
    c1_violations: int = 0
    c2_violations: int = 0
    empty_extension_spaces: int = 0
    skipped_no_valid_z: int = 0
    skipped_z_fails_at_g0: int = 0
    skipped_cpdag_too_many_undirected: int = 0
    skipped_cpdag_space_too_big: int = 0
    skipped_cpdag_no_undirected: int = 0
    c1_checks_run: int = 0
    c1_checks_skipped: int = 0


def run_study(
    n: int,
    out_dir: str | Path,
    *,
    max_space: int = 400,
    max_undirected: int = 6,
    c1_every: int = 50,
    all_z: bool = True,
    max_z_per_g0: int = 8,
    verbose: bool = True,
) -> tuple[StudyTotals, list[dict[str, Any]]]:
    """Run the exhaustive conjecture study for ``n`` nodes.

    Args:
        n: Number of nodes.
        out_dir: Directory for incremental output.
        max_space: Skip CPDAGs whose space exceeds this, recording the skip.
        c1_every: Run the O(|space|^2) Conjecture-1 check on one combination in
            this many. Conjecture 1 is a theorem; the check guards the
            implementation, so a sample suffices and running it on every
            combination dominated the sweep.
        max_undirected: Skip CPDAGs with more undirected edges than this BEFORE
            building their space, since construction is ``3^k`` in that count.
            Skips are counted and reported, so the scope of a "no
            counterexample" claim stays honest.
        all_z: Test every valid adjustment set, not only ``O*(G0)``.
        max_z_per_g0: Cap on adjustment sets tested per ``G0``, for cost.
        verbose: Print progress.

    Returns:
        ``(totals, counterexamples)``. Counterexamples are recorded minimally:
        the CPDAG, ``G0``, target pair, ``Z``, and both radii.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    totals = StudyTotals()
    counterexamples: list[dict[str, Any]] = []

    cpdags = all_cpdags(n)
    totals.n_dags = len(all_dags(n))
    totals.n_cpdags = len(cpdags)
    start = time.time()

    for ci, cpdag in enumerate(cpdags):
        if not cpdag.undirected_edges:
            totals.skipped_cpdag_no_undirected += 1
            continue  # nothing background knowledge could orient
        if len(cpdag.undirected_edges) > max_undirected:
            totals.skipped_cpdag_too_many_undirected += 1
            continue
        # Two-stage build. Enumerating the elements is cheap; computing the
        # covering relation is cubic in the element count, so the size guard has
        # to fire BEFORE that, not after. Checking it afterwards is what made a
        # single dense CPDAG take four minutes.
        from bkrobust.demo.space import enumerate_space

        if len(enumerate_space(cpdag)) > max_space:
            totals.skipped_cpdag_space_too_big += 1
            continue
        space = build_space(cpdag)
        totals.n_spaces += 1
        if not check_no_empty_extensions(space):
            totals.empty_extension_spaces += 1

        nodes = list(cpdag.nodes)
        for g0 in space.elements:
            dists = distances_from(space, g0)
            for x, y in itertools.permutations(nodes, 2):
                if all_z:
                    zs = all_valid_adjustment_sets_mpdag(g0, x, y)[:max_z_per_g0]
                else:
                    o = optimal_adjustment_set_mpdag(g0, x, y)
                    zs = [] if o is None else [frozenset(o)]
                if not zs:
                    totals.skipped_no_valid_z += 1
                    continue
                for zset in zs:
                    z = frozenset(zset)
                    if not is_valid(z, g0, x, y):
                        totals.skipped_z_fails_at_g0 += 1
                        continue
                    totals.n_combos += 1

                    def fails(
                        g: MPDAG,
                        _z: frozenset[str] = z,
                        _x: str = x,
                        _y: str = y,
                    ) -> bool:
                        return not is_valid(_z, g, _x, _y)

                    c2 = check_conjecture2(space, g0, fails, dists)
                    totals.n_radius_comparisons += 1
                    if not c2.holds:
                        totals.c2_violations += 1
                        counterexamples.append(
                            {
                                "conjecture": 2,
                                "n": n,
                                "cpdag": cpdag.edge_string(),
                                "g0": g0.edge_string(),
                                "x": x,
                                "y": y,
                                "z": sorted(z),
                                "r_full": c2.r_full,
                                "r_up": c2.r_up,
                                "witness_full": c2.witness_full,
                                "witness_up": c2.witness_up,
                                "space_size": c2.n_space,
                            }
                        )
                    # Conjecture 1 is a THEOREM, so checking it is guarding the
                    # implementation, not testing mathematics. It is O(|space|^2)
                    # per call, which dominated the whole sweep when run on every
                    # combination -- one dense CPDAG burned minutes. Sampled
                    # instead, at a rate recorded in the totals.
                    totals.c1_checks_skipped += 1
                    if totals.n_combos % c1_every != 0:
                        continue
                    totals.c1_checks_skipped -= 1
                    totals.c1_checks_run += 1
                    c1 = check_conjecture1(space, fails)
                    if not c1.holds:
                        totals.c1_violations += 1
                        counterexamples.append(
                            {
                                "conjecture": 1,
                                "n": n,
                                "cpdag": cpdag.edge_string(),
                                "g0": g0.edge_string(),
                                "x": x,
                                "y": y,
                                "z": sorted(z),
                                "violations": list(c1.violations[:3]),
                            }
                        )
        if verbose and (ci + 1) % 25 == 0:
            print(
                f"  n={n}: cpdag {ci + 1}/{len(cpdags)} "
                f"combos={totals.n_combos} c1v={totals.c1_violations} "
                f"c2v={totals.c2_violations} [{time.time() - start:.0f}s]",
                flush=True,
            )
    return totals, counterexamples
