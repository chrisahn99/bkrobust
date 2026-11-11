"""Vendored from experiments/ledger_sweep.py.

commit ba7ddc12c446f1087e58dc428720c84d33b1b622
git status of that file at copy time: A  experiments/ledger_sweep.py (staged, uncommitted; working
tree as it stood on 2026-11-06).

The constants and functions below are copied unchanged (extracted by name with ast), so E12 scores
rows exactly as the released ledger did. The imports are rewritten to the submitted tree plus the
two vendored modules; NO_RETRACTABLE_EDGES is absent from the submitted conventions and is defined
here with the working tree's value.
"""

from __future__ import annotations

import itertools
import math

import numpy as np

from bkrobust.benchmarks.measure import component_of
from bkrobust.demo.evaluate import is_valid_adjustment_set_dag, optimal_adjustment_set_dag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag
from bkrobust.mpdag_criterion.criterion import is_amenable

from vendor_audit import is_valid_adjustment_set_gac_dag

NO_RETRACTABLE_EDGES = "no_retractable_edges"

MAX_DEPTH = 3


SUBSET_BUDGET = 300


RUNGS = ("B2_free_rule", "B4_kg0", "B5_k", "B6_hop", "B7_claim")


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
