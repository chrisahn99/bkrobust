"""Closed-form optimal adjustment set on an amenable MPDAG.

Henckel, Perkovic & Maathuis (2022, JRSS-B), Section 5, and Witte et al.
(2020): for an MPDAG G amenable relative to (X, Y),

    cn(X, Y, G)  = nodes other than X on proper causal paths X -> ... -> Y
                   (in an amenable MPDAG these are the directed paths),
    forb(X, Y, G) = possde(cn(X, Y, G), G) u {X},
    O(X, Y, G)   = pa(cn(X, Y, G), G) \\ forb(X, Y, G),

with pa the definite (directed) parents. The submitted code reads O off by
enumerating DAG extensions and returns None when extensions disagree; this
module gives the graphical set and is checked against that enumeration
(E13 b).
"""

from __future__ import annotations

from bkrobust.demo.graph import MPDAG
from bkrobust.mpdag_criterion import is_amenable, possible_descendants


def _directed_desc(g: MPDAG, x: str) -> set[str]:
    out, stack = set(), [x]
    while stack:
        v = stack.pop()
        for c in g.children(v):
            if c not in out:
                out.add(c)
                stack.append(c)
    return out


def _directed_anc(g: MPDAG, y: str) -> set[str]:
    out, stack = set(), [y]
    while stack:
        v = stack.pop()
        for p in g.parents(v):
            if p not in out:
                out.add(p)
                stack.append(p)
    return out


def causal_nodes_directed(g: MPDAG, x: str, y: str) -> set[str]:
    """Nodes (other than x) on directed paths x -> ... -> y that avoid x after the start."""
    sub_nodes = [v for v in g.nodes if v != x]
    # proper paths: only the first node is x, so compute reachability in G minus x
    # from the children of x.
    desc = set()
    stack = list(g.children(x))
    desc.update(stack)
    while stack:
        v = stack.pop()
        for c in g.children(v):
            if c != x and c not in desc:
                desc.add(c)
                stack.append(c)
    if y not in desc:
        return set()
    anc = _directed_anc(g, y) | {y}
    del sub_nodes
    return {v for v in desc if v in anc}


def optimal_set_closed_form(g: MPDAG, x: str, y: str) -> tuple[frozenset[str] | None, str]:
    """(O, status). status in {"ok", "not_amenable", "y_not_descendant"}."""
    if not is_amenable(g, x, y):
        return None, "not_amenable"
    cn = causal_nodes_directed(g, x, y)
    if not cn:
        return None, "y_not_descendant"
    forb = set(possible_descendants(g, cn)) | {x}
    pa = set()
    for c in cn:
        pa |= g.parents(c)
    return frozenset(pa - forb), "ok"
