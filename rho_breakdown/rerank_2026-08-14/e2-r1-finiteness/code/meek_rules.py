"""
Rule-SELECTABLE Meek closure with a per-rule firing counter.

The pilot's graphs.meek_closure() hard-codes R1 -> R2 -> R3 -> R4 behind a single
boolean `orient` flag with no rule selector and no counter, so it can measure the
SIZE of a cascade but never ATTRIBUTE it. Attribution is exactly what the
R1-closed finiteness claim is about (Bang & Didelez, arXiv:2306.01638: tiered
knowledge needs only R1).

This module reproduces graphs.meek_closure() rule-for-rule -- the R1..R4 blocks
below are copied verbatim from graphs.py:122-157 -- and adds (a) a `rules` filter
and (b) a firing counter. Equality with the original on rules=('R1','R2','R3','R4')
is asserted by test_meek_rules.py.
"""
from itertools import combinations

import numpy as np

from graphs import (MeekFail, adjacent, has_directed_cycle, is_directed,
                    is_undirected, v_structures)

ALL_RULES = ("R1", "R2", "R3", "R4")


def meek_closure_rules(G, rules=ALL_RULES, counter=None):
    """Apply the SELECTED Meek rules to fixpoint. Returns a new amat.

    counter: optional dict; incremented per rule that fires.
    """
    G = G.copy()
    p = G.shape[0]
    use = set(rules)
    changed = True
    while changed:
        changed = False
        for a, b in [(a, b) for a in range(p) for b in range(p) if is_undirected(G, a, b)]:
            if not is_undirected(G, a, b):
                continue
            fired = None
            # R1: c -> a, a -- b, c not adj b  =>  a -> b
            if "R1" in use:
                for c in range(p):
                    if c in (a, b):
                        continue
                    if is_directed(G, c, a) and not adjacent(G, c, b):
                        fired = "R1"
                        break
            # R2: a -> c -> b, a -- b  =>  a -> b
            if fired is None and "R2" in use:
                for c in range(p):
                    if c in (a, b):
                        continue
                    if is_directed(G, a, c) and is_directed(G, c, b):
                        fired = "R2"
                        break
            # R3: a--c, a--d, c->b, d->b, c,d non-adjacent, a--b  =>  a -> b
            if fired is None and "R3" in use:
                cand = [c for c in range(p) if c not in (a, b)
                        and is_undirected(G, a, c) and is_directed(G, c, b)]
                for c, d in combinations(cand, 2):
                    if not adjacent(G, c, d):
                        fired = "R3"
                        break
            # R4: a--d, d->c, c->b, a--c(or adjacent), b,d non-adjacent, a--b => a -> b
            if fired is None and "R4" in use:
                for d in range(p):
                    if d in (a, b) or not is_undirected(G, a, d) or adjacent(G, b, d):
                        continue
                    for c in range(p):
                        if c in (a, b, d):
                            continue
                        if is_directed(G, d, c) and is_directed(G, c, b) and adjacent(G, a, c):
                            fired = "R4"
                            break
                    if fired is not None:
                        break
            if fired is not None:
                G[b, a] = 0
                changed = True
                if counter is not None:
                    counter[fired] = counter.get(fired, 0) + 1
    return G


def apply_background_knowledge_rules(C, K, rules=ALL_RULES, counter=None):
    """Perkovic et al. (UAI'17) Algorithm 1, with a selectable rule set.

    Identical to graphs.apply_background_knowledge() when rules == ALL_RULES.
    """
    G = C.copy()
    for (i, j) in K:
        if is_directed(G, j, i):
            raise MeekFail(f"conflict: {i}->{j} but graph has {j}->{i}")
        if G[i, j] == 0 and G[j, i] == 0:
            raise MeekFail(f"{i}-{j} not an edge of the CPDAG")
        G[j, i] = 0
        G = meek_closure_rules(G, rules=rules, counter=counter)
    if has_directed_cycle(G):
        raise MeekFail("directed cycle")
    if v_structures(G) != v_structures(C):
        raise MeekFail("new v-structure")
    return G


# --------------------------------------------------------------- orientation bookkeeping
def orientation_status(G, i, j):
    """1 if i->j, -1 if j->i, 0 if i--j, None if no edge."""
    if G[i, j] == 1 and G[j, i] == 0:
        return 1
    if G[j, i] == 1 and G[i, j] == 0:
        return -1
    if G[i, j] == 1 and G[j, i] == 1:
        return 0
    return None


def n_orientation_changes(G0, G1, U):
    """# edges of U whose orientation status differs between G0 and G1."""
    return sum(1 for (i, j) in U
               if orientation_status(G0, i, j) != orientation_status(G1, i, j))


def n_oriented(G, U):
    """# edges of U that are directed in G."""
    return sum(1 for (i, j) in U if orientation_status(G, i, j) != 0)


def cascade_size(G, U, K):
    """# edges of U oriented in G BEYOND those directly named in K.

    Non-negative by construction, unlike the pilot's
    n_extra_orientations = |dir(G)| - |dir(C)| - |K|, which double-subtracts when
    a statement of K names an edge that an earlier statement's closure already
    oriented -- a real case for tiered K, where K is the FULL cross-tier set.
    """
    named = {frozenset(e) for e in K}
    named_in_U = sum(1 for (i, j) in U if frozenset((i, j)) in named)
    return n_oriented(G, U) - named_in_U


def sym_diff_statements(K, Kp):
    """|K delta K'| over ORIENTED pairs, counting a reversed edge as 1 change,
    an added statement as 1, a removed statement as 1."""
    a = {tuple(e) for e in K}
    b = {tuple(e) for e in Kp}
    edges_a = {frozenset(e): e for e in a}
    edges_b = {frozenset(e): e for e in b}
    n = 0
    for k in set(edges_a) | set(edges_b):
        if k not in edges_a or k not in edges_b:
            n += 1                       # added or removed
        elif edges_a[k] != edges_b[k]:
            n += 1                       # reversed
    return n
