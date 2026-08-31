"""
X1 — the spurious required-edge operator, the Dor-Tarsi extendability oracle,
and the three draw laws.

THE DESIGN PROBLEM, stated where the code lives.
------------------------------------------------
Meek's background-knowledge framework (Meek 1995; Perkovic, Kalisch & Maathuis
UAI'17 Alg. 1) decides "is K consistent with C?" = "does SOME DAG in MEC(C)
satisfy K?".  Every member of MEC(C) has skeleton(C).  A required edge on a pair
NON-ADJACENT in C is therefore satisfied by NO member of MEC(C): the predicate is
not merely broken for this class, it is TRIVIALLY FALSE for every such statement,
which is the same as saying it is uninformative.  `graphs.apply_background_knowledge`
encodes exactly that by raising `MeekFail("i-j not an edge of the CPDAG")`.

So the class S7 names cannot be expressed in the pilot's operator at all, and the
question "what fraction of spurious statements are caught by the consistency
check?" has a definitional answer (all of them, trivially) that says nothing
about the object a practitioner's tool actually builds.

`bk_assert` implements operator (ii) of the design brief — b-LOAD's semantics:
write the required edge into the working adjacency and re-stamp it after every
Meek closure so it is never overwritten (`mb_by_mb.py:495`, `g[mask] = bk[mask]`).
Nothing raises.  ARM S1 (operator (i), "a careful tool") is then a *derived view*:
the same graph, labelled by an INTRINSIC coherence predicate that does not
reference C's skeleton --

    consistent_S1 := (no conflict) and pdag_extendable(G) and not has_directed_cycle(G)

`pdag_extendable` is Dor & Tarsi (1992); AUDIT M4 requires it because the design's
original predicate (`v_structures(G) != v_structures(C)`) is contaminated by a
purely skeletal artefact: `v_structures` demands non-adjacent parents, so ADDING
the edge {a,b} DELETES every unshielded collider of C whose parent pair is exactly
{a,b}, with no incoherence in it at all (measured floor: 10.94% of N(C) draws on
`original`).  That test is retained under the name `vstruct_vs_C` as a labelled
diagnostic, decomposed four ways.
"""
from itertools import combinations

import numpy as np

from graphs import (adjacent, has_directed_cycle, is_directed, is_undirected,
                    meek_closure, skeleton, v_structures)

UNREACH = 1 << 20   # hop_dist's sentinel (run_arm1.py:49). NEVER pooled, NEVER averaged.


# ------------------------------------------------------------ Dor & Tarsi (1992)
def pdag_extendable(G):
    """Does the PDAG G admit a consistent DAG extension?

    Dor & Tarsi (1992), "A simple algorithm to construct a consistent extension
    of a partially oriented graph".  Repeatedly find a vertex x such that
      (1) x has no outgoing DIRECTED edge, and
      (2) every undirected neighbour y of x is adjacent to all other neighbours
          of x,
    orient all of x's undirected edges into x, and delete x.  Succeeds iff a
    consistent extension exists (same skeleton, all directed edges preserved, no
    new unshielded collider, acyclic).

    Skeleton-agnostic and intrinsic: it never looks at C.
    """
    p = G.shape[0]
    H = G.copy()
    alive = list(range(p))
    while alive:
        found = -1
        for x in alive:
            out = False
            for w in alive:
                if w != x and is_directed(H, x, w):
                    out = True
                    break
            if out:
                continue
            nbrs = [w for w in alive if w != x and adjacent(H, x, w)]
            ok = True
            for y in nbrs:
                if not is_undirected(H, x, y):
                    continue
                for z in nbrs:
                    if z != y and not adjacent(H, y, z):
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                found = x
                break
        if found < 0:
            return False
        alive.remove(found)
        H[found, :] = 0
        H[:, found] = 0
    return True


# ------------------------------------------------------------ operator S
def bk_assert(C, Kp, restamp="prefix"):
    """Assert every (i,j) in Kp as a REQUIRED edge i->j on C.

    b-LOAD's 'protect' semantics: NOTHING raises.  The pair may be non-adjacent,
    in which case the edge is ADDED and skeleton(G) grows.  Closure runs after
    each statement (pilot-matched); the required edges applied so far are then
    re-stamped, mirroring `g[mask] = background_knowledge[mask]`.

    restamp: "prefix" (b-LOAD: re-stamp the whole applied prefix, AUDIT A13#4)
             "current" (re-stamp only the statement just applied)

    Returns (G, info).  `info` records WHICH branch fired per statement, whether
    the re-stamp was a no-op under both scopes (gate G7), and the three coherence
    flags -- recorded, never used to reject.
    """
    G = C.copy()
    kinds = []
    n_add = n_ovr = n_ori = n_noop = 0
    conflict = False
    restamp_changed_prefix = False
    restamp_changed_current = False
    applied = []
    for (i, j) in Kp:
        if G[i, j] == 0 and G[j, i] == 0:
            G[i, j] = 1
            G[j, i] = 0
            n_add += 1
            kinds.append("added")
        elif is_directed(G, j, i):
            G[i, j] = 1
            G[j, i] = 0
            n_ovr += 1
            conflict = True
            kinds.append("overwrote")
        elif is_undirected(G, i, j):
            G[j, i] = 0
            n_ori += 1
            kinds.append("oriented")
        else:                                   # already i -> j
            n_noop += 1
            kinds.append("noop")
        applied.append((i, j))
        Gc = meek_closure(G)
        Gpre = Gc.copy()
        for (a, b) in applied:
            Gpre[a, b] = 1
            Gpre[b, a] = 0
        Gcur = Gc.copy()
        Gcur[i, j] = 1
        Gcur[j, i] = 0
        if not np.array_equal(Gpre, Gc):
            restamp_changed_prefix = True
        if not np.array_equal(Gcur, Gc):
            restamp_changed_current = True
        G = Gpre if restamp == "prefix" else Gcur
    info = dict(kinds=kinds, n_added=n_add, n_overwrote=n_ovr,
                n_oriented=n_ori, n_noop=n_noop, conflict=conflict,
                restamp_changed_prefix=restamp_changed_prefix,
                restamp_changed_current=restamp_changed_current,
                cycle=bool(has_directed_cycle(G)),
                vstruct_vs_C=bool(v_structures(G) != v_structures(C)),
                skel_changed=bool(not np.array_equal(skeleton(G), skeleton(C))))
    return G, info


def bk_assert_stampall(C, Kp):
    """Gate G8 comparator: stamp ALL of K' first, then close once, then re-stamp.
    b-LOAD's `initialize_background_knowledge` builds the whole matrix up front."""
    G = C.copy()
    for (i, j) in Kp:
        G[i, j] = 1
        G[j, i] = 0
    G = meek_closure(G)
    for (i, j) in Kp:
        G[i, j] = 1
        G[j, i] = 0
    return G


def vstruct_decomposition(C, G, added_pairs):
    """AUDIT M4.3: decompose the `vstruct_vs_C` diagnostic.

    `collider_of_C_shielded` -- an unshielded collider (a,m,b) of C whose parent
    pair {a,b} is exactly a pair the operator made adjacent.  Pure skeletal
    artefact, no incoherence.
    `new_unshielded` -- a v-structure of G that is not one of C.
    """
    vc, vg = v_structures(C), v_structures(G)
    lost = vc - vg
    gained = vg - vc
    addset = {frozenset(pr) for pr in added_pairs}
    shielded = sum(1 for (a, m, b) in lost if frozenset((a, b)) in addset)
    return dict(n_lost=len(lost), n_lost_shielded_by_added=shielded,
                n_lost_other=len(lost) - shielded, n_gained=len(gained))


# ------------------------------------------------------------ draw laws
def nonadjacent_pairs(C):
    """N(C) = unordered pairs non-adjacent in skeleton(C).  skeleton(C) == skeleton(D)."""
    S = skeleton(C)
    p = C.shape[0]
    return [(a, b) for a in range(p) for b in range(a + 1, p) if S[a, b] == 0]


def draw_suni(C, k, rng):
    """AUDIT M9.1: k pairs WITHOUT REPLACEMENT over UNORDERED pairs of N(C),
    each oriented by an independent fair coin from the same stream.

    Fixes the design's "ordered pair without replacement", which (a) made the
    eligibility rule ambiguous (|N(C)|>=4 vs >=2, disagreeing on 7.6% of
    `original`) and (b) permitted drawing both (a,b) and (b,a), so the second
    assertion silently overwrote the first."""
    N = nonadjacent_pairs(C)
    if len(N) < k:
        return None
    idx = rng.choice(len(N), size=k, replace=False)
    out = []
    for t in idx:
        a, b = N[int(t)]
        out.append((int(a), int(b)) if rng.random() < 0.5 else (int(b), int(a)))
    return out


def bload_pool(D, x, y):
    """b-LOAD's own candidate pool, verbatim (background_knowledge.py:166-170):
        [(u,v) for u in target_set for v in nodes if u != v and (u,v) not in true_edges]
    ORACLE-guaranteed-false sampler (AUDIT M14i): it is defined by (a,b) not in E(D)."""
    p = D.shape[0]
    return [(int(a), int(b)) for a in (x, y) for b in range(p)
            if b != a and D[a, b] != 1]


def draw_sloc(D, k, x, y, rng):
    pool = bload_pool(D, x, y)
    if len(pool) < k:
        return None
    idx = rng.choice(len(pool), size=k, replace=False)
    return [pool[int(t)] for t in idx]


def label_statement(C, D, a, b, x, y):
    """ORACLE labels (AUDIT M14ii): the analyst cannot perform this partition, so
    a stratified rate may NEVER be presented as a triage rule."""
    S = skeleton(C)
    if S[a, b] == 0:
        kind = "pure_spurious"
    elif D[b, a] == 1:
        kind = "reversal_of_true_edge"
    elif D[a, b] == 1:
        kind = "restatement_of_true_edge"
    else:
        kind = "reorient_cpdag_edge"
    return kind, bool({a, b} == {x, y})


# ------------------------------------------------------------ hop
def hop_dist_from(Gsk, srcs, p):
    """BFS on a skeleton.  Sentinel UNREACH for unreachable vertices."""
    dist = np.full(p, UNREACH, dtype=int)
    frontier = list(srcs)
    for v in frontier:
        dist[v] = 0
    while frontier:
        nxt = []
        for v in frontier:
            for w in np.flatnonzero(Gsk[v]):
                if dist[w] > dist[v] + 1:
                    dist[w] = dist[v] + 1
                    nxt.append(int(w))
        frontier = nxt
    return dist


def stmt_hop(dist, a, b):
    return int(min(dist[a], dist[b]))
