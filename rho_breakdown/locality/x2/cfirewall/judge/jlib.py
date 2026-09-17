"""Judge's independent harness. Uses ONLY the project's own graph/adjust primitives."""
import sys, importlib.util
BASE = "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious"
sys.path.insert(0, BASE + "/x2/code")
import numpy as np
from itertools import product, combinations
from graphs import (consistent_dag_extensions, dag_to_cpdag, directed_edges,
                    has_directed_cycle, meek_closure, random_dag, skeleton,
                    undirected_edges, v_structures, adjacent, is_directed,
                    is_undirected, common_orientation_graph)
import adjust
from adjust import (poss_de, causal_nodes, forb, parents_of_set,
                    possibly_causal_paths, is_valid_adjustment_set,
                    optimal_adjustment_set, possibly_directed_neighbors)

_spec = importlib.util.spec_from_file_location("x1_ops", BASE + "/x1/code/x1_ops.py")
x1_ops = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(x1_ops)
bk_assert = x1_ops.bk_assert
pdag_extendable = x1_ops.pdag_extendable
nonadjacent_pairs = x1_ops.nonadjacent_pairs

sys.path.insert(0, BASE + "/x2/code")
from run_lemma import all_dags, reachable_mpdags


def ostar(G, x, y):
    """(O*, npaths, amenable, cn, fb). Re-implemented from the stated definitions."""
    paths = possibly_causal_paths(G, x, y)
    npaths = len(paths)
    if npaths == 0:
        return None, 0, False, set(), {x}
    if not all(is_directed(G, p[0], p[1]) for p in paths):
        cn = set()
        for p in paths: cn.update(p[1:])
        fb = (poss_de(G, cn) | {x}) if cn else {x}
        return None, npaths, False, cn, fb
    cn = set()
    for p in paths: cn.update(p[1:])
    fb = (poss_de(G, cn) | {x}) if cn else {x}
    O = parents_of_set(G, cn) - fb
    return frozenset(O), npaths, True, cn, fb


def cpdags(p):
    seen = {}
    for D in all_dags(p):
        C = dag_to_cpdag(D)
        seen.setdefault(C.tobytes(), C)
    return list(seen.values())


def dplus(D, a, b):
    Dp = D.copy(); Dp[a, b] = 1; Dp[b, a] = 0
    return Dp


def lifts(D, a, b, H):
    """Does D+ = D + (a->b) lie in [H]?  (bridge hypothesis of the LIFT lemma)"""
    Dp = dplus(D, a, b)
    if has_directed_cycle(Dp):
        return False, Dp, "cycle"
    for (i, j) in directed_edges(H):
        if not is_directed(Dp, i, j):
            return False, Dp, "orient"
    if v_structures(Dp) != v_structures(H):
        return False, Dp, "vstruct"
    return True, Dp, "ok"


def pbd_graph(D, x, y, cnD):
    Dp = D.copy()
    for w in cnD:
        if is_directed(Dp, x, w):
            Dp[x, w] = 0
    return Dp


def reach_outside_forb(D, x, y, cnD, fbD):
    """Nodes reachable from x in skeleton(D^pbd), never entering forb_D (x itself allowed)."""
    Dp = pbd_graph(D, x, y, cnD)
    S = skeleton(Dp)
    p = D.shape[0]
    seen = {x}; st = [x]
    while st:
        v = st.pop()
        for w in range(p):
            if S[v, w] and w not in seen and w not in fbD:
                seen.add(w); st.append(w)
    return seen


def certL(D, x, y, Z):
    """Z n forb = {} and Z contains every pa(cn)\forb node REACHABLE from x outside forb."""
    cnD = causal_nodes(D, x, y)
    fbD = (poss_de(D, cnD) | {x}) if cnD else {x}
    if set(Z) & fbD: return False
    if not cnD: return False
    R = reach_outside_forb(D, x, y, cnD, fbD)
    need = (parents_of_set(D, cnD) - fbD) & R
    return need <= set(Z)


def certC(D, x, y, Z):
    p = D.shape[0]
    paX = {u for u in range(p) if is_directed(D, u, x)}
    deX = poss_de(D, {x})
    return paX <= set(Z) and not (set(Z) & deX)
