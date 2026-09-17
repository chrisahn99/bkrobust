"""Exploration harness: census over p, recording set-relations between
G0 / H / D objects, to find the invariant that carries the proof."""
import sys, json, importlib.util
from collections import defaultdict
import numpy as np

BASE = "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious"
sys.path.insert(0, BASE + "/x2/code")
from graphs import (consistent_dag_extensions, dag_to_cpdag, has_directed_cycle,
                    skeleton, v_structures, is_directed, is_undirected, adjacent,
                    meek_closure, random_dag, undirected_edges)
import adjust
import x2lib as X
from run_lemma import all_dags, reachable_mpdags

_spec = importlib.util.spec_from_file_location("x1_ops", BASE + "/x1/code/x1_ops.py")
x1_ops = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(x1_ops)
bk_assert = x1_ops.bk_assert; pdag_extendable = x1_ops.pdag_extendable
nonadjacent_pairs = x1_ops.nonadjacent_pairs


def parts(G, x, y):
    """(cn, forb, pa(cn), O*) or None if not amenable."""
    paths = X.pcp_capped(G, x, y)
    if not paths: return None
    if not all(is_directed(G, q[0], q[1]) for q in paths): return None
    cn = set()
    for q in paths: cn.update(q[1:])
    fb = adjust.poss_de(G, cn) | {x}
    pa = adjust.parents_of_set(G, cn)
    return frozenset(cn), frozenset(fb), frozenset(pa), frozenset(pa - fb), paths


def rel(A, B):
    if A == B: return "eq"
    if A > B: return "sup"
    if A < B: return "sub"
    return "inc"


def cpdags(p):
    d = {}
    for D in all_dags(p):
        C = dag_to_cpdag(D)
        d.setdefault(C.tobytes(), C)
    return list(d.values())
