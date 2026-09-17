import sys, os
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/code")
import numpy as np
from itertools import product
from graphs import (dag_to_cpdag, meek_closure, v_structures, has_directed_cycle,
                    undirected_edges, skeleton, is_directed, dag_agrees_with,
                    directed_edges)
import adjust
import x2lib as X
from run_lemma import all_dags, reachable_mpdags

# ---------------- report state ----------------
POS, ZERO, REFUSE = "POS", "ZERO", "REFUSE"

def report_state(G, x, y):
    """(state, O) where state in {POS,ZERO,REFUSE}."""
    O, npaths, amen = X.ostar_and_paths(G, x, y)
    if npaths == 0:
        return ZERO, None
    if not amen:
        return REFUSE, None
    return POS, O

_valid_cache = {}
def valid_in_D(D, x, y, Z):
    k = (D.tobytes(), x, y, Z)
    v = _valid_cache.get(k)
    if v is None:
        v = adjust.is_valid_adjustment_set(D, x, y, set(Z))
        _valid_cache[k] = v
    return v

_cn_cache = {}
def cnD(D, x, y):
    k = (D.tobytes(), x, y)
    v = _cn_cache.get(k)
    if v is None:
        v = frozenset(adjust.causal_nodes(D, x, y))
        _cn_cache[k] = v
    return v

def sound(D, x, y, state, O):
    """Is the report issued from a graph SOUND against the true DAG D?
       POS  -> O* must be a valid adjustment set in D
       ZERO -> D must really have no causal path x~>y
       REFUSE -> a refusal is never unsound"""
    if state == REFUSE:
        return None
    if state == ZERO:
        return len(cnD(D, x, y)) == 0
    return valid_in_D(D, x, y, O)

def ident(state):
    return state in (POS, ZERO)
