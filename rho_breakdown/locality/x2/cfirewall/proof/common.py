import sys, importlib.util
BASE = "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious"
sys.path.insert(0, BASE + "/x2/code")
import numpy as np
from graphs import (dag_agrees_with, consistent_dag_extensions, dag_to_cpdag, directed_edges,
                    undirected_edges, has_directed_cycle, meek_closure,
                    random_dag, skeleton, v_structures, is_directed, is_undirected,
                    adjacent)
import adjust
import x2lib as X
from run_lemma import all_dags, reachable_mpdags

_spec = importlib.util.spec_from_file_location("x1_ops", BASE + "/x1/code/x1_ops.py")
x1_ops = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(x1_ops)
bk_assert = x1_ops.bk_assert
pdag_extendable = x1_ops.pdag_extendable
nonadjacent_pairs = x1_ops.nonadjacent_pairs

def acyclic(D):
    return not has_directed_cycle(D)

def preserves(W, H):
    """every directed edge of H is directed the same way in W"""
    for (i,j) in directed_edges(H):
        if not is_directed(W, i, j):
            return False
    return True
