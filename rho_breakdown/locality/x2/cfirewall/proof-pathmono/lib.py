import sys, importlib.util
import numpy as np
BASE = "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious"
sys.path.insert(0, BASE + "/x2/code")
from graphs import (consistent_dag_extensions, dag_to_cpdag, directed_edges,
                    has_directed_cycle, meek_closure, random_dag, skeleton,
                    undirected_edges, v_structures, is_directed, is_undirected,
                    adjacent)
import adjust
import x2lib as X
from run_lemma import all_dags, reachable_mpdags
_spec = importlib.util.spec_from_file_location("x1_ops", BASE + "/x1/code/x1_ops.py")
x1_ops = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(x1_ops)
bk_assert = x1_ops.bk_assert
pdag_extendable = x1_ops.pdag_extendable
nonadjacent_pairs = x1_ops.nonadjacent_pairs

def dir_set(G):
    return set(directed_edges(G))

def in_class(D, H):
    """Is DAG D a consistent extension of PDAG H? same skeleton, dir(H) subset, same vstructs."""
    if not np.array_equal(skeleton(D), skeleton(H)): return False
    if has_directed_cycle(D): return False
    for (i,j) in directed_edges(H):
        if not is_directed(D,i,j): return False
    return v_structures(D) == v_structures(H)
