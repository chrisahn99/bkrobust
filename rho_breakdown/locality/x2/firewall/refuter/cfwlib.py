"""Strong-form refutation harness for C-FIREWALL (agent A)."""
import sys, json, time, importlib.util
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/code")
import numpy as np
from itertools import product
from collections import defaultdict

from graphs import (dag_to_cpdag, meek_closure, v_structures, has_directed_cycle,
                    undirected_edges, skeleton, is_directed, is_undirected,
                    random_dag, consistent_dag_extensions, dag_agrees_with)
import adjust
import x2lib as X

_spec = importlib.util.spec_from_file_location(
    "x1_ops_A", "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x1/code/x1_ops.py")
x1_ops = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(x1_ops)


def enum_dags(p):
    idx = [(i, j) for i in range(p) for j in range(i + 1, p)]
    for state in product([0, 1, 2], repeat=len(idx)):
        D = np.zeros((p, p), dtype=np.int8)
        for (i, j), s in zip(idx, state):
            if s == 1: D[i, j] = 1
            elif s == 2: D[j, i] = 1
        A = D.astype(bool); R = A.copy()
        for _ in range(p): R = R | (R @ A)
        if np.any(np.diag(R)): continue
        yield D


def reachable_mpdags(C, cap=0):
    ref = v_structures(C)
    seen = {C.tobytes(): C}; frontier = [C]
    while frontier:
        if cap and len(seen) > cap: return None
        nxt = []
        for G in frontier:
            for (u, v) in undirected_edges(G):
                for (a, b) in ((u, v), (v, u)):
                    H = G.copy(); H[b, a] = 0; H = meek_closure(H)
                    if has_directed_cycle(H) or v_structures(H) != ref: continue
                    kb = H.tobytes()
                    if kb not in seen:
                        seen[kb] = H; nxt.append(H)
        frontier = nxt
    return list(seen.values())


_vcache = {}
def valid_in(D, x, y, Z):
    k = (D.tobytes(), x, y, Z)
    v = _vcache.get(k)
    if v is None:
        v = adjust.is_valid_adjustment_set(D, x, y, set(Z))
        _vcache[k] = v
    return v


def scan_G0(G0, st, hits, max_hits=12, ref=None):
    """All queries x all non-adjacent stamped edges x ALL D in [G0]."""
    p = G0.shape[0]
    S = skeleton(G0)
    NA = [(a, b) for a in range(p) for b in range(a + 1, p) if S[a, b] == 0]
    if not NA: return
    if ref is None: ref = v_structures(G0)
    dags = consistent_dag_extensions(G0, ref_vstructs=ref)
    if not dags: return
    st["n_G0"] += 1
    st["sum_class_size"] += len(dags)
    for x in range(p):
        for y in range(p):
            if x == y: continue
            O0, np0, amen0 = X.ostar_and_paths(G0, x, y)
            if not amen0: continue
            st["n_query_pos"] += 1
            for (a0, b0) in NA:
                for (a, b) in ((a0, b0), (b0, a0)):
                    H, info = x1_ops.bk_assert(G0, [(a, b)])
                    ok = ((not info["conflict"]) and (not has_directed_cycle(H))
                          and x1_ops.pdag_extendable(H))
                    st["n_trial"] += 1
                    if not ok:
                        st["n_rejected"] += 1
                        continue
                    st["n_stamped"] += 1
                    O1, np1, amen1 = X.ostar_and_paths(H, x, y)
                    if not amen1:
                        st["n_abort"] += 1
                        st["n_abort_zero" if np1 == 0 else "n_abort_unamen"] += 1
                        continue
                    st["n_pos"] += 1
                    moved = (O1 != O0)
                    if moved: st["n_moved"] += 1
                    bad = [D for D in dags if not valid_in(D, x, y, O1)]
                    st["n_Dcheck"] += len(dags)
                    if bad:
                        st["n_violation"] += 1
                        if moved: st["n_violation_moved"] += 1
                        if len(hits) < max_hits:
                            hits.append(dict(G0=G0.tolist(), x=int(x), y=int(y),
                                             edge=[int(a), int(b)],
                                             O0=sorted(int(z) for z in O0),
                                             O1=sorted(int(z) for z in O1),
                                             H=H.tolist(), moved=bool(moved),
                                             n_bad=len(bad), n_class=len(dags),
                                             Dbad=[D.tolist() for D in bad[:3]]))


def new_stats():
    return defaultdict(int)
