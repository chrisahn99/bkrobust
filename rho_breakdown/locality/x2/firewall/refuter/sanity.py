"""Pipeline sanity: (a) [G0] members really are consistent extensions of G0 and
Markov-equivalent to C; (b) HPM baseline -- O*(G0) is valid in EVERY D in [G0]."""
import sys, time
sys.path.insert(0, "/private/tmp/claude-501/-Users-josecosta-mugango/bceac0a4-93a7-4757-b3b3-a3d4e6c42d55/scratchpad/cfw_ref_A")
from cfwlib import *
n_ext = n_bad_ext = n_pair = n_bad_base = 0
for p in (4, 5):
    cp = {}
    for D in enum_dags(p):
        C = dag_to_cpdag(D); cp[C.tobytes()] = C
    for i, C in enumerate(list(cp.values())[:400]):
        ref = v_structures(C)
        for G0 in reachable_mpdags(C):
            dags = consistent_dag_extensions(G0, ref_vstructs=ref)
            for D in dags:
                n_ext += 1
                if not dag_agrees_with(D, G0) or v_structures(D) != ref or \
                   not np.array_equal(dag_to_cpdag(D), C):
                    n_bad_ext += 1
            for x in range(p):
                for y in range(p):
                    if x == y: continue
                    O0, np0, a0 = X.ostar_and_paths(G0, x, y)
                    if not a0: continue
                    for D in dags:
                        n_pair += 1
                        if not adjust.is_valid_adjustment_set(D, x, y, set(O0)):
                            n_bad_base += 1
print("class members checked", n_ext, "not-a-consistent-extension", n_bad_ext)
print("HPM baseline (O*(G0), D) pairs", n_pair, "INVALID", n_bad_base)
