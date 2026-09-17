import sys
sys.path.insert(0, "/private/tmp/claude-501/-Users-josecosta-mugango/bceac0a4-93a7-4757-b3b3-a3d4e6c42d55/scratchpad/cfw_ref_A")
from cfwlib import *
from collections import Counter
for p in (4, 5):
    cp = {}
    for D in enum_dags(p):
        C = dag_to_cpdag(D); cp[C.tobytes()] = C
    c = Counter(); n = 0
    lim = None if p == 4 else 1200
    for i, C in enumerate(list(cp.values())[:lim] if lim else cp.values()):
        for G0 in reachable_mpdags(C):
            S = skeleton(G0)
            NA = [(a, b) for a in range(p) for b in range(a+1, p) if S[a, b] == 0]
            for (a0, b0) in NA:
                for (a, b) in ((a0, b0), (b0, a0)):
                    H, info = x1_ops.bk_assert(G0, [(a, b)])
                    c[(bool(has_directed_cycle(H)), bool(x1_ops.pdag_extendable(H)))] += 1
                    n += 1
    print("p", p, "n_cpdag_scanned", len(cp) if not lim else min(lim, len(cp)), "trials", n)
    for k, v in sorted(c.items()): print("   (cycle,extendable)", k, v)
