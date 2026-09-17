"""Is the Dor-Tarsi check doing anything the cycle check does not, at k=1?
Count the joint distribution of (cycle, extendable) over single non-adjacent
stamps on Meek-closed MPDAGs."""
import sys, json, time
sys.path.insert(0, "/private/tmp/claude-501/-Users-josecosta-mugango/bceac0a4-93a7-4757-b3b3-a3d4e6c42d55/scratchpad/cfw_ref_A")
from cfwlib import *
from arm2 import gen
from collections import Counter
from multiprocessing import Pool

def one(args):
    seed, mode, kk = args
    rng = np.random.default_rng(seed); c = Counter()
    D0 = gen(rng, mode)
    if D0 is None: return c
    p = D0.shape[0]
    C = dag_to_cpdag(D0); ref = v_structures(C); G0 = C.copy()
    U = undirected_edges(C)
    if len(U) and rng.random() < 0.5:
        k = int(rng.integers(0, len(U)))
        for t in rng.permutation(len(U))[:k]:
            u, v = U[t]
            if not is_undirected(G0, u, v): continue
            a, b = (u, v) if rng.random() < 0.5 else (v, u)
            Hc = G0.copy(); Hc[b, a] = 0; Hc = meek_closure(Hc)
            if has_directed_cycle(Hc) or v_structures(Hc) != ref: continue
            G0 = Hc
    S = skeleton(G0)
    NA = [(a, b) for a in range(p) for b in range(a + 1, p) if S[a, b] == 0]
    if len(NA) < kk: return c
    Ks = ([[(a, b)] for (a0, b0) in NA for (a, b) in ((a0, b0), (b0, a0))] if kk == 1 else
          [[((int(NA[int(t)][0]), int(NA[int(t)][1])) if rng.random() < .5
             else (int(NA[int(t)][1]), int(NA[int(t)][0])))
            for t in rng.choice(len(NA), size=kk, replace=False)] for _ in range(40)])
    for K in Ks:
        H, info = x1_ops.bk_assert(G0, K)
        c[(bool(has_directed_cycle(H)), bool(x1_ops.pdag_extendable(H)),
           bool(info["conflict"]))] += 1
    return c

if __name__ == "__main__":
    n = int(sys.argv[1]); s0 = int(sys.argv[2]); kk = int(sys.argv[3])
    modes = ["er","sf","chain","dense"]
    tot = Counter()
    with Pool(4) as pool:
        for c in pool.imap_unordered(one, [(s0+i, modes[i%4], kk) for i in range(n)], chunksize=8):
            tot += c
    print("k =", kk, " (cycle, extendable, conflict) -> count")
    for k2, v in sorted(tot.items()): print("  ", k2, v)
    print("  total", sum(tot.values()))
