"""How much of the census does the FULLY PROVED Theorem 1 cover?
   hypothesis (M): poss_de_{G0}(X) subset poss_de_H(X)
   sufficient condition (M0): the Meek cascade is empty (H = G0 + a->b)."""
import sys
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof")
from explore import *
from collections import defaultdict


def scan_C(C):
    p = C.shape[0]
    acc = defaultdict(int)
    for G0 in reachable_mpdags(C):
        Q = []
        for x in range(p):
            for y in range(p):
                if x != y:
                    z = parts(G0, x, y)
                    if z is not None: Q.append((x, y, z))
        if not Q: continue
        pd0 = {}
        for (a0, b0) in nonadjacent_pairs(G0):
            for (a, b) in ((a0, b0), (b0, a0)):
                H, info = bk_assert(G0, [(a, b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                raw = G0.copy(); raw[a, b] = 1; raw[b, a] = 0
                nocasc = np.array_equal(raw, H)
                pdH = {}
                for (x, y, (cn0, fb0, pa0, O0, _)) in Q:
                    z = parts(H, x, y)
                    if z is None: continue
                    acc["n"] += 1
                    if nocasc: acc["nocascade"] += 1
                    if x not in pd0: pd0[x] = frozenset(adjust.poss_de(G0, {x}))
                    if x not in pdH: pdH[x] = frozenset(adjust.poss_de(H, {x}))
                    if pd0[x] <= pdH[x]:
                        acc["M_holds"] += 1
                        if z[3] != O0: acc["M_holds_moved"] += 1
                    if z[3] != O0: acc["moved"] += 1
    return dict(acc)


def _job(args):
    kb, p = args
    return scan_C(np.frombuffer(kb, dtype=np.int8).reshape(p, p).copy())


def main():
    p = int(sys.argv[1]); nw = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    from multiprocessing import Pool
    Cs = cpdags(p); tot = defaultdict(int)
    with Pool(nw) as pool:
        for acc in pool.imap_unordered(_job, [(C.tobytes(), p) for C in Cs], chunksize=4):
            for k, v in acc.items(): tot[k] += v
    print("p=%d" % p, json.dumps(dict(sorted(tot.items())), indent=1))


if __name__ == "__main__":
    main()
