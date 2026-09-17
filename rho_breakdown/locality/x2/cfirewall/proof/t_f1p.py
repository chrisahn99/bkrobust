"""Test F1' : O*(H) cap poss_de_{G0}(X) = empty   (stronger than F1).
Also: O*(G) cap poss_de_G(X) = empty for any amenable MPDAG G."""
import sys
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof")
from explore import *
from collections import defaultdict


def scan_C(C):
    p = C.shape[0]
    acc = defaultdict(int); wit = []
    for G0 in reachable_mpdags(C):
        pdX = {x: frozenset(adjust.poss_de(G0, {x})) for x in range(p)}
        NA = nonadjacent_pairs(G0)
        Q = []
        for x in range(p):
            for y in range(p):
                if x != y:
                    z = parts(G0, x, y)
                    if z is not None:
                        Q.append((x, y, z))
                        acc["G_n"] += 1
                        if z[3] & pdX[x]: acc["G_ownPD_viol"] += 1
        if not NA or not Q: continue
        for (a0, b0) in NA:
            for (a, b) in ((a0, b0), (b0, a0)):
                H, info = bk_assert(G0, [(a, b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                pdXH = {}
                for (x, y, (cn0, fb0, pa0, O0, _)) in Q:
                    z = parts(H, x, y)
                    if z is None: continue
                    cn1, fb1, pa1, O1, _ = z
                    acc["n"] += 1
                    if O1 & pdX[x]:
                        acc["F1p_viol"] += 1
                        if len(wit) < 6:
                            wit.append((G0.tolist(), H.tolist(), a, b, x, y,
                                        sorted(O1), sorted(pdX[x])))
                    if x not in pdXH: pdXH[x] = frozenset(adjust.poss_de(H, {x}))
                    if O1 & pdXH[x]: acc["F1pH_viol"] += 1
    return dict(acc), wit


def _job(args):
    kb, p = args
    return scan_C(np.frombuffer(kb, dtype=np.int8).reshape(p, p).copy())


def main():
    p = int(sys.argv[1]); nw = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    from multiprocessing import Pool
    Cs = cpdags(p); print("cpdags", len(Cs), flush=True)
    tot = defaultdict(int); W = []
    with Pool(nw) as pool:
        for acc, w in pool.imap_unordered(_job, [(C.tobytes(), p) for C in Cs], chunksize=4):
            for k, v in acc.items(): tot[k] += v
            if len(W) < 6: W.extend(w[:6 - len(W)])
    print(json.dumps(dict(sorted(tot.items())), indent=1))
    json.dump({"counts": dict(tot), "wit": W}, open("t_f1p_p%d.json" % p, "w"), indent=1)


if __name__ == "__main__":
    main()
