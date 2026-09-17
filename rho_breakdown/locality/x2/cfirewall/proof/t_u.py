"""Claim U: pa_H(cn_H) cap poss_de_{G0}(X)  subset  poss_de_H(X).
Also claim U2: cn_H cap ... ; and the localisation of the failure."""
import sys
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof")
from explore import *
from collections import defaultdict


def scan_C(C):
    p = C.shape[0]
    acc = defaultdict(int); wit = defaultdict(list)
    for G0 in reachable_mpdags(C):
        pd0 = {x: frozenset(adjust.poss_de(G0, {x})) for x in range(p)}
        NA = nonadjacent_pairs(G0)
        Q = []
        for x in range(p):
            for y in range(p):
                if x != y:
                    z = parts(G0, x, y)
                    if z is not None: Q.append((x, y, z))
        if not NA or not Q: continue
        for (a0, b0) in NA:
            for (a, b) in ((a0, b0), (b0, a0)):
                H, info = bk_assert(G0, [(a, b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                pdH = {}
                for (x, y, (cn0, fb0, pa0, O0, _)) in Q:
                    z = parts(H, x, y)
                    if z is None: continue
                    cn1, fb1, pa1, O1, _ = z
                    acc["n"] += 1
                    if x not in pdH: pdH[x] = frozenset(adjust.poss_de(H, {x}))
                    U = (pa1 & pd0[x]) - pdH[x]
                    if U:
                        acc["U_viol"] += 1
                        acc["U_nodes"] += len(U)
                        if not (U <= fb1): acc["U_escapes_fb1"] += 1
                        if len(wit["U"]) < 8:
                            wit["U"].append((G0.tolist(), H.tolist(), a, b, x, y, sorted(U),
                                             sorted(cn1), sorted(fb1), sorted(pd0[x]), sorted(pdH[x])))
                    # U3: cn_H subset poss_de_H(X) trivially; check cn_0 vs pdH
                    if not ((cn1 & pd0[x]) <= pdH[x]): acc["cnH_notsub"] += 1
    return dict(acc), {k: v for k, v in wit.items()}


def _job(args):
    kb, p = args
    return scan_C(np.frombuffer(kb, dtype=np.int8).reshape(p, p).copy())


def main():
    p = int(sys.argv[1]); nw = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    from multiprocessing import Pool
    Cs = cpdags(p); print("cpdags", len(Cs), flush=True)
    tot = defaultdict(int); W = defaultdict(list)
    with Pool(nw) as pool:
        for acc, w in pool.imap_unordered(_job, [(C.tobytes(), p) for C in Cs], chunksize=4):
            for k, v in acc.items(): tot[k] += v
            for k, v in w.items():
                if len(W[k]) < 8: W[k].extend(v[:8 - len(W[k])])
    print(json.dumps(dict(sorted(tot.items())), indent=1))
    json.dump({"counts": dict(tot), "wit": {k: v for k, v in W.items()}},
              open("t_u_p%d.json" % p, "w"), indent=1)


if __name__ == "__main__":
    main()
