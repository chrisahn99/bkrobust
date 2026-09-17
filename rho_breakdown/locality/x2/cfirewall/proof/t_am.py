"""Lemma AM:  G Meek-closed and AMENABLE rel (X,Y)  =>
     cn_G(X,Y) = (poss_de_G(X) cap poss_an_G(Y)) \\ {X}."""
import sys
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof")
from explore import *
from collections import defaultdict


def poss_an(G, S):
    return adjust.poss_de(np.ascontiguousarray(G.T), S)


def check(G, tag, acc, wit):
    p = G.shape[0]
    for x in range(p):
        pdx = frozenset(adjust.poss_de(G, {x}))
        for y in range(p):
            if x == y: continue
            paths = X.pcp_capped(G, x, y)
            if not paths: continue
            amen = all(is_directed(G, q[0], q[1]) for q in paths)
            cn = set()
            for q in paths: cn.update(q[1:])
            any_ = frozenset(poss_an(G, {y}))
            cand = (pdx & any_) - {x}
            t = tag + ("A" if amen else "N")
            acc[t + "_n"] += 1
            if cn != cand:
                acc[t + "_AM_viol"] += 1
                if not (cn <= cand): acc[t + "_cn_notsub"] += 1
                if not (cand <= cn): acc[t + "_cand_notsub"] += 1
                if amen and len(wit) < 8:
                    wit.append((tag, G.tolist(), x, y, sorted(cn), sorted(cand)))


def scan_C(C):
    p = C.shape[0]
    acc = defaultdict(int); wit = []
    for G0 in reachable_mpdags(C):
        check(G0, "mpdag", acc, wit)
        for (a0, b0) in nonadjacent_pairs(G0):
            for (a, b) in ((a0, b0), (b0, a0)):
                H, info = bk_assert(G0, [(a, b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                check(H, "H", acc, wit)
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
            if len(W) < 8: W.extend(w[:8 - len(W)])
    print(json.dumps(dict(sorted(tot.items())), indent=1))
    json.dump({"counts": dict(tot), "wit": W}, open("t_am_p%d.json" % p, "w"), indent=1)


if __name__ == "__main__":
    main()
