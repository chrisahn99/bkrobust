"""Single-graph claim PD2:  pa_G(cn_G) cap poss_de_G(X)  subset  cn_G u {X},
for every reachable MPDAG G (and for every H = MeekClose(G0+a->b) too)."""
import sys
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof")
from explore import *
from collections import defaultdict


def check(G, tag, acc, wit, require_amen=False):
    p = G.shape[0]
    for x in range(p):
        pdx = frozenset(adjust.poss_de(G, {x}))
        for y in range(p):
            if x == y: continue
            paths = X.pcp_capped(G, x, y)
            if not paths: continue
            amen = all(is_directed(G, q[0], q[1]) for q in paths)
            if require_amen and not amen: continue
            cn = set()
            for q in paths: cn.update(q[1:])
            pa = adjust.parents_of_set(G, cn)
            acc[tag + "_n"] += 1
            bad = (pa & pdx) - cn - {x}
            if bad:
                acc[tag + "_PD2_viol"] += 1
                if len(wit) < 6:
                    wit.append((tag, G.tolist(), x, y, sorted(bad), sorted(cn), amen))


def scan_C(C):
    p = C.shape[0]
    acc = defaultdict(int); wit = []
    for G0 in reachable_mpdags(C):
        check(G0, "mpdag", acc, wit)
        check(G0, "mpdagA", acc, wit, True)
        for (a0, b0) in nonadjacent_pairs(G0):
            for (a, b) in ((a0, b0), (b0, a0)):
                H, info = bk_assert(G0, [(a, b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                check(H, "H", acc, wit)
                check(H, "HA", acc, wit, True)
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
    json.dump({"counts": dict(tot), "wit": W}, open("t_pd2_p%d.json" % p, "w"), indent=1)


if __name__ == "__main__":
    main()
