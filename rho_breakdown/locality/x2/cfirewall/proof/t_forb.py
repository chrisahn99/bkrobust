"""Census: is forb(X,Y,G0) subset forb(X,Y,H)?  And the F1 invariant."""
import sys
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof")
from explore import *
from collections import defaultdict


def scan_C(C):
    p = C.shape[0]
    acc = defaultdict(int); wit = []
    for G0 in reachable_mpdags(C):
        NA = nonadjacent_pairs(G0)
        if not NA: continue
        Q = []
        for x in range(p):
            for y in range(p):
                if x != y:
                    z = parts(G0, x, y)
                    if z is not None: Q.append((x, y, z))
        if not Q: continue
        for (a0, b0) in NA:
            for (a, b) in ((a0, b0), (b0, a0)):
                H, info = bk_assert(G0, [(a, b)])
                coh = not (info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H))
                for (x, y, (cn0, fb0, pa0, O0, _)) in Q:
                    z = parts(H, x, y)
                    tag = "coh" if coh else "ctl"
                    if z is None:
                        acc[tag + "_abort"] += 1; continue
                    cn1, fb1, pa1, O1, _ = z
                    acc[tag + "_n"] += 1
                    if not (fb0 <= fb1):
                        acc[tag + "_FORBGROW_viol"] += 1
                        if coh and len(wit) < 6:
                            wit.append((G0.tolist(), H.tolist(), a, b, x, y,
                                        sorted(fb0), sorted(fb1)))
                    if fb1 != fb0: acc[tag + "_fb_ne"] += 1
                    if O1 & fb0: acc[tag + "_F1_viol"] += 1
                    if not (cn0 <= fb1): acc[tag + "_cn0_notsub_fb1"] += 1
                    if not (cn1 <= fb0): acc[tag + "_cn1_notsub_fb0"] += 1
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
    json.dump({"counts": dict(tot), "wit": W}, open("t_forb_p%d.json" % p, "w"), indent=1)


if __name__ == "__main__":
    main()
