"""Sharpness of the CLAUSE-(1) theorem: is Meek-closure of H needed for
   O*(H) cap forb_D = empty ?   (stamp a->b, do NOT close)"""
import sys
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof")
from explore import *
from collections import defaultdict


def scan_C(C):
    p = C.shape[0]; ref = v_structures(C)
    acc = defaultdict(int); wit = []
    for G0 in reachable_mpdags(C):
        ext = consistent_dag_extensions(G0, ref_vstructs=ref)
        if not ext: continue
        Q = []
        for x in range(p):
            for y in range(p):
                if x != y:
                    z = parts(G0, x, y)
                    if z is not None: Q.append((x, y, z))
        if not Q: continue
        for (a0, b0) in nonadjacent_pairs(G0):
            for (a, b) in ((a0, b0), (b0, a0)):
                Hn = G0.copy(); Hn[a, b] = 1; Hn[b, a] = 0
                H, info = bk_assert(G0, [(a, b)])
                if np.array_equal(Hn, H): continue
                closed_ok = not (info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H))
                nc_ok = (not has_directed_cycle(Hn)) and pdag_extendable(Hn)
                for (x, y, (cn0, fb0, pa0, O0, _)) in Q:
                    z = parts(Hn, x, y)
                    if z is None: continue
                    O1 = z[3]
                    acc["nc_report"] += 1
                    if nc_ok: acc["nc_report_coh"] += 1
                    if O1 & fb0:
                        acc["nc_F1_viol"] += 1
                        if nc_ok: acc["nc_F1_viol_coh"] += 1
                    for D in ext:
                        fbD = parts(D, x, y)[1]
                        if O1 & fbD:
                            acc["nc_clause1_viol"] += 1
                            if nc_ok:
                                acc["nc_clause1_viol_coh"] += 1
                                if len(wit) < 3:
                                    wit.append((G0.tolist(), a, b, x, y, sorted(O1), sorted(fbD), D.tolist()))
                            break
    return dict(acc), wit


def _job(args):
    kb, p = args
    return scan_C(np.frombuffer(kb, dtype=np.int8).reshape(p, p).copy())


def main():
    p = int(sys.argv[1]); nw = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    from multiprocessing import Pool
    Cs = cpdags(p); tot = defaultdict(int); W = []
    with Pool(nw) as pool:
        for acc, w in pool.imap_unordered(_job, [(C.tobytes(), p) for C in Cs], chunksize=4):
            for k, v in acc.items(): tot[k] += v
            if len(W) < 3: W.extend(w[:3 - len(W)])
    print("p=%d" % p, json.dumps(dict(sorted(tot.items())), indent=1))
    def es(M):
        M = np.array(M); o = []
        for i in range(len(M)):
            for j in range(len(M)):
                if M[i, j] and not M[j, i]: o.append('%d->%d' % (i, j))
                elif M[i, j] and M[j, i] and i < j: o.append('%d-%d' % (i, j))
        return o
    for w in W[:2]:
        print(" G0", es(w[0]), "assert %d->%d q=(%d,%d)" % (w[1], w[2], w[3], w[4]),
              "O*=", w[5], "forb_D=", w[6], "D=", es(w[7]))


if __name__ == "__main__":
    main()
