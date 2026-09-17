import sys
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof")
from explore import *
from collections import defaultdict
def scan_C(C):
    p = C.shape[0]; ref = v_structures(C); acc = defaultdict(int)
    for G0 in reachable_mpdags(C):
        ext = consistent_dag_extensions(G0, ref_vstructs=ref)
        if not ext: continue
        QN = [(x, y) for x in range(p) for y in range(p)
              if x != y and parts(G0, x, y) is None and X.pcp_capped(G0, x, y)]
        if not QN: continue
        for (a0, b0) in nonadjacent_pairs(G0):
            for (a, b) in ((a0, b0), (b0, a0)):
                H, info = bk_assert(G0, [(a, b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                pd0 = {}; pdH = {}
                for (x, y) in QN:
                    z = parts(H, x, y)
                    if z is None: continue
                    O1 = z[3]; acc["S4_report"] += 1
                    if x not in pd0: pd0[x] = frozenset(adjust.poss_de(G0, {x}))
                    if x not in pdH: pdH[x] = frozenset(adjust.poss_de(H, {x}))
                    M = pd0[x] <= pdH[x]
                    if M: acc["S4_report_M"] += 1
                    if O1 & pd0[x]: acc["S4_F1p_viol"] += 1
                    for D in ext:
                        acc["S4_checks"] += 1
                        fbD = parts(D, x, y)[1] if parts(D, x, y) else ({x} | set(adjust.forb(D, x, y)))
                        if O1 & frozenset(adjust.forb(D, x, y)):
                            acc["S4_clause1_viol"] += 1
                            if M: acc["S4_clause1_viol_M"] += 1
                        if not adjust.is_valid_adjustment_set(D, x, y, set(O1)): acc["S4_invalid"] += 1
    return dict(acc)
def _job(args):
    kb, p = args
    return scan_C(np.frombuffer(kb, dtype=np.int8).reshape(p, p).copy())
def main():
    p = int(sys.argv[1]); nw = int(sys.argv[2])
    from multiprocessing import Pool
    Cs = cpdags(p); tot = defaultdict(int)
    with Pool(nw) as pool:
        for acc in pool.imap_unordered(_job, [(C.tobytes(), p) for C in Cs], chunksize=4):
            for k, v in acc.items(): tot[k] += v
    print("p=%d" % p, json.dumps(dict(sorted(tot.items()))))


if __name__ == "__main__":
    main()
