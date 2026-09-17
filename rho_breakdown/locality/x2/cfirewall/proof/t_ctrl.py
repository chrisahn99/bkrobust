"""In the CONTROL arm (coherence check failed, edge stamped anyway), which CLAUSE
of the adjustment criterion does the reported O*(H) break?"""
import sys
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof")
from explore import *
from collections import defaultdict


def scan_C(C):
    p = C.shape[0]; ref = v_structures(C)
    acc = defaultdict(int)
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
                H, info = bk_assert(G0, [(a, b)])
                if not (info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H)):
                    continue
                for (x, y, (cn0, fb0, pa0, O0, _)) in Q:
                    z = parts(H, x, y)
                    if z is None: continue
                    O1 = z[3]
                    acc["ctrl_report"] += 1
                    if O1 & fb0: acc["ctrl_F1_viol"] += 1
                    for D in ext:
                        cnD, fbD, paD, OD, _ = parts(D, x, y)
                        ok = adjust.is_valid_adjustment_set(D, x, y, set(O1))
                        acc["ctrl_checks"] += 1
                        if not ok:
                            acc["ctrl_invalid"] += 1
                            if O1 & fbD: acc["ctrl_invalid_clause1"] += 1
                            else: acc["ctrl_invalid_clause2_only"] += 1
    return dict(acc)


def _job(args):
    kb, p = args
    return scan_C(np.frombuffer(kb, dtype=np.int8).reshape(p, p).copy())


def main():
    p = int(sys.argv[1]); nw = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    from multiprocessing import Pool
    Cs = cpdags(p); tot = defaultdict(int)
    with Pool(nw) as pool:
        for acc in pool.imap_unordered(_job, [(C.tobytes(), p) for C in Cs], chunksize=4):
            for k, v in acc.items(): tot[k] += v
    print("p=%d" % p, json.dumps(dict(sorted(tot.items())), indent=1))


if __name__ == "__main__":
    main()
