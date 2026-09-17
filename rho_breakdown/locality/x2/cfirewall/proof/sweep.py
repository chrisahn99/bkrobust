"""Census sweep recording relations among G0 / H / D objects."""
import sys
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof")
from explore import *
from collections import defaultdict


def scan_C(C):
    p = C.shape[0]
    ref = v_structures(C)
    acc = defaultdict(int)
    wit = defaultdict(list)
    for G0 in reachable_mpdags(C):
        NA = nonadjacent_pairs(G0)
        if not NA: continue
        Q = []
        for x in range(p):
            for y in range(p):
                if x == y: continue
                z = parts(G0, x, y)
                if z is not None: Q.append((x, y, z))
        if not Q: continue
        ext = consistent_dag_extensions(G0, ref_vstructs=ref)
        if not ext: continue
        extp = {}
        for (a0, b0) in NA:
            for (a, b) in ((a0, b0), (b0, a0)):
                H, info = bk_assert(G0, [(a, b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H):
                    continue
                for (x, y, (cn0, fb0, pa0, O0, _)) in Q:
                    z = parts(H, x, y)
                    if z is None:
                        acc["abort"] += 1
                        continue
                    cn1, fb1, pa1, O1, _ = z
                    acc["n"] += 1
                    if O1 & fb0: acc["F1_viol"] += 1
                    for di, D in enumerate(ext):
                        k = (di, x, y)
                        if k not in extp: extp[k] = parts(D, x, y)
                        zd = extp[k]
                        assert zd is not None
                        cnD, fbD, paD, OD, _ = zd
                        acc["nD"] += 1
                        # P2: O*(D) subset O*(H) ?
                        if not (OD <= O1): acc["P2_viol"] += 1
                        if not (O0 <= OD): acc["O0_subset_OD_viol"] += 1
                        if not (O0 <= O1): acc["O0_subset_O1_viol"] += 1
                        if not (cnD <= cn1): acc["cnD_sub_cn1_viol"] += 1
                        if cnD != cn0: acc["cnD_ne_cn0"] += 1
                        if fbD != fb0: acc["fbD_ne_fb0"] += 1
                        if not (O1 & fbD) : pass
                        else: acc["O1_meets_fbD"] += 1
                        if not adjust.is_valid_adjustment_set(D, x, y, set(O1)):
                            acc["INVALID"] += 1
                            if len(wit["inv"]) < 3:
                                wit["inv"].append((C.tolist(), G0.tolist(), a, b, x, y, sorted(O1), D.tolist()))
                        # is O1 union OD valid?
                        if not (OD <= O1) and len(wit["p2"]) < 6:
                            wit["p2"].append((G0.tolist(), H.tolist(), a, b, x, y,
                                              sorted(O1), sorted(OD), sorted(O0), D.tolist()))
    return dict(acc), {k: v for k, v in wit.items()}


def main():
    p = int(sys.argv[1])
    nw = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    from multiprocessing import Pool
    Cs = cpdags(p)
    print("cpdags", len(Cs), flush=True)
    tot = defaultdict(int); W = defaultdict(list)
    with Pool(nw) as pool:
        for acc, w in pool.imap_unordered(_job, [(C.tobytes(), p) for C in Cs], chunksize=4):
            for k, v in acc.items(): tot[k] += v
            for k, v in w.items():
                if len(W[k]) < 6: W[k].extend(v[:6 - len(W[k])])
    print(json.dumps(dict(sorted(tot.items())), indent=1))
    json.dump({"counts": dict(tot), "wit": {k: v for k, v in W.items()}},
              open("sweep_p%d.json" % p, "w"), indent=1)


def _job(args):
    kb, p = args
    return scan_C(np.frombuffer(kb, dtype=np.int8).reshape(p, p).copy())


if __name__ == "__main__":
    main()
