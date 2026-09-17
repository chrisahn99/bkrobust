"""S5 done properly (coherence checked on the UNCLOSED graph), plus the
UNIFIED conjecture: H = MeekClose(G0 + K) for K a mix of non-adjacent pairs and
re-orientations of undirected G0 edges."""
import sys
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof")
from explore import *
from collections import defaultdict
from itertools import combinations


def coh(G, conflict=False):
    return (not conflict) and (not has_directed_cycle(G)) and pdag_extendable(G)


def scan_C(C):
    p = C.shape[0]
    ref = v_structures(C)
    acc = defaultdict(int); wit = defaultdict(list)
    for G0 in reachable_mpdags(C):
        ext = consistent_dag_extensions(G0, ref_vstructs=ref)
        if not ext: continue
        QA = []
        for x in range(p):
            for y in range(p):
                if x == y: continue
                z = parts(G0, x, y)
                if z is not None: QA.append((x, y, z))
        if not QA: continue
        NA = nonadjacent_pairs(G0)
        UE = undirected_edges(G0)
        # candidate single statements: non-adjacent (both dirs) + undirected edges (both dirs)
        stmts = [((a, b), "NA") for (a0, b0) in NA for (a, b) in ((a0, b0), (b0, a0))]
        stmts += [((u, v), "UE") for (u0, v0) in UE for (u, v) in ((u0, v0), (v0, u0))]
        for (s, kind) in stmts:
            a, b = s
            Hn = G0.copy(); Hn[a, b] = 1; Hn[b, a] = 0          # stamp, NO closure
            H, info = bk_assert(G0, [(a, b)])
            for (x, y, (cn0, fb0, pa0, O0, _)) in QA:
                # S5 proper: unclosed graph passes the 3 checks, and is amenable
                if coh(Hn) and not np.array_equal(Hn, H):
                    z5 = parts(Hn, x, y)
                    if z5 is not None:
                        acc["S5_report"] += 1
                        bad = [D for D in ext if not adjust.is_valid_adjustment_set(D, x, y, set(z5[3]))]
                        if bad:
                            acc["S5_INVALID"] += 1
                            if len(wit["S5"]) < 3:
                                wit["S5"].append((C.tolist(), G0.tolist(), a, b, x, y,
                                                  sorted(z5[3]), sorted(O0), es_(H), es_(bad[0])))
                # unified single-statement arm
                if coh(H, info["conflict"]):
                    z = parts(H, x, y)
                    if z is not None:
                        acc[kind + "_report"] += 1
                        if z[3] != O0: acc[kind + "_moved"] += 1
                        if any(not adjust.is_valid_adjustment_set(D, x, y, set(z[3])) for D in ext):
                            acc[kind + "_INVALID"] += 1
                            if len(wit[kind]) < 3:
                                wit[kind].append((C.tolist(), G0.tolist(), a, b, x, y, sorted(z[3])))
        # MIXED pairs (k=2: one NA + one UE)
        for (s1, k1) in stmts:
            if k1 != "NA": continue
            for (s2, k2) in stmts:
                if k2 != "UE": continue
                H, info = bk_assert(G0, [s1, s2])
                if not coh(H, info["conflict"]): continue
                for (x, y, (cn0, fb0, pa0, O0, _)) in QA:
                    z = parts(H, x, y)
                    if z is None: continue
                    acc["MIX_report"] += 1
                    if z[3] != O0: acc["MIX_moved"] += 1
                    if any(not adjust.is_valid_adjustment_set(D, x, y, set(z[3])) for D in ext):
                        acc["MIX_INVALID"] += 1
                        if len(wit["MIX"]) < 3:
                            wit["MIX"].append((C.tolist(), G0.tolist(), s1, s2, x, y, sorted(z[3])))
    return dict(acc), {k: v for k, v in wit.items()}


def es_(M):
    o = []
    for i in range(len(M)):
        for j in range(len(M)):
            if M[i, j] and not M[j, i]: o.append('%d->%d' % (i, j))
            elif M[i, j] and M[j, i] and i < j: o.append('%d-%d' % (i, j))
    return o


def _job(args):
    kb, p = args
    return scan_C(np.frombuffer(kb, dtype=np.int8).reshape(p, p).copy())


def main():
    p = int(sys.argv[1]); nw = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    from multiprocessing import Pool
    Cs = cpdags(p); print("cpdags", len(Cs), flush=True)
    tot = defaultdict(int); W = defaultdict(list)
    with Pool(nw) as pool:
        for acc, w in pool.imap_unordered(_job, [(C.tobytes(), p) for C in Cs], chunksize=4):
            for k, v in acc.items(): tot[k] += v
            for k, v in w.items():
                if len(W[k]) < 3: W[k].extend(v[:3 - len(W[k])])
    print(json.dumps(dict(sorted(tot.items())), indent=1))
    json.dump({"counts": dict(tot), "wit": {k: v for k, v in W.items()}},
              open("sharp2_p%d.json" % p, "w"), indent=1)
    for k in ("S5", "UE", "MIX"):
        for w in W.get(k, [])[:2]: print(k, w[2:])


if __name__ == "__main__":
    main()
