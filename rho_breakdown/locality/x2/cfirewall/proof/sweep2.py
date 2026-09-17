"""Sweep 2: test refined sufficient conditions."""
import sys
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof")
from explore import *
from collections import defaultdict


def pbd(D, x, cnD):
    Dp = D.copy()
    for w in cnD:
        if is_directed(Dp, x, w): Dp[x, w] = 0
    return Dp


def reach_outside_forb(Dp, x, fbD, p):
    """Nodes reachable from x in skeleton(Dp) restricted to (V \\ fbD) u {x}."""
    S = ((Dp + Dp.T) > 0)
    allowed = set(v for v in range(p) if v not in fbD) | {x}
    seen = {x}; st = [x]
    while st:
        v = st.pop()
        for w in np.flatnonzero(S[v]):
            w = int(w)
            if w in allowed and w not in seen:
                seen.add(w); st.append(w)
    return seen


def scan_C(C):
    p = C.shape[0]
    ref = v_structures(C)
    acc = defaultdict(int); wit = defaultdict(list)
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
        cache = {}
        for (a0, b0) in NA:
            for (a, b) in ((a0, b0), (b0, a0)):
                H, info = bk_assert(G0, [(a, b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H):
                    continue
                for (x, y, (cn0, fb0, pa0, O0, _)) in Q:
                    z = parts(H, x, y)
                    if z is None: continue
                    cn1, fb1, pa1, O1, _ = z
                    acc["n"] += 1
                    for di, D in enumerate(ext):
                        k = (di, x, y)
                        if k not in cache:
                            cnD, fbD, paD, OD, _ = parts(D, x, y)
                            Dp = pbd(D, x, cnD)
                            R = reach_outside_forb(Dp, x, fbD, p)
                            cache[k] = (cnD, fbD, paD, OD, R)
                        cnD, fbD, paD, OD, R = cache[k]
                        acc["nD"] += 1
                        need = OD & R
                        if not (need <= O1):
                            acc["HypA_viol"] += 1
                            if len(wit["A"]) < 6:
                                wit["A"].append((G0.tolist(), H.tolist(), D.tolist(), a, b, x, y,
                                                 sorted(O1), sorted(OD), sorted(need)))
                        # cnD vs cn1
                        if not (cnD <= cn1): acc["cnD_notsub_cn1"] += 1
                        if not (cn1 <= fb0): acc["cn1_notsub_fb0"] += 1
                        if not (cn0 <= fb1): acc["cn0_notsub_fb1"] += 1
                        if fb1 != fb0: acc["fb1_ne_fb0"] += 1
                        if not (fb0 <= fb1): acc["fb0_notsub_fb1"] += 1
                        if not (fb1 <= fb0): acc["fb1_notsub_fb0"] += 1
    return dict(acc), {k: v for k, v in wit.items()}


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
                if len(W[k]) < 6: W[k].extend(v[:6 - len(W[k])])
    print(json.dumps(dict(sorted(tot.items())), indent=1))
    json.dump({"counts": dict(tot), "wit": {k: v for k, v in W.items()}},
              open("sweep2_p%d.json" % p, "w"), indent=1)


if __name__ == "__main__":
    main()
