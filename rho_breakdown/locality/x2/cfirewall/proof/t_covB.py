"""Coverage of the fully-proved clause-(2) condition SC-refined:
   O*(H) must contain  (pa_D(cn_D) \\ forb_D) cap Reach_D .
Sampled at p=5,6."""
import sys
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof")
from explore import *
from collections import defaultdict
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/code")
import cfirewall as CF


def blk(args):
    p, n, seed = args
    rng = np.random.default_rng(seed)
    acc = defaultdict(int)
    for _ in range(n):
        deg = float(rng.choice([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0]))
        D0 = random_dag(p, deg, rng); C = dag_to_cpdag(D0); ref = v_structures(C)
        G0 = CF.random_mpdag(C, ref, rng)
        if len(undirected_edges(G0)) > 12: continue
        ext = consistent_dag_extensions(G0, ref_vstructs=ref)
        if not ext: continue
        Q = []
        for x in range(p):
            for y in range(p):
                if x != y:
                    z = parts(G0, x, y)
                    if z is not None: Q.append((x, y, z))
        if not Q: continue
        pd0 = {}; cache = {}
        for (a0, b0) in nonadjacent_pairs(G0):
            for (a, b) in ((a0, b0), (b0, a0)):
                H, info = bk_assert(G0, [(a, b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                pdH = {}
                for (x, y, (cn0, fb0, pa0, O0, _)) in Q:
                    z = parts(H, x, y)
                    if z is None: continue
                    O1 = z[3]
                    if x not in pd0: pd0[x] = frozenset(adjust.poss_de(G0, {x}))
                    if x not in pdH: pdH[x] = frozenset(adjust.poss_de(H, {x}))
                    M = pd0[x] <= pdH[x]
                    moved = (O1 != O0)
                    for di, D in enumerate(ext):
                        k = (di, x, y)
                        if k not in cache:
                            cnD, fbD, paD, OD, _ = parts(D, x, y)
                            Dp = D.copy()
                            for wv in cnD:
                                if is_directed(Dp, x, wv): Dp[x, wv] = 0
                            S = ((Dp + Dp.T) > 0)
                            allowed = set(vv for vv in range(p) if vv not in fbD) | {x}
                            seen = {x}; st = [x]
                            while st:
                                vv = st.pop()
                                for ww in np.flatnonzero(S[vv]):
                                    ww = int(ww)
                                    if ww in allowed and ww not in seen:
                                        seen.add(ww); st.append(ww)
                            cache[k] = (fbD, OD, frozenset(seen))
                        fbD, OD, R = cache[k]
                        acc["nD"] += 1
                        if moved: acc["nD_moved"] += 1
                        c1 = not (O1 & fbD)
                        c2full = OD <= O1
                        c2ref = (OD & R) <= O1
                        if c1 and c2ref: acc["proved_ref"] += 1
                        if c1 and c2full: acc["proved_full"] += 1
                        if M and c1 and c2ref: acc["proved_ref_M"] += 1
                        if moved and c1 and c2ref: acc["proved_ref_moved"] += 1
    return dict(acc)


def main():
    p = int(sys.argv[1]); nd = int(sys.argv[2]); nw = int(sys.argv[3])
    from multiprocessing import Pool
    nb = 200
    jobs = [(p, max(1, nd // nb), 31337 + 7919 * i) for i in range(nb)]
    tot = defaultdict(int)
    with Pool(nw) as pool:
        for acc in pool.imap_unordered(blk, jobs):
            for k, v in acc.items(): tot[k] += v
    print("p=%d" % p, json.dumps(dict(sorted(tot.items()))))
    if tot["nD"]:
        print("  SC-refined coverage: %.4f%%   SC-full: %.4f%%   moved-only: %.4f%%" %
              (100 * tot["proved_ref"] / tot["nD"], 100 * tot["proved_full"] / tot["nD"],
               100 * tot["proved_ref_moved"] / max(tot["nD_moved"], 1)))


if __name__ == "__main__":
    main()
