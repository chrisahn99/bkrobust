"""p=6/7 SAMPLE check of the three lemmas that are not yet proved end-to-end:
   PD2 (in H), Claim U, F1', F1."""
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
        D = random_dag(p, deg, rng)
        C = dag_to_cpdag(D)
        ref = v_structures(C)
        G0 = CF.random_mpdag(C, ref, rng)
        pd0 = {}
        NA = nonadjacent_pairs(G0)
        Q = []
        for x in range(p):
            for y in range(p):
                if x != y:
                    z = parts(G0, x, y)
                    if z is not None: Q.append((x, y, z))
        if not NA or not Q: continue
        for (a0, b0) in NA:
            for (a, b) in ((a0, b0), (b0, a0)):
                H, info = bk_assert(G0, [(a, b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                pdH = {}
                for (x, y, (cn0, fb0, pa0, O0, _)) in Q:
                    z = parts(H, x, y)
                    if z is None: continue
                    cn1, fb1, pa1, O1, _ = z
                    acc["n"] += 1
                    if x not in pd0: pd0[x] = frozenset(adjust.poss_de(G0, {x}))
                    if x not in pdH: pdH[x] = frozenset(adjust.poss_de(H, {x}))
                    if not ((pa1 & pdH[x]) <= (cn1 | {x})): acc["PD2_viol"] += 1
                    S = (pa1 & pd0[x]) - {x}
                    if S:
                        acc["U_nonvac"] += 1
                        if not (S <= pdH[x]): acc["U_viol"] += 1
                        if not (S <= cn1): acc["U0_viol"] += 1
                    if O1 & pd0[x]: acc["F1p_viol"] += 1
                    if O1 & fb0: acc["F1_viol"] += 1
    return dict(acc)


def main():
    p = int(sys.argv[1]); ndraw = int(sys.argv[2]); nw = int(sys.argv[3])
    from multiprocessing import Pool
    nb = 200
    jobs = [(p, max(1, ndraw // nb), 900000 + 7919 * i) for i in range(nb)]
    tot = defaultdict(int)
    with Pool(nw) as pool:
        for acc in pool.imap_unordered(blk, jobs):
            for k, v in acc.items(): tot[k] += v
    print("p=%d n_draws=%d" % (p, ndraw), json.dumps(dict(sorted(tot.items()))))


if __name__ == "__main__":
    main()
