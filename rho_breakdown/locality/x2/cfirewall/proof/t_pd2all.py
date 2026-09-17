"""PD2 over ALL Meek-closed acyclic PDAGs on p nodes (not just reachable MPDAGs)."""
import sys
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof")
from explore import *
from itertools import product
from collections import defaultdict

p = int(sys.argv[1])
nw = int(sys.argv[2]) if len(sys.argv) > 2 else 10
idx = [(i, j) for i in range(p) for j in range(i + 1, p)]


def blk(args):
    lo, hi = args
    acc = defaultdict(int); wit = []
    for n in range(lo, hi):
        st = []
        t = n
        for _ in range(len(idx)):
            st.append(t & 3); t >>= 2
        G = np.zeros((p, p), dtype=np.int8)
        for (i, j), s in zip(idx, st):
            if s == 1: G[i, j] = 1
            elif s == 2: G[j, i] = 1
            elif s == 3: G[i, j] = G[j, i] = 1
        if has_directed_cycle(G): continue
        if not np.array_equal(meek_closure(G), G): continue
        acc["graphs"] += 1
        for x in range(p):
            pdx = frozenset(adjust.poss_de(G, {x}))
            for y in range(p):
                if x == y: continue
                paths = X.pcp_capped(G, x, y)
                if not paths: continue
                if not all(is_directed(G, q[0], q[1]) for q in paths): continue
                cn = set()
                for q in paths: cn.update(q[1:])
                pa = adjust.parents_of_set(G, cn)
                acc["amen_q"] += 1
                bad = (pa & pdx) - cn - {x}
                if bad:
                    acc["PD2_viol"] += 1
                    if len(wit) < 4: wit.append((G.tolist(), x, y, sorted(bad), sorted(cn)))
                # also the back-edge lemma NB
                for v in pdx:
                    if v != x and is_directed(G, v, x):
                        acc["NB_viol"] += 1
    return dict(acc), wit


def main():
    from multiprocessing import Pool
    N = 4 ** len(idx)
    nb = 400
    jobs = [(i * N // nb, (i + 1) * N // nb) for i in range(nb)]
    tot = defaultdict(int); W = []
    with Pool(nw) as pool:
        for acc, w in pool.imap_unordered(blk, jobs):
            for k, v in acc.items(): tot[k] += v
            if len(W) < 4: W.extend(w[:4 - len(W)])
    print("p=%d ALL Meek-closed acyclic PDAGs:" % p, json.dumps(dict(sorted(tot.items()))))
    print(W[:2])


if __name__ == "__main__":
    main()
