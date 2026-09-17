"""SC-refined:  Z cap forb_D = {}  and  Z superset (O*(D) cap Reach_D)  =>  Z valid.
Reach_D = nodes reachable from X in skeleton(D_pbd) restricted to (V \\ forb_D) u {X}.
Exhaustive over all DAGs on p nodes, all queries, all such Z."""
import sys
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof")
from explore import *
from itertools import combinations

p = int(sys.argv[1])
n_t = n_bad = 0; bad = []
for D in all_dags(p):
    for x in range(p):
        for y in range(p):
            if x == y: continue
            z = parts(D, x, y)
            if z is None: continue
            cn, fb, pa, O, _ = z
            Dp = D.copy()
            for w in cn:
                if is_directed(Dp, x, w): Dp[x, w] = 0
            S = ((Dp + Dp.T) > 0)
            allowed = set(v for v in range(p) if v not in fb) | {x}
            seen = {x}; st = [x]
            while st:
                v = st.pop()
                for w in np.flatnonzero(S[v]):
                    w = int(w)
                    if w in allowed and w not in seen: seen.add(w); st.append(w)
            need = O & frozenset(seen)
            free = [v for v in range(p) if v not in fb and v not in need]
            for r in range(len(free) + 1):
                for extra in combinations(free, r):
                    Z = set(need) | set(extra)
                    n_t += 1
                    if not adjust.is_valid_adjustment_set(D, x, y, Z):
                        n_bad += 1
                        if len(bad) < 4: bad.append((D.tolist(), x, y, sorted(Z), sorted(O), sorted(need)))
print("SC-refined p=%d: tested=%d violations=%d" % (p, n_t, n_bad))
for b in bad: print(b)
