"""TEST SC: for a DAG D and any Z with O*(D) subset Z and Z cap forb_D = empty,
is Z a valid adjustment set for (X,Y)?  Exhaustive over all DAGs on p nodes,
all ordered queries with a causal path, all such Z."""
import sys
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof")
from explore import *
from itertools import combinations

p = int(sys.argv[1]) if len(sys.argv) > 1 else 4
n_tested = n_bad = 0
bad = []
for D in all_dags(p):
    for x in range(p):
        for y in range(p):
            if x == y: continue
            z = parts(D, x, y)
            if z is None: continue
            cn, fb, pa, O, _ = z
            free = [v for v in range(p) if v not in fb and v not in O]
            for r in range(len(free) + 1):
                for extra in combinations(free, r):
                    Z = set(O) | set(extra)
                    n_tested += 1
                    if not adjust.is_valid_adjustment_set(D, x, y, Z):
                        n_bad += 1
                        if len(bad) < 5:
                            bad.append((D.tolist(), x, y, sorted(Z), sorted(O)))
print("SC p=%d: tested=%d violations=%d" % (p, n_tested, n_bad))
for b in bad: print(b)
