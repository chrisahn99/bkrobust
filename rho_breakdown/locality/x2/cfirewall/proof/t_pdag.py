"""Is PD2-under-amenability a fact about ARBITRARY acyclic PDAGs, or does it
need Meek-closedness?  Exhaustive over all PDAGs on p nodes."""
import sys
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof")
from explore import *
from itertools import product
from collections import defaultdict

p = int(sys.argv[1]) if len(sys.argv) > 1 else 4
idx = [(i, j) for i in range(p) for j in range(i + 1, p)]
acc = defaultdict(int); wit = []
for state in product([0, 1, 2, 3], repeat=len(idx)):
    G = np.zeros((p, p), dtype=np.int8)
    for (i, j), s in zip(idx, state):
        if s == 1: G[i, j] = 1
        elif s == 2: G[j, i] = 1
        elif s == 3: G[i, j] = G[j, i] = 1
    if has_directed_cycle(G): continue
    mc = np.array_equal(meek_closure(G), G)
    for x in range(p):
        pdx = frozenset(adjust.poss_de(G, {x}))
        for y in range(p):
            if x == y: continue
            paths = X.pcp_capped(G, x, y)
            if not paths: continue
            amen = all(is_directed(G, q[0], q[1]) for q in paths)
            if not amen: continue
            cn = set()
            for q in paths: cn.update(q[1:])
            pa = adjust.parents_of_set(G, cn)
            bad = (pa & pdx) - cn - {x}
            tag = "closed" if mc else "open"
            acc[tag + "_n"] += 1
            if bad:
                acc[tag + "_viol"] += 1
                if len(wit) < 6: wit.append((tag, G.tolist(), x, y, sorted(bad), sorted(cn)))
print(json.dumps(dict(sorted(acc.items())), indent=1))
def es(M):
    M = np.array(M); o = []
    for i in range(len(M)):
        for j in range(len(M)):
            if M[i, j] and not M[j, i]: o.append('%d->%d' % (i, j))
            elif M[i, j] and M[j, i] and i < j: o.append('%d-%d' % (i, j))
    return o
for w in wit: print(w[0], 'q=(%d,%d)' % (w[2], w[3]), 'bad', w[4], 'cn', w[5], es(w[1]))
