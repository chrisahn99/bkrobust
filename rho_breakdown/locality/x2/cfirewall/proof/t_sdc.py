"""Lemma SDC:  G Meek-closed, acyclic, amenable rel (X,Y)
  =>  there is NO directed edge v->w and node c in cn_G with  w ~> c ~> v  (p.d.).
(i.e. no semi-directed cycle through cn_G)."""
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
    if not np.array_equal(meek_closure(G), G): continue
    pd = {u: frozenset(adjust.poss_de(G, {u})) for u in range(p)}
    de_edges = [(u, w) for u in range(p) for w in range(p) if is_directed(G, u, w)]
    for x in range(p):
        for y in range(p):
            if x == y: continue
            paths = X.pcp_capped(G, x, y)
            if not paths: continue
            amen = all(is_directed(G, q[0], q[1]) for q in paths)
            cn = set()
            for q in paths: cn.update(q[1:])
            hit = False
            for (v, w) in de_edges:
                for c in cn:
                    if c in pd[w] and v in pd[c]:
                        hit = True; break
                if hit: break
            tag = "A" if amen else "N"
            acc[tag + "_n"] += 1
            if hit:
                acc[tag + "_SDC"] += 1
                if amen and len(wit) < 6:
                    wit.append((G.tolist(), x, y, sorted(cn)))
print("p=%d" % p, json.dumps(dict(sorted(acc.items()))))
def es(M):
    M = np.array(M); o = []
    for i in range(len(M)):
        for j in range(len(M)):
            if M[i, j] and not M[j, i]: o.append('%d->%d' % (i, j))
            elif M[i, j] and M[j, i] and i < j: o.append('%d-%d' % (i, j))
    return o
for w in wit: print(' q=(%d,%d) cn=%s' % (w[1], w[2], w[3]), es(w[0]))
