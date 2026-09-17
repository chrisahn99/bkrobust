import sys
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof")
from explore import *

def es(M):
    o = []
    for i in range(len(M)):
        for j in range(len(M)):
            if M[i, j] and not M[j, i]: o.append('%d->%d' % (i, j))
            elif M[i, j] and M[j, i] and i < j: o.append('%d-%d' % (i, j))
    return o

def mk(p, edges):
    G = np.zeros((p, p), dtype=np.int8)
    for e in edges:
        if '->' in e:
            i, j = e.split('->'); G[int(i), int(j)] = 1
        else:
            i, j = e.split('-'); G[int(i), int(j)] = G[int(j), int(i)] = 1
    return G

print("=== S4 witness: G0 NOT amenable (p=4) ===")
G0 = mk(4, ['2-3'])
x, y = 2, 3
print(" G0", es(G0), " amenable(G0,2,3)?", adjust.is_amenable(G0, x, y))
H, info = bk_assert(G0, [(0, 2)])
print(" assert 0->2   H", es(H), " coherent:",
      (not info['conflict'], not has_directed_cycle(H), pdag_extendable(H)))
O1, npaths, amen = X.ostar_and_paths(H, x, y)
print(" amenable(H)?", amen, " O*(H) =", sorted(O1) if O1 is not None else None)
ext = consistent_dag_extensions(G0, ref_vstructs=v_structures(G0))
for D in ext:
    print("   D =", es(D), " valid(O*_H) =", adjust.is_valid_adjustment_set(D, x, y, set(O1)))

print()
print("=== S2 witness: acyclicity/extendability check dropped (p=4) ===")
G0 = mk(4, ['1->3', '2->0', '2->1', '3->0'])
x, y = 1, 0
print(" G0", es(G0), " amen:", adjust.is_amenable(G0, x, y), " O*(G0)=", sorted(adjust.optimal_adjustment_set(G0, x, y)))
H, info = bk_assert(G0, [(3, 2)])
print(" assert 3->2   H", es(H), " conflict", info['conflict'],
      " cycle", has_directed_cycle(H), " extendable", pdag_extendable(H))
O1, _, amen = X.ostar_and_paths(H, x, y)
print(" amenable(H)?", amen, " O*(H) =", sorted(O1))
C = mk(4, ['1-2', '1-3', '2->0', '3->0'])
ext = consistent_dag_extensions(G0, ref_vstructs=v_structures(C))
for D in ext:
    print("   D =", es(D), " valid =", adjust.is_valid_adjustment_set(D, x, y, set(O1)))

print()
print("=== S5 witness: Meek closure omitted (p=4) ===")
G0 = mk(4, ['1->2', '1-3', '2->0', '3->0'])
x, y = 2, 0
C = mk(4, ['1-2', '1-3', '2->0', '3->0'])
print(" G0", es(G0), " O*(G0)=", sorted(adjust.optimal_adjustment_set(G0, x, y)))
Hn = G0.copy(); Hn[0, 1] = 1; Hn[1, 0] = 0
Hc, _ = bk_assert(G0, [(0, 1)])
print(" stamp 0->1 WITHOUT closure:", es(Hn), "   WITH closure:", es(Hc))
for nm, G in (("no-closure", Hn), ("closed", Hc)):
    O1, _, amen = X.ostar_and_paths(G, x, y)
    if not amen: print("  ", nm, "not amenable"); continue
    ext = consistent_dag_extensions(G0, ref_vstructs=v_structures(C))
    bad = [D for D in ext if not adjust.is_valid_adjustment_set(D, x, y, set(O1))]
    print("  ", nm, "O* =", sorted(O1), " invalid in", len(bad), "of", len(ext), "D",
          ([es(D) for D in bad[:1]] if bad else ""))
