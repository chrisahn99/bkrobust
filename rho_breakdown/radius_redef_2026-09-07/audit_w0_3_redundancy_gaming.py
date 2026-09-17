import sys; sys.path.insert(0,'.')
from audit_w0_rk import *

# GAMING VECTOR: add Meek-IMPLIED orientations to K.
# Those leave G0 (and therefore Z, and therefore the set of damaging DAGs)
# bit-identical, but can only raise m(d) for every d.  So r_K is non-decreasing
# under logically redundant restatement of one's own claims.
rng = np.random.default_rng(4242)
probs,_ = gen(rng, 400, n=7, p=0.40, k_claims=3, n_false=0)
raised = same = 0; deltas=[]; rows=[]
for pr in probs:
    cp, g0, K = pr['cpdag'], pr['g0'], pr['K']
    Kset = set(K)
    # every edge Meek DERIVED from K (directed in G0, undirected in Chat, not asserted)
    derived = [(a,b) for (a,b) in g0.directed_edges
               if tuple(sorted((a,b))) in {tuple(sorted(e)) for e in cp.undirected_edges}
               and (a,b) not in Kset]
    if not derived: continue
    K2 = list(K) + derived
    g0b = apply_orientations(cp, K2)
    assert g0b == g0, "graph changed"                      # identical knowledge state
    z2 = optimal_adjustment_set_mpdag(g0b, pr['x'], pr['y'])
    assert frozenset(z2) == pr['z'], "Z changed"           # identical adjustment set
    a = r_K(pr['fibre'], K,  pr['x'], pr['y'], pr['z'])
    b = r_K(pr['fibre'], K2, pr['x'], pr['y'], pr['z'])
    if math.isinf(a): continue
    rows.append((len(K), len(K2), a, b))
    if b > a: raised += 1; deltas.append(b-a)
    else: same += 1
print(f"problems with a Meek-derived orientation available and finite r_K: {len(rows)}")
print(f"  r_K STRICTLY INFLATED by restating implications: {raised}  ({raised/max(1,len(rows)):.3f})")
print(f"  unchanged: {same}")
if deltas: print(f"  mean inflation +{np.mean(deltas):.2f}, max +{max(deltas)}")
print(f"  |K| {np.mean([r[0] for r in rows]):.2f} -> {np.mean([r[1] for r in rows]):.2f}")
print(f"  mean r_K {np.mean([r[2] for r in rows]):.3f} -> {np.mean([r[3] for r in rows]):.3f}")
print(f"  fraction at r_K==1: {np.mean([r[2]==1 for r in rows]):.3f} -> {np.mean([r[3]==1 for r in rows]):.3f}")
print("\n  worked examples (|K|, |K'|, r_K, r_K'):")
for r in sorted(rows, key=lambda t:-(t[3]-t[2]))[:6]: print("   ", r)
