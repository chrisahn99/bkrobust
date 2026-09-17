import sys; sys.path.insert(0,'.')
from audit_w0_rk import *
def show(tag, probs):
    rk, rv = [], []
    for pr in probs:
        rk.append(r_K(pr['fibre'], pr['K'], pr['x'], pr['y'], pr['z']))
        b,_ = r_val(pr['cpdag'], pr['g0'], pr['x'], pr['y'], pr['z']); rv.append(b)
    rk=np.array([np.inf if x==math.inf else x for x in rk],float)
    rv=np.array([np.inf if x==math.inf else x for x in rv],float)
    fin=~np.isinf(rk); finv=~np.isinf(rv)
    print(f"\n### {tag}  n={len(rk)}  |fibre| mean {np.mean([len(p['fibre']) for p in probs]):.1f}")
    print(f"  r_K   inf {np.mean(np.isinf(rk)):.3f} | of FINITE: ==1 {np.mean(rk[fin]==1):.3f} ==2 {np.mean(rk[fin]==2):.3f} >=3 {np.mean(rk[fin]>=3):.3f}  (n_fin={fin.sum()})")
    print(f"  r_val inf {np.mean(np.isinf(rv)):.3f} | of FINITE: ==1 {np.mean(rv[finv]==1):.3f} ==2 {np.mean(rv[finv]==2):.3f} >=3 {np.mean(rv[finv]>=3):.3f}  (n_fin={finv.sum()})")
    agree = np.mean([(np.isinf(a) and np.isinf(b)) or a==b for a,b in zip(rk,rv)])
    print(f"  r_K == r_val on {agree:.4f} of problems")
    d = [(('inf' if np.isinf(a) else int(a)), ('inf' if np.isinf(b) else int(b))) for a,b in zip(rv,rk) if not ((np.isinf(a) and np.isinf(b)) or a==b)]
    print(f"  disagreements ({len(d)}): {Counter(d).most_common()}")
    return rk, rv
for (n,p,k) in [(6,0.35,3),(7,0.40,4),(8,0.35,4),(7,0.30,2)]:
    rng = np.random.default_rng(11+n*10+int(p*100))
    probs,_ = gen(rng, 250, n=n, p=p, k_claims=k, n_false=0)
    show(f"n={n} p={p} |K|={k}", probs)
