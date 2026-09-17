import sys; sys.path.insert(0,'.')
from audit_w0_rk import *
rng = np.random.default_rng(20260907)
t0=time.time()
probs, tried = gen(rng, 300, n=6, p=0.35, k_claims=3, n_false=0)
print(f"generated {len(probs)} problems from {tried} draws in {time.time()-t0:.1f}s")
rk, rv, fib, sp = [], [], [], []
for pr in probs:
    a = r_K(pr['fibre'], pr['K'], pr['x'], pr['y'], pr['z'])
    b, ns = r_val(pr['cpdag'], pr['g0'], pr['x'], pr['y'], pr['z'])
    rk.append(a); rv.append(b); fib.append(len(pr['fibre'])); sp.append(ns)
rk=np.array([np.inf if x==math.inf else x for x in rk],float)
rv=np.array([np.inf if x==math.inf else x for x in rv],float)
n=len(rk)
print(f"\n|[Chat]| mean {np.mean(fib):.1f} max {max(fib)}   |G_Chat| mean {np.mean(sp):.1f}")
print(f"r_K   : ==inf {np.mean(np.isinf(rk)):.3f}  ==1 {np.mean(rk==1):.3f}  ==2 {np.mean(rk==2):.3f}  >=3 {np.mean((rk>=3)&~np.isinf(rk)):.3f}")
print(f"r_val : ==inf {np.mean(np.isinf(rv)):.3f}  ==1 {np.mean(rv==1):.3f}  ==2 {np.mean(rv==2):.3f}  >=3 {np.mean((rv>=3)&~np.isinf(rv)):.3f}")
fin = ~np.isinf(rk)
print(f"\nAmong FINITE r_K (n={fin.sum()}): ==1 {np.mean(rk[fin]==1):.3f}  ==2 {np.mean(rk[fin]==2):.3f}  >=3 {np.mean(rk[fin]>=3):.3f}")
finv = ~np.isinf(rv)
print(f"Among FINITE r_val (n={finv.sum()}): ==1 {np.mean(rv[finv]==1):.3f}  ==2 {np.mean(rv[finv]==2):.3f}")
print(f"\nLemma 1 (r_K >= 1 under P): min r_K = {np.nanmin(rk):.0f}  violations(r_K==0) = {int((rk==0).sum())}")
print("\njoint (r_val,r_K) counts:")
for k,v in sorted(Counter(zip([('inf' if np.isinf(a) else int(a)) for a in rv],
                              [('inf' if np.isinf(a) else int(a)) for a in rk])).items(), key=lambda t:-t[1]):
    print("  r_val=%-4s r_K=%-4s : %d"%(k[0],k[1],v))
# how much entropy does r_K add?
import collections
print("\ndistinct r_K values realised:", sorted(set('inf' if np.isinf(a) else int(a) for a in rk)))
