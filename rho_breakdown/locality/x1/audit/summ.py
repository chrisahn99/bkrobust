import json, sys, numpy as np
from collections import Counter
res = json.load(open(sys.argv[1])); lab = sys.argv[2]
n = len(res)
def W(k, N):
    if N == 0: return (float('nan'),)*2
    z=1.959963985; ph=k/N; d=1+z*z/N
    c=(ph+z*z/(2*N))/d; h=z*np.sqrt(ph*(1-ph)/N + z*z/(4*N*N))/d
    return (c-h, c+h)
print(f"=== {lab}  n_scm_analysed = {n} ===")
nN = np.array([r['n_N'] for r in res]); nU=np.array([r['n_U'] for r in res])
nA = np.array([r['n_A'] for r in res]); P=np.array([r['p'] for r in res])
print(f"|N(C)| mean {nN.mean():.2f} median {np.median(nN):.0f} min {nN.min()} max {nN.max()}")
print(f"|U(C)| mean {nU.mean():.2f} ; |A(C)| mean {nA.mean():.2f}")
for thr in (0,1,2,3,4):
    k=int((nN<=thr).sum()); print(f"  frac |N(C)| <= {thr}: {k}/{n} = {k/n:.4f}  Wilson {W(k,n)[0]:.4f},{W(k,n)[1]:.4f}")
print(f"  frac |N(C)|>=4 (S-uni eligible, unordered-w/o-repl): {(nN>=4).sum()}/{n} = {(nN>=4).mean():.4f}")
print(f"  frac |N(C)|>=2 (S-uni eligible, ordered-w/o-repl)  : {(nN>=2).sum()}/{n} = {(nN>=2).mean():.4f}")
xy=np.array([r['xy_nonadj'] for r in res])
print(f"  frac X,Y NON-adjacent in C: {xy.sum()}/{n} = {xy.mean():.4f}")
# probability a uniform S-uni draw is the query pair
pq = np.array([ (2.0/(2*r['n_N'])) if (r['xy_nonadj'] and r['n_N']>0) else 0.0 for r in res])
print(f"  E[P(one S-uni draw == the query pair {{X,Y}})] = {pq.mean():.4f}")
print(f"  E[expected # of 4 S-uni slots that are the query pair] = {4*pq.mean():.4f}")
print(f"  frac SCMs where >=1 of 4 slots hits the query pair (approx 1-(1-q)^4) = {np.mean(1-(1-pq)**4):.4f}")
# hop composition
allN = Counter(); allK=Counter()
for r in res:
    allN.update(r['hopsN']); allK.update(r['hopsK'])
tN=sum(allN.values()); tK=sum(allK.values())
print(f"  hop distribution of N(C) pairs (pool, total {tN}): " +
      ", ".join(f"{h}:{c/tN:.4f}" for h,c in sorted(allN.items())))
print(f"  hop distribution of K stmts  (pool, total {tK}): " +
      ", ".join(f"{h}:{c/tK:.4f}" for h,c in sorted(allK.items())))
sh1 = sum(c for h,c in allN.items() if h>=1)/tN
print(f"  share of N(C) pool at hop>=1 : {sh1:.4f}  -> expected hop>=1 stmts at rho=1 = {4*n*sh1:.0f}")
shK1 = sum(c for h,c in allK.items() if h>=1)/tK
print(f"  share of K   pool at hop>=1 : {shK1:.4f}  -> arm R hop>=1 stmts at rho=1 = {4*n*shK1:.0f}")
# cycle lower bound
oc=np.array([r['n_ord_cyc'] for r in res]); od=np.array([r['n_ord'] for r in res])
m=od>0
print(f"  LOWER BOUND on P(added a->b closes a directed cycle already present in C): "
      f"{oc[m].sum()}/{od[m].sum()} = {oc[m].sum()/od[m].sum():.4f}")
# vstruct destruction
kv=np.array([r['n_N_kills_vstruct'] for r in res]); nv=np.array([r['n_vstruct_C'] for r in res])
print(f"  P(uniform N(C) pair is the endpoint pair of an existing unshielded collider of C) "
      f"= {kv.sum()}/{nN.sum()} = {kv.sum()/nN.sum():.4f}   (=> v_structures(G)!=v_structures(C) fires for a PURELY SKELETAL reason)")
print(f"  mean #v-structures in C: {nv.mean():.2f}")
# S-loc composition
lp=np.array([r['loc_pool'] for r in res]); lpu=np.array([r['loc_pure'] for r in res])
lr=np.array([r['loc_rev'] for r in res]); lo=np.array([r['loc_reo'] for r in res])
lq=np.array([r['loc_query'] for r in res])
print(f"  S-loc pool: mean size {lp.mean():.2f}; composition over the POOL: "
      f"pure_spurious {lpu.sum()/lp.sum():.4f}, reversal_of_true_edge {lr.sum()/lp.sum():.4f}, "
      f"reorient_cpdag_edge {lo.sum()/lp.sum():.4f}")
print(f"  S-loc: E[share of a draw that is the QUERY pair (X,Y) or (Y,X)] = {(lq/lp).mean():.4f}"
      f"  -> P(>=1 of 4 slots is a query-pair statement) approx {np.mean(1-(1-lq/lp)**4):.4f}")
print(f"  frac SCMs whose S-loc pool < 4 (cannot fill 4 slots w/o replacement): {(lp<4).mean():.4f}")
print(f"  O0 size mean {np.mean([r['O0'] for r in res]):.2f}")
